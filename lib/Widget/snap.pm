#------------------------------------------------------------------------
#----                        Widget::snap                            ---- 
#------------------------------------------------------------------------
package Widget::snap;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use Bio::Search::Hit::PhatHit::snap;
use Bio::Search::HSP::PhatHSP::snap;
use PhatHit_utils;
use IPC::Open3;
use Symbol;
use FastaSeq;
use Error qw(:try);
use Error::Simple;

@ISA = qw(
	Widget
       );

my $OPT_F; # global varible for maker clean up
my $LOG; # global varible for maker runlog
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
    my $class  = shift;
    my @args   = @_;
    
    my $self = $class->SUPER::new(@args);
    
    bless ($self, $class);
    return $self;
}

#------------------------------------------------------------------------
sub prep_for_genefinder {
    my $seq    = shift;
    my $mia    = shift;
    my $set    = shift;
    my $flank  = shift;
    my $seq_id = shift;
    
    my ($span,
	$strand,
	$p_set_coors,
	$n_set_coors,
	$i_set_coors) = process_hints($seq, $mia, $set);
    
    my $p = Shadower::getPieces($seq, $span, $flank);
    my $final_seq  = $p->[0]->{piece};
    my $offset = $p->[0]->{b} - 1;
    my $i_flank    = 2;
    
    my $xdef = get_xdef($seq,
			$p_set_coors,
			$n_set_coors,
			$i_set_coors,
			$strand == 1 ? '+' : '-',
			0.2, # Coding
			1000, # intron
			$offset,
			$i_flank,
			);
    
    return (\$final_seq, $strand, $offset, $xdef);
}
#------------------------------------------------------------------------
sub process_hints {
    my $seq    = shift;
    my $mia    = shift;
    my $set    = shift;
    
    my $gomiph = $set->{gomiph}    || [];
    my $models = $set->{gomod}     || [];
    my $alts   = $set->{alt_ests}  || [];
    my $preds  = $set->{preds}     || [];
    my $all_p  = $set->{all_preds} || [];
    my $ests   = $set->{ests}      || [];
    my $fusion = $set->{fusion}    || [];
    my @t_data;
    
    push(@t_data, @{$gomiph});
    push(@t_data, @{$preds});
    push(@t_data, @{$alts});
    push(@t_data, @{$models});
    push(@t_data, @{$ests});
    
    #get span
    my $plus  = 0;
    my $minus = 0;
    my $least;
    my $most;
    foreach my $hit (@t_data){
	foreach my $hsp ($hit->hsps()){
	    my $s = $hsp->start('query');
	    my $e = $hsp->end('query');
	    
	    $least = $s if !defined($least) || $s < $least;
	    $most  = $e if !defined($most)  || $e > $most;
	    
	    if ($hsp->strand('query') == 1) {
		$plus++;
	    }
	    else {
		$minus++;
	    }
	}
    }
    
    #instantiate seq position array
    # 0 => empty,
    # 1 => exon,
    # 2 => reserved_exon,
    # 3 => CDS,
    # 4 => reserved_CDS,
    # -1 => intron,
    # -2 => reserved_intron
    my @b_seq = map {0} ($least..$most);
    my @span_coors = [$least, $most];
    my $strand = $plus > $minus ? 1 : -1;
    my $offset = $least;
    
    #identify all possible confirmed splice sites (compared later to all preds)
    my %splices;
    foreach my $e (@$ests, @$gomiph, @$alts){
	next if($e->num_hsps() == 1);
	my @hsps = sort {$a->start('query') <=> $b->start('query')} $e->hsps;
	for(my $i = 0; $i < @hsps - 1; $i++){
	    my $hsp1 = $hsps[$i];
	    my $hsp2 = $hsps[$i+1];
	    
	    next if($hsp1->end('query') >= $hsp2->start('query'));
	    
	    $splices{start}{$hsp2->start('query')}++; #exon end
	    $splices{end}{$hsp1->end('query')}++; #exon start
	    $splices{exact}{$hsp1->end('query')}{$hsp2->start('query')}++; #intron (bordering exons)
	}
    }

    #use $mia to infer alt splicing (keep all structure from mia)
    my $fusion_skip = grep {$mia->algorithm eq $_->algorithm} @$fusion if($mia && $fusion);
    if($mia && !$fusion_skip){
	my @hsps = sort {$a->start('query') <=> $b->start('query')} $mia->hsps;
	for(my $i = 0; $i < @hsps - 1; $i++){
	    my $hsp1 = $hsps[$i];
	    my $hsp2 = $hsps[$i+1];
	    
	    my $s1 = $hsp1->start('query');
	    my $e1 = $hsp1->end('query');
	    my $s2 = $hsp2->start('query');
	    my $e2 = $hsp2->end('query');
	    
	    #first exon given lower precedence
	    my $val = ($i == 0) ? 1 : 2;
	    my $s = $s1 - $offset;
	    my $e = $e1 - $offset;
	    @b_seq[$s..$e] = map {$val} ($s..$e); #overrides all

	    #last exon given lower precedence
	    if($i + 1 == @hsps - 1){
		$s = $s2 - $offset;
		$e = $e2 - $offset;
		@b_seq[$s..$e] = map {1} ($s..$e); #overrides all  
	    }
	    
	    next if($e1 >= $s2); #no gap so skip
	    
	    #intron verify
	    my $si = $e1 + 1; #fix for inron coordiantes
	    my $ei = $s2 - 1; #fix for intron coordiantes
	    
	    $si -= $offset;
	    $ei -= $offset;
	    @b_seq[$si..$ei] = map {-2} ($si..$ei); #overrides all
	}
    }
    
    #map introns (and unknown coding status exons)
    my %seen; #skip redundant
    foreach my $e (@{$ests}){
	#map entire hit space (defines intron once exons are mapped)
	if($e->num_hsps() > 1){
	    my $s = $e->start('query') - $offset;
	    my $e = $e->end('query') - $offset;
	    
	    foreach my $i ($s..$e){
		$b_seq[$i] = -1 if($b_seq[$i] == 0);
	    }
	}
	
	#map exons
	foreach my $hsp ($e->hsps){
	    my $s = $hsp->start('query');
	    my $e = $hsp->end('query');
	    
	    next if($seen{$s}{$e}); #skip redundant HSPs
	    $seen{$s}{$e}++;
	    
	    $s -= $offset;
	    $e -= $offset;
	    
	    foreach my $i ($s..$e){
		$b_seq[$i] = 1 if($b_seq[$i] == 0 || $b_seq[$i] == -1);
	    }
	}
    }
    
    #map protein HSPs
    %seen = (); #reset, skip redundant
    foreach my $c (@{PhatHit_utils::get_hsp_coors($gomiph, 'query')}){
	my $s = $c->[0];
	my $e = $c->[1];
	($s, $e) = ($e, $s) if($s > $e);
	
	next if($seen{$s}{$e}); #skip redundant HSPs
	$seen{$s}{$e}++;
	
	$s -= $offset;
	$e -= $offset;
	foreach my $i ($s..$e){
	    if($b_seq[$i] == 0 || $b_seq[$i] == 1){ #cannot overide EST inferred intron
		$b_seq[$i] = 3;
	    }
	    elsif($b_seq[$i] == 2){
		$b_seq[$i] = 4;
	    }
	}
    }
    
    #add EST to protein CDS if spliced and all orf
    my $tM = new CGL::TranslationMachine();
    foreach my $e (@$ests){
	next if($e->num_hsps() == 1);
	my ($lAq) = $e->getLengths();
	next if($lAq < 300); #orf of 300 required (same as most prokaryotic gene finders)
	
	my $t_seq  = maker::auto_annotator::get_transcript_seq($e, $seq);
	my ($p_seq, $poffset) = $tM->longest_translation($t_seq);
	my $end = $poffset + (length($p_seq) * 3 + 1);
	
	if($poffset < 3 && length($t_seq) - (3 * length($p_seq) + $poffset) < 3 && $p_seq !~ /X{10}/){
	    foreach my $hsp ($e->hsps){
		my $s = $hsp->start('query');
		my $e = $hsp->end('query');
		
		next if($seen{$s}{$e}); #skip redundant HSPs
		$seen{$s}{$e}++;
		
		$s -= $offset;
		$e -= $offset;
		
		foreach my $i ($s..$e){
		    if($b_seq[$i] == 0 || $b_seq[$i] == 1 || $b_seq[$i] == -1){
			$b_seq[$i] = 3;
		    }
		    elsif($b_seq[$i] == 2){
			$b_seq[$i] = 4;
		    }
		}
	    }	
	    next;
	}
	
	#if not all orf test internal HSPs
	next if($e->num_hsps() < 2);
	next if((3 * length($p_seq)) / length($t_seq) < 0.5); #skip unless orf at least 50% of EST
	my $pos = 0;
	foreach my $hsp (@{PhatHit_utils::sort_hits($e, 'query')}){
	    my $B = $hsp->start('query');
	    my $E = $hsp->end('query');
	    my $L = abs($E - $B) + 1;
	    $pos += $L;
	    
	    next unless(($offset < ($pos - $L) + 1) && ($pos < $end)); #skip not fully orf exons
	    next if($hsp->start('query') == $e->start('query')); #skip first HSP
	    next if($hsp->end('query') == $e->end('query')); #skip last HSP
	    next if($L < 3); #can't be translated
	    next if($seen{$B}{$E}); #skip redundant HSPs
	    $seen{$B}{$E}++;
	    
	    my $piece = ($e->strand('query') == 1) ?
		substr_o($seq, $B-1, $L) : Fasta::revComp(substr_o($seq, $B-1, $L));	    
	    
	    my ($p_seq, $poffset) = $tM->longest_translation($piece);
	    
	    if($poffset < 3 && $L - (3 * length($p_seq) + $poffset) < 3 && $p_seq !~ /X{10}/){
		my $s = $B - $offset;
		my $e = $E - $offset;
		
		foreach my $i ($s..$e){
		    if($b_seq[$i] == 0 || $b_seq[$i] == 1 || $b_seq[$i] == -1){
			$b_seq[$i] = 3;
		    }
		    elsif($b_seq[$i] == 2){
			$b_seq[$i] = 4;
		    }
		}
	    }
	}
    }
    
    #get likely CDS from splice site crossing ESTs
    foreach my $c (@{PhatHit_utils::splice_infer_exon_coors($ests, $seq)}){
	my $s = $c->[0];
	my $e = $c->[1];
	
	next if($seen{$s}{$e}); #skip redundant HSPs
	$seen{$s}{$e}++;
	
	$s -= $offset;
	$e -= $offset;
	
	foreach my $i ($s..$e){
	    if($b_seq[$i] == 0 || $b_seq[$i] == 1 || $b_seq[$i] == -1){
		$b_seq[$i] = 3;
	    }
	    elsif($b_seq[$i] == 2){
		$b_seq[$i] = 4;
	    }
	}
    }
    
    #get intron info from protein2genome hits
    foreach my $p (@{maker::auto_annotator::get_selected_types($gomiph,'protein2genome')}){
	my $s = $p->start('query') - $offset;
	my $e = $p->end('query') - $offset;
	
	foreach my $i ($s..$e){
	    $b_seq[$i] = -1 if($b_seq[$i] == 0);
	}
    }
    
    #get dually confirmed ab-initio and EST CDS/introns
    my @intron_set;
    my @exon_set;
    foreach my $p (@$all_p){
	next if($p->num_hsps() == 1);
	my @hsps = sort {$a->start('query') <=> $b->start('query')} $p->hsps;
	for(my $i = 0; $i < @hsps - 1; $i++){
	    my $hsp1 = $hsps[$i];
	    my $hsp2 = $hsps[$i+1];
	    
	    my $s1 = $hsp1->start('query');
	    my $e1 = $hsp1->end('query');
	    my $s2 = $hsp2->start('query');
	    my $e2 = $hsp2->end('query');
	    
	    next if($e1 >= $s2); #no gap so skip
	    
	    #intron verify
	    if($splices{exact}{$e1}{$s2}){
		my $si = $e1 + 1; #fix for inron coordiantes
		my $ei = $s2 - 1; #fix for intron coordiantes
		push(@intron_set, [$si, $ei]);
		
		$si -= $offset;
		$ei -= $offset;
		
		#make intron (won't override dually verified CDS)
		foreach my $i ($si..$ei){
		    $b_seq[$i] = -1 if($b_seq[$i] != 2 && $b_seq[$i] != 4 && $b_seq[$i] != -2);
		}
	    }
	    
	    #CDS verify
	    if($i == 0 && $splices{start}{$s1} && $splices{end}{$e1}){
		my $s = $s1 - $offset;
		my $e = $e1 - $offset;
		push(@exon_set, [$s1, $e1]);
		
		#overides most other hints to force dually verified CDS
		foreach my $i ($s..$e){
		    $b_seq[$i] = 4 if($b_seq[$i] != -2);
		}
	    }
	    elsif($splices{end}{$e1} && ! $splices{start}{$s1}){
		my $s = $s1 - $offset;
		my $e = $e1 - $offset;
		foreach my $i (reverse($s..$e)){
		    if($b_seq[$i] == 1){ #upgrade exon regions to CDS
			$b_seq[$i] = 3;
		    }
		    elsif($b_seq[$i] == 2){
			$b_seq[$i] = 4;
		    }
		    else{ #stop on empty or intron
			push(@exon_set, [$i+1+$offset,$e1]) if($i != $e);
			last;
		    }
		}
	    }
	    
	    if($splices{start}{$s2} && $splices{end}{$e2}){
		my $s = $s2 - $offset;
		my $e = $e2 - $offset;
		push(@exon_set, [$s2, $e2]);
		
		#overides most other hints to force dually verified CDS
		foreach my $i ($s..$e){
		    $b_seq[$i] = 4 if($b_seq[$i] != -2);
		}
	    }
	    elsif($splices{start}{$s2} && ! $splices{end}{$e2}){
		my $s = $s2 - $offset;
		my $e = $e2 - $offset;
		foreach my $i ($s..$e){
		    if($b_seq[$i] == 1){ #upgrade exon regions to CDS
			$b_seq[$i] = 3;
		    }
		    elsif($b_seq[$i] == 2){
			$b_seq[$i] = 4;
		    }
		    else{ #stop on empty or intron
			push(@exon_set, [$s2,($i-1)+$offset]) if($i != $s);
			last;
		    }
		}
	    }
	}
    }
    
    #get coordinates
    my @n_set_coors;
    my @p_set_coors;
    my @i_set_coors;
    my $nfirst;
    my $pfirst;
    my $ifirst;
    for (my $i = 0; $i < @b_seq; $i++){
	$nfirst = $i if(!defined($nfirst) && $b_seq[$i] > 0);
	$pfirst = $i if(!defined($pfirst) && $b_seq[$i] > 2);
	$ifirst = $i if(!defined($ifirst) && $b_seq[$i] < 0);
	
	if($b_seq[$i] <= 0 && defined($nfirst)){
	    push(@n_set_coors, [$nfirst+$offset, $i-1+$offset]);
	    $nfirst = undef;
	}
	if($b_seq[$i] <= 2 && defined($pfirst)){
	    push(@p_set_coors, [$pfirst+$offset, $i-1+$offset]);
	    $pfirst = undef;
	}
	if($b_seq[$i] >= 0 && defined($ifirst)){
	    push(@i_set_coors, [$ifirst+$offset, $i-1+$offset]);
	    $ifirst = undef;
	}
    }
    
    my $j = @b_seq;
    push(@n_set_coors, [$nfirst+$offset, $j-1+$offset]) if(defined($nfirst));
    push(@p_set_coors, [$pfirst+$offset, $j-1+$offset]) if(defined($pfirst));
    push(@i_set_coors, [$ifirst+$offset, $j-1+$offset]) if(defined($ifirst));
    
    return (\@span_coors, $strand, \@p_set_coors, \@n_set_coors, \@i_set_coors);
}
#------------------------------------------------------------------------
sub get_pred_shot {
        my $seq           = shift;
        my $def           = shift;
        my $the_void      = shift;
        my $mia           = shift;
        my $set           = shift;
	my $set_id        = shift;
        my $snap_flank    = shift;
        my $snap_command  = shift;
        my $hmm           = shift;
	   $OPT_F         = shift;
	   $LOG           = shift;


        my ($shadow_seq, $strand, $offset, $xdef) =
            prep_for_genefinder($seq, $mia, $set, $snap_flank);
	my $id = $set->{c_id}.'_'.$set_id;
	my $end = $offset + length($$shadow_seq);
        my $shadow_fasta = Fasta::toFasta($def." $id offset:$offset",
                                          $shadow_seq,
                                         );

        my ($exe, $param) = $snap_command =~ /(\S+)\s+(\S+)/;

        my $alt_snap_command;
        if ($strand == 1){
                $alt_snap_command = $exe.' -plus ';
        }
        else {
                 $alt_snap_command = $exe.' -minus ';
        }

        $alt_snap_command .= $param;

        my $gene_preds = snap($shadow_fasta,
			      $the_void,
			      $id,
			      $strand,
			      $offset,
			      $end,
			      $xdef,
			      $alt_snap_command,
			      $hmm
			     );


        return ($gene_preds, $strand);

}
#------------------------------------------------------------------------
sub snap {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
        my $strand     = shift;
        my $offset     = shift;
	my $end        = shift;
        my $xdef       = shift;
        my $command    = shift;
        my $hmm        = shift;

        my $snap_keepers = [];
	my ($hmm_name) = $hmm =~ /([^\:\/]+)(\:[^\:\/]+)?$/;

	my $tmp = GI::get_global_temp();
	my $rank = GI::RANK();
	my $t_dir = "$tmp/$rank";

        my $file_name = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.snap.fasta";
        my $xdef_file = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.snap";
        my $o_file    = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.snap";
        my $backup    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.snap";
	
	my $run = $command;
	$run .= " -xdef $xdef_file ";
	$run .= " $file_name";
	$run .= " > $o_file";

	$LOG->add_entry("STARTED", $backup, "") if(defined $LOG);

        if (-f $backup && ! $OPT_F){
                print STDERR "re reading snap report.\n" unless $main::quiet;
                print STDERR "$backup\n" unless $main::quiet;
		$o_file = $backup;
        }
        else {
                print STDERR "running  snap.\n" unless $main::quiet;
		write_xdef_file($xdef, $xdef_file);
		FastaFile::writeFile(\$fasta, $file_name);
		my $w = new Widget::snap();
                $w->run($run);
        }
	unlink($xdef_file) if(-f $xdef_file);
	unlink($file_name) if(-f $file_name);

	my %params;
        $params{min_exon_score}  = -100;
        $params{min_gene_score}  = -20;
	
	my $keepers;
	try { #make sure it parse correctly
	    $keepers = parse($o_file,
			     \%params,
			     $fasta);
	    #File::Copy::copy($o_file, $backup) unless(-f $backup); #temp
	    unlink($o_file); #temp
	}
	catch Error::Simple with {
	    my $E = shift;

	    unlink($o_file);

	    #retry predictor in different location and parse again
	    $file_name = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.snap.fasta";
	    $xdef_file = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.snap";
	    $o_file    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.snap";
	    my $run = $command;
	    $run .= " -xdef $xdef_file ";
	    $run .= " $file_name";
	    $run .= " > $o_file";
	    write_xdef_file($xdef, $xdef_file);
	    FastaFile::writeFile(\$fasta, $file_name);
	    my $w = new Widget::snap();
	    $w->run($run);
	    unlink($xdef_file) if(-f $xdef_file);
	    unlink($file_name) if(-f $file_name);
	    $keepers = parse($o_file,
			     \%params,
			     $fasta);
	    unlink($o_file);
	};

	$LOG->add_entry("FINISHED", $backup, "") if(defined $LOG);
	
	PhatHit_utils::add_offset($keepers,
				  $offset,
                                  );
	
        return $keepers;

}
#------------------------------------------------------------------------------
sub run {
	my $self    = shift;
	my $command = shift;

	if (defined($command)){
	    $self->print_command($command);
	    my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
	    my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	    {
		local $/ = \1;
		while (my $line = <$CHLD_ERR>){
		    print STDERR $line unless($main::quiet);
		}
	    }
	    waitpid $pid, 0;
	    die "ERROR: Snap failed\n" if $? != 0;
	}
	else {
	    die "you must give Widget::snap a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
        shift @_ if($_[0] =~ /^Widget\:\:/); #ignore 
        my $report = shift;
        my $params = shift;
	my $q_file = shift;
	
	my $def;
	my $q_seq;

	if(ref($q_file) eq 'FastaChunk'){ #object not scalar
	    $q_seq = $q_file->seq_w_flank;
	    $def = $q_file->def_w_flank;
	}
	elsif($q_file =~ /^>/){
            $def     = Fasta::getDef(\$q_file);
            $q_seq   = Fasta::getSeqRef(\$q_file);
        }
        else{
	    my $index = GI::build_fasta_index($q_file);
	    my ($id) = $index->get_all_ids;
	    $def = ">$id";
            $q_seq   = $index->get_Seq_by_id($id);
        }

        my ($q_name)  = $def =~ /^>(.+)/;

	my $fh = new FileHandle();
	   $fh->open($report);

	my $header = <$fh>; # shift off the header

        #checks if snap really finished
        unless($header){
            unlink($report);
            die "ERROR: The file $report appears to be incomplete\n".
		"MAKER will need to delete the file, before trying again\n\n";
        }

	my %g;
	my $i = -1;
	while (my $line = <$fh>){
		chomp($line);
		my @stuff = split(/\s+/, $line);

		$i = 0  if  !defined($g{$stuff[8]});

		my $b = $stuff[1];
		my $e = $stuff[2];
		#($b, $e) = ($e, $b) if $stuff[3] eq '-';

		$g{$stuff[8]}[$i]{type}     = $stuff[0];
		$g{$stuff[8]}[$i]{b}        = $stuff[1];
		$g{$stuff[8]}[$i]{e}        = $stuff[2];
		$g{$stuff[8]}[$i]{strand}   = $stuff[3] eq '+' ? 1 : -1;
		$g{$stuff[8]}[$i]{score}    = $stuff[4];
		$g{$stuff[8]}[$i]{five}     = $stuff[5];
		$g{$stuff[8]}[$i]{three}    = $stuff[6];
		$g{$stuff[8]}[$i]{frame}    = $stuff[7];	
		$i++;
	}
	$fh->close();
	my $phat_hits = load_phat_hits($q_name, $q_seq, \%g, $params);
	return $phat_hits;
}
#-------------------------------------------------------------------------------
sub get_xdef {
	my $seq     = shift;
	my $p_coors = shift;
	my $n_coors = shift;
	my $i_coors = shift;
	my $s       = shift;
	my $c_u     = shift;
	my $i_u     = shift;
	my $offset  = shift;
	my $i_flank = shift;	

	my @index;

        my @xdef;
	push(@xdef, ">xdef file offset: $offset i_flank:$i_flank");
        foreach my $p (@$p_coors){
                my $c_b = $p->[0] - $offset;
                my $c_e = $p->[1] - $offset;
		my $l = "Coding\t".$c_b."\t".$c_e;
		$l .= "\t$s\t$c_u\t\.\t\.\t\.\tADJ";
                push(@xdef, $l);
        }

	return \@xdef if(!@$i_coors);

        my @intron;

	foreach my $i (@$i_coors){
                my $i_b = ($i->[0] - $offset) + ($i_flank-1);
                my $i_e = ($i->[1] - $offset) - ($i_flank-1);
	
		next if abs($i_b - $i_e) < 2*$i_flank;
		next if abs($i_b - $i_e) < 25;

		my $l = "Intron\t".$i_b."\t".$i_e;
                $l .= "\t$s\t$i_u\t\.\t\.\t\.\tADJ";

		push(@xdef, $l);

		my $c = "Coding\t".$i_b."\t".$i_e;
		$c .= "\t$s\t-100\t\.\t\.\t\.\tADJ";
		push(@xdef, $c);
        }

	return \@xdef;
}
#-------------------------------------------------------------------------------
sub adjust_xdef {
	my $xdef       = shift;
	my $offset     = shift;
	my $chunk_size = shift;
	my @new_xdef;
	my $i = 0;

	my $end = $offset + $chunk_size;

	foreach my $l (@{$xdef}){
		if ($i == 0){
			push(@new_xdef, $l);
		}
		else {
			my @fields = split (/\s+/, $l);
			if ($fields[1] >=  $offset && $fields[2] <= $end){
				push(@new_xdef, $l);
			}
			else {
			} 
		}
			
		$i++;
	}
	return \@new_xdef;
}
#-------------------------------------------------------------------------------
sub write_xdef_file {
	my $xdef = shift;
	my $loc  = shift;

	my $fh = new FileHandle();
	   $fh->open(">$loc") || die "couldn't open $loc";	

	foreach my $l (@{$xdef}){
		print $fh $l."\n";
	}
	$fh->close();

}
#-------------------------------------------------------------------------------
sub get_exon_seq {
	my $exon  = shift;
	my $q_seq = shift;

       my $e_b = $exon->{b};
       my $e_e = $exon->{e};
        
        my $length = abs($e_e - $e_b) + 1;
        
        ($e_b, $e_e) = ($e_e, $e_b) if($e_b > $e_e);
	
        my $e_seq = substr_o($q_seq, $e_b - 1, $length);

        $e_seq = Fasta::revComp($e_seq) if $exon->{strand} == -1;

	return $e_seq;
}
#-------------------------------------------------------------------------------
sub total_score {
	my $g = shift;

	my $total = 0;
	foreach my $exon (@{$g}){
		$total += $exon->{score};
	}
	return $total;
}
#-------------------------------------------------------------------------------
sub load_phat_hits {
	my $q_name  = shift;
	my $q_seq   = shift;
	my $g       = shift;
	my $params  = shift;
	
	my $q_len = length_o($q_seq);
	
	my @keepers;
	foreach my $gene (keys %{$g}){
	
		my $total_score = total_score($g->{$gene});

		#build hit
		my %hsps;
		my $i = 0;
		my $f = new Bio::Search::Hit::PhatHit::snap('-name'         => $gene,
							    '-description'  => 'NA',
							    '-algorithm'    => 'snap',
							    '-length'       => $q_len,
							    '-score'        => $total_score,
							    );

                $f->queryLength($q_len);

		#added 3/19/2009
		#check for single and double base pair overhangs
		@{$g->{$gene}} = sort {$a->{b} <=> $b->{b}} @{$g->{$gene}};
		my $length = 0;
		foreach my $exon (@{$g->{$gene}}){
		    $length += abs($exon->{e} - $exon->{b}) + 1;
		}

		my $overhang = $length % 3;
		if($overhang != 0){
		    if($g->{$gene}->[0]->{strand} == 1){
			my $last = $g->{$gene}->[-1];
			my $l_length = abs($last->{e} - $last->{b}) + 1;

			while($l_length <= $overhang){
			    pop(@{$g->{$gene}});
			    $overhang -= $l_length;
			    $last = $g->{$gene}->[-1];
			    $l_length = abs($last->{e} - $last->{b}) + 1;
			}

			$last->{e} -= $overhang;
		    }
		    elsif($g->{$gene}->[0]->{strand} == -1){
			my $last = $g->{$gene}->[0];
                        my $l_length = abs($last->{e} - $last->{b}) + 1;

                        while($l_length <= $overhang){
                            shift(@{$g->{$gene}});
                            $overhang -= $l_length;
                            $last = $g->{$gene}->[0];
                            $l_length = abs($last->{e} - $last->{b}) + 1;
                        }

			$last->{b} += $overhang;
		    }
		    else{
			die "FATAL: No exon strand in Widget::snap\n";
		    }
		}
		
		#build hsps
		my $hit_start = 1;
		foreach my $exon (@{$g->{$gene}}){
			my @args;
			my $exon_seq = get_exon_seq($exon, $q_seq); # revcomped!
			my $hit_end = abs($exon->{e} - $exon->{b}) + $hit_start;

 			push(@args, '-query_start');
        		push(@args, $exon->{b});

       			push(@args, '-query_seq');
       			push(@args, $exon_seq);

        		push(@args, '-score');
        		push(@args, $exon->{score});
       			push(@args, '-homology_seq');
       			push(@args, $exon_seq);

        		push(@args, '-hit_start');
        		push(@args, $hit_start);

       			push(@args, '-hit_seq');
       			push(@args, $exon_seq);

        		push(@args, '-hsp_length');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-identical');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_length');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-query_name');
        		push(@args, $q_name);

        		push(@args, '-algorithm');
        		push(@args, 'SNAP');

        		push(@args, '-bits');
        		push(@args, 'NA');

        		push(@args, '-evalue');
        		push(@args, 'NA');

        		push(@args, '-pvalue');
        		push(@args, 'NA');

        		push(@args, '-query_length');
        		push(@args, $q_len);

        		push(@args, '-query_end');
        		push(@args, $exon->{e});

        		push(@args, '-conserved');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_name');
        		push(@args, $gene.":".$i.":".$exon->{type});

        		push(@args, '-hit_end');
        		push(@args, $hit_end);

        		push(@args, '-query_gaps');
        		push(@args, 0);

        		push(@args, '-hit_gaps');
        		push(@args, 0);

        		my $hsp = new Bio::Search::HSP::PhatHSP::snap(@args);
        		   $hsp->queryName($q_name);
        		#-------------------------------------------------
        		# setting strand because bioperl is all messed up!
        		#------------------------------------------------
			if ($exon->{strand} == 1 ){
				$hsp->{_strand_hack}->{query} = 1;
				$hsp->{_strand_hack}->{hit}   = 1;
			}
			else {
				$hsp->{_strand_hack}->{query} = -1;
				$hsp->{_strand_hack}->{hit}   =  1;
			}
			$f->add_hsp($hsp);

			$hit_start += abs($exon->{e} - $exon->{b}) + 1; 
			
			$i++;
		}

              	push(@keepers, $f);

        }

        my $final = keepers(\@keepers, $params);

	#set start and stop coordinates for codons
	foreach my $f (@$final){
	    my $seq = maker::auto_annotator::get_transcript_seq($f,$q_seq);
	    maker::auto_annotator::get_translation_seq($seq, $f, 1);
	}

	return $final;
}
#-------------------------------------------------------------------------------
sub keepers {
    my $genes   = shift;
    my $params = shift; 
    
    my $start = @{$genes};
    my @keepers;
    foreach my $phat_hit (@{$genes}){
	my @hsps;
	my $total_score = 0;
	while(my $hsp = $phat_hit->next_hsp) {
	    push(@hsps, $hsp)
		if(!defined($params->{min_exon_score}) ||
		   $hsp->score() > $params->{min_exon_score});
	    
	    $total_score += $hsp->score();
	}
	$phat_hit->hsps(\@hsps);
	push(@keepers, $phat_hit) 
	    if ($phat_hit->hsps() &&
		( !defined($params->{min_gene_score}) || 
		  $total_score > $params->{min_gene_score}));	
	
    }
    my $end     = @keepers;
    my $deleted = $start - $end;
    print STDERR "deleted:$deleted genes\n" unless $main::quiet;
    
    return \@keepers;
}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
    my ($self, $arg) = @_;
    
    my $caller = caller();
    use vars qw($AUTOLOAD);
    my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
    $call =~/DESTROY/ && return;
    
    print STDERR "Widget::snap::AutoLoader called for: ",
    "\$self->$call","()\n";
    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
    
    if (defined($arg)){
	$self->{$call} = $arg;
    }
    else {
	return $self->{$call};
    }
}
#------------------------------------------------------------------------

1;


