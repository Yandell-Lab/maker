#------------------------------------------------------------------------
#----             Widget::fgenesh			             ---- 
#------------------------------------------------------------------------
package Widget::fgenesh;
use strict;
use lib '~/maker/lib';
use lib '/data1/hao/projects/MAKER-fgenesh/lib';

use vars qw(@ISA $TMP);
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use PhatHit_utils;
use IPC::Open3;
use Bio::Search::Hit::PhatHit::fgenesh;
use Bio::Search::HSP::PhatHSP::fgenesh;

@ISA = qw(
	Widget
       );
my $OPT_F;
my $LOG;
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

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
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
sub prep_for_genefinder {
        my $seq    = shift;
        my $set    = shift;
        my $flank  = shift;

        my $gomiph = $set->{gomiph};
        my $mia    = $set->{mia};
        my $models = $set->{gomod};
	my $alts   = $set->{alt_ests};
        my $preds  = $set->{preds};
        my $ests   = $set->{ests};
        my @t_data;

        push(@t_data, @{$gomiph})  if defined($gomiph);
        push(@t_data, @{$preds})   if defined($preds);
        push(@t_data, $mia)        if defined($mia);
        push(@t_data, @{$alts})     if defined($alts);
        push(@t_data, @{$models})   if defined($models);
        push(@t_data, @{$ests})    if defined($ests);

        my $p_set_coors = PhatHit_utils::get_hsp_coors($gomiph, 'query');

        my $n_set_coors =
        defined($ests) ? PhatHit_utils::get_hsp_coors($ests, 'query')
                      : [];

        my @coors;
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
        my @span_coors = [$least, $most];

        my $p = Shadower::getPieces($seq, \@span_coors, $flank);

        my $final_seq = $p->[0]->{piece};

        my $offset    = $p->[0]->{b} - 1;

        my $strand = $plus > $minus ? 1 : -1;

        my $i_flank = 2;

        my $xdef = get_xdef($seq,
                            $p_set_coors,
                            $n_set_coors,
                            $strand == 1 ? '+' : '-',
                            $offset,
                            $i_flank,
                           );


        return (\$final_seq, $strand, $offset, $xdef);
}

#-------------------------------------------------------------------------------
sub run {
        my $self    = shift;
        my $command = shift;
        
	if (defined($command)){
	   $self->print_command($command);
	   my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, $command);
	   local $/ = \1; #read in everytime a byte becomes available
	   while (my $line = <CHLD_ERR>){
	      print STDERR $line unless($main::quiet);
	   }
	   waitpid $pid, 0;
	   #die "ERROR: FgenesH failed\n" if $? != 0;
	}
	else {
	   die "you must give Widget::fgenesh a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
sub parse {
        my $report = shift;
        my $params = shift;
        my $q_file = shift;

        my $iterator = new Iterator::Fasta($q_file);
        my $fasta = $iterator->nextFasta();

        my $def     = Fasta::getDef($fasta);
        my $q_seq   = Fasta::getSeqRef($fasta);
        my $q_name  = Fasta::def2SeqID($def);;
        
        my $fh = new FileHandle();
           $fh->open($report);

	my %g;
	my @content;

	## get rid of the first nine lines of the output.
	foreach (1..6) {
		<$fh>;
	}

	my $line = <$fh>;

	#($g{official_score}) = $line =~ /Score:(.+)\n/;
	# the score that fgenesh assigned to the whole seq.

	<$fh>;

	while (my $line=<$fh>) {
		next if $line =~ /^\n/;
		push @content, $line;
		last if $line =~ /Predicted protein/;
	}
	pop @content;


	foreach my $line (@content) {
		chomp $line;

		my %item;

		my (@stuff) = $line =~ /(\S+)/g;
		my $gene_num = $stuff[0];
		my $gene_def = $q_name."-fgenesh.$gene_num";

		if ($stuff[2] eq 'TSS' ||$stuff[2] eq 'PolA') {
			$item{gene_num} = $stuff[0];
			$item{strand} = $stuff[1] eq '+' ? '1' : '-1';
			$item{type}   = $stuff[2];
			$item{start}  = $stuff[3];
			$item{score}  = $stuff[4];
		
		}

		else {
			$item{gene_num} = $stuff[0];
			$item{strand}   = $stuff[1] eq '+' ? '1' : '-1';
			$item{type}     = $stuff[3];
			$item{b}    = $stuff[4];
			$item{e}      = $stuff[6];
			$item{score}    = $stuff[7];
			$item{ORF_b}    = $stuff[8];
			$item{ORF_e}    = $stuff[10];
			$item{length}   = $stuff[11];
		}
		push @{$g{$gene_def}}, \%item;
	}
	
	$fh->close;
	my $phat_hits = load_phat_hits($q_name, $q_seq, \%g, $params);
	return $phat_hits;

}


#-------------------------------------------------------------------------------
sub total_score {
        my $g = shift;

        my $total = 0;
        foreach my $item (@{$g}){
                $total += $item->{score};
        }
        return $total;
}
#-------------------------------------------------------------------------------
sub get_pred_shot {
	my $seq           = shift;
        my $def           = shift;
        my $id            = shift;
        my $the_void      = shift;
        my $set           = shift;
        my $pred_flank    = shift;
        my $pred_command  = shift;
           $OPT_F         = shift;
	   $LOG           = shift;

        my ($shadow_seq, $strand, $offset, $xdef) =
            prep_for_genefinder($seq, $set, $pred_flank);

        my $shadow_fasta = Fasta::toFastaRef($def." $id offset:$offset",
                                          $shadow_seq,
                                         );

        my $gene_preds = fgenesh ($$shadow_fasta,
                                  $the_void,
                                  $id,
                                  $strand,
                                  $offset,
                                  $xdef,
                                  $pred_command,
                                 );


        return ($gene_preds, $strand);
}

#-------------------------------------------------------------------------------
sub fgenesh {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
        my $strand     = shift;
        my $offset     = shift;
        my $xdef       = shift;
        my $command    = shift;

	
	my $wrap = "$FindBin::Bin/../lib/Widget/fgenesh/fgenesh_wrap"; #fgenesh wrapper

	my $file_name = "$the_void/$seq_id\.auto_annotator\.$offset\.fgenesh.fasta";
        my $o_file    = "$the_void/$seq_id\.$offset\.auto_annotator\.fgenesh";

        my $xdef_file = "$the_void/$seq_id\.$offset\.auto_annotator\.xdef\.fgenesh";
                            
        write_xdef_file($xdef, $xdef_file) if defined $xdef;
                            
        FastaFile::writeFile(\$fasta, $file_name);
        $command = $wrap . " $command"; #prepend wrapper
	$command .= " $file_name";
	$command .= " -tmp $TMP";                        
        $command .= ' -exon_table:'.$xdef_file if -e $xdef_file;
                           
        $command .= " > $o_file";
        
	$LOG->add_entry("STARTED", $o_file, "") if(defined $LOG);

        if (-e $o_file && ! $OPT_F){
                print STDERR "re reading fgenesh report.\n"
                        unless $main::quiet;
                print STDERR "$o_file\n"
                        unless $main::quiet;
        } 
        else { 
                print STDERR "running fgenesh.\n" unless $main::quiet;
                my $w = new Widget::fgenesh();
                $w->run($command);
        }

	$LOG->add_entry("FINISHED", $o_file, "") if(defined $LOG);
        
        my %params;
           $params{min_exon_score}  = -100;
           $params{min_gene_score}  = -20;
        
        my $keepers = parse($o_file,
                           \%params,
                            $file_name,
                           );
        
        PhatHit_utils::add_offset($keepers,
                                  $offset,
                                  );
        
        return $keepers; 
}
#-------------------------------------------------------------------------------
sub get_exon_seq {                        
        my $exon  = shift;                
        my $q_seq = shift;                

        my $e_b = $exon->{b};
        my $e_e = $exon->{e};
                
        my $length = abs($e_e - $e_b) + 1;
                        
        ($e_b, $e_e) = ($e_e, $e_b) if $e_b > $e_e;
                
        my $e_seq = substr($$q_seq, $e_b - 1, $length);
                        
        $e_seq = Fasta::revCompRef($e_seq) if $exon->{strand} == -1;

        return $e_seq;  
}
#-------------------------------------------------------------------------------
sub load_phat_hits {
        my $q_name  = shift;
        my $q_seq   = shift;
        my $g       = shift;
        my $params  = shift;
        
        my $q_len = length($$q_seq);

        my @keepers; 
        foreach my $gene (keys %{$g}){

                my $total_score = total_score($g->{$gene});
                
                my %hsps;
                my $i = 0;
                my $f = Bio::Search::Hit::PhatHit::fgenesh->new('-name'        => $gene,
								'-description' => 'NA',
								'-algorithm'   => 'fgenesh',
								'-length'      => $q_len,
								'-score'       => $total_score,
							       );

                $f->queryLength($q_len);

		#added 3/19/2009
		#check for single and double base pair overhangs
                @{$g->{$gene}} = grep {$_->{type} !~ /TSS|PolA/} @{$g->{$gene}};
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

			## skip if the item is promoter or poly A sequence.
			next if $exon->{type} eq 'TSS' 
				|| $exon->{type} eq 'PolA';

                        my @args;
                        my $exon_seq = get_exon_seq($exon, $q_seq);
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
                        push(@args, 'fgenesh');

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
                        my $hsp = new Bio::Search::HSP::PhatHSP::fgenesh(@args);
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
        return keepers(\@keepers, $params);

}

#------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_xdef {
       my $seq       = shift;
        my $p_coors  = shift;
        my $n_coors = shift;
        my $s       = shift;
        my $offset  = shift;
        my $i_flank = shift;

        my $group = 'gb|bogus;';

        my $p_pieces = Shadower::getPieces($seq, $p_coors, 0);
        my $n_pieces = Shadower::getPieces($seq, $n_coors, 0);

        my @xdef;

	my $num = (scalar @{$p_pieces} ) + (scalar @{$n_pieces});
	push @xdef, "$num $s";
        for (my $i = 0; $i < @{$p_pieces}; $i++){
                my $p = $p_pieces->[$i];
                my $c_b = $p->{b} - $offset;
                my $c_e = $p->{e} - $offset;

		my $l = $c_b.' '.$c_e.' '.'1000';
                
                push(@xdef, $l);
        }
        
        
        return \@xdef if !defined($n_pieces->[1]);
        
        for (my $i = @{$n_pieces} - 1; $i > 0;  $i--){
                my $p_r = $n_pieces->[$i]; 
                my $p_l = $n_pieces->[$i - 1];
                my $i_b = ($p_l->{e} - $offset) + $i_flank;
                my $i_e = ($p_r->{b} - $offset) - $i_flank;
                
                next if abs($i_b - $i_e) < 2*$i_flank;
                next if abs($i_b - $i_e) < 25;

		my $l  = $i_b.' '.$i_e.' -1000';

                push(@xdef, $l);
        }

	# The next part also put ESTs as arbitrary exons.
=pod
        for (my $i = 0; $i < @{$n_pieces}; $i++){
                my $p = $n_pieces->[$i];
                my $e_b = $p->{b} - $offset;
                my $e_e = $p->{e} - $offset;

                my $l  = "$q_id\tEST-EXON\texonpart";
                   $l .= "\t".$e_b."\t".$e_e."\t"."1e-1000"."\t".$s;
                   $l .= "\t".'.'."\t".'group=n_peices;'."source=E";

                push(@xdef, $l);
        }
=cut
        return \@xdef;

}
#------------------------------------------------------------------------
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
                        if $hsp->score() > $params->{min_exon_score};
                        
                        $total_score += $hsp->score();
             }
             $phat_hit->hsps(\@hsps);
             push(@keepers, $phat_hit) 
             if ($phat_hit->hsps() && $total_score > $params->{min_gene_score});
        
        }
        my $end     = @keepers;
        my $deleted = $start - $end;
        print STDERR "deleted:$deleted genes\n" unless $main::quiet;
        
        return \@keepers;
}
#------------------------------------------------------------------------
1;


