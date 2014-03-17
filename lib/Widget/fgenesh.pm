#------------------------------------------------------------------------
#----             Widget::fgenesh			             ---- 
#------------------------------------------------------------------------
package Widget::fgenesh;
use strict;
use vars qw(@ISA);
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
use Symbol;
use FastaSeq;
use Error qw(:try);
use Error::Simple;

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
        my $mia    = shift;
        my $set    = shift;
        my $flank  = shift;

	my ($span,
            $strand,
            $p_set_coors,
            $n_set_coors,
	    $i_set_coors) = Widget::snap::process_hints($seq, $mia, $set);
	
        my $p = Shadower::getPieces($seq, $span, $flank);
        my $final_seq  = $p->[0]->{piece};
        my $offset = $p->[0]->{b} - 1;
        my $i_flank    = 2;

        my $xdef = get_xdef($seq,
                            $p_set_coors,
                            $n_set_coors,
                            $i_set_coors,
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
	   my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
	   my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	   {
	       local $/ = \1; #read in everytime a byte becomes available
	       while (my $line = <$CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
	       }
	   }
	   waitpid $pid, 0;
	   die "ERROR: FgenesH failed\n" if $? != 0;
	}
	else {
	   die "you must give Widget::fgenesh a command to run!\n";
	}
}
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
        my $the_void      = shift;
        my $mia           = shift;
        my $set           = shift;
	my $set_id        = shift;
        my $pred_flank    = shift;
        my $pred_command  = shift;
        my $hmm           = shift;
           $OPT_F         = shift;
	   $LOG           = shift;

        my ($shadow_seq,
	    $strand,
	    $offset,
	    $xdef) = prep_for_genefinder($seq, $mia, $set, $pred_flank);
        my $id = $set->{c_id}.'_'.$set_id;
        my $end = $offset + length($$shadow_seq);
        my $shadow_fasta = Fasta::toFastaRef($def." $id offset:$offset",
                                          $shadow_seq,
                                         );

        my $gene_preds = fgenesh ($$shadow_fasta,
                                  $the_void,
                                  $id,
                                  $strand,
                                  $offset,
				  $end,
                                  $xdef,
                                  $pred_command,
				  $hmm
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
	my $end        = shift;
        my $xdef       = shift;
        my $command    = shift;
        my $hmm        = shift;

	my ($hmm_name) = $hmm =~ /([^\:\/]+)(\:[^\:\/]+)?$/;
	my $wrap = "$FindBin::RealBin/../lib/Widget/fgenesh/fgenesh_wrap"; #fgenesh wrapper

        my $tmp = GI::get_global_temp();
        my $rank = GI::RANK();
        my $t_dir = "$tmp/$rank";

	my $file_name = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.fgenesh.fasta";
        my $xdef_file = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.fgenesh";
        my $o_file    = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.fgenesh";
        my $backup    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.fgenesh";
                            
        my $run = $wrap . " $command"; #prepend wrapper
	$run .= " $file_name";
	#$run .= " -tmp $TMP";                        
        $run .= ' -exon_table:'.$xdef_file if(defined $xdef);
        $run .= " > $o_file";
        
	$LOG->add_entry("STARTED", $backup, "") if(defined $LOG);

        if (-f $backup && ! $OPT_F){
                print STDERR "re reading fgenesh report.\n" unless $main::quiet;
                print STDERR "$backup\n" unless $main::quiet;
		$o_file = $backup;
        } 
        else { 
                print STDERR "running fgenesh.\n" unless $main::quiet;
		write_xdef_file($xdef, $xdef_file) if(defined $xdef);
		FastaFile::writeFile(\$fasta, $file_name);
                my $w = new Widget::fgenesh();
                $w->run($run);
		#File::Copy::copy($o_file, $backup) unless();
        }

	unlink($xdef_file) if(-f $xdef_file);
        unlink($file_name) if(-f $file_name);
 
        my %params;
           $params{min_exon_score}  = -100;
           $params{min_gene_score}  = -20;
        
	my $keepers;
        try { #make sure it parses correctly
            $keepers = parse($o_file,
                             \%params,
                             $fasta);
            unlink($o_file);
            #File::Copy::copy($o_file, $backup) unless();
        }
        catch Error::Simple with {
            my $E = shift;

            unlink($o_file);

            #retry predictor in different location and parse again
	    $file_name = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.fgenesh.fasta";
	    $xdef_file = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.fgenesh";
	    $o_file    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.fgenesh";
	    my $run = $wrap . " $command"; #prepend wrapper
	    $run .= " $file_name";
	    #$run .= " -tmp $TMP";
	    $run .= ' -exon_table:'.$xdef_file if(defined $xdef);
	    $run .= " > $o_file";
            write_xdef_file($xdef, $xdef_file) if(defined $xdef);
	    FastaFile::writeFile(\$fasta, $file_name);
            my $w = new Widget::fgenesh();
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
#-------------------------------------------------------------------------------
sub get_exon_seq {                        
        my $exon  = shift;                
        my $q_seq = shift;                

        my $e_b = $exon->{b};
        my $e_e = $exon->{e};
                
        my $length = abs($e_e - $e_b) + 1;
                        
        ($e_b, $e_e) = ($e_e, $e_b) if $e_b > $e_e;

        my $e_seq = substr_o($q_seq, $e_b - 1, $length);
                        
        $e_seq = Fasta::revComp($e_seq) if $exon->{strand} == -1;

        return $e_seq;  
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
                        die "FATAL: No exon strand in Widget::fgenesh\n";
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

	my $final = keepers(\@keepers, $params);

        #set start and stop coordinates for codons
        foreach my $f (@$final){
            my $seq = maker::auto_annotator::get_transcript_seq($f,$q_seq);
	    maker::auto_annotator::get_translation_seq($seq, $f, 1);
	}

	return $final;
}

#------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_xdef {
    my $seq       = shift;
    my $p_coors  = shift;
    my $n_coors = shift;
    my $i_coors = shift;
    my $s       = shift;
    my $offset  = shift;
    my $i_flank = shift;
    
    my $group = 'gb|bogus;';
    
    my @xdef;
    foreach my $p (@$p_coors){
	my $c_b = $p->[0] - $offset;
	my $c_e = $p->[1] - $offset;
	
	my $l = "$c_b $c_e 100";
	
	push(@xdef, $l);
    }
    
    foreach my $i (@$i_coors){
	my $i_b = ($i->[0] - $offset) + ($i_flank-1);
	my $i_e = ($i->[1] - $offset) - ($i_flank-1);
	
	next if abs($i_b - $i_e) < 2*$i_flank;
	next if abs($i_b - $i_e) < 25;
	
	my $l  = "$i_b $i_e -1000";
	
	push(@xdef, $l);
    }
    
    my $num = @xdef;
    unshift(@xdef, "$num $s");
    

=pod
    # The next part also put ESTs as arbitrary exons.
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
		    if(!defined($params->{min_exon_score}) ||
		       $hsp->score() > $params->{min_exon_score});
		
		$total_score += $hsp->score();
	    }
	    $phat_hit->hsps(\@hsps);
	    push(@keepers, $phat_hit)
		if ($phat_hit->hsps() &&
		    (!defined($params->{min_gene_score}) ||
		     $total_score > $params->{min_gene_score}));
	    
	}
        my $end     = @keepers;
        my $deleted = $start - $end;
        print STDERR "deleted:$deleted genes\n" unless $main::quiet;
        
        return \@keepers;
}
#------------------------------------------------------------------------
1;


