#------------------------------------------------------------------------
#----                        Widget::augustus                        ---- 
#------------------------------------------------------------------------
package Widget::augustus;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use Bio::Search::Hit::PhatHit::augustus;
use Bio::Search::HSP::PhatHSP::augustus;
use PhatHit_utils;
use IPC::Open3;
use FindBin;
use Symbol;
use FastaSeq;
use Error qw(:try);
use Error::Simple;

@ISA = qw(
	Widget
       );

my $OPT_F; # global option -f to force cleanup
my $LOG; #global varible for maker runlog
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

	my $q_id = Fasta::def2SeqID($def);
        my $id = $set->{c_id}.'_'.$set_id;
        my ($shadow_seq, $strand, $offset, $xdef) =
            prep_for_genefinder($seq, $mia, $set, $pred_flank, $id, $q_id);
        my $end = $offset + length($$shadow_seq);
        my $shadow_fasta = Fasta::toFasta($def." $id offset:$offset",
                                          $shadow_seq,
                                         );

        my $alt_aug_command;
        if ($strand == 1){
                $alt_aug_command = $pred_command.' --strand=forward';
        }
        else {
                 $alt_aug_command = $pred_command.' --strand=backward';
        }

        my $gene_preds = augustus($shadow_fasta,
                                  $the_void,
                                  $id,
                                  $strand,
                                  $offset,
				  $end,
                                  $xdef,
                                  $alt_aug_command,
				  $hmm
                                 );


        return ($gene_preds, $strand);
}
#------------------------------------------------------------------------
sub prep_for_genefinder {
        my $seq         = shift;
        my $mia         = shift;
        my $set         = shift;
        my $flank       = shift;
        my $seq_id      = shift;
        my $q_id        = shift;

	my ($span,
	    $strand,
	    $p_set_coors,
	    $n_set_coors,
	    $i_set_coors) = Widget::snap::process_hints($seq, $mia, $set);

        my $p = Shadower::getPieces($seq, $span, $flank);
        my $final_seq = $p->[0]->{piece};
        my $offset = $p->[0]->{b} - 1;
        my $i_flank = 2;

        my $xdef = get_xdef($seq,
			    $p_set_coors,
                            $n_set_coors,
                            $i_set_coors,
                            $strand == 1 ? '+' : '-',
                            $offset,
                            $i_flank,
                            $q_id,
                           );


        return (\$final_seq, $strand, $offset, $xdef);
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
#------------------------------------------------------------------------
sub augustus {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
        my $strand     = shift;
        my $offset     = shift;
	my $end        = shift;
        my $xdef       = shift;
        my $command    = shift;
        my $hmm        = shift;

        my $aug_keepers = [];
	my ($hmm_name) = $hmm =~ /([^\:\/]+)(\:[^\:\/]+)?$/;
        my $cfg_file = "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg";

        my $tmp = GI::get_global_temp();
        my $rank = GI::RANK();
        my $t_dir = "$tmp/$rank";

        my $file_name = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.augustus.fasta";
        my $xdef_file = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.augustus";
        my $o_file    = "$t_dir/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.augustus";
        my $backup    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.augustus";

	my $run = $command;
	$run .= ' --UTR=off';
	$run .= ' --hintsfile='.$xdef_file if(defined $xdef);
	$run .= ' --extrinsicCfgFile='.$cfg_file if(-e $cfg_file);
	$run .= " $file_name";
	$run .= " > $o_file";

	$LOG->add_entry("STARTED", $backup, "") if(defined $LOG);

        if (-f $backup && ! $OPT_F){
                print STDERR "re reading augustus report.\n" unless $main::quiet;
                print STDERR "$backup\n" unless $main::quiet;
		$o_file = $backup;
        }
        else {
                print STDERR "running  augustus.\n" unless $main::quiet;
		write_xdef_file($xdef, $xdef_file) if defined $xdef;
		FastaFile::writeFile(\$fasta, $file_name);

		my $w = new Widget::augustus();
                $w->run($run);
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
            #File::Copy::copy($o_file, $backup) unless(! -f $backup); #temp
	    unlink($o_file); #temp
        }
        catch Error::Simple with {
            my $E = shift;

	    unlink($o_file);

            #retry predictor in different location and parse again
	    my $file_name = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.augustus.fasta";
	    my $xdef_file = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.xdef\.augustus";
	    my $o_file    = "$the_void/$seq_id\.$offset-$end\.$hmm_name\.auto_annotator\.augustus";
	    my $run = $command;
	    $run .= ' --UTR=off';
	    $run .= ' --hintsfile='.$xdef_file if defined $xdef;
	    $run .= ' --extrinsicCfgFile='.$cfg_file if -e $cfg_file;
	    $run .= " $file_name";
	    $run .= " > $o_file";
            write_xdef_file($xdef, $xdef_file) if(defined $xdef);
	    FastaFile::writeFile(\$fasta, $file_name);
            my $w = new Widget::augustus();
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
	   
	   #try one more time
	   if ($? != 0){
	       ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
	       $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	       {
		   local $/ = \1;
		   while (my $line = <$CHLD_ERR>){
		       print STDERR $line unless($main::quiet);
		   }
	       }
	       waitpid $pid, 0;
	   }
	   die "ERROR: Augustus failed\n" if $? != 0;
	}
	else {
	   die "you must give Widget::augustus a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_xdef {
       my $seq       = shift;
       my $p_coors  = shift;
       my $n_coors = shift;
       my $i_coors = shift;
       my $s       = shift;
       my $offset  = shift;
       my $i_flank = shift;
       my $q_id    = shift;

	my $group = 'gb|bogus;';

        my @xdef;

       foreach my $p (@$p_coors){
                my $c_b = $p->[0] - $offset;
                my $c_e = $p->[1] - $offset;

		my $l  = "$q_id\tPROTEIN\tCDSpart";
		   $l .= "\t".$c_b."\t".$c_e."\t"."1e-100"."\t".$s;
                   $l .= "\t".'.'."\t".'group=p_pieces;'."source=P";

                push(@xdef, $l);
        }

        return \@xdef if(!@$i_coors);

        foreach my $i (@$i_coors){
	   my $i_b = ($i->[0] - $offset);# + ($i_flank-1);
	   my $i_e = ($i->[1] - $offset);# - ($i_flank-1);

	   #next if abs($i_b - $i_e) < 2*$i_flank;
	   next if abs($i_b - $i_e) < 25;
	   
	   my $l  = "$q_id\tEST-INTRON\tintronpart";
	   $l .= "\t".$i_b."\t".$i_e."\t"."1e-1000"."\t".$s;
	   $l .= "\t".'.'."\t".'group=n_pieces;'."source=E";
	   
	   push(@xdef, $l);
        }

       foreach my $n (@$n_coors){
                my $e_b = $n->[0] - $offset;
                my $e_e = $n->[1] - $offset;

                my $l  = "$q_id\tEST-EXON\texonpart";
                   $l .= "\t".$e_b."\t".$e_e."\t"."1e-1000"."\t".$s;
                   $l .= "\t".'.'."\t".'group=n_peices;'."source=E";

                push(@xdef, $l);
        }

        return \@xdef;

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
sub parse_gene {
	my $stuff = shift;

	shift(@{$stuff});

	return if(!@{$stuff});

	my %this_gene;
	my $i = 0;
	my $gene_name;
	my $score;
	foreach my $datum (@{$stuff}){
		last if $datum =~ /\# end gene/;

		my @fields = split(/\s+/, $datum);

		if (! @fields){
		}
		elsif ($fields[0] =~ /^\#/){
                }
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2] eq 'gene'){
			$gene_name = $fields[0].'-'.$fields[8];
			$score     = $fields[5];

			$this_gene{name}  = $gene_name;
			$this_gene{score} = $score;
		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'CDS' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

			push(@{$this_gene{CDS}}, {'b'      => $fields[3],
			                          'e'      => $fields[4],
					          'strand' => $strand,
					          'frame'  => $fields[7],
						  'score'  => $score,
						  'type'   => 'CDS',
				                  });

		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'start_codon' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

			$this_gene{start_codon} = {'b'      => $fields[3],
						   'e'      => $fields[4], 
		                                  'strand' => $strand,
			                          'frame'  => $fields[7], 
		                                  'score'  => $score,
						  'type'   => 'start_codon',
                                                  };

		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'stop_codon' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

                        $this_gene{stop_codon} = {'b'      => $fields[3],
			                          'e'      => $fields[4],
                                                  'strand' => $strand,
                                                  'frame'  => $fields[7],
                                                  'score'  => $score,
						  'type'   => 'stop_codon',
                                                  };

		}
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'transcript'){
                }
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'intron'){
                }
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  =~ /^exon$|^internal$|^initial$|^terminal$|UTR$|^single$/){
		    my $strand = $fields[6] eq '+' ? 1 : -1;

		    push(@{$this_gene{exons}}, {'b'      => $fields[3],
					        'e'      => $fields[4],
					        'strand' => $strand,
					        'frame'  => $fields[7],
					        'score'  => $score,
					        'type'   => 'exon',
					       });
                }
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'tts'){
                }
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'tss'){
                }
                elsif ($fields[1] eq 'protein' && $fields[2]  eq 'sequence'){
                }
		else {
			PostData(\@fields);	
			die "dead in Widget::augustus::parse_gene unknown feature type:".$fields[2]."\n";
		}
		
	}

	return \%this_gene; 
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

        my ($q_name)  = $def =~ /^>(.+)/;

	local $/ = '### gene';
	my $is_new = `grep -c \"\# start gene\" $report`;
	$is_new =~ s/\n+//g;
	local $/ = '# start gene' if($is_new);

	my $fh = new FileHandle();
	   $fh->open($report);

	my $header = <$fh>; # shift off the header

        #checks if Augustus really finished
        unless($header){
            unlink($report);
            die "ERROR: The file $report appears to be incomplete\n".
		"MAKER will need to delete the file, before trying again\n\n";
        }

	my @genes;
	while (my $line = <$fh>){
		chomp($line);

		my @stuff = split(/\n/, $line);

		my $t_count = grep {/^[^\#]/ && /^[^\s]+\s+AUGUSTUS\s+transcript\s+/} @stuff;

		next unless(@stuff > 1); #empty report

		#augustus can call multiple transcripts for inconsistent hints
		if($t_count > 1){
		    my @multi_stuff;
		    my $count;
		    foreach my $l (@stuff){
			if($l =~ /^\#/){ #maybe important someday
			    for (my $i = 0; $i < $t_count; $i++){
				push(@{$multi_stuff[$i]}, $l);
			    }
			}
			elsif($l =~ /^[^\s]+\s+AUGUSTUS\s+gene\s+/){
			    for (my $i = 0; $i < $t_count; $i++){
				push(@{$multi_stuff[$i]}, $l);
			    }
			}
			elsif($l =~ /^[^\s]+\s+AUGUSTUS\s+transcript\s+/){
			    $count++;
			    push(@{$multi_stuff[$count-1]}, $l);
                        }
			elsif($l =~ /^[^\s]+\s+AUGUSTUS\s+[^\s]+\s+/){
			    push(@{$multi_stuff[$count-1]}, $l);
			}
			else{
			    for (my $i = 0; $i < $t_count; $i++){
                                push(@{$multi_stuff[$i]}, $l);
                            }
			}
		    }

		    foreach my $s (@multi_stuff){
			my $gene = parse_gene($s);
			push(@genes, $gene);
		    }
		    next;
		}

		my $gene = parse_gene(\@stuff);
		push(@genes, $gene);
	}

	local $/= "\n";

	$fh->close();

	refactor(\@genes);
	my $phat_hits = load_phat_hits($q_name, $q_seq, \@genes, $params);
	return $phat_hits;
}
#-------------------------------------------------------------------------------
sub refactor {
	my $genes = shift;

	foreach my $g (@{$genes}){
	    my $CDS = $g->{CDS};
	    my $start_codon = $g->{start_codon} || undef;
	    my $stop_codon  = $g->{stop_codon}  || undef;

	    next if(! $stop_codon);

	    my ($last) = @{[sort sort_cdss @{$CDS}]}[-1];

	    #stop codon may be it's own exon for weird predictions
	    my $class= compare::compare($stop_codon->{b}, $stop_codon->{e}, $last->{b}, $last->{e});
            if($class eq '0' && abs($last->{e} - $stop_codon->{b}) != 1 && abs($last->{b} - $stop_codon->{e}) != 1){
		my %CDS = ( 'b' => $stop_codon->{b},
			    'e' => $stop_codon->{e},
			    'frame' => $stop_codon->{frame},
			    'score' => $stop_codon->{score},
			    'strand' => $stop_codon->{strand},
			    'type' => 'CDS',
			    );
		push (@{$g->{CDS}}, \%CDS);
		$last = \%CDS;
	    }
	    
	    if ($last->{strand} == 1){
		$last->{e} = $stop_codon->{e} if($stop_codon);
	    }
	    else {
		$last->{b} = $stop_codon->{b} if ($stop_codon);
	    }
	}
}
#-------------------------------------------------------------------------------
sub sort_cdss {
	if ($a->{strand} == 1){
		$a->{b} <=> $b->{b};
	}
	else {
		$b->{b} <=> $a->{b}
	}
}
#-------------------------------------------------------------------------------
sub load_phat_hits {
	my $q_name  = shift;
	my $q_seq   = shift;
	my $genes   = shift;
	my $params  = shift;

	my $q_len = length_o($q_seq);

	my @keepers;
	foreach my $g (@{$genes}){
	
		my $total_score = $g->{score};

		my %hsps;
		my $i = 0;
		my $f = new Bio::Search::Hit::PhatHit::augustus('-name'         => $g->{name},
								'-description'  => 'NA',
								'-algorithm'    => 'augustus',
								'-length'       => $q_len,
								'-score'        => $total_score,
								);
		
                $f->queryLength($q_len);

		my $features = $g->{CDS};
		#$features = $g->{exons} if (defined $g->{exons}); #removed 3/19/2009

		#added 3/19/2009
		#check for single and double base pair overhangs
                @{$features} = sort {$a->{b} <=> $b->{b}} @{$features};
                my $length = 0;
                foreach my $exon (@{$features}){
                    $length += abs($exon->{e} - $exon->{b}) + 1;
                }

                my $overhang = $length % 3;
                if($overhang != 0){
                    if(($features->[0]->{strand} == 1 && ! defined $g->{stop_codon}) ||
		       ($features->[0]->{strand} == -1 && ! defined $g->{start_codon})
		      ){
                        my $last = $features->[-1];
                        my $l_length = abs($last->{e} - $last->{b}) + 1;

                        while($l_length <= $overhang){
                            pop(@{$features});
                            $overhang -= $l_length;
                            $last = $features->[-1];
                            $l_length = abs($last->{e} - $last->{b}) + 1;
                        }

                        $last->{e} -= $overhang;
                    }
                    elsif(($features->[0]->{strand} == -1 && ! defined $g->{stop_codon}) ||
			  ($features->[0]->{strand} == 1 && ! defined $g->{start_codon})
			 ){
                        my $last = $features->[0];
                        my $l_length = abs($last->{e} - $last->{b}) + 1;

                        while($l_length <= $overhang){
                            shift(@{$features});
                            $overhang -= $l_length;
                            $last = $features->[0];
                            $l_length = abs($last->{e} - $last->{b}) + 1;
                        }

                        $last->{b} += $overhang;
                    }
                    else{
                        die "FATAL: Can not fix overhang in Widget::augustus\n";
                    }
                }

		#build hsps
		my $hit_start = 1;
		foreach my $exon (@{$features}){
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
        		push(@args, 'AUGUSTUS');

        		push(@args, '-bits');
        		push(@args, 'NA');

        		push(@args, '-evalue');
        		push(@args, 'NA');

        		push(@args, '-pvalue');
        		push(@args, $exon->{score});

        		push(@args, '-query_length');
        		push(@args, $q_len);

        		push(@args, '-query_end');
        		push(@args, $exon->{e});

        		push(@args, '-conserved');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_name');
        		push(@args, $g->{name}.":".$i.":".$exon->{type});

        		push(@args, '-hit_end');
        		push(@args, $hit_end);

        		push(@args, '-query_gaps');
        		push(@args, 0);

        		push(@args, '-hit_gaps');
        		push(@args, 0);

        		my $hsp = new Bio::Search::HSP::PhatHSP::augustus(@args);
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

        print STDERR "Widget::augustus::AutoLoader called for: ",
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
