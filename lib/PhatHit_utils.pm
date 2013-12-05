#------------------------------------------------------------------------
#----                          PhatHit_utils                         ---- 
#------------------------------------------------------------------------
package PhatHit_utils;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use Fasta;
use SimpleCluster;
use FastaSeq;
use Carp;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- CLASS FUNCTIONS ----------------------------
#------------------------------------------------------------------------
sub sort_hits {
	my $hit  = shift;
	my $what = shift;

        my $sorted = [];
        if    ($hit->strand($what) == 1 || $hit->strand($what) == 0) {
                $sorted = $hit->sortFeatures($what);
        }
        elsif ($hit->strand($what) == -1 ){
                $sorted = $hit->revSortFeatures($what);
        }     
        else {
	        $hit->show();
		PostData($hit);
		print STDERR "caller:".caller()."STRAND:".$hit->strand('query')."\n";

                confess "ERROR: not yet supported in PhatHit_utils::sort_hits\n";
        }

	return $sorted;
}
#------------------------------------------------------------------------
sub separate_by_strand {
	my $what = shift;
	my $hits = shift;
	my $exonerate_flag = shift || 0;
	
	#The exonerate flag specifies whether to take into account
	#exonerate realignment alignment flip.  In other words should
	#a blastn hit be considered to belong to the opposite strand
	#if its exonerate counterpart was flipped

	my @p;
	my @m;
	my @x;
	my @z;
	foreach my $hit (@{$hits}){
	    my $obj = (ref($hit) eq 'ARRAY') ? $hit->[0] : $hit;
	    my $strand = $obj->strand($what);
	    
	    unless($strand =~ /^\-?\d$/){
		confess "FATAL: Could not get stand correctly. Perhaps".
		        "your using the wrong version of BioPerl.\n\n";
	    }	    

	    $strand *= -1 if($exonerate_flag && $obj->{_exonerate_flipped});
	    
	    if ($strand == 1) {
		push(@p, $hit);
	    }        
	    elsif ($strand == -1 ){
		push(@m, $hit);
	    }      
	    elsif ($strand == 0 ){
		push(@z, $hit);
	    }

	}
	return (\@p, \@m, \@x, \@z);
}
#------------------------------------------------------------------------
sub sort_hsps_by_score {
	my $hit = shift;
	
	my @hsps;
	foreach my $hsp ($hit->hsps){
		push(@hsps, $hsp);
	}
	my @sorted = sort {$b->score() <=> $a->score()} @hsps;
	return \@sorted;

}
#------------------------------------------------------------------------
#this method fills in the space between splice site crossing reads to see if
#there is an ORF that infers the location of an exon
sub splice_infer_exon_coors {
    my $ests = shift;
    my $seq = shift;

    #map likely CDS space based on filling in space between splice site crossing reads
    my @sorted;
    my $tM = new CGL::TranslationMachine();
    foreach my $e (@$ests){
	next if($e->num_hsps() <= 1);
	my $fC = $e->start('query');
	my $lC = $e->end('query');
	my $first;
	my $last;
	
	foreach my $hsp ($e->hsps){
	    $first = $hsp if($hsp->start('query') == $fC);
	    $last = $hsp if($hsp->end('query') == $lC);
	    last if($first && $last);
	}
	push(@sorted, [$first, $last]);
    }
    @sorted = sort {$a->[0]->start <=> $b->[0]->start} @sorted; #sort on start coordinate of first HSP
    
    my @coors; #coordinates to return
    my %done; #keep tabs on coordinates already checked
    for(my $i = 0; $i < @sorted-1; $i++){
	my $iB = $sorted[$i]->[0]->start; #entire EST start
	my $iE = $sorted[$i]->[1]->end; #entire EST end
	my $B = $sorted[$i]->[1]->start; #bridging exon space begin

	next if($done{B}); #skip if this starting coor already checked for ORF
	$done{$B} = {};

	my $bad = 0; #anything longer than this already been checked and doesn't work
	my $ok; #anything shorter than this has already been checked and is ok
	for(my $j = $i+1; $j < @sorted; $j++){
	    my $jB = $sorted[$j]->[0]->start; #entire EST start
	    my $jE = $sorted[$j]->[1]->end; #entire EST end
	    my $E = $sorted[$j]->[0]->end; #bridging exon space end
	
	    next if($ok && $E <= $ok); #already verified up to this length
	    next if($bad && $E >= $bad); #over this length is already known to be bad
	    next if($E - $B < 300); #min orf of 300 required (same as most prokaryotic gene finders)
	    next if(compare::compare($iB, $iE, $jB, $jE) ne '0'); #ESTs overlap so there is no in between space
	    next if($done{$B}{$E}); #skip if these coors already checked
	    confess "ERROR: Logic error in Widget:snap::get_xdef\n" if($B > $E);
	    $done{$B}{$E}++;
	    
	    my $L = abs($E - $B) + 1;
	    my $piece = ($sorted[$j]->[0]->strand('query') == 1) ?
		substr_o($seq, $B-1, $L) : Fasta::revComp(substr_o($seq, $B-1, $L));

	    my ($p_seq, $poffset) = $tM->longest_translation($piece);
	    if($poffset < 3 && $L - (3 * length($p_seq) + $poffset) < 3 && $p_seq !~ /X{10}/){
		push(@coors, [$B, $E]);
		$ok = $E if(! $ok || $E > $ok);
	    }
	    else{
		$bad = $E if(!$bad || $E < $bad);
	    }
	}
    }

    return \@coors;
}
#------------------------------------------------------------------------
sub to_begin_and_end_coors {
	my $hits = shift;
	my $what = shift;
	my @coors;
	foreach my $hit (@{$hits}){
		push(@coors, [get_span_of_hit($hit, $what)]);
	}
	return \@coors;
}
#------------------------------------------------------------------------
sub get_hsp_coors_from_hit {
        my $hit = shift;
        my $what = shift;
        my @coors;
        foreach my $hsp ($hit->hsps){
        	push(@coors, [$hsp->nB($what), $hsp->nE($what)]);
        }
        return \@coors;
}
#------------------------------------------------------------------------
sub get_total_score_of_hit {
        my $hit = shift;
	my $total_score = 0;
        foreach my $hsp ($hit->hsps){
		$total_score += $hsp->score();
        }
	return $total_score;
}
#------------------------------------------------------------------------
sub get_hsp_coors {
        my $hits = shift;
        my $what = shift;
        my @coors;
        foreach my $hit (@{$hits}){
		foreach my $hsp ($hit->hsps){
                	push(@coors, [$hsp->nB($what), $hsp->nE($what)]);
		}
        }
        return \@coors;
}
#------------------------------------------------------------------------
sub get_hit_coors {
        my $hits = shift;
        my $what = shift;
        my @coors;
        foreach my $hit (@{$hits}){
	    push(@coors, [$hit->nB($what), $hit->nE($what)]);
        }
        return \@coors;
}
#------------------------------------------------------------------------
sub shadow_i_corrected_coors {
        my $hits = shift;
        my $what = shift;
	
	my $B;
	my $E;
	foreach my $hit (@{$hits}){
	    $B = $hit->start('query') if(!$B || $hit->start('query') < $B);
	    $E = $hit->end('query') if(!$E || $hit->end('query') > $E);
	}

	my $offset = $B;
	my $length = abs($E - $B) + 1;

	my @b_seq = map {0} ($B..$E);

        foreach my $hit (@{$hits}){
	    my $last;
	    my $s = $hit->start('query') - $offset;
	    my $e = $hit->end('query') - $offset;

	    #map hit space (defines intron once exons are mapped)
	    foreach my $i ($s..$e){
		$b_seq[$i] = 1 if($b_seq[$i] == 0);
	    }

	    foreach my $hsp ($hit->hsps){
		$s = $hsp->start('query') - $offset;
		$e = $hsp->end('query') - $offset;

		@b_seq[$s..$e] = map {2} ($s..$e); #map exon
	    }
        }

	my @coors;
	my $first;
	for (my $i = 0; $i < @b_seq; $i++){
	    if($first && $b_seq[$i] == 1){
		my $last = $i - 1 + $offset;	
		push(@coors, [$first, $last]);
		$first = undef;
	    }
	    elsif(! $first && $b_seq[$i] != 1){
		$first = $i + $offset;
	    }

	    #last bp
	    if($i == @b_seq - 1 && $first && $b_seq[$i] != 1){
		my $last = $i + $offset;
                push(@coors, [$first, $last]);
	    }
	}
	
        return \@coors;
}
#------------------------------------------------------------------------
sub split_hit_on_intron {
   my $hits        = shift;
   my $max_intron = shift;
   
   push (@{$hits}, $hits) if (ref($hits) ne 'ARRAY');
   
   my @new_hits;
   
   foreach my $hit (@{$hits}){
      my $ref = ref($hit);
      
      my $sorted = $hit->sortFeatures('query');
      
      unless (defined($sorted->[1])){
	 push(@new_hits, $hit);
	 next;
      }
      
      my %splits;
      my $k = 0;
      my $distance;
      for (my $i = 0; $i < @{$sorted} - 1; $i++){
	 my $hsp_i = $sorted->[$i];
	 my $hsp_j = $sorted->[$i + 1];

	 my $end_i = $hsp_i->end('query'); #end is normalized to be larger than start
	 my $beg_j = $hsp_j->start('query'); #start is normalized to be smaller than end
	 
	 $distance = $beg_j - $end_i;
	 
	 
	 push(@{$splits{$k}}, $hsp_i);
	 
	 $k++ if $distance > $max_intron;
      }
      
      push(@{$splits{$k}}, $sorted->[-1]);
      
      my $num = (keys %splits);
      
      foreach my $key (sort keys %splits){
	 
	 my $new_hit = 
	 new $ref('-name'         => $hit->name,
		  '-description'  => $hit->description,
		  '-algorithm'    => $hit->algorithm,
		  '-length'       => $hit->length,
		 );
	 
	 $new_hit->queryLength($hit->queryLength);
	 $new_hit->database_name($hit->database_name);
	 $new_hit->{is_split} = 1 if $num > 1;
	 
	 foreach my  $hsp (@{$splits{$key}}){
	    $new_hit->add_hsp($hsp);
	 }

	 $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
	 $new_hit->{_label} = $hit->{_label} if($hit->{_label});

	 push(@new_hits, $new_hit);
	 
	 $num++;
      }
   }
   
   return \@new_hits;
}
#------------------------------------------------------------------------
sub split_hit_by_strand {
   my $hits = shift;
   $hits = [$hits] if(ref($hits) ne 'ARRAY' &&
		      ref($hits) =~ /Bio\:\:Search\:\:PhatHit\:\:/);

   my @new_hits;
   foreach my $hit (@{$hits}){
      my $ref = ref($hit);
      
      my @hsps = $hit->hsps;
      
      my %pm_splits;
      
      foreach my $hsp (@hsps){
	 my $strand = $hsp->strand('query');
	 
	 if($strand == 1 || $strand == 0){
	    push(@{$pm_splits{plus}}, $hsp); #plus strand
	 }
	 elsif($strand == -1){
	    push(@{$pm_splits{minus}}, $hsp); #minus strand
	 }
	 else{
	    confess "ERROR: There is no strand for this HSP\n";
	 }
      }  
      
      my @keys = keys %pm_splits;
      
      if(@keys == 1){
	 push(@new_hits, $hit);
	 next;
      }
      
      foreach my $k (@keys){
	 my $new_hit = new $ref('-name'         => $hit->name,
				'-description'  => $hit->description,
				'-algorithm'    => $hit->algorithm,
				'-length'       => $hit->length,
			       );
	 
	 $new_hit->queryLength($hit->queryLength);
	 $new_hit->database_name($hit->database_name);
	 $new_hit->{is_split} = 1;
	 
	 foreach my  $hsp (@{$pm_splits{$k}}){
	    $new_hit->add_hsp($hsp);
	 }

	 $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
	 $new_hit->{_label} = $hit->{_label} if($hit->{_label});
	 
	 push(@new_hits, $new_hit);
      }
   }

   return \@new_hits;
}
#------------------------------------------------------------------------
sub shatter_hit {
	my $hit = shift;

	my $ref = ref($hit);

	my @new_hits;
	foreach my $hsp ($hit->hsps){
		
            my $new_hit = new $ref('-name'         => $hit->name,
                                   '-description'  => $hit->description,
                                   '-algorithm'    => $hit->algorithm,
                                   '-length'       => $hsp->length,
				   '-score'        => $hsp->score,
				   '-bits'         => $hsp->bits,
				   '-significance' => $hsp->significance
                                   );

	    $new_hit->queryLength($hit->queryLength);
	    $new_hit->database_name($hit->database_name);
	    $new_hit->{is_shattered} = 1;
	    
	    $new_hit->add_hsp($hsp);

	    $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
	    $new_hit->{_label} = $hit->{_label} if($hit->{_label});

	    push(@new_hits, $new_hit);
	}

	return \@new_hits;
}
#------------------------------------------------------------------------
sub shatter_all_hits {
	my $hits = shift;

	my $ref = ref($hits->[0]);

	my @new_hits;
	foreach my $hit (@$hits){
	    foreach my $hsp ($hit->hsps){
		
		my $new_hit = new $ref('-name'         => $hit->name,
				       '-description'  => $hit->description,
				       '-algorithm'    => $hit->algorithm,
				       '-length'       => $hsp->length,
				       '-score'        => $hsp->score,
				       '-bits'         => $hsp->bits,
				       '-significance' => $hsp->significance
				       );
		
		$new_hit->queryLength($hit->queryLength);
		$new_hit->database_name($hit->database_name);
		$new_hit->{is_shattered} = 1;
		
		$new_hit->add_hsp($hsp);

		$new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
		$new_hit->{_label} = $hit->{_label} if($hit->{_label});
		
		push(@new_hits, $new_hit);
	    }
	}

	return \@new_hits;
}
#------------------------------------------------------------------------
#this method assumes that the features are all on the same strand
#they must also be shattered hits
sub make_flat_hits {
	my $hits = shift;
	my $seq = shift;

	return [] if(! @$hits);

	my $ref = ref($hits->[0]);
	my $hsp_ref = ref($hits->[0]->{_hsps}->[0]);
	my $strand = $hits->[0]->strand;

	my $coors  = get_hit_coors($hits, 'query');
	my $pieces = Shadower::getVectorPieces($coors, 0);

	my @new_hits;

	foreach my $p (@$pieces){
	    my $nB = $p->{b};
	    my $nE = $p->{e};
	    my $length = abs($nE -$nB) + 1;
	    
	    #natural end and beginning are non-normalized values
	    ($nB, $nE) = ($nE, $nB) if($nB < $nE && $strand == -1);

	    my $pSeq = ($seq) ? substr_o($seq, $p->{b}-1, $p->{e}-$p->{b}+1) : 'N'x($p->{e}-$p->{b}+1);
	    $pSeq = Fasta::revComp(\$pSeq) if($strand == -1);

	    my $new_hit = new $ref('-name'         => "flat_hit_$nB\_$nE",
				   '-description'  => '',
				   '-algorithm'    => $hits->[0]->algorithm,
				   '-length'       => $length,
				   '-score'        => $length,
				   '-bits'         => $length*2,
				   '-significance' => 0
				   );
		
	    $new_hit->queryLength($length);
	    
	    my @args;
	    
	    push(@args, '-query_start');
	    push(@args, $nB);

	    push(@args, '-query_seq');
	    push(@args, $pSeq);
	    
	    push(@args, '-score');
	    push(@args, $length);
	    push(@args, '-homology_seq');
	    push(@args, $pSeq);
	    
	    push(@args, '-hit_start');
	    push(@args, 1);
	    
	    push(@args, '-hit_seq');
	    push(@args, $pSeq);
	    
	    push(@args, '-hsp_length');
	    push(@args, $length);
	    
	    push(@args, '-identical');
	    push(@args, $length);
	    
	    push(@args, '-hit_length');
	    push(@args, $length);
		 
	    push(@args, '-query_name');
	    push(@args, $hits->[0]->{_hsps}->[0]->query_name);

	    push(@args, '-algorithm');
	    push(@args, $hits->[0]->algorithm);
	    
	    push(@args, '-bits');
	    push(@args, $length*2);
	    
	    push(@args, '-evalue');
	    push(@args, 0);

	    push(@args, '-pvalue');
	    push(@args, 0);
	    
	    push(@args, '-query_length');
	    push(@args, $length);
	    
	    push(@args, '-query_end');
	    push(@args, $nE);
	    
	    push(@args, '-conserved');
	    push(@args, $length);
	    
	    push(@args, '-hit_name');
	    push(@args, "flat_hit_$nB\_$nE");

	    push(@args, '-hit_end');
	    push(@args, $length);

	    push(@args, '-query_gaps');
	    push(@args, 0);
	    
	    push(@args, '-hit_gaps');
	    push(@args, 0);
	    
	    my $hsp = new $hsp_ref(@args);
	    $hsp->queryName($hits->[0]->{_hsps}->[0]->query_name);
	    #-------------------------------------------------
	    # setting strand because bioperl is all messed up!
	    #------------------------------------------------
	    if ($strand == 1 ){
		$hsp->{_strand_hack}->{query} = 1;
		$hsp->{_strand_hack}->{hit}   = 1;
	    }
	    else {
		$hsp->{_strand_hack}->{query} = -1;
		$hsp->{_strand_hack}->{hit}   =  1;
	    }

	    $new_hit->add_hsp($hsp);
	    
	    push(@new_hits, $new_hit);
	}

	return \@new_hits;
}
#------------------------------------------------------------------------
sub clip_5_utr {return _clip(shift,shift,1,0);}
#------------------------------------------------------------------------
sub clip_3_utr {return _clip(shift,shift,0,1);}
#------------------------------------------------------------------------
sub clip_utr {return _clip(shift,shift,1,1);}
#------------------------------------------------------------------------
sub _clip {
    my $hit = shift;
    my $seq = shift;
    my $trim5 = shift;
    my $trim3 = shift;

    my $offset = $hit->{translation_offset};
    my $end = $hit->{translation_end};
    my $strand = $hit->strand('query');

    if(!$end || !defined($offset)){
	confess "ERROR: Need seq to determine translation in PhatHit_utils::_clip\n" if(!$seq);
	my $transcript_seq  = maker::auto_annotator::get_transcript_seq($hit, $seq);
	(undef, $offset, $end, undef, undef) = maker::auto_annotator::get_translation_seq($transcript_seq, $hit);
    }

    return undef if(!$end); #no CDS

    my @hsps = @{sort_hits($hit, 'query')};

    my $coorB = $hit->{_TSTART}{query};
    my $coorE = $hit->{_TEND}{query};

    my $tB = ($trim5) ? $coorB : $hit->nB('query');
    my $tE = ($trim3) ? $coorE : $hit->nE('query');

    my $hit_start = 1;
    my $hit_length = 0;
    my @new_hsps;
    ($tB, $tE) = ($tE, $tB) if($tB > $tE); #sort by value
    foreach my $hsp ($hit->hsps){
	my $B = $hsp->start();
	my $E = $hsp->end();

	#see if hsp even overlaps the current CDS
	my $class = compare::compare($B, $E, $tB, $tE);

	next if($class eq '0');

	#figure out change to trim off seq string
        my $change = 0; #only change from start of translation needed
        if($B < $tB){
            $change = abs($tB - $B) if($strand == 1);
            $B = $tB;
        }
        if($E > $tE){
            $change = abs($tE - $E) if($strand == -1);
            $E = $tE;
        }

	my $length = abs($E-$B)+1; #new length

	#set natural begining (important for correct strandedness)
	my ($nB, $nE) = ($strand == 1) ? ($B, $E) : ($E, $B) ;

	my $qrSeq = substr($hsp->query_string(),
			  $change,
			  $length);

	my $htSeq = substr($hsp->hit_string(),
			  $change,
			  $length);

	my $hoSeq = substr($hsp->homology_string(),
			  $change,
			  $length);

	my @args;
	
	push(@args, '-query_start');
	push(@args, $nB);
	
	push(@args, '-query_seq');
	push(@args, $qrSeq);
	
	push(@args, '-score');
	push(@args, $hsp->score);
	push(@args, '-homology_seq');
	push(@args, $hoSeq);
	
	push(@args, '-hit_start');
	push(@args, $hit_start);
	
	push(@args, '-hit_seq');
	push(@args, $htSeq);
	
	push(@args, '-hsp_length');
	push(@args, $length);
	
	push(@args, '-identical');
	push(@args, $length);
	    
	push(@args, '-hit_length');
	push(@args, $length);
	
	push(@args, '-query_name');
	push(@args, $hsp->query_name);
	
	push(@args, '-algorithm');
	push(@args, $hsp->algorithm);
	
	push(@args, '-bits');
	push(@args, $hsp->bits);
	
	push(@args, '-evalue');
	push(@args, $hsp->evalue);
	
	push(@args, '-pvalue');
	push(@args, $hsp->pvalue);
	
	push(@args, '-query_length');
	push(@args, $length);
	
	push(@args, '-query_end');
	push(@args, $nE);
	
	push(@args, '-conserved');
	push(@args, $length);
	
	push(@args, '-hit_name');
	push(@args, $hsp->name);
	
	push(@args, '-hit_end');
	push(@args, $hit_start + $length - 1);
	
	push(@args, '-query_gaps');
	push(@args, 0);
	
	push(@args, '-hit_gaps');
	push(@args, 0);

	my $hsp_ref = ref($hsp);	
	my $new_hsp = new $hsp_ref(@args);
	$new_hsp->queryName($hsp->query_name);
	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#------------------------------------------------
	if ($strand == 1 ){
	    $new_hsp->{_strand_hack}->{query} = 1;
	    $new_hsp->{_strand_hack}->{hit}   = 1;
	}
	else {
	    $new_hsp->{_strand_hack}->{query} = -1;
	    $new_hsp->{_strand_hack}->{hit}   =  1;
	}
	$new_hsp->{_label} = $hsp->{_label} if($hsp->{_label});
	
	$hit_start = $hit_start + $length;
	$hit_length += $length;

	push(@new_hsps, $new_hsp);
    }

    my $ref = ref($hit);
    my $new_hit = new $ref('-name'         => $hit->name,
			   '-description'  => $hit->description,
			   '-algorithm'    => $hit->algorithm,
			   '-length'       => $hit_length,
			   '-score'        => $hit_length,
			   '-bits'         => $hit->bits,
			   '-significance' => $hit->significance
			   );    

    foreach my $hsp (@new_hsps){
	$new_hit->add_hsp($hsp);
    }

    confess "ERROR: Logic error there are no HSPs left in PhatHit_utils::_clip\n"
	if($hit_start == 1 || $new_hit->num_hsps == 0); #should change with each HSP

    #adjust translation offset for clip
    if($trim5){
	$end -= $offset;
	$offset -= $offset;
    }

    #add hidden attributes that may be part of a gene prediction
    $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
    $new_hit->{_label} = $hit->{_label} if($hit->{_label});
    $new_hit->{_REMOVE} = $hit->{_REMOVE} if($hit->{_REMOVE});
    $new_hit->{_tran_name} = $hit->{_tran_name} if($hit->{_tran_name});
    $new_hit->{_tran_id} = $hit->{_tran_id} if($hit->{_tran_id});
    $new_hit->{gene_name} = $hit->{gene_name} if($hit->{gene_name});
    $new_hit->{gene_id} = $hit->{gene_id} if($hit->{gene_id});
    $new_hit->{-attrib} = $hit->{-attrib} if($hit->{-attrib});
    $new_hit->{gene_attrib} = $hit->{gene_attrib} if($hit->{gene_attrib});
    $new_hit->{_Alias} = $hit->{_Alias} if($hit->{_Alias});
    $new_hit->{_AED} = $hit->{_AED} if($hit->{_AED});
    $new_hit->{_eAED} = $hit->{_eAED} if($hit->{_eAED});
    $new_hit->{translation_offset} = $offset;
    $new_hit->{translation_end} = $end;

    return $new_hit;
}
#------------------------------------------------------------------------
sub adjust_start {return _adjust(shift,shift,1,0)}
#------------------------------------------------------------------------
sub adjust_stop {return _adjust(shift,shift,0,1)}
#------------------------------------------------------------------------
sub adjust_start_stop {return _adjust(shift,shift,1,1)}
#------------------------------------------------------------------------
sub _adjust {
    my $hit = shift;
    my $seq = shift;
    my $fixstart = shift;
    my $fixstop = shift;
	
    my $offset = $hit->{translation_offset};
    my $end = $hit->{translation_end};
    my $strand = $hit->strand('query');
    my $B = $hit->start;
    my $E = $hit->end;

    confess "ERROR: Need seq to determine translation adjustments in PhatHit_utils::_adjust\n" if(!$seq);

    my $transcript_seq  = maker::auto_annotator::get_transcript_seq($hit, $seq);
    my $slength = length_o($seq);
    my $tlength = length($transcript_seq); #edge

    if(!$end || !defined($offset)){
        (undef, $offset, $end, undef, undef) = maker::auto_annotator::get_translation_seq($transcript_seq, $hit);
    }

    my $tM = new CGL::TranslationMachine();

    #fix stop codon by walking downstream
    my $has_stop = $tM->is_ter_codon(substr($transcript_seq, $end-1-3, 3));
    my $repeat = 0; #count X codons
    if($fixstop && !$has_stop){
	my $z = $end;
	my $zseq = $transcript_seq;
	my $zlength = $tlength;
	my $edge = $hit->nE('query'); #edge of transcript in contig space
	while(!$has_stop){
	    $z += 3;
	    if($z > $zlength+1){
		my $last = $edge;
		$edge = ($strand == 1) ? $edge + 100 : $edge - 100;
		$edge = $slength if($edge > $slength); #cannot walk past end of contig
		$edge = 1 if($edge < 1); #can't walk apst start of contig 
		my $l = abs($last - $edge);
		my $add_seq = ($strand == 1) ?
			substr_o($seq, $last, $l) : Fasta::revComp(substr_o($seq, $last-$l-1, $l));

		$zseq .= $add_seq;
		$zlength += $l;
	    }
            last if($z > $zlength + 1); #can't go off edge of translatable seq (contig)
            my $codon = substr($zseq, $z-1-3, 3);
	    
	    #don't try and extend through string of N's beyond edge of existing transcript
	    if($tM->translate($codon) eq 'X'){
		$repeat++;
		last if($repeat > 5 && $z-$tlength+1 > 0);
	    }
	    else{
		$repeat = 0;
	    }
	    
	    if($tM->is_ter_codon($codon)){
                $has_stop = 1;
                $end = $z;
            }
	}
	
	if($has_stop){
	    my $diff = $end-($tlength+1);
	    if($strand == 1){
		$E += $diff if($diff > 0);
	    }
	    else{
		$B -= $diff if($diff > 0);
	    }    
	}
    }

    #fix start codon by first walking upsream and then walking into CDS
    my $has_start = $tM->is_start_codon(substr($transcript_seq, $offset, 3));
    $repeat = 0;
    if($fixstart && !$has_start){
        #step upstream to find longer ORF
	my $i = $offset;
	my $iseq = $transcript_seq;
	my $edge = $hit->nB('query');
	while(!$has_start){ # must be within contig seq
	    $i -= 3;
	    if($i < 0){
		my $last = $edge;
		$edge = ($strand == 1) ? $edge - 100 : $edge + 100;
		$edge = $slength if($edge > $slength); #cannot walk past contig end 
		$edge = 1 if($edge < 1);  #cannot walk past contig start
		my $l = abs($last - $edge);
		my $add_seq = ($strand == 1) ?
			substr_o($seq, $last-$l-1 , $l) : Fasta::revComp(substr_o($seq, $last, $l));

		$iseq = $add_seq . $iseq;
		$i += $l;
	    }
	    last if($i < 0);
	    my $codon = substr($iseq, $i, 3);
	    last if($tM->is_ter_codon($codon));
	    if($tM->is_start_codon($codon)){
		$has_start = 1;
		$offset = $i;
	    }

	    #don't try and extend through string of N's beyond edge of existing transcript
	    if($tM->translate($codon) eq 'X'){
		$repeat++;
		my $diff = (length($iseq) - $tlength) - $offset;
		last if($repeat > 5 &&  $diff > 0);
	    }
	    else{
		$repeat = 0;
	    }
	}

	#step downstream if no good upstream start
	my $j = $offset;
	while(!$has_start){
	    $j += 3;
	    last if($j > $tlength-3); #can't wlk downstream of contig end
	    my $codon = substr($transcript_seq, $j, 3);
	    last if($tM->is_ter_codon($codon));
	    if($tM->is_start_codon($codon)){
		$has_start = 1;
		$offset = $j;
	    }
	}

	if($has_start){
	    my $diff = (length($iseq) - $tlength) - $offset;
	    if($strand == 1){
		$B -= $diff if($diff > 0);
	    }
	    else{
		$E += $diff if($diff > 0);
	    }	    
	}
    }

    my $ref = ref($hit);
    my $hsp_ref = ref($hit->{_hsps}->[0]);

    my $new_hit = new $ref('-name'         => $hit->name,
			   '-description'  => $hit->description,
			   '-algorithm'    => $hit->algorithm,
			   '-length'       => $hit->length,
			   '-score'        => $hit->score,
			   '-bits'         => $hit->bits,
			   '-significance' => $hit->significance
			   );

    my @hsps = @{sort_hits($hit, 'query')};
    my $hit_start = 1;
    for (my $i=0; $i < @hsps; $i++){
	my $hsp = $hsps[$i];

	my $hB = $hsp->start('query');
	my $hE = $hsp->end('query');

	#fix first and last hsps
	if(($i == 0 && $strand == 1) || ($i == @hsps - 1 && $strand == -1)){
	    $hB = $B;
	}
	if(($i == 0 && $strand == -1) || ($i == @hsps - 1 && $strand == 1)){
	    $hE = $E;
	}

	my $length = abs($hE-$hB)+1; #new length

	#set natural begining (important for correct strandedness)
	my ($nB, $nE) = ($strand == 1) ? ($hB, $hE) : ($hE, $hB) ;

	my $qrSeq = substr_o($seq,
			     $hB-1,
			     $length);
	$qrSeq = Fasta::revComp($qrSeq) if($strand == -1);
	my $htSeq = $qrSeq;
	my $hoSeq = '|'x$length;

	my @args;
	
	push(@args, '-query_start');
	push(@args, $nB);
	
	push(@args, '-query_seq');
	push(@args, $qrSeq);
	
	push(@args, '-score');
	push(@args, $hsp->score);
	push(@args, '-homology_seq');
	push(@args, $hoSeq);
	
	push(@args, '-hit_start');
	push(@args, $hit_start);
	
	push(@args, '-hit_seq');
	push(@args, $htSeq);
	
	push(@args, '-hsp_length');
	push(@args, $length);
	
	push(@args, '-identical');
	push(@args, $hsp->identical);
	    
	push(@args, '-hit_length');
	push(@args, $length);
	
	push(@args, '-query_name');
	push(@args, $hsp->query_name);
	
	push(@args, '-algorithm');
	push(@args, $hsp->algorithm);
	
	push(@args, '-bits');
	push(@args, $hsp->bits);
	
	push(@args, '-evalue');
	push(@args, $hsp->evalue);
	
	push(@args, '-pvalue');
	push(@args, $hsp->pvalue);
	
	push(@args, '-query_length');
	push(@args, $length);
	
	push(@args, '-query_end');
	push(@args, $nE);
	
	push(@args, '-conserved');
	push(@args, $length);
	
	push(@args, '-hit_name');
	push(@args, $hsp->name);
	
	push(@args, '-hit_end');
	push(@args, $hit_start + $length - 1);
	
	push(@args, '-query_gaps');
	push(@args, 0);
	
	push(@args, '-hit_gaps');
	push(@args, 0);
	
	my $new_hsp = new $hsp_ref(@args);
	$new_hsp->queryName($hsp->query_name);
	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#------------------------------------------------
	if ($strand == 1 ){
	    $new_hsp->{_strand_hack}->{query} = 1;
	    $new_hsp->{_strand_hack}->{hit}   = 1;
	}
	else {
	    $new_hsp->{_strand_hack}->{query} = -1;
	    $new_hsp->{_strand_hack}->{hit}   =  1;
	}
	$new_hsp->{_label} = $hsp->{_label} if($hsp->{_label});

	$new_hit->add_hsp($new_hsp);
	
	$hit_start = $hit_start + $length;
    }    

    $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
    $new_hit->{_label} = $hit->{_label} if($hit->{_label});

    return $new_hit;
}
#------------------------------------------------------------------------
sub add_splice_data {
    my $new_hsp = shift;
    my $old_hsp = shift;
    my $what    = shift;
    
    my $d = $old_hsp->donor();
    my $a = $old_hsp->acceptor();
    if ($what eq 'revq' || $what eq 'both'){
	$new_hsp->acceptor(Fasta::revComp($d))
	    if defined($d);
	
	$new_hsp->donor(Fasta::revComp($a))
	    if defined($a);
	
    }
    else {
	$new_hsp->donor($d);
	$new_hsp->acceptor($a);
    }
}
#------------------------------------------------------------------------
sub load_args {
	my $hsp    = shift;
	my $action = shift;

        my @args;

	push(@args, '-query_start');
	push(@args, $hsp->start('query'));

	my $q_s = $hsp->query_string();
        push(@args, '-query_seq');

        if ($action eq 'copy'){
                push(@args, $q_s);
        }
        elsif ($action eq 'revq' || $action eq 'both'){
	   if (ref($hsp) =~ /blastp|tblastn|blastx/ && ref($hsp) !~ /tblastx/){
	      push(@args, $q_s);
	   }
	   else{
	      push(@args, Fasta::revComp($q_s));
	   }
        }

        push(@args, '-score');
        push(@args, $hsp->{SCORE});

        push(@args, '-homology_seq');
        push(@args, $hsp->homology_string);


        push(@args, '-hit_start');
        push(@args, $hsp->start('hit'));

        my $h_s = $hsp->hit_string();

	push(@args, '-hit_seq');

        if ($action eq 'copy'){
                push(@args, $h_s);
        }
        elsif ($action eq 'revh' || $action eq 'both'){
	   if(ref($hsp) =~ /blastp|tblastn|blastx/ && ref($hsp) !~ /tblastx/){
	       push(@args, $h_s);
	   }
	   else{
	      push(@args, Fasta::revComp($h_s))
           } 
        }

        push(@args, '-hsp_length');
        push(@args, $hsp->query->length);
	
	#my $spaces = $hsp->homology_string =~ tr/ / /;
	#my $midd = length($hsp->homology_string) - $spaces;

        push(@args, '-identical');
        push(@args, $hsp->num_identical);

        push(@args, '-hit_length');
        push(@args, $hsp->hit->length);

        push(@args, '-query_name');
        push(@args, $hsp->name());

        push(@args, '-algorithm');
        push(@args, $hsp->algorithm);

        push(@args, '-bits');
        push(@args, $hsp->bits);

        push(@args, '-evalue');
        push(@args, $hsp->evalue);

        push(@args, '-pvalue');
        push(@args, $hsp->signifcance);

        push(@args, '-query_length');
        push(@args, $hsp->query->length);

        push(@args, '-query_end');
        push(@args, $hsp->end('query'));

        push(@args, '-conserved');
        push(@args, $hsp->num_conserved);

        push(@args, '-hit_name');
        push(@args, $hsp->name());

        push(@args, '-hit_end');
        push(@args, $hsp->end('hit'));

        push(@args, '-query_gaps');
        push(@args, $hsp->gaps('query'));

        push(@args, '-hit_gaps');
        push(@args, $hsp->gaps('hit'));	

	return \@args;
}
#------------------------------------------------------------------------
sub copy {

	my $hit   = shift;
	my $what  = shift;

	confess "PhatHit_utils::copy have what arg revq, revh, both, copy!\n"
	unless defined($what);

	my $ref = ref($hit);

        my $new_hit = new $ref('-name'         => $hit->name,
                               '-description'  => $hit->description,
                               '-algorithm'    => $hit->algorithm,
                               '-length'       => $hit->length,
                              );

        $new_hit->queryLength($hit->queryLength);
	$new_hit->database_name($hit->database_name);
	$new_hit->{_label} = $hit->{_label};

	my @new_hsps;
	foreach my $hsp ($hit->hsps){
	    my $args  = load_args($hsp, $what);
	    
	    my $ref = ref($hsp);
	    my $new_hsp = new $ref(@{$args});
	    
	    
	    add_splice_data($new_hsp, $hsp, $what) if ref($hsp) =~ /2genome$/;
	    
	    $new_hsp->queryName($hsp->queryName)
		if defined($hsp->queryName);
	    
	    my $n_q_s = $hsp->strand('query');
	    my $n_h_s = $hsp->strand('hit');
	    
	    if ($what eq 'both') {
		$n_q_s = -1*$n_q_s;
		$n_h_s = -1*$n_h_s;		    
	    }
	    elsif ($what eq 'revq'){
		$n_q_s = -1*$n_q_s;
	    }
	    elsif ($what eq 'revh'){
		$n_h_s = -1*$n_h_s;
	    }
	    
	    $new_hsp->{_strand_hack}->{query} = $n_q_s;
	    $new_hsp->{_strand_hack}->{hit}   = $n_h_s;
	    $new_hsp->{_identical_hack}       = $hsp->frac_identical();
	    $new_hsp->{_label}                = $hsp->{_label};
	    $new_hsp->{_CIGAR}                = $hsp->{_CIGAR} if($hsp->{_CIGAR});
	    
	    push(@new_hsps, $new_hsp);
	}

	my @sorted;
	if ($what eq 'both' || $what eq 'rev'){
		@sorted = sort {$b->nB('query') <=> $a->nB('query') } @new_hsps;
		$new_hit->{_was_flipped} = 1;
	}
	elsif ($what eq 'revh'){
		@sorted = sort {$b->nB('hit') <=> $a->nB('hit') } @new_hsps;
	}

	foreach my $hsp (@sorted){
	    $new_hit->add_hsp($hsp);
	}

	$new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
	$new_hit->{_label} = $hit->{_label} if($hit->{_label});

	return $new_hit;
}
#------------------------------------------------------------------------
#flattens a hit with introns
sub normalize {
    my $hit = shift;
    my $seq = shift;

    my @hsps = @{sort_hits($hit, 'query')};

    my $bad = 0;
    for (my $i = 0; $i < @hsps-1; $i++){
	my $fi = $hsps[$i];
	my $fj = $hsps[$i+1];
	if(compare::compare($fi->start, $fi->end, $fj->start, $fj->end)){
	    $bad = 1;
	    last;
	}
    }
    return $hit if(!$bad);

    my $coors  = PhatHit_utils::get_hsp_coors_from_hit($hit, 'query');
    my $pieces = Shadower::getPieces($seq, $coors, 0);

    my $hit_start = 1;
    my $hit_length = 0;
    my $strand = $hit->strand('query');
    my @new_hsps;
    foreach my $p (@$pieces){
	my $B = $p->{b};
	my $E = $p->{e};
	my $hseq = ($strand == 1) ? $p->{piece} : Fasta::revComp($p->{piece});
	my $length = length($hseq);

	my @args;
	
	push(@args, '-query_start');
	push(@args, $B);
	
	push(@args, '-query_seq');
	push(@args, $hseq);
	
	push(@args, '-score');
	push(@args, $length);

	push(@args, '-homology_seq');
	push(@args, '|'x$length);
	
	push(@args, '-hit_start');
	push(@args, $hit_start);

	push(@args, '-hit_seq');
	push(@args, $hseq);
	
	push(@args, '-hsp_length');
	push(@args, $length);
	
	push(@args, '-identical');
	push(@args, $length);
	    
	push(@args, '-hit_length');
	push(@args, $length);
	
	push(@args, '-query_name');
	push(@args, $hsps[0]->query_name);
	
	push(@args, '-algorithm');
	push(@args, $hsps[0]->algorithm);
	
	push(@args, '-bits');
	push(@args, $length);
	
	push(@args, '-evalue');
	push(@args, 0);
	
	push(@args, '-pvalue');
	push(@args, 0);
	
	push(@args, '-query_length');
	push(@args, $length);
	
	push(@args, '-query_end');
	push(@args, $E);
	
	push(@args, '-conserved');
	push(@args, $length);
	
	push(@args, '-hit_name');
	push(@args, $hsps[0]->name);
	
	push(@args, '-hit_end');
	push(@args, $hit_start + $length - 1);
	
	push(@args, '-query_gaps');
	push(@args, 0);
	
	push(@args, '-hit_gaps');
	push(@args, 0);

	my $hsp_ref = ref($hsps[0]);
	my $new_hsp = new $hsp_ref(@args);
	$new_hsp->queryName($hsps[0]->query_name);
	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#------------------------------------------------
	if ($strand == 1 ){
	    $new_hsp->{_strand_hack}->{query} = 1;
	    $new_hsp->{_strand_hack}->{hit}   = 1;
	}
	else {
	    $new_hsp->{_strand_hack}->{query} = -1;
	    $new_hsp->{_strand_hack}->{hit}   =  1;
	}
	$new_hsp->{_label} = $hsps[0]->{_label} if($hsps[0]->{_label});
	
	$hit_start = $hit_start + $length;
	$hit_length += $length;

	push(@new_hsps, $new_hsp);
    }

    my $ref = ref($hit);
    my $new_hit = new $ref('-name'         => $hit->name,
			   '-description'  => $hit->description,
			   '-algorithm'    => $hit->algorithm,
			   '-length'       => $hit_length,
			   '-score'        => $hit_length,
			   '-bits'         => $hit->bits,
			   '-significance' => $hit->significance
			   );    

    foreach my $hsp (@new_hsps){
	$new_hit->add_hsp($hsp);
    }

    confess "ERROR: Logic error there are no HSPs left in PhatHit_utils::_clip\n"
	if($hit_start == 1 || $new_hit->num_hsps == 0); #should change with each HSP

    #add hidden attributes that may be part of a gene prediction
    $new_hit->{_HMM} = $hit->{_HMM} if($hit->{_HMM});
    $new_hit->{_label} = $hit->{_label} if($hit->{_label});
    $new_hit->{_REMOVE} = $hit->{_REMOVE} if($hit->{_REMOVE});
    $new_hit->{_tran_name} = $hit->{_tran_name} if($hit->{_tran_name});
    $new_hit->{_tran_id} = $hit->{_tran_id} if($hit->{_tran_id});
    $new_hit->{gene_name} = $hit->{gene_name} if($hit->{gene_name});
    $new_hit->{gene_id} = $hit->{gene_id} if($hit->{gene_id});
    $new_hit->{-attrib} = $hit->{-attrib} if($hit->{-attrib});
    $new_hit->{gene_attrib} = $hit->{gene_attrib} if($hit->{gene_attrib});
    $new_hit->{_Alias} = $hit->{_Alias} if($hit->{_Alias});
    $new_hit->{_AED} = $hit->{_AED} if($hit->{_AED});
    $new_hit->{_eAED} = $hit->{_eAED} if($hit->{_eAED});

    return $new_hit;
}
#------------------------------------------------------------------------
sub is_contigous {
	my $hit = shift;
	
	my $q_sorted = sort_hits($hit, 'query');
	my $h_sorted = sort_hits($hit, 'hit');

	#print $hit->name."\n";
	for (my $i = 0; $i < @{$q_sorted};$i++){
		my $q_hsp = $q_sorted->[$i];
		my $h_hsp = $h_sorted->[$i];

		my $q_nB = $q_hsp->nB('query');
		my $q_nE = $q_hsp->nE('query');

		my $h_nB = $h_hsp->nB('query');
                my $h_nE = $h_hsp->nE('query');

		#print "q_nB:$q_nB h_nB:$h_nB q_nE:$q_nE h_nE:$h_nE\n";

		return 0 if $q_nB != $h_nB;
		return 0 if $q_nE != $h_nE;
	}
	return 1;
}
#------------------------------------------------------------------------
sub get_span_of_hit {
	my $hit  = shift;
	my $what = shift || 'query';

	confess "ERROR: Incorrect reference type\n" if(ref($hit) !~ /Bio/);
	
	if ($hit->strand($what) eq '-1/1'){
		print STDERR " mixed strand feature in PhatHit_utils\n";	
	}
        elsif ($hit->strand($what) eq '1/-1'){
                print STDERR " mixed strand feature in PhatHit_utils\n";
        }

	return ($hit->nB($what), $hit->nE($what));
}
#------------------------------------------------------------------------
sub add_offset {
	my $lil_fish = shift;
	my $offset   = shift;

	foreach my $f (@{$lil_fish}){
		foreach my $hsp ($f->hsps){
			my $new_start = $offset + $hsp->start('query');
			my $new_end   = $offset + $hsp->end('query');

			$hsp->query->location->start($new_start);
			$hsp->query->location->end($new_end);
			$hsp->{'_sequenceschanged'} = 1;
		}
		$f->{nB}{query} += $offset if (defined($f->{nB}) && defined($f->{nB}{query}));
		$f->{nE}{query} += $offset if (defined($f->{nE}) && defined($f->{nE}{query})); 
		$f->{'_sequenceschanged'} = 1;
	} 

}
#------------------------------------------------------------------------
sub reset_query_lengths {
	my $features       = shift;
        my $query_length   = shift;

        foreach my $f (@{$features}){
		$f->queryLength($query_length);
                $f->{'_sequenceschanged'} = 1;
        }

}
#------------------------------------------------------------------------
sub merge_hits {
   my $big_fish = shift;
   my $lil_fish = shift;
   my $max_sep  = shift;   

   return if(! @{$lil_fish});

   print STDERR "merging blast reports...\n" unless($main::quiet);

   my %sets;
   foreach my $h (@$big_fish, @$lil_fish){
       push (@{$sets{$h->name}}, $h);
   }

   my @merged;
   while(my $key = each %sets){
       if(@{$sets{$key}} == 1){
	   push(@merged, @{$sets{$key}});
	   next;
       }

       my $clusters = SimpleCluster::cluster_hits($sets{$key}, $max_sep);

       foreach my $c (@$clusters){
	   if(@$c == 1){
	       push(@merged, @$c);
	   }
	   else{
	       my @new_hsps;
	       foreach my $h (@$c){
		   push(@new_hsps, $h->hsps);
	       }
	       
	       $c->[0]->hsps(\@new_hsps);
	       $c->[0]->{'_sequenceschanged'} = 1;
	       $c->[0]->{'_sequences_was_merged'} = 1;
	       $c->[0]->{'nB'} = undef; #force recompute of start
	       $c->[0]->{'nE'} = undef; #force recompute of end
	       $c->[0]->{'strand'} = undef; #force recompute of strand
	       $c->[0]->{'strand'} = undef; #force recompute of strand
	       $c->[0]->{'LAlnQ'} = undef; #reset this
	       $c->[0]->{'LAlnH'} = undef; #reset this
	       push(@merged, $c->[0]);
	   }
       }
   }

   print STDERR "...finished\n" unless($main::quiet);

  @{$big_fish} = @merged;
   @{$lil_fish} = ();
}
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
#sub AUTOLOAD {
#        my ($self, $arg) = @_;
#
#        my $caller = caller();
#        use vars qw($AUTOLOAD);
#        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
#        $call =~/DESTROY/ && return;
#
#        #print STDERR "PhatHit::AutoLoader called for: ",
#        #      "\$self->$call","()\n";
#        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";
#
#        if ($arg){
#                $self->{$call} = $arg;
#        }
#        else {
#                return $self->{$call};
#        }
#}
#------------------------------------------------------------------------
1;


