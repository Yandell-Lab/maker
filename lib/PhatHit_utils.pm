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
@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- CLASS FUNCTIONS ----------------------------
#------------------------------------------------------------------------
sub sort_hits {
	my $hit  = shift;
	my $what = shift;

        my $sorted;
        if    ($hit->strand($what) == 1 || $hit->strand($what) == 0) {
                $sorted = $hit->sortFeatures($what);
        }
        elsif ($hit->strand($what) == -1 ){
                $sorted = $hit->revSortFeatures($what);
        }     
        else {
		$hit->show();
		PostData($hit);
		print "caller:".caller()."STRAND:".$hit->strand('query')."\n";

                print "not yet supported in PhatHit_utils::sort_hits\n";
                die;
        }

	return $sorted;
}
#------------------------------------------------------------------------
sub seperate_by_strand {
	my $what = shift;
	my $hits = shift;

	my @p;
	my @m;
	my @x;
	my @z;
	foreach my $hit (@{$hits}){
		
        	if ($hit->strand($what) eq '-1/1'){
                	push(@x, $hit);
        	}
		elsif ($hit->strand($what) eq '1/-1'){
                	push(@x, $hit);
        	}
		elsif ($hit->strand($what) == 1) {
                	push(@p, $hit); 
        	}        
        	elsif ($hit->strand($what) == -1 ){
                	push(@m, $hit);
        	}      
		elsif ($hit->strand($what) == 0 ){
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
	 $new_hit->{is_split} = 1 if $num > 1;
	 
	 foreach my  $hsp (@{$splits{$key}}){
	    $new_hit->add_hsp($hsp);
	 }
	 
	 push(@new_hits, $new_hit);
	 
	 $num++;
      }
   }
   
   return \@new_hits;
}
#------------------------------------------------------------------------
sub split_hit_by_strand {
   my $hits = shift;

   push (@{$hits}, $hits) if (ref($hits) ne 'ARRAY');

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
	    die "ERROR: There is no strand for this HSP\n";
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
	 $new_hit->{is_split} = 1;
	 
	 foreach my  $hsp (@{$pm_splits{$k}}){
	    $new_hit->add_hsp($hsp);
	 }
	 
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
	my $i = 0;
	foreach my $hsp ($hit->hsps){
		
            my $new_hit = new $ref('-name'         => $hit->name,
                                   '-description'  => $hit->description,
                                   '-algorithm'    => $hit->algorithm,
                                   '-length'       => $hit->length,
                                   );

		$new_hit->queryLength($hit->queryLength);
		$new_hit->{is_shattered} = 1;

		$new_hit->add_hsp($hsp);

		push(@new_hits, $new_hit);
		$i++;
	}

	return \@new_hits;
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
                push(@args, Fasta::revComp($q_s))
		unless ref($hsp) =~ /blastp/ 
		    || ref($hsp) =~ /tblast/;
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
                push(@args, Fasta::revComp($h_s))
                unless ref($hsp) =~ /blastp/
                    || ref($hsp) =~ /tblast/
		    || ref($hsp) =~ /protein/;
        }

        push(@args, '-hsp_length');
        push(@args, $hsp->query->length);

	
	my $spaces = $hsp->homology_string =~ tr/ / /;
	my $midd = length($hsp->homology_string) - $spaces;

        push(@args, '-identical');
        push(@args, $midd);

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
        push(@args, $midd);

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

	die "PhatHit_utils::copy have what arg revq, revh, both, copy!\n"
	unless defined($what);

	my $ref = ref($hit);

        my $new_hit = new $ref('-name'         => $hit->name,
                               '-description'  => $hit->description,
                               '-algorithm'    => $hit->algorithm,
                               '-length'       => $hit->length,
                              );

        $new_hit->queryLength($hit->queryLength);


	my @new_hsps;
	foreach my $hsp ($hit->hsps){
		my $args  = load_args($hsp, $what);

		my $ref = ref($hsp);
		my $new_hsp = new $ref(@{$args});
		

		add_splice_data($new_hsp, $hsp, $what) if ref($hsp) =~ /est2genome$/;

                $new_hsp->queryName($hsp->queryName)
		if defined($hsp->queryName);


		my $n_q_s = $hsp->strand('query');
		my $n_h_s = $hsp->strand('hit');

		if    ($what eq 'both') {
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
               $new_hsp->{_indentical_hack}      = $hsp->frac_identical();

		push(@new_hsps, $new_hsp);
	}

	my @sorted;
	if ($what eq 'both' || $what eq 'rev'){
		@sorted = sort {$b->nB('query') <=> $a->nB('query') } @new_hsps;
	}
	elsif ($what eq 'revh'){
		@sorted = sort {$b->nB('hit') <=> $a->nB('hit') } @new_hsps;
	}

	foreach my $hsp (@sorted){
		$new_hit->add_hsp($hsp);
	}

	return $new_hit;
}
#------------------------------------------------------------------------
sub normalize {
        my $hit  = shift;
	my $what = shift;

        my $sorted = sort_hits($hit, $what);

	my %seen;
	my @keepers;
        for (my $i = 0; $i < @{$sorted} -1;$i++){
                my $hsp_i = $sorted->[$i];
		my $nbeg_i = $hsp_i->nB($what);
		my $nend_i = $hsp_i->nE($what);

		for (my $j = 1; $j < @{$sorted};$j++){
			my $hsp_j = $sorted->[$j];
			my $nbeg_j = $hsp_j->nB($what);
			my $nend_j = $hsp_j->nE($what);	

			if ($nbeg_i >= $nbeg_j && $nbeg_i <= $nend_j){
				if ($hsp_i->score > $hsp_j->score){
					push(@keepers, $hsp_i)
					unless $seen{$nbeg_i}{$nend_i};
				}
				else {
					push(@keepers, $hsp_j)
					unless $seen{$nbeg_j}{$nend_j};
				}
				$seen{$nbeg_i}{$nend_i}++;
				$seen{$nbeg_j}{$nend_j}++;

			}
			elsif ($nbeg_j >= $nbeg_i && $nbeg_j <= $nend_i){
                               if ($hsp_i->score > $hsp_j->score){
                                        push(@keepers, $hsp_i)
					unless $seen{$nbeg_i}{$nend_i};
                                }
                                else {
                                        push(@keepers, $hsp_j)
					unless $seen{$nbeg_j}{$nend_j};
                                }

				$seen{$nbeg_i}{$nend_i}++;
				$seen{$nbeg_j}{$nend_j}++;
			}
			else {
			}
			
		}
	}

	my $ref = ref($hit);

        my $new_hit = new $ref('-name'         => $hit->name,
                               '-description'  => $hit->description,
                               '-algorithm'    => $hit->algorithm,
                               '-length'       => $hit->length,
                              );

                $new_hit->queryLength($hit->queryLength);
                $new_hit->{is_normalized} = 1;

		my %crap;
		foreach my $hsp (@keepers){
			$crap{$hsp->nB('hit')}{$hsp->nE('hit')}++;
                	$new_hit->add_hsp($hsp);
		}

	my $s_size = @{$sorted};
	my $k_size = @keepers;

	#print "s_size:$s_size k_size:$k_size\n";
	#sleep 2;

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
	
        return unless @{$lil_fish};

        unless (@{$big_fish}){
                @{$big_fish} = @{$lil_fish};
                return;
        }

	print STDERR "merging blast reports...\n";

	#-- working
	my %big_names;
	my %little_names;
	my @merged;
        foreach my $b_hit (@{$big_fish}){
		$big_names{$b_hit->hsp(0)->hit->seq_id} = $b_hit;
	}
	foreach my $l_hit (@{$lil_fish}){
		$little_names{$l_hit->hsp(0)->hit->seq_id}++;
	}

	my %was_merged;
        foreach my $l_hit (@{$lil_fish}){

		next unless defined($big_names{$l_hit->hsp(0)->hit->seq_id});

		my $b_hit = $big_names{$l_hit->hsp(0)->hit->seq_id};

                my $b_start = $b_hit->nB('query');
                my $b_end   = $b_hit->nE('query');

                ($b_start, $b_end) = ($b_end, $b_start)
                if $b_start > $b_end;
	
                my @new_l_hsps;

                my $l_start = $l_hit->nB('query');
                my $l_end   = $l_hit->nE('query');

               ($l_start, $l_end) = ($l_end, $l_start)
                if $l_start > $l_end;

                #print STDERR "b_start:$b_start l_start:$l_start\n";
                #print STDERR "b_end:$b_end l_end:$l_end\n";

                my $distance =
                $b_start <= $l_start ? $l_start - $b_end : $b_start - $l_end;

                #print STDERR "distance:$distance\n";

                next unless $distance < $max_sep;

                next unless $b_hit->strand('query') eq $l_hit->strand('query');

                next unless $b_hit->strand('hit') eq $l_hit->strand('hit');

                print STDERR "adding new hsp to ".$b_hit->name." ".$b_hit->description."\n";

		$was_merged{$l_hit->hsp(0)->hit->seq_id}++;

                my @new_b_hsps;
                foreach my $b_hsp ($b_hit->hsps) {
                        push(@new_b_hsps, $b_hsp);
                }


               foreach my $l_hsp ($l_hit->hsps){
               		push(@new_b_hsps, $l_hsp);
               }

		$b_hit->hsps(\@new_b_hsps);
		$b_hit->{'_sequenceschanged'} = 1;
		push(@merged, $b_hit)	
	}

        foreach my $b_hit (@{$big_fish}){
		next if $was_merged{$b_hit->hsp(0)->hit->seq_id};
		push(@merged, $b_hit);
        }
        foreach my $l_hit (@{$lil_fish}){
              	next if $was_merged{$l_hit->hsp(0)->hit->seq_id};
		push(@merged, $l_hit);
        }

	print STDERR "...finished\n";
        @{$big_fish} = @merged;
}
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "PhatHit::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if ($arg){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------
1;


