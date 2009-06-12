#------------------------------------------------------------------------
#----             	AED.pm		             ---- 
#------------------------------------------------------------------------
package evaluator::AED;
use strict;
use FileHandle;
use PhatHit_utils;
use Shadower;
use cluster;
use compare;
use evaluator::pseudo_hit;
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub txnAED {
	my $box = shift;
	my $parameter = shift;

	my $t =		$box->{transcript};
	my $pol_e =	$box->{est};
	my $pol_p =	$box->{exonerate};

	my $pseudo_transcript = evaluator::pseudo_hit::convert_to_pseudo_hit
					($t);

	my $hit_bag = combine($pol_e, $pol_p);

	return (1,1) if scalar @$hit_bag == 0;

	my @ph_bag;
	foreach my $hit (@$hit_bag) {
		my $ph = evaluator::pseudo_hit::convert_to_pseudo_hit($hit);
		push @ph_bag, $ph;
	}

	foreach my $hit (@ph_bag) {
		evaluator::pseudo_hit::add_splice($hit);
	}

	my $extended_hits = extend_hits(\@ph_bag);

	my $non_redundant_hits = remove_redundancy($extended_hits);


	my $min_AED = 1;
	foreach my $ph (@$non_redundant_hits) {
		my $AED = ph_AED($pseudo_transcript, $ph, $parameter);
		$min_AED = $AED if $AED < $min_AED;
	}


	my $closestAED = cal_closest($pseudo_transcript, \@ph_bag, $parameter);

	return ($min_AED, $closestAED);
}

#------------------------------------------------------------------------------
sub gene_AED {
	my $c = shift; # The gene object;
	my $pol_e_hits = shift;
	my $pol_p_hits = shift;
	my $blastx_hits = shift;
	my $ab_inits  = shift;
	my $seq = shift;

	my @transcripts_in_gene;
	foreach my $t (@$c) {
		my $ph_transcript = evaluator::pseudo_hit::convert_to_pseudo_hit
				($t);
		push @transcripts_in_gene, $ph_transcript;
	}
	
	my $gene_hit = evaluator::pseudo_hit::combine_pseudo_hit (
				\@transcripts_in_gene);


	my @bag;
	foreach my $hit (@$pol_e_hits) {
		push @bag, $hit;
	}
	foreach my $hit (@$pol_p_hits) {
                push @bag, $hit;
        }
	foreach my $hit (@$blastx_hits) {
                push @bag, $hit;
        }
=pod
	foreach my $hit (@$ab_inits) {
		push @bag, $hit;
	}
=cut

	my @good_hits;
	foreach my $hit (@bag) {
		my $ph = evaluator::pseudo_hit::convert_to_pseudo_hit($hit);
		push @good_hits, $ph ;
				
	}

	return 1 if (scalar @good_hits == 0);
	my $evi_ph = evaluator::pseudo_hit::combine_pseudo_hit(\@good_hits);

	my $para = {'exon'=>1};

	my $AED = ph_AED($gene_hit, $evi_ph, $para);
	return $AED;
}
	
#------------------------------------------------------------------------------
sub ph_AED {
	my $a = shift;
	my $b = shift;
	my $para= shift;

	return 1 if $a->{strand} != $b->{strand};

        if (!defined $para->{exon}) { $para->{exon} = 0;}
        if (!defined $para->{intron}) { $para->{intron} = 0;}
        if (!defined $para->{donor}) { $para->{donor} = 0;}
        if (!defined $para->{acceptor}) { $para->{acceptor} = 0;}
        if (!defined $para->{start}) { $para->{start} = 0;}
        if (!defined $para->{stop}) { $para->{stop} = 0;}

	my ($score1, $score2) = (0, 0);
	my $overlap_score;

	foreach my $a_hsp (@{$a->{hsps}}) {
		$score1 += ($a_hsp->[1]-$a_hsp->[0] + 1) * $para->{exon};
		
		for (my $i = $a_hsp->[0]; $i <= $a_hsp->[1]; $i++) {
			if (within($i, $b->{hsps})) {
				$overlap_score += $para->{exon};
			}
		}
	}
	foreach my $b_hsp (@{$b->{hsps}}) {
		$score2 += ($b_hsp->[1]-$b_hsp->[0] + 1) * $para->{exon};
	}

	if ($para->{intron} != 0 ) {
		foreach my $a_intron (@{$a->{introns}}) {
			$score1 += ($a_intron->[1] - $a_intron->[0] + 1)
					* $para->{intron};
		
			for (my $i=$a_intron->[0]; $i<=$a_intron->[1];$i++) {
				if (within($i, $b->{introns})) {
					$overlap_score += $para->{intron};
				}
			}
		}
		foreach my $b_intron (@{$b->{introns}}) {
			$score2 += ($b_intron->[1] - $b_intron->[0] + 1)
					* $para->{intron};
		}
	}

	my %special2;
	$special2{$b->{b}} = 'start';
	$special2{$b->{e}} = 'stop';

	$score1 += $para->{start};
	$score1 += $para->{stop};
	$score2 += $para->{start};
	$score2 += $para->{stop};


	$score1 += (scalar @{$a->{introns}}) * $para->{donor};
	$score1 += (scalar @{$a->{introns}}) * $para->{acceptor};

	foreach my $b_intron (@{$b->{introns}}) {
		$special2{$b_intron->[0]} = 'donor';
		$special2{$b_intron->[1]} = 'acceptor';

		$score2 += $para->{donor};
		$score2 += $para->{acceptor};
	}

        foreach my $a_intron (@{$a->{introns}}) {
		my $donor = $a_intron->[0];
		my $acceptor = $a_intron->[1];

		if (defined $special2{$donor} &&
			$special2{$donor} eq 'donor') {
			$overlap_score += $para->{donor};
		}
		if (defined $special2{$acceptor} &&
			$special2{$acceptor} eq 'acceptor') {
			$overlap_score += $para->{acceptor};
		}
	}
		
	if (defined $special2{$a->{b}} &&
		$special2{$a->{b}} eq 'start') {
			$overlap_score += $para->{start};
	}
        if (defined $special2{$a->{e}} &&
                $special2{$a->{e}} eq 'stop') {
                        $overlap_score += $para->{stop};
        }

	my $congruence = 0.5 * ($overlap_score/$score1 + $overlap_score/$score2);

	return 1-$congruence;
}

#-------------------------------------------------------------------------------
sub cal_closest {
	my $transcript = shift;
	my $bag = shift;
	my $para = shift;

        if (!defined $para->{exon}) { $para->{exon} = 0;}
        if (!defined $para->{intron}) { $para->{intron} = 0;}
        if (!defined $para->{donor}) { $para->{donor} = 0;}
        if (!defined $para->{acceptor}) { $para->{acceptor} = 0;}
        if (!defined $para->{start}) { $para->{start} = 0;}
        if (!defined $para->{stop}) { $para->{stop} = 0;}

	my (@hsp_bag, @intron_bag);
	foreach my $ph (@$bag) {
		push @hsp_bag, @{$ph->{hsps}};
		push @intron_bag, @{$ph->{introns}};
	}

	my $overlap = 0;
	my ($score1, $score2) = (0,0);
	foreach my $hsp (@{$transcript->{hsps}}) {
		for (my $i = $hsp->[0]; $i<= $hsp->[1]; $i++) {
			$score1 += $para->{exon};

			if (within($i, \@hsp_bag)) {
				$overlap += $para->{exon};
				$score2 += $para->{exon};
			}
			elsif (within($i, \@intron_bag)) {
				$score2 += $para->{intron};
			}
		}
	}
	foreach my $intron (@{$transcript->{introns}}) {
                for (my $i = $intron->[0]; $i<= $intron->[1]; $i++) {
                        $score1 += $para->{intron};

                        if (within($i, \@intron_bag)) {
                                $overlap += $para->{intron};
                                $score2 += $para->{intron};
                        }
                        elsif (within($i, \@hsp_bag)) {
                                $score2 += $para->{exon};
                        }
                }
        }
	
	$score1 += $para->{start};
	$score1 += $para->{stop};

	foreach my $ph (@$bag) {
		if ($ph->{b} eq $transcript->{b}) {
			$overlap += $para->{start};
			$score2 += $para->{start};

			last;
		}
	}
        foreach my $ph (@$bag) {
                if ($ph->{e} eq $transcript->{e}) {
                        $overlap += $para->{stop};
			$score2 += $para->{stop};
                        last;
                }
        }
		

	foreach my $intron (@{$transcript->{introns}}) {
		my ($donor, $acceptor) = @$intron;
		$score1 += $para->{donor};	
		$score1 += $para->{acceptor};

		my ($find_donor, $find_acceptor) = (0,0);
		foreach my $coor (@intron_bag) {
			if ($coor->[0] == $intron->[0]) {
				$find_donor = 1;
			}
			if ($coor->[1] == $intron->[1]) {
				$find_acceptor = 1;
			}
		}

		if ($find_donor) {
			$overlap += $para->{donor};
			$score2 += $para->{donor};
		}
		if ($find_acceptor) {
			$overlap += $para->{acceptor};
			$score2 += $para->{acceptor};
		}
	}
	my $congruence = 0.5 * ($overlap/$score1 + $overlap/$score2);

	return 1-$congruence;
}
			
				

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub within {
	my $i      = shift;
	my $coors  = shift;

	foreach my $coor (@$coors) {
		return 1 if $i>= $coor->[0] && $i <= $coor->[1];
	}
	
	return 0;
}

#-------------------------------------------------------------------------------
sub overlap {
	my $gene_boundary = shift;
	my $gene_strand   = shift;
	my $ph = shift;

	return 0 unless $ph->{strand} eq $gene_strand;

	my $b = $ph->{b};
	my $e = $ph->{e};

	($b, $e) = ($e, $b) unless $b<=$e;

	return 0 if $b > $gene_boundary->[1] || $e < $gene_boundary->[0];
	return 1;
}
	
#-------------------------------------------------------------------------------
sub combine {
        my @bag;
        while (my $hits = shift(@_)){
                foreach my $hit (@{$hits}){
                        push(@bag, $hit);
                }
        }
        return \@bag;
}

#------------------------------------------------------------------------
sub ensure_unique {
	my @hits = @_;

	my @polished;

	foreach my $hit (@hits) {
		my $unique = 1;

		foreach my $one (@polished) {
			if (evaluator::pseudo_hit::are_same_pseudo_hits
				($one, $hit)) {
				$unique = 0;
				last; 
			}
		}
		push @polished, $hit if $unique == 1;
	}

	return @polished;
}
#------------------------------------------------------------------------
sub already_exist {
	my $hit = shift;
	my $bag = shift;

	foreach my $one (@$bag) {
		return 1 if evaluator::pseudo_hit::is_redundant_alt_form
				($one, $hit);
	}

	return 0;
}
#------------------------------------------------------------------------
sub remove_redundancy {
	my $hits = shift;

	my @non_redundant;

	foreach my $hit (@$hits) {
		my $status = 'non_redund';
		
		foreach my $good_hit (@non_redundant) {
			if (evaluator::pseudo_hit::is_redundant_alt_form($hit,
				$good_hit)) {
				$status = 'redund';
				last;
			}
		}

		push @non_redundant, $hit if $status eq 'non_redund';
	}

	return \@non_redundant;
}
#------------------------------------------------------------------------
sub extend_hits {
	my $bag = shift;

	my ($b, $e) = ($bag->[0]->{hsps}->[0]->[0], $bag->[0]->{hsps}->[0]->[1]);
	foreach my $hit (@$bag) {
		evaluator::pseudo_hit::add_splice($hit);

		if ($hit->{b} < $b) { $b = $hit->{b};}
		if ($hit->{e} > $e) { $e = $hit->{e};}
	}
	
	my $size = $e-$b+1;
	my @array;

	for (my $i = 0; $i<= $size-1; $i++) {
		$array[$i] = 0;
	}
	# 0 stands for unknonw region; 1 for coding; 2 for intron;
	# -1 for contradicting region;	

	foreach my $hit (@$bag) {
		foreach my $hsp (@{$hit->{hsps}}) {
			my $left = $hsp->[0] - $b;
			my $right = $hsp->[1] - $b;

			for (my $i = $left; $i <= $right; $i++) {
				if ($array[$i] == 0) { $array[$i] = 1;}
				elsif ($array[$i] == 2) {
					$array[$i] = -1;
				}
			}
		}

		foreach my $intron (@{$hit->{introns}}) {
			my $left = $intron->[0] - $b;
			my $right = $intron->[1] - $b;
			
			for (my $i = $left; $i <= $right; $i++) {
                                if ($array[$i] == 0) { $array[$i] = 2;}
                                elsif ($array[$i] == 1) {
                                        $array[$i] = -1;
                                }
                        }
                }
	}

	my @conserved;
	my $left = -1;
	for (my $i = 0; $i <= $size-1; $i ++) {
		next if $array[$i] == -1 || $array[$i] == 0;
		if ($left == -1) { $left = $i; }
		if (!defined $array[$i+1] || $array[$i+1] != $array[$i]) {
			push @conserved, { 'b'=> $left + $b, 'e'=>$i + $b, 
					   'type' => $array[$i],
					 };
			$left = -1;
		}
	}

	my @extended_hits;
	foreach my $hit (@$bag) {
		my @left_additional_hsps;
		my @right_additional_hsps;
		my @new_hsps;

		my $first_hsp = $hit->{hsps}->[0];
		my $last_hsp  = $hit->{hsps}->[scalar @{$hit->{hsps}}-1];

		my $new_left = $first_hsp->[0];
		my $new_right = $last_hsp->[1];

		if ($array[$first_hsp->[0] - $b] == 1) {	
			my $first_tag = 1;
			my $left_coor = $first_hsp->[0];

			for (my $k = scalar @conserved -1; $k>=0; $k--) {
				my $region = $conserved[$k];

				next unless $region->{b} <= $first_hsp->[0];
				if ($first_tag == 1) {

					$new_left = $region->{b};
					$left_coor = $region->{b} - 1;
					$first_tag = 0;
				}

				else {
					last if $left_coor < $b || $array[$left_coor-$b] == -1
						|| $array[$left_coor-$b] == 0;

					if ($array[$left_coor-$b] == 1) {
						unshift @left_additional_hsps, [$region->{b}, $region->{e}];
					}

					$left_coor = $region->{b} - 1;
				}
			}
		}

		if ($array[$last_hsp->[1] - $b] == 1) {
                        my $first_tag = 1;
                        my $right_coor = $last_hsp->[1];

                        for (my $k = 0; $k <= scalar @conserved -1; $k++) {
                                my $region = $conserved[$k];

                                next unless $region->{e} >= $last_hsp->[1];
                                if ($first_tag == 1) {

                                        $new_right = $region->{e};
                                        $right_coor = $region->{e} + 1;
                                        $first_tag = 0;
                                }       

                                else {
                                        last if  $right_coor > $e|| $array[$right_coor-$b] == -1
                                                || $array[$right_coor-$b] == 0;

                                        if ($array[$right_coor-$b] == 1) {
                                                push @right_additional_hsps, [$region->{b}, $region->{e}];
                                        }

                                        $right_coor = $region->{e} + 1;
                                }
                        }
                }

		if (scalar @{$hit->{hsps}} == 1) {
			push @new_hsps, @left_additional_hsps;
			push @new_hsps, [$new_left, $new_right];
			push @new_hsps, @right_additional_hsps;
		}
	
		else {	
			push @new_hsps, @left_additional_hsps;
			push @new_hsps, [$new_left, $first_hsp->[1]];

			if (scalar @{$hit->{hsps}} > 2) {
				for (my $l = 1; $l <= scalar @{$hit->{hsps}} - 2; $l ++) {
					push @new_hsps, $hit->{hsps}->[$l];
				}
			}

			push @new_hsps, [$last_hsp->[0], $new_right];
			push @new_hsps, @right_additional_hsps;
		}
			
		my $new_hit =  {  'b' => $new_hsps[0]->[0],
       				  'e' => $new_hsps[scalar @new_hsps -1]->[1],
			          'hsps' => \@new_hsps,
			       	  'strand' => $hit->{strand},
				};
		evaluator::pseudo_hit::add_splice($new_hit);
	
		push @extended_hits, $new_hit;

	}

	return \@extended_hits;
}	
				
#------------------------------------------------------------------------
1;
