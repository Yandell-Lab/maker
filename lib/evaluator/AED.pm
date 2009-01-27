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
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub txnAED {
	my $box = shift;
	my $parameter = shift;

	my $t =		$box->{transcript};
	my $pol_e =	$box->{est};
	my $pol_p =	$box->{exonerate};
	my $seq =	$box->{seq};

	my $bag = combine($pol_e, $pol_p);

	my $same_alt_clusters = combine_same_alt_forms($bag, $seq);

	my $min_AED = 1;
	foreach my $cluster (@$same_alt_clusters) {
		my $AED = cluster_AED($t, $cluster, $seq, $parameter);
		$min_AED = $AED if $AED < $min_AED;
	}

	return $min_AED;
}
#------------------------------------------------------------------------------
sub gene_AED {
	my $c = shift; # The gene object;
	my $pol_e_hits = shift;
	my $pol_p_hits = shift;
	my $blastx_hits = shift;
	my $seq = shift;

	my $annotation_coors = PhatHit_utils::get_hsp_coors($c, 'query');
	my $annotation_pieces = Shadower::getPieces($seq, 
				$annotation_coors, 0);
	my $gene_strand = $c->[0]->strand;
	my $gene_boundary = get_boundary($annotation_pieces);

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

	my @good_hits;
	foreach my $hit (@bag) {
		push @good_hits, $hit if overlap($gene_boundary, 
					$gene_strand, $hit);
	}

	my $evi_coors = PhatHit_utils::get_hsp_coors(\@good_hits, 'query');
	my $evi_pieces= Shadower::getPieces($seq, $evi_coors, 0);

	my $coorA = get_coor($annotation_pieces);
	my $coorB = get_coor($evi_pieces);

	my $para = {'exon'=>1};

	my $AED = cal($coorA, $coorB, $para);
	return $AED;
}
	
#------------------------------------------------------------------------------
sub cal {
	my $coor1 = shift; # Annotation coordinations
	my $coor2 = shift; # Evidence coordinations
	my $para  = shift;

	if (!defined $para->{exon}) { $para->{exon} = 0;}
	if (!defined $para->{intron}) { $para->{intron} = 0;}
	if (!defined $para->{donor}) { $para->{donor} = 0;}
	if (!defined $para->{acceptor}) { $para->{acceptor} = 0;}
	if (!defined $para->{start}) { $para->{start} = 0;}
	if (!defined $para->{stop}) { $para->{stop} = 0;}
	

	die "Annotation set cannot be empty!\n" 
		if @$coor1 == 0 ;

	return 1 if @$coor2 == 0;

	my (@array1, @array2);

	# Creating two arrays containing the two annotations, in which 
	# -1 means UTR, 0 means introns and 1 means exons.
	foreach my $coor (@$coor1) {
		my $start = $coor->[0] - 1;
		my $end   = $coor->[1] - 1;

		for(my $i = $start; $i <= $end; $i++) {
			$array1[$i] = 1;
		}
	}

        foreach my $coor (@$coor2) {
                my $start = $coor->[0] - 1;
                my $end   = $coor->[1] - 1;

                for(my $i = $start; $i <= $end; $i++) {
                        $array2[$i] = 1;
                }
        }

        my $start_point = 0; 
        for (my $i= 0; ; $i++) {
                last if defined $array1[$i] || defined $array2[$i];
                $start_point++; 
		$array1[$i] = -1;
		$array2[$i] = -1;
        }

	for (my $i= $start_point; $i<=$#array1; $i++) { 
		$array1[$i] = 0 unless defined $array1[$i]; }
	for (my $i= $start_point; $i<=$#array2; $i++) { 
		$array2[$i] = 0 unless defined $array2[$i]; }
	

	# Make UTRs
	for(my $i = $start_point; ; $i++) {
		last if $array1[$i] == 1;
		$array1[$i] = -1;
	}
	for(my $i = $#array1; ; $i--) {
		last if $array1[$i] == 1;
		$array1[$i] = -1;
	}
        for(my $i = $start_point; ; $i++) {
                last if $array2[$i] == 1;
                $array2[$i] = -1;
        }
        for(my $i = $#array2; ; $i--) {
                last if $array2[$i] == 1;
                $array2[$i] = -1;
        }


	# Get the coordination of start codon, stop codon, donor site,
	# and acceptor sites into two hashes: %special1 and %special2.
	my (%special1, %special2);
	my ($score1, $score2) = (0,0);
	for (my $i = $start_point; $i <= $#array1; $i++) {
		if ($array1[$i] == 1) { $score1 += $para->{exon}; }
		elsif ($array1[$i] == 0) {
			$score1 += $para->{intron};
		}
		
		if 	($i == 0 && $array1[$i] == 1) {	
			$special1{$i} = 'start'; }
		elsif ($i == $#array1 && $array1[$i] == 1) { 
			$special1{$i} = 'stop'; 
		}
		elsif ($array1[$i] == 1) {
			if ($array1[$i-1] == -1) { $special1{$i} = 'start';}
			elsif ($array1[$i-1] == 0) { $special1{$i} = 'acceptor';}
			elsif ($array1[$i+1] == -1) { $special1{$i} = 'stop';}
			elsif ($array1[$i+1] == 0) { $special1{$i} = 'donor';}
		}
	}
        for (my $i = $start_point; $i <= $#array2; $i++) {
                if ($array2[$i] == 1) { $score2 += $para->{exon}; }
                elsif ($array2[$i] == 0) {
                        $score2 += $para->{intron};
                }

                if      ($i == 0 && $array2[$i] == 1) {        
                        $special2{$i} = 'start'; }
                elsif ($i == $#array2 && $array2[$i] == 1) { 
                        $special2{$i} = 'stop'; 
                }       
                elsif ($array2[$i] == 1) {
                        if ($array2[$i-1] == -1) { $special2{$i} = 'start';}
                        elsif ($array2[$i-1] == 0) { $special2{$i} = 'acceptor';}
                        elsif ($array2[$i+1] == -1) { $special2{$i} = 'stop';}
                        elsif ($array2[$i+1] == 0) { $special2{$i} = 'donor';}
                }
        }

	while (my ($position, $type) = each %special1) {
		if ($type eq 'start') { $score1 += $para->{start};}
		elsif ($type eq 'stop') {
			$score1 += $para->{stop};
		}
		elsif ($type eq 'donor') {
			$score1 += $para->{donor};
		}
		elsif ($type eq 'acceptor') {
			$score1 += $para->{acceptor};
		}
	}

        while (my ($position, $type) = each %special2) {
                if ($type eq 'start') { $score2 += $para->{start};}
                elsif ($type eq 'stop') {
                        $score2 += $para->{stop};
                }
                elsif ($type eq 'donor') {
                        $score2 += $para->{donor};
                }
                elsif ($type eq 'acceptor') {
                        $score2 += $para->{acceptor};
                }
        }

	my ($o_exon, $o_intron, $o_start, $o_stop, $o_donor, $o_acceptor)=(0,0,
		0,0,0,0);

	for (my $i = $start_point; $i <= $#array1 && $i <= $#array2; $i++) {
		if ($array1[$i] == 0 && $array2[$i] == 0) {
			$o_intron ++;
		}
		elsif ($array1[$i] ==1 && $array2[$i] ==1) {
			$o_exon ++;
		}
	}
	
	while (my ($position, $type1) = each %special1) {
		if (defined $special2{$position}) {
			my $type2 = $special2{$position};
			if ($type1 eq 'start' && $type2 eq 'start') {
				$o_start ++;
			}
			elsif ($type1 eq 'stop' && $type2 eq 'stop') {
				$o_stop ++;
			}
			elsif ($type1 eq 'donor' && $type2 eq 'donor') {
				$o_donor ++;
			}
			elsif ($type1 eq 'acceptor' && $type2 eq 'acceptor') {
		 		$o_acceptor ++;
			}
		}
	}

	my $overlap_score = $o_exon * $para->{exon} + $o_intron * $para->{intron}
			+ $o_start * $para->{start} + $o_stop   * $para->{stop}
			+ $o_donor * $para->{donor} + $o_acceptor * $para->{acceptor};

	my $congruence = 0.5 * ($overlap_score/$score1 + $overlap_score/$score2);

	return 1- $congruence;
}
		

#-------------------------------------------------------------------------------
sub cluster_AED {
	my $t = shift;
	my $cluster = shift;
	my $seq = shift;
	my $para= shift;

	return 1 if (scalar @$cluster == 0) ||
			$t->strand ne $cluster->[0]->strand;

	my $annotation_coors = PhatHit_utils::get_hsp_coors([$t], 'query');
	my $annotation_pieces= Shadower::getPieces($seq, $annotation_coors, 0);

	my $evi_coors = PhatHit_utils::get_hsp_coors($cluster, 'query');
	my $evi_pieces= Shadower::getPieces($seq, $evi_coors, 0);

	my $coorA = get_coor($annotation_pieces);
	my $coorB = get_coor($evi_pieces);

	my $AED = cal($coorA, $coorB, $para);
	return $AED;
}

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub overlap {
	my $gene_boundary = shift;
	my $gene_strand   = shift;
	my $hit = shift;

	return 0 unless $hit->strand eq $gene_strand;

	my $b = $hit->nB('query');
	my $e = $hit->nE('query');

	($b, $e) = ($e, $b) unless $b<=$e;

	return 0 if $b > $gene_boundary->[1] || $e < $gene_boundary->[0];
	return 1;
}
	
#-------------------------------------------------------------------------------
sub get_coor {
	my $pieces = shift;

	my @coors;

	foreach my $piece (@$pieces) {
                my $b = $piece->{b};
                my $e = $piece->{e};
                
                ($b, $e) = ($e, $b) if $b > $e;
		push @coors, [$b, $e];
	}
	return \@coors;
}
#-------------------------------------------------------------------------------
sub get_boundary {
	my $pieces = shift;

	my ($left, $right);
	foreach my $piece (@$pieces) {
		my $b = $piece->{b};
		my $e = $piece->{e};

		($b, $e) = ($e, $b) if $b > $e;

		$left = $b if (not defined $left) || $left > $b;
		$right= $e if (not defined $right)|| $right< $e;
	}
	return [$left, $right];
}
#-------------------------------------------------------------------------------
#------------------------------------------------------------------------
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
sub combine_same_alt_forms {
	my $bag = shift;
	my $seq = shift;

	my @clusters;

	foreach my $hit (@$bag) {
		my @clusters_to_be_added;
		foreach my $cluster (@clusters) {
			if (hit_agree_hits_alt_form($hit, $cluster, $seq)) {
				my @new_cluster = @$cluster;
				push @new_cluster, $hit;
				push @clusters_to_be_added, \@new_cluster;
			}
		}
		push @clusters, @clusters_to_be_added;
		push @clusters, [$hit];
	}


	## Maybe I need to get rid of the redundant clusters and 
	## short clusters here as well!	

	return \@clusters;
}
	
#------------------------------------------------------------------------
sub hit_agree_hits_alt_form {
	my ($hit, $cluster, $seq ) = @_;

	if ($hit->strand ne $cluster->[0]->strand) { return 0; }
	
	foreach my $one_hit_in_cluster (@$cluster) {
		return 0 unless compare::is_same_alt_form($hit, 
				$one_hit_in_cluster, $seq, 0);
	}
	
	return 1;
}

	
#------------------------------------------------------------------------
1;


