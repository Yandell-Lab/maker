#------------------------------------------------------------------------
#----                 evaluator::pseudo_hit.pm                       ---- 
#------------------------------------------------------------------------
package evaluator::pseudo_hit;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use cluster;
use clean;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub convert_to_pseudo_hit {
	my $hit = shift;

	my $strand = $hit->strand('query');

	my $nB = $hit->nB('query');
	my $nE = $hit->nE('query');
	($nB, $nE) = ($nE, $nB) if ($nB > $nE);

	my @hsps;
	foreach my $hsp ($hit->hsps) {
		my $start = $hsp->nB('query');
		my $end	  = $hsp->nE('query');
		($start, $end) = ($end, $start) if $start > $end;
		
		push @hsps, [$start, $end];
	}

	my @sorted_hsps = sort {$a->[0] <=> $b->[0]} @hsps;
	
	my $combined_hit = { 'b'		=> $nB,
		    'e'		=> $nE,
		    'strand'	=> $strand,
		    'hsps'	=> \@sorted_hsps,
		  };

	return $combined_hit;
}
#------------------------------------------------------------------------
sub combine_pseudo_hit {
	my $pseu_hits = shift;

	my $strand = $pseu_hits->[0]->{strand};

	my @hsps;
	my ($b, $e) = ($pseu_hits->[0]->{b}, $pseu_hits->[0]->{e});;
		
	foreach my $ph (@$pseu_hits) {
		die "Cannot combine hits in different strands!" if
			$ph->{strand} != $strand;
		push @hsps, @{$ph->{hsps}};
		$b = $ph->{b} if $b > $ph->{b};
		$e = $ph->{e} if $e < $ph->{e};
	}

	my $combined_hsps = combine_pieces(\@hsps);

	my $combined_ph = {'b'	=> $b,
			   'e'  => $e,
			   'strand' => $strand,
			   'hsps' => $combined_hsps,
			  };

	return $combined_ph;
}


#------------------------------------------------------------------------
sub is_same_alt_form {
        my $a     = shift;
        my $b     = shift;
	my $flank = shift;

	die "only one feature in is_same_alt_form!\n"
	unless defined($a) && defined($b);

	if ($a->{strand} != $b->{strand}) {return 0;}
	if ($a->{b} >$b->{e} || $a->{e} < $b->{b}) {return 0;}


	my ($s_to_a_str, $s_to_b_str) = compare_by_shadow($a, $b, $flank);

	my $a_to_b_str = compare_phat_hits($a, $b, $flank);

	if     ($s_to_a_str eq $s_to_b_str && $a_to_b_str !~ /^0+$/){
		return 1;
	}

        elsif ($s_to_a_str =~ /[^0]0+[^0]/){
		#print STDERR "QAAAA s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
		#print STDERR $a->name." ".$b->name."\n";
		#sleep 3;
                return 0;
        }
        elsif ($s_to_b_str =~ /[^0]0+[^0]/){
                #print STDERR "QBBBB s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
		#print STDERR $a->name." ".$b->name."\n";
                #sleep 3;

                return 0;
        }
	else {
		#print STDERR "not caught s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str".$a->name." ".$b->name."\n";
		#sleep 3;
		return 1;
	}
}
#------------------------------------------------------------------------
sub is_redundant_alt_form {
        my $a     = shift;
        my $b     = shift;
        my $seq   = shift;
        my $flank = shift;

       die "only one feature in is_redundant_form!\n"
        unless defined($a) && defined($b);

	# note that b will always have fewer exons or be shorter...

	return 0 if $a->{strand} != $b->{strand};
        my ($s_to_a_str, $s_to_b_str) = compare_by_shadow($a, $b, $seq, $flank);

        my $a_to_b_str = compare_phat_hits($a, $b, 'query', $flank);
	my $b_to_a_str = compare_phat_hits($b, $a, 'query', $flank);
        if     ($s_to_a_str eq $s_to_b_str && $a_to_b_str !~ /^0*$/){
                return 1;
        }
        elsif ($b_to_a_str =~ /^B?1+b?$/ && $a_to_b_str =~ /^0*A?1+a?$/){
                #print STDERR "RAAAA s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
                #sleep 3;
                return 1;
        }
        elsif ($s_to_a_str =~ /^1+$/ && $s_to_b_str =~ /^0*A?1+a?$/){
                #print STDERR "RBBBB s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
                sleep 3;
                #return 1;
        }

        else {
                #print STDERR "redun not caught s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str a_to_b_str:$a_to_b_str b_to_a_str:$b_to_a_str".$a->name." ".$b->name."\n";
                #sleep 3;
                return 0;
        }
}
#------------------------------------------------------------------------
sub overlap {
       my $hit_a = shift;
       my $hit_b = shift;
       my $what  = shift;
       my $r     = shift || 0;


        my $sorted_a = PhatHit_utils::sort_hits($hit_a, $what);
        my $sorted_b = PhatHit_utils::sort_hits($hit_b, $what);

        foreach my $hsp_a (@{$sorted_a}){

                my $aB = $hsp_a->nB($what);
                my $aE = $hsp_a->nE($what);

		($aB, $aE) = ($aE, $aB) if $aB > $aE;

                foreach my $hsp_b (@{$sorted_b}){
                        my $bB = $hsp_b->nB($what);
                        my $bE = $hsp_b->nE($what);

			($bB, $bE) = ($bE, $bB) if $bB > $bE;
                        my $class = compare($aB, $aE, $bB, $bE, $r);

			return 1 if $class ne '0'
			
		}
	}
	return 0;
}
#------------------------------------------------------------------------
sub compare_by_shadow {
        my $a     = shift;
        my $b     = shift;
	my $flank = shift;

        my $a_coors  = $a->{hsps};
	my $b_coors  = $b->{hsps};

	my $combined_ph = combine_pseudo_hit([$a, $b]);

	my $s_code_a = get_shadow_code($combined_ph->{hsps}, $a_coors, $flank);
	my $s_code_b = get_shadow_code($combined_ph->{hsps}, $b_coors, $flank);

	return ($s_code_a, $s_code_b);
}
#------------------------------------------------------------------------
sub get_shadow_code {
	my $combined_coors = shift;
	my $coors  = shift;
	my $flank  = shift;

	my @codes;
	my $i = 0;
	foreach my $p (@{$combined_coors}){
		my $pB = $p->[0];
		my $pE = $p->[1];

		foreach my $pair (@{$coors}){
			my $cB = $pair->[0];
			my $cE = $pair->[1];

			($cB, $cE) = ($cE, $cB) if $cB > $cE;
			
			my $class = compare($pB, $pE, $cB, $cE, $flank);

			push(@{$codes[$i]}, $class);
		}
		$i++;
	}

	my $code = simplify_comparison_code(\@codes);


	return $code;
}
#------------------------------------------------------------------------
sub simplify_comparison_code {
	my $codes = shift;

	my $code = '';
	foreach my $exon (@{$codes}){
		my @sorted  = 
		sort {code_values($b) <=> code_values($a)} @{$exon};

		$code .= shift @sorted;
	}
	return $code;
}
#------------------------------------------------------------------------
sub code_values {
	my $v = shift;

	if    ($v eq '0'){
		return 0;
	}
	elsif ($v eq  '1'){
		return 10;
	}
	elsif ($v =~  /[ABab]/) {
		return 5;
	}
	else {
		return 1;
        }

}
#------------------------------------------------------------------------
sub same_strand {
        my $aB = shift;
        my $aE = shift;
        my $bB = shift;
        my $bE = shift;

	if    ($aB <= $aE && $bB <= $bE){
		return 1;
	}
	elsif ($aB >= $aE && $bB >= $bE){
		return 1;
	}
	else {
		return 0;
	}
}
#------------------------------------------------------------------------
sub compare_phat_hits {
        my $hit_a = shift;
        my $hit_b = shift;
	my $r     = shift || 0;

        my $sorted_a = sort_hits($hit_a);
        my $sorted_b = sort_hits($hit_b);

	my $alpha = 0;
	my $omega = @{$sorted_a} - 1;

	my @codes;
	my $i = 0;
        foreach my $hsp_a (@{$sorted_a}){
			
                my $aB = $hsp_a->[0];
                my $aE = $hsp_a->[1];


		my $j = 0;
                foreach my $hsp_b (@{$sorted_b}){
                        my $bB = $hsp_b->[0];
                        my $bE = $hsp_b->[1];
			
			my $class = compare($aB, $aE, $bB, $bE, $r);

			push(@{$codes[$i]}, $class);

                }
		$i++;
        }
	my $code = simplify_comparison_code(\@codes);

	return $code;
}
#------------------------------------------------------------------------
sub compare {
	my $aB = shift;
	my $aE = shift;
	my $bB = shift;
	my $bE = shift;
	my $r  = shift || 0;

	my $class;
	if(same_strand($aB, $aE, $bB, $bE)){

	}
	else {
		return 'R';
	}

	my $s = 1;
	if ($aB > $aE && $bB > $bE){
		($aB, $aE) = ($aE, $aB);
		($bB, $bE) = ($bE, $bB);
		$s = -1;
	}
		

	if     ($bE < $aB || $bB > $aE){
		# a        ----
		# b  ---- 
		$class = 0;
	}
	elsif (abs($aB - $bB) <= $r && abs($aE - $bE) <= $r){
		#  a ----- 
		#  b -----
		$class = 1;
	}
	elsif ($bB < $aB && abs($aE - $bE) <= $r){
		# a   ----
		# b ------
		$class = $s == 1 ? 'B' : 'b'; 
	} 
        elsif ($aB < $bB && abs($aE - $bE) <= $r){
                # a ------
                # b   ----
                $class = $s == 1 ? 'A' : 'a';
        }
	elsif (abs($aB - $bB) <= $r && $aE < $bE){
		# a ----
		# b ------
		$class = $s == 1 ? 'b' : 'B';
	}
        elsif (abs($aB - $bB) <= $r && $aE > $bE){
                # a ------
                # b ----
                $class = $s == 1 ? 'a' : 'A'; 
        }
	elsif ($aB > $bB && $aE < $bE){
                # a   ----
                # b --------
                $class = 'I'; # I for In
	}
        elsif ($aB < $bB && $aE > $bE){
                # a --------
                # b   ----
                $class = 'i'; # i for in
        }
        elsif ($aB < $bB && $aE < $bE){
                # a ------
                # b   -------
                $class = $s == 1 ? 'Z' : 'z'; ; # Z for zig-zag 
        }
        elsif ($aB > $bB && $aE > $bE){
                # a      ------
                # b -------
                $class = $s == 1 ? 'z' : 'Z'; # z for Zig-zag
        }

	return $class;	
}
#------------------------------------------------------------------------
sub sort_hits {
	my $hit = shift;
	
	my $hsps = $hit->{hsps};
	my @sorted = sort {$a->[0] <=> $b->[0] } 
			@$hsps;

	return \@sorted;
}
#------------------------------------------------------------------------
sub combine_pieces {
	my $coors = shift;

	my ($start, $end) = @{$coors->[0]};

	foreach my $coor (@$coors) {
		$start = $coor->[0] if $coor->[0] < $start;
		$end   = $coor->[1] if $coor->[1] > $end;
	}

	my $length = $end - $start + 1;
	
	my @array;
	for (my $i = 0; $i <= $length - 1; $i ++) {
		$array[$i] = 0;
	}

	foreach my $coor (@$coors) {
		for (my $i = $coor->[0]-$start; $i <= $coor->[1]-$start; $i ++) {
			$array[$i] = 1;
		}
	}

	my @result;
	my ($left, $right);
	for (my $i = 0; $i <= $length -1; $i ++) {
		next if $array[$i] == 0;

		if ($i == 0 || $array[$i-1] == 0) {$left = $i;}
		if ($i == $length-1 || $array[$i+1] == 0) {
			$right = $i;
			push @result, [$left, $right];
		}
	}

	foreach my $coor (@result) {
		$coor->[0] += $start;
		$coor->[1] += $start;
	}

	return \@result;
}

#------------------------------------------------------------------------
sub are_same_pseudo_hits {
	my $ph_a = shift;
	my $ph_b = shift;

	return 0 if $ph_a->{strand} != $ph_b->{strand};
	return 0 if $ph_a->{b} != $ph_b->{b};
	return 0 if $ph_a->{e} != $ph_b->{e};

	my @hsps_a = sort {$a->[0] <=> 
			   $b->[0] } @{$ph_a->{hsps}};

        my @hsps_b = sort {$a->[0] <=> 
                           $b->[0] } @{$ph_b->{hsps}};

	return 0 if (scalar @hsps_a) != (scalar @hsps_b);

	for(my $i = 0; $i <= $#hsps_a; $i ++) {
		my $coor_a = $hsps_a[$i];
		my $coor_b = $hsps_b[$i];

		return 0 if $coor_a->[0] != $coor_b->[0] 
			 || $coor_a->[1] != $coor_b->[1];
	}

	return 1;
}
#------------------------------------------------------------------------
1;


