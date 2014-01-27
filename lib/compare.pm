#------------------------------------------------------------------------
#----                            compare                             ---- 
#------------------------------------------------------------------------
package compare;
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
sub is_same_alt_form {
        my $a     = shift;
        my $b     = shift;
	my $flank = shift;

	die "only one feature in is_same_alt_form!\n"
	unless defined($a) && defined($b);

	return 0 if($a->num_hsps != $b->num_hsps);

	my ($s_to_a_str, $s_to_b_str) = compare_by_shadow($a, $b, $flank);

	my $a_to_b_str = compare_phat_hits($a, $b, 'query', $flank);

	if($s_to_a_str eq $s_to_b_str && $a_to_b_str !~ /^0+$/){
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
	elsif($a_to_b_str =~ /^(0+|Z+|z+)$/){
	   return 0;
	}
	else {
	        #print STDERR "not caught s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str ".$a->name." ".$b->name."\n";
		#sleep 3;

		return 1;
	}
}
#------------------------------------------------------------------------
sub is_redundant_alt_form {
        my $a     = shift;
        my $b     = shift;
        my $flank = shift;

       die "only one feature in is_redundant_form!\n"
        unless defined($a) && defined($b);

	# note that b will always have fewer exons or be shorter...

        my ($s_to_a_str, $s_to_b_str) = compare_by_shadow($a, $b, $flank);

        my $a_to_b_str = compare_phat_hits($a, $b, 'query', $flank);
	my $b_to_a_str = compare_phat_hits($b, $a, 'query', $flank);
        if     ($s_to_a_str eq $s_to_b_str && $a_to_b_str !~ /^0*$/){
                return 1;
        }
        elsif ($b_to_a_str =~ /^B?1+b?$/ && $a_to_b_str =~ /^0*A?1+a?0*$/){
                #print STDERR "RAAAA s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
                #sleep 3;
                return 1;
        }
        elsif ($s_to_a_str =~ /^1+$/ && $s_to_b_str =~ /^0*A?1+a?0*$/){
                #print STDERR "RBBBB s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str\n";
	        #sleep 3;
                #return 1;
        }
        else {
                #print STDERR "redun not caught s_to_a_str:$s_to_a_str s_to_b_str:$s_to_b_str a_to_b_str:$a_to_b_str b_to_a_str:$b_to_a_str".$a->name." ".$b->name."\n";
                #sleep 3;
                return 0;
        }
}
#------------------------------------------------------------------------
sub hsps_overlap {
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

		#fix for flank
		$aB -= $r;
		$aB = 1 if($aB < 1);
		$aE += $r;
                foreach my $hsp_b (@{$sorted_b}){
                        my $bB = $hsp_b->nB($what);
                        my $bE = $hsp_b->nE($what);

			($bB, $bE) = ($bE, $bB) if $bB > $bE;
                        my $class = compare($aB, $aE, $bB, $bE);

			return 1 if $class ne '0';			
		}
	}
	return 0;
}
#------------------------------------------------------------------------
sub overlap {
       my $hit_a = shift;
       my $hit_b = shift;
       my $what  = shift;
       my $r     = shift || 0;

       my $aB = $hit_a->nB($what);
       my $aE = $hit_a->nE($what);
       
       ($aB, $aE) = ($aE, $aB) if $aB > $aE;
       
       #fix for flank
       $aB -= $r;
       $aB = 1 if($aB < 1);
       $aE += $r;

       my $bB = $hit_b->nB($what);
       my $bE = $hit_b->nE($what);
       
       ($bB, $bE) = ($bE, $bB) if $bB > $bE;
       my $class = compare($aB, $aE, $bB, $bE);
       
       return 1 if $class ne '0';
	   
       return 0;
}
#------------------------------------------------------------------------
sub compare_by_shadow {
        my $a     = shift;
        my $b     = shift;
	my $flank = shift;

        my $a_coors  = PhatHit_utils::get_hsp_coors_from_hit($a, 'query');
	my $b_coors  = PhatHit_utils::get_hsp_coors_from_hit($b, 'query');

	my @t_coors = (@{$a_coors}, @{$b_coors});

        my $pieces  = Shadower::getVectorPieces(\@t_coors, 0);

	my $s_code_a = get_shadow_code($pieces, $a_coors, $flank);
	my $s_code_b = get_shadow_code($pieces, $b_coors, $flank);

	return ($s_code_a, $s_code_b);
}
#------------------------------------------------------------------------
sub get_shadow_code {
    my $pieces = shift;
    my $coors  = shift;
    my $r      = shift;

    @$coors = sort {$a->[0] <=> $b->[0]} @$coors;

    my @codes;
    my $i = 0;
    my $first = 0;
    foreach my $p (@{$pieces}){
	#assume correct order
	my $pB = $p->{b};
	my $pE = $p->{e}; 
	
	my $flag;
	$codes[$i] = 0; #initialize
	for(my $j = $first; $j < @{$coors}; $j++){
	    my $pair = $coors->[$j];
	    my $cB = $pair->[0];
	    my $cE = $pair->[1];
	    ($cB, $cE) = ($cE, $cB) if($cB > $cE);

	    last if($pE < $cB);
	    if(!$flag){ #itterate until the first match
		if($pB <= $cE){
		    $first = $j;
		    $flag++;
		}
		else{
		    $first = $j+1;
		    next;
		}
	    }

	    my $c = compare($pB, $pE, $cB, $cE, $r);
	    next if(!$c);
	    $codes[$i] = $c if(code_values($c) > code_values($codes[$i]));
	    last if($c eq '1');
	}
	$i++;
    }
    
    my $code = join('', @codes);
    
    return $code;
}
#------------------------------------------------------------------------
sub simplify_comparison_code {
	my $codes = shift;

	my $code = '';
	foreach my $exon (@{$codes}){
	    $code .= 0 if(!$exon || !@$exon);
	    my @sorted  = sort {code_values($b) <=> code_values($a)} @{$exon};
	    
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
		return 0;
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
#assumes same strand
sub compare_phat_hits {
    my $hit_a = shift;
    my $hit_b = shift;
    my $what  = shift;
    my $r     = shift || 0;
    
    my $sorted_a = $hit_a->sortFeatures($what);
    my $sorted_b = $hit_b->sortFeatures($what);
    
    my @codes;
    my $i = 0;
    my $first = 0;
    foreach my $hsp_a (@{$sorted_a}){
	my $aB = $hsp_a->start($what);
	my $aE = $hsp_a->end($what);

	my $flag;
	$codes[$i] = 0; #initialize
	for(my $j = $first; $j < @{$sorted_b}; $j++){
	    my $hsp_b = $sorted_b->[$j];
	    my $bB = $hsp_b->start($what);
	    my $bE = $hsp_b->end($what);

	    last if($aE < $bB);
	    if(!$flag){ #itterate until the first match
		if($aB <= $bE){
		    $first = $j;
		    $flag++;
		}
		else{
		    $first = $j+1;
		    next;
		}
	    }

	    my $c = compare($aB, $aE, $bB, $bE, $r);
	    next if(!$c);
	    $codes[$i] = $c if(code_values($c) > code_values($codes[$i]));
	    last if($c eq '1');
	}
	$i++;
    }
    
    my $code = join('', @codes);
    
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
	if(!same_strand($aB, $aE, $bB, $bE)){
	    return 'R';
	}

	my $s = 1;
	if ($aB > $aE && $bB > $bE){
		($aB, $aE) = ($aE, $aB);
		($bB, $bE) = ($bE, $bB);
		$s = -1;
	}
		

	if($bE < $aB || $bB > $aE){
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
1;


