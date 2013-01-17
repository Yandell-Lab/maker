#------------------------------------------------------------------------------
#----                            SimpleCluster.pm                          ---- 
#------------------------------------------------------------------------------
package SimpleCluster;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub cluster_hits {
    my $hits = shift;
    my $flank = shift;

    my $pairs = build_hit_pairs($hits, $flank);
    my $cMap = cluster_pairs($pairs);

    my @clusters;
    for(my $i = 0; $i < @$cMap; $i++){
        my $members = $cMap->[$i];
        next if(!$members || !@$members);

        my @array = map {$hits->[$_]} @{$members};

	push(@clusters, \@array);
    }
    
    return \@clusters;
}
#-------------------------------------------------------------------------------
#generic form - requires an array of hashes or other features with the keys
#start, end, and stand (optional) defined
sub cluster {
    my $array = shift;
    my $flank = shift;

    my $pairs = build_pairs($array, $flank);
    my $cMap = cluster_pairs($pairs);

    my @clusters;
    for(my $i = 0; $i < @$cMap; $i++){
        my $members = $cMap->[$i];
        next if(!$members || !@$members);

        my @set = map {$array->[$_]} @{$members};

	push(@clusters, \@set);
    }
    
    return \@clusters;
}
#-------------------------------------------------------------------------------
sub build_hit_pairs {
    my $hits = shift;
    my $flank = shift;
    my $ignore_strand = shift || 0;

    my @pairs;

    for (my $i = 0; $i < @$hits; $i++) {
	my $aB = $hits->[$i]->start('query');
	my $aE = $hits->[$i]->end('query');
	my $aSt = $hits->[$i]->strand('query');

	#fix for flank
	$aB -= $flank;
	$aB = 1 if($aB < 1);
	$aE += $flank;

	push(@pairs, [$i, $i]);

	for (my $j = $i+1; $j < @$hits; $j++) {
	    my $bB = $hits->[$j]->start('query');
	    my $bE = $hits->[$j]->end('query');
	    my $bSt = $hits->[$j]->strand('query');

	    next if($aSt ne $bSt && !$ignore_strand);

	    my $code = compare::compare($aB, $aE, $bB, $bE);
	    push(@pairs, [$i, $j]) if($code);
	}
    }

    return \@pairs;
}
#-------------------------------------------------------------------------------
sub build_pairs {
    my $array = shift;
    my $flank = shift;
    my $jaccard = shift;

    my @pairs;
    my @matrix_h;
    for (my $i = 0; $i < @$array; $i++) {
	my $aB = $array->[$i]->{start};
	my $aE = $array->[$i]->{end};
	my $aSt = $array->[$i]->{strand};
	push(@pairs, [$i, $i]);

	#fix for flank
	$aB -= $flank;
	$aB = 1 if($aB < 1);
	$aE += $flank;

	for (my $j = $i+1; $j < @$array; $j++) {
	    my $bB = $array->[$j]->{start};
	    my $bE = $array->[$j]->{end};
	    my $bSt = $array->[$j]->{strand};

	    next if(defined($aSt) && defined ($bSt) && ($aSt ne $bSt));

	    my $code = compare::compare($aB, $aE, $bB, $bE);
	    if($code){
		push(@pairs, [$i, $j]);
		if($jaccard){
		    $matrix_h[$i]->{$j}++;
		    $matrix_h[$j]->{$i}++;
		}
	    }
	}
    }

    if($jaccard){
	my @matrix;
	foreach my $m (@matrix_h){
	    push(@matrix, [keys %{$m}]);
	}

	my @n_pairs;
	for (my $i = 0; $i < @pairs; $i++){
	    my ($j, $k) = @{$pairs[$i]};

	    if($j == $k){
		push(@n_pairs, $pairs[$i]);
		next;
	    }

	    my $j_size = @{$matrix[$j]};
	    my $k_size = @{$matrix[$k]};

	    if($jaccard == 1 && $j_size != $k_size){
		next;
	    }

	    my $int = (grep {exists ($matrix_h[$j]{$_})} @{$matrix[$k]}) + 1; #see comment below
	    #you must add 1 in order to include the edge between the query pair
	    #i.e. edge A->B will not get counted as inersecting  because 'A' != 'B'

	    my $union = ($j_size + $k_size) - $int;

	    if ($int/$union >= $jaccard){
		push(@n_pairs, $pairs[$i]);
	    }
	}

	return \@n_pairs;
    }

    return \@pairs;
}
#-------------------------------------------------------------------------------
sub cluster_on_hsps {
    my $hits = shift;
    my $flank = shift;

    my $pairs = build_hsp_pairs($hits, $flank);
    my $cMap = cluster_pairs($pairs);

    my @clusters;
    for(my $i = 0; $i < @$cMap; $i++){
        my $members = $cMap->[$i];
        next if(!$members || !@$members);

        my @array = map {$hits->[$_]} @{$members};

	push(@clusters, \@array);
    }
    
    return \@clusters;
}
#-------------------------------------------------------------------------------
sub build_hsp_pairs {
    my $hits = shift;
    my $flank = shift;

    my @pairs;

    for (my $i = 0; $i < @$hits; $i++) {
	my $aSt = $hits->[$i]->strand('query');
	push(@pairs, [$i, $i]);

	for (my $j = $i+1; $j < @$hits; $j++) {
	    my $bSt = $hits->[$j]->strand('query');
	    next if($aSt ne $bSt);
	    my $code = compare::hsps_overlap($hits->[$i], $hits->[$j], 'query', $flank);
	    push(@pairs, [$i, $j]) if($code);
	}
    }

    return \@pairs;
}
#-------------------------------------------------------------------------------
sub cluster_pairs {
    my $pairs = shift;
    my $c_flag = shift;
    my $cId = 0;
    my @cId_index;
    my @cMap;
    my @cCount;

    foreach my $p (@$pairs){
        my ($mUidI, $mUidJ) = @{$p};
        if (!defined($cId_index[$mUidI]) && !defined($cId_index[$mUidJ])){
            $cId_index[$mUidI] = $cId;
            $cId_index[$mUidJ] = $cId;
	    if($mUidI == $mUidJ){
		push(@{$cMap[$cId]}, $mUidI);
	    }
	    else{
		push(@{$cMap[$cId]}, $mUidI, $mUidJ);
		$cCount[$cId]++;
	    }
            $cId++;
        }
        elsif (defined($cId_index[$mUidI]) && !defined($cId_index[$mUidJ])){
            my $cId = $cId_index[$mUidI];
            $cId_index[$mUidJ] = $cId;
            push(@{$cMap[$cId]}, $mUidJ);
	    $cCount[$cId]++;
        }
        elsif (!defined($cId_index[$mUidI]) && defined($cId_index[$mUidJ])){
            my $cId = $cId_index[$mUidJ];
            $cId_index[$mUidI] = $cId;
            push(@{$cMap[$cId]}, $mUidI);
	    $cCount[$cId]++;
        }
        elsif (defined($cId_index[$mUidI]) && defined($cId_index[$mUidJ])){
	    next if($mUidI == $mUidJ);

	    my $cIdI = $cId_index[$mUidI];
	    $cCount[$cIdI]++;
            next if ($cId_index[$mUidI] == $cId_index[$mUidJ]);

            my @copy;
            my $cIdJ = $cId_index[$mUidJ];
            foreach my $mUid (@{$cMap[$cIdJ]}){
                $cId_index[$mUid]= $cIdI;
                push(@{$cMap[$cIdI]}, $mUid);
            }

	    $cCount[$cIdI] += $cCount[$cIdJ];
	    undef $cCount[$cIdJ];
            undef $cMap[$cIdJ];
        }
    }

    return ($c_flag) ?  (\@cMap, \@cCount) : \@cMap;
}

1;


