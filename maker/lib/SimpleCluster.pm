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

    my @pairs;

    for (my $i = 0; $i < @$hits; $i++) {
	my $aB = $hits->[$i]->nB('query');
	my $aE = $hits->[$i]->nE('query');
	my $aSt = $hits->[$i]->strand('query');
	push(@pairs, [$i, $i]);

	for (my $j = $i+1; $j < @$hits; $j++) {
	    my $bB = $hits->[$j]->nB('query');
	    my $bE = $hits->[$j]->nE('query');
	    my $bSt = $hits->[$j]->strand('query');

	    next if($aSt ne $bSt);

	    my $code = (abs($aB - $bE) + 1 <= $flank ||
			abs($aB - $bB) + 1 <= $flank ||
			abs($aE - $bE) + 1 <= $flank ||
			abs($aE - $bB) + 1 <= $flank);

	    $code = compare::compare($aB, $aE, $bB, $bE, $flank) if(!$code);
	    push(@pairs, [$i, $j]) if($code);
	}
    }

    return \@pairs;
}
#-------------------------------------------------------------------------------
sub build_pairs {
    my $array = shift;
    my $flank = shift;

    my @pairs;
    for (my $i = 0; $i < @$array; $i++) {
	my $aB = $array->[$i]->{start};
	my $aE = $array->[$i]->{end};
	my $aSt = $array->[$i]->{strand};
	push(@pairs, [$i, $i]);

	for (my $j = $i+1; $j < @$array; $j++) {
	    my $bB = $array->[$j]->{start};
	    my $bE = $array->[$j]->{end};
	    my $bSt = $array->[$j]->{strand};

	    next if(defined($aSt) && defined ($bSt) && ($aSt ne $bSt));

	    my $code = (abs($aB - $bE) + 1 <= $flank ||
			abs($aB - $bB) + 1 <= $flank ||
			abs($aE - $bE) + 1 <= $flank ||
			abs($aE - $bB) + 1 <= $flank);

	    $code = compare::compare($aB, $aE, $bB, $bE, $flank) if(!$code);
	    push(@pairs, [$i, $j]) if($code);
	}
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
    my $cId = 0;
    my @cId_index;
    my @cMap;

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
	    }
            $cId++;
        }
        elsif (defined($cId_index[$mUidI]) && !defined($cId_index[$mUidJ])){
            my $cId = $cId_index[$mUidI];
            $cId_index[$mUidJ] = $cId;
            push(@{$cMap[$cId]}, $mUidJ);
        }
        elsif (!defined($cId_index[$mUidI]) && defined($cId_index[$mUidJ])){
            my $cId = $cId_index[$mUidJ];
            $cId_index[$mUidI] = $cId;
            push(@{$cMap[$cId]}, $mUidI);
        }
        elsif (defined($cId_index[$mUidI]) && defined($cId_index[$mUidJ])){
            next if ($cId_index[$mUidI] == $cId_index[$mUidJ]);

            my @copy;
            my $cIdI = $cId_index[$mUidI];
            my $cIdJ = $cId_index[$mUidJ];
            foreach my $mUid (@{$cMap[$cIdJ]}){
                $cId_index[$mUid]= $cIdI;
                push(@{$cMap[$cIdI]}, $mUid);
            }

            undef $cMap[$cIdJ];
        }
    }

    return \@cMap;
}

1;


