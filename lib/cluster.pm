#------------------------------------------------------------------------
#----                           cluster                              ---- 
#------------------------------------------------------------------------
package cluster;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use SimpleCluster;
use clean;
use compare;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub clean_and_cluster {
	my $keepers = shift;
	my $seq     = shift;
	my $depth   = shift;

	my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $keepers);

	my $p_clusters = shadow_cluster($depth, $seq, $p);
        my $m_clusters = shadow_cluster($depth, $seq, $m);

	my @clusters = (@{$p_clusters}, @{$m_clusters});

	my $num_c = @clusters;
	print STDERR "cleaning clusters....\n" unless $main::quiet;
	my $counter = 0;
	my @clean_clusters;

	#show_clusters(\@clusters);
	#die;

	foreach my $c (@clusters){
		print STDERR "total clusters:$num_c now processing $counter\n"
			unless $main::quiet;
		my $alts = clean::get_best_alt_splices($c, $seq, 10);
		my $i = 0;
		my @new_cluster;
		foreach my $a (@{$alts}){
			push(@new_cluster, $a);
			last if $i > $depth;
			$i++;
		}
		push(@{$clean_clusters[$counter]}, @new_cluster);
		$counter++;
	}
		

	return \@clean_clusters;

}
#------------------------------------------------------------------------
sub special_cluster_phat_hits {
        my $phat_hits = shift;        
	my $ests      = shift;
	my $seq       = shift;        
	my $flank     = shift || 10;

        my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $phat_hits);

        my $p_clusters = shadow_cluster(20, $seq, $p, $flank);
        my $m_clusters = shadow_cluster(20, $seq, $m, $flank);

        my @careful_clusters;
        foreach my $c (@{$p_clusters}){
                my $cares = careful_cluster($seq, $c, $flank);
                push(@careful_clusters, @{$cares});
        }

        foreach my $c (@{$m_clusters}){
                my $cares = careful_cluster($seq, $c, $flank);
                push(@careful_clusters, @{$cares});
        }

	my $temp_id = 0;
	foreach my $est (@{$ests}){
		$est->{temp_id} = $temp_id;
		$temp_id++;
	}
	my @special_clusters;
	my %done;
	my $i = 0;
	foreach my $c (@careful_clusters){
		my $coors  = PhatHit_utils::to_begin_and_end_coors($c, 'query');
		my $pieces = Shadower::getPieces($seq, $coors, $flank);

		push(@{$special_clusters[$i]}, @{$c});
		die "multiple pieces in cluster::special_cluster_phat_hits!\n"
		if defined($pieces->[1]);

		foreach my $est (@{$ests}){
			my ($eB, $eE) = PhatHit_utils::get_span_of_hit($est, 'query');
			($eB, $eE) = ($eE, $eB) if $eB > $eE; 
			my $s_code = compare::get_shadow_code($pieces, 
			                                      [[$eB, $eE]], 
			                                      $flank,
			                                      );


			if ($s_code ne '0'){
				push(@{$special_clusters[$i]}, $est);
				$done{$est->{temp_id}}++;
			}
		}
		$i++;
	}

	my @remains;
	foreach my $est (@{$ests}){
		push(@remains, $est) 
		unless defined($done{$est->{temp_id}});
	}
	my $new_e_clusters = shadow_cluster(20, $seq, \@remains, $flank);

	my $start = @special_clusters;

	foreach my $c (@{$new_e_clusters}){
		push(@{$special_clusters[$start]}, @{$c});
		$start++;
	}

        return \@special_clusters;
	
}
#------------------------------------------------------------------------
sub careful_cluster_phat_hits {
        my $phat_hits = shift;
	my $seq       = shift;
	my $flank     = shift || 10;
        my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $phat_hits);

        my $p_clusters = shadow_cluster(20, $seq, $p, $flank);
        my $m_clusters = shadow_cluster(20, $seq, $m, $flank);
        
        my @careful_clusters;
        foreach my $c (@{$p_clusters}){
                my $cares = careful_cluster($seq, $c, $flank);
                push(@careful_clusters, @{$cares});
        }

        foreach my $c (@{$m_clusters}){
                my $cares = careful_cluster($seq, $c, $flank);
                push(@careful_clusters, @{$cares});
        }

        return \@careful_clusters;
}
#------------------------------------------------------------------------
sub careful_cluster {
        my $seq       = shift;
        my $phat_hits = shift;
        my $flank     = shift;


	my @careful_clusters;
	if (!defined($phat_hits->[1])){
		push(@careful_clusters, $phat_hits);
		return \@careful_clusters;
	}
	my $temp_id = 0;
	foreach my $hit (@{$phat_hits}){
		$hit->{temp_id} = $temp_id;
		
		$temp_id++;
	}
	print STDERR "now careful_clustering....\n"
		unless $main::quiet;
	my %lookup;
	my %matrix;
	for (my $i = 0; $i < @{$phat_hits} - 1;$i++){
		my $hit_i = $phat_hits->[$i];
		$lookup{$hit_i->{temp_id}} = $hit_i;
		for (my $j = $i +1 ; $j < @{$phat_hits}; $j++){
			my $hit_j = $phat_hits->[$j];
			$lookup{$hit_j->{temp_id}} = $hit_j;
			if (compare::overlap($hit_i, $hit_j, 'query', $flank)){
				$matrix{$hit_i->{temp_id}}{$hit_j->{temp_id}}++;
			}
			else {
				#
			}

		}
	}

	my $pairs = SimpleCluster::pairs(\%matrix);
	my $map   = SimpleCluster::singleLinkageClusters($pairs);

	my %clustered;
	my $i = 0;
	foreach my $c (keys %{$map}){
		foreach my $m (@{$map->{$c}}){
			my $hit = $lookup{$m};
			die "name not found in careful cluster!\n"
			unless defined $m;
			push(@{$careful_clusters[$i]}, $hit);
			$clustered{$m}++;
		}
		$i++;
	}
	# get the left overs...
	foreach my $temp_id (keys %lookup){
		next if defined($clustered{$temp_id});
		push(@{$careful_clusters[$i]}, $lookup{$temp_id});
		$i++;
	}
	
	return \@careful_clusters;
}
#------------------------------------------------------------------------
sub mani_sort {
	criteria($b) <=> criteria($a) || criteria_2($b) <=> criteria_2($a) || criteria_3($a) <=> criteria_3($b);
}
#------------------------------------------------------------------------
sub criteria {
	my $hit = shift;

	my $ref = ref($hit);

	return  2 if $ref =~ /blast/;
	return  3 if $ref =~ /2genome/;
	return  1 if $ref =~ /snap/;
	return  1 if $ref =~ /augustus/;
	return  0;
}
#------------------------------------------------------------------------
sub criteria_2 {
        my $hit = shift;

        my $ref = ref($hit);

        return  $hit->hsp('best')->score  if $ref =~ /blast/;
	return  $hit->hsp('best')->score  if $ref =~ /repeatmasker/;
	return  $hit->hsp('best')->bits   if $ref =~ /2genome/;
	return $hit->score()              if $ref =~/snap/;	
	return $hit->score()              if $ref =~/augustus/;
	die "UNKNOWN CLASS(".ref($hit).") in cluster::criteria_2\n";
}
#------------------------------------------------------------------------
#third citeria used to address order issue with mpi vs standard maker
sub criteria_3 {
    my $hit = shift;

    my $ref = ref($hit);

    return  $hit->hsp('best')->evalue if $ref =~ /blast/;
    return  1                         if $ref =~ /repeatmasker/;
    return  1                         if $ref =~ /2genome/;
    return  1                         if $ref =~/snap/;
    return  1                         if $ref =~/augustus/;
    die "UNKNOWN CLASS(".ref($hit).") in cluster::criteria_3\n";
}
#------------------------------------------------------------------------
sub shadow_cluster {
        my $depth     = shift;
        my $seq       = shift;
        my $phat_hits = shift;
        my $flank     = shift;

        my $coors  = PhatHit_utils::to_begin_and_end_coors($phat_hits, 'query');
        my $pieces = Shadower::getPieces($seq, $coors, $flank);
 
        my $temp_id = 0;
        foreach my $hit (@{$phat_hits}){
                $hit->{temp_id} = $temp_id;
                $temp_id++;
        }

        my @clusters;
        my $i_size = @{$phat_hits};
        my $j_size = @{$pieces};

        print STDERR " in cluster:shadow cluster...\n" unless $main::quiet;
        print STDERR "    i_size:$i_size j_size:$j_size\n"
		unless $main::quiet;

	my $i = 0;
	my %c_size;

	print STDERR " sorting hits in shadow cluster...\n"
		unless $main::quiet;

	my @sorted = sort mani_sort @{$phat_hits};
 
	print STDERR "... finished.\n" unless $main::quiet;

	foreach my $hit (@sorted){

		my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit, 'query');
                   ($nB, $nE) = ($nE, $nB) if $nB > $nE;

		my $j = 0;
        	foreach my $s (@{$pieces}){

		        if (defined($depth) &&
			    defined($c_size{$j}) &&
			    $c_size{$j} >  $depth
			   ){
			      $j++;
			      next;
			}

                	my $sB = $s->{b};
                	my $sE = $s->{e};

                        my $class = compare::compare($sB, $sE, $nB, $nE);

                        if ($class ne '0'){
                                push(@{$clusters[$j]}, $hit); 

				$c_size{$j}++;
                        }
			
			$j++;
		}

		print STDERR " i_size:$i_size   current i:$i\n" unless $main::quiet;
		$i++;
	}

	return \@clusters;
}
#------------------------------------------------------------------------
sub shadow_cluster_old {
	my $depth     = shift;
	my $seq       = shift;
	my $phat_hits = shift;
	my $flank     = shift;

	my $coors  = PhatHit_utils::to_begin_and_end_coors($phat_hits, 'query');
	my $pieces = Shadower::getPieces($seq, $coors, $flank);

	my $temp_id = 0;
	foreach my $hit (@{$phat_hits}){
		$hit->{temp_id} = $temp_id;
		$temp_id++;
	}
	my @clusters;
	my $i = 0;
	my $i_size = @{$pieces};
	my $j_size = @{$phat_hits};
	print STDERR " in cluster:shadow cluster...\n" unless $main::quiet;
	print STDERR "    i_size:$i_size j_size:$j_size\n" unless $main::quiet;
	my %clustered;
        foreach my $s (@{$pieces}){
        	my $sB = $s->{b};
                my $sE = $s->{e};

		my $n = 0;
                foreach my $hit (@{$phat_hits}){

			next if defined($clustered{$hit->{temp_id}});

			print STDERR "     i:$i   i_size:$i_size j_size:$j_size\n"
				unless $main::quiet;
			print STDERR "          n:$n\n" unless $main::quiet;
		        my ($nB, $nE) = 
			PhatHit_utils::get_span_of_hit($hit, 'query');

			($nB, $nE) = ($nE, $nB) if $nB > $nE;

                 	my $class = compare::compare($sB, $sE, $nB, $nE);

                        if ($class ne '0'){
				#print "i:$i n:$n\n";
                                push(@{$clusters[$i]}, $hit) if $n < $depth;
                                $clustered{$hit->{temp_id}}++;
                                $n++;
                        }

                }
		$i++;
        }

	#show_clusters(\@clusters);
	return \@clusters;
}
#------------------------------------------------------------------------
sub show_clusters {
	my $c = shift;

	for (my $i = 0; $i < @{$c}; $i++){
		print "cluster:$i\n";
			foreach my $m (@{$c->[$i]}){
				my $num_hsps = $m->hsps();

				print "   name:".$m->name;
				print " type:".ref($m);
				print " num hsps:$num_hsps";
				print " length:".$m->length();
				print " nB:".$m->nB('query')."\n";
			}
	}
}
#------------------------------------------------------------------------
1;


