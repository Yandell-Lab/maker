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
use clean;
use compare;
use SimpleCluster;
use Shadower;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub clean_and_cluster {
	my $keepers = shift;
	my $seq     = shift;
	my $depth   = shift;

	$depth = 0 if (! defined($depth) || $depth < 0);

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
		my $clean = clean::complexity_filter($c, $seq);
		my $alts = clean::get_best_alt_splices($clean, $seq, 10);
		my $i = 0;
		my @new_cluster;
		foreach my $a (@{$alts}){
			push(@new_cluster, $a);
			last if ($i > $depth && $depth > 0);
			$i++;
		}
		next if(! @new_cluster);
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
sub mani_sort {
	criteria($b) <=> criteria($a) || criteria_2($b) <=> criteria_2($a) || criteria_3($b) <=> criteria_3($a);
}
#------------------------------------------------------------------------
sub criteria {
	my $hit = shift;

	my $ref = $hit->algorithm;

	return  0 if $ref =~ /^repeat/i;
	return  0 if $ref =~ /^blastx\:repeatmask$/i;
	return  0 if $ref =~ /^repeat_gff\:/i;
	return  1 if $ref =~ /^snap$/i;
	return  1 if $ref =~ /^augustus$/i;
	return  1 if $ref =~ /^fgenesh$/i;
	return  1 if $ref =~ /^twinscan$/i;
	return  1 if $ref =~ /^jigsaw$/i;
	return  1 if $ref =~ /^pred_gff\:/i;
	return  2 if $ref =~ /^blastn$/i;
	return  2 if $ref =~ /^tblastx$/i;
	return  2 if $ref =~ /^altest_gff\:/i;
	return  3 if $ref =~ /^blastx$/i;
	return  3 if $ref =~ /^protein_gff\:/i;
	return  4 if $ref =~ /2genome$/i;
	return  4 if $ref =~ /^est_gff\:/i;
	return  5 if $ref =~ /^maker$/i;
	return  5 if $ref =~ /^model_gff\:/i;
	die "UNKNOWN CLASS(".ref($hit)."), ALGORITHM($ref), in cluster::criteria\n";
}
#------------------------------------------------------------------------
sub criteria_2 {
        my $hit = shift;

	my $score =  $hit->score();
	$score = '' if(! defined $score);
        #check if value is numerical ad adjust if not
	if($score !~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/){
	    $score = $hit->bits();
	    $score = '' if(! defined $score);
	}
	if($score !~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/){
	    $score = $hit->significance;
	    $score = '' if(! defined $score);
	    if($score =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/){
		$score = 10000 - (10000 * $score);
		$score = 0 if($score < 0);
	    }
	    else{
		$score = 0;
	    }
	}

	return $score;
}
#------------------------------------------------------------------------
#third citeria used to address order issue with mpi vs standard maker
sub criteria_3 {
    my $hit = shift;

    my ($lAq) = $hit->getLengths();
    return $lAq;
}
#------------------------------------------------------------------------
sub shadow_cluster {
        my $depth     = shift;
        my $seq       = shift;
        my $phat_hits = shift;
        my $flank     = shift;

	$depth = 0 if (! defined($depth) || $depth < 0);

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

	my @sorted = ($depth == 0) ? @{$phat_hits} : sort mani_sort @{$phat_hits};
 
	print STDERR "... finished.\n" unless $main::quiet;

	foreach my $hit (@sorted){

		my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit, 'query');
                   ($nB, $nE) = ($nE, $nB) if $nB > $nE;

		my $j = 0;
        	foreach my $s (@{$pieces}){

		        if ($depth != 0 &&
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
				last;
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
#for clustering multiple external_gff::PhatHit objects using the gene_id
sub gene_cluster {
   my $gff_hits = shift;

   my %index;
   my @clusters;

   foreach my $trans (@{$gff_hits}){
      my $gene_id = $trans->gene_id();

      #external_gff::PhatHit hits must have a gene_id
      die "ERROR: No gene id in hit\n" if (! defined $gene_id);

      if(! exists $index{$gene_id}){
	 push(@clusters, [$trans]);
	 $index{$gene_id} = @clusters - 1;
      }
      else{
	 my $c_id = $index{$gene_id};

	 push(@{$clusters[$c_id]}, $trans);
      }
   }

   return \@clusters;
}
#------------------------------------------------------------------------
sub get_overlap_evidence {
	my $p_hit     = shift;
        my $phat_hits = shift;
        my $depth     = shift;
        my $flank     = shift;
	
	$flank = 0 if (!$flank || $flank < 0);
	$depth = 0 if (!$depth || $depth < 0);

	my ($sB, $sE) = PhatHit_utils::get_span_of_hit($p_hit, 'query');
	($sB, $sE) = ($sE, $sB) if $sB > $sE;
	
	$sB = $sB - $flank;
	$sE = $sE + $flank;    
	
	my @cluster;
	push(@cluster, $p_hit); 
	
	my $i = 0;
	foreach my $hit (@{$phat_hits}){	      
	   if ($depth > 0 && $i >= $depth){
	      last;
	   }

	   my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit, 'query');
	   ($nB, $nE) = ($nE, $nB) if $nB > $nE;
	   
	   my $class = compare::compare($sB, $sE, $nB, $nE);
	   
	   if ($class ne '0'){
	      push(@cluster, $hit);
	      $i++;
	   }
	}

	return \@cluster;
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
	    if (compare::hsps_overlap($hit_i, $hit_j, 'query', $flank)){
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
	next if(! defined $map->{$c});
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


