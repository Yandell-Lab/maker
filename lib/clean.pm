#------------------------------------------------------------------------
#----                                clean                           ---- 
#------------------------------------------------------------------------
package clean;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use exonerate::splice_info;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub purge_short_single_exons {
    my $phat_hits = shift;
    my $min = shift;

    my @keepers;
    foreach my $hit (@{$phat_hits}){
	if ($hit->num_hsps == 1){
	    push(@keepers, $hit) if($hit->{_hsps}->[0]->length('query') >= $min);
	}
	else{
	    push(@keepers, $hit);
	}
    }
    
    return \@keepers
}
#------------------------------------------------------------------------
sub throw_out_bad_splicers {
	my $phat_hits = shift;
	my $seq = shift;

        my @keepers;
        foreach my $hit (@{$phat_hits}){
                if ($hit->num_hsps == 1){
		    push(@keepers, $hit);
		    next;
		}

		exonerate::splice_info::set_donors_acceptors($hit, $seq);
		my $str = exonerate::splice_info::get_splice_str($hit);

		my $num = exonerate::splice_info::get_numbers($str);

		my $s_ratio = $num->{p}/$num->{l};
		next if $s_ratio < 0.75;

                push(@keepers, $hit);
        }
        return \@keepers;
}
#------------------------------------------------------------------------
sub purge_single_exon_hits {
	my $phat_hits = shift;

	my @keepers;
	foreach my $hit (@{$phat_hits}){
		my $num_hsps = $hit->num_hsps();
		next unless $num_hsps > 1;
		push(@keepers, $hit);
	}
	return \@keepers;
}
#------------------------------------------------------------------------
sub alt_sort {
        $b->num_hsps <=> $a->num_hsps || $b->length  <=>  $a->length;
}
#------------------------------------------------------------------------
sub is_identical_form {
    my $hit1 = shift;
    my $hit2 = shift;

    my @hsps1 = sort {$a->start('query') <=> $b->start('query')} $hit1->hsps;
    my @hsps2 = sort {$a->start('query') <=> $b->start('query')} $hit2->hsps;

    return if (@hsps1 != @hsps2);

    for(my $i = 0; $i < @hsps1; $i++){
	return if($hsps1[$i]->start('query') != $hsps2[$i]->start('query'));
	return if($hsps1[$i]->end('query') != $hsps2[$i]->end('query'));
    }

    return 1;
}
#------------------------------------------------------------------------
sub remove_redundant_alt_splices {
        my $candidates = shift;
        my $flank      = shift;
	my $depth      = shift;

        return [] unless defined($candidates);
	return $candidates if(@{$candidates} == 1);

        my @sorted =  sort alt_sort @{$candidates};

        my $num = @sorted;
        my %dead;
        my @keepers;
	for(my $i = 0; $i < @sorted; $i++){
	    print STDERR " ...processing $i of $num\n" unless($main::quiet);
	    
	    my $hit_i = $sorted[$i];
	    for(my $j = 0; $j < @keepers; $j++){
		my $hit_j = $keepers[$j];
		
		if (compare::is_redundant_alt_form($hit_i, $hit_j, $flank)){
		    $dead{$i}++;
		    last;
		}
	    }
	    next if ($dead{$i});
	    
	    push(@keepers, $hit_i);
	    if($depth && @keepers == $depth){
		print STDERR " ...trimming the rest\n" unless($main::quiet);
		last;
	    }
        }

        return \@keepers;
}
#------------------------------------------------------------------------
sub get_best_alt_splices {
        my $candidates = shift;
	my $flank      = shift;

	return [] unless defined($candidates);

        my @sorted =  sort alt_sort @{$candidates};

	my $num = @sorted;
        my %dead;
        my @transcripts;
        for(my $i = 0; $i < @sorted;$i++){
	    print STDERR " ...processing $i of $num\n" unless($main::quiet);
	    my $hit_i = $sorted[$i];
	    
	    for(my $j = 0; $j < @transcripts; $j++){
		my $hit_j = $transcripts[$j];
		
		if (compare::is_same_alt_form($hit_i, $hit_j, $flank)){
		    $dead{$i}++;
		    last;
		}
	    }
	    next if ($dead{$i});

	    push(@transcripts, $sorted[$i]);
        }

        return \@transcripts;
}
#------------------------------------------------------------------------
#this function throws out hits that result from multiple overlapping HSPs
#from the same region (caused by low complexity repeats) 
sub complexity_filter {
    my $candidates = shift;
    
    my @keepers;
    
    foreach my $f (@$candidates){
	#handle cluster of clusters
	if(ref $f eq 'ARRAY'){
	    my $set = complexity_filter($f);
	    push(@keepers, $f) if(@$set);
	    next;
	}
	
	my @q_space; #space based coordinates
	my @h_space; #space based coordinates
	my @q_set; #hsp regions on query
	my @h_set; #hsp regions pn hit
	foreach my $hsp ($f->hsps){
	    my $qS = $hsp->start('query');
	    my $qE = $hsp->end('query');
	    
	    my $hS = $hsp->start('hit');
	    my $hE = $hsp->end('hit');

	    push(@q_space, ($qS-1, $qE));
	    push(@h_space, ($hS-1, $hE));
	    push(@q_set, [$qS, $qE]);
	    push(@h_set, [$hS, $hE]);
	    
	    #@q_cov[$qS..$qE] = map {(defined $_) ? $_ + 1 : 1} (@q_cov[$qS..$qE]);
	    #@h_cov[$hS..$hE] = map {(defined $_) ? $_ + 1 : 1} (@h_cov[$hS..$hE]);
	}

	@q_set = sort {$a->[0] <=> $b->[0]} @q_set;
	@h_set = sort {$a->[0] <=> $b->[0]} @h_set;

	my %uniq;
	@q_space = grep {!$uniq{$_}++} sort {$a <=> $b} @q_space; 
	undef %uniq;
	@h_space = grep {!$uniq{$_}++} sort {$a <=> $b} @h_space;

	#now count coverage on query regions
	my $count = 0;
	my $total = 0;
	for(my $i = 0; $i < @q_space - 1; $i++){
	    my $j = $i+1;

	    my $sB = $q_space[$i]; #these are space based
	    my $sE = $q_space[$j]; #these are space based

	    my $l = $sE - $sB;
	    $count += $l;
	    foreach my $h (@q_set){
		last if($h->[0] > $sE);
		my $B = $h->[0]; #these are coordinate based
		my $E = $h->[1]; #these are coordinate based
		
		$total += $l  if($B-1 <= $sB && $sE <= $E); #should be fully contained
	    }
	}
	my $q_ave = ($count) ? $total/$count : 0; #average coverage for hit

	#now count coverage on hit regions
	$count = 0;
	$total = 0;
	for(my $i = 0; $i < @h_space - 1; $i++){
	    my $j = $i+1;

	    my $sB = $h_space[$i]; #these are space based
	    my $sE = $h_space[$j]; #these are space based

	    my $l = $sE - $sB;
	    $count += $l;
	    foreach my $h (@h_set){
		last if($h->[0] > $sE);
		my $B = $h_set[$i]; #these are coordinate based
		my $E = $h_set[$j]; #these are coordinate based
		
		$total += $l  if($B-1 <= $sB && $sE <= $E); #should be fully contained
	    }
	}
	my $h_ave = ($count) ? $total/$count : 0; #average coverage for hit

	#if average basepair hits 1.5 times on both query and subject, skip
	push(@keepers, $f) unless(($q_ave + $h_ave)/2 > 1.5);
    }
    
    return \@keepers;
}
#------------------------------------------------------------------------
1;


