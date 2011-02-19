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

        return [] unless defined($candidates);

        my @sorted =  sort alt_sort @{$candidates};

        my $num = @sorted;
        my %dead;
        for (my $i = 0; $i < @sorted -1 ;$i++){
                print STDERR " ...processing $i of $num\n" unless($main::quiet);
                my $hit_i = $sorted[$i];
                if (defined($dead{$i})){
                        next;
                }
                for (my $j = $i+1; $j < @sorted; $j++){
                        my $hit_j = $sorted[$j];
                        if (defined($dead{$j})){
                                next;
                        }
                        if (compare::is_redundant_alt_form($hit_i, $hit_j, $flank)){
                                $dead{$j}++;
                        }
                }
        }

        my @transcripts;
        for (my $i = 0; $i < @sorted; $i++){
                push(@transcripts, $sorted[$i])
                unless defined($dead{$i});
        }
        return \@transcripts;
}
#------------------------------------------------------------------------
sub get_best_alt_splices {
        my $candidates = shift;
	my $flank      = shift;

	return [] unless defined($candidates);

        my @sorted =  sort alt_sort @{$candidates};

	my $num = @sorted;
        my %dead;
        for (my $i = 0; $i < @sorted -1 ;$i++){
		print STDERR " ...processing $i of $num\n" unless($main::quiet);
		my $hit_i = $sorted[$i];
		if (defined($dead{$i})){
			next;
		}
                for (my $j = $i+1; $j < @sorted; $j++){
			my $hit_j = $sorted[$j];
			if (defined($dead{$j})){
				next;
			}
			if (compare::is_same_alt_form($hit_i, $hit_j, $flank)){
				$dead{$j}++;

			}
                }
        }

        my @transcripts;
        for (my $i = 0; $i < @sorted; $i++){
                push(@transcripts, $sorted[$i])
                unless defined($dead{$i});
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
	
	my $qOff = $f->start('query');
	my $hOff = $f->start('hit');
	
	my @q_cov;
	my @h_cov;
	foreach my $hsp ($f->hsps){
	    my $qS = $hsp->start('query') - $qOff;
	    my $qE = $hsp->end('query') - $qOff;
	    
	    my $hS = $hsp->start('hit') - $hOff;
	    my $hE = $hsp->end('hit') - $hOff;
	    
	    @q_cov[$qS..$qE] = map {(defined $_) ? $_ + 1 : 1} (@q_cov[$qS..$qE]);
		@h_cov[$hS..$hE] = map {(defined $_) ? $_ + 1 : 1} (@h_cov[$hS..$hE]);
	}
	
	#average coverage for qeury
	my $count = 0;
	my $total = 0;
	foreach my $c (@q_cov){
	    next if(! defined $c);
	    $count++;
	    $total += $c;		
	}
	
	my $q_ave = ($count) ? $total/$count : 0;
	
	#average coverage for hit
	$count = 0;
	$total = 0;
	foreach my $c (@h_cov){
	    next if(! defined $c);
	    $count++;
	    $total += $c;		
	}

	my $h_ave = ($count) ? $total/$count : 0;

	#if average basepair hits twice on both query and subject, skip
	push(@keepers, $f) unless(($q_ave + $h_ave)/2 > 1.5);
    }
    
    return \@keepers;
}
#------------------------------------------------------------------------
1;


