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
use SimpleCluster;
use compare;
use cluster;
use exonerate::splice_info;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub throw_out_bad_splicers {
	my $phat_hits = shift;
        my @keepers;
        foreach my $hit (@{$phat_hits}){
		die " clean::throw_out_bad_slicers only works on est2genome hits!\n"
		unless ref($hit) =~ /est2genome$/;

                die " clean::throw_out_bad_slicers only works on spliced ests
		call clean::purge_single_exon_hits first!\n"
                if $hit->num_hsps == 1;

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
sub remove_redundant_alt_splices {
        my $candidates = shift;
        my $seq        = shift;
        my $flank      = shift;

        return [] unless defined($candidates);

        my @sorted =  sort alt_sort @{$candidates};

        my $num = @sorted;
        my %dead;
        for (my $i = 0; $i < @sorted -1 ;$i++){
                print STDERR " ...processing $i of $num\n";
                my $hit_i = $sorted[$i];
                if (defined($dead{$i})){
                        next;
                }
                for (my $j = $i+1; $j < @sorted; $j++){
                        my $hit_j = $sorted[$j];
                        if (defined($dead{$j})){
                                next;
                        }
                        if (compare::is_redundant_alt_form($hit_i, $hit_j, $seq, $flank)){
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
	my $seq        = shift;
	my $flank      = shift;

	return [] unless defined($candidates);

        my @sorted =  sort alt_sort @{$candidates};

	my $num = @sorted;
        my %dead;
        for (my $i = 0; $i < @sorted -1 ;$i++){
		print STDERR " ...processing $i of $num\n";
		my $hit_i = $sorted[$i];
		if (defined($dead{$i})){
			next;
		}
                for (my $j = $i+1; $j < @sorted; $j++){
			my $hit_j = $sorted[$j];
			if (defined($dead{$j})){
				next;
			}
			if (compare::is_same_alt_form($hit_i, $hit_j, $seq, $flank)){
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
1;


