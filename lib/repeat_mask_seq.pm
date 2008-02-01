#------------------------------------------------------------------------
#----                         repeat_mask_seq                        ---- 
#------------------------------------------------------------------------
package repeat_mask_seq;
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
use clean;
@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub process {
	my $query_seq             = shift;
	my $rm_keepers            = shift;
	my $repeat_blastx_keepers = shift;

	my @features = (@{$rm_keepers}, @{$repeat_blastx_keepers});

	my ($tes, $lcs) = seperate_types(\@features);

	my $masked_seq = mask_seq($query_seq, $tes, $lcs);

	my $shattered_lcs = shatter_hits($lcs);
	my $tes_keepers   = clean_tes($query_seq, $tes);

	my @best_keepers  = (@{$shattered_lcs}, @{$tes_keepers});

	return ($masked_seq, \@best_keepers);

}
#-----------------------------------------------------------------------------
sub shatter_hits {
	my $hits = shift;
        my @shattered_hits;
        foreach my $hit (@{$hits}){
                my $shatter_hits = PhatHit_utils::shatter_hit($hit);
                push(@shattered_hits, @{$shatter_hits});
        }

	return \@shattered_hits;
}
#------------------------------------------------------------------------
sub span_length {
	my $hit = shift;

	my ($b, $e) = PhatHit_utils::get_span_of_hit($hit, 'query');

	return abs($b - $e);

}
#------------------------------------------------------------------------
sub l_sort {
        span_length($b) <=> span_length($a)
}
#------------------------------------------------------------------------
sub clean_tes {
	my $seq = shift;
	my $tes = shift;

	my $shattered_hits = shatter_hits($tes);

	my $clusters = cluster::shadow_cluster(10, $seq, $shattered_hits, 20);

	my @keepers;
	foreach my $c (@{$clusters}){
		my @sorted = sort l_sort @{$c};
		push(@keepers, $sorted[0]);
	}
	return \@keepers;
}
#-----------------------------------------------------------------------------
sub mask_seq {
	my $seq = shift;
	my $tes = shift;
	my $lcs = shift;

	my $tes_coors = get_coors($tes);
	my $lcs_coors = get_coors($lcs);

	my $masked_seq = Shadower::maskSequence($seq, $tes_coors, 50, 'N');
    	   $masked_seq = Shadower::softMaskSequence($masked_seq,
                                                    $lcs_coors,
                                                    0,
                                                   );
	return $masked_seq;
}
#-----------------------------------------------------------------------------
sub get_coors {
	my $hits = shift;

	my @coors;
	foreach my $hit (@{$hits}){
        	foreach my $hsp ($hit->hsps()){
			push(@coors, [$hsp->nB('query'), $hsp->nE('query')]);
		}
	}

	return \@coors;
}
#-----------------------------------------------------------------------------
sub seperate_types {
	my $features = shift;

	my @tes;
	my @lcs;
	foreach my $f (@{$features}){
		if (ref($f) =~ /repeatmasker/){
			my $n = $f->hsp('best')->hit->seq_id();
			if ($n =~ /Simple_repeat/ || $n =~ /Low_complexity/){
				push(@lcs, $f);
			}
			else {
				push(@tes, $f);
			}
		}
		elsif (ref($f) =~ /blastx/) {
			push(@tes, $f);
		}
	}
	return (\@tes, \@lcs);
}
#-----------------------------------------------------------------------------
sub gff {
    my $qseq = shift @_;
    my $qname = shift @_;
    my $gff_file = shift @_;
    my %repeat_regions;
    my $parser = new Bio::Tools::GFF(-gff_version => 3,
				     -file        => "$gff_file");
    
    while (my $feature = $parser->next_feature()) {
	my $tag = $feature->primary_tag;
	my $chromosome = $feature->seq_id;
	my $start = $feature->start;
	my $end = $feature->end;
	my $offset = int($start/1e6)*1e6;
	my $chend = $offset + 1e6;
	my $sname = join(":",$chromosome, ($offset+1));
	my $seqname = join("..", $sname, $chend);
	my $b = $start - $offset - 1;
	my $e = $end - $offset - 1;
	push @{$repeat_regions{$seqname}}, [$b, $e];
    }
    my @seq = split(//, $qseq);
    foreach my $coor (@{$repeat_regions{$qname}}) {
	my ($b, $e) = @{$coor};
	foreach my $i ($b..$e) {
	    $seq[$i] = "N";
	}
    }
    my $mseq = join("", @seq);
    return $mseq;
}
#-----------------------------------------------------------------------------
1;


