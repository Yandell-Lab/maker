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
use compare;
use cluster;
use clean;
use Bio::Search::HSP::PhatHSP::blastx;
use Bio::Search::HSP::PhatHSP::blastn;
use Bio::Search::HSP::PhatHSP::tblastx;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub process {
	my @features;  
	while (my $f = shift){
	   push(@features, @{$f});
	}

	my ($tes, $lcs) = separate_types(\@features);

	my $shattered_lcs = shatter_hits($lcs);
	my $tes_keepers   = clean_tes($tes);

	my @best_keepers  = (@{$shattered_lcs}, @{$tes_keepers});

	return (\@best_keepers);
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
	my $tes = shift;

	my $shattered_hits = shatter_hits($tes);
	my $clusters = cluster::shadow_cluster(0, $shattered_hits, 20);

	my @keepers;
	foreach my $c (@{$clusters}){
		my @sorted = sort l_sort @{$c};
		push(@keepers, colapse_hit(\@sorted));
	}
	return \@keepers;
}
#------------------------------------------------------------------------
sub colapse_hit {
    my $hits = shift;
    my $best = $hits->[0];

    my $ref = ref($best);

    my $new_hit = new $ref('-name'         => $best->name,
			   '-description'  => $best->description,
			   '-algorithm'    => $best->algorithm,
			   '-length'       => $best->length,
	);
    
    $new_hit->queryLength($best->queryLength);
    $new_hit->database_name($best->database_name);
    $new_hit->{_label} = $best->{_label};

    my $B;
    my $E;
    foreach my $h (@$hits){
	$B = $h->start('query') if(!$B || $h->start('query') < $B);
	$E = $h->end('query') if(!$E || $h->end('query') > $E);
    }
    my $l = $E-$B+1;

    my @args;
    (my $hsp) = $best->hsps; #just need one hsp


    push(@args, '-query_start');
    push(@args, $B);

    push(@args, '-query_end');
    push(@args, $E);

    push(@args, '-score');
    push(@args, $hsp->{SCORE});

    #push(@args, '-query_seq');
    #push(@args, 'N'x$l);

    #push(@args, '-homology_seq');
    #push(@args, '|'x$l);

    #push(@args, '-hit_seq');
    #push(@args, 'N'x$l);

    push(@args, '-hit_start');
    push(@args, $hsp->start('hit'));

    push(@args, '-hsp_length');
    push(@args, $l);

    push(@args, '-identical');
    push(@args, $hsp->num_identical);

    push(@args, '-hit_length');
    push(@args, $l);

    push(@args, '-query_name');
    push(@args, $hsp->name());

    push(@args, '-algorithm');
    push(@args, $hsp->algorithm);

    push(@args, '-bits');
    push(@args, $hsp->bits);

    push(@args, '-evalue');
    push(@args, $hsp->evalue);

    push(@args, '-pvalue');
    push(@args, $hsp->signifcance);

    push(@args, '-query_length');
    push(@args, $hsp->query->length);

    push(@args, '-conserved');
    push(@args, $hsp->num_conserved);

    push(@args, '-hit_name');
    push(@args, $hsp->name());

    push(@args, '-hit_end');
    push(@args, $hsp->end('hit'));

    push(@args, '-query_gaps');
    push(@args, $hsp->gaps('query'));

    push(@args, '-hit_gaps');
    push(@args, $hsp->gaps('hit'));

    $ref = ref($hsp);
    my $new_hsp = new $ref(@args);
    
    $new_hsp->queryName($hsp->queryName) if defined($hsp->queryName);
    
    $new_hsp->{_strand_hack}->{query} = $hsp->strand('query');
    $new_hsp->{_strand_hack}->{hit}   = $hsp->strand('hit');
    $new_hsp->{_identical_hack}       = $hsp->frac_identical();
    $new_hsp->{_label}                = $hsp->{_label};
    
    $new_hit->add_hsp($new_hsp);
    $new_hit->{_HMM} = $best->{_HMM} if($best->{_HMM});
    $new_hit->{_label} = $best->{_label} if($best->{_label});

    return $new_hit;
}
#-----------------------------------------------------------------------------
sub mask_chunk {
	my $chunk = shift;
	my $features = shift;

	my ($tes, $lcs) = separate_types($features);

	my $chunk_offset = $chunk->offset();

	my $tes_coors = get_hsp_coors($tes, $chunk_offset);
	my $lcs_coors = get_hsp_coors($lcs, $chunk_offset);

	my $seq = $chunk->seq();

	_hard_mask_seq ($seq, $tes_coors, 50, 'N');
    	_soft_mask_seq($seq, $lcs_coors, 0);

	$chunk->seq($seq);
	
	return $chunk;
}
#-----------------------------------------------------------------------------
sub mask_seq {
	my $seq = shift;
	my $features = shift;

	my ($tes, $lcs) = separate_types($features);

	my $tes_coors = get_hsp_coors($tes, 0);
	my $lcs_coors = get_hsp_coors($lcs, 0);

	_hard_mask_seq (\$seq, $tes_coors, 50, 'N');
    	_soft_mask_seq(\$seq, $lcs_coors, 0);
	
	return $seq;
}
#-----------------------------------------------------------------------------
sub _hard_mask_seq {
   my $seq = shift;
   my $features = shift;
   my $flank = shift || 0;
   my $replace = shift || 'N';
   
   foreach my $p (@{$features}){
      my $b = $p->[0];
      my $e = $p->[1];
      
      ($b, $e) = ($e, $b) if $e < $b;
      
      my $first = $b - $flank;
      my $last = $e + $flank;
      $b = ($first > 0) ? $b - $flank : 1;
      $e = ($last <= length($$seq)) ? $e + $flank : length($$seq);
   
      my $l = $e - $b + 1;
      next if($l <= 0);

      #my $replace_string = substr($$seq, $b -1 , $l);
      #$replace_string =~ s/./$replace/g;
      
      #substr($$seq, $b -1 , $l, $replace_string);
      substr($$seq, $b -1 , $l, "$replace"x$l);
   }  
}
#-----------------------------------------------------------------------------
sub _soft_mask_seq {
   my $seq = shift;
   my $features = shift;
   my $flank = shift || 0;
   
   foreach my $p (@{$features}){
      my $b = $p->[0];
      my $e = $p->[1];
      
      ($b, $e) = ($e, $b) if $e < $b;

      my $first = $b - $flank;
      my $last = $e + $flank;
      $b = ($first > 0) ? $b - $flank : 1;
      $e = ($last <= length($$seq)) ? $e + $flank : length($$seq);      
   
      my $l = $e - $b + 1;
      next if($l <=0);
      
      substr($$seq, $b -1 , $l, lc(substr($$seq, $b -1 , $l)));
   }
}
#-----------------------------------------------------------------------------
sub get_hsp_coors {
	my $hits = shift;
	my $offset = shift || 0;

	my @coors;
	foreach my $hit (@{$hits}){
        	foreach my $hsp ($hit->hsps()){
			push(@coors, [$hsp->nB('query') - $offset, $hsp->nE('query') - $offset]);
		}
	}

	return \@coors;
}
#-----------------------------------------------------------------------------
sub separate_types {
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
		elsif (ref($f) =~ /gff/){
		    my $n = $f->hsp('best')->hit->seq_id();
		    my $s = $f->algorithm;
		    if ($n =~ /Simple_repeat/ || $n =~ /Low_complexity/){
			push(@lcs, $f);
		    }
		    elsif($s =~ /blastx/i){
			push(@tes, $f);
		    }
		    else {
			push(@tes, $f);
		    }
		}
		else{
		    push(@tes, $f);
		}
	}
	return (\@tes, \@lcs);
}
#-----------------------------------------------------------------------------
#sub gff {
#    my $chunk = shift @_;
#    my $gff_file = shift @_;

#     my %repeat_regions;
#     my $parser = new Bio::Tools::GFF(-gff_version => 3,
# 				     -file        => "$gff_file");
    
#     while (my $feature = $parser->next_feature()) {
# 	my $tag = $feature->primary_tag;
# 	my $id = $feature->seq_id;
# 	my $b = $feature->start;
# 	my $e = $feature->end;
# 	push @{$repeat_regions{$id}}, [$b, $e];
#     }
#     my @seq = split(//, $qseq);
#     foreach my $coor (@{$repeat_regions{$qname}}) {
# 	my ($b, $e) = @{$coor};
# 	foreach my $i ($b..$e) {
# 	    $seq[$i] = "N";
# 	}
#     }
#     my $mseq = join("", @seq);
#    return $chunk;
#}
#-----------------------------------------------------------------------------
1;


