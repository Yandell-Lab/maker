#------------------------------------------------------------------------
#----                          maker::auto_annotator                 ----
#------------------------------------------------------------------------
package maker::auto_annotator;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use FastaChunker;
use Widget::snap;
use PhatHit_utils;
use CGL::TranslationMachine;
use compare;
use cluster;
use clean;
use maker::join;
use maker::quality_index;
use Widget::snap;
use Widget::augustus;
use GFFDB;
use shadow_AED;

$Storable::forgive_me = 1;

@ISA = qw(
       );

my $LOG; #GLOBAL VARIABLE
my $SEEN;

#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
#prepares hits by combining evidence into appropriate clusters
#returns clusters for abinits, gff_models, and standard evidence
sub prep_hits {
	my $prot_hits        = shift;
	my $est_hits         = shift;
	my $altest_hits      = shift;
	my $predictions      = shift;
	my $models           = shift;
	my $seq              = shift;
	my $single_exon      = shift;
	my $single_length    = shift;
	my $pred_flank       = shift;
	my $organism_type    = shift;
	my $est_forward      = shift;

	#---build clusters for basic evidence

	# combine puts type in order they are given
	# important for later culstering, as hits in
	# first arg more likey to make in into to cluster
	# than second arg, etc
	my $c_bag = combine($prot_hits,
			    $est_hits,
			    $altest_hits
			   );

	my ($p, $m, $x, $z) = PhatHit_utils::separate_by_strand('query', $c_bag);
	my $p_clusters = cluster::shadow_cluster(100, $seq, $p, $pred_flank, 1);
	my $m_clusters = cluster::shadow_cluster(100, $seq, $m, $pred_flank, 1);
	
	#this method will cause clusters that are near each other and are connected by an orf to merge.
	#this solves issues with mRNAseq splice site crossing reads and other EST partial exon coverage
	$p_clusters = join_clusters_around_orf($p_clusters, $seq);
	$m_clusters = join_clusters_around_orf($m_clusters, $seq);

	my $careful_clusters = [];
	push(@{$careful_clusters}, @{$p_clusters}, @{$m_clusters});

	#===purge ESTs after clustering so as to still have the effect of evidence joining ESTs
	# don't use unpliced single exon ESTs-- may be genomic contamination
	if($single_exon != 1 && $organism_type eq 'eukaryotic' && !$est_forward) {
	    $careful_clusters = purge_single_exon_hits_in_cluster($careful_clusters);
	}
	# throw out the exonerate est hits with weird splice sites
	if(!$est_forward){
	    $careful_clusters = throw_out_bad_splicers_in_cluster($careful_clusters, $seq);
	}
	#throw out short ESTs
	if($single_exon == 1 || $organism_type eq 'prokaryotic') {
	    $careful_clusters = purge_short_ESTs_in_clusters($careful_clusters, $single_length);
	}

	#===now start preparing data for different types of input

	#--model_gff3 input
	# identify the models that fall within and between basic clusters
	my ($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($models,
								      $careful_clusters,
								      $seq,
								      );

	#join the clusters on the models
	my $model_clusters = join_clusters_on_pred($models, $careful_clusters, $c_index);
	
	# identify the abinits that fall within and between clusters
	($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
								   $model_clusters,
								   $seq,
								   );
	
	merge_into_cluster($hit_one, $model_clusters, $c_index);
	merge_into_cluster($hit_mult, $model_clusters, $c_index); #these have an internal tag

	#==prep model_gff data
	my @all_data;
	my $c_id = 0;
	#my @gf_data;
	foreach my $c (@{$model_clusters}){
	    my $gf = prep_gff_data($c, $c_id, $seq);
	    push(@all_data, @{$gf}) if defined $gf;

	    $c_id++;
	}

	#--abinit input
	($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
								   $careful_clusters,
								   $seq,
								   );

	#join clusters on the ab-inits
	my $pred_clusters = join_clusters_on_pred($predictions, $careful_clusters, $c_index);

	# identify the abinits that fall within and between clusters
	($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
								   $pred_clusters,
								   $seq,
								   );
	
	merge_into_cluster($hit_one, $pred_clusters, $c_index);
	merge_into_cluster($hit_mult, $pred_clusters, $c_index); #these have an internal tag

	#==prep abinit data
	#my @pr_data;
	foreach my $c (@{$pred_clusters}){
	    my $pr = prep_pred_data($c, $c_id, $seq);
	    push(@all_data, @{$pr}) if defined $pr;

	    $c_id++;
	}

	#--build clusters joined together by models (for hint based annotations)
	my $hint_clusters = [];
	if(@$models){
	    ($p, $m, $x, $z) = PhatHit_utils::separate_by_strand('query', $models);
	    $p_clusters = cluster::shadow_cluster(0, $seq, [@$p,@$p_clusters], $pred_flank);
	    $m_clusters = cluster::shadow_cluster(0, $seq, [@$m,@$m_clusters], $pred_flank);

	    #this method will cause clusters that are near each other and are connected by an orf to merge.
	    #this solves issues with mRNAseq splice site crossing reads and other EST partial exon coverage
	    $p_clusters = join_clusters_around_orf($p_clusters, $seq);
	    $m_clusters = join_clusters_around_orf($m_clusters, $seq);
	    
	    push(@{$hint_clusters}, @{$p_clusters}, @{$m_clusters});
	}
	else{
	    $hint_clusters = $careful_clusters;
	}

	# identify the abinits that fall within and between clusters
	($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
								   $hint_clusters,
								   $seq,
								   );

	merge_into_cluster($hit_one, $hint_clusters, $c_index);
	merge_into_cluster($hit_mult, $hint_clusters, $c_index); #these have an internal tag

	#==prep hint data
	#my @bx_data;
	foreach my $c (@{$hint_clusters}){
	   my $bx = prep_blastx_data($c, $c_id, $seq, $organism_type, $est_forward);
	   push(@all_data, @{$bx}) if defined $bx;

	   $c_id++;
	}

	#add index value to ESTs (corresponds to order in array)
	my $index = 0;
	foreach my $f (@$est_hits){
	    $f->{_index}{est} = $index++;
	}

	#add index value to proteins (corresponds to order in array)
	$index = 0;
	foreach my $f (@$prot_hits){
	    $f->{_index}{prot} = $index++;
	}

	#add index value to preds (corresponds to order in array)
	$index = 0;
	foreach my $f (@$predictions){
	    $f->{_index}{pred} = $index++;
	}

	#add index value to clusters (corresponds to order in array)
	$index = 0;
	foreach my $d (@all_data){
	    $d->{index} = $index++;
	}

	return (\@all_data);
	#return (\@bx_data, \@gf_data, \@pr_data);
}
#------------------------------------------------------------------------
#called in prep_hits to classify abinits as overlaping one, many, or no
#clusters of other evidence
sub segment_preds {
	my $preds            = shift;
	my $careful_clusters = shift;
	my $seq              = shift;

	#identify span of preds
	my @index;
	my $temp_id =0;
	my %p_ab_spans;
	my %m_ab_spans;
	foreach my $s (@{$preds}){
	    $s->{temp_id} = $temp_id;
	    $index[$temp_id] = $s;

	    my ($sB, $sE) = PhatHit_utils::get_span_of_hit($s, 'query');
	    ($sB, $sE) = ($sE, $sB) if $sB > $sE;

	    if($s->strand('query') == 1){
		$p_ab_spans{$temp_id} = [$sB, $sE];
	    }
	    elsif($s->strand('query') == -1){
		$m_ab_spans{$temp_id} = [$sB, $sE];
	    }
	    else{
		die "FATAL: Logic error in segmenting preds\n";
	    }

	    $temp_id++;
	}

	#identify span of clusters
	my $cid = 0;
	my %p_c_spans;
	my %m_c_spans;
	foreach my $c (@{$careful_clusters}){
	    my $cB;
	    my $cE;
	    foreach my $h (@$c){
		my ($hB, $hE) = PhatHit_utils::get_span_of_hit($h, 'query');
		($hB, $hE) = ($hE, $hB) if $hB > $hE;
		$cB = $hB if(!$cB || $hB < $cB);
		$cE = $hE if(!$cE ||$hE > $cE);
	    }

	    if($c->[0]->strand('query') == 1){
		$p_c_spans{$cid} = [$cB, $cE];
	    }
	    elsif($c->[0]->strand('query') == -1){
		$m_c_spans{$cid} = [$cB, $cE];
	    }
	    else{
		die "FATAL: Logic error in segmenting preds\n";
	    }

	    $cid++;
	}

	#identify who each pred overlaps and assign to an overlap group
	my %clusters_hit;
	my @hit_none;
	my @hit_one;
	my @hit_mult;

	my $i_size = $temp_id;
	my $current = 0;
	#on plus strand
	while(my $i = each %p_ab_spans){
	    my $sB = $p_ab_spans{$i}->[0];
	    my $sE = $p_ab_spans{$i}->[1];

	    my $count = 0;
	    while(my $j = each %p_c_spans){
		my $cB = $p_c_spans{$j}->[0];
		my $cE = $p_c_spans{$j}->[1];
		my $class = compare::compare($sB, $sE, $cB, $cE);

		if ($class ne '0'){
		    $clusters_hit{$i}{$j}++;
		    $count++;
		}
	    }

	    if ($count == 0){
		push(@hit_none, $index[$i]);
		$index[$i]->{_hit_multi} = 0;
	    }
	    elsif ($count == 1){
		push(@hit_one, $index[$i]);
		$index[$i]->{_hit_multi} = 0;
	    }
	    else {
		push(@hit_mult, $index[$i]);
		$index[$i]->{_hit_multi} = 1;
	    }
	}

	#on minus strand
	while(my $i = each %m_ab_spans){
	    my $sB = $m_ab_spans{$i}->[0];
	    my $sE = $m_ab_spans{$i}->[1];

	    my $count = 0;
	    while(my $j = each %m_c_spans){
		my $cB = $m_c_spans{$j}->[0];
		my $cE = $m_c_spans{$j}->[1];
		my $class = compare::compare($sB, $sE, $cB, $cE);

		if ($class ne '0'){
		    $clusters_hit{$i}{$j}++;
		    $count++;
		}
	    }

	    if ($count == 0){
		push(@hit_none, $index[$i]);
		$index[$i]->{_hit_multi} = 0;
	    }
	    elsif ($count == 1){
		push(@hit_one, $index[$i]);
		$index[$i]->{_hit_multi} = 0;
	    }
	    else {
		push(@hit_mult, $index[$i]);
		$index[$i]->{_hit_multi} = 1;
	    }
	}

	return (\%clusters_hit, \@hit_one, \@hit_none, \@hit_mult);
}
#------------------------------------------------------------------------
#takes abinits identified by segment_preds and merges them back into
#all the clusters they overlap
sub merge_into_cluster {
   my $hits = shift;
   my $clusters = shift;
   my $c_index = shift;

   foreach my $s (@{$hits}, ){
      my @keys = keys %{$c_index->{$s->{temp_id}}};
      foreach my $i (@keys){
	 push(@{$clusters->[$i]}, $s)
	     unless(grep {$_ eq $s} @{$clusters->[$i]});
      }
   }
}
#------------------------------------------------------------------------
#takes individual clusters and merges them if an abinit identified by
#segment_preds overlaps both clusters, this returns the set of new clusters
sub join_clusters_on_pred {
   my $hits = shift;
   my $clusters = shift;
   my $c_index = shift;

   my @p_clusters;

   foreach my $s (@{$hits}, ){
      my @keys = keys %{$c_index->{$s->{temp_id}}};
      my @new_c = ($s);
      foreach my $i (@keys){
	 push(@new_c, @{$clusters->[$i]});
      }
      push(@p_clusters, \@new_c);
   }

   return \@p_clusters;
}
#------------------------------------------------------------------------
#this method will cause clusters that are near each other and are connected by an orf
#to merge. this solves issues with mRNAseq splice site crossing reads and other ESTs,
#where reads on both extremes of the exon cluster, but there is no joining overage
sub join_clusters_around_orf {
    my $clusters = shift;
    my $seq = shift;

    my @cs = sort {$a->[0]->start('query') <=> $b->[0]->start('query')} @$clusters;
    for(my $i = 0; $i < @cs - 1; $i++){
	#get cluster end
	my $iS = $cs[$i][0]->end('query'); #already will be array space equivilent
	foreach my $h (@{$cs[$i]}){
	    $iS = $h->end('query') if($h->end('query') > $iS);
	}

	#get neighbor cluster start
	my $iE = $cs[$i+1][0]->start('query') - 2;
	foreach my $h (@{$cs[$i+1]}){
	    $iE = $h->start('query') if($h->start('query') < $iE);
	}
	$iE -= 2; #fix for array space

	my $iL = abs($iE - $iS) + 1; #space-in-between length

	next if($iL < 300); #min orf of 300 required (same as most prokaryotic gene finders)

	my $piece = ($cs[$i][0]->strand('query') == 1) ? substr($$seq, $iS, $iL) : Fasta::revComp(substr($$seq, $iS, $iL));
	my $tM = new CGL::TranslationMachine();
	my ($p_seq , $poffset) = $tM->longest_translation($piece);

	#all orf, so merge clusters (merges forward onto next cluster)
	if($poffset < 3 && $iL - (3 * length($p_seq) + $poffset) < 3 && $p_seq !~ /X{10}/){
	    $cs[$i+1] = [@{$cs[$i]}, @{$cs[$i+1]}];
	    $cs[$i] = [];
	}
    }

    @cs = grep {scalar @{$_}} @cs; #remove empty clusters
    
    return \@cs;
}
#------------------------------------------------------------------------
#used to identify hits related to ESTs in a cluster and then remove those
#that are single exon
sub purge_single_exon_hits_in_cluster{
    my $clusters = shift;

    my @c_keepers;
    foreach my $c (@$clusters){
	my @new_c;
	my $ests = [];
	my $preds = [];
	foreach my $h (@$c){
	    if(grep {$h->algorithm =~ /^(.*exonerate\:\:)?$_(\:.*)?$/i} qw(est2genome est_gff blastn tblastx altest_gff)){
		push(@$ests, $h);
	    }
	    elsif(grep {$h->algorithm =~ /^$_(_masked)?(\:.*)?$/i} qw(snap augustus fgenesh genemark pred_gff)){
		push(@$preds, $h);
	    }
	    else{
		push(@new_c, $h);
	    }
	}

	$ests = clean::purge_single_exon_hits($ests);

	push(@new_c, @$ests);

	if(@new_c){
	    push(@new_c, @$preds);
	    push(@c_keepers, \@new_c);
	}
    }

    return \@c_keepers;
}
#------------------------------------------------------------------------
#used to identify hits related to ESTs in a cluster and then remove those
#that have weird splice sites
sub throw_out_bad_splicers_in_cluster{
    my $clusters = shift;
    my $seq = shift;

    my @c_keepers;
    foreach my $c (@$clusters){
	my @new_c;
	my $ests = [];
	my $preds = [];
	foreach my $h (@$c){
	    if(grep {$h->algorithm =~ /^(.*exonerate\:\:)?$_(\:.*)?$/i} qw(est2genome est_gff blastn tblastx altest_gff)){
		push(@$ests, $h);
	    }
	    elsif(grep {$h->algorithm =~ /^$_(_masked)?(\:.*)?$/i} qw(snap augustus fgenesh genemark pred_gff)){
		push(@$preds, $h);
	    }
	    else{
		push(@new_c, $h);
	    }
	}

	$ests = clean::throw_out_bad_splicers($ests, $seq);

	push(@new_c, @$ests);

	if(@new_c){
	    push(@new_c, @$preds);
	    push(@c_keepers, \@new_c);
	}
    }

    return \@c_keepers;
}
#------------------------------------------------------------------------
#used to identify hits related to ESTs in a cluster and then remove those
#below a certain threshold, returns the revised cluster
sub purge_short_ESTs_in_clusters{
    my $clusters = shift;
    my $min = shift;

    my @c_keepers;
    foreach my $c (@$clusters){
	my @new_c;
	my $ests = [];
	my $preds = [];
	foreach my $h (@$c){
	    if(grep {$h->algorithm =~ /^(.*exonerate\:\:)?$_(\:.*)?$/i} qw(est2genome est_gff blastn tblastx altest_gff)){
		push(@$ests, $h);
	    }
	    elsif(grep {$h->algorithm =~ /^$_(_masked)?(\:.*)?$/i} qw(snap augustus fgenesh genemark pred_gff)){
		push(@$preds, $h);
	    }
	    else{
		push(@new_c, $h);
	    }
	}

	$ests = clean::purge_short_single_exons($ests, $min);

	push(@new_c, @$ests);

	if(@new_c){
	    push(@new_c, @$preds);
	    push(@c_keepers, \@new_c);
	}
    }

    return \@c_keepers;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#called by prep_hits for standard evidence clusters
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => each est from best ests
#snap predictions included
sub prep_blastx_data {
	my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;
	my $org_type = shift || 'eukaryotic';
	my $est_forward = shift;

	my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff', 'blastn');
	my $ps_in_cluster    = get_selected_types($c,'protein2genome');
	my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
	my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $models_in_cluster = get_selected_types($c,'model_gff', 'maker');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark', 'pred_gff');
	my @uniq_preds = grep {$_->{_hit_multi} == 0} @$preds_in_cluster;
	
	# groups of most informative protein hits
	# go ahead and inclde the proteion2genome data as well... why not?
	my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	# group of most informative alt splices
	my $gomias = [];

	if($est_forward){
	    $gomias = $ests_in_cluster;
	}
	elsif($org_type eq 'eukaryotic'){
	    $gomias = clean::purge_single_exon_hits($ests_in_cluster);
	    $gomias = clean::get_best_alt_splices($gomias, $seq, 10);
	}
	else{
	    $gomias = PhatHit_utils::make_flat_hits($ests_in_cluster, $seq);
	}

	my @data;
	my $i = 0;
	$gomias->[0] = undef if(! @{$gomias});
	foreach my $mia (@{$gomias}){
	    push(@data, {'gomiph'    => $gomiph,
			 'preds'     => \@uniq_preds,
			 'all_preds' => $preds_in_cluster,
			 'ests'      => $ests_in_cluster,
			 'alt_ests'  => $alt_ests_in_cluster,
			 'mia'       => $mia,
			 'model'     => undef,
			 'gomod'     => $models_in_cluster,
			 'c_id'      => $c_id,
			 'set_id'    => $i++,
			 'type'      => 'bx',
			 }
		 );
	}

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#called by prep_hits for model_gff based clusters
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => none
#predictions included

sub prep_gff_data {
	my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;

	my $models_in_cluster = get_selected_types($c,'model_gff', 'maker');
	my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff', 'blastn');
	my $ps_in_cluster    = get_selected_types($c,'protein2genome');
	my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
	my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark',  'pred_gff');
	my @uniq_preds = grep {$_->{_hit_multi} == 0} @$preds_in_cluster;

	# groups of most informative protein hits
	my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;
	my $i = 0;
	foreach my $model (@{$models_in_cluster}){
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => \@uniq_preds,
			'all_preds' => $preds_in_cluster,
			'ests'      => $ests_in_cluster,
			'alt_ests'  => $alt_ests_in_cluster,
			'mia'       => undef,
			'model'     => $model,
			'gomod'     => undef,
			'c_id'      => $c_id,
			'set_id'    => $i++,
			'type'      => 'gf'
		       }
	       );
	}

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#called by prep_hits for abinit based  clusters
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => none
#predictions included

sub prep_pred_data {
	my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;

	my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff', 'blastn');
	my $ps_in_cluster    = get_selected_types($c,'protein2genome');
	my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
	my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark', 'pred_gff');
	my @uniq_preds = grep {$_->{_hit_multi} == 0} @$preds_in_cluster;

	#model is always first cluster entry
	my $abinits = get_selected_types([$c->[0]],'snap', 'augustus', 'fgenesh',
					 'twinscan', 'genemark', 'pred_gff');
	die "ERROR: Logic problem in maker::auto_annotator::prep_pred_data\n"
	    if(!@$abinits);

	# groups of most informative protein hits
	my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;
	my $i = 0;
	foreach my $pred (@$abinits){
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => \@uniq_preds,
			'all_preds' => $preds_in_cluster,
			'model'     => $pred,
			'gomod'     => undef,
			'ests'      => $ests_in_cluster,
			'alt_ests'  => $alt_ests_in_cluster,
			'mia'       => undef,
			'c_id'      => $c_id,
			'set_id'    => $i++,
			'type'      => 'pr'
			}
		);
        }

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#called by prep_hits for standard evidence clusters
#ests => set of best ests from all ests
#protein homology => each protein exonerate structure
#alternative splice forms => each protein exonerate structure
#no snap predictions included
#sub prep_polpro_data {
#	my $c    = shift;
#	my $c_id = shift;
#	my $seq  = shift;
#
#	my $ests_in_cluster = get_selected_types($c, 'est2genome', 'est_gff', 'blastn');
#	my $ps_in_cluster   = get_selected_types($c,'protein2genome');
#
#	my $possible_ext_sources = combine($ests_in_cluster, $ps_in_cluster);
#
#	my $best_exts = clean::get_best_alt_splices($possible_ext_sources, $seq, 10);
#
#	# group of most informative alt splices
#	my $gomias = clean::get_best_alt_splices($ps_in_cluster, $seq, 10);
#
#	my @data;
#	my $i = 0;
#	foreach my $mia (@{$gomias}){
#		push(@data, {'gomiph' => [$mia],
#			     'ests'   => $best_exts,
#			     'mia'    => $mia,
#			     'c_id'   => $c_id,
#			     'set_id'    => $i++,
#			     'type'   => 'pp'
#			     });
#	}
#
#	return \@data;
#}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#called by prep_hits for standard evidence clusters
#ests => each best est from all ests
#protein homology =>  best set from combined protein exonerate and blastx data
#alternative splice forms => based on each best ests from all ests
#no snap predictions included
#sub prep_polest_data {
#	my $c    = shift;
#	my $c_id = shift;
#	my $seq  = shift;
#
#
#	my $ests_in_cluster = get_selected_types($c, 'est2genome', 'est_gff', 'blastn');
#	my $ps_in_cluster   = get_selected_types($c,'protein2genome');
#	my $bx_in_cluster   = get_selected_types($c,'blastx', 'protein_gff');
#
#	my $i_set      = combine($ps_in_cluster, $bx_in_cluster);
#	my $best_p_set = clean::remove_redundant_alt_splices($i_set, $seq, 10);
#
#	my $best_exts  = clean::get_best_alt_splices($ests_in_cluster, $seq, 10);
#
#	# group of most informative alt splices
#	my $gomias = clean::get_best_alt_splices($ps_in_cluster, $seq, 10);
#
#	my @data;
#	my $i = 0;
#	foreach my $alt_splice_est (@{$best_exts}){
#		push(@data, {'gomiph' => $best_p_set,
#			     'ests'   => [$alt_splice_est],
#			     'mia'    => $alt_splice_est,
#			     'c_id'   => $c_id,
#			     'set_id'    => $i++,
#			     'type'   => 'pe'
#			     });
#	}
#	return \@data;
#}
#------------------------------------------------------------------------
#this subroutine returns finished transcripts for all predictors.
#called outside of package by maker
sub annotate_trans {
    my $v_seq_ref        = shift;
    my $m_seq_ref        = shift;
    my $def              = shift;
    my $seq_id           = shift;
    my $all_data         = shift; #evidence clusters
    my $the_void         = shift;
    my $CTL_OPT          = shift;
    $LOG                 = shift;

    my %transcripts;
    foreach my $dc (@$all_data) {
	#---model passthrough here
	if($dc->{type} eq 'gf'){
	    my $trans = run_it([$dc],
			       $the_void,
			       $m_seq_ref,
			       $v_seq_ref,
			       $def,
			       'model_gff',
			       $CTL_OPT
			       );

	    push(@{$transcripts{'model_gff'}}, @$trans);
	}
	#---hint based gene prediction here (includes est2genome)
	elsif($dc->{type} eq 'bx'){
	    foreach my $prdr (@{$CTL_OPT->{_predictor}}){
		next if($prdr eq 'model_gff' || $prdr eq 'pred_gff');
		
		my $trans = run_it([$dc],
				   $the_void,
				   $m_seq_ref,
				   $v_seq_ref,
				   $def,
				   $prdr,
				   $CTL_OPT
				   );
		
		push(@{$transcripts{$prdr}}, @$trans); 
	    }
	}
	#---abinit scoring here
	elsif($dc->{type} eq 'pr'){
	    my $trans = run_it([$dc],
			       $the_void,
			       $m_seq_ref,
			       $v_seq_ref,
			       $def,
			       'abinit', #all abinits not just pred_gff
			       $CTL_OPT
			       );

	    #add abinits to their predictor type after they have been proccessed
	    #remeber they get treated differently so you want to add them as a
	    #separate step from hint based predictions.
	    foreach my $s (@$trans){
		my $t = $s->[0];
		my $al = lc($t->algorithm);
		$al =~ s/_masked$//;
		$al =~ s/(pred_gff).*$/$1/;

		if($al =~ /^snap$|^genemark$|^augustus$|^fgenesh$/){
		    push(@{$transcripts{"$al\_abinit"}}, $s);
		}
		elsif($al =~ /^pred_gff$/){
		    push(@{$transcripts{$al}}, $s);
		}
		else{
		    die "ERROR: Not a supported algorithm: ".$t->algorithm."\n";
		}
	    }
	}
    }
    
    return \%transcripts;
}

#------------------------------------------------------------------------
#this subroutine returns seperates and groups transcripts for all predictors.
#called outside of package by maker
sub annotate_genes {
    my $transcripts      = shift;
    my $all_data         = shift; #evidence clusters
    my $v_seq_ref        = shift;
    my $seq_id           = shift;
    my $chunk_number     = shift; #required to name genes for each chunk
    my $the_void         = shift;
    my $CTL_OPT          = shift;

    #reset gene names
    #my $GFF_DB = new GFFDB($CTL_OPT);
    $SEEN = {};#$GFF_DB->get_existing_gene_names($seq_id);

    my %annotations;
    while (my $key = each %{$transcripts}){
	$annotations{$key} = group_transcripts($transcripts->{$key},
					       $all_data,
					       $v_seq_ref,
					       $seq_id,
					       $chunk_number,
					       $key,
					       $the_void,
					       $CTL_OPT
					       );
    }

    return \%annotations;
}
#------------------------------------------------------------------------
sub annotate_stats {
    my $annots    = shift;
    my $seq       = shift;
    my $CTL_OPT   = shift;

    my %annotations;
    while (my $key = each %{$annots}){
	foreach my $g (@{$annots->{$key}}){
	    my $evidence = $g->{g_evidence};
	    my $g_name = $g->{g_name};

	    #load transcript stats
	    my $AED = 1;
	    my $eAED = 1;
	    my $i = 1;
	    my @t_structs;
	    foreach my $s (@{$g->{t_structs}}) {
		my $t_struct = load_transcript_stats($s, $g_name, $i, $evidence, $seq, $CTL_OPT);
		
		my $bad = 1 if($t_struct->{p_length} <= $CTL_OPT->{min_protein} ||
			       $t_struct->{eAED} > $CTL_OPT->{AED_threshold}
			       );

		next if($bad && $key !~ /_abinit$/);

		push(@t_structs, $t_struct);

		$t_struct->{hit}->{_REMOVE} = 1 if($bad);
		$AED = $t_struct->{AED} if($t_struct->{AED} < $AED);
		$eAED = $t_struct->{eAED} if($t_struct->{eAED} < $eAED);
		$i++;
	    }

	    next if(!@t_structs);

	    delete($g->{g_evidence}); #remove evidence (compact for transmission)
	    $g->{t_structs} = \@t_structs;
	    $g->{AED}       = $AED;
	    $g->{eAED}      = $eAED;

	    push(@{$annotations{$key}}, $g);
	}
    }
    add_abAED(\%annotations);    

    return \%annotations;
}
#------------------------------------------------------------------------
sub annotate_finalize {
    my $annotations = shift;

}
#------------------------------------------------------------------------
#counts the number of nucleotides in a hit on the query sequence
sub get_n_count{
    my $g = shift;

    my ($start, $end) = ($g->{g_start}, $g->{g_end});
    my $offset = $start;
    my $length = $end - $start + 1;
    my @b_seq = map {0} (1..$length);

    foreach my $s (@{$g->{t_structs}}){
	my $hit = $s->{hit};
	my @hsps = $hit->hsps() if defined $hit->hsps();
	foreach my $hsp (@hsps){
	    my $s = $hsp->start('query') - $offset;
	    my $e = $hsp->end('query') - $offset;

	    #array space coors
	    die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	    die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

	    @b_seq[$s..$e] = map {1} ($s..$e);
	}
    }

    my @index = (0, 0);

    foreach my $i (@b_seq){
	$index[$i]++;
    }

    return $index[1];
}
#------------------------------------------------------------------------
#returns all phathits inside a gene annotation
sub gene2all_phathits {
    my $g = shift;

    my @hits;
    foreach my $t (@{$g->{t_structs}}){
	push (@hits, $t->{hit});
    }

    return \@hits;
}
#------------------------------------------------------------------------
#returns all p_bases (transcript phat_hits without UTRs added by pneu)
#for a gene annotation, called for abAED calculation
sub gene2allPbases {
    my $g = shift;

    my @hits;
    foreach my $t (@{$g->{t_structs}}){
	my $f = $t->{p_base};
	$f = $t->{hit} if(! $f);
	push (@hits, $f);
    }

    return \@hits;
}
#------------------------------------------------------------------------
#adds the abAED to all transcripts within a gene annotation.
#abAED is the average AED across alternate models.
#It is  caculated by comaring each of the p_bases or CDSs
#of alternate gene models back to the current gene model.
#lower abAED means that this model is closer to the consensus of all the gene
#models. est2genome is not considered as an alternate model when calculating
#abAED; however, abinits, and model_gff annotations are.  This means abAED is
#calculated for est2genome, but not against est2genome.
sub add_abAED{
    my $annotations = shift;

    my @p_bag;
    my @m_bag;
    my @p_genes;
    my @m_genes;

    #collect all no UTR base models to get abAED
    foreach my $p (keys %$annotations){
	foreach my $g (@{$annotations->{$p}}){
	    if($g->{g_strand} == 1){
		push(@p_genes, $g);#calculate for est2genome
		next if($p eq 'est2genome');#but not against est2genome
		my $hits = gene2allPbases($g);
		push(@p_bag, @$hits);
	    }
	    elsif($g->{g_strand} == -1) {
		push(@m_genes, $g);#calculate for est2genome
		next if($p eq 'est2genome');#but not against est2genome
		my $hits = gene2allPbases($g);
		push(@m_bag, @$hits);
	    }
	    else{
		die "ERROR: Logic error in auto_annotator::best_annotations\n";
	    }
	}
    }

    #calculate abAED plus strand
    foreach my $g (@p_genes){
	my $hits = get_hits_overlapping_gene($g, \@p_bag);
	my $abAED = 1;
	foreach my $t (@{$g->{t_structs}}){
	    my $f = $t->{p_base};
	    $f = $t->{hit} if(! $f);
	    my $ab = shadow_AED::get_abAED($hits, $f);
	    $t->{abAED} = $ab;
	    $abAED = $ab if($ab < $abAED);
	}
	$g->{abAED} = $abAED;
    }

    #calculate abAED minus strand
    foreach my $g (@m_genes){
	my $hits = get_hits_overlapping_gene($g, \@m_bag);
	my $abAED = 1;
	foreach my $t (@{$g->{t_structs}}){
	    my $f = $t->{p_base};
	    $f = $t->{hit} if(! $f);
	    my $ab = shadow_AED::get_abAED($hits, $f);
	    $t->{abAED} = $ab;
	    $abAED = $ab if($ab < $abAED);
	}
	$g->{abAED} = $abAED;
    }
}
#------------------------------------------------------------------------
#this subrutine returns finished MAKER annotations
#it looks through all sets of annotations and returns the best subset based
#on AED and abAED.  Only returns for predictors specified. Predictor order
#is important if AED and abAED are the same for an overlapping annotation!

sub best_annotations {
    my $annotations = shift;
    my $CTL_OPT = shift;

    my @predictors = @{$CTL_OPT->{_predictor}};

    my @p_keepers;
    my @m_keepers;

    #keep all gff3 passthrough if there's nothing else
    if(@predictors == 1 && $predictors[0] eq 'model_gff'){
	my @final;
	foreach my $g (@{$annotations->{'model_gff'}}){
	    if($g->{g_strand} == 1){
		push(@final, $g);
	    }
	    elsif($g->{g_strand} == -1){
		push(@final, $g);
	    }
	}

	return \@final;
    }
    #keep all est2genome genes if mapping forward onto a new assembly
    elsif($CTL_OPT->{est_forward} &&
	  @predictors == 1 &&
	  $predictors[0] eq 'est2genome'
	  ){
	my @final;
	foreach my $g (@{$annotations->{'est2genome'}}){
	    if($g->{g_strand} == 1){
		push(@final, $g);
	    }
	    elsif($g->{g_strand} == -1){
		push(@final, $g);
	    }
	}

	return \@final;
    }
    elsif(@predictors){
	#set up lists for plus and minus strands as well as possible mergers
	#predictor types are processed in the order given by control files
	my $p_list = [];
	my $m_list = [];
	my @p_est2g;
	my @m_est2g;
	foreach my $p (@predictors){
	    my $pa = "$p\_abinit";
	    my @hints = @{$annotations->{$p}} if($annotations->{$p});
	    my @abinits = @{$annotations->{$pa}} if($annotations->{$pa} && $p !~ /2genome$/);

	    foreach my $g (@hints, @abinits){
		next if($g->{t_structs}->[0]->{hit}->{_REMOVE}); #added to filter low support abinits

		if($p ne 'est2genome' && $p ne 'protein2genome' && $g->{g_strand} == 1){
		    push(@$p_list, $g) if(($g->{eAED} < 1 && $g->{eAED} <= $CTL_OPT->{AED_threshold})
					  || $p eq 'model_gff'
					  );
		}
		elsif($p ne 'est2genome' && $p ne 'protein2genome' && $g->{g_strand} == -1) {
		    push(@$m_list, $g) if(($g->{eAED} < 1 && $g->{eAED} <= $CTL_OPT->{AED_threshold})
					  || $p eq 'model_gff'
					  );
		}
		elsif($g->{g_strand} == 1){
		    #separate est2genome and protein2genome genes
		    push(@p_est2g, $g) if($g->{eAED} < 1 && $g->{eAED} <= $CTL_OPT->{AED_threshold});
		}
		elsif($g->{g_strand} == -1){
		    #separate est2genome and protein2genome genes
		    push(@m_est2g, $g) if($g->{eAED} < 1 && $g->{eAED} <= $CTL_OPT->{AED_threshold});
		}
		else{
		    die "ERROR: Logic error in auto_annotator::best_annotations\n";
		}
	    }
	}

	#remove low scoring overlaping genes
	@$p_list  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @$p_list;
	@$m_list  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @$m_list;
	push(@$p_list, @p_est2g); #est2genome added to end, will only appear if nothing else overlaps
	push(@$m_list, @m_est2g); #est2genome added to end, will only appear if nothing else overlaps
	$p_list = _best($p_list);
	$m_list = _best($m_list);

	#almost final  sets
	push(@p_keepers, @$p_list);
	push(@m_keepers, @$m_list);
    }

    #remove CDS competition on opposite strand
    my $final = remove_CDS_competitors(\@p_keepers, \@m_keepers);

    return $final;
}
#------------------------------------------------------------------------
#filter for CDS competition on opposite strands
#thow out genes that compete on opposite strands and
#don't have abinit, exonerate, or spliced EST confirmation
sub remove_CDS_competitors {
    my $plus = shift;
    my $minus = shift;

    my @p_final;
    my @m_final;
    my @p_maybe;
    my @m_maybe;

    #--find maybe suspects on plus strand
    foreach my $g (@$plus){
	my $add_flag = 0;

	#always keep gff passthrough annotations
	if($g->{predictor} eq 'model_gff'){
	    push(@p_final, $g);
	    $add_flag = 1;
	    next;
	}

	#test each gene for abinit, spliced EST, or exonerate protein support
	foreach my $t (@{$g->{t_structs}}){
	    #check EST splice sites (fast)
	    my @qi = split(/\|/, $t->{t_qi});
	    if($qi[1] > 0){
		push(@p_final,$g);
		$add_flag = 1;
		last;
	    }

	    #then check exact exon overlap (exonerate est/protein)
	    if($qi[2] > 0){
                push(@p_final,$g);
                $add_flag = 1;
                last;
	    }

	    #then check abinits
	    if($qi[4] > 0 || $qi[5] > 0){
		push(@p_final,$g);
		$add_flag = 1;
		last;
	    }
	}

	#add suspects to maybe list
	push(@p_maybe, $g) if(! $add_flag);
    }

    #--find maybe suspects on minus strand
    foreach my $g (@$minus){
	my $add_flag = 0;

	#always keep gff passthrough annotations
	if($g->{predictor} eq 'model_gff'){
	    push(@m_final, $g);
	    $add_flag = 1;
	    next;
	}

	#test each gene for abinit, spliced EST, or exonerate protein support
	foreach my $t (@{$g->{t_structs}}){
	    #check EST splice sites (fast)
	    my @qi = split(/\|/, $t->{t_qi});
	    if($qi[1] > 0){
		push(@m_final,$g);
		$add_flag = 1;
		last;
	    }

	    #then check exact exon overlap (exonerate est/protein)
	    if($qi[2] > 0){
                push(@m_final,$g);
                $add_flag = 1;
                last;
	    }

	    #then check abinits
	    if($qi[4] > 0 || $qi[5] > 0){
		push(@m_final,$g);
		$add_flag = 1;
		last;
	    }
	}

	#add suspects to maybe list
	push(@m_maybe, $g) if(! $add_flag);
    }


    #--check if plus strand maybes overlap final set minus
    my @suspects;
    if(@m_final){ #skip if nothing to compare to
	foreach my $s (@p_maybe){
	    my $sB = $s->{g_start};
	    my $sE = $s->{g_end};
	    ($sB, $sE) = ($sE, $sB) if($sE < $sB);

	    my $bad;
	    foreach my $g (@m_final){
		my $gB = $g->{g_start};
		my $gE = $g->{g_end};
		($gB, $gE) = ($gE, $gB) if($gE < $gB);

		my $comp = compare::compare ($sB,
					     $sE,
					     $gB,
					     $gE
					    );
		if($comp ne '0'){
		    $bad = 1;
		    last;
		}
	    }

	    #keep as suspect if not overlapping final set
	    push(@suspects, $s) if(!$bad);
	}
    }
    else{
	push(@suspects, @p_maybe);#everyones a supsect
    }

    #--check if minus strand maybes overlap final set plus
    if(@p_final){ #skip if nothing to compare to
	foreach my $s (@m_maybe){
	    my $sB = $s->{g_start};
	    my $sE = $s->{g_end};
	    ($sB, $sE) = ($sE, $sB) if($sE < $sB);

	    my $bad;
	    foreach my $g (@p_final){
		my $gB = $g->{g_start};
		my $gE = $g->{g_end};
		($gB, $gE) = ($gE, $gB) if($gE < $gB);

		my $comp = compare::compare ($sB,
					     $sE,
					     $gB,
					     $gE
					    );
		if($comp ne '0'){
		    $bad = 1;
		    last;
		}
	    }

	    #keep as suspect if not overlapping final set
	    push(@suspects, $s) if(!$bad);
	}
    }
    else{
	    push(@suspects, @m_maybe);#everyones a supsect
    }

    #remove suspects that overlap other better suspects
    @suspects = sort {crit4($b) <=> crit4($a)} @suspects;
    @suspects = @{_best(\@suspects)};

    #return final annotations
    my $final;
    push(@$final, @p_final);
    push(@$final, @m_final);
    push(@$final, @suspects);

    return $final;
}
#------------------------------------------------------------------------
#called by best_annotations, returns non_overlapping annotion set.
#expects the input array to already be sorted so that the first annotations
#are more likely to be returned than the last. Expects all annotations to
#given to be from the same strand and will not double check this.
sub _best{
    my $list = shift;

    my @keepers;
    foreach my $g (@$list){
	my $g_B = $g->{g_start};
	my $g_E = $g->{g_end};

	my $bad;
	foreach my $k (@keepers){
	    my $k_B = $k->{g_start};
	    my $k_E = $k->{g_end};

	    my $class = compare::compare($g_B, $g_E, $k_B, $k_E);

	    if($class ne '0'){
		$bad = 1;
		last;
	    }
	}

	push(@keepers, $g) if(! $bad);
    }

    return \@keepers;
}
#------------------------------------------------------------------------
#sort by combined abinit-evidence AED score
sub crit1 {
   my $g = shift;

   return $g->{eAED} + ($g->{abAED} * 1/3);
}
#------------------------------------------------------------------------
#sort by evidence AED score
sub crit2 {
    my $g = shift;

    return $g->{AED};
}
#------------------------------------------------------------------------
#sort by abAED score
sub crit3 {
    my $g = shift;

    return $g->{abAED};
}
#------------------------------------------------------------------------
#sort by gene length (note replace this with fathom score some day)
sub crit4 {
    my $g = shift;

    my $length = length($g->{t_structs}->[0]->{p_seq});

    return $length;
}
#------------------------------------------------------------------------
#called by subroutine annotate.
#takes evidence clusters prepared by prep_hits and runs gene predictons
#using hints from that evidence.  For model_gff and abinit clusters it
#just returns the abinits with UTR added or model_gffs as they are.
#It also returns p_bases or the models before UTR was added.
#est2genome predictions are also built here by adding UTR to the best
#spliced ESTs

sub run_it {
    my $data         = shift;
    my $the_void     = shift;
    my $seq          = shift;
    my $v_seq        = shift;
    my $def          = shift;
    my $predictor    = shift;
    my $CTL_OPT      = shift;

    my $q_id = Fasta::def2SeqID($def);
    my @transcripts;
    foreach my $set (@{$data}) {
	my $mia      = $set->{mia};
	my $ests     = $set->{ests};
	my $model    = $set->{model};
	my $gomiph   = $set->{gomiph};
	my $blastx   = get_selected_types($gomiph,'blastx', 'protein_gff');
	my $pol_p    = get_selected_types($gomiph,'protein2genome');
	my $alt_ests = $set->{alt_ests};
	my $preds    = $set->{preds};
	my $all_preds = $set->{all_preds};

	#------gff passthrough
	if ($predictor eq 'model_gff') {
	    next if(! defined $model);
	    my $transcript = $model;

	    push(@transcripts, [$transcript, $set->{index}, undef]);

	    next;
	}

	#------ab-init passthrough
	if ($predictor eq 'abinit') {
	    next if(! defined $model);

	    #added 2/23/2009 to reduce spurious gene predictions with only single exon blastx support
	    my $remove;
	    if($CTL_OPT->{organism_type} eq 'eukaryotic'){
		$remove = 1;
		#make sure the spliced EST evidence actually overlaps
		if($remove && defined $mia){
		    my $mAED = shadow_AED::get_eAED([$mia],$model); #verifies at least a splice
		    $remove = 0 if($mAED < 1);
		}

		#check all ESTs for splice support if $mia does not exist
		if($remove && @$ests){
		    my $mAED = shadow_AED::get_eAED($ests, $model);
		    $remove = 0 if($mAED < 1);
		}

		#make sure the polished protein evidence actually overlaps
		if($remove && (@$pol_p > 1 || (@$pol_p == 1 && $pol_p->[0]->hsps > 1))){
		    my $pAED = shadow_AED::get_eAED($pol_p, $model); #veifies reading frame
		    $remove = 0 if($pAED < 1);
		}

		#make sure the alt est evidence is not single exon and actually overlaps
		if($remove && @$alt_ests){
		    my $clean  = clean::purge_single_exon_hits($alt_ests);
		    my $aAED = (! @$clean) ? 1 : shadow_AED::get_eAED($clean,$model); #splice site and overlap check
		    $remove = 0 if ($aAED < 1);
		}

		#make sure blastx evidence is sufficient
		if($remove && @$blastx){
		    my $coors  = PhatHit_utils::get_hsp_coors($blastx, 'query');
		    my $pieces = Shadower::getPieces($seq, $coors, 0);

		    if(@$pieces <= 1 && $model->hsps <= 2){ # if single exon evidence model should be close
			my $abAED = shadow_AED::get_abAED($all_preds, $model);

			if($abAED <= 0.25){
			    my $bAED = shadow_AED::get_eAED($blastx, $model); #also verifies reading frame
			    $remove = 0 if($bAED <= 0.25);
			}
		    }
		    elsif(@$pieces > 1){
			$remove = 0;
		    }
		}
	    }

	    #add UTR to ab-inits
	    my $select = $model;
	    my $transcript = pneu($ests, $select, $seq); #helps tile ESTs
	    while(! compare::is_same_alt_form($select, $transcript, $seq, 0)){
		$select = $transcript;
		$transcript = pneu($ests, $select, $seq); #helps tile ESTs
	    }

	    #don't filter imediately just mark for downstream filtering
	    $transcript->{_REMOVE} = $remove;

	    push(@transcripts, [$transcript, $set->{index}, $model]);

	    next;
	}

	#------est2genome
	if ($predictor eq 'est2genome') {
	    next if (! defined $mia);

	    my $transcript = $mia;
	    if($CTL_OPT->{est_forward}){
		my $select = $mia;
		$transcript = pneu($ests, $select, $seq); #helps tile ESTs
		while(! compare::is_same_alt_form($select, $transcript, $seq, 0)){
		    $select = $transcript;
		    $transcript = pneu($ests, $select, $seq); #helps tile ESTs
		}
	    }
	    $transcript->{_HMM} = 'est2genome';

	    if(! $CTL_OPT->{est_forward} && $CTL_OPT->{organism_type} eq 'prokaryotic'){
		my $transcript_seq  = get_transcript_seq($transcript, $seq);
		my ($translation_seq, $offset, $end, $has_start, $has_stop) = get_translation_seq($transcript_seq, $transcript);
		#at least 60% of EST must be CDS to make a gene prediction
		next if((length($translation_seq)+1) * 3 / length($transcript_seq) < .60 || ! $has_stop);
	    }

	    next if !$transcript;

	    $transcript->{_tran_name} = $mia->name if($CTL_OPT->{est_forward});

	    push(@transcripts, [$transcript, $set->{index}, $mia]);

	    next;
	}

	#------protein2genome
	if ($predictor eq 'protein2genome') {
	    next if($CTL_OPT->{organism_type} eq 'eukaryotic');
	    next if(! @$gomiph);

	    my $utr = [];
	    push(@$utr, $mia) if($mia);
	    my $miphs = PhatHit_utils::make_flat_hits($gomiph, $seq);

	    foreach my $miph (@$miphs){
		my $transcript_seq  = get_transcript_seq($miph, $seq);
		my ($translation_seq, $offset, $end, $has_start, $has_stop) = get_translation_seq($transcript_seq, $miph);

		#at least 80% of protein must be CDS to make a gene prediction
		next if(length($translation_seq) * 3 / length($transcript_seq) < .80);

		#adjust CDS pre-UTR identification
		my $copy = PhatHit_utils::adjust_start_stop($miph, $v_seq);
		$copy = PhatHit_utils::clip_utr($copy, $v_seq);

		my $select = $copy;
		my $transcript = pneu($utr, $select, $seq); #helps tile ESTs
		while(! compare::is_same_alt_form($select, $transcript, $seq, 0)){
		    $select = $transcript;
		    $transcript = pneu($ests, $select, $seq); #helps tile ESTs
		}
		$transcript->{_HMM} = 'protein2genome';

		next if(! $transcript);

		push(@transcripts, [$transcript, $set->{index}, $miph]);
	    }

	    next;
	}

	#------genemark does not have hints enabled
	return [] if ($predictor eq 'genemark');

	#------default hint based behavior
	my ($pred_shots, $strand) = get_pred_shot($seq,
						  $def,
						  $the_void,
						  $set,
						  $predictor,
						  $CTL_OPT,
						  $LOG
						  );

	my $on_right_strand = get_best_pred_shots($strand, $pred_shots);
	
	#only keep multi-exon hint based predictions single exon prediction
	#are more likely to be spurious if hint based, these are better
	#derived from the ab initio predictions
	@$on_right_strand = grep {$_->hsps > 1} @$on_right_strand if($CTL_OPT->{organism_type} eq 'eukaryotic');

	#added 2/23/2009 to reduce spurious gene predictions with only single exon blastx suport
	if($CTL_OPT->{organism_type} eq 'eukaryotic' && @$on_right_strand){
	    my @keepers;
	    
	    my $clean;
	    my $pieces;
	    foreach my $h (@$on_right_strand){
		my $remove = 1;

		#make sure the spliced EST evidence actually overlaps
		if($remove && defined $mia){
		    my $mAED = shadow_AED::get_eAED([$mia],$h);
		    $remove = 0 if($mAED < 1);
		}

		#check all ESTs for splice support if $mia does not exist
		if($remove && @$ests){
		    my $mAED = shadow_AED::get_eAED($ests,$h);
		    $remove = 0 if($mAED < 1);
		}
		
		#make sure the polished protein evidence actually overlaps
		if($remove && (@$pol_p > 1 || (@$pol_p == 1 && $pol_p->[0]->hsps > 1))){
		    my $pAED = shadow_AED::get_eAED($pol_p, $h);
		    $remove = 0 if($pAED < 1);
		}
		
		#make sure the alt est evidence is not single exon and actually overlaps
		if($remove && @$alt_ests){
		    $clean  = clean::purge_single_exon_hits($alt_ests) if(! $clean); #only calculate once
		    my $aAED = (! @$clean) ? 1 : shadow_AED::get_eAED($clean,$h);
		    $remove = 0 if ($aAED < 1);
		}
		
		#make sure blastx evidence is sufficient
		if($remove && @$blastx){
		    if(! $pieces){ # only calculate once
			my $coors  = PhatHit_utils::get_hsp_coors($blastx, 'query');
			$pieces = Shadower::getPieces($seq, $coors, 0);
		    }
		    
		    if(@$pieces <= 1 && $h->hsps <= 2){ # if single exon evidence then model should be close
			#make sure ab initio evidence can support a single exon alignment
			#this step not needed in ab inits because the test model is an abinit model
			if(grep {$_->hsps <= 2} @$all_preds){
			    my $abAED = shadow_AED::get_abAED($all_preds, $h);
			    
			    if($abAED <= 0.25){
				my $bAED = shadow_AED::get_eAED($blastx, $h);
				$remove = 0 if($bAED <= 0.25);
			    }
			}
		    }
		    elsif(@$pieces > 1){
			$remove = 0;
		    }
		}
		push(@keepers, $h) if(! $remove);
	    }

	    @$on_right_strand = @keepers;
	}

	#add transcripts
	next if(!@{$on_right_strand});
	foreach my $pred_shot (@{$on_right_strand}) {
	    if (defined($pred_shot)){
		my $transcript = $pred_shot;
		if(defined($mia)){
		    my $select = ($CTL_OPT->{alt_splice}) ? pneu([$mia], $transcript, $seq) : $transcript;
		    $transcript = pneu($ests, $select, $seq); #helps tile ESTs
		    while(! compare::is_same_alt_form($select, $transcript, $seq, 0)){
			$select = $transcript;
			$transcript = pneu($ests, $select, $seq); #helps tile ESTs
		    }
		}

		push(@transcripts, [$transcript, $set->{index}, $pred_shot]);
	    }
	}
    }

    return \@transcripts;
}
#------------------------------------------------------------------------
#runs the gene prdictors with hints.  Called by run_it.
sub get_pred_shot {
   my $seq         = shift;
   my $def         = shift;
   my $the_void    = shift;
   my $set         = shift;
   my $predictor   = shift;
   my $CTL_OPT    = shift;
   my $LOG         = shift;

   my $strand;
   my @all_preds;

   if($predictor eq 'snap'){
       foreach my $entry (split(',', $CTL_OPT->{snaphmm})){
	   my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
	   my $pred_command = $CTL_OPT->{snap}.' '.$hmm;
	   (my $preds, $strand) = Widget::snap::get_pred_shot($seq,
							      $def,
							      $the_void,
							      $set,
							      $CTL_OPT->{pred_flank},
							      $pred_command,
							      $hmm,
							      $CTL_OPT->{alt_splice},
							      $CTL_OPT->{force},
							      $LOG
							      );
	   foreach my $p (@$preds){
	       $p->{_HMM} = $hmm;
	       $p->{_label} = $label if($label);
	       map {$_->{_label} = $label} $p->hsps if($label);
	   }

	   push(@all_preds, @$preds);
       }
   }
   elsif($predictor eq 'augustus'){
       foreach my $entry (split(',', $CTL_OPT->{augustus_species})){
	   my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
	   my $pred_command = $CTL_OPT->{augustus}.' --species='.$hmm;
	   (my $preds, $strand) = Widget::augustus::get_pred_shot($seq,
								  $def,
								  $the_void,
								  $set,
								  $CTL_OPT->{pred_flank},
								  $pred_command,
								  $hmm,
								  $CTL_OPT->{alt_splice},
								  $CTL_OPT->{force},
								  $LOG
								  );
	   foreach my $p (@$preds){
	       $p->{_HMM} = $hmm;
	       $p->{_label} = $label if($label);
	       map {$_->{_label} = $label} $p->hsps if($label);
	   }

	   push(@all_preds, @$preds);
       }
   }
   elsif($predictor eq 'fgenesh'){
       foreach my $entry (split(',', $CTL_OPT->{fgenesh_par_file})){
	   my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
	   my $pred_command = $CTL_OPT->{fgenesh}.' '.$hmm;
	   (my $preds, $strand) = Widget::fgenesh::get_pred_shot($seq,
								 $def,
								 $the_void,
								 $set,
								 $CTL_OPT->{pred_flank},
								 $pred_command,
								 $hmm,
								 $CTL_OPT->{alt_splice},
								 $CTL_OPT->{force},
								 $LOG
								 );
	   foreach my $p (@$preds){
	       $p->{_HMM} = $hmm;
	       $p->{_label} = $label if($label);
	       map {$_->{_label} = $label} $p->hsps if($label);
	   }

	   push(@all_preds, @$preds);
       }
   }
   else{
      die "ERROR: Not a valid predictor in auto_annoator::get_pred_shot\n";
   }

   return (\@all_preds, $strand);
}
#------------------------------------------------------------------------
#takes the gene predictions and evidence and builds transcript name,
#QI, AED and gets protein and mRNA sequence then puts it al into a
#HASH.  Called by group transcripts. Transcript name is set here, as
#well as final UTR boudaries
sub load_transcript_struct {
	my $f            = shift;
	my $g_name       = shift;
	my $i            = shift;
	my $seq          = shift;
	my $p_base       = shift;
	my $predictor    = shift;
	my $CTL_OPT      = shift;

	my $transcript_seq  = get_transcript_seq($f, $seq);
	my ($translation_seq, $offset, $end, $has_start, $has_stop) = get_translation_seq($transcript_seq, $f);

	if($predictor !~ /model_gff/ && !$CTL_OPT->{est_forward}){
	    #fix for non-canonical (almost certainly bad) 5' and 3' UTR
	    my $trim5 = (!$has_stop && $end != length($transcript_seq) + 1);
	    my $trim3 = (!$has_start && $offset != 0);
	    if($trim5 || $trim3){
		$f = PhatHit_utils::_clip($f, $seq, $trim5, $trim3); #WARNING: this removes any non-standard values added to the object hash
		$transcript_seq  = get_transcript_seq($f, $seq);
		($translation_seq, $offset, $end, $has_start, $has_stop) = get_translation_seq($transcript_seq, $f);
	    }

	    #walk out edges to force completion
	    if($CTL_OPT->{always_complete} && (!$has_start || !$has_stop)){
		$f = PhatHit_utils::adjust_start_stop($f, $seq);
		$transcript_seq  = get_transcript_seq($f, $seq);
		($translation_seq, $offset, $end, $has_start, $has_stop) = get_translation_seq($transcript_seq, $f);
	    }
	}

	my $len_3_utr = length($transcript_seq) - $end + 1;
	my $l_trans =  length($translation_seq);

	my $t_name = ($f->{_tran_name}) ? $f->{_tran_name} : "$g_name-mRNA-$i"; #affects GFFV3.pm
	my $t_id = ($f->{_tran_id}) ? $f->{_tran_id} : $t_name; #affects GFFV3.pm
	$f->name($t_name);

	my $t_struct = {'hit'       => $f,
			't_seq'     => $transcript_seq,
			'p_seq'     => $translation_seq,
			't_offset'  => $offset,
			't_end'     => $end,
			't_name'    => $t_name,
			't_id'      => $t_id,
			'has_start' => $has_start,
			'has_stop'  => $has_stop,
			'p_base'    => $p_base,
			'p_length'  => length($translation_seq)
		    };

	return $t_struct;
}
#------------------------------------------------------------------------
#takes the gene predictions and evidence and builds transcript name,
#QI, AED and gets protein and mRNA sequence then puts it al into a
#HASH.  Called by group transcripts. Transcript name is set here, as
#well as final UTR boudaries
sub load_transcript_stats {
	my $struct       = shift;
	my $g_name       = shift;
	my $i            = shift;
	my $evi          = shift;
	my $seq          = shift;
	my $CTL_OPT      = shift;

	my $f               = $struct->{hit};
	my $p_base          = $struct->{p_base};
	my $offset          = $struct->{t_offset};
	my $end             = $struct->{t_end};
	my $transcript_seq  = $struct->{t_seq};
	my $translation_seq = $struct->{p_seq};

	my $len_3_utr = length($transcript_seq) - $end + 1;
	my $l_trans   = length($translation_seq);

	my $pol_p_hits   = get_selected_types($evi->{gomiph}, 'protein2genome');
	my $pol_e_hits   = get_selected_types($evi->{ests}, 'est2genome', 'est_gff', 'blastn');
	my $blastx_hits  = get_selected_types($evi->{gomiph},'blastx', 'protein_gff');
	my $tblastx_hits = get_selected_types($evi->{alt_ests},'tblastx', 'altest_gff');
	my $abinits      = $evi->{all_preds};

	my @bag = (@$pol_p_hits,
		   @$pol_e_hits,
		   @$blastx_hits,
		   @$tblastx_hits
		  );

	#evidence AED
	my $AED  = shadow_AED::get_AED(\@bag, $f);
	my $eAED = shadow_AED::get_eAED(\@bag, $f, $seq);
	my $qi   = maker::quality_index::get_transcript_qi($f,$evi,$offset,$len_3_utr,$l_trans);
	$f->{_AED}  = $AED;
	$f->{_eAED} = $eAED;

	if($p_base && $p_base->algorithm !~ /est2genome|est_gff|protein2genome|protein_gff/){
	    $p_base->name($f->name);
	    $p_base->{_AED} = shadow_AED::get_AED(\@bag, $p_base);
	    $p_base->{_eAED} = shadow_AED::get_eAED(\@bag, $p_base, $seq);
	}

	my $t_name = ($f->{_tran_name}) ? $f->{_tran_name} : "$g_name-mRNA-$i"; #affects GFFV3.pm
	my $t_id = ($f->{_tran_id}) ? $f->{_tran_id} : $t_name; #affects GFFV3.pm

	#double check name (mRNA count may have changed)
	if($f->name ne $t_name){
	    $struct->{t_name} = $t_name;
	    $struct->{t_id}   = $t_id;
	    $f->name($t_name);
	}

	#add statistics to existing data structure
	$struct->{hit}  = $f;
	$struct->{t_qi} = $qi;
	$struct->{AED}  = $AED;
	$struct->{eAED} = $eAED;

	return $struct;
}
#------------------------------------------------------------------------
#takes an array of annotations and only returns those that overlap a maker
#annotation (doesn't care if on opposite strand so be careful)
sub get_genes_overlapping_gene {
   my $gene  = shift;
   my $genes = shift;

   my $B = $gene->{g_start};
   my $E = $gene->{g_end};

   ($B, $E) = ($E, $B) if($E < $B);

   my @keepers;
   foreach my $g (@{$genes}){
       my $gB = $g->{g_start};
       my $gE = $g->{g_end};

       ($gB, $gE) = ($gE, $gB) if($gE < $gB);
       my $comp = compare::compare ($B,
				    $E,
				    $gB,
				    $gE
				  );
       if($comp ne '0'){
	   push(@keepers, $g);
       }
   }

   return \@keepers;
}
#------------------------------------------------------------------------
#takes an array of phathits and only returns those that overlap a maker
#annotation (will not return hits on opposite strand)
sub get_hits_overlapping_gene {
   my $g  = shift;
   my $hits = shift;

   my $B = $g->{g_start};
   my $E = $g->{g_end};

   my @keepers;
   foreach my $hit (@{$hits}){
      next unless $hit->strand eq $g->{g_strand};

      my $comp = compare::compare ($B,
				   $E,
				   $hit->start,
				   $hit->end
				  );
      if($comp ne '0'){
	  push(@keepers, $hit);
      }
   }

   return \@keepers;
}
#------------------------------------------------------------------------
#takes two arrays of phathits and returns all members of the second
#array that do not overlap those of the first array
sub get_non_overlapping {
   my $hits_a = shift;
   my $hits_b = shift;

   my @keepers;
   foreach my $b (@{$hits_b}){
      my $hit_one = 0;
      foreach my $a (@{$hits_a}){
	 $hit_one++ if compare::overlap($a, $b, 'query', 5);
      }

      push(@keepers, $b) unless $hit_one;
   }
   return \@keepers;
}
#------------------------------------------------------------------------
#called by subroutine annotate.  Actual maker annotation structure built here.
#groups overlapping transcripts from given predictor type as isoforms of the
#same gene.  Gene name is set here. load_transcript_struct is called here
#to set up individual transcript structure. model_gff transcripts are grouped
#based on the gene_id given in the phat_hit structure and not by overlap.

sub group_transcripts {
   my $transcripts  = shift;
   my $data         = shift;
   my $seq          = shift;
   my $seq_id       = shift;
   my $chunk_number = shift;
   my $predictor    = shift;
   my $the_void     = shift;
   my $CTL_OPT      = shift;

   #fix weird sequence names
   $seq_id = Fasta::seqID2SafeID($seq_id);

   #place evidence and p_bases in index for easy retrieval
   my @transcripts;
   my %lookup;
   my %p_bases;
   my %snap_lookup;
   my $temp_id = 0;
   foreach my $s (@$transcripts){
       my $tra    = $s->[0];
       my $set    = $data->[$s->[1]];
       my $p_base = $s->[2];

       #may overlap more predictions than seen in original cluster
       #my $all_preds = get_overlapping_hits($tra, $predictions);
       #$set->{all_preds} = $all_preds;

       $tra->{set_id} = $temp_id;
       $lookup{$temp_id} = $set;
       $p_bases{$temp_id} = $p_base;
       push(@transcripts, $tra);
      $temp_id++;
   }

   #cluster the transcripts to get genes
   my $careful_clusters = [];

   if ($predictor =~ /^model_gff$|_abinit$|^pred_gff$/ ||
       ($predictor =~ /^est2genome$/ && $CTL_OPT->{est_forward})
       ) {
       my @to_do;
       my %index;
       my $i = @$careful_clusters;
       foreach my $t (@transcripts) {
	   my $j;
	   if($predictor =~ /^est2genome$/ && ! exists $t->{gene_id}){
	       push(@to_do, $t);
	       next;
	   }
	   elsif(! exists $t->{gene_id}){
	       $j = $i;
	       $i++;
	   }
	   elsif (exists $index{$t->{gene_id}}) {
	       $j = $index{$t->{gene_id}};
	   }
	   else {
	       $j = $i;
	       $index{$t->{gene_id}} = $j;
	       $i++;
	   }
	   push(@{$careful_clusters->[$j]}, $t);
       }
       @transcripts = @to_do;
   }

   #seperate out when multiple HMM's are provided (comma seperated list)
   my %sources;
   foreach my $t (@transcripts){
       push(@{$sources{$t->{_HMM}}}, $t);
       die "ERROR: No hit source {_HMM} in maker::auto_annotator\n" if(! $t->{_HMM});
   }
   #now cluster each list seperately
   foreach my $set (values %sources){
       my $clusters = cluster::careful_cluster_phat_hits($set, $seq);
       
       #remove redundant transcripts in gene
       foreach my $c (@{$clusters}) {
	   my $best_alt_forms = ($CTL_OPT->{est_forward}) ?
	       $c : clean::remove_redundant_alt_splices($c, $seq, 10);
	   push(@$careful_clusters, $best_alt_forms);
       }
   }

   #process clusters into genes
   my $c_id = 0;
   my @annotations;
   foreach my $c (@$careful_clusters) {
      #build gene name here
      my %pred_sources;
      foreach my $f (@{$c}) {
	 my $source = $f->algorithm;
	 $source =~ s/\:+/_/;
	 $pred_sources{$source}++;
      }
      my $sources = join ('-', keys %pred_sources);

      my $g_name;
      my $g_id;
      if ($predictor eq 'model_gff') {
	 $g_name = $c->[0]->{gene_name} || $c->[0]->{gene_id}; #affects GFFV3.pm
	 $g_id = $c->[0]->{gene_id} || $c->[0]->{gene_name}; #affects GFFV3.pm
	 $SEEN->{$g_name}++;
	 $SEEN->{$g_id}++;
	 if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
	     $SEEN->{$1}++;
	 }
	 if($g_id =~ /(\d+\.\d+)\-mRNA\-\d+/){
	     $SEEN->{$1}++;
	 }
      }
      elsif ($predictor =~ /_abinit$/) {
	  #now check for preexisting name
	  if ($c->[0]->{gene_name} || $c->[0]->{gene_id}){
	      $g_name = $c->[0]->{gene_name} || $c->[0]->{gene_id}; #affects GFFV3.pm
	      $g_id = $c->[0]->{gene_id} || $c->[0]->{gene_name}; #affects GFFV3.pm
	      $SEEN->{$g_name}++;
	      $SEEN->{$g_id}++;
	      if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
		  $SEEN->{$1}++;
	      }
	      if($g_id =~ /(\d+\.\d+)\-mRNA\-\d+/){
		  $SEEN->{$1}++;
	      }
	  }
	  elsif ($c->[0]->name =~ /^maker-$seq_id|$seq_id-abinit/) {
	      $g_name = $c->[0]->name;
	      $g_name =~ s/-mRNA-\d.*//;
	      $g_id = $g_name;
	      $SEEN->{$g_name}++;
	      if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
		  $SEEN->{$1}++;
	      }
	  }
	  else{
	      $g_name = "$sources-$seq_id-abinit-gene-$chunk_number"; #affects GFFV3.pm
	      $c_id++ while(exists $SEEN->{"$chunk_number\.$c_id"} || exists $SEEN->{"$g_name.$c_id"});
	      $g_name = "$g_name.$c_id";
	      $g_id = $g_name;
	      $SEEN->{$g_name}++;
	      $SEEN->{"$chunk_number\.$c_id"}++;
	  }
      }
      elsif ($predictor eq 'est2genome' && $CTL_OPT->{est_forward} && ($c->[0]->{gene_id} || $c->[0]->{gene_name})) {
	  $g_name = $c->[0]->{gene_name} || $c->[0]->{gene_id}; #affects GFFV3.pm
	  $g_id = $c->[0]->{gene_id} || $c->[0]->{gene_name}; #affects GFFV3.pm
	  $SEEN->{$g_name}++;
	  $SEEN->{$g_id}++;
	 if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
	     $SEEN->{$1}++;
	 }
	 if($g_id =~ /(\d+\.\d+)\-mRNA\-\d+/){
	     $SEEN->{$1}++;
	 }
      }
      else{
	 $g_name = "maker-$seq_id-$sources-gene-$chunk_number"; #affects GFFV3.pm
	 $c_id++ while(exists $SEEN->{"$chunk_number\.$c_id"} || exists $SEEN->{"$g_name.$c_id"});
	 $g_name = "$g_name.$c_id";
	 $g_id = $g_name;
	 $SEEN->{$g_name}++;
	 $SEEN->{"$chunk_number\.$c_id"}++;
      }

      #combine evidence for all transcripts
      my $evidence = {};
      foreach my $f (@{$c}) {
	  my $evi = defined($f->{set_id}) ? $lookup{$f->{set_id}} : {};
	  merge_evidence($evidence, $evi);
	  push(@{$evidence->{index}}, $evi->{index});
      }

      #load transcript structs
      my $i = 1;
      my @t_structs;
      foreach my $f (@{$c}) {
	  my $p_base = defined($f->{set_id}) ? $p_bases{$f->{set_id}} : undef;
		
	  my $t_struct = load_transcript_struct($f, $g_name, $i, $seq, $p_base, $predictor, $CTL_OPT);
		
	  push(@t_structs, $t_struct);
	  $i++;
      }

      #infer gene name for pass-through transcripts (no gene name but yes tran name)
      if($CTL_OPT->{est_forward} && !$t_structs[0]->{hit}->{gene_name} &&
	 !$t_structs[0]->{hit}->{gene_id} && $t_structs[0]->{hit}->{_tran_name}
	 ){
	  $g_name = $t_structs[0]->{t_name}."-gene";
	  $g_id = $t_structs[0]->{t_id}."-gene";
	  $SEEN->{$g_name};
      }

      my ($g_start, $g_end, $g_strand) = get_start_and_end_on_seq(\@t_structs);
      my $g_attrib = (exists $t_structs[0]->{hit}->{gene_attrib}) ? $t_structs[0]->{hit}->{gene_attrib} : undef;

      my $annotation = { 't_structs'  => \@t_structs, 
			 'g_name'     => $g_name,
			 'g_id'       => $g_id,
			 'g_start'    => $g_start,
			 'g_end'      => $g_end,
			 'g_strand'   => $g_strand,
			 'g_evidence' => $evidence,
			 'g_evi_index'=> $evidence->{index},
			 'predictor'  => $predictor,
			 'algorithm'  => $sources,
			 'g_attrib'   => $g_attrib
		       };

      push(@annotations, $annotation);
      $c_id++;
   }

   return \@annotations;
}
#------------------------------------------------------------------------
#merges evidence of data supporting annotations.  Merges them so that
#all members are uniq and not repeated.  This is important for combining
#evidence from mutliple transcripts as a single gene.

sub merge_evidence {
    my $evi1 = shift;
    my $evi2 = shift;

    while(my $key = each %$evi2){
	next if(ref($evi2->{$key}) ne 'ARRAY');

	#reset uniq for set to be added
	foreach my $f (@{$evi2->{$key}}){
	    $f->{_uniq_set} = 0;
	}

        #set uniq for existing set
	$evi1->{$key} = [] if(! $evi1->{$key});
	foreach my $f (@{$evi1->{$key}}){
	    $f->{_uniq_set} = 1;
	}

	#now only add hits where uniq is not set
	foreach my $f (@{$evi2->{$key}}){
	    push(@{$evi1->{$key}}, $f) if(! $f->{_uniq_set});
	    $f->{_uniq_set} = 1;
	}
    }
}
#------------------------------------------------------------------------
#attempts to map names from model_gff annotions forward onto new annotations.
#this is done using AED to decide who best represents the previous model,
#best decribed as reciprical best AED. called outside of package by maker.

sub map_forward {
    my $ann_set = shift;
    my $gff_set = shift;

    return $ann_set if(! $gff_set || ! @$gff_set);

    #make temp IDs and indexes to keep track with
    my $temp = 0;
    my @index;
    my @g_index;

    #separate annotations by strand
    my @p_ann;
    my @m_ann;
    foreach my $g (@$ann_set){
	my $hits = gene2all_phathits($g);
	if($g->{g_strand} == 1){
	    push(@p_ann, @{$hits});
	}
	elsif($g->{g_strand} == -1){
	    push(@m_ann, @{$hits});
	}
	else{
	    die "ERROR: Logic error in auto_annotator::verify_old_forms\n";
	}

	#add temp ID
	foreach my $hit (@$hits){
	    $hit->{_temp_id} = $temp;
	    $index[$temp] = $hit;
	    $g_index[$temp] = $g;
	    $temp++;
	}
    }

    #separate model_gffs hits by strand
    my @p_gff;
    my @m_gff;
    foreach my $g (@$gff_set){
	my $hits = gene2all_phathits($g);
	if($g->{g_strand} == 1){
	    push(@p_gff, @{$hits});
	}
	elsif($g->{g_strand} == -1){
	    push(@m_gff, @{$hits});
	}
	else{
	    die "ERROR: Logic error in auto_annotator::verify_old_forms\n";
	}

	#add temp ID
	foreach my $hit (@$hits){
	    $hit->{_temp_id} = $temp;
	    $index[$temp] = $hit;
	    $g_index[$temp] = $g;
	    $temp++;
	}
   }

    #identify closest forms on plus strand
    my @closest;
    my @AEDs;
    foreach my $ann (@p_ann){
	my $bag = get_overlapping_hits($ann, \@p_gff);
	my $a_id = $ann->{_temp_id};
	foreach my $gff (@$bag){
	    my $g_id = $gff->{_temp_id};
	    my $AED = shadow_AED::get_AED([$gff], $ann);
	    next if($AED > 0.5); #must be closer than 0.5
	    push(@{$ann->{_Alias}}, $gff->name);
	    push(@{$ann->{_Alias}}, $ann->name); 
	    push(@{$gff->{_Alias}}, $ann->name); 

	    if (! defined($AEDs[$a_id]) || $AEDs[$a_id] > $AED){
		$closest[$a_id] = $g_id;
		$AEDs[$a_id] = $AED;
	    }
	    if (! defined($AEDs[$g_id]) || $AEDs[$g_id] > $AED){
		$closest[$g_id] = $a_id;
		$AEDs[$g_id] = $AED;
	    }
	}
    }

    #identify closest forms on minus strand
    foreach my $ann (@m_ann){
	my $bag = get_overlapping_hits($ann, \@m_gff);
	my $a_id = $ann->{_temp_id};
	foreach my $gff (@$bag){
	    my $g_id = $gff->{_temp_id};
	    my $AED = shadow_AED::get_AED([$gff], $ann);
	    next if($AED > 0.5); #must be closer than 0.5
	    push(@{$ann->{_Alias}}, $gff->name); 
	    push(@{$ann->{_Alias}}, $ann->name); 
	    push(@{$gff->{_Alias}}, $ann->name); 
	    if (! defined($AEDs[$a_id]) || $AEDs[$a_id] > $AED){
		$closest[$a_id] = $g_id;
		$AEDs[$a_id] = $AED;
	    }
	    if (! defined($AEDs[$g_id]) || $AEDs[$g_id] > $AED){
		$closest[$g_id] = $a_id;
		$AEDs[$g_id] = $AED;
	    }
	}
    }

    #map names forward
    foreach my $g (@$ann_set){
	my $new; #new gene name
	my $nid; #new gene name
	my $AED = 1; #for comparing which isoform is closest
	my $fg;
	foreach my $t (@{$g->{t_structs}}){
	    my $id = $t->{hit}->{_temp_id}; #get id to see who to map forward
	    my $cl = $closest[$id]; #get closest model

	    #nothing to change so skip
	    next if(! defined($cl) || $closest[$cl] != $id || $cl == $id);

	    my $f = $index[$cl]; #closet gff hit

	    #decide which isoform to get gene name from
	    #this affects the gene not the transcript
	    if($AEDs[$id] < $AED){
		$AED = $AEDs[$id];
		$new = $f->{gene_name}; #get new gene name
		$nid = $f->{gene_id}; #get new gene name
		$fg = $g_index[$cl]; #get annotation for that gene
	    }
	    $t->{is_changed} = ($AEDs[$id] == 0) ? 0 : 1; #not a true change if identical
	    $t->{is_changed} = 1 if(defined $t->{-attrib}); #changed if new transcript has own attributes (pred_gff?)
	    $g->{is_changed} = 1 if($t->{is_changed}); #gene changed only if trans changed
	    $t->{t_name} = $f->{_tran_name}; #set transcript name
	    $t->{t_id} = $f->{_tran_id}; #set transcript id

	    #this helps with attribute passthrough
	    $t->{hit} = $f if(! $t->{is_changed}); #use gff hit if they are identical

	    #!!remember hit name has not been changed on non identical!! ;-)
	}

	next if(! $new); #no change so skip

	$g->{g_name} = $new; #set new gene name
	$g->{g_id}   = $nid; #set new gene name

	#see if gene changed by altered transcript content
	if(@{$g->{t_structs}} == @{$fg->{t_structs}} && ! $g->{is_changed}){
	    @{$g->{t_structs}} = sort {$a->{t_name} cmp $b->{t_name}} @{$g->{t_structs}};
	    @{$fg->{t_structs}} = sort {$a->{t_name} cmp $b->{t_name}} @{$fg->{t_structs}};

	    my $changed;
	    for (my $i = 0; $i < @{$g->{t_structs}}; $i++){
		my $ni = $g->{t_structs}->[$i]->{t_name};
		my $nj = $fg->{t_structs}->[$i]->{t_name};
		$changed = 1 if($ni ne $nj);
	    }

	    $g->{is_changed} = 1 if($changed);
	}
	else{
	    $g->{is_changed} = 1; #gene changed if trans count changed
	}

	#if there was no chage in gene just use the old gff gene
	#this helps with attribute passthrough
	$g = $fg if(! $g->{is_changed});;
    }

    return $ann_set;
}
#------------------------------------------------------------------------
#called outside of package by maker to identify abinit maker annotations
#that don't overlap maker's final annotation set. The non-overlapping
#set is filtered using abAED so that they don't overlap each other.
#Only masked ab inits are considered (otherwise I get a lot of revere
#transcriptase genes).  Unmasked are allowed when there is no masking
#performed.

sub get_non_overlaping_abinits {
   my $ann_set = shift;
   my $all_set = shift;
   my $CTL_OPT = shift;

   my @overlap;
   my @none;

   my $abin_set = [];
   my @ab_keys = grep {/_abinit$|^pred_gff/} keys %$all_set;

   foreach my $key (@ab_keys){
       push(@$abin_set, @{$all_set->{$key}});
   }

   #separate annotations by strand
   my @p_ann;
   my @m_ann;
   foreach my $g (@$ann_set){
       if($g->{g_strand} == 1){
	   push(@p_ann, $g);
       }
       elsif($g->{g_strand} == -1){
	   push(@m_ann, $g);
       }
       else{
	   die "ERROR: Logic error in auto_annotator::best_annotations\n";
       }
   }

   #separate abinits by strand
   my @p_ab;
   my @m_ab;
   foreach my $g (@$abin_set){
       #only accept masked predictions unless I'm not masking or the predictor is genemark
       my $src = $g->{algorithm};
       unless($src =~ /_masked$|^pred_gff/ || $CTL_OPT->{_no_mask} || $CTL_OPT->{predictor} eq 'genemark'){
	   next;
       }

       if($g->{g_strand} == 1){
	   push(@p_ab, $g);
       }
       elsif($g->{g_strand} == -1){
	   push(@m_ab, $g);
       }
       else{
	   die "ERROR: Logic error in auto_annotator::best_annotations\n";
       }
   }

   #identify non-overlapping in plus strand
   my @p_keepers;
   foreach my $ab (@p_ab){
      my $bad = 0;
      foreach my $ann (@p_ann){
	 my $comp = compare::compare ($ab->{g_start},
				      $ab->{g_end},
				      $ann->{g_start},
				      $ann->{g_end}
				     );
	 if($comp ne '0'){
	    $bad = 1;
	    last;
	 }
      }

      if(! $bad){
	 push (@p_keepers, $ab);
      }
   }

   #identify non-overlapping in minus strand
   my @m_keepers;
   foreach my $ab (@m_ab){
      my $bad = 0;
      foreach my $ann (@m_ann){
	 my $comp = compare::compare ($ab->{g_start},
				      $ab->{g_end},
				      $ann->{g_start},
				      $ann->{g_end}
				     );
	 if($comp ne '0'){
	    $bad = 1;
	    last;
	 }
      }

      if(! $bad){
	 push (@m_keepers, $ab);
      }
   }

   #get best non-overlapping set
   my @keepers;
   @p_keepers  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @p_keepers;
   @m_keepers  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @m_keepers;
   push(@keepers, @{_best(\@p_keepers)});
   push(@keepers, @{_best(\@m_keepers)});

   return (\@keepers);
}
#------------------------------------------------------------------------
#called by group_transcripts to get gene start and gene end
sub get_start_and_end_on_seq {
    my $transcripts = shift;
    
    my @exons;
    foreach my $t (@{$transcripts}){
	my $phat_hit = $t->{hit};
	
	foreach my $hsp ($phat_hit->hsps){
	    push(@exons, $hsp);
	}
    }

    if(! @exons){
	return;
    }

    my @sorted_b = sort {$a->start('query') <=> $b->start('query')} @exons;
    my @sorted_e = sort {$b->end('query')   <=> $a->end('query')}   @exons;
    
    my $ref = ref($sorted_b[0]);
    
    my $g_start = $sorted_b[0]->start('query');
    my $g_end   = $sorted_e[0]->end('query');
    
    my $strand  = $sorted_b[0]->strand('query');
    
    return ($g_start, $g_end, $strand);    
}
#------------------------------------------------------------------------
#takes an array of phathits and returns only those of a given algorithm
#i.e. used to separate est2genome or blastx hits from a mixed array
sub get_selected_types {
	my $c = shift;

	my @keepers;
	foreach my $f (@{$c}){
	    foreach my $type (@_){
		if ($f->algorithm =~ /^$type$/i){
		    push(@keepers, $f);
		    last;
		}
		elsif ($f->algorithm =~ /^$type\:/i){
		    push(@keepers, $f);
		    last;
		}
		elsif ($f->algorithm =~ /^$type\_masked$/i){
		    push(@keepers, $f);
		    last;
		}
		elsif ($f->algorithm =~ /exonerate\:\:$type$/i){
		    push(@keepers, $f);
		    last;
		}
	    }
	}

	return \@keepers;
}
#------------------------------------------------------------------------
#gets the transcript sequence for a phathit using the query sequence
sub get_transcript_seq {
	my $anno = shift;
	my $seq  = shift;

	my $sorted = PhatHit_utils::sort_hits($anno, 'query');

	if(!$seq){
	    return join('', map {$_->hit_string()} @{$sorted});
	}

	my $transcript = '';

	foreach my $hsp (@{$sorted}){
		my $e_b = $hsp->nB('query');
		my $e_e = $hsp->nE('query');

		my $length = abs($e_e - $e_b) + 1;

	       ($e_b, $e_e) = ($e_e, $e_b)
		if $hsp->strand('query') == -1;

		my $exon_seq = substr($$seq, $e_b - 1, $length);

		$exon_seq = Fasta::revComp($exon_seq)
		if $hsp->strand('query') == -1;

		$transcript .= $exon_seq;
	}

	return $transcript;
}
#------------------------------------------------------------------------
#finds longest translatable sequence beginning with Methionine.
#returns sequence and offset. Called by get_translation_seq
sub get_longest_m_seq {
	my $seq = shift;

	my ($off_0, $p_seq_0) = get_off_and_str($seq, 0);
	my ($off_1, $p_seq_1) = get_off_and_str($seq, 1);
	my ($off_2, $p_seq_2) = get_off_and_str($seq, 2);

	my @data;
	push(@data, [$p_seq_0, $off_0]) if defined($p_seq_0);
	push(@data, [$p_seq_1, $off_1]) if defined($p_seq_1);
	push(@data, [$p_seq_2, $off_2]) if defined($p_seq_2);

	my @sorted = sort {length($b->[0]) <=> length($a->[0])} @data;

	my $best = shift(@sorted);

	return ($best->[1], $best->[0]) if($best);
}
#------------------------------------------------------------------------
#takes a hit and an array of phathits.  Returns all hits in the array
#overlapping the first hit. Strand is tested before returning.
sub get_overlapping_hits {
    my $eat  = shift;
    my $hits = shift;

    my @keepers;
    foreach my $hit (@{$hits}){
	next unless $hit->strand('query') eq $eat->strand('query');
	push(@keepers, $hit)
	    if compare::overlap($hit, $eat, 'query', 3);
    }
    return \@keepers;
}
#------------------------------------------------------------------------
#returns longest translation starting with M and the offset for a
#given reading frame. called by get_longest_m_seq
sub get_off_and_str {
	my $seq  = shift;
	my $offset = shift;

	my $tM = new CGL::TranslationMachine();

	my $p_seq = $tM->translate_from_offset($seq, $offset);

	return (undef, undef) unless $p_seq =~ /M/;

	my $best_pos;
	my $best_len;
	my $pos = -1;
	while(($pos = index($p_seq, 'M', $pos)) > -1){
		my ($open_run) = substr($p_seq, $pos) =~ /(^[^\*]+)\*?/;
		my  $length    = length($open_run);
		if (!defined($best_len) || $length > $best_len){
			$best_len = $length;
			$best_pos = $pos;
		}

		$pos++;
	}

	$offset = $offset + 3*$best_pos;

	my $p_seq_2 = $tM->translate_from_offset($seq, $offset);

	my ($t_seq) = $p_seq_2 =~ /(^[^\*]+\*?)/;

	die "logic error in auto_annotate::get_off_and_str!\n"
	unless $t_seq  =~ /^M/;

	return ($offset, $t_seq);
}
#------------------------------------------------------------------------
#returns what maker believes to be the best translation seq and offset
sub get_translation_seq {
    my $seq    = shift;
    my $f      = shift;
    
    my $tM = new CGL::TranslationMachine();
    my $p_seq;
    my $offset;

    #use offset and end already in model to guide seq selection
    if(defined($f->{translation_offset}) && $f->{translation_end}){
	$offset = $f->{translation_offset};
	my $has_start = $tM->is_start_codon(substr($seq, $offset, 3));

	#step upstream to find longer ORF
	for(my $i = $offset - 3; $i >= 0; $i -= 3){
	    my $codon = substr($seq, $i, 3);
	    last if($tM->is_ter_codon($codon));
	    
	    if(my $s = $tM->is_start_codon($codon) || !$has_start){
		$offset = $i;
		$has_start = $s;
	    }
	}

	#translate to stop
        $p_seq = $tM->translate_from_offset($seq, $offset);
        $p_seq =~ s/(^[^\*]+\*?).*/$1/;
    }
    
    #build a new translation
    if(! defined($offset) || ($p_seq !~ /^M/ && $offset >= 3)){
	($p_seq , $offset) = $tM->longest_translation_plus_stop($seq);

	#does not begin with M....
	if ($p_seq !~ /^M/){
	    my ($off_new, $p_seq_new) = get_longest_m_seq($seq);

	    #take M start sequence in most cases
            unless(!$p_seq_new || ($offset < 3 && $off_new - $offset > 90)){
		$offset = $off_new;
		$p_seq = $p_seq_new;
            }
	}
    }

    #fix for translations longer than seq (because of third codon ambiguity)
    if(length($p_seq)*3 + $offset > length($seq)){
	$p_seq =~ s/.$//; #remove trailing peptide
    }

    #get end, and see if there is a stop, and remove it
    my $end = length($p_seq)*3 + $offset + 1;
    my $has_start = 1 if($p_seq =~ /^M/);
    my $has_stop = 1 if($p_seq =~ s/\*$//);

    #if very small CDS try again with longest ORF rather than internal offest
    if(($f->{translation_offset} || $f->{translation_end}) && ($end - $offset - 1)/(length($p_seq)*3) < 40){
	$f->{translation_offset} = undef;
	$f->{translation_end} = undef;
	return get_translation_seq($seq, $f);
    }

    #set CDS internally in hit (CDS includes stop)
    $f->{translation_offset} = $offset;
    $f->{translation_end} = $end;

    #find correct spacial coordinates for easy access
    my $coorB;
    my $coorE;
    my ($toffset, $tend) = ($offset, $end);
    foreach my $hsp (@{PhatHit_utils::sort_hits($f, 'query')}){
	my $l = abs($hsp->nE('query') - $hsp->nB('query')) + 1;
        #find first bp coordinate
	if($toffset && $l <= $toffset){
	    $toffset -= $l;
	    $tend -= $l;
	    next;
	}
	elsif(!$coorB){
	    $coorB = ($hsp->strand == 1) ? $hsp->nB('query') + $toffset : $hsp->nB('query') - $toffset;
	    $tend -= $l;
	    $toffset = 0;
	}
	else{
	    $tend -= $l;
	}

        #find last bp coordinate
	if($tend <= 1){ #end is always bp after last translated bp
	    $coorE = ($hsp->strand == 1) ? $hsp->nE('query') + ($tend - 1) : $hsp->nE('query') - ($tend - 1);
	    last;
	}
    }

    die "ERROR: Logic problem in maker::auto_annotator::get_translation_seq\n" if(!defined($coorE));

    $f->{_TSTART}{query} = $coorB;
    $f->{_TSTART}{hit} = $offset + 1;
    $f->{_TEND}{query} = $coorE;
    $f->{_TEND}{hit} = $end - 1;

    #return
    return ($p_seq , $offset, $end, $has_start, $has_stop);
}
#------------------------------------------------------------------------
#takes the gene predictions produce by get_pred_shot and finds the one
#with the higest score.  Called by run_it
#sub get_best_pred_shot {
#   my $wanted_strand = shift;
#   my $gene_preds    = shift;
#
#   my @gs;
#   foreach my $g (@{$gene_preds}){
#      next unless defined($g);
#      next unless defined($g->strand('query'));
#      next unless $g->strand('query') == $wanted_strand;
#      my $total_score = PhatHit_utils::get_total_score_of_hit($g);
#      push(@gs, [$total_score, $g]);
#   }
#   my @sorted = sort {$b->[0] <=> $a->[0]} @gs;
#
#   my $best  = $sorted[0]->[1];
#
#   return $best;
#}
#------------------------------------------------------------------------
#returns all pred shots on the right strand, so best_pred_shots might be
#a misnomer.  Called by run_it
sub get_best_pred_shots {
   my $wanted_strand = shift;
   my $gene_preds    = shift;

   my @gs;
   foreach my $g (@{$gene_preds}){
      next unless defined($g);
      next unless defined($g->strand('query'));
      next unless $g->strand('query') == $wanted_strand;
      push(@gs, $g);
   }

   return \@gs;
}
#------------------------------------------------------------------------
#takes a transcript phathit and a set of ESTs then tries and determine
#UTRs from those.  The new transcript plus UTR transcript is returned.
sub pneu {
	my $ests       = shift;
	my $g          = shift;
	my $q_seq      = shift;

	my @transcripts;

	my $num_hsps = $g->hsps;

	my $b_5;
	my $b_3;

	my @hsps = $g->hsps;

	if ($num_hsps == 1){
		($b_5, $b_3)  = maker::join::find_best_one($g, $ests);
	}
	else {
	   $b_5 = maker::join::find_best_five($g, $ests);
	   $b_3 = maker::join::find_best_three($g, $ests);
	}

	my $anno_transcript =
	maker::join::join_f($b_5, $g, $b_3, $q_seq);

	return $anno_transcript;
}
#------------------------------------------------------------------------
#combines two phathit arrays into one
sub combine {
	my @bag;
	while (my $hits = shift(@_)){
		foreach my $hit (@{$hits}){
			push(@bag, $hit);
		}
	}
	return \@bag;
}
#------------------------------------------------------------------------
1;
