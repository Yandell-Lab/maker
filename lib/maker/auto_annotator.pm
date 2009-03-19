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
#use evaluator::so_classifier;
#use evaluator::evaluate;
#use evaluator::funs;
#use evaluator::AED;
use shadow_AED;

$Storable::forgive_me = 1; 

@ISA = qw(
       );

my $LOG; #GLOBAL VARIABLE
my $SEEN;

#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub prep_hits {
	my $prot_hits        = shift;
	my $est_hits         = shift;
	my $alt_est_hits     = shift;
	my $predictions      = shift;
	my $models           = shift;
	my $seq              = shift;
	my $single_exon      = shift;

	#only ESTs from splice concious algorithms i.e. no blastn
	my $clean_est = get_selected_types($est_hits,'est2genome', 'est_gff');
	my $clean_altest = $alt_est_hits;
	if ($single_exon == 1) {
	    #do nothing
	}else {
	    # don't use unpliced single exon ESTs-- usually genomic contamination
	    $clean_est = clean::purge_single_exon_hits($clean_est);
	    $clean_altest = clean::purge_single_exon_hits($clean_altest);
	}

	# throw out the exonerate est hits with weird splice sites
	$clean_est = clean::throw_out_bad_splicers($clean_est, $seq);
	
	# combine puts type in order they are given
	# important for later culstering, as hits in
	# first arg more likey to make in into to cluster
	# than second arg, etc
	
	my $c_bag = combine($models,
			    $prot_hits,
	                    $clean_est, 
			    $clean_altest
			   );

	#--- c_bag processing for gene prediction and gff3 models
        my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $c_bag);
        my $p_clusters = cluster::shadow_cluster(0, $seq, $p, 10);
        my $m_clusters = cluster::shadow_cluster(0, $seq, $m, 10);

	#purge after clustering so as to still have the effect of evidence joining ESTs
	$p_clusters = purge_short_ESTs_in_clusters($p_clusters, 250);
	$m_clusters = purge_short_ESTs_in_clusters($m_clusters, 250);

	my $careful_clusters = [];
	push(@{$careful_clusters}, @{$p_clusters});
	push(@{$careful_clusters}, @{$m_clusters});

	# identify the ab-inits that fall within and between clusters
	my ($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions, 
	                                                              $careful_clusters,
	                                                              $seq,
								     );
	foreach my $h (@{$hit_mult}){
	    $h->{_hit_multi} = 1;
	}

	#--new clusters built by merging clusters that overlap preds
	#--used for ab-initio coring
	#--must be ran before merge as merge alters cluster content
	my $pred_clusters = join_clusters_on_pred($predictions, $careful_clusters, $c_index);

	#--split preds across clusters that they overlap
	merge_into_cluster($hit_one, $careful_clusters, $c_index);
	merge_into_cluster($hit_mult, $careful_clusters, $c_index);

	#--data prep
	my $c_id = 0;
	my @bx_data;
	my @gf_data;
	foreach my $c (@{$careful_clusters}){
	   my $bx = prep_blastx_data($c, $c_id, $seq);
	   push(@bx_data, @{$bx}) if defined $bx;
	   
	   if(@{$models}){
	       my $gf = prep_gff_data($c, $c_id, $seq);
	       push(@gf_data, @{$gf}) if defined $gf;
	   }

	   $c_id++;
	}

	#--- processing for scoring ab-initio predictions
	my @pr_data;
	foreach my $c (@{$pred_clusters}){
	    my $pr = prep_pred_data($c, $c_id, $seq);
	    push(@pr_data, @{$pr}) if defined $pr;
	    
	    $c_id++;
	}

	return (\@bx_data, \@gf_data, \@pr_data);
}
#------------------------------------------------------------------------
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
	    }
	    elsif ($count == 1){
		push(@hit_one, $index[$i]);
	    }
	    else {
		push(@hit_mult, $index[$i]);
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
	    }
	    elsif ($count == 1){
		push(@hit_one, $index[$i]);
	    }
	    else {
		push(@hit_mult, $index[$i]);
	    }
	}

        return (\%clusters_hit, \@hit_one, \@hit_none, \@hit_mult);
}
#------------------------------------------------------------------------
sub merge_into_cluster {
   my $hits = shift;
   my $clusters = shift;
   my $c_index = shift;

   foreach my $s (@{$hits}, ){
      my @keys = keys %{$c_index->{$s->{temp_id}}};
      foreach my $i (@keys){
	 push(@{$clusters->[$i]}, $s);
      }
   }
}
#------------------------------------------------------------------------
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
sub purge_short_ESTs_in_clusters{
    my $clusters = shift;
    my $min = shift;
    
    my @c_keepers;
    foreach my $c (@$clusters){
        my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
        my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
        my $models_in_cluster = get_selected_types($c,'model_gff', 'maker');
        my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh', 'twinscan', 'pred_gff');

	$ests_in_cluster = clean::purge_short_single_exons($ests_in_cluster, $min);
	$alt_ests_in_cluster = clean::purge_short_single_exons($alt_ests_in_cluster, $min);

	my @new_c;
	push(@new_c, @$ests_in_cluster);
	push(@new_c, @$alt_ests_in_cluster);
	push(@new_c, @$ps_in_cluster);
	push(@new_c, @$bx_in_cluster);
	push(@new_c, @$models_in_cluster);

	if(@new_c){
	    push(@new_c, @$preds_in_cluster);
	    push(@c_keepers, \@new_c);
	}
    }

    return \@c_keepers;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => each est from best ests
#snap predictions included
sub prep_blastx_data {
        my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;

        my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
        my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
        my $models_in_cluster = get_selected_types($c,'model_gff', 'maker');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh', 'twinscan', 'pred_gff');

	my @single = grep {! $_->{_hit_multi}} @{$preds_in_cluster};

        # groups of most informative protein hits
	# go ahead and inclde the proteion2genome data as well... why not?
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

        # group of most informative alt splices
        my $gomias = clean::purge_single_exon_hits($ests_in_cluster);
	$gomias = clean::get_best_alt_splices($gomias, $seq, 10);

	my @data;
        if (defined($gomias->[0])){
	   foreach my $mia (@{$gomias}){
	      push(@data, {'gomiph'    => $gomiph,
			   'preds'     => \@single,
			   'all_preds' => $preds_in_cluster,
			   'ests'      => $ests_in_cluster,
			   'alt_ests'  => $alt_ests_in_cluster,
			   'mia'       => $mia,
			   'model'     => undef,
			   'gomod'     => $models_in_cluster,
			   'c_id'      => $c_id
			  }
		  );
	   }
        }
        else {
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => \@single,
			'all_preds' => $preds_in_cluster,
			'ests'      => $ests_in_cluster,
			'alt_ests'  => $alt_ests_in_cluster,
			'mia'       => undef,
			'model'     => undef,
			'gomod'     => $models_in_cluster,
			'c_id'      => $c_id
		       }
	       );
	}

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => none
#predictions included

sub prep_gff_data {
        my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;

        my $models_in_cluster = get_selected_types($c,'model_gff', 'maker');
        my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
        my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh', 'twinscan', 'pred_gff');

	my @single = grep {! $_->{_hit_multi}} @{$preds_in_cluster};

        # groups of most informative protein hits
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;

	foreach my $model (@{$models_in_cluster}){
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => \@single,
			'all_preds' => $preds_in_cluster,
			'ests'      => $ests_in_cluster,
			'alt_ests'  => $alt_ests_in_cluster,
			'mia'       => undef,
			'model'     => $model,
			'gomod'     => undef,
			'c_id'      => $c_id
		       }
	       );
	}

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#ests => set of all ests
#protein homology =>  set of combined protein exonerate and blastx data
#alternative splice form => none
#predictions included

sub prep_pred_data {
        my $c    = shift;
	my $c_id = shift;
	my $seq  = shift;

        my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
        my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh', 'twinscan', 'pred_gff');

        # groups of most informative protein hits
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;

	foreach my $pred (@{$preds_in_cluster}){
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => $preds_in_cluster,
			'all_preds' => $preds_in_cluster,
			'model'     => $pred,
			'gomod'     => undef,
			'ests'      => $ests_in_cluster,
			'alt_ests'  => $alt_ests_in_cluster,
			'mia'       => undef,
			'c_id'      => $c_id
		       }
	       );
	}

	return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#ests => set of best ests from all ests
#protein homology => each protein exonerate structure
#alternative splice forms => each protein exonerate structure
#no snap predictions included
sub prep_polpro_data {
        my $c    = shift;
        my $c_id = shift;
        my $seq  = shift;

        my $ests_in_cluster = get_selected_types($c, 'est2genome', 'est_gff');
        my $ps_in_cluster   = get_selected_types($c,'protein2genome');

	my $possible_ext_sources = combine($ests_in_cluster, $ps_in_cluster);

	my $best_exts = clean::get_best_alt_splices($possible_ext_sources, $seq, 10);

        # group of most informative alt splices
        my $gomias = clean::get_best_alt_splices($ps_in_cluster, $seq, 10);

        my @data;
        foreach my $mia (@{$gomias}){
        	push(@data, {'gomiph' => [$mia],
                             'ests'   => $best_exts,
                             'mia'    => $mia,
                             'c_id'   => $c_id});
        }

        return \@data;
}
#------------------------------------------------------------------------
#returns an array of hashes with the following atributes
#ests => each best est from all ests
#protein homology =>  best set from combined protein exonerate and blastx data
#alternative splice forms => based on each best ests from all ests
#no snap predictions included
sub prep_polest_data {
        my $c    = shift;
        my $c_id = shift;
        my $seq  = shift;


        my $ests_in_cluster = get_selected_types($c, 'est2genome', 'est_gff');
        my $ps_in_cluster   = get_selected_types($c,'protein2genome');
        my $bx_in_cluster   = get_selected_types($c,'blastx', 'protein_gff');

        my $i_set      = combine($ps_in_cluster, $bx_in_cluster);
	my $best_p_set = clean::remove_redundant_alt_splices($i_set, $seq, 10); 

        my $best_exts  = clean::get_best_alt_splices($ests_in_cluster, $seq, 10);

        # group of most informative alt splices
        my $gomias = clean::get_best_alt_splices($ps_in_cluster, $seq, 10);

        my @data;
	foreach my $alt_splice_est (@{$best_exts}){
               	push(@data, {'gomiph' => $best_p_set,
                       	     'ests'   => [$alt_splice_est],
                             'mia'    => $alt_splice_est,
                             'c_id'   => $c_id});
	}
        return \@data;
}
#------------------------------------------------------------------------
#this subrutine returns finished annoations for all predictors
sub annotate {
    my $virgin_fasta     = shift;
    my $masked_fasta     = shift;
    my $chunk_number     = shift; #required to name genes for each chunk
    my $prot_evidence    = shift;
    my $est_evidence     = shift;
    my $alt_est_evidence = shift;
    my $predictions      = shift;
    my $models           = shift;
    my $the_void         = shift;
    my $build            = shift;
    my $CTL_OPTIONS      = shift;
    $LOG                 = shift;
    
    #process fasta files
    my $def   = Fasta::getDef($masked_fasta);
    my $seq_id  = Fasta::def2SeqID($def);
    my $seq_ref   = Fasta::getSeqRef($masked_fasta);
    my $v_seq_ref = Fasta::getSeqRef($virgin_fasta);
    
    #reset gene names
    my $GFF_DB = new GFFDB($CTL_OPTIONS);
    $SEEN = $GFF_DB->get_existing_gene_names($seq_id);
    
    #group evidence and predictions
    my ($bx_data, $gf_data, $pr_data) = prep_hits($prot_evidence,
						  $est_evidence,
						  $alt_est_evidence,
						  $predictions,
						  $models,
						  $v_seq_ref,
						  $CTL_OPTIONS->{single_exon}
						  );
    
    my %annotations;
    
    #---model passthrough here
    if(grep {/^model_gff$/} @{$CTL_OPTIONS->{_predictor}}){
	print STDERR "Processing GFF3 passthrough annotations\n" unless($main::quiet);
	my $model_trans = run_it($gf_data,
				 $the_void,
				 $seq_ref,
				 $v_seq_ref,
				 $def,
				 'model_gff',
				 $the_void,
				 $CTL_OPTIONS
				 );
	
	$annotations{'model_gff'} = group_transcripts($model_trans,
						      $v_seq_ref,
						      $seq_id,
						      $chunk_number,
						      $build,
						      'model_gff',
						      $the_void,
						      $CTL_OPTIONS
						      );
    }
    
    #---hint based gene prediction here
    foreach my $prdr (@{$CTL_OPTIONS->{_predictor}}){
	next if($prdr eq 'model_gff' || $prdr eq 'abinit');
	print STDERR "Producing $prdr hint based annotations\n" unless($main::quiet);
	
	my $transcripts = run_it($bx_data,
				 $the_void,
				 $seq_ref,
				 $v_seq_ref,
				 $def,
				 $prdr,
				 $CTL_OPTIONS
				 );
	
	my $annot = group_transcripts($transcripts,
				      $v_seq_ref,
				      $seq_id,
				      $chunk_number,
				      $build,
				      $prdr,
				      $the_void,
				      $CTL_OPTIONS
				      );
	
	   $annotations{$prdr} =  $annot;
    }
    
    #---abinit scoring here
    print STDERR "Processing ab-initio predictions\n" if(@$pr_data && ! $main::quiet);
    
    my $pred_trans = run_it($pr_data,
			    $the_void,
			    $seq_ref,
			    $v_seq_ref,
			    $def,
			    'abinit',
			    $CTL_OPTIONS
			    );
    
    my $all_ab = group_transcripts($pred_trans,
				   $v_seq_ref,
				   $seq_id,
				   $chunk_number,
				   $build,
				   'abinit',
				   $the_void,
				   $CTL_OPTIONS
				   );
    
    $annotations{'abinit'} = $all_ab;
    
    add_abAED(\%annotations);
    
    return \%annotations;
}
#------------------------------------------------------------------------
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
sub gene2all_phathits {
    my $g = shift;
    
    my @hits;
    foreach my $t (@{$g->{t_structs}}){
	push (@hits, $t->{hit});
    }

    return \@hits;
}
#------------------------------------------------------------------------
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
		push(@p_genes, $g);
		next if($p eq 'est2genome');
		my $hits = gene2allPbases($g);
		push(@p_bag, @$hits);
	    }
	    elsif($g->{g_strand} == -1) {
		push(@m_genes, $g);
		next if($p eq 'est2genome');
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
sub best_annotations {
    my $annotations = shift;
    my $out_base = shift;
    my $CTL_OPTIONS = shift;

    print STDERR "Choosing best annotations\n" unless($main::quiet);
    
    my @keepers;
    if(@{$CTL_OPTIONS->{_predictor}} > 1 || grep {/^abinit$/} @{$CTL_OPTIONS->{_predictor}}){
	#set up lists for plus and minus strands as well as possible mergers
	#predictor types are processed in the order given by control files
	my $p_list = [];
	my $m_list = [];
	my @p_est2g;
	my @m_est2g;
	foreach my $p (@{$CTL_OPTIONS->{_predictor}}){
	    foreach my $g (@{$annotations->{$p}}){
		if($p ne 'est2genome' && $g->{g_strand} == 1){
		    push(@$p_list, $g) if($g->{AED} < 1 || $p eq 'model_gff');
		}
		elsif($p ne 'est2genome' && $g->{g_strand} == -1) {
		    push(@$m_list, $g) if($g->{AED} < 1 || $p eq 'model_gff');
		}
		elsif($g->{g_strand} == 1){
		    push(@p_est2g, $g); #seperate est2genome genes
		}
		elsif($g->{g_strand} == -1){
		    push(@m_est2g, $g); #seperate est2genome genes
		}
		else{
		    die "ERROR: Logic error in auto_annotator::best_annotations\n";
		}
	    }
	}

	#remove low scoring overlaping genes
	@$p_list  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b)} @$p_list;
	@$m_list  = sort {crit1($a) <=> crit1($b) || crit2($a) <=> crit2($b) || crit3($a) <=> crit3($b)} @$m_list;
	push(@$p_list, @p_est2g); #est2genome added to end, will only appear if nothing else overlaps
	push(@$m_list, @m_est2g); #est2genome added to end, will only appear if nothing else overlaps 
	$p_list = _best($p_list);
	$m_list = _best($m_list);
	
	#final set
	push(@keepers, @$p_list);
	push(@keepers, @$m_list);
    }
    elsif(@{$CTL_OPTIONS->{_predictor}}){
	my $key = $CTL_OPTIONS->{_predictor}->[0];
	foreach my $g (@{$annotations->{$key}}){
	    push(@keepers, $g) if($g->{AED} < 1  || $key eq 'model_gff');
	}
    }

    return \@keepers;
}
#------------------------------------------------------------------------
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
   
   return ($g->{AED} + $g->{abAED})/2;
}
#------------------------------------------------------------------------
#sort by evidence AED score
sub crit2 {
    my $g = shift;

    return $g->{AED};
}
#------------------------------------------------------------------------
#sort by abinit AED score
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
sub run_it {
    my $data         = shift;
    my $the_void     = shift;
    my $seq          = shift;
    my $v_seq        = shift;
    my $def          = shift;
    my $predictor    = shift;
    my $CTL_OPTIONS  = shift;
    
    my $q_id = Fasta::def2SeqID($def); 
    my @transcripts;
    my $i = 0;
    foreach my $set (@{$data}) {
	my $mia      = $set->{mia};
	my $ests     = $set->{ests};
	my $model    = $set->{model};

	#------gff passthrough
	if ($predictor eq 'model_gff') {
	    next if(! defined $model);
	    my $transcript = $model;
	    	    
	    push(@transcripts, [$transcript, $set, undef]);
	    
	    $i++;
	    next;
	}

	#------ab-init passthrough
	if ($predictor eq 'abinit') {
	    next if(! defined $model);
	    my $transcript = $model;
	    
	    #add UTR to ab-inits
	    my $utr_trans = pneu($ests, $transcript, $seq);
	    push(@transcripts, [$utr_trans, $set, $transcript]);
	    
	    $i++;
	    next;
	}
	
	#------est2genome
	if ($predictor eq 'est2genome') {
	    next if (! defined $mia);
	    my $transcript = pneu($ests, $mia, $seq);
	    push(@transcripts, [$transcript, $set, $mia]) if $transcript;
	    $i++;
	    next;
	}
	
	#------default hint based behavior
	my $gomiph   = $set->{gomiph};
	my $blastx   = get_selected_types($gomiph,'blastx', 'protein_gff');
	my $pol_p    = get_selected_types($gomiph,'protein2genome');
	my $alt_ests = $set->{alt_ests};
	my $abinits  = $set->{all_preds};
	
	my ($pred_shots, $strand) = get_pred_shot($seq, 
						  $def, 
						  $i, 
						  $the_void, 
						  $set,
						  $predictor,
						  $CTL_OPTIONS,
						  $LOG
						  );
	
	my $on_right_strand = get_best_pred_shots($strand, $pred_shots);
	
	#added 2/23/2009 to reduce spurious gene predictions with only blastx
	if(! @$abinits && ! defined $mia && ! @$pol_p && @$on_right_strand > 1){
	    my $clean  = clean::purge_single_exon_hits($alt_ests);
	    my $coors  = PhatHit_utils::get_hsp_coors($blastx, 'query');
	    my $pieces = Shadower::getPieces($seq, $coors, 0);
	    
	    if(! @$clean && @$pieces <= 1){
		next; #skip these spurious predictions
	    }
	}
	
	#add transcripts
	foreach my $pred_shot (@{$on_right_strand}) {
	    my $copy = $pred_shot;
	    if (defined($pred_shot) && defined($mia)) {
		my $transcript = pneu($ests, $copy, $seq);
		if (defined($transcript)) {
		    push(@transcripts, [$transcript, $set, $pred_shot]);
		}
		else {
		    push(@transcripts, [$copy, $set, $pred_shot]);
		}
	    }
	    elsif (defined($pred_shot) && !defined($mia)) {
		push(@transcripts, [$copy, $set, $pred_shot]);
	    }
	}
	$i++;
    }
    
    return \@transcripts;
}
#------------------------------------------------------------------------
sub get_pred_shot {
   my $seq         = shift; 
   my $def         = shift;
   my $set_id      = shift;
   my $the_void    = shift;
   my $set         = shift;
   my $predictor   = shift;
   my $CTL_OPTIONS = shift;
   my $LOG         = shift;

   if($predictor eq 'snap'){
      my $pred_command = $CTL_OPTIONS->{snap}.' '.$CTL_OPTIONS->{snaphmm};
      return Widget::snap::get_pred_shot($seq, 
					 $def, 
					 $set_id, 
					 $the_void, 
					 $set, 
					 $CTL_OPTIONS->{pred_flank}, 
					 $pred_command,
					 $CTL_OPTIONS->{force},
					 $LOG
					);
   }
   elsif($predictor eq 'augustus'){
      my $pred_command = $CTL_OPTIONS->{augustus}.' --species='.$CTL_OPTIONS->{augustus_species};
      return Widget::augustus::get_pred_shot($seq, 
					     $def, 
					     $set_id, 
					     $the_void, 
					     $set, 
					     $CTL_OPTIONS->{pred_flank}, 
					     $pred_command,
					     $CTL_OPTIONS->{force},
					     $LOG
					    );
   }
   elsif($predictor eq 'fgenesh'){
      my $pred_command = $CTL_OPTIONS->{fgenesh}.' '.$CTL_OPTIONS->{fgenesh_par_file};
      return Widget::fgenesh::get_pred_shot($seq, 
					    $def, 
					    $set_id, 
					    $the_void, 
					    $set, 
					    $CTL_OPTIONS->{pred_flank}, 
					    $pred_command,
					    $CTL_OPTIONS->{force},
					    $LOG
					   );
   }
   elsif($predictor eq 'twinscan'){
      my $pred_command = $CTL_OPTIONS->{twinscan};
      return Widget::twinscan::get_pred_shot($seq, 
					     $def, 
					     $set_id, 
					     $the_void, 
					     $set, 
					     $CTL_OPTIONS->{pred_flank}, 
					     $pred_command,
					     $CTL_OPTIONS->{force},
					     $LOG
					    );
   }
   else{
      die "ERROR: Not a valid predictor in auto_annoator::get_pred_shot\n";
   }
}
#------------------------------------------------------------------------
sub load_transcript_struct {
	my $f            = shift;
	my $g_name       = shift;
	my $i            = shift;
	my $seq          = shift;
	my $evi          = shift;
	my $p_base       = shift;
	my $the_void     = shift;
	my $CTL_OPTIONS  = shift;

	my $transcript_seq  = get_transcript_seq($f, $seq);

        my ($translation_seq, $offset, $end) = get_translation_seq($transcript_seq);
	
	my $len_3_utr = length($transcript_seq) - $end + 1;
	my $l_trans =  length($translation_seq);
	
        my $t_name = ($f->{_tran_name}) ? $f->{_tran_name} : "$g_name-mRNA-$i"; #affects GFFV3.pm

	my $pol_p_hits  = get_selected_types($evi->{gomiph}, 'protein2genome');
	my $pol_e_hits  = get_selected_types($evi->{ests}, 'est2genome', 'est_gff');
	my $blastx_hits = get_selected_types($evi->{gomiph},'blastx', 'protein_gff');
	my $tblastx_hits = get_selected_types($evi->{alt_ests},'tblastx', 'altest_gff');
	my $abinits = $evi->{all_preds};

	my @bag = (@$pol_p_hits,
		   @$pol_e_hits,
		   @$blastx_hits,
		   @$tblastx_hits
		  );
	
	#evidence AED
        my $AED = shadow_AED::get_AED(\@bag, $f);
        my $qi    = maker::quality_index::get_transcript_qi($f,$evi,$offset,$len_3_utr,$l_trans);
	$f->name($t_name);

	if($p_base && $p_base->algorithm !~ /est2genome|est_gff/){
	    $p_base->name($t_name);
	    $p_base->{_AED} = shadow_AED::get_AED(\@bag, $p_base);
	}

        my $t_struct = {'hit'      => $f,
                        't_seq'    => $transcript_seq,
                        'p_seq'    => $translation_seq,
                        't_offset' => $offset,
                        't_end'    => $end,
		        't_name'   => $t_name,
			't_qi'     => $qi,
			'AED'      => $AED,
			'evi'      => $evi,
			'p_base'   => $p_base
                    };

	return $t_struct;
}
#------------------------------------------------------------------------                                                           
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
sub group_transcripts {
   my $data         = shift;
   my $seq          = shift;
   my $seq_id       = shift;
   my $chunk_number = shift;
   my $build        = shift;
   my $predictor    = shift;
   my $the_void     = shift;
   my $CTL_OPTIONS  = shift;

   #place evidence and p_bases in index for easy retrieval
   my @transcripts;
   my %lookup;
   my %p_bases;
   my %snap_lookup;
   my $temp_id = 0;
   foreach my $datum (@{$data}) {
      my $tra       = $datum->[0];
      my $set       = $datum->[1];
      my $p_base    = $datum->[2]; # added 12-15-2006
      
      $tra->{set_id} = $temp_id;
      
      push(@transcripts, $tra);
      
      $lookup{$temp_id} = $set;
      $p_bases{$temp_id} = $p_base;
      $temp_id++;
   }
   
   #cluster the transcripts to get genes
   my $careful_clusters = [];
   if ($predictor eq 'model_gff' ) {
      my %index;
      my $i = 0;
      foreach my $t (@transcripts) {
	 my $j;
	 if (exists $index{$t->{gene_id}}) {
	    $j = $index{$t->{gene_id}};
	 }
	 else {
	    $j = $i;
	    $index{$t->{gene_id}} = $j;
	    $i++;
	 }
	 push(@{$careful_clusters->[$j]}, $t);
      }
   }
   elsif ($predictor eq 'abinit' ) {
      foreach my $t (@transcripts) {
	  push(@{$careful_clusters}, [$t]);
      }
   }
   else {
      $careful_clusters = cluster::careful_cluster_phat_hits(\@transcripts, $seq);
   }
   
   #remove redundant tanscripts in gene
   my @keepers;
   foreach my $c (@{$careful_clusters}) {
      my $best_alt_forms = 
      clean::remove_redundant_alt_splices($c, $seq, 10);
      push(@keepers, $best_alt_forms);
   }	

   #process clusters into genes
   my $c_id = 0;   
   my @annotations;
   foreach my $c (@keepers) {
      my @t_structs;
      
      #build gene name here
      my %pred_sources;
      foreach my $f (@{$c}) {
	 my $source = $f->algorithm;
	 $source =~ s/\:+/_/;
	 $pred_sources{$source}++;
      }
      my $sources = join ('-', keys %pred_sources);
      
      my $g_name;
      if ($predictor eq 'model_gff') {
	 $g_name = $c->[0]->{gene_name}; #affects GFFV3.pm
      }
      elsif ($predictor eq 'abinit') {
	  if ($c->[0]->name =~ /^maker-$seq_id|$seq_id-abinit/) {
	      $g_name = $c->[0]->name;
	      $g_name =~ s/-mRNA-\d.*//;
	  }
	  else{
	      $g_name = "$sources-$seq_id-abinit-gene-$chunk_number"; #affects GFFV3.pm	   
	      $c_id++ while(exists $SEEN->{$c_id} || exists $SEEN->{"$g_name.$c_id"});
	      $g_name = "$g_name.$c_id";
	  }
      }
      else {
	 $g_name = "maker-$seq_id-$sources-gene-$chunk_number"; #affects GFFV3.pm	   
	 $c_id++ while(exists $SEEN->{$c_id} || exists $SEEN->{"$g_name.$c_id"});
	 $g_name = "$g_name.$c_id";
      }
      
      #combine evidence for all transcripts
      my $evidence = {};
      foreach my $f (@{$c}) {
	  my $evi = defined($f->{set_id}) ? $lookup{$f->{set_id}} : {};
	  merge_evidence($evidence, $evi);
      }

      #load transcript structs
      my $AED = 1;
      my $i = 1;
      foreach my $f (@{$c}) {
	 my $p_base = defined($f->{set_id}) ? $p_bases{$f->{set_id}} : undef;

	 my $t_struct = load_transcript_struct($f, $g_name, $i, $seq, $evidence, $p_base, $the_void, $CTL_OPTIONS);

	 push(@t_structs, $t_struct);
	 $AED = $t_struct->{AED} if($t_struct->{AED} < $AED);
	 $i++;
      }
      
      my ($g_start, $g_end, $g_strand) = get_start_and_end_on_seq(\@t_structs);
      
      my $annotation = { 't_structs' => \@t_structs,
			 'g_name'    => $g_name,
			 'g_start'   => $g_start,
			 'g_end'     => $g_end,
			 'g_strand'  => $g_strand,
			 'AED'       => $AED,
			 'predictor' => $predictor
		       };
      
      push(@annotations, $annotation);
      $c_id++;
   }
   
   return \@annotations;
}
#------------------------------------------------------------------------
sub merge_evidence {
    my $evi1 = shift;
    my $evi2 = shift;

    #reset uniq for set to be added
    while(my $key = each %$evi2){
	next if(ref($evi2->{$key}) ne 'ARRAY');
	foreach my $f (@{$evi2->{$key}}){
	    $f->{_uniq_set} = 0;
	}
    }
    #set uniq for existing set
    while(my $key = each %$evi1){
	next if(ref($evi1->{$key}) ne 'ARRAY');
	foreach my $f (@{$evi1->{$key}}){
	    $f->{_uniq_set} = 1;
	}
    }

    #now only add hits where uniq is not set
    while(my $key = each %$evi2){
	next if(ref($evi2->{$key}) ne 'ARRAY');
	$evi1->{$key} = [] if(! $evi1->{$key});
	foreach my $f (@{$evi2->{$key}}){
	    push(@{$evi1->{$key}}, $f) if(! $f->{_uniq_set});
	    $f->{_uniq_set} = 1;
	}
    }
}
#------------------------------------------------------------------------
sub get_non_overlaping_abinits {
   my $ann_set = shift;
   my $abin_set = shift;

   my @overlap;
   my @none;

   #seperate annotations by strand
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
   @p_keepers  = sort {crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @p_keepers;
   @m_keepers  = sort {crit3($a) <=> crit3($b) || crit4($b) <=> crit4($a)} @m_keepers;
   push(@keepers, @{_best(\@p_keepers)});
   push(@keepers, @{_best(\@m_keepers)});

   return (\@keepers);
}
#------------------------------------------------------------------------
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
		elsif ($f->algorithm =~ /exonerate\:\:$type$/i){
		    push(@keepers, $f);
		    last;
		}
	    }
	}

	return \@keepers;
}
#------------------------------------------------------------------------
sub get_transcript_seq {
        my $anno = shift;
        my $seq  = shift;

        my $sorted = PhatHit_utils::sort_hits($anno, 'query');

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
sub get_overlapping_hits {
    my $eat  = shift;
    my $hits = shift;

    my @keepers;
    foreach my $hit (@{$hits}){
	next unless $hit->strand eq $eat->strand;
      push(@keepers, $hit)
	  if compare::overlap($hit, $eat, 'query', 3);
    }
    return \@keepers;
}
#------------------------------------------------------------------------
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
sub get_translation_seq {
	my $seq    = shift;

	my $tM = new CGL::TranslationMachine();

	my ($p_seq , $offset) = $tM->longest_translation_plus_stop($seq);

	my $end = length($p_seq)*3 + $offset + 1;

	$p_seq =~ s/\*$//; # added 11/21/06

	if ($p_seq =~ /^M/){
		# easy  it begins with an M....
		return ($p_seq , $offset, $end );
	}
	else{
		# get the longest internal prot that begins with an M....
		my ($off_new, $p_seq_new) = get_longest_m_seq($seq);

		my $n_end = length($p_seq_new)*3 + $off_new + 1
		            if defined($p_seq_new);

		$p_seq_new =~ s/\*$// if($p_seq_new); # added 11/21/06

		if (!defined($p_seq_new)){
			return ($p_seq , $offset, $end);
		}
		elsif ($offset < 3 && length($p_seq) - length($p_seq_new) > 30){
			return ($p_seq , $offset, $end);
		}
		else {
			return ($p_seq_new, $off_new, $n_end);
		}
	}
}
#------------------------------------------------------------------------
sub get_best_pred_shot {
   my $wanted_strand = shift;
   my $gene_preds    = shift;
   
   my @gs;
   foreach my $g (@{$gene_preds}){
      next unless defined($g);
      next unless defined($g->strand('query'));
      next unless $g->strand('query') == $wanted_strand;
      my $total_score = PhatHit_utils::get_total_score_of_hit($g);
      push(@gs, [$total_score, $g]);
   }
   my @sorted = sort {$b->[0] <=> $a->[0]} @gs;
   
   my $best  = $sorted[0]->[1];
   
   return $best;
}
#------------------------------------------------------------------------
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

	#return $g;
	return $anno_transcript;
}
#------------------------------------------------------------------------
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

