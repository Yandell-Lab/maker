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
	my $alt_est_hits     = shift;
	my $predictions      = shift;
	my $models           = shift;
	my $seq              = shift;
	my $single_exon      = shift;
	my $single_length    = shift;

	#only ESTs from splice concious algorithms i.e. no blastn
	my $clean_est = get_selected_types($est_hits,'est2genome', 'est_gff');
	my $clean_altest = $alt_est_hits;
	if ($single_exon == 1) {
	    #do nothing
	}else {
	    # don't use unpliced single exon ESTs-- may be genomic contamination
	    $clean_est = clean::purge_single_exon_hits($clean_est);
	    $clean_altest = clean::purge_single_exon_hits($clean_altest);
	}

	# throw out the exonerate est hits with weird splice sites
	$clean_est = clean::throw_out_bad_splicers($clean_est, $seq);
	

	#---build clusters for basic evidence

	# combine puts type in order they are given
	# important for later culstering, as hits in
	# first arg more likey to make in into to cluster
	# than second arg, etc	
	my $c_bag = combine($prot_hits,
	                    $clean_est, 
			    $clean_altest
			   );

        my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $c_bag);
        my $p_clusters = cluster::shadow_cluster(0, $seq, $p, 10);
        my $m_clusters = cluster::shadow_cluster(0, $seq, $m, 10);

	#purge after clustering so as to still have the effect of evidence joining ESTs
	$p_clusters = purge_short_ESTs_in_clusters($p_clusters, $single_length);
	$m_clusters = purge_short_ESTs_in_clusters($m_clusters, $single_length);

	my $careful_clusters = [];
	push(@{$careful_clusters}, @{$p_clusters});
	push(@{$careful_clusters}, @{$m_clusters});

	#---now start preparing data for different types of input

	#--model_gff3 input
        # identify the models that fall within and between basic clusters
        my ($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($models,
                                                                      $careful_clusters,
                                                                      $seq,
								      );

	#join the clusters on the models
	my $model_clusters = join_clusters_on_pred($models, $careful_clusters, $c_index);



	#--abinit input
        ($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
								   $careful_clusters,
								   $seq,
								   );

	#join clusters on the ab-inits
	my $pred_clusters = join_clusters_on_pred($predictions, $careful_clusters, $c_index);


	#--build clusters joined together by models (for hint based annotations)
        my $hint_clusters = [];
	if(@$models){
	    my $m_bag = combine($models,
				$prot_hits,
				$clean_est,
				$clean_altest
				);
	    
	    ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $m_bag);
	    $p_clusters = cluster::shadow_cluster(0, $seq, $p, 10);
	    $m_clusters = cluster::shadow_cluster(0, $seq, $m, 10);

	    #purge after clustering so as to still have the effect of evidence joining ESTs
	    $p_clusters = purge_short_ESTs_in_clusters($p_clusters, $single_length);
	    $m_clusters = purge_short_ESTs_in_clusters($m_clusters, $single_length);
	    
	    
	    push(@{$hint_clusters}, @{$p_clusters});
	    push(@{$hint_clusters}, @{$m_clusters});
	}
	else{
	    $hint_clusters = $careful_clusters;
	}

	# identify the abinits that fall within and between clusters
	($c_index, $hit_one, $hit_none, $hit_mult) = segment_preds($predictions,
                                                                   $hint_clusters,
                                                                   $seq,
								   );

	#only add as hints those that don't cross multiple clusters
	merge_into_cluster($hit_one, $hint_clusters, $c_index);

	#--prep hint data
	my $c_id = 0;
	my @bx_data;
	foreach my $c (@{$hint_clusters}){
	   my $bx = prep_blastx_data($c, $c_id, $seq);
	   push(@bx_data, @{$bx}) if defined $bx;

	   $c_id++;
	}

	#--prep model_gff data
	my @gf_data;
	foreach my $c (@{$model_clusters}){
	    my $gf = prep_gff_data($c, $c_id, $seq);
	    push(@gf_data, @{$gf}) if defined $gf;
	    
	    $c_id++;
	}

	#--prep abinit data
	my @pr_data;
	foreach my $c (@{$pred_clusters}){
	    my $pr = prep_pred_data($c, $c_id, $seq);
	    push(@pr_data, @{$pr}) if defined $pr;
	    
	    $c_id++;
	}

	return (\@bx_data, \@gf_data, \@pr_data);
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
#takes abinits identified by segment_preds and merges them back into
#all the clusters they overlap
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
#used to identify hits related to ESTs in a cluster and then remove those
#below a certain threshold, returns the revised cluster
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
        my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark', 'pred_gff');

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
#called by prep_hits for standard evidence clusters
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
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark', 'pred_gff');

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
			   'preds'     => $preds_in_cluster,
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
			'preds'     => $preds_in_cluster,
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
        my $ests_in_cluster  = get_selected_types($c,'est2genome', 'est_gff');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx', 'protein_gff');
        my $alt_ests_in_cluster = get_selected_types($c,'tblastx', 'altest_gff');
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark',  'pred_gff');

        # groups of most informative protein hits
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;

	foreach my $model (@{$models_in_cluster}){
	   push(@data, {'gomiph'    => $gomiph,
			'preds' => $preds_in_cluster,
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
#called by prep_hits for abinit based  clusters
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
	my $preds_in_cluster = get_selected_types($c,'snap', 'augustus', 'fgenesh',
						  'twinscan', 'genemark', 'pred_gff');

        # groups of most informative protein hits
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

	my @data;

	foreach my $pred (@{$preds_in_cluster}){
	   push(@data, {'gomiph'    => $gomiph,
			'preds'     => $preds_in_cluster,
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
#called by prep_hits for standard evidence clusters
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
#called by prep_hits for standard evidence clusters
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
#this subrutine returns finished annotations for all predictors.
#called outside of package by maker
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
    #my $GFF_DB = new GFFDB($CTL_OPTIONS);
    $SEEN = {};#$GFF_DB->get_existing_gene_names($seq_id);
    
    #group evidence and predictions
    my ($bx_data, $gf_data, $pr_data) = prep_hits($prot_evidence,
						  $est_evidence,
						  $alt_est_evidence,
						  $predictions,
						  $models,
						  $v_seq_ref,
						  $CTL_OPTIONS->{single_exon},
						  $CTL_OPTIONS->{single_length}
						  );
    
    my %annotations;
    
    #---model passthrough here
    if(@$gf_data){
	print STDERR "Processing GFF3 passthrough annotations\n" unless($main::quiet);
	my $model_trans = run_it($gf_data,
				 $the_void,
				 $seq_ref,
				 $v_seq_ref,
				 $def,
				 'model_gff',
				 $predictions,
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
    
    #---hint based gene prediction here (includes est2genome)
    foreach my $prdr (@{$CTL_OPTIONS->{_predictor}}){
	next if($prdr eq 'model_gff' || $prdr eq 'pred_gff');
	print STDERR "Producing $prdr hint based annotations\n" unless($main::quiet);
	
	my $transcripts = run_it($bx_data,
				 $the_void,
				 $seq_ref,
				 $v_seq_ref,
				 $def,
				 $prdr,
				 $predictions,
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
			    'abinit', #all abinits not just pref_gff
			    $predictions,
			    $CTL_OPTIONS
			    );
    
    my $all_ab = group_transcripts($pred_trans,
				   $v_seq_ref,
				   $seq_id,
				   $chunk_number,
				   $build,
				   'abinit', #all abinits not_just pred_gff
				   $the_void,
				   $CTL_OPTIONS
				   );

    $annotations{'abinit'} = $all_ab; #all abinit

    #add abinits to their predictor type after they have been proccessed
    #remeber they get treated differently so you want to add them as a
    #seperate step.
    foreach my $g (@$all_ab){
	if($g->{algorithm} =~ /^snap/i){
            push(@{$annotations{'snap'}}, $g);
	}
	elsif($g->{algorithm} =~ /^genemark/i){
	    push(@{$annotations{'genemark'}}, $g);
	}
	elsif($g->{algorithm} =~ /^augustus/i){
            push(@{$annotations{'augustus'}}, $g);
	}
	elsif($g->{algorithm} =~ /^fgenesh/i){
            push(@{$annotations{'fgenesh'}}, $g);
	}
	elsif($g->{algorithm} =~ /^pred_gff/i){
            push(@{$annotations{'pred_gff'}}, $g);
	}
	else{
	    die "ERROR: Not a supported algorithm: ".$g->{algorithm}."\n";
	}
    }
    
    add_abAED(\%annotations);

    return \%annotations;
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
	next if ($p eq 'abinit'); #these are redundant inside predictor type
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
    my $out_base = shift;
    my $CTL_OPTIONS = shift;

    print STDERR "Choosing best annotations\n" unless($main::quiet);
    
    my @p_keepers;
    my @m_keepers;

    #keep all gff3 passthrough if there's nothing else
    if(@{$CTL_OPTIONS->{_predictor}} == 1 && $CTL_OPTIONS->{_predictor}->[0] eq 'model_gff'){
	foreach my $g (@{$annotations->{'model_gff'}}){
	    if($g->{g_strand} == 1){
		push(@p_keepers, $g);
	    }
	    elsif($g->{g_strand} == -1){
		push(@m_keepers, $g);
	    }
	}
    }
    elsif(@{$CTL_OPTIONS->{_predictor}}){
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

            #then check abinits (slow)
	    my $preds = $g->{g_evidence}->{all_preds};
            my $pAED = shadow_AED::get_AED($preds, $t->{hit});
            if($pAED < 1){
                push(@p_final,$g);
                $add_flag = 1;
                last;
            }

	    #then check exonerate proteins (slow)
	    my $gomiph = $g->{g_evidence}->{gomiph};
	    my $p_ex = [];
	    $p_ex = get_selected_types($gomiph,'protein2genome') if($gomiph);
	    my $eAED = shadow_AED::get_AED($p_ex, $t->{hit});
	    if($eAED < 1){
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

            #then check abinits (slow)
	    my $preds = $g->{g_evidence}->{all_preds};
            my $pAED = shadow_AED::get_AED($preds, $t->{hit});
            if($pAED < 1){
                push(@m_final,$g);
                $add_flag = 1;
                last;
            }

	    #then check exonerate proteins (slow)
	    my $gomiph = $g->{g_evidence}->{gomiph};
	    my $p_ex = [];
	    $p_ex = get_selected_types($gomiph,'protein2genome') if($gomiph);
	    my $eAED = shadow_AED::get_AED($p_ex, $t->{hit});
	    if($eAED < 1){
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
   
   return ($g->{AED} + $g->{abAED})/2;
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
    my $predictions  = shift;
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
	
    	    my $all_preds = get_overlapping_hits($transcript, $predictions);
	    $set->{all_preds} = $all_preds;
	    push(@transcripts, [$transcript, $set, undef]);
	    
	    $i++;
	    next;
	}

	#------ab-init passthrough
	if ($predictor eq 'abinit') {
	    next if(! defined $model);
	    
	    #add UTR to ab-inits
	    my $transcript = pneu($ests, $model, $seq);

	    my $all_preds = get_overlapping_hits($transcript, $predictions);
            $set->{all_preds} = $all_preds;

	    push(@transcripts, [$transcript, $set, $model]);
	    
	    $i++;
	    next;
	}
	
	#------est2genome
	if ($predictor eq 'est2genome') {
	    next if (! defined $mia);
	    my $transcript = pneu($ests, $mia, $seq);
	    
	    next if(! $transcript);

            my $all_preds = get_overlapping_hits($transcript, $predictions);
            $set->{all_preds} = $all_preds;

	    push(@transcripts, [$transcript, $set, $mia]);
	    $i++;
	    next;
	}

        #------genemark does not have hints enabled 
        return [] if ($predictor eq 'genemark');
	
	#------default hint based behavior
	my $gomiph   = $set->{gomiph};
	my $blastx   = get_selected_types($gomiph,'blastx', 'protein_gff');
	my $pol_p    = get_selected_types($gomiph,'protein2genome');
	my $alt_ests = $set->{alt_ests};
	my $preds    = $set->{preds};
	
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
	if(@$on_right_strand > 1 && ! defined $mia && ! @$pol_p && ! @$preds){
	    my $clean  = clean::purge_single_exon_hits($alt_ests);
	    my $coors  = PhatHit_utils::get_hsp_coors($blastx, 'query');
	    my $pieces = Shadower::getPieces($seq, $coors, 0);
	    
	    if(! @$clean && @$pieces <= 1){
		my $keep;
		foreach my $h (@$on_right_strand){
		    my $set = get_overlapping_hits($h, $predictions);
		    $keep = 1 if(@$set);
		}

		next if(! $keep); #skip these spurious predictions
	    }
	}
	
	#add transcripts
	foreach my $pred_shot (@{$on_right_strand}) {
	    if (defined($pred_shot)){
		my $transcript = (defined($mia)) ? pneu($ests, $pred_shot, $seq) : $pred_shot;
		
		my $all_preds = get_overlapping_hits($transcript, $predictions);
		$set->{all_preds} = $all_preds;
		
		push(@transcripts, [$transcript, $set, $pred_shot]);
	    }
	}
	$i++;
    }
    
    return \@transcripts;
}
#------------------------------------------------------------------------
#runs the gene prdictors with hints.  Called by run_it.
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
#takes the gene predictions and evidence and builds transcript name,
#QI, AED and gets protein and mRNA sequence then puts it al into a
#HASH.  Called by group transcripts. Transcript name is set here, as
#well as UTR boudaries
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
			'p_base'   => $p_base,
			'p_length' => length($translation_seq)
                    };

	return $t_struct;
}
#------------------------------------------------------------------------                                                           
#takes an array of annotations and only returns those that overlap a maker
#annotation (doesn't care if on opposite strand so be careful)
sub get_geness_overlapping_gene {
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
   
   #remove redundant transcripts in gene
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
	 $SEEN->{$g_name}++;
	 if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
	     $SEEN->{$1}++;
	 }
      }
      elsif ($predictor eq 'abinit') {
	  if ($c->[0]->name =~ /^maker-$seq_id|$seq_id-abinit/) {
	      $g_name = $c->[0]->name;
	      $g_name =~ s/-mRNA-\d.*//;
	      $SEEN->{$g_name}++; 
	      if($g_name =~ /(\d+\.\d+)\-mRNA\-\d+/){
		  $SEEN->{$1}++;
	      }
	  }
	  else{
	      $g_name = "$sources-$seq_id-abinit-gene-$chunk_number"; #affects GFFV3.pm	   
	      $c_id++ while(exists $SEEN->{"$chunk_number\.$c_id"} || exists $SEEN->{"$g_name.$c_id"});
	      $g_name = "$g_name.$c_id";
	      $SEEN->{$g_name}++; 
	      $SEEN->{"$chunk_number\.$c_id"}++;
	  }
      }
      else {
	 $g_name = "maker-$seq_id-$sources-gene-$chunk_number"; #affects GFFV3.pm
	 $c_id++ while(exists $SEEN->{"$chunk_number\.$c_id"} || exists $SEEN->{"$g_name.$c_id"});
	 $g_name = "$g_name.$c_id";
	 $SEEN->{$g_name}++; 
	 $SEEN->{"$chunk_number\.$c_id"}++; 
      }
      
      #combine evidence for all transcripts
      my $evidence = {};
      foreach my $f (@{$c}) {
	  my $evi = defined($f->{set_id}) ? $lookup{$f->{set_id}} : {};
	  merge_evidence($evidence, $evi);
      }

      #load transcript structs
      my $AED = 1;
      my $p_length = 0;
      my $i = 1;
      foreach my $f (@{$c}) {
	 my $p_base = defined($f->{set_id}) ? $p_bases{$f->{set_id}} : undef;

	 my $t_struct = load_transcript_struct($f, $g_name, $i, $seq, $evidence, $p_base, $the_void, $CTL_OPTIONS);

	 push(@t_structs, $t_struct) unless ($t_struct->{p_length} < $CTL_OPTIONS->{min_protein});

	 $AED = $t_struct->{AED} if($t_struct->{AED} < $AED);
	 $i++;
      }
      
      if(! @t_structs){
	  $c_id++;
	  next;
      }

      my ($g_start, $g_end, $g_strand) = get_start_and_end_on_seq(\@t_structs);
      my $g_attrib = (exists $t_structs[0]->{hit}->{gene_attrib}) ? $t_structs[0]->{hit}->{gene_attrib} : undef;

      my $annotation = { 't_structs'  => \@t_structs,
			 'g_name'     => $g_name,
			 'g_start'    => $g_start,
			 'g_end'      => $g_end,
			 'g_strand'   => $g_strand,
			 'g_evidence' => $evidence,
			 'AED'        => $AED,
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
#merges evidence of data supporting annotaitons.  Merges them so that
#all members are uniq and not repeated.  This is important for combining
#evidence from mutliple transcripts as a single gene.

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

    #seperate annotations by strand
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
		$fg = $g_index[$cl]; #get annotation for that gene
	    }
	    $t->{is_changed} = ($AEDs[$id] == 0) ? 0 : 1; #not a true change if identical
	    $g->{is_changed} = 1 if($t->{is_changed}); #gene changed only if trans changed
	    $t->{t_name} = $f->{_tran_name}; #set transcript name

	    #this helps with attribute passthrough
	    $t->{hit} = $f if(! $t->{is_changed}); #use gff hit if they are identical

	    #!!remember hit name has not been changed on non identical!! ;-)
	}

	next if(! $new); #no change so skip

	$g->{g_name} = $new; #set new gene name

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
   my $abin_set = shift;
   my $CTL_OPT = shift;

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

   #separate abinits by strand (only masked unless I'm not masking)
   my @p_ab;
   my @m_ab;
   foreach my $g (@$abin_set){
       if($g->{g_strand} == 1){
	   push(@p_ab, $g) if($g->{algorithm} =~ /_masked$|^pred_gff/ || $CTL_OPT->{_no_mask});
       }
       elsif($g->{g_strand} == -1){
	   push(@m_ab, $g) if($g->{algorithm} =~ /_masked$|^pred_gff/ || $CTL_OPT->{_no_mask});
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
#i.e. used to seperate est2genome or blastx hits from a mixed array
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
	next unless $hit->strand eq $eat->strand;
      push(@keepers, $hit)
	  if compare::overlap($hit, $eat, 'query', 3);
    }
    return \@keepers;
}
#------------------------------------------------------------------------
#returns longest translation starting with M and the offet for a
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
#takes the gene predictions produce by get_pred_shot and finds the one
#with the higest score.  Called by run_it
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

