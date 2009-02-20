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
use evaluator::so_classifier;
use evaluator::evaluate;
use evaluator::funs;
use evaluator::AED;
use shadow_AED;
use Storable qw(dclone);

$Storable::forgive_me = 1; 

@ISA = qw(
       );

my $OPT_PREDS; #GLOBAL VARIABLE
my $LOG; #GLOBAL VARIABLE
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

        my $temp_id = 0;
        foreach my $s (@{$preds}){
                $s->{temp_id} = $temp_id;
                $temp_id++;
        }

        my %clusters_hit;
	my %counts;
	my %spans;
	# initialize...

	foreach my $s (@{$preds}){
		$counts{$s->{temp_id}} = 0;

		my ($sB, $sE) = PhatHit_utils::get_span_of_hit($s, 'query');
		   ($sB, $sE) = ($sE, $sB) if $sB > $sE;

		$spans{$s->{temp_id}} = [$sB, $sE];
	}

        print STDERR "Segmenting predictions\n" unless ($main::quiet);

	my $i = 0;
        foreach my $c (@{$careful_clusters}){
        	my $c_strand = $c->[0]->strand('query');
                
                my $coors  = PhatHit_utils::to_begin_and_end_coors($c, 'query');
                my $pieces = Shadower::getPieces($seq, $coors, 10);
                
		my $cB;
		my $cE;

                if (defined($pieces->[1])){	   
		   #die "error in auto_annotator::segment_preds!\n";
		   $pieces = [sort {$a->{b} <=> $b->{b}} @{$pieces}];
		   $cB = $pieces->[0]->{b};
		   $cE = $pieces->[-1]->{e};
		}
		else{
		   $cB = $pieces->[0]->{b};
		   $cE = $pieces->[0]->{e};
		}

        	foreach my $s (@{$preds}){
                	my $s_strand = $s->strand('query');

			next if $s_strand != $c_strand;

			my $pairs = $spans{$s->{temp_id}};

			my $sB = $pairs->[0];
			my $sE = $pairs->[1];

                 	my $class = compare::compare($sB, $sE, $cB, $cE);

                        $clusters_hit{$s->{temp_id}}{$i}++
                        if $class ne '0';

                        $counts{$s->{temp_id}}++ if $class ne '0';
		}
		$i++;
	}

        my @hit_none;
        my @hit_one;
        my @hit_mult;

        foreach my $s (@{$preds}){
                if ($counts{$s->{temp_id}} == 0){
                        push(@hit_none, $s);
                }
                elsif ($counts{$s->{temp_id}} == 1){
                        push(@hit_one, $s);
                }
                else {
                        push(@hit_mult, $s);
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
# sub purge_single_preds {
#    my $careful_clusters = shift;

#    my @c_keepers;
#    foreach my $c (@{$careful_clusters}){
#       my $ests_in_cluster  = get_selected_types($c, 'est2genome');
#       my $ps_in_cluster    = get_selected_types($c,'protein2genome');
#       my $bx_in_cluster    = get_selected_types($c,'blastx');
#       my $snaps_in_cluster = get_selected_types($c,'snap');
#       my $augs_in_cluster = get_selected_types($c,'augustus');
#       if (defined($ests_in_cluster->[0]) ||
# 	  defined($ps_in_cluster->[0]) ||
# 	  defined($bx_in_cluster->[0])
# 	 ){
# 	    push (@c_keepers, $c);
#       }
#    }
   
#    return (\@c_keepers);
# }
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

        my $def   = Fasta::getDef($masked_fasta);
        my $seq_id  = Fasta::def2SeqID($def);
        my $seq_ref   = Fasta::getSeqRef($masked_fasta);
        my $v_seq_ref = Fasta::getSeqRef($virgin_fasta);

	my ($bx_data, $gf_data, $pr_data) = prep_hits($prot_evidence,
						      $est_evidence,
						      $alt_est_evidence,
						      $predictions,
						      $models,
						      $v_seq_ref,
						      $CTL_OPTIONS->{single_exon}
						     );

        my %annotations;
	my @bag;

	#---model passthrough here
	if(grep {/^gff$/} @{$CTL_OPTIONS->{_predictor}}){
	   my $model_trans = run_it($gf_data,
				    $the_void,
				    $seq_ref,
				    $v_seq_ref,
				    $def,
				    'gf',
				    'gff',
				    $the_void,
				    $CTL_OPTIONS
				   );
	   
	   $annotations{'gff'} = group_transcripts($model_trans,
						   $v_seq_ref,
						   $seq_id,
						   $chunk_number,
						   $build,
						   'gff',
						   $the_void,
						   $CTL_OPTIONS
						  );
	   push(@bag, @{$annotations{'gff'}});
	}

	#---hint based gene prediction here
	foreach my $prdr (@{$CTL_OPTIONS->{_predictor}}){
	   next if($prdr eq 'gff' || $prdr eq 'abinit');
	   my $transcripts = run_it($bx_data,
				    $the_void,
				    $seq_ref,
				    $v_seq_ref,
				    $def,
				    'bx',
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
	my $pred_trans = run_it($pr_data,
				$the_void,
				$seq_ref,
				$v_seq_ref,
				$def,
				'pr',
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
{
my %SCORE; #used for sort in function crit2
#this subrutine returns finished MAKER annotations
sub best_annotations {
   my $annotations = shift;
   my $out_base = shift;
   my $CTL_OPTIONS = shift;

   return $annotations->{gff} 
      if(@{$CTL_OPTIONS->{_predictor}} == 1 && grep {/^gff$/} @{$CTL_OPTIONS->{_predictor}});

   my %order; #this is important for sorting
   my $p_list = [];
   my $m_list = [];
   my $p_merge = [];
   my $m_merge = [];
   my $i = @{$CTL_OPTIONS->{_predictor}};
   foreach my $p (@{$CTL_OPTIONS->{_predictor}}){
      $order{$p} = $i;
      $i--;
      foreach my $g (@{$annotations->{$p}}){
	 if($g->{g_strand} == 1){
	     next unless($g->{AED} < 1 || $p eq 'gff');
	     push(@$p_list, $g);
	     push(@$p_merge, $g) if($g->{t_structs}->[0]->{hit}->{_hit_multi});
	 }
	 elsif($g->{g_strand} == -1){
	     next unless($g->{AED} < 1 || $p eq 'gff');
	     push(@$m_list, $g);
	     push(@$m_merge, $g) if($g->{t_structs}->[0]->{hit}->{_hit_multi});
	 }
	 else{
	    die "ERROR: Logic error in auto_annotator::best_annotations\n";
	 }
      }
   }
   %SCORE = %order; #must be set for sort to work correctly

   #remove low scoring overlaping genes
   @$p_list  = sort {crit1($a) <=> crit1($b) || crit2($b) <=> crit2($a)|| crit3($b) cmp crit3($a)} @$p_list;
   @$m_list  = sort {crit1($a) <=> crit1($b) || crit2($b) <=> crit2($a)|| crit3($b) cmp crit3($a)} @$m_list;
   $p_list = _best($p_list);
   $m_list = _best($m_list);

   #use abinit that merges clusters if average AED improves
   my $p_m_list = [];
   foreach my $g (@$p_merge){
       my $g_B = $g->{g_start};
       my $g_E = $g->{g_end};

       my $tot_length;
       my $tot_AED;

       foreach my $k (@$p_list){
           my $k_B = $k->{g_start};
           my $k_E = $k->{g_end};

           my $class = compare::compare($g_B, $g_E, $k_B, $k_E);

           if($class ne '0'){
	       my $nc = get_n_count($k);
	       $tot_length += $nc;
	       $tot_AED += $nc * $k->{AED}
           }
       }

       next if(!$tot_length);

       my $ave_AED = $tot_AED/$tot_length;
 
       push(@$p_m_list, $g) if($g->{AED} < $ave_AED);
   }

   my $m_m_list = [];
   foreach my $g (@$m_merge){
       my $g_B = $g->{g_start};
       my $g_E = $g->{g_end};

       my $tot_length;
       my $tot_AED;

       foreach my $k (@$m_list){
           my $k_B = $k->{g_start};
           my $k_E = $k->{g_end};

           my $class = compare::compare($g_B, $g_E, $k_B, $k_E);

           if($class ne '0'){
	       my $nc = get_n_count($g);
	       $tot_length += $nc;
	       $tot_AED += $nc * $g->{AED}
           }
       }

       next if(!$tot_length);

       my $ave_AED = $tot_AED/$tot_length;

       push(@$m_m_list, $g) if($g->{AED} < $ave_AED);
   }

   @$p_m_list = sort {crit1($a) <=> crit1($b) || crit2($b) <=> crit2($a)|| crit3($b) cmp crit3($a)} @$p_m_list;
   @$m_m_list  = sort {crit1($a) <=> crit1($b) || crit2($b) <=> crit2($a)|| crit3($b) cmp crit3($a)} @$m_m_list;
   my $p_keepers = _best($p_m_list);
   my $m_keepers = _best($m_m_list);

   #now combine mergers with other keepers
   push(@$p_keepers, @$p_list);
   push(@$m_keepers, @$m_list);
   $p_keepers = _best($p_keepers);
   $m_keepers = _best($m_keepers);
   push(@$p_keepers, @$m_keepers);

   write evaluator reports
   foreach my $ann (@$p_keepers){
       foreach my $t (@{$ann->{t_structs}}){
   	   my $dir = "$out_base/evaluator";
   	   mkdir($dir) if(! -e $dir);
   	   my $file = "$dir/".Fasta::seqID2SafeID($t->{hit}->name).".eva";
   	   open(my $FH, "> $file");
   	   print $FH $t->{report};
   	   close($FH);
       }       
   }

   return $p_keepers;
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
#sort by evaluator score
sub crit1 {
   my $g = shift;
   
   return $g->{AED};
}
#------------------------------------------------------------------------
#sort by order given in control files
sub crit2 {
   my $g = shift;

   my $score = $SCORE{$g->{predictor}} || 0;
   return $score;
}
#------------------------------------------------------------------------

sub crit3 {
   my $g = shift;
   
   return $g->{predictor};
}
}
#------------------------------------------------------------------------
sub run_it {
   my $data         = shift;
   my $the_void     = shift;
   my $seq          = shift;
   my $v_seq        = shift;
   my $def          = shift;
   my $id           = shift;
   my $predictor    = shift;
   my $CTL_OPTIONS  = shift;

   my $q_id = Fasta::def2SeqID($def); 
   my @transcripts;
   my $i = 0;
   foreach my $set (@{$data}) {
      my $goimph   = $set->{gomiph};
      my $mia      = $set->{mia};
      my $c_id     = $set->{c_id};
      my $ests     = $set->{ests};
      my $alt_ests = $set->{alt_ests};
      my $model    = $set->{model};
	   
      if ($predictor eq 'gff' || $predictor eq 'abinit') {
	 next if(! defined $model);
	 my $transcript = $model;

	 my $utr_trans; 
	 if($predictor eq 'abinit' && defined($ests)){
	     my $copy = dclone($transcript);
	     $utr_trans = pneu($ests, $copy, $seq);
	     push(@transcripts, [$utr_trans, $set, $transcript]);
	     if (! clean::is_identical_form($utr_trans, $transcript)){
		 push(@transcripts, [$utr_trans, $set, undef]);
	     }
	 }
	 push(@transcripts, [$transcript, $set, undef]);

	 $i++;
	 next;
      }

      if ($predictor eq 'est2genome') {
	 next if (! defined $mia);
	 my $transcript = pneu($ests, $mia, $seq);
	 push(@transcripts, [$transcript, $set, undef]) if $transcript;
	 $i++;
	 next;
      }

      my ($pred_shots, $strand) = get_pred_shot($seq, 
						$def, 
						$id.'.'.$i, 
						$the_void, 
						$set,
						$predictor,
						$CTL_OPTIONS,
						$LOG
					       );

      my $on_right_strand = get_best_pred_shots($strand, $pred_shots);

      #commented out 12/11/2008 to remove redundancy when est2genome is also a predictor
      #my $best_pred       = get_best_pred_shot($strand, $pred_shots);      
      #if (!defined($best_pred) && defined($mia)) {
      #   my $transcript = pneu($ests, $mia, $seq);
      #   push(@transcripts, [$transcript, $set, undef]);
      #}	  
	   
      foreach my $pred_shot (@{$on_right_strand}) {
	 my $copy = $pred_shot;
	 if (defined($pred_shot) && defined($mia)) {
	    my $transcript = pneu($ests, $copy, $seq);
	    if (defined($transcript)) {
	       push(@transcripts, [$transcript, $set, $pred_shot])
	    }
	    else {
	       push(@transcripts, [$copy, $set, $pred_shot])
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
	my $so_code      = shift;
	my $geneAED      = shift;
	my $alt_spli_sup = shift;
	my $the_void     = shift;
	my $CTL_OPTIONS  = shift;

	my $transcript_seq  = get_transcript_seq($f, $seq);

        my ($translation_seq, $offset, $end) = get_translation_seq($transcript_seq);
	

	my $len_3_utr = length($transcript_seq) - $end + 1;
	my $l_trans =  length($translation_seq);
	
        my $t_name = "$g_name-mRNA-$i"; #affects GFFV3.pm

	#----evaluator here
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

	#holds evalutor struct
        my $eva = evaluator::evaluate::power_evaluate($f,
						      $seq,
						      $pol_p_hits,
						      $pol_e_hits,
						      $blastx_hits,
						      $abinits,
						      $so_code,
						      $geneAED,
						      $alt_spli_sup,
						      $t_name,
						      $CTL_OPTIONS,
						      $the_void,
						     );
	
        my $AED   = shadow_AED::get_AED(\@bag, $f);#$eva->{txnAED};

        my $score = $eva->{score};
        my $qi    = $eva->{qi};#maker::quality_index::get_transcript_qi($f,$evi,$offset,$len_3_utr,$l_trans);#
	my $report = $eva->{report};
	#----
	$t_name .= " AED:";
	$t_name .= sprintf '%.2f', $AED; # two decimal places
	$t_name .= " $qi";
	$f->name($t_name); #give name to hit

        my $t_struct = {'hit'      => $f,
                        't_seq'    => $transcript_seq,
                        'p_seq'    => $translation_seq,
                        't_offset' => $offset,
                        't_end'    => $end,
		        't_name'   => $t_name,
			't_qi'     => $qi,
			'AED'      => $AED,
			'score'    => $score,
			'report'   => $report
                    };

	return $t_struct;
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
{
my %SEEN; #holds names that have already been seen from gff passthrough
sub group_transcripts {
   my $data         = shift;
   my $seq          = shift;
   my $seq_id       = shift;
   my $chunk_number = shift;
   my $build        = shift;
   my $predictor    = shift;
   my $the_void     = shift;
   my $CTL_OPTIONS  = shift;
   
   my @transcripts;
   my %lookup;
   my %snap_lookup;
   my $temp_id = 0;
   foreach my $datum (@{$data}) {
      my $tra       = $datum->[0];
      my $set       = $datum->[1];
      my $o_snap    = $datum->[2]; # added 12-15-2006
      
      $tra->{set_id} = $temp_id;
      
      push(@transcripts, $tra);
      
      $lookup{$temp_id} = $set;
      $temp_id++;
   }
   
   my $careful_clusters = [];
   if ($predictor eq 'gff' ) {
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
	 push(@{$careful_clusters}, [$t])
      }
   }
   else {
      $careful_clusters = cluster::careful_cluster_phat_hits(\@transcripts, $seq);
   }
   
   my @keepers;
   foreach my $c (@{$careful_clusters}) {
      my $best_alt_forms = 
      clean::remove_redundant_alt_splices($c, $seq, 10);
      push(@keepers, $best_alt_forms);
   }	
   my $c_id = 0;
   
   my @annotations;
   foreach my $c (@keepers) {
      my @t_structs;
      
      my %pred_sources;
      foreach my $f (@{$c}) {
	 my $source = $f->algorithm;
	 $source =~ s/\:+/_/;
	 $pred_sources{$source}++;
      }
      my $sources = join ('-', keys %pred_sources);
      
      my $g_name;
      if ($predictor eq 'gff') {
	 $g_name = $c->[0]->{gene_name}; #affects GFFV3.pm
	 if ($g_name =~ /^maker-$seq_id/) {
	    $SEEN{$g_name}++;
	 }
      }
      elsif ($predictor eq 'abinit') {
	 $g_name = "$sources-$seq_id-abinit-gene-$chunk_number"; #affects GFFV3.pm	   
	 $c_id++ while(exists $SEEN{"$g_name.$c_id"});
	 $g_name = "$g_name.$c_id";
      }
      else {
	 $g_name = "maker-$seq_id-$sources-gene-$chunk_number"; #affects GFFV3.pm	   
	 $c_id++ while(exists $SEEN{"$g_name.$c_id"});
	 $g_name = "$g_name.$c_id";
      }
      
      #----evaluator here
      my @pol_e_hits;
      my @pol_p_hits;
      my @blastx_hits;
      my @ab_inits;

      foreach my $f (@{$c}) {
	 my $evi = defined($f->{set_id}) ? $lookup{$f->{set_id}} : [];
	 my $ests = get_selected_types($evi->{ests}, 'est2genome', 'est_gff');
	 my $prot2gen = get_selected_types($evi->{gomiph}, 'protein2genome');
	 my $prot = get_selected_types($evi->{gomiph}, 'blastx', 'protein_gff');
	 my $ab_init = $evi->{all_preds};
	 
	 #remove redundant evidence (some added > 1 time)
	 foreach my $e (@{$ests}, @{$prot2gen}, @{$prot}, @{$ab_init}){
	    $e->{_uniq_set} = 0;  #reset _uniq_set before beginning
	 }
	 foreach my $e (@{$ests}){
	    if(! $e->{_uniq_set}){
	       $e->{_uniq_set} = 1;
	       push(@pol_e_hits, $e);
	    }
	 }
         foreach my $e (@{$prot2gen}){
	     if(! $e->{_uniq_set}){
		 $e->{_uniq_set} = 1;
		 push(@pol_p_hits, $e);
	     }
         }
         foreach my $e (@{$prot}){
	     if(! $e->{_uniq_set}){
		 $e->{_uniq_set} = 1;
		 push(@blastx_hits, $e);
	     }
         }
	 foreach my $e (@{$ab_init}) {
             if(! $e->{_uniq_set}){
                 $e->{_uniq_set} = 1;
		 push(@ab_inits, $e) ;
		}
	 }

      }
      
      my $so_code = evaluator::so_classifier::so_code($c);
      my $alt_spli_sup = evaluator::funs::alt_spli($c, \@pol_e_hits, $seq);
      my $geneAED = evaluator::AED::gene_AED($c, \@pol_e_hits, \@pol_p_hits, \@blastx_hits, \@ab_inits, $seq);
      #----
      
      my $score = 0;
      my $AED = 1;
      my $i = 1;
      foreach my $f (@{$c}) {
	 my $evidence = defined($f->{set_id}) ? $lookup{$f->{set_id}} : undef;
	 
	 my $t_struct = 
	 load_transcript_struct($f, $g_name, $i, $seq, $evidence, $so_code,
				$geneAED, $alt_spli_sup, $the_void, $CTL_OPTIONS);

	 if(!$f || ! $t_struct){
	     sleep 1;
	 }

	 push(@t_structs, $t_struct);
	 $score = $t_struct->{score} if($t_struct->{score} > $score);
	 $AED = $t_struct->{AED} if($t_struct->{AED} < $AED);
	 $i++;
      }
      
      my ($g_start, $g_end, $g_strand) = get_start_and_end_on_seq(\@t_structs);
      
      my $annotation = { 't_structs' => \@t_structs,
			 'g_name'    => $g_name,
			 'g_start'   => $g_start,
			 'g_end'     => $g_end,
			 'g_strand'  => $g_strand,
			 'score'     => $score,
			 'AED'       => $geneAED,
			 'predictor' => $predictor,
			 'so_code'   => $so_code
		       };
      
      push(@annotations, $annotation);
      $c_id++;
   }
   
   return \@annotations;
}
}
#------------------------------------------------------------------------
sub abinit_annotation_overlap {
   my $abin_set = shift;
   my $ann_set = shift;

   my @overlap;
   my @none;

   foreach my $abin (@{$abin_set}){
      my $over = 0;
      foreach my $ann (@{$ann_set}){
	 next if ($abin->{g_strand} != $ann->{g_strand});
	 my $comp = compare::compare ($abin->{g_start},
				      $abin->{g_end},
				      $ann->{g_start},
				      $ann->{g_end}
				     );
	 if($comp ne '0'){
	    $over = 1;
	    last;
	 }
      }
      if($over){
	 push (@overlap, $abin);
      }
      else{
	 push (@none, $abin);
      }
   }
   
   return (\@overlap, \@none);
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
		   }
		   elsif ($f->algorithm =~ /^$type\:/i){
		      push(@keepers, $f);
		   }
		   elsif ($f->algorithm =~ /exonerate\:\:$type$/i){
		      push(@keepers, $f);
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


