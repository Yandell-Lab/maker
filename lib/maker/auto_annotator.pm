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
use SimpleCluster;
use CGL::TranslationMachine;
use compare;
use cluster;
use clean;
use maker::join;
use maker::quality_index;
@ISA = qw(
       );

my $OPT_F; #GLOBAL VARIABLE
my $OPT_SNAPS; #GLOBAL VARIABLE
my $OPT_PRED; #GLOBAL VARIABLE
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub prep_hits {
	my $exonerate_p_hits = shift;
	my $exonerate_e_hits = shift;
	my $blastx_keepers   = shift;
	my $snap_predictions = shift;
	my $seq              = shift;
	my $single_exon      = shift;
	
	my $clean_exonerate_e;
	
	if ($single_exon == 1) {
	    $clean_exonerate_e = $exonerate_e_hits;
	}else {
	    # don't use unpliced ESTs-- usually genomic contamination
	    my $spliced_exonerate_e = clean::purge_single_exon_hits($exonerate_e_hits);
	    # throw out the exonerate est hits with weird splice sites
	    $clean_exonerate_e = clean::throw_out_bad_splicers($spliced_exonerate_e);
	}
	
	# combine puts type in order they are given
	# important for later culstering, as hits in
	# first arg more likey to make in into to cluster
	# than second arg, etc

	my $c_bag = combine($exonerate_p_hits,
	                    $clean_exonerate_e, 
		             $blastx_keepers,
			    );

        #my $s_bag = combine($exonerate_p_hits,
        #                    $blastx_keepers
        #                    );

	# preferred method if the est data is good. This was standard untill 12/11/06
	#my $careful_clusters = cluster::careful_cluster_phat_hits($c_bag, $seq, 50);

	# -- BEGIN nGASP modification
	# I  set this up to deal with nGASP
        my ($p, $m, $x, $z) = PhatHit_utils::seperate_by_strand('query', $c_bag);

        my $p_clusters = cluster::shadow_cluster(30, $seq, $p, 10);
        my $m_clusters = cluster::shadow_cluster(30, $seq, $m, 10);

	my $careful_clusters = [];

	push(@{$careful_clusters}, @{$p_clusters}) if defined $p_clusters->[0];
	push(@{$careful_clusters}, @{$m_clusters}) if defined $m_clusters->[0];

        my $temp_id = 0;
        foreach my $s (@{$snap_predictions}){
                $s->{temp_id} = $temp_id;
                $temp_id++;
        }

	
	# identify the snap ab-inits that fall within and between clusters
	my ($c_index, $hit_one, $hit_none, $hit_mult) = segment_snaps($snap_predictions, 
	                                                              $careful_clusters,
	                                                              $seq,
	                                                              );
	
	# add the ab-initio snaps that hit only a single cluster
        foreach my $s (@{$hit_one}){
                my @keys = keys %{$c_index->{$s->{temp_id}}};
                my $i = $keys[0];
                die "logic error in segment_snaps\n" if defined($keys[1]);
	
                push(@{$careful_clusters->[$i]}, $s);
        }
	#--- end nGASP modification;

	#compromise method becase we don't have reliable strand info for S. med.
	# no longer needs because of exonerate::splice_info.pm
	#my $careful_clusters = cluster::special_cluster_phat_hits($s_bag, 
	#                                                          $clean_exonerate_e, 
	#                                                          $seq, 
	#                                                          50,
	#                                                          );
	
	my @bx_data;
	my @pp_data;
	my @pe_data;
	my $c_id = 0;
	foreach my $c (@{$careful_clusters}){
		my $bx = prep_blastx_data($c, $c_id, $seq);
		push(@bx_data, @{$bx}) if defined $bx;

		#my $pp = prep_polpro_data($c, $c_id, $seq);
		#push(@pp_data, @{$pp}) if defined $pp;

                #my $pe = prep_polest_data($c, $c_id, $seq);
                #push(@pe_data, @{$pe}) if defined $pe;

		$c_id++;
	}
	return (\@bx_data, \@pp_data, \@pe_data);

}
#------------------------------------------------------------------------
sub segment_snaps {
        my $snaps            = shift;
        my $careful_clusters = shift;
        my $seq              = shift;


        my %clusters_hit;
	my %counts;
	my %spans;
	# initialize...

	foreach my $s (@{$snaps}){
		$counts{$s->{temp_id}} = 0;

		my ($sB, $sE) = PhatHit_utils::get_span_of_hit($s, 'query');
		   ($sB, $sE) = ($sE, $sB) if $sB > $sE;

		$spans{$s->{temp_id}} = [$sB, $sE];
	}

        print STDERR "Segmenting snap predctions\n";

	my $i = 0;
        foreach my $c (@{$careful_clusters}){
        	my $c_strand = $c->[0]->strand('query');
                
                my $coors  = PhatHit_utils::to_begin_and_end_coors($c, 'query');
                my $pieces = Shadower::getPieces($seq, $coors, 10);
                
		my $cB;
		my $cE;

                if (defined($pieces->[1])){	   
		   #die "error in auto_annotator::segment_snaps!\n";
		   $pieces = [sort {$a->{b} <=> $b->{b}} @{$pieces}];
		   $cB = $pieces->[0]->{b};
		   $cE = $pieces->[-1]->{e};
		}
		else{
		   $cB = $pieces->[0]->{b};
		   $cE = $pieces->[0]->{e};
		}

        	foreach my $s (@{$snaps}){

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

        foreach my $s (@{$snaps}){
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
sub perge_single_snaps {
   my $careful_clusters = shift;

   my @c_keepers;
   foreach my $c (@{$careful_clusters}){
      my $ests_in_cluster  = get_selected_types($c, 'est2genome');
      my $ps_in_cluster    = get_selected_types($c,'protein2genome');
      my $bx_in_cluster    = get_selected_types($c,'blastx');
      my $snaps_in_cluster = get_selected_types($c,'snap');

      if (defined($ests_in_cluster->[0]) ||
	  defined($ps_in_cluster->[0]) ||
	  defined($bx_in_cluster->[0])
	 ){
	    push (@c_keepers, $c);
      }
   }
   
   return (\@c_keepers);
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

        my $ests_in_cluster  = get_selected_types($c, 'est2genome');
        my $ps_in_cluster    = get_selected_types($c,'protein2genome');
        my $bx_in_cluster    = get_selected_types($c,'blastx');
	my $snaps_in_cluster = get_selected_types($c,'snap');

        # groups of most informative protein hits
	# go ahead and inclde the proteion2genome data as well... why not?
        my $gomiph = combine($ps_in_cluster, $bx_in_cluster);

        # group of most informative alt splices
        my $gomias = clean::get_best_alt_splices($ests_in_cluster, $seq, 10);

	my @data;
        if (defined($gomias->[0])){
        	foreach my $mia (@{$gomias}){
        	        push(@data, {'gomiph' => $gomiph,
                                     'ests'   => $ests_in_cluster,
				     'snaps'  => $snaps_in_cluster,
                                     'mia'    => $mia,
                                     'c_id'   => $c_id});

               	}
        }
        else {
        	push(@data, {'gomiph' => $gomiph,
			     'snaps'  => $snaps_in_cluster,
                             'ests'   => undef,
                             'mia'    => undef,
                             'c_id'   => $c_id});

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

        my $ests_in_cluster = get_selected_types($c, 'est2genome');
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


        my $ests_in_cluster = get_selected_types($c, 'est2genome');
        my $ps_in_cluster   = get_selected_types($c,'protein2genome');
        my $bx_in_cluster   = get_selected_types($c,'blastx');

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
#this subrutine returns finished MAKER annoations
sub annotate {
        my $virgin_fasta     = shift;
        my $masked_fasta     = shift;
	my $chunk_number     = shift; #required to name genes for each chunk
        my $exonerate_p_hits = shift;
        my $exonerate_e_hits = shift;
        my $blastx_hits      = shift;
	my $predictions      = shift;
        my $the_void         = shift;
	my $pred_command     = shift;
	my $pred_flank       = shift;
	my $single_exon      = shift;
	$OPT_F               = shift;
	$OPT_SNAPS           = shift;
	$OPT_PRED            = shift;

        my $def   = Fasta::getDef($masked_fasta);
        my $seq   = Fasta::getSeq($masked_fasta);
        my $v_seq = Fasta::getSeq($virgin_fasta);

	my ($bx_data, $pp_data, $pe_data) = prep_hits($exonerate_p_hits,
	                                              $exonerate_e_hits,
	                                              $blastx_hits,
						      $predictions,
	                                              $v_seq,
						      $single_exon
						     );

        my ($seq_id)  = $def =~ /^>(\S+)/;

        my @annotations;
        my $i = 0;
	my $g_name;

	my $bx_transcripts = 
	run_it($bx_data, $the_void, $seq, $v_seq, $def, 'bx', $pred_command, $pred_flank);


	#my $pp_transcripts = 
	#run_it($pp_data, $the_void, $seq, $v_seq, $def, 'pp', $snap_command, $snap_flank);

	#my $pe_transcripts = 
	#run_it($pe_data, $the_void, $seq, $v_seq, $def, 'pe', $snap_command, $snap_flank);
	
	#my @transcripts = (@{$bx_transcripts}, @{$pp_transcripts}, @{$pe_transcripts});

	my @transcripts = (@{$bx_transcripts});

        my $annotations = group_transcripts(\@transcripts, $v_seq, $seq_id, $chunk_number, $predictions);

        return $annotations;
}
#------------------------------------------------------------------------
sub run_it {
	my $data         = shift;
	my $the_void     = shift;
	my $seq          = shift;
	my $v_seq        = shift;
	my $def          = shift;
	my $id           = shift;
	my $pred_command = shift;
	my $pred_flank   = shift;

	my @transcripts;
	my $i = 0;
	foreach my $set (@{$data}){

		my $goimph = $set->{gomiph};
		my $mia    = $set->{mia};
		my $c_id   = $set->{c_id};
		my $ests   = $set->{ests};

		my ($snap_shots, $strand);

		($snap_shots, $strand)   = get_snap_shot($seq, 
		  				         $def, 
							 $id.'.'.$i, 
							 $the_void, 
							 $set, 
							 $pred_flank, 
							 $pred_command,
							 ) if $OPT_PRED eq 'snap';

                ($snap_shots, $strand)   = get_aug_shot($seq,
                                                        $def,
                                                        $id.'.'.$i,
                                                        $the_void,
                                                        $set,
                                                        $pred_flank,
                                                        $pred_command,
                                                        ) if $OPT_PRED eq 'augustus';


		my $best_pred       = get_best_snap_shot($strand, $snap_shots);
		my $on_right_strand = get_best_snap_shots($strand, $snap_shots);


		if (!defined($best_pred) && defined($mia)){
		    my $transcript = pneu($ests, $mia, $seq);
		    push(@transcripts, [$transcript, $set, undef]);
		}
		
		foreach my $snap_shot (@{$on_right_strand}){
		    my $copy = $snap_shot;
		    if    (defined($snap_shot) && defined($mia)){
			my $transcript = pneu($ests, $copy, $seq);
			if (defined($transcript)){
			    push(@transcripts, [$transcript, $set, $snap_shot])
			}
			else {
			    push(@transcripts, [$copy, $set, $snap_shot])
			}
		    }
		    elsif (defined($snap_shot) && !defined($mia)){
			push(@transcripts, [$copy, $set, $snap_shot])
			    if defined($copy);
		    }
		}
		$i++;
	}

	return \@transcripts;
}
#------------------------------------------------------------------------
sub load_transcript_struct {
	my $f            = shift;
	my $c_id         = shift;
	my $i            = shift;
	my $seq          = shift;
	my $seq_id       = shift;
	my $evi          = shift;
	my $snap_abinits = shift;

	my $transcript_seq  = get_transcript_seq($f, $seq);

        my ($translation_seq, $offset, $end) =
           get_translation_seq($transcript_seq);
	

	my $len_3_utr = length($transcript_seq) - $end + 1;
	my $l_trans =  length($translation_seq);
	
	my $qi = defined($evi)
	? maker::quality_index::get_transcript_qi($f,$evi,$offset,$len_3_utr,$snap_abinits, $l_trans)
	: 'non-overlapping-snap';

	my ($source) = ref($f) =~ /(\S+)\:\:PhatHit/;;

        my $t_name = "maker-$seq_id-gene-$source-$c_id-mRNA-$i"; #affects GFFV3.pm
           $t_name .= " $qi";

        my $t_struct = {'hit'      => $f,
                        't_seq'    => $transcript_seq,
                        'p_seq'    => $translation_seq,
                        't_offset' => $offset,
                        't_end'    => $end,
                        'c_id'     => $c_id,
		        't_name'   => $t_name,
			't_qi'       => $qi,
                    };

	return $t_struct;
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
	my $snap_abinits = shift;

	my @transcripts;
	my %lookup;
	my %snap_lookup;
	my $temp_id = 0;
	foreach my $datum (@{$data}){
		my $tra       = $datum->[0];
		my $set       = $datum->[1];
		my $o_snap    = $datum->[2]; # added 12-15
		
		$tra->{set_id} = $temp_id;

		push(@transcripts, $tra);

		$lookup{$temp_id} = $set;
		$temp_id++;
	}

	#-- add the non overlapping ab initio snap predictions
	if ($OPT_SNAPS){
	   my $non_overlapping_snap_abinits = get_non_overlapping(\@transcripts,
								  $snap_abinits
								 );
	   
	   push (@transcripts, @{$non_overlapping_snap_abinits});
	}
	#---
	
	my $careful_clusters = 
	cluster::careful_cluster_phat_hits(\@transcripts, $seq);

	my @keepers;
	foreach my $c (@{$careful_clusters}){
		my $best_alt_forms = 
		clean::remove_redundant_alt_splices($c, $seq, 10);
		push(@keepers, $best_alt_forms);
	}	
	my $c_id = 0;

	my %pred_sources;
	my @annotations;
	foreach my $c (@keepers){
		my @t_structs;
		my $i = 1;
		foreach my $f (@{$c}){

			my ($source) = ref($f) =~ /(\S+)\:\:PhatHit/;
			$pred_sources{$source}++;

			my $evidence = defined($f->{set_id}) ? $lookup{$f->{set_id}} : undef;

			my $t_struct = 
			load_transcript_struct($f, "$chunk_number.$c_id", $i, $seq, $seq_id, $evidence, $snap_abinits);
			push(@t_structs, $t_struct);

			$i++;
		}

		my ($g_start, $g_end, $g_strand) = get_start_and_end_on_seq(\@t_structs);

		my @sources = keys %pred_sources;
		die "ERROR: There is more than one source for annotations which may break GFF3V.pm\n" if (@sources > 1);
		my $source = $sources[0];

		my $annotation = { 't_structs' => \@t_structs,
		                   'g_name'    => "maker-gene-$source-.$chunk_number.$c_id", #affects GFFV3.pm
				   'g_start'   => $g_start,
			           'g_end'     => $g_end,
				   'g_strand'  => $g_strand,
				 };

		push(@annotations, $annotation);
		$c_id++;
	}

	return \@annotations;
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
	my @sorted_b = sort {$a->start('query') <=> $b->start('query')} @exons;
	my @sorted_e = sort {$b->end('query')   <=> $a->end('query')}   @exons;

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
			if (ref($f) =~ /$type$/){
				push(@keepers, $f);
			}
			elsif (ref($f) eq 'snap::PhatHit' && $type eq 'snap'){
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

	return ($best->[1], $best->[0]);

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

	my ($t_seq) = $p_seq_2 =~ /(^[^\*]+)\*?/;

        die "logic error in auto_annotate::get_off_and_str!\n"
        unless $t_seq  =~ /^M/;

        return ($offset, $t_seq);
}
#------------------------------------------------------------------------
sub get_translation_seq {
	my $seq    = shift;

	my $tM = new CGL::TranslationMachine();

	my ($p_seq , $offset) = $tM->longest_translation($seq);

	$p_seq =~ s/\*$//; # added 11/21/06

	my $end = length($p_seq)*3 + $offset + 1 + 3;

	my ($trim, $leader);
	if ($p_seq =~ /^M/){
		# easy  it begins with an M....
		return ($p_seq , $offset, $end );
	}
	else{
		# get the longest internal prot that begins with an M....
		my ($off_new, $p_seq_new) = get_longest_m_seq($seq);

		my $n_end = length($p_seq_new)*3 + $off_new + 1 + 3
		            if defined($p_seq_new);


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
sub get_best_snap_shot {
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
sub get_best_snap_shots {
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
	maker::join::join_f($b_5, $g, $b_3, $q_seq, $OPT_PRED);

	#return $g;
	return $anno_transcript;
}
#------------------------------------------------------------------------
sub get_snap_shot {
	my $seq           = shift;
	my $def           = shift;
	my $id            = shift;
	my $the_void      = shift;
	my $set           = shift;
	my $snap_flank    = shift;
	my $snap_command  = shift;
	

	my ($shadow_seq, $strand, $offset, $xdef, $use_full) =  
	    prep_for_genefinder($seq, $set, $snap_flank);

	my $shadow_fasta = Fasta::toFasta($def." $id offset:$offset", 
                                          $shadow_seq,
                                         );


	my ($exe, $hmm) = $snap_command =~ /(\S+)\s+(\S+)/;

	my $alt_snap_command;
        if ($strand == 1){
                $alt_snap_command = $exe.' -plus ';
        }
        else {
                 $alt_snap_command = $exe.' -minus ';
        }

	if ($use_full && -e $hmm.'.full_only'){
		$alt_snap_command .= $hmm.'.full_only';
	}
	else {
		$alt_snap_command .= $hmm;
	}
	
	my $gene_preds;

	$gene_preds = snap($$shadow_fasta,
                           $the_void,
                           $id,
	                   $strand,
		           $offset,
		           $xdef,
		           $alt_snap_command,
                          );


       $gene_preds = snap($$shadow_fasta,
                           $the_void,
                           $id,
                           $strand,
                           $offset,
                           $xdef,
                           $snap_command,
                          ) if $use_full && !defined($gene_preds->[0]);

	return ($gene_preds, $strand);

}
#------------------------------------------------------------------------
sub get_aug_shot {
        my $seq           = shift;
        my $def           = shift;
        my $id            = shift;
        my $the_void      = shift;
        my $set           = shift;
        my $pred_flank    = shift;
        my $pred_command  = shift;


        my ($shadow_seq, $strand, $offset, $xdef, $use_full) =
            prep_for_genefinder($seq, $set, $pred_flank);

        my $shadow_fasta = Fasta::toFasta($def." $id offset:$offset",
                                          $shadow_seq,
                                         );


        my $alt_snap_command;
        if ($strand == 1){
                $alt_snap_command = $pred_command.' --strand=forward';
        }
        else {
                 $alt_snap_command = $pred_command.' --strand=backward';
        }

        if ($use_full){
                $alt_snap_command .= $pred_command. '--genemodel=complete';
        }
        my $gene_preds;

        $gene_preds = augustus($$shadow_fasta,
                               $the_void,
                               $id,
                               $strand,
                               $offset,
                               $xdef,
                               $alt_snap_command,
                              );


       $gene_preds = augustus($$shadow_fasta,
                              $the_void,
                              $id,
                              $strand,
                              $offset,
                              $xdef,
                             $pred_command,
                             ) if $use_full && !defined($gene_preds->[0]);

        return ($gene_preds, $strand);

}
#------------------------------------------------------------------------
sub snap {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
	my $strand     = shift;
	my $offset     = shift;
	my $xdef       = shift;
	my $command    = shift;

        my $snap_keepers = [];

	my $file_name = "$the_void/$seq_id\.auto_annotator\.$offset\.snap.fasta";

	my $o_file    = "$the_void/$seq_id\.$offset\.auto_annotator\.snap";

	my $xdef_file = "$the_void/$seq_id\.$offset\.auto_annotator\.xdef\.snap";

	Widget::snap::write_xdef_file($xdef, $xdef_file);
	
                FastaFile::writeFile(\$fasta, $file_name);

                runSnap($file_name, $o_file, $xdef_file, $command);

                my %params;
                   $params{min_exon_score}  = -100;
                   $params{min_gene_score}  = -20;

                my $keepers =
                Widget::snap::parse($o_file,
                                    \%params,
                                    $file_name,
                                     );

                PhatHit_utils::add_offset($keepers,
                                          $offset,
                                          );

        return $keepers;

}
#------------------------------------------------------------------------
sub augustus {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
        my $strand     = shift;
        my $offset     = shift;
        my $xdef       = shift;
        my $command    = shift;

        my $aug_keepers = [];

        my $file_name = "$the_void/$seq_id\.auto_annotator\.$offset\.augustus.fasta";

        my $o_file    = "$the_void/$seq_id\.$offset\.auto_annotator\.augustus";

        my $xdef_file = "$the_void/$seq_id\.$offset\.auto_annotator\.xdef\.augustus";

        Widget::augustus::write_xdef_file($xdef, $xdef_file) if defined $xdef;

                FastaFile::writeFile(\$fasta, $file_name);

                runAugustus($file_name, $o_file, $xdef_file, $command);

                my %params;
                   $params{min_exon_score}  = -100;
                   $params{min_gene_score}  = -20;

                my $keepers =
                Widget::augustus::parse($o_file,
                                       \%params,
                                        $file_name,
                                        );

                PhatHit_utils::add_offset($keepers,
                                          $offset,
                                          );

        return $keepers;

}
#------------------------------------------------------------------------
sub runSnap {
        my $q_file    = shift;
        my $out_file  = shift;
	my $xdef_file = shift;
	my $command   = shift;

	$command .= " -xdef $xdef_file ";

    #$command .= " /usr/local/SNAP/HMM/planaria.f1.hmm -xdef $xdef_file ";
    #$command .= " /usr/local/SNAP/HMM/planaria.f1.hmm -xdef $xdef_file -ACoding 0.2 -AIntron 0.2 ";

            $command .= " $q_file";
            $command .= " > $out_file";

        my $w = new Widget::snap();

        if (-e $out_file && ! $OPT_F){
                print STDERR "re reading snap report.\n";
                print STDERR "$out_file\n";
        }
        else {
                print STDERR "running  snap.\n";
                $w->run($command);
        }

}
#------------------------------------------------------------------------
sub runAugustus {
        my $q_file    = shift;
        my $out_file  = shift;
        my $xdef_file = shift;
        my $command   = shift;

        $command .= " -xdef $xdef_file " if -e $xdef_file;

            $command .= " $q_file";
            $command .= " > $out_file";

        my $w = new Widget::augustus();

        if (-e $out_file && ! $OPT_F){
                print STDERR "re reading augustus report.\n";
                print STDERR "$out_file\n";
        }
        else {
                print STDERR "running  augustus.\n";
                $w->run($command);
        }

}
#------------------------------------------------------------------------
sub prep_for_genefinder {
	my $seq    = shift;
	my $set    = shift;
	my $flank  = shift;

	my $gomiph = $set->{gomiph};
	my $mia    = $set->{mia};
	my $snaps  = $set->{snaps};
	my $ests   = $set->{ests};

	my @t_data;

	push(@t_data, @{$gomiph})  if defined($gomiph);
	push(@t_data, @{$snaps})   if defined($snaps);
	push(@t_data, $mia)        if defined($mia);
	push(@t_data, @{$ests})    if defined($ests);

	my $p_set_coors = PhatHit_utils::get_hsp_coors($gomiph, 'query');

	my $n_set_coors = 
	defined($ests) ? PhatHit_utils::get_hsp_coors($ests, 'query')
                      : [];

	my @coors;
	my $plus  = 0;
	my $minus = 0;


	my $least;
	my $most;
        foreach my $hit (@t_data){
		foreach my $hsp ($hit->hsps()){
			my $s = $hsp->start('query');
			my $e = $hsp->end('query');   

			$least = $s if !defined($least) || $s < $least;
			$most  = $e if !defined($most)  || $e > $most;

			if ($hsp->strand('query') == 1) {
				$plus++;
			}
			else {
				$minus++;
			}
		}
        }

	my @span_coors = [$least, $most];

        my $p = Shadower::getPieces($seq, \@span_coors, $flank);

        my $final_seq = $p->[0]->{piece}; 

	my $offset    = $p->[0]->{b} - 1;
	#my $end_d     = length($$seq) - $p->[0]->{e};
	#my $size      = $p->[0]->{e} - $p->[0]->{b} + 1;

	#my $min_gene_size = 2000;
	#my $max_gene_size = 10000;

	my $use_full = 0;
	
	# 	$use_full = 1 if $size > $min_gene_size 
	# 	     && $size < $max_gene_size 
	# 	     && $offset > $min_gene_size 
	# 	     && $end_d > $min_gene_size; 
	
	
	my $strand = $plus > $minus ? 1 : -1;

	my $i_flank = 2;
	my $xdef;

        $xdef = Widget::snap::get_xdef($seq,
                                       $p_set_coors,
				       $n_set_coors,
                                       $strand == 1 ? '+' : '-',
                                       0.2, # Coding
                                       1000, # intron 
                                       $offset,
                                       $i_flank,
                                      ) if $OPT_PRED eq 'snap';

        $xdef = Widget::augustus::get_xdef($seq,
                                           $p_set_coors,
                                           $n_set_coors,
                                           $strand == 1 ? '+' : '-',
                                           0.2, # Coding
                                           1000, # intron
                                           $offset,
                                           $i_flank,
                                           ) if $OPT_PRED eq 'augustus';


	return (\$final_seq, $strand, $offset, $xdef, $use_full);
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


