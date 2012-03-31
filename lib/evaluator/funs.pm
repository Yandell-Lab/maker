#------------------------------------------------------------------------
#----                      evaluator::funs	                  ----
#------------------------------------------------------------------------
package evaluator::funs;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use clean;
use evaluator::quality_seq;
use maker::auto_annotator;
use compare;
use evaluator::Widget::blast;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub prepare_box_in_maker {
        my $t            = shift;
        my $set          = shift;
        my $length_3     = shift;
        my $snap_abinits = shift;
        my $l_trans      = shift;
	my $seq		 = shift;
	my $transcription_seq = shift;
	my $translation_seq  = shift;
	my $translation_offset = shift;
	my $translation_end    = shift;

	my $length_5     = $translation_offset;

	my @bag;
	push(@bag, @{$set->{gomiph}}) if defined($set->{gomiph});
        push(@bag, @{$set->{ests}})   if defined($set->{ests});
        push(@bag, @{$set->{fusion}}) if defined($set->{fusion});

        my $pol_est_hits = maker::auto_annotator::get_selected_type(\@bag, 'est2genome', 'est_gff');
        my $pol_pro_hits = maker::auto_annotator::get_selected_type(\@bag, 'protein2genome');
	my $blastx_hits  = maker::auto_annotator::get_selected_type(\@bag, 'blastx', 'protein_gff');
	
	my @good_splicers;
        push(@good_splicers, @{$pol_est_hits});
        push(@good_splicers, @{$pol_pro_hits}); 

	my $transcript_type = 'mRNA'; ## Right now maker can only 
				      ## mRNA genes!

	my %box = (
		'transcript' => $t,
		'evidence'=> $set,
		'bag'	=> \@bag,
		'good_splicers' => \@good_splicers,
		'abinits'  => $snap_abinits,
		'5_len'		=> $length_5,
		'3_len'		=> $length_3,
		'seq'		=> $seq,
		'translational_length' => $l_trans,
		'gomiph'	=> $set->{gomiph},
		'ests'		=> $pol_est_hits,
		'exonerate'	=> $pol_pro_hits,
		'blastx'	=> $blastx_hits,
		'transcription_seq'=> $transcription_seq,
		'translation_seq' => $translation_seq,
		'transcript_type' => $transcript_type,
		'translation_offset' => $translation_offset,
		'translation_end'    => $translation_end,
	);

	return \%box;
}

#------------------------------------------------------------------------
sub prepare_box_for_gff {
	my $eat = shift;
	my $seq = shift;
	my $exonerate_p_hits = shift;
	my $exonerate_e_hits = shift;
	my $blastx_hits = shift;
	my $abinits_hits = shift;
	my $t_name = shift;

	my $good_pol_p_hits = maker::auto_annotator::get_overlapping_hits
		($eat, $exonerate_p_hits);
	my $good_pol_e_hits = maker::auto_annotator::get_overlapping_hits
		($eat, $exonerate_e_hits);
	my $good_blastx_hits= maker::auto_annotator::get_overlapping_hits
		($eat, $blastx_hits);
	my $good_abinits_hits = maker::auto_annotator::get_overlapping_hits
		($eat, $abinits_hits);

	my @bag;
	push @bag, @$good_pol_p_hits;
	push @bag, @$good_pol_e_hits;
	push @bag, @$good_blastx_hits;

	my @good_splicers;
	push @good_splicers, @$good_pol_p_hits;
	push @good_splicers, @$good_pol_e_hits;

	my $transcription_seq = 
		maker::auto_annotator::get_transcript_seq($eat, $seq);
	my ($translation_seq, $offset, $end) =
		maker::auto_annotator::get_translation_seq($transcription_seq);

	my $len_3_utr = length($transcription_seq) - $end +1;
	my $l_trans = length($translation_seq);

	my $transcript_type = $eat->{type}; ## to be implemented

	my %box = (
		't_name'		=> $t_name,
		'transcript'		=> $eat,
		'bag'			=> \@bag,
		'good_splicers'		=> \@good_splicers,
		'5_len'			=> $offset,
		'3_len'			=> $len_3_utr,
		'seq'			=> $seq,
		'translational_length'	=> $l_trans,
		'est'			=> $good_pol_e_hits,
		'exonerate'		=> $good_pol_p_hits,
		'blastx'		=> $good_blastx_hits,
		'abinits'		=> $good_abinits_hits,
		'transcription_seq'	=> \$transcription_seq,
		'translation_seq'	=> \$translation_seq,
		'transcript_type'	=> $transcript_type,
		'translation_offset'	=> $offset,
		'translation_end'	=> $end,
	);
	
	return \%box;
}


#------------------------------------------------------------------------
sub solexa_support {
	my $CTL = shift;
	my $box = shift;
	my $splice_sites = shift;

	return undef unless ($CTL->{est_reads} &&
		scalar @$splice_sites > 0);

	my $transcript_seq = $box->{transcription_seq};
	my $path = $CTL->{current_tmp_path};
	my $cpus = $CTL->{cpus};
	my $side_thre = $CTL->{side_thre};
	my $est_file = $CTL->{est_reads};
	my $percov = $CTL->{eva_percov_blastn};
	my $percid = $CTL->{eva_percid_blastn};
	my $expection = $CTL->{eva_eval_blastn};
	my $bit = $CTL->{eva_bit_blastn};
	my $t_name = $box->{t_name};
	my $window_size = $CTL->{eva_window_size};	
	my $blastn = $CTL->{blastn};
	my $hspmax = $CTL->{eva_hspmax};
	my $gspmax = $CTL->{eva_gspmax};
	my $split_hit = $CTL->{eva_split_hit};
	my $xdformat = $CTL->{xdformat};

	$path = $path.'/'.'evalutor';
	unless (-d $path) {
		`mkdir $path`;
	}

	unless ((-e $est_file.'.xns') && (-e $est_file.'.xnt') 
			&& (-e $est_file.'.xnd') ) {
		`$xdformat -n $est_file`;
	}

        unless ((-e $est_file.'.xns') && (-e $est_file.'.xnt')
                        && (-e $est_file.'.xnd') ) {
		my $new_est_file = $path.'/'.'evaluator_solexa_reads_master_db.fasta';
		`cp $est_file $new_est_file`;
		$CTL->{est_reads} = $new_est_file;
		$est_file = $new_est_file;
		`$xdformat -n $est_file`;
	}

	my $solexa_support = evaluator::Widget::blast::short_est_support(
		'path' => $path,
		'cpus' => $cpus,
		'seq'  => $transcript_seq,
		'splice_sites' => $splice_sites,
		'db' => $est_file,
		't_name' => $t_name,
		'window_size' => $window_size,
		'expection' => $expection,
		'blastn' => $blastn,
		'side_thre' => $side_thre,
		'hsp_bit_min' => $bit,
		'percov' => $percov,
		'percid' => $percid,
		'split_hit' => $split_hit,
		'hspmax' => $hspmax,
		'gspmax' => $gspmax,
	);
	
	return $solexa_support;
}

#------------------------------------------------------------------------
sub get_code_info {
        my $t    = shift;
        my $hits = shift;

        my @codes;
        foreach my $hit (@{$hits}){
                my $code = compare::compare_phat_hits($t, $hit, 'query');
                push(@codes, $code);
        }

        return \@codes;
}
#------------------------------------------------------------------------
sub get_exon_coors {
        my $hit = shift;

        my $sorted = PhatHit_utils::sort_hits($hit, 'query');

        my @exon_coors;
        for (my $i = 0; $i < @{$sorted}; $i++){

                my $e_b = $sorted->[$i]->nB('query');
                my $e_e = $sorted->[$i]->nE('query');

                push(@exon_coors, [$e_b, $e_e]);
        }

        return \@exon_coors;
}
#------------------------------------------------------------------------
sub get_intron_coors {
	my $hit = shift;

	my $sorted = PhatHit_utils::sort_hits($hit, 'query');

	my @intron_coors;
	for (my $i = 1; $i < @{$sorted}; $i++){

		my $i_b = $sorted->[$i-1]->nE('query');
		my $i_e = $sorted->[$i]->nB('query');

		push(@intron_coors, [$i_b, $i_e]);
	}

	return \@intron_coors;
}
#------------------------------------------------------------------------
sub confirm_overlaps {
        my $t    = shift;
        my $hits = shift;
	my $expr = shift;

        my $exon_coors_t = get_exon_coors($t);

        my @all_supporting_data_exon_coors;
        foreach my $hit (@{$hits}){
                my $exon_coors_h = get_exon_coors($hit);
                push(@all_supporting_data_exon_coors, $exon_coors_h);
        }

        my $num_confirmed = 0;
        foreach my $t_e (@{$exon_coors_t}){
                $num_confirmed++
                if overlaps($t_e, \@all_supporting_data_exon_coors, $expr);
        }
        return $num_confirmed;
}
#------------------------------------------------------------------------
sub confirm_snaps {
        my $t    = shift;
        my $hits = shift;

        my $exon_coors_t = get_exon_coors($t);

        my @all_supporting_data_exon_coors;
        foreach my $hit (@{$hits}){
                my $exon_coors_h = get_exon_coors($hit);
                push(@all_supporting_data_exon_coors, $exon_coors_h);
        }

        my $num_confirmed = 0;

	if ($t->num_hsps() > 1){
		my $alpha = shift(@{$exon_coors_t});
        	my $omega =   pop(@{$exon_coors_t});

        	$num_confirmed++
        	if overlaps($alpha, \@all_supporting_data_exon_coors, '[A1]');

        	$num_confirmed++
        	if overlaps($omega, \@all_supporting_data_exon_coors, '[a1]');

	}
	elsif ($t->num_hsps() == 1){
		my $a_o = shift(@{$exon_coors_t});

		$num_confirmed++
                if overlaps($a_o, \@all_supporting_data_exon_coors, '[Ai1a]');
	}

	foreach my $t_e (@{$exon_coors_t}){
                $num_confirmed++
                if confirm($t_e, \@all_supporting_data_exon_coors);
        }
        return $num_confirmed;
}
#------------------------------------------------------------------------
sub confirm_exons {
        my $t    = shift;
        my $hits = shift;

        my $exon_coors_t = get_exon_coors($t);

        my @all_supporting_data_exon_coors;
        foreach my $hit (@{$hits}){
                my $exon_coors_h = get_exon_coors($hit);
                push(@all_supporting_data_exon_coors, $exon_coors_h);
        }

        my $num_confirmed = 0;
        foreach my $t_e (@{$exon_coors_t}){
                $num_confirmed++
                if confirm($t_e, \@all_supporting_data_exon_coors);
        }
        return $num_confirmed;
}
#------------------------------------------------------------------------
sub confirm {
        my $t_i   = shift;
        my $asdic = shift;

        my $t_b = $t_i->[0];
        my $t_e = $t_i->[1];
        foreach my $i_set (@{$asdic}){
                foreach my $s_i (@{$i_set}){
                        my $s_b = $s_i->[0];
                        my $s_e = $s_i->[1];

                        return 1 if $t_b == $s_b && $t_e == $s_e;
                }
        }
        return 0;
}
#------------------------------------------------------------------------
sub confirm_splices {
	my $t    = shift;
	my $hits = shift;

	my $intron_coors_t = get_intron_coors($t);

	my @all_supporting_data_intron_coors;
	foreach my $hit (@{$hits}){
		next unless $hit->num_hsps() > 1;
		my $intron_coors_h = get_intron_coors($hit);
		push(@all_supporting_data_intron_coors, $intron_coors_h);
	}

	my $num_confirmed = 0;
	foreach my $t_i (@{$intron_coors_t}){
		$num_confirmed++ 
		if confirm($t_i, \@all_supporting_data_intron_coors);
	}
	return $num_confirmed;
}
#------------------------------------------------------------------------
sub overlaps {
        my $t_i   = shift;
        my $asdic = shift;     
	my $expr  = shift;
        my $t_b = $t_i->[0];
        my $t_e = $t_i->[1];
        foreach my $i_set (@{$asdic}){
                foreach my $s_i (@{$i_set}){
                        my $s_b = $s_i->[0];
                        my $s_e = $s_i->[1];

			my $code = compare::compare($t_b, $t_e, $s_b, $s_e);
                        return 1 if $code =~ /$expr/;
                }
        }
        return 0;
}
#------------------------------------------------------------------------
sub get_percent {
	my $codes       = shift;
	my $num_t_exons = shift;
	my $wanted      = shift;

	my $largest = 0;
	foreach my $code (@{$codes}){
		my @stuff = split('', $code);

		my $num = 0;
		foreach my $e (@stuff){
			$num++ if $e =~ /$wanted/;
		}
		
		$largest = $num if $num > $largest;

	}
	return substr($largest/$num_t_exons, 0, 5);
}
#------------------------------------------------------------------------
sub get_transcript_qi {
	my $box = shift;

	my $t = 		$box->{transcript};
	my $good_splicers =     $box->{good_splicers};
	my $length_5	 = 	$box->{'5_len'};
	my $length_3	 = 	$box->{'3_len'};
	my $abinits = 		$box->{abinits};
	my $l_trans =		$box->{translational_length};
	my $bag		=       $box->{bag};

	my $qi = 'QI:';
	# length 5-prime utr
	$qi .= $length_5;

	my $num_e = $t->num_hsps();
	my $num_i = $num_e - 1;

	#  fra. of splices confirmed
	my $num_confirmed  = $num_i == 0 ? -1 : confirm_splices($t, $good_splicers);
	my $frac_confirmed = $num_i == 0 ? -1 : substr($num_confirmed/$num_i, 0, 4);

	$qi .= "|$frac_confirmed";

	# frac. exons perfectly tagged by an est2genome or protein2genome hit
	$num_confirmed  = confirm_exons($t, $good_splicers);
	$frac_confirmed = substr($num_confirmed/$num_e, 0, 4);

	$qi .= "|$frac_confirmed";

	# frac. exons overlapped by anything
        $num_confirmed = confirm_overlaps($t, $bag, '[1AaBbZzIi]');
	$frac_confirmed = substr($num_confirmed/$num_e, 0, 4);
	
	$qi .= "|$frac_confirmed";

        #  fra. of splices confirmed confirmed by a snap abinitio prediction
        $num_confirmed  = $num_i == 0 ? -1 : confirm_splices($t, $abinits);
        $frac_confirmed = $num_i == 0 ? -1 : substr($num_confirmed/$num_i, 0, 4);

	$qi .= "|$frac_confirmed";

        # frac. exons perfectly tagged by a snap abinitio prediction 
        $num_confirmed  = confirm_snaps($t, $abinits);
        $frac_confirmed = substr($num_confirmed/$num_e, 0, 4);

	$qi .= "|$frac_confirmed";

        # num HSPs
        $qi .= "|$num_e";

	# length 3-prime utr
	$qi .= "|$length_3";
	#length of translation
	
	$qi .= "|$l_trans";
	return $qi;
}
#------------------------------------------------------------------------
sub get_quality_seq {
	my $box = shift;

	my $ests   = $box->{ests};
	my $abinits= $box->{abinits};
	my $t	   = $box->{transcript};
	my $seq	   = $box->{seq};
	my $exonerate = $box->{exonerate};
	my $blastx = $box->{blastx};

	my $quality_seq = evaluator::quality_seq::prep($t,
						       $ests,
						       $exonerate,
						       $blastx,
						       $abinits,
						       $seq,);

	return $quality_seq;
}
#------------------------------------------------------------------------
sub completion {
	my $box = shift;

	my $seq = $box->{transcription_seq};
	my $length = length($$seq);

	my $first_codon = uc(substr($$seq, $box->{translation_offset}, 3));
	my $last_codon  = uc(substr($$seq, $box->{translation_end}-4, 3));

	my ($s_codon, $e_codon)=('0','0');

	$s_codon = '1' if $first_codon eq 'ATG';

	foreach (qw/TAG TAA TGA/) {
                $e_codon = '1' if $last_codon eq $_;
        }

        return join '', ($s_codon, $e_codon);
}
		
#------------------------------------------------------------------------
sub get_splice_j_locations {
        my $eat = shift;

        my @hsps = $eat->hsps();

        shift(@hsps); # don't want the first one....
        
        my @locs;
        foreach my $hsp (@hsps){
                push(@locs, $hsp->nB('hit'));
        }
        return \@locs;
}
#------------------------------------------------------------------------
sub alt_spli {
    my $eats           = shift;
    my $pol_e_hits     = shift;
    my $seq            = shift;

    
    #pre-label alt splice forms
    my $i = 1;
    foreach my $eat (@{$eats}){
	$eat->{_splice_form} = $i;
	$i++;
    }
    
    #find same alt_splice forms
    my %how_many;
    my %which_one;
    my $temp_id = 0;
    foreach my $est (@{$pol_e_hits}) {
	$est->temp_id($temp_id);
	
	foreach my $eat (@{$eats}){
	    if(compare::hsps_overlap($eat, $est, 'query', 0) &&
	       compare::is_same_alt_form($eat, $est, $seq, 0)
	       ){
		$how_many{$est->temp_id}++;
		push(@{$which_one{$eat->{_splice_form}}}, $est->temp_id); 
	    }
	}
    }
    
    #quantify the alt_splice support
    my $num_eats = @{$eats};    
    my %is_supported;
    foreach my $eat (@{$eats}){
	if($num_eats == 1){
	    $is_supported{$eat->{_splice_form}} = 'NA';
	    next;
	}
	else{
	    $is_supported{$eat->{_splice_form}} = 0;
	}

	foreach my $id (@{$which_one{$eat->{_splice_form}}}){
	    my $how_many = $how_many{$id};
	    
	    $is_supported{$eat->{_splice_form}} += 
		($num_eats - $how_many)/($num_eats - 1);
	}
    } 
    
    return \%is_supported;
}

#------------------------------------------------------------------------
1;


