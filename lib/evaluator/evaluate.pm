#------------------------------------------------------------------------
#----             	evaluator::evaluate.pm		             ---- 
#------------------------------------------------------------------------
package evaluator::evaluate;
use strict;
use vars qw/@ISA/;
use FileHandle;
use FastaChunker;
use PhatHit_utils;
use maker::auto_annotator;
use evaluator::funs;
use evaluator::scoring;
use evaluator::so_classifier;
use evaluator::AED;
use evaluator::fathom_utils;
use Fasta;


use vars qw/$OPT_F $OPT_PREDS $OPT_PREDICTOR $LOG $CTL/;
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
=pod
sub evaluate_for_gff {
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
        $OPT_PREDS           = shift;
        $OPT_PREDICTOR       = shift;
        $LOG                 = shift;

	my $seq = Fasta::getSeq($masked_fasta);
	my $def = Fasta::getDef($masked_fasta);
	my ($seq_id) = $def =~ /^>(\S+)/;

	##Get gene predictions here;
	my @genes;

	foreach my $gene (@genes) {
		my $temp_id = 1;

		my $so_code = evaluator::so_classifier::so_code($gene);
		my $alt_spli_sup = evaluator::funs::alt_spli($gene, 
			$exonerate_e_hits, $seq, $seq_id, $temp_id);
		$temp_id ++;
		
		foreach my $eat (@$gene) {
			
			my $evaluation = power_evaluate($eat, $seq,$exonerate_p_hits,
				$exonerate_e_hits, $blastx_hits, $predictions,$so_code, 
				$alt_spli_sup);
			
			##Output the evaluation here, or store it somewhere.
		}
	}
}
=cut
#-------------------------------------------------------------------------------
sub evaluate_maker_annotations {
    my $annotations = shift;
    my $seq         = shift;
    my $out_base    = shift;
    my $the_void    = shift;
    my $CTL_OPTIONS = shift;

    #iterate over each gene annotation
    foreach my $ann (@$annotations){
	my @t_cluster;
	my @pol_p_hits;
	my @pol_e_hits;
	my @blastx_hits;
	my @tblastx_hits;
	my @ab_inits;

	#collect all evidence and transcripts hit for the gene
	foreach my $t_struct (@{$ann->{t_structs}}){
	    my $hit = $t_struct->{hit};
	    my $evi = $t_struct->{evi};

	    my $pol_p   = maker::auto_annotator::get_selected_types($evi->{gomiph}, 'protein2genome');
	    my $pol_e   = maker::auto_annotator::get_selected_types($evi->{ests}, 'est2genome', 'est_gff');
	    my $pol_f   = maker::auto_annotator::get_selected_types($evi->{fusion}, 'est2genome', 'est_gff');
	    my $blastx  = maker::auto_annotator::get_selected_types($evi->{gomiph},'blastx', 'protein_gff');
	    my $tblastx = maker::auto_annotator::get_selected_types($evi->{alt_ests},'tblastx', 'altest_gff');
	    my $ab      = $evi->{all_preds};

	    push(@t_cluster, $hit);
	    push(@pol_p_hits, @$pol_p);
	    push(@pol_e_hits, @$pol_e, @$pol_f);
	    push(@blastx_hits, @$blastx);
	    push(@tblastx_hits, @$tblastx);
	    push(@ab_inits, @$ab);
	}

	#unique the set of evidence
	foreach my $hit (@pol_p_hits, @pol_e_hits, @blastx_hits, @tblastx_hits, @ab_inits){
	    $hit->{_uniq_set} = 0;
	}
	@pol_p_hits   = grep {$_->{_uniq_set}++ == 0} @pol_p_hits;
        @pol_e_hits   = grep {$_->{_uniq_set}++ == 0} @pol_e_hits;
        @blastx_hits  = grep {$_->{_uniq_set}++ == 0} @blastx_hits;
        @tblastx_hits = grep {$_->{_uniq_set}++ == 0} @tblastx_hits;
        @ab_inits     = grep {$_->{_uniq_set}++ == 0} @ab_inits;

	#get gene level statistics
	my $so_code = evaluator::so_classifier::so_code(\@t_cluster);
	my $alt_spli_sup = evaluator::funs::alt_spli(\@t_cluster, \@pol_e_hits, $seq);
	my $geneAED = evaluator::AED::gene_AED(\@t_cluster, \@pol_e_hits, \@pol_p_hits, \@blastx_hits, \@ab_inits, $seq);

	$ann->{geneAED} = $geneAED;
	$ann->{so_code} = $so_code;

	#get transcript level statistics;
	foreach my $t_struct (@{$ann->{t_structs}}){
	    my $t_name = $t_struct->{t_name};
	    my $f = $t_struct->{hit};

	    #run evaluator
	    my $eva = power_evaluate($f,
				     $seq,
				     \@pol_p_hits,
				     \@pol_e_hits,
				     \@blastx_hits,
				     \@ab_inits,
				     $so_code,
				     $geneAED,
				     $alt_spli_sup,
				     $t_name,
				     $ann->{g_name},
				     $CTL_OPTIONS,
				     $the_void,
				     );

	    $t_struct->{report} = $eva->{report};

	    #print report
	    my $dir = "$out_base/evaluator";
	    mkdir($dir) if(! -e $dir);	    
	    my $file = "$dir/".Fasta::seqID2SafeID($t_name);
	    ($file) = $file =~ /^([^\s\t\n]+)/;
	    $file .= ".eva";
	    open(my $FH, "> $file");
	    print $FH $eva->{report};
	    close($FH);
	}
    }

}
#-------------------------------------------------------------------------------
sub power_evaluate {
	my $eat = shift;
	my $seq = shift;
	my $exonerate_p_hits = shift;
	my $exonerate_e_hits = shift;
	my $blastx_hits	= shift;
	my $abinits_hits = shift;
	my $so_code = shift;
	my $geneAED = shift;
	my $alt_spli_sup = shift;
	my $t_name = shift;
	my $g_name = shift;
	$CTL = shift;
	my $the_void = shift;

	print STDERR "\nEVALUATing transcript $t_name...\n" unless $main::quiet;

	my $gff3_identity = get_identity($eat);

	$eat->{g_name} = $g_name;
	$eat->{t_name} = $t_name;

	my $box = evaluator::funs::prepare_box_for_gff($eat, $seq, 
						       $exonerate_p_hits, $exonerate_e_hits, 
						       $blastx_hits, $abinits_hits, $t_name);

	my ($txnAED, $closestAED) = evaluator::AED::txnAED($box, {'exon'=>1 });
						
	my $gene_length = abs($eat->nB('query') - $eat->nE('query')) + 1;

	my $snap_backwards = {overall_score => 'NA'};
	$snap_backwards = evaluator::fathom_utils::snap_backwards($box, $CTL
			, $the_void) if 
		$CTL->{enable_fathom} == 1;
	
	my $qi = evaluator::funs::get_transcript_qi($box);

	my $quality_seq	= evaluator::funs::get_quality_seq($box);

	my $splice_sites = evaluator::funs::get_splice_j_locations($eat);

	my $transcript_type = $box->{transcript_type};
	
	my $completion = evaluator::funs::completion($box);

	my $alt = $alt_spli_sup->{$eat->{_splice_form}};

	my $solexa_for_splices = evaluator::funs::solexa_support($CTL, $box,
			$splice_sites);

	my $score = scoring($qi, $quality_seq, $so_code, $transcript_type,
			    $completion, $alt, $solexa_for_splices);


	my $report = generate_report($eat, $box, $qi, $quality_seq, $splice_sites,
				     $transcript_type, $completion, $alt, $score, 
					$so_code, $geneAED, $txnAED, $closestAED, 
					$solexa_for_splices, $gff3_identity,
					$snap_backwards, $gene_length);

	print STDERR "Finished.\n\n" unless $main::quiet;

        my $eva = { 'qi'                => $qi,
                    'score'             => $score,
                    'quality_seq'       => $quality_seq,
                    'completion'        => $completion,
                    'alt'               => $alt,
                    'transcript_type'   => $transcript_type,
		    'report'		=> $report,
		    'so_code'		=> $so_code,
		    'txnAED'		=> $txnAED,
		    'closestAED'	=> $closestAED,
		    'gene_AED'          => $geneAED,
		    'snap_backwards'    => $snap_backwards->{overall_score},
		    'gff3_identity'	=> $gff3_identity,
		    'gene_length'	=> $gene_length,
                  };

	return $eva;
}


#-------------------------------------------------------------------------------
sub generate_report {
	my $eat 		= shift;
	my $box			= shift;
	my $qi  		= shift;
	my $quality_seq		= shift;
	my $splice_sites	= shift;
	my $type		= shift;
	my $completion		= shift;
	my $alt			= shift;
	my $score		= shift;
	my $so_code		= shift;
	my $geneAED             = shift;
	my $txnAED		= shift;
	my $closestAED		= shift;
	my $solexa		= shift;
	my $gff3_identity	= shift;
	my $snap_backwards	= shift;
	my $gene_length		= shift;

	my $g_name = $eat->{g_name};
	my $t_name = $eat->{t_name};
	my $prefix = $t_name.'@'.$g_name;
	($qi) = $qi =~ /QI:(.+)$/;

	my $report;

	$report .=
	"##################################################################\n";
	$report .= "#####          $t_name";
	my $space_number = (46- length($t_name));
	for(my $temp =1; $temp<= $space_number; $temp++) { $report .= ' ';}

	$report .= "#####\n";

	$report .=
	"##################################################################\n";

	$report .= $prefix."\t"."Quality_Index"."\t";
	$report .= $qi."\n";

	$report .= $prefix."\t"."translation_offset"."\t";
	$report .= $box->{translation_offset}."\n";
	
	$report .= $prefix."\t"."translation_end"."\t";
	$report .= $box->{translation_end}."\n";

	$report .= $prefix."\t"."translation_length"."\t";
	$report .= $box->{translational_length}."\n";

	$report .= $prefix."\t"."gene_length"."\t";
	$report .= $gene_length."\n";

	$report .= $prefix."\t"."splice_sites"."\t";
	foreach my $junction (@$splice_sites) {
		$report .= "$junction,";
	}
	chop $report if $report =~ /,$/;
	$report .= "\n";

	$report .= $prefix."\t"."completion_of_protein"."\t";
	$report .= $completion."\n";

	$report .= $prefix."\t"."transcript_type"."\t";
	$report .= (defined $type)? ($type."\n"):('unknown'."\n");

	$report .= $prefix."\t"."SO_CODE"."\t";
	$report .= $so_code."\n";

	$report .= $prefix."\t"."geneAED"."\t";
	$report .= $geneAED."\n";

        $report .= $prefix."\t"."txnAED"."\t";
        $report .= $txnAED."\n";

        $report .= $prefix."\t"."closestAED"."\t";
        $report .= $closestAED."\n";
	
	$report .= $prefix."\t"."support_for_altsplicing"."\t";
	$report .= $alt."\n";

	$report .= $prefix."\t"."splice_support_by_short_ests"."\t";
	if (!defined $solexa) 	{ $report .= "NA\n";}
	else {
		foreach my $site (sort {$a <=> $b} (keys %$solexa)) {
			$report .= ($site.'='.$solexa->{$site}.',');
		}
		chop $report;
		$report .= "\n";
	}


	$report .= $prefix."\t"."gff3_identity"."\t";
	$report .= $gff3_identity->[0].';'.$gff3_identity->[1].';'.
			$gff3_identity->[2].';'.$gff3_identity->[3]."\n";

	$report .= $prefix."\t"."snap_backwards_score"."\t";
	$report .= $snap_backwards->{overall_score}."\n";

	$report .= $prefix."\t"."quality_sequence"."\t";
	$report .= $$quality_seq."\n";

	return $report;
}
	
#-------------------------------------------------------------------------------
sub scoring {
	my $qi 		= shift;
	my $quality_seq = shift;
	my $so_code	= shift;
	my $type	= shift;
	my $completion	= shift;
	my $alt 	= shift;
	my $solexa_sup  = shift;  ## If doesn't allow to calculate solexa 
				  ## support, this will be undef.

	($qi) = $qi =~ /QI:(.+)$/;
	my @qi = split /\|/, $qi;
	my @so_code = split /:/, $so_code;
	my @quality_seq = split /\|/, $$quality_seq;

	my $junk_score = 0;
	my $count = 0;

	foreach my $one (@quality_seq) {
		next unless $one =~ /^(\d|-|\.)+$/;
		$junk_score += $one;
		$count ++;
	}
	$junk_score = ($count==0) ? 0 : ($junk_score/$count);

	my %arg;
	$arg{QI} = \@qi;
	$arg{SO_CODE} = \@so_code;
	$arg{type} = $type;
	$arg{junk} = $junk_score;
	$arg{alt_splicing} = $alt eq 'NA' ? 0: $alt ;
	$arg{completion} = $completion;

	my $score = evaluator::scoring::valuecal(%arg);

	return $score;
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_identity {
	my $eat = shift;

	my $gff3_g_name = 'UNKNOWN';
	my $gff3_g_id = 'UNKNOWN';
	my $gff3_t_name = 'UNKNOWN';
	my $gff3_t_id = 'UNKNOWN';

	if (defined $eat->{_description}) {
		my $content = $eat->{_description};

		my @items = split /;/, $content;
		my %hash;

		foreach my $item (@items) {
			my ($first, $second) = split /=/, $item;
			$hash{$first} = $second if defined $second;
		}

		$gff3_g_name = $hash{g_name} if defined $hash{g_name};
		$gff3_g_id = $hash{g_id} if defined $hash{g_id};
		$gff3_t_name = $hash{t_name} if defined $hash{t_name};
		$gff3_t_id = $hash{t_id} if defined $hash{t_id};
	}

	return [$gff3_g_name, $gff3_g_id, $gff3_t_name, $gff3_t_id];
}
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
1;


