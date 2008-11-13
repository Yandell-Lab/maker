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
use Fasta;
@ISA = qw(
       );
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------

sub evaluate_in_maker {

	my  ($f, $c_id, $i, $seq,
             $seq_id, $evi, $abinits, $so_code, $OPT_PREDICTOR, 
             $transcript_seq, $translation_seq, $offset, $end,
             $len_3_utr, $l_trans, $alt_spli_sup)  = @_;

	my $box = evaluator::funs::prepare_box_in_maker($f, $evi, 
		$len_3_utr, $abinits, $l_trans, $seq, \$transcript_seq,
		\$translation_seq, $offset, $end);

	## Get QI; 
        my $qi = defined($evi)
        ? evaluator::funs::get_transcript_qi($box)
        : 'non-overlapping-'.$OPT_PREDICTOR;

	## Get quality sequence.
	my $quality_seq = evaluator::funs::get_quality_seq($box);

	## Get splice site location.
	my $splice_sites = evaluator::funs::get_splice_j_locations($f);

	## Get transcript type (currently always 'mRNA' for maker).
	my $transcript_type = $box->{transcript_type};

	## Get the completion info for the transcript.
	my $completion = evaluator::funs::completion($box);

	## Get the alternative splicing support.
	my $alt = $alt_spli_sup->{$f->{_tran_id}};

	## Get the score for this struct.
	my $score = scoring($qi, $quality_seq, $so_code, $transcript_type, 
		$completion, $alt);


	my $eva = { 'qi'		=> $qi,
		    'score'		=> $score,
		    'quality_seq'	=> $quality_seq,
		    'completion'	=> $completion,
		    'alt'		=> $alt,
		    'transcript_type'	=> $transcript_type,
		    
		  };

        return $eva;
}


#-------------------------------------------------------------------------------
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
				$exonerate_e_hits, $blastx_hits, $so_code, 
				$alt_spli_sup);
			
			##Output the evaluation here, or store it somewhere.
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
	my $so_code = shift;
	my $alt_spli_sup = shift;

	my $box = evaluator::funs::prepare_box_for_gff($eat, $seq, 
		$exonerate_p_hits, $exonerate_e_hits, $blastx_hits);


=cut
#-------------------------------------------------------------------------------
sub scoring {
	my $qi 		= shift;
	my $quality_seq = shift;
	my $so_code	= shift;
	my $type	= shift;
	my $completion	= shift;
	my $alt 	= shift;

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

#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
1;


