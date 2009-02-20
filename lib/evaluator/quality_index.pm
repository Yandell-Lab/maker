#------------------------------------------------------------------------
#----                          evaluator::quality_index              ----
#------------------------------------------------------------------------
package evaluator::quality_index;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use clean;
use maker::auto_annotator;
use CGL::TranslationMachine;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
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
sub get_length_protein {
        my $eat = shift;

        my $seq      = $eat->seq();
        my $t_offset = $eat->translation_offset();
        my $t_end    = $eat->translation_end();


	my $length = int(($t_end - $t_offset)/3);

	return $length;
}
#------------------------------------------------------------------------
sub get_length_of_utrs {
	my $eat = shift;

	my $seq      = $eat->seq();
	my $t_offset = $eat->translation_offset();
	my $t_end    = $eat->translation_end();	

	my $length_5 = $t_offset;
	my $length_3 = length($seq) - $t_end;

	return ($length_5, $length_3);
}
#------------------------------------------------------------------------
sub get_transcript_qi {
        my $t            = shift;
	my $blastx       = shift;
	my $pol_ests     = shift;
	my $pol_pro      = shift; 
	my $snap_abinits = shift;


	my $l_trans = get_length_protein($t);

	my ($length_5, $length_3) = get_length_of_utrs($t);

	my @bag;
	push(@bag, @{$blastx})     if defined($blastx->[0]);
	push(@bag, @{$pol_ests})   if defined($pol_ests->[0]);
	push(@bag, @{$pol_pro})    if defined($pol_pro->[0]);

	my $pol_est_hits = maker::auto_annotator::get_selected_type(\@bag, 'est2genome', 'est_gff');
	my $pol_pro_hits = maker::auto_annotator::get_selected_type(\@bag, 'protein2genome');

	my @good_splicers;
	push(@good_splicers, @{$pol_est_hits});
	push(@good_splicers, @{$pol_pro_hits});	

	my $qi = 'QI:';
	# length 5-prime utr
	$qi .= $length_5;

	my $num_e = $t->num_hsps();
	my $num_i = $num_e - 1;

	#  fra. of splices confirmed
	my $num_confirmed  = $num_i == 0 ? -1 : confirm_splices($t, \@good_splicers);
	my $frac_confirmed = $num_i == 0 ? -1 : substr($num_confirmed/$num_i, 0, 4);

	$qi .= "|$frac_confirmed";

	# frac. exons perfectly tagged by an est2genome or protein2genome hit
	$num_confirmed  = confirm_exons($t, \@good_splicers);
	$frac_confirmed = substr($num_confirmed/$num_e, 0, 4);

	$qi .= "|$frac_confirmed";

	# frac. exons overlapped by anything
        $num_confirmed = confirm_overlaps($t, \@bag, '[1AaBbZzIi]');
	$frac_confirmed = substr($num_confirmed/$num_e, 0, 4);
	
	$qi .= "|$frac_confirmed";

        #  fra. of splices confirmed confirmed by a snap abinitio prediction
        $num_confirmed  = $num_i == 0 ? -1 : confirm_splices($t, $snap_abinits);
        $frac_confirmed = $num_i == 0 ? -1 : substr($num_confirmed/$num_i, 0, 4);

	$qi .= "|$frac_confirmed";

        # frac. exons perfectly tagged by a snap abinitio prediction 
        $num_confirmed  = confirm_snaps($t, $snap_abinits);
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
1;


