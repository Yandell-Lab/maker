#------------------------------------------------------------------------
#----                          evaluator::so_classifier              ----
#------------------------------------------------------------------------
package evaluator::so_classifier;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use clean;
use CGL::TranslationMachine;
use Shadower;
@ISA = qw(
       );

our %MATRIX;
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub so_code {
	my $eats     = shift;

	my $num = @{$eats};

	return '0:0:0' unless $num > 1;

	my %code;
	$code{X} = 0;
	$code{Y} = 0;
	$code{Z} = 0;

	I: for (my $i = 0; $i < @{$eats} - 1; $i++){
		my @hsps_i   = $eats->[$i]->hsps();
		my $strand_i = $eats->[$i]->strand('query');

		J: for (my $j = $i + 1; $j < @{$eats}; $j++){
			my @hsps_j   = $eats->[$j]->hsps();
			my $strand_j = $eats->[$j]->strand('query');

			if ($strand_i != $strand_j){
				$code{X}++;
				next J;
			}
			

			my $type = compare_hsps(\@hsps_i, \@hsps_j);

			if    ($type eq 'identical_exons'){
				$code{Z}++;
			}
			elsif ($type eq 'overlapping_exons'){
				$code{Y}++;
			}
			elsif ($type eq 'no_overlap'){
				$code{X}++;
			} 
		}
	}

	return $code{X}.':'.$code{Y}.':'.$code{Z};
}
#------------------------------------------------------------------------
sub compare_hsps {
	my $hsps_i = shift;
	my $hsps_j = shift;

	my $overlap = 0;
	foreach my $hsp_i (@{$hsps_i}){
		my $iB = $hsp_i->nB('query');
		my $iE = $hsp_i->nE('query');
		foreach my $hsp_j (@{$hsps_j}){
			my $jB = $hsp_j->nB('query');
			my $jE = $hsp_j->nE('query');

			my $class = compare::compare($iB, $iE, $jB, $jE);

			return 'identical_exons' if $class eq '1';
			
			if ($class ne '0'){
				$overlap = 1;
			}
		}
	}
	if ($overlap == 1){
		return 'overlapping_exons';
	}
	else {
		return 'no_overlap';
	}
}
#------------------------------------------------------------------------
1;


