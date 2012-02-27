#------------------------------------------------------------------------
#----                        evaluator::gff3_to_phatHit              ----
#------------------------------------------------------------------------
package evaluator::gff3_to_phatHit::gff3_to_phatHit;
use strict;
use warnings;
use evaluator::gff3_to_phatHit::FlyBase;
use evaluator::gff3_to_phatHit::WormBase;

#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub new {

	my $class = shift;
	my $gff3_file  = shift;
	my $gff3_type  = shift;
	my $fasta_file = shift;
	my $ids	       = shift;

	my $format;


	if ($gff3_type eq 'fly') {
		$format = 'evaluator::gff3_to_phatHit::FlyBase';
	}
	elsif ($gff3_type eq 'worm') {
		$format = 'evaluator::gff3_to_phatHit::WormBase';
	}
	else {
		die "Not a valid gff3 file type!\n";
	}	

	my $object = "$format"->new($gff3_file, $fasta_file, $ids);
	return $object;
}
#------------------------------------------------------------------------



