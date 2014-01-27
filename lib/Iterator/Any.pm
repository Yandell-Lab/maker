#-------------------------------------------------------------------------------
#------                               Iterator::Any                     -------
#-------------------------------------------------------------------------------
package Iterator::Any;
use strict;
use Iterator;
use Iterator::Fasta;
use Iterator::GFF3;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use Exporter;

@ISA = qw(Iterator);

#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub new {
        my $class      = shift;
	my %args = @_;
	my $iter;
	
	my $fasta_file = $args{-fasta} || undef;
	my $gff_file = $args{-gff} || undef;

	if($gff_file){
	   $iter = new Iterator::GFF3($gff_file, $fasta_file);
	}
	elsif($fasta_file){
	   $iter = new Iterator::Fasta($fasta_file);
	}
	else{
	   die "ERROR: Failure to identify input arguments in Iterator::Any\n";
	}

	return $iter;
}
#------------------------------------------------------------------------
1;
