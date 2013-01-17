#------------------------------------------------------------------------
#----                            polisher                            ---- 
#------------------------------------------------------------------------
package polisher;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use Fasta;
use FastaFile;
use PhatHit_utils;
use Iterator::Fasta;
@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub prep {
	my $file = shift;

	my $iterator = new Iterator::Fasta($file);
	my $fasta = $iterator->nextEntry();
	my $seq   = Fasta::getSeqRef($fasta);
	my $len = length($$seq);

	return $len;
}
#------------------------------------------------------------------------
1;


