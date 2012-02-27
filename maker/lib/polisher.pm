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
sub add_offset {
	my $offset = shift;
	my $f      = shift;

       foreach my $hsp ($f->hsps){
       		my $new_start = $offset + $hsp->start('query');
                my $new_end   = $offset + $hsp->end('query');
        
                $hsp->query->location->start($new_start);
                $hsp->query->location->end($new_end);
                $hsp->{'_sequenceschanged'} = 1;
      }
      $f->{'_sequenceschanged'} = 1;

}
#------------------------------------------------------------------------
1;


