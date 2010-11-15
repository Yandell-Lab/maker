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
	my $g_file = shift;
	my $p_file = shift;

	my $g_iterator = new Iterator::Fasta($g_file);
	my $g_fasta = $g_iterator->nextEntry();
	my $g_seq   = Fasta::getSeqRef($g_fasta);

	my $p_iterator = new Iterator::Fasta($p_file);
	my $p_fasta = $p_iterator->nextEntry();
	my $p_seq   = Fasta::getSeqRef($p_fasta);

	my $g_len = length($$g_seq);
	my $p_len = length($$p_seq);

	return ($p_len, $g_len);	
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


