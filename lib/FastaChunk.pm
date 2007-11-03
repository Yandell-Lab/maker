#----------------------------------------------------------------------------
#----                            FastaChunk                              ---- 
#----------------------------------------------------------------------------
package FastaChunk;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Fasta;
use FastaFile;

@ISA = qw(
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;

        my $self = {};
        bless $self;

	return $self;
}
#-------------------------------------------------------------------------------
sub write_file {
	my $self      = shift;
	my $file_name = shift;


	$self->fasta_file_location($file_name);

	FastaFile::writeFile($self->fasta, $file_name);
}
#-------------------------------------------------------------------------------
sub erase_fasta_file {
	my $self     = shift;
	my $location = shift;

	if    (defined($location)){
		unlink($location);
	}
	elsif (defined($self->fasta_file_location)){
		unlink($self->fasta_file_location);
	}
	else {
		print STDERR "cant find a file to erase!\n";
	}
}
#-------------------------------------------------------------------------------
sub fasta {
	my $self = shift;

	my $def = $self->def();
	my $seq = $self->seq();

	return Fasta::toFasta($def, \$seq);
}
#-------------------------------------------------------------------------------
#---------------------------  CLASS FUNCTIONS ----------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "FastaChunk::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------
1;


