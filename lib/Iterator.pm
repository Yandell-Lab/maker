#----------------------------------------------------------------------------
#----                             Iterator                               ---- 
#----------------------------------------------------------------------------
package Iterator;
use strict;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use Exporter;
use PostData;
use FileHandle;
use Carp;

@ISA = qw();

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;
	my $arg   = shift;

        my $self = {};
        bless $self;

	$self->fileHandle($arg);

	return $self;
}
#-------------------------------------------------------------------------------
sub nextEntry {
        my $self = shift;

        my $fh = $self->fileHandle();
	my $line;
        while($line = <$fh>){
	    return $line;
        }
        $fh->close();
        return undef;
}
#-------------------------------------------------------------------------------
sub fileHandle {
	my $self = shift;
	my $arg  = shift;

	if    (defined($arg) && ref($arg) eq 'FileHandle'){
	    die "ERROR: You must provide a file name and not a file handle Iterator::fileHandle\n";
	}
	elsif (defined($arg) && -e $arg){
		my $fh = new FileHandle();
		$fh->open("$arg") or die "ERROR: Could not open file: $!\n";

		$self->{fileHandle} = $fh;
		$self->startPos($self->fileHandle()->getpos());
	}
	elsif (!defined($arg) && defined($self->{fileHandle})){
		return $self->{fileHandle};
	}
	else {
		confess "unknown event in Iterator::fileHandle\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Iterator::AutoLoader called for: ",
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


