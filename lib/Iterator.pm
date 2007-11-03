#----------------------------------------------------------------------------
#----                             Iterator                               ---- 
#----------------------------------------------------------------------------
package Iterator;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
@ISA = qw(
       );

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
        while(my $line = <$fh>){
                $self->offsetInFile(1);
                return $line;
        }
        $fh->close();
        return 0;
}
#-------------------------------------------------------------------------------
sub offsetInFile {
	my $self      = shift;
	my $increment = shift;

	if (defined($increment) && defined($self->{offsetInFile})){
		$self->{offsetInFile} += $increment;
	}
	elsif (defined($increment) && !defined($self->{offsetInFile})) {
		$self->{offsetInFile} = $increment -1;
		
	}
	else {
		return $self->{offsetInFile};
	}
}
#-------------------------------------------------------------------------------
sub fileHandle {
	my $self = shift;
	my $arg  = shift;
	

	if    (defined($arg) && ref($arg) eq 'FileHandle'){
		 $self->{fileHandle} = $arg;
	}
	elsif (defined($arg) && -e $arg){
		my $fh = new FileHandle();
		   $fh->open("$arg");

		$self->{fileHandle} = $fh;
	}
	elsif (!defined($arg) && defined($self->{fileHandle})){
		return $self->{fileHandle};
	}
	else {
		die "unknown event in Iterator::fileHandle\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
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


