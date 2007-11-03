#----------------------------------------------------------------------------
#----                             Iterator::Fasta                        ---- 
#----------------------------------------------------------------------------
package Iterator::Fasta;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Iterator;
use Fasta;

@ISA = qw(
		Iterator
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;
	my $arg   = shift;

        my $self = {};
        bless $self;

	$self->set_number_of_entries($arg);

	$self->fileHandle($arg);

	return $self;
}
#-------------------------------------------------------------------------------
sub set_number_of_entries {
	my $self = shift;
	my $arg  = shift;

	my $fh = new FileHandle();
	   $fh->open($arg);

	my $i = 0;
	while (my $line = <$fh>){
		$i++ if $line =~ /^>/;
	}
	$fh->close();

	$self->number_of_entries($i);
}
#-------------------------------------------------------------------------------
{
my %hash;
sub find {
	my $self = shift;
	my $id   = shift;

	return $hash{$id} if defined($hash{$id});

	while (my $query = $self->nextEntry()){
	        my $query_def = Fasta::getDef($query);

		my ($id) = $query_def =~ />(\S+)/;

		$hash{$id} = $query;
	}

	return $hash{$id};

}
}
#-------------------------------------------------------------------------------
sub nextEntry {
	my $self = shift;

	my $fh = $self->fileHandle();
	$/ = "\n>";
	while(my $line = <$fh>){
                $line =~ s/>//;
		$line =~ s/>$//;
                $line = ">".$line;
		$self->offsetInFile(1);
		$/ = "\n";
		return $line;
	}
	$/ = "\n";
	
	$fh->close();
	return 0;
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

        print STDERR "Iterator::Fasta::AutoLoader called for: ",
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


