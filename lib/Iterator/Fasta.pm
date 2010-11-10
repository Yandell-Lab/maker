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
use Scalar::Util qw(openhandle);

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

	$self->_set_number_of_entries($arg);

	$self->fileHandle($arg);

	return $self;
}
#-------------------------------------------------------------------------------
sub _set_number_of_entries {
	my $self = shift;
	my $arg  = shift;

	my $fh = new FileHandle();
	   $fh->open($arg);

	my $i = 0;
	my $line;
	while ($line = <$fh>){
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
		my ($id) = Fasta::getSeqID(\$query);
		$hash{$id} = $query;
	}

	return $hash{$id};

}
}
#-------------------------------------------------------------------------------
{
my @BUF; #buffer for pushing back meta character contamination
sub nextEntry {
	my $self = shift;

	#return buffered entry first
	if(@BUF){
	    return shift @BUF;
	}

	my $fh = $self->fileHandle();
	local $/ = "\n>";

	if (! openhandle($fh)){ #checks to see if file handle is open
	    return undef; 
	}
	my $line;
	while($line = <$fh>){
                $line =~ s/>//;
		$line =~ s/>$//;
                $line = ">".$line;
		$self->offsetInFile(1);
		$/ = "\n";

		if($line =~ /^M\n?|\cM\n?/){
		    $line =~ s/^M\n?|\cM\n?/\n/g;
		    my @set = grep {$_ ne "\n" } split(/\n>/, $line);
		    foreach my $s (@set){
			$s = ">".$s if($s !~ /^>/);
		    }
		    $line = shift @set;
		    push(@BUF, @set);
		}

		return $line;
	}
	
	$fh->close();
	return undef;
}
}
#-------------------------------------------------------------------------------
sub nextFasta {#alias to nextEntry
   my $self = shift;
   return $self->nextEntry;
}
#-------------------------------------------------------------------------------
sub finished {
    my $self = shift;

    my $fh = $self->fileHandle();

    if (openhandle($fh)){ #checks to see if file handle is open                                                              
	return 0;
    }
    else{
	return 1;
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

        #print STDERR "Iterator::Fasta::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


