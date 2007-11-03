#----------------------------------------------------------------------------
#----                             Widget                                 ---- 
#----------------------------------------------------------------------------
package Widget;
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
	my @args   = @_;

        my $self = {};
        bless $self;

	return $self;
}
#------------------------------------------------------------------------------
sub run {
	my $self  = shift;
	my $params = shift;
	die "run is an abstract method\n";
}
#-------------------------------------------------------------------------------
sub print_command {
	my $self    = shift;
	my $command = shift;

	print STDERR "#--------- command -------------#\n";
	print STDERR ref($self).":\n";
	if (defined($command)){
		print STDERR $command."\n";
	}
	else {
		print STDERR "executing default command...\n";
	}
	print STDERR "#-------------------------------#\n";
}
#-------------------------------------------------------------------------------
sub queryName {
	my $self = shift;
	my $name = shift;

	if    (defined($name)){
		$self->{queryName} = $name;
	}
	elsif (defined($self->{queryName})){
		return $self->{queryName};
	}
	else {
		my $file = $self->queryFastaFile();
		my ($name) = $file =~ /.*\/(\S+)$/;
		$self->{queryName} = $name;
		return $self->{queryName};
	}
}
#-------------------------------------------------------------------------------
sub queryFasta {
	my $self = shift;

	my $fh = new FileHandle();
	   $fh->open($self->queryFastaFile);

	my $fasta = '';
	while (my $line .= <$fh>){ $fasta .= $line};
	$fh->close;
	return $fasta;	
}
#-------------------------------------------------------------------------------
sub parse {
	my $self   = shift;
	my $params = shift;

	die "parse is an abstract method\n";
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub toFasta {
        my $def = shift;
        my $seq = shift;

        my $fasta = $def."\n";
        for (my $i=0; $i< length($seq);$i+=60){
                $fasta .= substr($seq, $i, 60). "\n";
        }
        return $fasta;
}

#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::AutoLoader called for: ",
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


