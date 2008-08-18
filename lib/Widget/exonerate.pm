#------------------------------------------------------------------------
#----                        Widget::exonerate                       ---- 
#------------------------------------------------------------------------
package Widget::exonerate;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use IPC::Open3;

@ISA = qw(
	Widget
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class  = shift;
        my @args   = @_;

        my $self = $class->SUPER::new(@args);

	bless ($self, $class);
        return $self;
}
#------------------------------------------------------------------------------
sub run {
	my $self    = shift;
	my $command = shift;

	if (defined($command)){
		$self->print_command($command);
		my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, $command);
		local $/ = \1;
		while (my $line = <CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
		}
		waitpid $pid, 0;
	}
	else {
		die " Widget::exonerate::run needs a command!\n";

	}

}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::exonerate::AutoLoader called for: ",
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


