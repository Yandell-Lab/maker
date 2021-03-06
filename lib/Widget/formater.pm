#------------------------------------------------------------------------
#----                        Widget::formater                        ---- 
#------------------------------------------------------------------------
package Widget::formater;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use PhatHit_utils;
use IPC::Open3;
use Symbol;

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
	        my ($exe) = $command =~ /^([^\s]+)/;
		$self->print_command($command);
		my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
		my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
		local $/ = \1;
		while (my $line = <$CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
		}
		waitpid $pid, 0;
		die "ERROR: $exe failed in Widget::formater\n" if $? != 0;
		unlink("formatdb.log") if(-e "formatdb.log");
		unlink("setdb.log") if(-e "setdb.log");
	}
	else {
		die "you must give Widget::formater a command to run!\n";
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

        #print STDERR "Widget::blastn::AutoLoader called for: ",
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


