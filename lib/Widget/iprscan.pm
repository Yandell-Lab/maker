#------------------------------------------------------------------------
#----                        Widget::iprscan                         ---- 
#------------------------------------------------------------------------
package Widget::iprscan;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Widget;
use IPC::Open3;

@ISA = qw(
	Widget
       );

my $OPT_F; # global option -f to force cleanup
my $LOG; #global varible for maker runlog
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
#------------------------------------------------------------------------
sub run {
	my $self    = shift;
	my $command = shift;

	if (defined($command)){
	   $self->print_command($command);
	   my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, $command);
	   local $/ = \1;
	   my $err;
	   while (my $line = <CHLD_ERR>){
	      print STDERR $line unless($main::quiet);
	      $err .= $line;
	   }
	   waitpid $pid, 0;
	   
	   my $fail;
	   if($err !~ /^SUBMITTED iprscan-\d+-\d+\n*$/){
	       $fail = 1;
	   }

	   die "ERROR: Iprscan failed\n" if ($? != 0 || $fail);
	}
	else {
	   die "you must give Widget::iprscan a command to run!\n";
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

        print STDERR "Widget::iprscan::AutoLoader called for: ",
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
