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
	my $o_file  = shift;

	if (defined($command)){
		$self->print_command($command);
		my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
		my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
		{
		    local $/ = \1;
		    while (my $line = <$CHLD_ERR>){
			print STDERR $line unless($main::quiet);
		    }
		}
		waitpid $pid, 0;

		#try a second time
		if($? != 0){
		   $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
		   {
		       local $/ = \1;
		       while (my $line = <$CHLD_ERR>){
			   print STDERR $line unless($main::quiet);
		       }
		   }
		   waitpid $pid, 0;
		   
		   if($? != 0){
		      unlink($o_file) if($o_file && -f $o_file);
		      die "ERROR: Exonerate failed\n";
		   }
		}
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


