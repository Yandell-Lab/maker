#------------------------------------------------------------------------
#----                        Widget::iprscan                         ---- 
#------------------------------------------------------------------------
package Widget::iprscan;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Widget;
use IPC::Open3;
use POSIX ":sys_wait_h";
use File::Path;
use Cwd;

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
	   my $fail;

	   #run with alarm to correct for program hanging
	   eval{
	       local $SIG{ALRM} = sub { die "ERROR: The iprscan instance is frozen\n" };
	       alarm 600;
	       while (my $line = <CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
		   $err .= $line;
	       }
	       alarm 0;
	   };

	   if($@){
	       #warn $@; #removed - don't report first round errors
	       $fail = 1
	   }

	   #reap process
	   my $stat = waitpid($pid, WNOHANG);
	   $fail = 1 if($?);

	   #reap failed so kill
	   for(my $i = 0; $i < 10 && ! $stat; $i++){
	       sleep 1;
	       kill(9, $pid);
	       $stat = waitpid($pid, WNOHANG);
	       $fail = 1 if($?);
	   }
	   $fail = 1 if(! $stat);

	   #close handles
	   close(CHLD_IN);
	   close(CHLD_OUT);
	   close(CHLD_ERR);

	   #check for correct STDERR and cleanup tmpdir
	   if($err =~ /^SUBMITTED iprscan-(\d+)-(\d+)\n*$/){
	       my $dir = "$1/iprscan-$1-$2";
	       my ($exe) = $command =~ /^(.*iprscan) -cli .* -appl /;
	       $exe = Cwd::abs_path($exe);
	       my ($base) = $exe =~ /^(.*\/)bin\/iprscan/;
	       my $tmpdir = File::Spec->tmpdir();
	       if(-d "$base/tmp/$dir"){
		   eval{File::Path::rmtree("$base/tmp/$dir");}; #ignores failure
	       }
	       elsif(-d "$tmpdir/$dir"){
		   eval{File::Path::rmtree("$tmpdir/$dir");}; #ignores failure
	       }
	   }
	   else{
	       $fail = 1;
	   }

	   #stop if everything is ok
	   return if (! $fail);

	   #always try twice because iprscan is unstable
	   $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, $command);
	   undef $err;
	   undef $fail;

	   #run with alarm to correct for program hanging
	   eval{
	       local $SIG{ALRM} = sub { die "ERROR: The iprscan instance is frozen\n" };
	       alarm 600;
	       while (my $line = <CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
		   $err .= $line;
	       }
	       alarm 0;
	   };

	   if($@){
	       warn $@;
	       $fail = 1
	   }

	   #reap process
	   $stat = waitpid($pid, WNOHANG);
	   $fail = 1 if($?);

	   #reap failed so kill
	   for(my $i = 0; $i < 10 && ! $stat; $i++){
	       sleep 1;
	       kill(9, $pid);
	       $stat = waitpid($pid, WNOHANG);
	       $fail = 1 if($?);
	   }
	   $fail = 1 if(! $stat);

	   #close handles
	   close(CHLD_IN);
	   close(CHLD_OUT);
	   close(CHLD_ERR);

	   #check for correct STDERR and cleanup tmpdir
	   if($err =~ /^SUBMITTED iprscan-(\d+)-(\d+)\n*$/){
	       my $dir = "$1/iprscan-$1-$2";
	       my ($exe) = $command =~ /^(.*iprscan) -cli .* -appl /;
	       $exe = Cwd::abs_path($exe);
	       my ($base) = $exe =~ /^(.*\/)bin\/iprscan/;
	       my $tmpdir = File::Spec->tmpdir();
	       if(-d "$base/tmp/$dir"){
		   eval{File::Path::rmtree("$base/tmp/$dir");}; #ignores failure
	       }
	       elsif(-d "$tmpdir/$dir"){
		   eval{File::Path::rmtree("$tmpdir/$dir");}; #ignores failure
	       }
	   }
	   else{
	       $fail = 1;
	   }

	   die "ERROR: Iprscan failed\n" if ($fail);
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
