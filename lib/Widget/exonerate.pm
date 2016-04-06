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
sub get_model_order {
    my $v = shift;

    my $str = '';
    foreach my $o (@{$v->{operations}}){
	$str .= $o->{state};
    }

    my $type;
    if ($str =~ /3I5/ && $str =~ /5I3/){
	die "ERROR: Mixed model in Widget::exonerate!\n";
    }
    elsif ($str =~ /5I3/){
	$type = '5I3';
    }
    elsif ($str =~ /3I5/){
	$type = '3I5';
    }

    if(!$type){
            if($str =~ /[A-HJ-Z]53|53[A-HJ-Z]/ && $str =~ /[A-HJ-Z]35|35[A-HJ-Z]/){
                die "ERROR: Mixed zero intron model in Widget::exonerate!\n";
            }
            elsif($str =~ /[A-HJ-Z]53|53[A-HJ-Z]/){
                $type = '5I3';
            }
            elsif($str =~ /[A-HJ-Z]35|35[A-HJ-Z]/){
                $type = '3I5';
            }
            elsif($str =~ /[35]/){
                die "ERROR: Could not determine model type in Widget::exonerate!\n";
            }
        }

        return $type;
}
#------------------------------------------------------------------------
1;


