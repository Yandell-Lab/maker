#------------------------------------------------------------------------
#----                   repeatmasker::PhatHsp                        ---- 
#------------------------------------------------------------------------
package repeatmasker::PhatHsp;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Bio::Search::HSP::PhatHSP::Base;
@ISA = qw(
	Bio::Search::HSP::PhatHSP::Base
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
#-------------------------------------------------------------------------------
sub strand {
	my $self = shift;
	my $w    = shift;

	$w = 'hit' if $w eq 'sbjct';

	return $self->{_strand_hack}->{$w};
}
#-------------------------------------------------------------------------------
sub equivalent_pos_in_alignment_partner {

	print " method equivalent_pos_in_alignment_partner is not available for ";
	print " repeatmasker::PhatHsps\n";

	die;
}
#-------------------------------------------------------------------------------
sub whatIsThere {

       print " method whatIsThere is not available for ";
       print " repeatmasker::PhatHsps\n";

        die;

}
#-------------------------------------------------------------------------------
sub whatIsInTheMiddle {
	

	print " method whatIsInTheMiddle is not available for ";
        print " repeatmasker::PhatHsps\n";

        die;

}
#-------------------------------------------------------------------------------
sub show {
	my $self = shift;

        print "--------------------------------------\n";
        print "repeatmasker::PhatHSP\n";
	print "   ".$self->hit->seq_id()."\n";
        print "--------------------------------------\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDq:".$self->frac_identical('query')."\n";
	print "fIDt:".$self->frac_identical('total')."\n";
}
#-------------------------------------------------------------------------------
sub nB {
        my $self = shift;
        my $w    = shift;

        die unless defined($w);

        if ($self->strand($w) == 1){
                return $self->start($w);
        }
        else {
                return $self->end($w);
        }
}
#-------------------------------------------------------------------------------
sub nE {
        my $self = shift;
        my $w    = shift;

        if ($self->strand($w) == 1){
                return $self->end($w);
        }
        elsif ($w eq 'query' && $self->strand($w) == 0){
                return $self->end($w);
        }

        else {
                return $self->start($w);
        }
}
#-------------------------------------------------------------------------------
sub hasRun {

        print " method hasRun is not available for ";
        print " repeatmasker::PhatHsps\n";

        die;

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

        print STDERR "repeatmasker::PhatHsp::AutoLoader called for: ",
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


