#------------------------------------------------------------------------
#----                      repeatmasker::PhatHit                     ---- 
#------------------------------------------------------------------------
package repeatmasker::PhatHit;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Bio::Search::Hit::PhatHit::Base;
@ISA = qw(
	Bio::Search::Hit::PhatHit::Base
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
#------------------------------------------------------------------------
sub nB {
        my $self = shift;
        my $w    = shift;

        my $sorted;
        my $hsp;

        if    ($self->strand($w) == -1){
                $sorted = $self->revSortFeatures($w);
                $hsp = shift(@{$sorted});
                return  $hsp->end($w);
        }
        elsif ($self->strand($w) == 1){
                $sorted = $self->sortFeatures($w);
                $hsp = shift(@{$sorted});
                return  $hsp->start($w);
        }
        else {
                die "dead in repeatmasker::PhatHit::nB\n";
        }
}
#------------------------------------------------------------------------
sub nE {
        my $self = shift;
        my $w    = shift;
        my $sorted;
        my $hsp;
        if    ($self->strand($w) == -1){
                $sorted = $self->revSortFeatures($w);
                $hsp = pop(@{$sorted});
                return  $hsp->start($w);
        }
        elsif ($self->strand($w) == 1){
                $sorted = $self->sortFeatures($w);
                $hsp = pop(@{$sorted});
                return  $hsp->end($w);
        }
        else {
                die "dead in repeatmasker::PhatHit::nE\n";
        }

}
#------------------------------------------------------------------------
sub show {
	my $self = shift;

	print "---------------------------------------------------------\n";
	print "------------- repeatmasker::PhatHit ---------------------\n";
	print "description:".$self->hsp(0)->hit->seq_id()."\n";
	print "---------------------------------------------------------\n";
	print $self->name()."\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
	print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "pAq:". $self->pAq."\n";
	print "pAH:". $self->pAh."\n";
	print "E/P:". $self->significance()."\n";
	print "queryLength:".$self->queryLength()."\n";
	print "sbjctLength:".$self->length()."\n";
	my $i = 0;
	foreach my $hsp ($self->hsps){
		print "-------------- HSP:$i -----------------\n";
		$hsp->show();
		$i++;
	}
	print "---------------------------------------------------------\n";
}
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "repeatmasker::PhatHit::AutoLoader called for: ",
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


