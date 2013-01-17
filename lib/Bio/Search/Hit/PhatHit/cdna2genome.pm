#------------------------------------------------------------------------
#----                 Bio::Search::Hit::PhatHit::cdna2genome                ---- 
#------------------------------------------------------------------------
package Bio::Search::Hit::PhatHit::cdna2genome;
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
sub show {
	my $self = shift;

	print "---------------------------------------------------------\n";
	print "-------------io::Search::Hit::PhatHit::cdna2genome--------------\n";
	print "description:".$self->name()."\n";
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
	my $evi = $self->evidence || [];
	foreach my $e (@{$evi}){
		print "evidence:".substr($e, 0, 60)."\n";
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

        #print STDERR "Bio::Search::Hit::PhatHit::cdna2genome::AutoLoader called for: ",
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
