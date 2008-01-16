package exonerate::PhatHit::protein2genome;

use strict;
use warnings;

use Bio::Search::Hit::PhatHit::Base;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
	    Bio::Search::Hit::PhatHit::Base
	   );
}

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  isa_ok($hit, "Bio::Search::Hit::PhatHit::blastp", "check type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub new
 {
        my $class  = shift;
	my @args   = @_;


	my $self = $class->SUPER::new(@args);

	bless ($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 nB

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_begin = $hit->nB('query');
 my $h_natural_begin = $hit->nB('hit');
=for example end

=for example_testing
  is($q_natural_begin, 23, "Check the query's natural begin.");
  is($h_natural_begin, 16596, "Check the hit's natural begin.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nB
 {
        my $self = shift;
        my $w    = shift;

        my $sorted;
        my $hsp;

        if ($self->strand($w) == -1){
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
                die "dead in blastn::PhatHit::nB\n";
        }
}

################################################ subroutine header begin ##

=head2 nE

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastp',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_end = $hit->nE('query');
 my $h_natural_end = $hit->nE('hit');
=for example end

=for example_testing
  is($q_natural_end, 742, "Check the query's natural end.");
  is($h_natural_end, 17492, "Check the hit's natural end.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nE
 {
        my $self = shift;
        my $w    = shift;
        my $sorted;
        my $hsp;

        if ($self->strand($w) == -1){
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
                die "dead in blastp::PhatHit::nE\n";
        }

}

################################################ subroutine header begin ##

=head2 show

 Usage     : How to use this function/method

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub show
 {
	my $self = shift;

	print "---------------------------------------------------------\n";
	print "-- ", ref($self), " -------------\n";
	print "---------------------------------------------------------\n";
	print $self->name()."\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
	print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	#print "pAq:". $self->pAq."\n";
	#print "pAH:". $self->pAh."\n";
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

################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub AUTOLOAD
 {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "blastn::PhatHit::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

