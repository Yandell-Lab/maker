###################################################### main header begin ##

=head1 Bio::Search::Hit::PhatHit::rapsearch

Bio::Search::Hit::PhatHit::rapsearch - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('rapsearch',
              'sample_data/rapsearch.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  isa_ok($hit, "Bio::Search::Hit::PhatHit::rapsearch", "check type.");

=head1 DESCRIPTION

Stub documentation for this module was created by
ExtUtils::ModuleMaker.  And, then it was poked, prodded, and otherwise
massaged into it's current form by George.

Hopefully the module author wasn't negligent enough to leave the stub
unedited.

Blah blah blah.

=head1 USAGE

Expand on the examples from the SYNOPSIS.

=head1 BUGS

Not yet.

=head1 AUTHOR

 Mark Yandell
 myandell@fruitfly.org
 http://www.yandell-lab.org

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

List other relevant resources.

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package Bio::Search::Hit::PhatHit::rapsearch;

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
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('rapsearch',
              'sample_data/rapsearch.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  isa_ok($hit, "Bio::Search::Hit::PhatHit::rapsearch", "check type.");

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
	#print "queryLength:".$self->queryLength()."\n";
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

        if($ENV{CGL_CHATTER}) {
	    print STDERR "rapsearch::PhatHit::AutoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

