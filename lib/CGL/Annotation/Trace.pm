###################################################### main header begin ##

=head1 CGL::Annotation::Trace

CGL::Annotation::Trace - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use CGL::Annotation::Trace;
  my $foo = new CGL::Annotation::Trace;

=for example end

=for example_testing
  isa_ok($foo, "CGL::Annotation::Trace", "Check if it's the right type.");

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

package CGL::Annotation::Trace;

use strict;
use warnings;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub new {
	my $class    = shift;

	my $self = {};
	
	$self->{features} = [];

	bless($self, $class);
	return $self;
}

################################################ subroutine header begin ##

=head2 depth

 Usage     :

=for example begin

  # this is an un-natural example, features aren't really
  # strings....
  use CGL::Annotation::Trace;
  my $t;			# a trace object
  my $count_zero;		# a counter
  my $count_two;		# another counter

  $t = new CGL::Annotation::Trace;

  $count_zero = $t->depth();

  $t->add_feature("Dummy Feature", 0); # strings don't really make sense...
  $t->add_feature("A dummy feature at level 1", 1);
  $t->add_feature("Another dummy feature at level 1", 1);

  $count_two = $t->depth();

=for example end

=for example_testing
  is($count_zero, 0, "Did it start off empty.");
  is($count_two, 2, "Does it count correctly.");

 Purpose   : Counts the number of levels of data added to the trace.
 Returns   : The number of levels, an integer >= 0.
 Argument  : None.
 Throws    : None.
 Comments  :
 See Also  :

=cut

################################################## subroutine header end ##

sub depth {
  my $self = shift;

  return scalar(@{$self->{features}}); # return the number of levels.
}

################################################ subroutine header begin ##

=head2 terminii

 Usage     :

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  : Appears to be unused.
 See Also  :

=cut

################################################## subroutine header end ##

sub terminii {
	my $self = shift;
	
	return $self->features($self->depth -1);
}

################################################ subroutine header begin ##

=head2 features

 Usage     :

=for example begin

  # this is an un-natural example, features aren't really
  # strings....
  use CGL::Annotation::Trace;
  my $t;			# a trace object

  $t = new CGL::Annotation::Trace;

  $t->add_feature("Dummy Feature", 0); # strings don't really make sense...
  $t->add_feature("A dummy feature at level 1", 1);
  $t->add_feature("Another dummy feature at level 1", 1);

=for example end

=for example_testing
  like(($t->features(0))[0], qr/Dummy.*/, "Check lvl 0 structure");
  like(($t->features(1))[0], qr/A dummy.*/, "Check lvl 0, index 0 structure");
  like(($t->features(1))[1], qr/Another dummy.*/, "Check lvl 1, index 1");

 Purpose   : Returns the array of features at a specified level,
           : or an array of all of the levels.
 Returns   : see above.
 Argument  : An optional level.
 Throws    :
 Comments  : the example code above is contrived, in reality you'd be adding
           : features to the trace, not strings, but features are too hard to
           : just cons up.
 See Also  : CGL::Annotation::Feature.

=cut

################################################## subroutine header end ##

sub features {
	my $self  = shift;
	my $level = shift;

	if (defined($level)){
		return @{$self->{features}->[$level]|| []};
	}
	else {
		return @{$self->{features}|| []};
	}
}

################################################ subroutine header begin ##

=head2 add_feature

 Usage     :

=for example begin

  # this is an un-natural example, features aren't really
  # strings....
  use CGL::Annotation::Trace;
  my $t;			# a trace object

  $t = new CGL::Annotation::Trace;

  $t->add_feature("Dummy Feature", 0); # strings don't really make sense...
  $t->add_feature("A dummy feature at level 1", 1);
  $t->add_feature("Another dummy feature at level 1", 1);

=for example end

=for example_testing
  like(($t->features(0))[0], qr/Dummy.*/, "Check lvl 0 structure");
  like(($t->features(1))[0], qr/A dummy.*/, "Check lvl 0, index 0 structure");
  like(($t->features(1))[1], qr/Another dummy.*/, "Check lvl 1, index 1");

 Purpose   : Add a feature to the list at a specified level.
 Returns   : Nothing.
 Argument  : A feature and a level at which to add it.
 Throws    : Nothing.
 Comments  : the example code above is contrived, in reality you'd be adding
           : features to the trace, not strings, but features are too hard to
           : just cons up.
 See Also  : CGL::Annotation::Feature.

=cut

################################################## subroutine header end ##

sub add_feature {
        my $self = shift;
        my $f    = shift;
        my $l    = shift;

        push(@{$self->{features}->[$l]}, $f);
}

################################################ subroutine header begin ##

=head2 _show

 Usage     : *private*

 Purpose   : prints a debugging dump of a trace object.
 Returns   :
 Argument  :
 Throws    :
 Comments  : unused.
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _show {
        my $self = shift;

        my $depth = $self->depth();
        print "CGL::Annotation::Trace depth:$depth\n";
        for (my $d = 0; $d < $depth; $d++){
                print "CGL::Annotation::Trace\n";
                print "   features level:$d\n";
                foreach my $f ($self->features($d)){
                        print "   "." "x $d.$f->id()."\n";
                }
        }
}

################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
 See Also  :

=cut

################################################## subroutine header end ##

sub AUTOLOAD {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::Trace::AutotoLoader called for: ",
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

