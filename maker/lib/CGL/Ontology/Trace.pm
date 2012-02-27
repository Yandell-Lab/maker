###################################################### main header begin ##

=head1 CGL::Ontology::Trace

CGL::Ontology::Trace - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use CGL::Ontology::Trace;
  my $foo = new CGL::Ontology::Trace;

=for example end

=for example_testing
  isa_ok($foo, "CGL::Ontology::Trace", "Check if it's the right type.");

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

package CGL::Ontology::Trace;

use strict;
use warnings;

use CGL::Ontology::Node;

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

sub new
 {
	my $class    = shift;

	my $self = {};

	bless($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 depth

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

sub depth
 {
	my $self = shift;

	return $#{$self->{features}} + 1;
}

################################################ subroutine header begin ##

=head2 terminii

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

sub terminii
 {
	my $self = shift;
	
	return $self->features($self->depth -1);
}

################################################ subroutine header begin ##

=head2 show

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

sub show
 {
	my $self = shift;

	my $depth = $self->depth();
	print "CGL::Ontology::Trace depth:$depth\n";
	for (my $d = 0; $d < $depth; $d++){
		print "CGL::Ontology::Trace\n";
		print "   features level:$d\n";
		foreach my $f ($self->features($d)){
			print "   "." "x $d.$f->id()."\n";
		}
	}
}

################################################ subroutine header begin ##

=head2 features

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

sub features
 {
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

sub add_feature
 {
        my $self = shift;
        my $f    = shift;
        my $l    = shift;

        push(@{$self->{features}->[$l]}, $f);
}
#-------------------------------------------------------------------------------
#------------------------------- PRIVATE ---------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------

################################################ subroutine header begin ##

=head2 AUTOLOAD

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

sub AUTOLOAD
 {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Ontology::Trace::AutotoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}

1; #this line is important and will help the module return a true value
__END__

