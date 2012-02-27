###################################################### main header begin ##

=head1 NAME

PROTO - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

=head1 DESCRIPTION

Stub documentation for this module was created by
ExtUtils::ModuleMaker.  And, then it was poked, prodded, and otherwise
massaged into it's current form by George.

Hopefully the module author wasn't negligent enough to leave the stub
unedited.

Blah blah blah.

=head1 USAGE

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



=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package PROTO;
use strict;
use warnings;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}

# subroutine exporting version
#use vars qw( $VERSION @EXPORT_OK %EXPORT_TAGS );
#BEGIN {
#    $VERSION     = 0.01;
#    @EXPORT_OK   = qw( ... );
#    %EXPORT_TAGS = ( ... );
#    # don't pull in Way Too Much from Exporter;
#    require Exporter;
#    *import = \&Exporter::import;
#}


################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end
=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : 

=cut

################################################## subroutine header end ##


sub new
{
	my ($class, %parameters) = @_;

	my $self = bless ({}, ref ($class) || $class);

	return ($self);
}

################################################ subroutine header begin ##

=head2 attr

 Usage     :

=for example begin

  use PROTO;
  my $foo;
  my $attribute;

  $foo = new PROTO;
  $foo->attr("testing");
  $attribute = $foo->attr();

=for example end
=for example_testing


 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : 

=cut

################################################## subroutine header end ##

#
# can call this three ways.
#  $m->foo(42)    to assign the value 42.
#  $m->foo(undef) to assign the value undef which frees whatever
#                   was previously there.
#  $m->foo()      to just retrieve the current value.
#
sub foo {
  my $self = shift;

  if (@_) {
    $self->{_foo} = shift;
  }
  return $self->{_foo};
}

################################################ subroutine header begin ##
# ....
# ....
# ....
################################################## subroutine header end ##

1; #this line is important and will help the module return a true value
__END__
