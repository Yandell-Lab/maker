###################################################### main header begin ##

=head1 CGL::Annotation::FeatureLocation

CGL::Annotation::FeatureLocation - A class to keep track of feature
locations.

=head1 SYNOPSIS

=head1 DESCRIPTION

Objects in this class are created from a hash that is in turn
populated by the xml parsing routine in CGL::Annotation.pm.

Other than new, all of the methods are explicit calls that have the
form of the CGL generic getter/setter.  Calling one of these methods
without any arguments will return the corresponding field's current
value.  Calling it with an argument will set the corresponding field
to the value of that argument.  Calling a method with an explicit
value of undef will undefine that field in the object's underlying
hash.

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

CGL::Annotation.pm
CGL::Annotation::Feature.pm

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::Annotation::FeatureLocation;

use strict;
use warnings;

use CGL::Clone qw(clone);

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );
}

################################################ subroutine header begin ##

=head2 new

 Usage     :

=for example begin

  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

=for example end

=for example_testing
  isa_ok($f_loc, "CGL::Annotation::FeatureLocation",
         "Check if it's the right type.");

 Purpose   : Create a new FeatureLocation object.
 Returns   : A reference to the new object.
 Argument  : A hash full of info derived from the underlying chaos-xml.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub new
 {
	my $class = shift;
	my $loc = shift;

	my $self = clone($loc);
	bless($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 group

 Usage     :

=for example begin

  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  $f_loc->group("Testing");
  my $group_value = $f_loc->group();
  $f_loc->group(undef);
  my $group_value2 = $f_loc->group();

=for example end

=for example_testing
  is($group_value, "Testing", "Check if it's the right value.");
  is($group_value2, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the object's group field.
 Returns   : The current value of the group field.
 Argument  : Optionally, a scalar value to set for the field
 Throws    :
 Comments  : Explicitly passing in 'undef' will undefine the value
           : in the hash.
 See Also  :

=cut

################################################## subroutine header end ##

sub group {
  my $self  = shift;

 if (@_) {
    $self->{group} = shift;
  }
  return $self->{group};
}

################################################ subroutine header begin ##

=head2 nbeg

 Usage     : How to use this function/method

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  my $nbeg = $f_loc->nbeg();
  $f_loc->nbeg(505);
  my $nbeg2 = $f_loc->nbeg();
  $f_loc->nbeg(undef);
  my $nbeg3 = $f_loc->nbeg();

=for example end

=for example_testing
  is($nbeg, 500, "Check if it's the right value. (500)");
  is($nbeg2, 505, "Check if it's the right value.(505)");
  is($nbeg3, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the feature's natural begin.
 Returns   : The current natural begin for the feature.
 Argument  : Optionally, a new value for the natural begin.
 Throws    :
 Comments  : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  :

=cut

################################################## subroutine header end ##

sub nbeg {
  my $self  = shift;

  if (@_) {
    $self->{nbeg} = shift;
  }
  return $self->{nbeg};

}

################################################ subroutine header begin ##

=head2 nend

 Usage     : How to use this function/method

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  my $nend = $f_loc->nend();
  $f_loc->nend(505);
  my $nend2 = $f_loc->nend();
  $f_loc->nend(undef);
  my $nend3 = $f_loc->nend();

=for example end

=for example_testing
  is($nend, 26220, "Check if it's the right value. (26220)");
  is($nend2, 505, "Check if it's the right value.(505)");
  is($nend3, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the feature's natural end.
 Returns   : The current natural end for the feature.
 Argument  : Optionally, a new value for the natural end.
 Throws    :
 Comments  : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  :

=cut

################################################## subroutine header end ##

sub nend {
  my $self  = shift;

  if (@_) {
    $self->{nend} = shift;
  }
  return $self->{nend};
}

################################################ subroutine header begin ##

=head2 srcfeature

 Usage     :

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  $f_loc->srcfeature(505);
  my $srcfeature = $f_loc->srcfeature();
  $f_loc->srcfeature(undef);
  my $srcfeature2 = $f_loc->srcfeature();

=for example end

=for example_testing
  is($srcfeature, 505, "Check if it's the right value.(505)");
  is($srcfeature2, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the srcfeature field.
 Returns   : The value of the srcfeature field.
 Argument  : An optional argument to set.
 Throws    :
 Comments  : In the real world, this would store a reference to
           : a feature object.
           : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub srcfeature {
  my $self  = shift;

  if (@_) {
    $self->{srcfeature} = shift;
  }
  return $self->{srcfeature};
}

################################################ subroutine header begin ##

=head2 srcfeature_id

 Usage     : How to use this function/method

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  my $srcfeature_id = $f_loc->srcfeature_id();
  $f_loc->srcfeature_id(undef);
  my $srcfeature_id2 = $f_loc->srcfeature();

=for example end

=for example_testing
  is($srcfeature_id, 'contig-3088', "Check if it's the right value.(contig-3088)");
  is($srcfeature_id2, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the srcfeature_id field.
 Returns   : The value of the srcfeature_id field.
 Argument  : An optional argument to set.
 Throws    :
 Comments  : In the real world, this would store a reference to
           : a feature object.
           : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub srcfeature_id {
  my $self  = shift;

  if (@_) {
    $self->{srcfeature_id} = shift;
  }
  return $self->{srcfeature_id};
}

################################################ subroutine header begin ##

=head2 strand

 Usage     :

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  $f_loc->strand("+");
  my $strand = $f_loc->strand();
  $f_loc->strand(undef);
  my $strand2 = $f_loc->strand();

=for example end

=for example_testing
  is($strand, '+', "Check if it's the right value.(+)");
  is($strand2, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the strand field.
 Returns   : The value of the strand field.
 Argument  : An optional argument to set.
 Throws    :
 Comments  : In the real world, this would store a reference to
           : a feature object.
           : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub strand {
  my $self  = shift;

  if (@_) {
    $self->{strand} = shift;
  }
  return $self->{strand};
}

################################################ subroutine header begin ##

=head2 rank

Usage:

=for example begin
  use CGL::Annotation::FeatureLocation;
  my $location = {'srcfeature_id' => 'contig-3088',
                  'nbeg' => 500,
                  'nend' => 26220,
                  'rank' => 0,
                 };
  my $f_loc = new CGL::Annotation::FeatureLocation ($location);

  my $rank = $f_loc->rank();
  $f_loc->rank(undef);
  my $rank2 = $f_loc->rank();

=for example end

=for example_testing
  is($rank, 0, "Check if it's the right value.(0)");
  is($rank2, undef, "Check if it's the right value (undef).");

 Purpose   : Get/set the rank field.
 Returns   : The value of the rank field.
 Argument  : An optional argument to set.
 Throws    :
 Comments  : In the real world, this would store a reference to
           : a feature object.
           : Explicitly passing in 'undef' will cause the field in the
           : underlying hash to be undefined.
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub rank {
  my $self  = shift;

  if (@_) {
    $self->{rank} = shift;
  }
  return $self->{rank};
}

################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : *private*

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
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

	if ($ENV{CGL_CHATTER}){
        	print STDERR "CGL::Annotation::FeatureLocation::AutoLoader called for: ",
              "\$self->$call","()\n";
        	print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (@_) {
	  $self->{$call} = shift;
        }
	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

