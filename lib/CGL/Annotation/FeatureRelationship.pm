###################################################### main header begin ##

=head1 CGL::Annotation::FeatureRelationship

CGL::Annotation::FeatureRelationship - The platonic ideal of a CGL module.

=head1 SYNOPSIS

  use CGL::Annotation::FeatureRelationship;
  # fake up some test data, real hash would come from chaos xml....
  my $fr_info = {'feature_relationship' =>
                   {'object_id' => 'gene-290267',
                    'subject_id' => 'STS-290319',
                    'type' => 'part_of'}
                };
  my $relationship = new CGL::Annotation::FeatureRelationship($fr_info);

  my $object = $relationship->oF();
  my $subject = $relationship->sF();
  my $logus = $relationship->logus();

=head1 DESCRIPTION

This class manages relationships between chaos features.

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

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::Annotation::FeatureRelationship;

use strict;
use warnings;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}

################################################ subroutine header begin ##

=head2 new

 Usage     :

=for example begin

  use CGL::Annotation::FeatureRelationship;
  # fake up some test data, real hash would come from chaos xml....
  my $fr_info = {'feature_relationship' =>
                   {'object_id' => 'gene-290267',
                    'subject_id' => 'STS-290319',
                    'type' => 'part_of'}
                };

  my $relationship = new CGL::Annotation::FeatureRelationship($fr_info);

=for example end

=for example_testing
  isa_ok($relationship, "CGL::Annotation::FeatureRelationship",
         "Check if it's the right type.");

 Purpose   : Create a new FeatureRelationship object and populate it.
 Returns   : A reference to a CGL::Annotatation::FeatureRelationship object.
 Argument  : The class and a hash of feature relationship information.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub new
 {
	my $class = shift;
	my $hash  = shift;

	my $self ={};
	bless($self, $class);

	foreach my $k (keys %{$hash->{feature_relationship}}){
		my $v = $hash->{feature_relationship}->{$k};
		
		if    ($k eq 'object_id'){
			$self->oF($v);
		}
		 elsif ($k eq 'subject_id'){
			$self->sF($v);
		}	
                elsif ($k eq 'type'){
                        $self->logus($v);
                }
		else {
			$self->{$k} = $v;
		}
	}

	return $self;
}

################################################ subroutine header begin ##

=head2 oF

 Usage     :

=for example begin

  use CGL::Annotation::FeatureRelationship;
  # fake up some test data, real hash would come from chaos xml....
  my $fr_info = {'feature_relationship' =>
                   {'object_id' => 'gene-290267',
                    'subject_id' => 'STS-290319',
                    'type' => 'part_of'}
                };
  my $relationship = new CGL::Annotation::FeatureRelationship($fr_info);
  my $object_feature = $relationship->oF();

=for example end

=for example_testing
  is($object_feature, 'gene-290267', "Check the object feature.");

 Purpose   : get/set the value of the object of the relationship.
 Returns   : The value of relationship's "object"
 Argument  : optionally, the value to set the object to.
 Throws    :
 Comments  : Explicitly passing 'undef' will undefine the field in the
           : underlying object.
 See Also  :

=cut

################################################## subroutine header end ##

sub oF {
  my $self  = shift;

  if (@_) {
    $self->{oF} = shift;
  }
  return $self->{oF};
}

################################################ subroutine header begin ##

=head2 sF

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation::FeatureRelationship;
  # fake up some test data, real hash would come from chaos xml....
  my $fr_info = {'feature_relationship' =>
                   {'object_id' => 'gene-290267',
                    'subject_id' => 'STS-290319',
                    'type' => 'part_of'}
                };
  my $relationship = new CGL::Annotation::FeatureRelationship($fr_info);
  my $subject_feature = $relationship->sF();

=for example end

=for example_testing
  is($subject_feature, 'STS-290319', "Check the subject of the feature.");

 Purpose   : Get/set the subject of the relationship.
 Returns   : The value of the subject of the relationship.
 Argument  : Optionally, a value to set.
 Throws    :
 Comments  : Explicitly passing 'undef' will undefine the field in the
           : underlying object.
 See Also  :

=cut

################################################## subroutine header end ##

sub sF {
  my $self  = shift;

  if (@_) {
    $self->{sF} = shift;
  }
  return $self->{sF};
}

################################################ subroutine header begin ##

=head2 logus

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation::FeatureRelationship;
  # fake up some test data, real hash would come from chaos xml....
  my $fr_info = {'feature_relationship' =>
                   {'object_id' => 'gene-290267',
                    'subject_id' => 'STS-290319',
                    'type' => 'part_of'}
                };
  my $relationship = new CGL::Annotation::FeatureRelationship($fr_info);
  my $logus = $relationship->logus();

=for example end

=for example_testing
  is($logus, 'part_of', "Check the relationship's.");

 Purpose   : get/set the logus value for the relationship.
 Returns   : the value of the field.
 Argument  : optionally, a value to set.
 Throws    :
 Comments  : Explicitly passing 'undef' will undefine the field in the
           : underlying object.
           : The term is a reference to the invisible forces that bind
           : the universe together in Zelazny's Amber novels.
 See Also  :

=cut

################################################## subroutine header end ##

sub logus {
  my $self  = shift;

  if (@_) {
    $self->{logus} = shift;
  }
  return $self->{logus};
}

################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : This routine should never be called directly.

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
	    print STDERR "CGL::Annotation::FeatureRelationship::AutoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (@_){
	  $self->{$call} = shift
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

