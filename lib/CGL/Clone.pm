###################################################### main header begin ##

=head1 NAME

Clone.pm - Make a "deep" clone of an object.

=head1 DESCRIPTION

The constructors for many of the CGL classes create a new object based
on another object.  The package encapsulates the clone() routine that
they all use.

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

  CGL::Annotation::Feature::FeatureLocation.pm
  CGL::Annotation::Feature::Intron.pm
  CGL::Annotation::Feature::Exon.pm
  CGL::Annotation::Feature::Transcript.pm
  CGL::Annotation::Feature::Gene.pm
  CGL::Annotation::Feature::Protein.pm
  CGL::Annotation::Feature::Contig.pm

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::Clone;

use strict;
use warnings;

use Data::Dumper;

# subroutine exporting version
use vars qw( $VERSION @EXPORT_OK %EXPORT_TAGS );
BEGIN {
    $VERSION     = 0.01;
    @EXPORT_OK   = qw( clone );
    # %EXPORT_TAGS = ( ... );
    # don't pull in Way Too Much from Exporter;
    require Exporter;
    *import = \&Exporter::import;
}

################################################ subroutine header begin ##

=head2 clone

 Usage     : This is a toy example of how to use this function.  See
           : one of the CGL classes that uses CGL::Clone for a more
           : concrete example.

=for example begin

  use CGL::Clone qw(clone);
  use CGL::Annotation::Feature::Exon;

  my $object;
  my $newObject;

  # In real life, a real CGL::Annotation::Feature::Exon...
  $object = bless {}, CGL::Annotation::Feature::Exon;
  $newObject = clone($object);

=for example end

=for example_testing
  isa_ok($newObject, "CGL::Annotation::Feature::Exon",
         "Check if it's the right type.");

 Purpose   : Clone an object, using Data::Dumper;
 Returns   : A reference to the newly created object.
 Argument  : A reference to the object to clone.
 Throws    :
 Comments  :
 See Also  :

=cut

################################################## subroutine header end ##

sub clone
 {
	my $feature = shift;

	my $d = Data::Dumper->new([$feature], ['new_obj']);
	   $d->Purity(1);

	my $s = $d->Dump;

	my $new_obj;
	eval "$s";

	return $new_obj;


}

1; #this line is important and will help the module return a true value
__END__
