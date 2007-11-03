###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Contig

CGL::Annotation::Feature::Contig - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

=for example end

=for example_testing


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

package CGL::Annotation::Feature::Contig;

use strict;
use warnings;

use CGL::Annotation::Feature;
use CGL::Clone qw(clone);

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
	    CGL::Annotation::Feature
	   );
}

################################################ subroutine header begin ##

=head2 new

 Usage     :

=for example begin

  use CGL::Annotation::Feature::Contig;

  my $feature;			# the raw CGL Feature
  my $c;			# the contig object

  $feature = {};		# It's actually a CGL::Annotation::Feature
  $c = new CGL::Annotation::Feature::Contig($feature);

=for example end

=for example_testing
  isa_ok($c, "CGL::Annotation::Feature::Contig", "Check its type.");

 Purpose   : Create a CGL::Annotation::Feature::Contig object.
 Returns   : The object.
 Argument  : A CGL::Annotation::Feature object that is actually a contig.
 Throws    : Nothing.
 Comments  : The example of silly, in real life, $feature would be
           : a CGL::Annotation::Feature object.
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub new {
  my $class    = shift;
  my $feature  = shift;

  my $self = clone($feature);

  bless($self, $class);

  return $self;
}

################################################ subroutine header begin ##

=head2 location

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $c;			# reference to a contig feature object
  my $loc;
  my $loc10;

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $c = $a->contig(0);
  $loc = $c->location(0);
  $loc10 = $c->location(10);

=for example end

=for example_testing
  isa_ok($loc, "CGL::Annotation::FeatureLocation", "Check it's type.");
  is($loc->nbeg, 8174345, "Is it the location we expected?");
  is($loc10, undef, "Test for a crazy location.");

 Purpose   : Retrieve one of the contig's FeatureLocations.
 Returns   : The requested FeatureLocation object, undef if it doesn't exist.
 Argument  : The index of the location.
 Throws    : Nothing.
 Comments  :
 See Also  : CGL::Annotation::FeatureLocation

=cut

################################################## subroutine header end ##

sub location {
	my $self = shift;
	my $i    = shift;

	return $self->{locations}->[$i];
}

################################################ subroutine header begin ##

=head2 metaPos

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $c;			# reference to a contig feature object
  my $e_ref;			# ref to list of exons.
  my $e;			# an exon
  my $t;			# a transcript
  my $pos_in_e;
  my $pos_in_t;

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $c = $a->contig(0);
  $e_ref = $a->exons();
  $e = $e_ref->[0];		# snag the first one.
  $pos_in_e = $c->metaPos($e, 501);

  $a = new CGL::Annotation("sample_data/cint.sample.chaos.xml");
  $c = $a->contig(0);
  $t = $a->transcript(0);
  $pos_in_t = $c->metaPos($t, 8000);

=for example end

=for example_testing
  is($pos_in_e, 1, "Check metaPos in exon.");
  is($pos_in_t, 921, "Check metaPos in transcript.");

 Purpose   : Get the position in another object of a position in this object.
 Returns   : The position, or undef it it isn't a reasonabel request.
 Argument  : The other feature of interest, and the position in this object.
 Throws    : Nothing.
 Comments  : Currently only implemented for exons and transcripts.
           : XXXX (remember hyponym stuff...
 See Also  :

=cut

################################################## subroutine header end ##

sub metaPos {
        my $self          = shift;
        my $what_f        = shift;
        my $where_in_self = shift;

	my $so = $self->_so();

        my $where_in_feature;
	if ($so->a_is_hyponym_of_b($what_f->type, 'exon')){
		my $pos_of_e_begin_on_s = $what_f->metaPos($self, 0);
		my $pos_of_e_end_on_s   = $what_f->metaPos($self, $what_f->length);

		if ($what_f->strand() == 1){
			return undef if $where_in_self < $pos_of_e_begin_on_s;
			return undef if $where_in_self > $pos_of_e_end_on_s;

			return $where_in_self - $pos_of_e_begin_on_s;
		}
		else {
			return undef if $where_in_self > $pos_of_e_begin_on_s;
			return undef if $where_in_self < $pos_of_e_end_on_s;
		
			return $pos_of_e_begin_on_s - $where_in_self;
		}
        }
	if ($so->a_is_hyponym_of_b($what_f->type, 'transcript')){
		# caller was $s->metaPos($t, some_pos_on_s);
		my $i = 0;
		my $t_pos = 0;
		while (my $e = $what_f->exon($i)){
			my $e_pos = $self->metaPos($e, $where_in_self);

			return $e->metaPos($what_f, $e_pos) if defined($e_pos);

			$i++;
		}
		return undef;
	}
	else {
		die "this metaPos not supported yet!\n";
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
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub AUTOLOAD {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;
	if ($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::Feature::Contig AutoLoader called for: ",
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

