###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Exon

CGL::Annotation::Feature::Exon - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use CGL::Annotation::Feature::Exon;

  my $feature;			# the raw CGL Feature
  my $e;			# the exon object.

  $feature = {};		# In real life, a CGL::Annotation::Feature...
  $e = new CGL::Annotation::Feature::Exon($feature);

=for example end

=for example_testing
  isa_ok($e, "CGL::Annotation::Feature::Exon", "Check its type.");

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

package CGL::Annotation::Feature::Exon;

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

=head1 new

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation::Feature::Exon;

  my $feature;			# the raw CGL Feature
  my $e;			# the exon object.

  $feature = {};		# In real life, a CGL::Annotation::Feature...
  $e = new CGL::Annotation::Feature::Exon($feature);

=for example end

=for example_testing
  isa_ok($e, "CGL::Annotation::Feature::Exon", "Check its type.");

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
	my $feature  = shift;

	my $self = clone($feature);

	bless($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head1 length

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $l;			# it's length.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $e = $t->exon(0);
  $l = $e->length();

=for example end

=for example_testing
  is($l, 1788, "Check an exon's length.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub length {
	my $self = shift;

	return abs($self->nbeg - $self->nend);
}

################################################ subroutine header begin ##

=head1 isCoding

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $p;			# the transcript's translation.
  my $e;			# reference to an exon feature
  my @coding;			# true if it's coding....

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $p = $t->translation(0);
  $e = $t->exon(0);
  $coding[0] = $e->isCoding($p, $t);
  $e = $t->exon(1);
  $coding[1] = $e->isCoding($p, $t);

=for example end

=for example_testing
  is($coding[0], 1, "Test if the first exon is coding.");
  is($coding[1], 1, "Test if the second exon is coding.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub isCoding {
	my $self  = shift;
	my $p     = shift;
	my $t     = shift;

	my $transEinTranscript = $t->translationEndInTranscript($p->id);
	my $transSinTranscript = $t->translationStartInTranscript($p->id);
	my $eSinTranscript     = $t->exonStartInTranscript($self);
	my $eEinTranscript     = $t->exonEndInTranscript($self);

	return 0 unless defined($eSinTranscript);
	return 0 unless defined($eEinTranscript);

	if ($eEinTranscript <= $transSinTranscript
	||  $eSinTranscript > $transEinTranscript){
		return 0;
	}
	else {
		return 1;
	}
}

################################################ subroutine header begin ##

=head1 metaPos

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# a transcript
  my $e;			# an exon
  my $c;			# a contig

  my $pos_in_t;
  my $pos_in_c;

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $e = $t->exon(0);
  $pos_in_t = $e->metaPos($t, 0);

  $a = new CGL::Annotation("sample_data/cint.sample.chaos.xml");
  $c = $a->contig(0);
  $t = $a->transcript(0);
  $e = $t->exon(0);
  $pos_in_c = $e->metaPos($c, 0);

=for example end

=for example_testing
  is($pos_in_t, 0, "Check an exon's position in its transcript.");
  is($pos_in_c, 9550, "Check an exon's position in its contig.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub metaPos {
  my $self          = shift;
  my $what_f        = shift;
  my $where_in_self = shift;

  my $where_in_feature;
  my $so = $self->_so();

  if ($so->a_is_hyponym_of_b($what_f->type, 'transcript') ||
      $so->a_is_hyponym_of_b($what_f->type, 'pseudogene')){

    my $self_start_in_transcript = $what_f->exonStartInTranscript($self);

    $where_in_feature = $self_start_in_transcript + $where_in_self;
  }
  elsif ($so->a_is_hyponym_of_b($what_f->type, 'contig')){
    if ($self->strand() == 1){
      return $self->nbeg + $where_in_self;
    }
    else {
      return $self->nbeg - $where_in_self;
    }
  }
  else {
    die "this metaPos not yet supported!\n";
  }

  return $where_in_feature;
}

################################################ subroutine header begin ##

=head1 order

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t_dmel;			# a transcript in the dmel annot.
  my $t_atha;			# a transcript in the atha annot.
  my $e;			# reference to an exon feature
  my $i_dmel;			# position of dmel exon in dmel transcript.
  my $i_atha;			# position of atha exon in dmel transcript.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t_dmel = $a->transcript(0);
  $e = $t_dmel->exon(1);
  $i_dmel = $e->order($t_dmel);

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t_atha = $a->transcript(0);
  $e = $t_atha->exon(0);
  $i_atha = $e->order($t_dmel);	# shouldn't find it here!

=for example end

=for example_testing
  is($i_dmel, 1, "Check that it finds the correct position in dmel.");
  is($i_atha, undef, "Check that it's not found in atha.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub order {
	my $self = shift;
	my $t    = shift;

	my $i = 0;
	while (my $e = $t->exon($i)){
		return $i if $self->id() eq $e->id();
		$i++;
	}
	return undef;
}

################################################ subroutine header begin ##

=head1 AUTOLOAD

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
            print STDERR "CGL::Annotation::Feature::Exon::AutoLoader called for: ",
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

