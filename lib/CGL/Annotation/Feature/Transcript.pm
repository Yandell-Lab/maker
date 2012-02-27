###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Transcript

CGL::Annotation::Feature::Transcript - The platonic ideal of a CGL module.

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

package CGL::Annotation::Feature::Transcript;

use strict;
use warnings;

use CGL::Annotation::Feature;
use CGL::Annotation::Feature::Intron;
use CGL::Annotation::FeatureLocation;
use CGL::Annotation::FeatureRelationship;
use CGL::Clone qw(clone);
use CGL::Revcomp qw(revcomp);

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
	    CGL::Annotation::Feature
	   );
}

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation::Feature::Transcript;

  my $feature;			# the raw CGL Feature
  my $t;			# the transcript object.

  $feature = {};		# In real life, a CGL::Annotation::Feature...
  $t = new CGL::Annotation::Feature::Transcript($feature);

=for example end

=for example_testing
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check its type.");

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

=head2 exon_with_id

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $e = $t->exon_with_id('exon-107701');

=for example end

=for example_testing
  isa_ok($e, "CGL::Annotation::Feature::Exon", "Did we get an exon?");
  is($e->id(), 'exon-107701', "Did we get the exon we expected?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exon_with_id {
        my $self = shift;
        my $id   = shift;

        my $i = 0;
        while (my $e = $self->exon($i)){
                return $e if $e->id eq $id;
                $i++;
        }
}

################################################ subroutine header begin ##

=head2 intron_with_id

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $i;			# reference to an intron feature

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $i = $t->intron_with_id("CG4852-region-ext1.3:588:647");

=for example end

=for example_testing
  isa_ok($i, "CGL::Annotation::Feature::Intron", "Did we get back an intron?");
  is($i->id(), 'CG4852-region-ext1.3:588:647', "Is it the expected intron?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub intron_with_id {
        my $self = shift;
        my $id   = shift;

        my $i = 0;
        while (my $e = $self->intron($i)){
                return $e if $e->id eq $id;
                $i++;
        }
}

################################################ subroutine header begin ##

=head2 residues

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $r;			# the transcript's residue

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $r = $t->residues();

=for example end

=for example_testing
  like($r, qr/^GTGACC.*/, "Do the residues look correct?");
  is(length($r), 1365, "Are the residues the correct length?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub residues {
        my $self = shift;
        my $arg  = shift;

        if (defined($arg)){
                $self->{residues} = $arg;
        }
        else {
                return $self->{residues};
        }
}

################################################ subroutine header begin ##

=head2 nbeg

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $i0;			# a location
  my $i1;			# a location

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $i0 = $t->nbeg();
  $t = $a->transcript(1);
  $i1 = $t->nbeg();

=for example end

=for example_testing
  is($i0, 1, "Does the first transcript begin at the right place?");
  is($i1, 6, "Does the second transcript begin at the right place?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nbeg {
	my $self = shift;
	return $self->exons->[0]->featureLocation(0)->nbeg();
}

################################################ subroutine header begin ##

=head2 nend

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $i0;			# a location
  my $i1;			# a location

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $i0 = $t->nend();
  $t = $a->transcript(1);
  $i1 = $t->nend();

=for example end

=for example_testing
  is($i0, 1425, "Does the first transcript begin at the right place?");
  is($i1, -1990, "XXXX Does the second transcript begin at the right place?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nend {
	my $self = shift;
	return $self->exons->[-1]->featureLocation(0)->nend();
}

################################################ subroutine header begin ##

=head2 strand

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $s;			# strand

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $s = $t->strand();

=for example end

=for example_testing
  is($s, 1, "Is the strand correct?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub strand {
	my $self = shift;

	# XXXX what's this, what to do with it?
	#die if $self->_mixedStrand();

	my $nbeg = $self->{exons}->[0]->{locations}->[0]->{nbeg};
	my $nend = $self->{exons}->[0]->{locations}->[0]->{nend};
	my $strand = $nbeg < $nend ? 1 : -1;
	return $strand;
}

################################################ subroutine header begin ##

=head2 length

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $l;			# its length.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $t = $a->transcript(0);
  $l = $t->length();

=for example end

=for example_testing
  is($l, 1365, "Check the transcript's length.");

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

	my $length = 0;
	foreach my $e (@{$self->exons}){
		$length += $e->length();
	}
	return $length;
}

################################################ subroutine header begin ##

=head2 acceptor

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# refrence to a gene feature
  my $t;			# reference to a transcript feature
  my $acc;			# the acceptor

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $acc = $t->acceptor(1, 4, 4);	# grab four bases on either side.

=for example end

=for example_testing
  is($acc, 'GCAG|ATCT', "Check the first acceptor.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub acceptor {
        my $self         = shift;
        my $exonNumber   = shift;
        my $upStream     = shift;
        my $dnStream     = shift;

        if ($exonNumber == 0){
	  # XXXX should this die?
                print STDERR "No acceptor for first exon($exonNumber)\n";
                return;
        }
        my $e = $self->exon($exonNumber);
        my $i = $self->intron($exonNumber - 1);

	return undef unless defined($i);
        return undef unless (defined($e->residues) && defined($i->residues));

	my $iOffset  = $i->length() - $upStream;
        my $prefix   = substr($i->residues(), $iOffset);
        my $postfix  = substr($e->residues(), 0, $dnStream);

        return $prefix."|".$postfix;
}

################################################ subroutine header begin ##

=head2 donor

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# refrence to a gene feature
  my $t;			# reference to a transcript feature
  my $d;			# the donor

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $d = $t->donor(0, 4, 4);	# grab four bases on either side.

=for example end

=for example_testing
  is($d, 'TTCG|GTAA', "Check the first donor.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub donor {
	my $self         = shift;
	my $exonNumber   = shift;
	my $upStream     = shift;
	my $dnStream     = shift;

	if ($exonNumber == $#{$self->exons}){
		print STDERR "No donor for last exon($exonNumber)\n";
		return;
	}
	my $e = $self->exon($exonNumber);
	my $i = $self->intron($exonNumber);

	return undef unless defined($i);
	return undef unless (defined($e->residues) && defined($i->residues));

	my $prefix   = substr($e->residues(), $e->length() - $upStream, $upStream);
	my $postfix  = substr($i->residues(), 0, $dnStream);

	return $prefix."|".$postfix;
}

################################################ subroutine header begin ##

=head2 exonNumber

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $e_num;			# how many exons are there?

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $e = $t->exon(1);
  $e_num = $t->exonNumber($e);

=for example end

=for example_testing
  is($e_num, 1, "Can we retrieve an exon's position in transcript.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonNumber {
	my $self = shift;
	my $exon = shift;

	my $i = 0;
	foreach my $e (@{$self->exons}){
		return $i if $e->id eq $exon->id;
		$i++;
	}
	# XXXX yikes, what if it's not there?
}

################################################ subroutine header begin ##

=head2 isFirstExon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $yes;			# boolean
  my $no;			# bolean

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $e = $t->exon(0);
  $yes = $t->isFirstExon($e);
  $e = $t->exon(1);
  $no = $t->isFirstExon($e);

=for example end

=for example_testing
  is($yes, 1, "Did we recognize the first exon?");
  is($no, 0, "Did we recognize NOT the first exon?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub isFirstExon {
	my $self = shift;
	my $exon = shift;

	return 1 if $self->exon(0)->id eq $exon->id;
	return 0;
}

################################################ subroutine header begin ##

=head2 isLastExon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $yes;			# boolean
  my $no;			# bolean

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->exon(2);
  $yes = $t->isLastExon($e);
  $e = $t->exon(0);
  $no = $t->isLastExon($e);

=for example end

=for example_testing
  is($yes, 1, "Did we recognize the last exon?");
  is($no, 0, "Did we recognize NOT the last exon?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub isLastExon {
        my $self = shift;
        my $exon = shift;

        return 1 if $self->exon(-1)->id eq $exon->id;
        return 0;
}

################################################ subroutine header begin ##

=head2 exonStartInTranscript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $start;			# exon start in transcript

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->exon(2);
  $start = $t->exonStartInTranscript($e);

=for example end

=for example_testing
  is($start, 322, "Can we find an exon's start?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonStartInTranscript {
	my $self = shift;
	my $exon = shift;

	my $offset_on_transcript = 0;
	foreach my $e (@{$self->exons}){
		return $offset_on_transcript if $exon->id eq $e->id;
		$offset_on_transcript += $e->length();
	}
	return undef;
}

################################################ subroutine header begin ##

=head2 exonEndInTranscript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature
  my $end;			# exon start in transcript

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->exon(2);
  $end = $t->exonEndInTranscript($e);

=for example end

=for example_testing
  is($end, 569, "Can we find an exon's end?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonEndInTranscript {
        my $self = shift;
        my $exon = shift;

        my $offset_on_transcript = 0;
        foreach my $e (@{$self->exons}){
		$offset_on_transcript += $e->length();
                return $offset_on_transcript if $exon->id eq $e->id;
        }
	return undef;
}

################################################ subroutine header begin ##

=head2 fetchTranslationById

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $p;			# reference to a translation feature

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $p = $t->fetchTranslationById('385099');
  $id = $p->id();

=for example end

=for example_testing
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check its type.");
  is($id, "385099", "Check the protein's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub fetchTranslationById  {
	my $self = shift;
	my $id   = shift;

	foreach my $t (@{$self->translations}){
		return $t if $t->id eq $id;
	}
	return undef;
}

################################################ subroutine header begin ##

=head2 translation

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $p = $t->translation(0);
  $id = $p->id();

=for example end

=for example_testing
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check its type.");
  is($id, "385099", "Check the protein's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub translation {
	my $self = shift;
	my $i    = shift;

	return $self->{translations}->[$i];
}

################################################ subroutine header begin ##

=head2 firstCodingExon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $p;			# reference to a translation feature

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->firstCodingExon('385099');

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

sub firstCodingExon {
	my $self          = shift;
	my $translationId = shift;

	my $translation = $self->fetchTranslationById($translationId);
	
	foreach my $e (@{$self->exons}){
		if ($self->strand == 1){
			return $e if ($e->nbeg <= $translation->nbeg
			          &&  $e->nend >= $translation->nbeg);
		}
		else {
			return $e if ($e->nbeg >= $translation->nbeg
                                  &&  $e->nend <= $translation->nbeg);
		}
	}
	die "logic error in Chaos::Feature::Transcript::firstCodingExon\n";
	# XXXX is there a better way to handle weird case?
}

################################################ subroutine header begin ##

=head2 lastCodingExon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $p;			# reference to a translation feature

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->lastCodingExon('385099');

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

sub lastCodingExon {
        my $self          = shift;
        my $translationId = shift;

        my $translation = $self->fetchTranslationById($translationId);

        foreach my $e (@{$self->exons}){
                if ($self->strand == 1){
                        return $e if ($e->nbeg <= $translation->nend
                                  &&  $e->nend >= $translation->nend);
                }
                else {
                        return $e if ($e->nbeg >= $translation->nend
                                  &&  $e->nend <= $translation->nend);
                }
        }
        die "logic error in CGL::Annotation::Feature::Transcript::lastCodingExon\n";
	# XXXX is there a better way to handle error case?
}

################################################ subroutine header begin ##

=head2 translationStartInTranscript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $i;			# a location in the transcript.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $i = $t->translationStartInTranscript('385099');

=for example end

=for example_testing
  is($i, 141, "Check start in transcript.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub translationStartInTranscript {
	my $self = shift;
	my $id   = shift;

	return $self->{translationStartInTranscript}->{$id}
	  if defined($self->{translationStartInTranscript}->{$id});

	my $translation = $self->fetchTranslationById($id);

	my $first_coding_exon = $self->firstCodingExon($id);

	my $translation_start_offset_in_first_coding_exon =
	abs($first_coding_exon->nbeg - $translation->nbeg);

	my $translationStartInTranscript =
	$first_coding_exon->metaPos($self,
	                            $translation_start_offset_in_first_coding_exon,
	                           );

	$self->{translationStartInTranscript}->{$id} = $translationStartInTranscript;

	return $self->{translationStartInTranscript}->{$id};
}

################################################ subroutine header begin ##

=head2 translationEndInTranscript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $i;			# a location in the transcript.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $i = $t->translationEndInTranscript('385099');

=for example end

=for example_testing
  is($i, 477, "Check end in transcript.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub translationEndInTranscript {
        my $self = shift;
        my $id   = shift;

        return $self->{translationEndInTranscript}->{$id}
        if defined($self->{translationEndInTranscript}->{$id});

        my $translation = $self->fetchTranslationById($id);

        my $last_coding_exon = $self->lastCodingExon($id);


        my $translation_end_offset_in_last_coding_exon =
        abs($last_coding_exon->nbeg - $translation->nend);

        my $translationEndInTranscript =
        $last_coding_exon->metaPos($self,
                                    $translation_end_offset_in_last_coding_exon,
                                   );

        $self->{translationEndInTranscript}->{$id} = $translationEndInTranscript;

        return $self->{translationEndInTranscript}->{$id};
}

################################################ subroutine header begin ##

=head2 metaPos

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $c;			# reference to a contig feature
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $i_e;			# a position in an exon
  my $i_p;			# a position in a protein
  my $i_c;			# a position in a contig

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $c = $a->contig(0);
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->lastCodingExon('385099');
  $p = $t->translation(0);
  
  $i_e = $t->metaPos($e, 322);
  $i_p = $t->metaPos($p, 322);
  $i_c = $t->metaPos($c, 0);

=for example end

=for example_testing
  is($i_e, 0, "Check position in the exon.");
  is(int($i_p), 60, "Check position in the translation (handle fractional frame).");
  is($i_c, 6, "Check position in the contig.");

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

	my $so = $self->_so();

	my $where_in_feature;
	if ($so->a_is_hyponym_of_b($what_f->type, 'exon')){
		return undef
		  unless $where_in_self >= $self->exonStartInTranscript($what_f);

		return undef
		  unless $where_in_self <= $self->exonEndInTranscript($what_f);

		return $where_in_self - $self->exonStartInTranscript($what_f);
	}
	elsif (($what_f->type() eq 'protein') ||
	       ($what_f->type() eq 'polypeptide') ||
	       ($so->a_is_hypomeronym_of_b($what_f->type, 'protein')) ||
	       ($so->a_is_hypomeronym_of_b($what_f->type, 'polypeptide'))){
		my $translation_start_in_transcript =
		$self->translationStartInTranscript($what_f->id);

		my $translation_end_in_transcript =
		$self->translationEndInTranscript($what_f->id);

		return undef
		unless $translation_start_in_transcript <= $where_in_self;

		return undef
		unless $translation_end_in_transcript >= $where_in_self;
	
		$where_in_feature =
		($where_in_self - $self->translationStartInTranscript($what_f->id))/3;
		
		return $where_in_feature;
	}
	elsif ($so->a_is_hyponym_of_b($what_f->type, 'contig')){
                my $i = 0;
                while (my $e = $self->exon($i)){
                        my $e_pos = $self->metaPos($e, $where_in_self);
			return $e->metaPos($what_f, $e_pos)
			if defined($e_pos);
                        $i++;
                }
	}
	else {
		die "this metaPos not supported yet!\n";
	}

}

################################################ subroutine header begin ##

=head2 startCodon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $codon;			# the start codon

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $codon = $t->startCodon('385099');

=for example end

=for example_testing
  is($codon, "ATG", "Check the protein's start codon.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub startCodon {
	my $self = shift;
	my $id   = shift;

	return substr($self->residues, $self->translationStartInTranscript($id), 3);
}

################################################ subroutine header begin ##

=head2 stopCodon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $codon;			# the start codon

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $codon = $t->stopCodon('385099');

=for example end

=for example_testing
  is($codon, "TAA", "Check the protein's stop codon.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub stopCodon {
        my $self = shift;
        my $id   = shift;

        return substr($self->residues, $self->translationEndInTranscript($id), 3);
}

################################################ subroutine header begin ##

=head2 spliceJunction

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $splice;			# the splice junction.


  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $splice = $t->spliceJunction(0);

=for example end

=for example_testing
  is($splice, 192, "Check the protein's stop codon.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub spliceJunction {
	my $self = shift;
	my $i    = shift;

	if (defined($self->exon($i+1))){
		my $e = $self->exon($i);
		
		return $e->metaPos($self, $e->length);
	}
	else {
		return undef;
	}
}

################################################ subroutine header begin ##

=head2 exonJunction

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e_0;			# reference to an exon
  my $e_1;			# ditto...
  my $junction;			# an exon junction

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e_0 = $t->exon(0);
  $e_1 = $t->exon(1);
  $junction = $t->exonJunction($e_0, $e_1);

=for example end

=for example_testing
  is($junction, 192, "Check an exon junction.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonJunction {
	my $self = shift;
	my $e_0  = shift;
	my $e_1  = shift;

	return undef unless defined($e_0) && defined($e_1);
	$self->_add_exon_junctions() unless defined($self->{exonJunctions});
	
	return $self->{exonJunctions}->{$e_0->id}->{$e_1->id};
}

################################################ subroutine header begin ##

=head2 exon

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to an exon feature.
  my $id;			# the exon id.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e = $t->exon(0);
  $id = $e->id();

=for example end

=for example_testing
  isa_ok($e, "CGL::Annotation::Feature::Exon", "Check its type.");
  is($id, "1277106", "Check the exon's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exon {
	my $self = shift;
	my $i    = shift;

	return $self->exons->[$i];
}

################################################ subroutine header begin ##

=head2 exons

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e_ref;			# reference to an array of exon features.
  my $id;			# the exon id.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e_ref = $t->exons();
  $id = $e_ref->[0]->id();

=for example end

=for example_testing
  isa_ok($e_ref->[0], "CGL::Annotation::Feature::Exon", "Check its type.");
  is($id, "1277106", "Check the exon's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exons {
	my $self = shift;

	if (defined($self->{_sorted_exons})) {
        	return($self->{_sorted_exons});
	}

	my @exons = @{$self->{exons}};

	_bubbleSort(\@exons);

	if ($self->strand == 1){
        	$self->{_sorted_exons} = \@exons || [];
	}
	else {
		$self->{_sorted_exons} = [reverse @exons] || [];
	}
	return $self->{_sorted_exons};
}

################################################ subroutine header begin ##

=head2 _bubbleSort

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

sub _bubbleSort {

	my $array = shift;
	
	my $i;
	my $j;

	for ($i = $#{$array}; $i; $i--){
		for ($j = 1; $j <= $i; $j++){
			if ($array->[$j-1]->nbeg()  > $array->[$j]->nbeg()){
				@{$array}[$j, $j-1] = @{$array}[$j-1, $j];
			}
		}
	}
}

################################################ subroutine header begin ##

=head2 translations

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $p_ref;			# reference to an array of protein features
  my $id;			# reference to an array of protein features

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $p_ref = $t->translations();
  $id = $p_ref->[0]->id();

=for example end

=for example_testing
  isa_ok($p_ref, "ARRAY", "Check its type.");
  isa_ok($p_ref->[0], "CGL::Annotation::Feature::Protein", "Check an element's type.");
  is($id, "385099", "Check the protein's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub translations {
        my $self = shift;

       return $self->{translations};
}

################################################ subroutine header begin ##

=head2 intron

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $i;			# reference to an intron feature.
  my $id;			# the exon id.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $i = $t->intron(0);
  $id = $i->id();

=for example end

=for example_testing
  isa_ok($i, "CGL::Annotation::Feature::Intron", "Check its type.");
  is($id, "CG4852-region-ext1.3:-186:-283", "Check the exon's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub intron {
        my $self = shift;
        my $i    = shift;

        return $self->introns->[$i];
}

################################################ subroutine header begin ##

=head2 introns

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $i_ref;			# reference to an array of intron features.
  my $id;			# the exon id.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $i_ref = $t->introns();
  $id = $i_ref->[0]->id();

=for example end

=for example_testing
  isa_ok($i_ref, "ARRAY", "Check its type.");
  isa_ok($i_ref->[0], "CGL::Annotation::Feature::Intron", "Check its type.");
  is($id, "CG4852-region-ext1.3:-186:-283", "Check the exon's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub introns {
        my $self = shift;

	return [] unless defined($self->{introns});

        my @introns = sort {$a->nbeg <=> $b->nbeg} @{$self->{introns}};

        if ($self->strand == 1){
                return \@introns || []
        }
        else {
                return [reverse @introns] || [];
        }
}

################################################ subroutine header begin ##

=head2 exonPairs

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature
  my $t;			# reference to a transcript feature
  my $e_pairs_ref;		# reference to an array of refs to exon pairs.
  my $e_pair_ref;		# reference to a pair of exons
  my ($e1, $e2);		# the exon pair.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(1);
  $t = $g->transcript(0);
  $e_pairs_ref = $t->exonPairs();
  $e_pair_ref = $e_pairs_ref->[0]; 
  ($e1, $e2) = @{$e_pair_ref};

=for example end

=for example_testing
  isa_ok($e1, "CGL::Annotation::Feature::Exon", "Check its type.");
  isa_ok($e2, "CGL::Annotation::Feature::Exon", "Check its type.");
  is($e1->id(), "1277106", "Check the exon's id.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonPairs {
	my $self = shift;

	my @exons = @{$self->exons};

	my @pairs;
	for (my $i = 0; $i < @{$self->exons} -1; $i++){
		my $j = $i+1;
		push(@pairs, [$self->exon($i), $self->exon($j)]);
	}
	return \@pairs || [];
}
#-------------------------------------------------------------------------------
#------------------------------- PRIVATE ---------------------------------------
#-------------------------------------------------------------------------------

################################################ subroutine header begin ##

=head2 _add_residues_2

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

sub _add_residues_2 {
        my $self    = shift;

	return if defined($self->residues);

        my $nB = $self->nbeg();
        my $nE = $self->nend();

	if ($nB < 0 || $nE < 0){
		$self->{_residues} = undef;

	}
        my $residues  ='';
        foreach my $e (@{$self->exons}){
		unless (defined($e->residues)){
			$self->{residues} = undef;
			return;
		}
                $residues .= $e->residues();
        }

	if (defined($self->residues) && $residues ne $self->residues){
		print "stated residues:".$self->residues."\n";
		print "residues:".$residues."\n";
		die "dead in Transcript::_add_residues_2\n";
	}

        $self->{residues} = $residues;

}

################################################ subroutine header begin ##

=head2 _add_residues

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

sub _add_residues {
        my $self    = shift;
        my $contig = shift;

        my $nB = $self->nbeg();
        my $nE = $self->nend();

        return undef if ($nB < 0 || $nE < 0);
        return undef if ($nB > $contig->length || $nE > $contig->length);

	my $residues  ='';
	foreach my $e (@{$self->exons}){
		$residues .= $e->residues();
	}

        if    ($contig->strand() == 1 && $self->strand() == 1){
                $self->residues($residues);
        }
        elsif ($contig->strand() == 1 && $self->strand() == -1){
                $self->residues(revcomp($residues));
        }
        elsif ($contig->strand() == -1 && $self->strand() == 1){
                $self->residues($residues);
        }
        elsif ($contig->strand() == -1 && $self->strand() == -1){
                $self->residues(revcomp($residues));
        }

        else {
                die "contig strand is -1 in Feature::Transcript::_add_residues!\n";
        }

}

################################################ subroutine header begin ##

=head2 _add_exon_junctions

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

sub _add_exon_junctions {
        my $self = shift;
	
	my %j;
	my $t_natural_offset = 0;
	foreach my $p (@{$self->exonPairs}){
		my $e_0 = $p->[0];
		my $e_1 = $p->[1];
		
		$t_natural_offset += $e_0->length();

		$j{$e_0->id}{$e_1->id} = $t_natural_offset;
	}

	$self->{exonJunctions} = \%j;
}

################################################ subroutine header begin ##

=head2 _add_exon

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

sub _add_exon {
        my $self = shift;
        my $exon = shift;

        push(@{$self->{exons}}, $exon);

	# see sub exons()
	$self->{_sorted_exons} = undef;
}

################################################ subroutine header begin ##

=head2 _add_translation

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

sub _add_translation {
	my $self        = shift;
	my $translation = shift;

	push(@{$self->{translations}}, $translation);
}

################################################ subroutine header begin ##

=head2 _add_introns

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

sub _add_introns {
	my $self    = shift;
	my $contig  = shift;

	return if defined($self->{introns});

	my $pairs = $self->exonPairs();
	my $i = 0;
	my $k = 1;
	foreach my $pair (@{$pairs}){
		my $f = {};
		my $eI = $pair->[0];
		my $eJ = $pair->[1];

		my ($name) = $eI->name() =~ /(\S+)\:\d+/;
		    $name .= ":intron:$k";
		
		$f->{dbxref}     =  $contig->uniquename.":".$eI->metaPos($contig, $eI->length).":".$eJ->metaPos($contig, 0);
		$f->{name}       = $name;
		$f->{residues}   = '';
		$f->{seqlen}     = '';
		$f->{strand}     = $eI->strand();
		$f->{type}       = 'intron';

		my %f_location;
		
		$f_location{srcfeature_id} = $eI->featureLocation(0)->srcfeature_id();
		$f_location{strand}     = $eI->strand();
		$f_location{rank}       = undef;
		
		if ($eI->strand == 1){
			$f_location{nbeg} = $eI->nend;
			$f_location{nend} = $eJ->nbeg;
		}
		else {
			$f_location{nbeg} = $eI->nend;
			$f_location{nend} = $eJ->nbeg;	
		}

		push(@{$f->{locations}}, new CGL::Annotation::FeatureLocation(\%f_location));

		my $intron = new CGL::Annotation::Feature::Intron($f);
		   $intron->_add_residues($eI, $eJ, $contig) if defined($contig);

		foreach my $r ($eI->relationships){
			my %hash;
			$hash{feature_relationship}{objfeature}  = $r->oF();
			$hash{feature_relationship}{subjfeature} = $f->{dbxref};
			$hash{feature_relationship}{type}        = $r->logus();

			my $rel = new CGL::Annotation::FeatureRelationship(\%hash);
			$intron->_add_relationship($rel);
		}

		$intron->_add_src_id($contig, $eI, $eJ);

		$self->_add_intron($intron);
		$i++;
		$k++;
	}
}

################################################ subroutine header begin ##

=head2 _add_intron

 Usage     : *private

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _add_intron {
        my $self   = shift;
        my $intron = shift;

        push(@{$self->{introns}}, $intron);
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
	if ($ENV{CGL_CHATTER}){
	    print STDERR "CGL::Annotation::Feature::Transcript AutoLoader called for: ",
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

