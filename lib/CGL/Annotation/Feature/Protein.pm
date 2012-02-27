###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Protein

CGL::Annotation::Feature::Protein - The platonic ideal of a CGL module.

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

package CGL::Annotation::Feature::Protein;

use strict;
use warnings;

use CGL::Annotation::Feature;
use CGL::Clone qw(clone);
use CGL::TranslationMachine();

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

  use CGL::Annotation::Feature::Contig;

  my $feature;			# the raw CGL Feature
  my $p;			# the protein object

  $feature = {};		# It's actually a CGL::Annotation::Feature
  $p = new CGL::Annotation::Feature::Protein($feature);

=for example end

=for example_testing
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check its type.");

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
	my $feature  = shift;

	
	my $self = clone($feature);

	bless($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 residues

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature
  my $r;			# the protein's residues

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $p = $t->translation(0);
  $r = $p->residues();

=for example end

=for example_testing
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check the type.");
  like($r, qr/^MRLWDND.*/, "Are they the residues we expect.");
  is(length($r), 595, "Check the length of the residues.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub residues
 {
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

=head2 length

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $p = $t->translation(0);

=for example end

=for example_testing
  is($p->length, 595, "Check the length() method.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub length
 {
	my $self = shift;

	return length($self->residues);

}

################################################ subroutine header begin ##

=head2 metaPos

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature
  my $begin_in_t;			# location in the transcript.
  my $end_in_t;			# location in the transcript.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcript(0);
  $p = $t->translation(0);
  $begin_in_t = $p->metaPos($t, 0);
  $end_in_t = $p->metaPos($t, $p->length);

=for example end

=for example_testing
  is($begin_in_t, 0, "Check the start of p in t.");
  is($end_in_t, 1785, "Check the end of p in t.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub metaPos
 {
        my $self          = shift;
        my $what_f        = shift;
        my $where_in_self = shift;

	my $so = $self->_so();

	my $where_in_feature;
	if ($so->a_is_hyponym_of_b($what_f->type, 'transcript')){	
		$where_in_feature =
		$what_f->translationStartInTranscript($self->id)+$where_in_self*3;
	}
	return $where_in_feature;
}

################################################ subroutine header begin ##

=head2 triplet

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# a gene.
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature
  my $start_triplet;		# the proteins first triplet.
  my $final_triplet;		# the proteins last triplet.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);		# XXXX need to go through gene to work!!!
  $t = $g->transcript(0);
  $p = $t->translation(0);
  $start_triplet = $p->triplet($t, 0);
  $final_triplet = $p->triplet($t, $p->length());

=for example end

=for example_testing
  is($start_triplet, "ATG", "Check the first three bases.");
  is($final_triplet, "TAG", "Check the last three bases.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub triplet
 {
	my $self                  = shift;
	my $transcript            = shift;
	my $translation_offset    = shift;

	my $offset_of_codon_in_transcript =
	$self->metaPos($transcript,int($translation_offset));

	return substr($transcript->residues,
		      $offset_of_codon_in_transcript,
		      3);
}

################################################ subroutine header begin ##

=head2 aa

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# a gene.
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature
  my $aa;			# an amino acid

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);		# XXXX need to go through gene to work!!!
  $t = $g->transcript(0);
  $p = $t->translation(0);
  $aa = $p->aa(0);

=for example end

=for example_testing
  is($aa, "M", "Check the first amino acid.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub aa
 {
	my $self              = shift;
	my $offset_on_protein = shift;

	
	return substr($self->residues(), int($offset_on_protein), 1);
}

################################################ subroutine header begin ##

=head2 exonJunction

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# a gene.
  my $t;			# reference to a transcript feature
  my $p;			# reference to a protein feature
  my $offset;			# offset of exon junction on protein.

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);		# XXXX need to go through gene to work!!!
  $t = $g->transcript(0);
  $p = $t->translation(0);
  $offset = $p->exonJunction($t, 0, 1);

=for example end

=for example_testing
  is(int($offset), int(135.333), "Check location of the first exon junction.");
  # deal with is() and floating point numbers....

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub exonJunction
 {
	my $self        = shift;
	my $transcript  = shift;
	my $exon_num_0  = shift;
	my $exon_num_1  = shift;

	my $e_0 = $transcript->exon($exon_num_0);
	my $e_1 = $transcript->exon($exon_num_1);

	my $offset_of_junction_on_transcript =
	$transcript->exonJunction($e_0, $e_1);
	
	my $offset_of_junction_on_protein =
	$transcript->metaPos($self, $offset_of_junction_on_transcript);

	return $offset_of_junction_on_protein;
}

################################################ subroutine header begin ##

=head2 _add_residues_2

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  : XXXX why?
           : XXXX returns various things, depending on how far it goes....
 See Also  :

=cut

################################################## subroutine header end ##

sub _add_residues_2
 {
        my $self    = shift;
	my $t       = shift;

	return if defined($self->residues);

	my $translation_offset = $self->metaPos($t, 0);

        unless (defined($t->residues)  && defined($translation_offset)){
                $self->{residues} = undef;

        }

        my $tM = new CGL::TranslationMachine();

        my $translation = $tM->translate_from_offset($t->residues,
						     $translation_offset);

        ($translation) = $translation =~ /([A-Z]+)/;

	if (defined($self->residues) && $translation  ne $self->residues){
		warn  "problems in Protein::_add_residues_2\n";
	}

        $self->{residues} = $translation;
}


################################################ subroutine header begin ##

=head2 _add_transcript

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

sub _add_transcript
 {
	my $self       = shift;
	my $transcript = shift;

	push(@{$self->{transcripts}}, $transcript);
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

sub AUTOLOAD
 {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::Feature::Protein AutoLoader called for: ",
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

