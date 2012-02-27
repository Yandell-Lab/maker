###################################################### main header begin ##

=head1 CGL::Annotation::Iterator

CGL::Annotation::Iterator - A mechanism for iterating across annotations.

=head1 SYNOPSIS

=for example

=for example begin

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $ref;			# an array reference
  my $g;			# a gene
  my $t;			# a transcript
  my $p;			# a translation (protein...)

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);
  $g = $iter->next_by_gene();
  $ref = $iter->next_by_transcript;
  ($t, $g) = @{$ref};
  $ref = $iter->next_by_translation;
  ($p, $t, $g) = @{$ref};

=for example end

=for example_testing
  isa_ok($iter, "CGL::Annotation::Iterator", "Check if it's the right type.");
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Check its type.");
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check its type.");
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check its type.");

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

package CGL::Annotation::Iterator;

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

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $g;			# a gene
  my $t;			# a transcript
  my $p;			# a translation (protein...)

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);

=for example end

=for example_testing
  isa_ok($iter, "CGL::Annotation::Iterator", "Check if it's the right type.");

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
  my $class      = shift;
  my $annotation = shift;

  my $self = {};
  bless $self, $class;

  $self->_annotation($annotation);
  $self->reset();

  return $self;
}

################################################ subroutine header begin ##

=head2 reset

 Usage     :

=for example begin

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $g_1;			# a gene
  my $g_2;			# another gene (should be the same one...)

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);

  $g_1 = $iter->next_by_gene();
  $iter->reset();
  $g_2 = $iter->next_by_gene();

=for example end

=for example_testing
  is($g_1, $g_2, "Check reset.");

 Purpose   : Reset the iterator to the beginning of the annotation.
 Returns   : Nothing.
 Argument  : None.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub reset {
  my $self = shift;

  $self->{_i} = 0;		# iterates over genes in an annotation.
  $self->{_j} = 0;		# iterates over transcripts in a gene in...
  $self->{_k} = 0;		# iterates over translations in a transcript...
}

################################################ subroutine header begin ##

=head2 next_by_gene

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $g;			# a gene
  my $count = 0;		# counts the genes

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);
  $g = $iter->next_by_gene();

  $iter->reset();
  while ($iter->next_by_gene()) {
    $count++;
  }

=for example end

=for example_testing
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Check that returned a gene.");
  is($count, 2, "Check that it counted right.");

 Purpose   : Return a reference to the next gene in the annotation.
 Returns   : A reference to a CGL::Annotation::Feature::Gene, or undef
           : if there are no more genes.
 Argument  : None.
 Throws    : None.
 Comments  : Updates the iterator so that the next transcript or translation
           : will come from this gene.
 See Also  : CGL::Annotation::Feature::Gene

=cut

################################################## subroutine header end ##

sub next_by_gene {
  my $self = shift;

  my $g = undef;		# the next gene.

  goto done if (! $self->_valid());

  $g = $self->_annotation->gene($self->{_i});

  if (!defined($g)) {		# no more genes...
    $self->_invalidate();
    goto done;
  }

  $self->{_i}++;		# bump the gene counter
  $self->{_j} = 0;		# start at the beginning of *this* gene
  $self->{_k} = 0;		# start at the beginning of *this* gene
	
done:
  return $g;
}

################################################ subroutine header begin ##

=head2 next_by_transcript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $ref;			# results ref.
  my $g;			# a gene
  my $t;			# a transcript
  my $count = 0;		# counts the transcripts

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);
  $ref = $iter->next_by_transcript();
  ($t, $g) = @{$ref};

  $iter->reset();
  while ($iter->next_by_transcript()) {
    $count++;
  }

=for example end

=for example_testing
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check its type.");
  is($count, 3, "Check that it counts right.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub next_by_transcript {
  my $self = shift;

  my $g;			# a gene
  my $t;			# a transcript
  my $result;			# the result...

  $g = $self->_annotation->gene($self->{_i});

  if (!defined($g)) {		# no more genes, give up.
    $self->{_i} = undef;
    $self->{_j} = undef;
    $self->{_k} = undef;
    $result = undef;
    goto done;
  }

  $t = $g->transcript($self->{_j});

  if (!defined($t)){		# no more transcripts in this gene.
    $self->{_i}++;		# try the next one...
    $self->{_j} = 0;
    $self->{_k} = 0;
    $result = $self->next_by_transcript();
    goto done;
  }

  $self->{_j}++;		# phew, found one.
  $self->{_k} = 0;
  $result = [$t, $g];

 done:
  return $result;
}

################################################ subroutine header begin ##

=head2 next_by_translation

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  use CGL::Annotation::Iterator;

  my $a;			# reference to a CGL::Annotation
  my $iter;			# iterator
  my $ref;			# results ref.
  my $g;			# a gene
  my $t;			# a transcript
  my $p;			# a translation (protein)
  my $count = 0;		# counts the transcripts

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $iter = new CGL::Annotation::Iterator($a);
  $ref = $iter->next_by_translation();
  ($p, $t, $g) = @{$ref};

  $iter->reset();
  while ($iter->next_by_transcript()) {
    $count++;
  }

=for example end

=for example_testing
  isa_ok($p, "CGL::Annotation::Feature::Protein", "Check the type.");
  is($count, 3, "Check that it counts right.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub next_by_translation {
  my $self = shift;

  my $g;			# a gene
  my $t;			# a transcript
  my $p;			# a protein (translation...)
  my $result;			# the answer

  $g = $self->_annotation->gene($self->{_i});

  if (!defined($g)) {		# no more genes, give up.
    $self->{_i} = undef;
    $self->{_j} = undef;
    $self->{_k} = undef;
    $result = undef;
    goto done;
  }

  $t = $g->transcript($self->{_j});

  if (!defined($t)){		# no more transcripts in this gene.
    $self->{_i}++;		# try the next one.
    $self->{_j} = 0;
    $self->{_k} = 0;
    $result = $self->next_by_translation();
    goto done;
  }

  $p = $t->translation($self->{_k});

  if (!defined($p)){		# no more translations for this transcript
    $self->{_j}++;		# try the next one
    $self->{_k} = 0;
    $result = $self->next_by_translation();	
    goto done;
  }

  $self->{_k}++;		# phew, found one...
  $result = [$p, $t, $g];

 done:
  return $result;
}

################################################ subroutine header begin ##

=head2 _annotation

 Usage     : *private*

 Purpose   : get/set the iterator's annotation.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub _annotation {
  my $self  = shift;

  if (@_) {
    $self->{annotation} = shift;
  }
  return $self->{annotation};
}

################################################ subroutine header begin ##

=head2 _invalidate

 Usage     : *private*

 Purpose   : set the internal counters to invalid values
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _invalidate {
  my $self = shift;

  $self->{_i} = undef;	# iterates over genes in an annotation.
  $self->{_j} = undef;	# iterates over transcripts in a gene in...
  $self->{_k} = undef;	# iterates over translations in a transcript...
}

################################################ subroutine header begin ##

=head2 _valid

 Usage     : *private*

 Purpose   : returns true if the iterator is hasn't run off the end of the
           : annotation.
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _valid {
  my $self = shift;

  return defined($self->{_i});
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

sub AUTOLOAD {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

	if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::Iterator::AutoLoader called for: ",
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

