###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Gene

CGL::Annotation::Feature::Gene - The platonic ideal of a CGL module.

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

package CGL::Annotation::Feature::Gene;

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

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation::Feature::Exon;

  my $feature;			# the raw CGL Feature
  my $g;			# the gene object.

  $feature = {};		# In real life, a CGL::Annotation::Feature...
  $g = new CGL::Annotation::Feature::Gene($feature);

=for example end

=for example_testing
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Check its type.");

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

=head2 transcript

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an exon feature
  my $t;			# reference to an exon feature

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);

=for example end

=for example_testing
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Check the gene's type.");
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check transcript type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub transcript
 {
        my $self = shift;
        my $i    = shift;

        return $self->transcripts->[$i];
}

################################################ subroutine header begin ##

=head2 nbeg

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an exon feature
  my $t;			# reference to an transcript feature
  my $begin;			# the beginning

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $begin = $t->nbeg();

=for example end

=for example_testing
  is($begin, 500, "Check the gene's natural begin.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nbeg
 {
	my $self = shift;

	return $self->transcripts->[0]->nbeg;
}

################################################ subroutine header begin ##

=head2 nend

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an exon feature
  my $t;			# reference to an transcript feature
  my $end;			# the end.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $end = $t->nend();

=for example end

=for example_testing
  is($end, 2288, "Check the gene's natural end.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nend
 {
	my $self = shift;

        return $self->transcripts->[-1]->nend;

}

################################################ subroutine header begin ##

=head2 transcripts

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an exon feature
  my $t;			# ref to transcript array

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcripts();

=for example end

=for example_testing
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Did we get a gene.");
  isa_ok($t->[0], "CGL::Annotation::Feature::Transcript", "Did we get transcripts.");
  is(scalar(@{$t}), 1, "Did we get the right number of transcripts?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub transcripts
 {
        my $self = shift;


	unless (defined($self->{transcripts})){
		print STDERR "\n";
		print STDERR "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
		print STDERR "WARNING WARNING NO TRANSCRPTS FOR THIS GENE! ".$self->name."\n";
		print STDERR "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
		print STDERR "\n";

		#sleep(5);

		return [];

	}

        my @trans = sort {$a->nbeg <=> $b->nbeg} @{$self->{transcripts}};

        if ($self->{transcripts}->[0]->strand == 1){
                return \@trans
        }
        else {
                return [reverse @trans];
        }
}

################################################ subroutine header begin ##

=head2 transcript_with_id

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an exon feature
  my $t;			# reference to an transcript feature
  my $begin;			# the beginning

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript_with_id("mRNA-107699");

=for example end

=for example_testing
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check the type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub transcript_with_id
 {
	my $self = shift;
	my $id   = shift;

	my $i = 0;
	while (my $t = $self->transcript($i)){
		return $t if $t->id eq $id;
		$i++;
	}
}

################################################ subroutine header begin ##

=head2 bubbleSort

 Usage     : unused.

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
 See Also  :

=cut

################################################## subroutine header end ##

sub bubbleSort
 {

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
#-------------------------------------------------------------------------------
#------------------------------- PRIVATE ---------------------------------------

################################################ subroutine header begin ##

=head2 _add_transcript

 Usage     : *private*

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

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
	    print STDERR "CGL::Annotation::Feature::Gene AutoLoader called for: ",
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

