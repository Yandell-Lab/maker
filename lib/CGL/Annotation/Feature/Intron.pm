###################################################### main header begin ##

=head1 CGL::Annotation::Feature::Intron

CGL::Annotation::Feature::Intron - The platonic ideal of a CGL module.

XXXX NOT FINISHED, NO KNOWN EXAMPLES OF CREATING THIS....

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

package CGL::Annotation::Feature::Intron;

use strict;
use warnings;

use CGL::Annotation::Feature;
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

  use CGL::Annotation::Feature::Intron;

  my $feature;			# the raw CGL Feature
  my $i;			# the intron object

  $feature = {};		# It's actually a CGL::Annotation::Feature
  $i = new CGL::Annotation::Feature::Intron($feature);

=for example end

=for example_testing
  isa_ok($i, "CGL::Annotation::Feature::Intron", "Check its type.");

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

	$self->_add_id('dbxref');

	return $self;
}

################################################ subroutine header begin ##

=head2 get_begin_end_on_src

 Usage     : How to use this function/method

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an gene feature
  my $t;			# reference to a transcript feature
  my $e;			# reference to array of exon features

  $a = new CGL::Annotation("sample_data/dmel.sample.chaos.xml");
  $g = $a->gene(0);
  $t = $g->transcript(0);
  $e = $t->exons();

=for example end

=for example_testing
  isa_ok($g, "CGL::Annotation::Feature::Gene", "Check the gene's type.");
  isa_ok($t, "CGL::Annotation::Feature::Transcript", "Check transcript type.");
  isa_ok($e->[0], "CGL::Annotation::Feature::Exon", "Check the exon's type.");
  is(scalar(@{$e}), 2, "Check the number of exons.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub get_begin_end_on_src
 {
        my $self = shift;
        my $c    = shift;
	my $e_i  = shift;
	my $e_j  = shift;

        my $b_on_c = $e_i->metaPos($c, $e_i->length());
        my $e_on_c = $e_j->metaPos($c, 0);

        my $src_b = $c->location(0)->nbeg();
        my $src_e = $c->location(0)->nend();

        my $b_on_s;
        my $e_on_s;

        if ($src_b < $src_e){
                $b_on_s = $b_on_c + $src_b;
                $e_on_s = $e_on_c + $src_b;
        }
        else {
                $b_on_s = $src_b - $b_on_c;
                $e_on_s = $src_b - $e_on_c;
        }
        return ($b_on_s, $e_on_s);

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

# XXXX common code w/ CGL::Annotation::Feature::_add_residues()

sub _add_residues
 {
        my $self    = shift;
	my $eI      = shift;
	my $eJ      = shift;
        my $contig = shift;

        my $nB = $self->nbeg();
        my $nE = $self->nend();

	if ( $nB eq "" || $nE eq "" || $nB < 0 || $nE < 0 ){
		$self->inScope(0);
		return undef;
	}
	if ($nB > $contig->length || $nE > $contig->length){
		 $self->inScope(0);
		 return undef;
	}

        if    ($contig->strand() == 1 && $self->strand() == 1){
		my $l = $eJ->nbeg - $eI->nend;
		my $r = substr($contig->residues(), $eI->nend(), $l);
                $self->residues($r);
		$self->inScope(1);
        }
        elsif ($contig->strand() == 1 && $self->strand() == -1){
		my $l = $eI->nend - $eJ->nbeg;
		my $r = substr($contig->residues(), $eJ->nbeg, $l);
                $self->residues(revcomp($r));
		$self->inScope(1);
        }
        elsif ($contig->strand() == -1 && $self->strand() == 1){
                my $l = $eJ->nbeg - $eI->nend;
                my $r = substr($contig->residues(), $eI->nend(), $l);
                $self->residues($r);
                $self->inScope(1);
        }
        elsif ($contig->strand() == -1 && $self->strand() == -1){
                my $l = $eI->nend - $eJ->nbeg;
                my $r = substr($contig->residues(), $eJ->nbeg, $l);
                $self->residues(revcomp($r));
                $self->inScope(1);
        }

        else {
                die "contig strand is -1 in Feature::Intron::_add_residues!\n";
        }

}

################################################ subroutine header begin ##

=head2 _add_src_id

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

sub _add_src_id
 {
        my $self = shift;
        my $c    = shift;
	my $e_i  = shift;
	my $e_j  = shift;

        my ($src_b, $src_e) = $self->get_begin_end_on_src($c, $e_i, $e_j);

        $c->location(0)->srcfeature_id();

        my $src_id = $c->location(0)->srcfeature_id().":".$src_b.":".$src_e;

        $self->src_id($src_id);

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
            print STDERR "CGL::Annotation::Feature::Intron::AutoLoader called for: ",
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

