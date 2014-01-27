###################################################### main header begin ##

=head1 NAME

Revcomp.pm - A simple module to reverse complement a sequence.

=head1 SYNOPSIS

=for example begin

  use CGL::Revcomp qw(revcomp);
  my $s1 = revcomp("acgt");
  my $s2 = revcomp("AttackAKayak");

=for example end

=for example_testing
  is(revcomp("acgt"), "acgt", "Check simple dna string");
  is(revcomp("AttackAKayak"), "mtrtMTmgtaaT",
     "Check funky ambiguity and case.");

=head1 DESCRIPTION

This module contains a single routine that reverse complements a
sequence.

=head1 USAGE

See the SYNOPSIS section.

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

package CGL::Revcomp;

use strict;
use warnings;

# subroutine exporting version
use vars qw( $VERSION @EXPORT_OK %EXPORT_TAGS );
BEGIN {
    $VERSION     = 0.01;
    @EXPORT_OK   = qw( revcomp );
    # %EXPORT_TAGS = ( ... );
    # don't pull in Way Too Much from Exporter;
    require Exporter;
    *import = \&Exporter::import;
}

################################################ subroutine header begin ##

=head2 revcomp

 Usage     :

=for example begin

  use CGL::Revcomp qw(revcomp);
  my $s1 = revcomp("acgt");
  my $s2 = revcomp("AttackAKayak");

=for example end

=for example_testing
  is(revcomp("acgt"), "acgt", "Check simple dna string");
  is(revcomp("AttackAKayak"), "mtrtMTmgtaaT",
     "Check funky ambiguity and case.");

 Purpose   : Reverse complements a sequence.
 Returns   : The reverse complement of its argument, in a scalar.
 Argument  : A scalar variable containing the sequence.
 Throws    : 
 Comments  : 
 See Also  : 

=cut

################################################## subroutine header end ##

sub revcomp
 {
   my $seq = shift;

   # XXXX I added the final V->B since there was only B->V before....
   $seq =~ tr/AaCcGgTtYyRrKkMmBbVvUuSsWwDdHh/TtGgCcAaRrYyMmKkVvBbAaWwSsHhDd/;

   return reverse $seq;
}

1; #this line is important and will help the module return a true value
__END__
