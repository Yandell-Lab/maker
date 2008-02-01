###################################################### main header begin ##

=head1 NAME

CGL::Version

=head1 SYNOPSIS

=for example begin

  use CGL::Version;
  print "CGL Version " . CGL::Version::version();

=for example end

=for example_testing
  is( $_STDOUT_, "CGL Version 0.08", "Checking explicit call.");

=head1 DESCRIPTION

Provide the current version for the CGL package.

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

package CGL::Version;

use strict;
use warnings;

use vars qw( $VERSION @EXPORT_OK );
BEGIN {
    $VERSION     = 0.08;
    @EXPORT_OK   = qw( version );
    # don't pull in Way Too Much from Exporter;
    require Exporter;
    *import = \&Exporter::import;
}

################################################ subroutine header begin ##

=head2 version

 Usage     :

=for example begin

  # refer to the version function explicitly.
  use CGL::Version;
  print "CGL Version " . CGL::Version::version();

=for example end

=for example_testing
  is( $_STDOUT_, "CGL Version 0.08", "Checking for correct version number.");

=for example begin

  # import the version function into your namespace (yikes!).
  use CGL::Version qw( version );
  print "CGL Version " . version();

=for example end

=for example_testing
  is( $_STDOUT_, "CGL Version 0.08", "Checking exported call.");



 Purpose   : Provide access to the version number for the CGL library.
 Returns   : The CGL version.
 Argument  : none.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub version
{
	return ($CGL::Version::VERSION);
}

################################################ subroutine header begin ##
# ....
# ....
# ....
################################################## subroutine header end ##

1; #this line is important and will help the module return a true value
__END__

