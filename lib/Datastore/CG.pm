package Datastore::CG;

use 5.006;
use strict;
use warnings;

use Datastore::Base;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter Datastore::Base);

our $VERSION = '0.01';

# given an identifier, transform it into a directory name.
sub id_to_dir {
  my($self, $id) = @_;
  my($dir);

  # snag the first three digits of the CG name.
  ($dir) = ($id =~ m|CG(\d\d\d).*|);

  # build up the directory path.
  $dir = $self->{_root} . $dir . "/" . $id . "/";

  return($dir);
}

# given an directory name, transform it into an identifier.
sub dir_to_id {
  my($self, $dir) = @_;
  my($id);
  
  $id = $dir;
  $id =~ s|(.*)/$|$1|;
  $id =~ s|.*/(.*)$|$1|;

  return($id);
}

1;
__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Datastore::CG - foo

=head1 SYNOPSIS

  use Datastore::CG;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Datastore::CG, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.


=head1 AUTHOR

A. U. Thor, E<lt>a.u.thor@a.galaxy.far.far.awayE<gt>

=head1 SEE ALSO

L<perl>.

=cut
