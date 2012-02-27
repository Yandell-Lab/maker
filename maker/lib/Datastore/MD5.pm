package Datastore::MD5;

use 5.006;
use strict;
use warnings;
use Digest::MD5 qw(md5_hex);
use Datastore::Base;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter Datastore::Base);

our $VERSION = '0.01';

sub id_to_dir {
  my($self, $id) = @_;
  my($dir);
  my($i);
  my($digest) = uc md5_hex($id); # the hex string, uppercased

  # make sure we have enough hex digits to build the dir string.
  if ($self->depth() > 16) {
    return(undef);
  }

  $dir = $self->{_root};
  for($i = 0; $i < $self->{_depth}; $i++) {
    $dir .= substr($digest, $i*2, 2) . "/";
  }
  $dir .= $id . "/";

  return($dir);
}

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

Datastore::MD5 - A datastore implementation with an MD5 digest dir structure.

=head1 SYNOPSIS

  use Datastore::MD5;
  
  # create/user a new datastore object, rooted at /scratch/space
  # and using a two level directory hiearachy.
  $ds = new Datastore::MD5("root" => "/scratch/space", "depth" => 2);

  # convert an identifier into a datastore directory.
  $dir = $ds->id_to_dir("CG0666");

  # and convert the name of a directory in a datastore back into an 
  # identifier
  $id = $ds->dir_to_id($dir);

  # temporarily change to the directory for "IDENTIFIER" and run 
  # a command.
  $ds->system("CG9876", "sim4 input.fa output.fa > sim4.output");

  # change the process's current working directory to the one for
  # "CG9876"
  $ds->chdir("CG9876");

=head1 DESCRIPTION

Create and manage a hierarchy of directories that make up a datastore.
This module uses pairs of digits from an MD5 digest of an identifier
(e.g. "CG1234") to create a hierarchy of directories that avoid the
problems caused by a attempting to create a very large number of
directories in a single directory.

For example, given a datastore rooted in "/tmp", and an IDENTIFIER of
"CG9876", the corresponding directory would be "/tmp/E2/96/CG9876/".

Most of the actual operations are performed by functions in the parent
class Datastore::Base.

=head2 METHODS

=head3 id_to_dir

  Convert an identifier to the corresponding name of a directory in the datastore.
  Returns undef on failure and the directory name on success.

=head3 dir_to_id

  Convert the name of a directory in a datastore to the corresponding identifier.
  Should never fail (famous last words, so no error info available...).

=head1 AUTHOR

George Hartzell, E<lt>hartzell @ fruitfly.orgE<gt>

=head1 SEE ALSO

L<Datastore::Base>, L<perl>.

=cut
