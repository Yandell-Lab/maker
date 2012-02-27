package Datastore::Base;

use 5.006;
use strict;
use warnings;

require Exporter;

use AutoLoader qw(AUTOLOAD);
use Cwd;
use File::Path;

our @ISA = qw(Exporter);

our $VERSION = '0.11';

sub new {
  my($class) = shift;
  my(%args) = @_;
  my($self) = bless {}, $class;

  # the Base class shouldn't be instantiated.
  if($class eq "Datastore::Base") {
    die "Don't use Datastore::Base directly, use a subclass.\n";
  }

  $self->{_root} = "";
  $self->{_depth} = 2;
  $self->{_readonly} = 0;

  if (defined($args{root})) {
    $self->root($args{"root"});
  }

  if (defined($args{"depth"})) {
    $self->depth($args{"depth"});
  }

  $self->{_dStack} = [];

  return ($self);
}

# set and/or get the root directory of the datastore.
# take care of adding a trailing "/" if necessary.
sub root {
  my($self, $root) = @_;

  if(defined($root)) {
    $self->{_root} = $root;
    # make sure it has a trailing "/"
    if ($self->{_root} !~ m|.*/$|) {
      $self->{_root} .= "/";
    }
  }
  
  return($self->{_root});
}

# set and/or get the root directory of the datastore.
# take care of adding a trailing "/" if necessary.
sub depth {
  my($self, $depth) = @_;

  if(defined($depth)) {
    $self->{_depth} = $depth;
  }
  
  return($self->{_depth});
}

# set and/or get the root directory of the datastore.
# take care of adding a trailing "/" if necessary.
sub readonly {
  my($self, $readonly) = @_;

  if(defined($readonly)) {
    $self->{_readonly} = $readonly;
  }
  
  return($self->{_readonly});
}

sub version {
  my($self) = @_;
  return($VERSION);
}

# recursively make the path to a datastore directory.
# redirect mkdir output to the bit bucket and return
# true on success, false on failure.
sub mkdir {
  my($self, $id) = @_;
  my($dir);
  my($retval);
  my($return);

  $dir = $self->id_to_dir($id);
  if (! $self->{_readonly}) {
    $retval = system("mkdir -p $dir > /dev/null 2>&1");
    $return = ($retval == 0 ? 1 : 0); # convert to boolean...
  }
  else {
    $return = 0;
  }

  return($return);
}

# make a directory if necessary, and chdir to it.
# return true on success, false on failure.
sub chdir {
  my($self, $id) = @_;
  my($dir);
  my($retval) = 1;
  
  $dir = $self->id_to_dir($id);

  if (! -x $dir) {
    $retval = $self->mkdir($id);
  }

  if ($retval) {
    $retval = CORE::chdir($dir);
  }
  
  return($retval);
}

# handy alias for chdir.
sub cd {
  my($self, $id) = @_;
  return($self->chdir());
}

# push the current working directory onto self's stack and cd to $id's datastore.
# return chdir's return value.
sub pushd {
  my($self, $id) = @_;
  
  push @{$self->{_dStack}}, cwd();
  return($self->chdir($id));
}

# pop a directory from the objects stash and cd to it.
# return CORE::chdir's value.
# return chdir's return value.
sub popd {
  my($self, $id) = @_;
  my($newdir);

  $newdir = pop @{$self->{_dStack}};
  
  return(CORE::chdir($newdir));
}


# run a command in the datastore directory that belongs to $id.
# return true on success, zero on failure.
sub system {
  my($self, $id, $command) = @_;
  my($retval);
  
  $self->pushd($id) || return undef;
  $retval = system($command);
  $self->popd();
  return($retval == 0 ? 1 : 0);
}

# Iterate over all of the data directories in the datastore,
# calling a user supplied function for each of them.
# Returns false if there was a problem walking the directory tree, true
# if it was able to walk the dir tree sucessfully.
# The code snippet should return a zero exit status on success, and a
# non-zero exit code of it's choice to indicate a problem.
# If the caller includes a reference to a list for the status information,
# then the code will stuff lists of ("IDENTIFIER", exit_status) pairs into the
# list iff the error status is non-zero.
# The code snippet gets pass a reference to the datastore object and the
# current id.
sub foreach {
  my($self, $code_ref, $status_list_ref) = @_;
  my($cwd);
  my($root);
  my($datadir_depth, $curr_depth);
  my(@dir_list);
  my($dir);
  my($status);

  $cwd = getcwd();

  # trim the trailing "/", for compat w/ while loop code below.
  ($root = $self->root()) =~ s|(.*)/|$1|;
  
  # set up the initial stuff.
  $datadir_depth=$self->depth() + 1;
  $curr_depth = 0;
  push @dir_list, [$curr_depth, $root];

  while(@dir_list) {
    ($curr_depth, $dir) = @{pop(@dir_list)};
    next if (! -d $dir);
    if ($curr_depth == $datadir_depth) {
      CORE::chdir($dir);
      $status = &$code_ref($self, $self->dir_to_id($dir));
      if ($status != 0 && defined($status_list_ref)) {
	push @{$status_list_ref}, [$self->dir_to_id($dir), $status];
      }
      CORE::chdir($cwd);
    }
    else {
      opendir(DIR, $dir) || return(0);
      $curr_depth++;
      push @dir_list,
	map {[$curr_depth, $_]}
	  map {$dir . "/". $_}
	    grep {!/^\./} readdir(DIR);
      closedir(DIR);
    }
  }
  return(1);
}

#
# read a list of identifiers from $ioh, an IO::Handle object
# for each identifier, call the code specified by the reference in $code_ref AND
# if a reference for a list of status values has been provided, stuff the id and that
# code_ref's return status into it.
#

sub ioh_foreach  {
  my($self, $ioh, $code_ref, $status_list_ref) = @_;
  my($cwd);
  my($id);
  my($dir);
  my($status);

  $cwd = getcwd();

  while($id = $ioh->getline()) {
    chomp $id;
    $dir = $self->id_to_dir($id);
    $status = CORE::chdir($dir);
    if (!$status) {		# if the chdir didn't work, status gests
      $status = 2;		# ENOENT (no such file or directory).
    }
    else {
      $status = &$code_ref($self, $self->dir_to_id($dir));
    }
    if ($status != 0 && defined($status_list_ref)) {
      push @{$status_list_ref}, [$id, $status];
    }
    CORE::chdir($cwd);
  }

  return(1);
}

sub iterate {
  my($self, $code_ref, $status_list_ref) = @_;

  return ($self->foreach($code_ref, $status_list_ref));
}

sub rmdir {
  my($self, $id) = @_;
  my($dir);
  my($retval);

  $dir = $self->id_to_dir($id);
  $retval = CORE::system("rm -r $dir > /dev/null 2>&1");

  return($retval == 0 ? 1 : 0);
}

1;

__END__

=head1 NAME

Datastore::Base - A base class for the BDGP's datastore directory
hierarchy tools.

=head1 SYNOPSIS

This class isn't used directly, use a subclass that has the complete
set of bells and whistles.

=head1 DESCRIPTION

Stub documentation for Datastore::Base, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

=head2 METHODS

=head3 C<new>

  $ds = new Datastore::Base("root" => "/tmp/foo", "depth" => 2);
  
  This creates a new datastore object.  It has two optional arguments,
  the name of the root directory used by the datastore and the depth 
  of the datastore subtrees.

=head3 C<root>

  $root = $ds->root("/tmp/the/other/directory");
  $root = $ds->root();

  This set's the name of the root directory used by the datastore
  object.

=head3 C<depth>

  $depth = $ds->depth(1);
  $depth = $ds->depth();

  This set's the depth of the datastore's directory hiearchy.

=head3 C<version>
  
  $version = $ds->version();

  Returns the version number for the entire Datastore library.

=head3 C<mkdir>

  $success = $ds->mkdir("IDENTIFIER");

  This creates a directory in the datastore that corresponds to the
  datastore object $ds.  It runs the system's "mkdir" command w/ the
  "-p" switch via perl's system() function.  Returns a true value if
  the call succeeds and a false value on failure.

  mkdir() respects the datastore's readonly() attribute, and will 
  fail if readonly() is true.

=head3 C<chdir>

  $success = $ds->chdir("IDENTIFIER");

  Creates a datastore directory for "IDENTIFIER" and change the
  process's current working directory to it.  Returns true if the call
  succeeds and false on failure.

=head3 C<cd>

  $success = $ds->cd("IDENTIFIER");

  This is simply an alias for chdir().

=head3 C<pushd>

  $success = $ds->pushd("IDENTIFIER");

  Saves the current working directory in the datastore object and
  chdir() to the directory in the "IDENTIFIER" in the datastore.
  Returns a true value on success and a false value on failure.

=head3 C<popd>

  $success = $ds->popd();

  Change the process's current working directory to the one most
  recently saved by a pushd().  Returns a true value on success and a

=head3 C<readonly>

  $boolean = $ds->readonly();

  Set's the datastore's notion of whether it is read/write or readonly.
  The default is currently read/write ($ds->readonly defaults to 0).

=head3 C<system>

  $success = $ds->system("IDENTIFIER", "command")

  Uses $ds->pushd to change the process's current working directory to
  the location for "IDENTIFIER" in the datastore, runs "command" using
  perl's system() command, and uses $ds->popd change back to the
  original directory.  Returns a true value on success and a false
  value on failure.

=head3 C<foreach>

  Iterate over all of the data directories in the datastore, calling a
  user supplied reference to a function for each of them.  Returns
  false if there was a problem walking the directory tree, true if it
  was able to walk the dir tree sucessfully.  The user's function is
  passed two arguments, a reference to the datastore being traversed
  and the identifier that corresponds to the current directory.  The
  function should return a zero exit status on success, and a non-zero
  exit code of it's choice to indicate a problem.  If the caller
  includes a reference to a list for the status information, then the
  code will stuff lists of ["IDENTIFIER", exit_status] pairs into the
  list iff the error status is non-zero.  The code snippet gets pass a
  reference to the datastore object and the current id.

  Here's a little sample program that makes some directories in a
  datastore and then iterates over them calling a little function at
    a) prints out some information about where it is 
    b) exits with a successful code unless the ID ends in "9"
    c) prints out the list of errors that were found.

  #!/usr/bin/env perl
  
  use Datastore::MD5;
  use Cwd;
  
  $ds = new Datastore::MD5("root" => "/tmp/test-me", "depth" => 2);
  
  $ds->mkdir("CG1239");
  $ds->mkdir("CG1234");
  $ds->mkdir("CG0670");
  $ds->mkdir("CG0669");
  $ds->mkdir("CG0668");
  $ds->mkdir("CG0667");
  $ds->mkdir("CG0666");
  
  @status = ();	
  $success = $ds->foreach(\&doit_toit, \@status);
  
  print "Errors occurred (ID -- status)\n";
  print join "\n", map {"  " . join " -- ",  @{$_}} @status;
  print "\n";
  
  # silly little demo sub that returns an error for any id ending in 9
  sub doit_toit {
    my($datastore) = shift;
    my($id) = shift;
  
    print "Working on $id: ", getcwd(), "\n";
    if ($id =~ m|.*9|) {
      return 1;
    }
    else {
      return 0;
    }
  }
    
=head3 C<ioh_foreach>

  Read a list of identifiers from $ioh, an IO::Handle object, and for
  each identifier, call the code specified by the reference in
  $code_ref, and if a reference for a list of status values has been
  provided, stuff the id and that code_ref's return status into it.

  The return status value "2" is reserved for ioh_foreach()'s use, it
  corresponds to the errno error "ENOENT" which is used when there is
  no such file or directory.  ioh_foreach() will save this status
  value in the status array when it is unable to chdir to the
  directory that corresponds to an id.

  use Datastore::MD5;
  use Cwd;
  use IO::File;

  $ds = new Datastore::MD5("root" => "/tmp/test-me", "depth" => 2);

  $fh = new IO::File;
  $fh->open("< CG_LIST") || die "Unable to open CG_LIST";

  $ds->ioh_foreach($fh, \&doit_toit, \@status);
    
  $fh->close() || die "Unable to close CG_LIST";

  print "Errors occured for (ID -- status)\n";
  print join "\n", map {"  " . join " -- ",  @{$_}} @status;
  print "\n";
  
  # silly little demo sub that returns an error for any id ending in 9
  sub doit_toit {
    my($id) = shift;
  
    print "Working on $id: ", getcwd(), "\n";
    if ($id =~ m|.*9|) {
      return 1;
    }
    else {
      return 0;
    }
  }

  Similar results can be achieved for identifier lists that are stored
  in various "in-core" datastructures using IO handles created with
  the perl objects described in the IO::Stringy man pages.

=head3 C<iterate>

  This is the deprecated name for C<foreach>.  It's only supported for
  backward compatibility.

=head3 C<rmdir>

  $success = $ds->rmdir("IDENTIFIER");

  This removes a directory in the datastore that corresponds to the
  datastore object $ds.  It runs the system's "rm" command w/ the "-r"
  switch via perl's system() function.  Returns a true value if the call
  succeeds and a false value on failure.

=head1 AUTHOR

George Hartzell, E<lt>hartzell @ fruitfly.orgE<gt>

=head1 SEE ALSO

L<Datastore::MD5>, L<Datastore::CG>, L<perl>.

=cut
