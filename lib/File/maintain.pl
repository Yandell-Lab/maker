#! /usr/bin/env perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/../../perl/lib";
use File::NFSLock;
use Storable;
use vars qw($LOCK);
use Proc::Signal;
use URI::Escape;

BEGIN {
    foreach my $sig (keys %SIG){
	$SIG{$sig} = sub{$LOCK->unlock if($LOCK); exit(0)};
    }

    $SIG{QUIT} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{KILL} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{TERM} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{STOP} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{INT}  = sub{$LOCK->unlock if($LOCK); exit(0)};
}

my $pid = shift;
my $time = shift;
my $serial = shift;

select((select(STDOUT), $|=1)[0]);

die "ERROR: Lock serialization does not exist\n\n"
    if(! defined($serial));
die "ERROR:  Lacking input for lock maintainer\n\n"
    if(! defined($time) || ! defined($pid));

$serial = uri_unescape($serial);
$LOCK = Storable::thaw($serial);

die "ERROR: Could not retrieve lock" if(! $LOCK);

while(-f $LOCK->{lock_file}){
    $LOCK->refresh;
    sleep $time;
    my $name = Proc::Signal::get_name_by_id($pid);
    last if(! $name || $name !~ /^maker$|^mpi_maker$/);
}

$LOCK->unlock;

exit(0);
