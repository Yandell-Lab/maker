#!/usr/bin/env perl
use warnings;
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
    if(! Proc::Signal->exists_proc_by_id($pid)){
	$LOCK->unlock if($LOCK);
	exit(0);
    }

    $LOCK->refresh;
    sleep $time;
}

$LOCK->unlock;

exit(0);

DESTROY {
    $LOCK->unlock if($LOCK);
}
