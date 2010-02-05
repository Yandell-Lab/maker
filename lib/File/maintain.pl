#! /usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use File::NFSLock;
use Storable;
use vars qw($LOCK);
use Proc::Signal;

BEGIN {
    $SIG{QUIT} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{INT} = sub{$LOCK->unlock if($LOCK); exit(0)};
}

my $file = shift;
my $time = shift;
my $pid = shift;

select((select(STDOUT), $|=1)[0]);

die "ERROR: Lock serialization file does not exist\n\n"
    if(! defined($file) || ! -f $file);
die "ERROR:  Lacking input for lock maintainer\n\n"
    if(! defined($time) || ! defined($pid));

$LOCK = Storable::retrieve($file);
unlink($file);

#attempted fix of file not being deleted error
for(my $i = 0; $i < 20 && -f $file; $i++){
    sleep 1;
    unlink($file);
}

die "ERROR: Could not retrieve lock" if(! $LOCK);

my $ok = (-f $file) ? 0 : 1;
print $ok."\n";

while(-f $LOCK->{lock_file}){
    $LOCK->refresh;
    sleep $time;
    my $name = Proc::Signal::get_name_by_id($pid);
    last if(! $name || $name !~ /^maker$|^mpi_maker$/);
}

$LOCK->unlock;

