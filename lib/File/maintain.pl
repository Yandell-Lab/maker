#! /usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use File::NFSLock;
use Storable;
use vars qw($LOCK);

BEGIN {
    $SIG{QUIT} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{INT} = sub{$LOCK->unlock if($LOCK); exit(0)};
}

my $file = shift;
my $time = shift;
my $pid = shift;

die "ERROR: Lock serialization file does not exist\n\n"
    if(! defined($file) || ! -f $file);
die "ERROR:  Lacking input for lock maintainer\n\n"
    if(! defined($time) || ! defined($pid));

my $LOCK = Storable::retrieve($file);
unlink($file);

#attempted fix of file not being deleted error
for(my $i = 0; $i < 20 && -f $file; $i++){
    unlink($file);
    sleep 1 if(-f $file);
}

die "ERROR: Could not retrieve lock" if(! $LOCK);

while(-f $LOCK->{lock_file}){
    $LOCK->refresh;
    sleep $time;
    my $p = `ps -p $pid`;
    last if(@{[$p =~ /(\n)/g]} < 2); #pid not running
}

$LOCK->unlock;
