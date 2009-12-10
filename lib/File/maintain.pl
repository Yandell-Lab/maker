#! /usr/bin/perl -w
BEGIN {
    $SIG{QUIT} = sub{exit(0)};
    $SIG{INT} = sub{exit(0)};
}

use strict;
use FindBin;
use lib "$FindBin::Bin/../";
use File::NFSLock;
use Storable;

my $file = shift;
my $time = shift;
my $pid = shift;

die "ERROR: Lock serialization file does not exist\n\n"
    if(! defined($file) || ! -f $file);
die "ERROR:  Lacking input for lock maintainer\n\n"
    if(! defined($time) || ! defined($pid));

my $lock = Storable::retrieve($file);
unlink($file);

die "ERROR: Could not retrieve lock" if(! $lock);

while(-f $lock->{lock_file}){
    $lock->refresh;
    sleep $time;
    my $p = `ps -p $pid`;
    last if(@{[$p =~ /(\n)/g]} < 2); #pid not running
}

$lock->unlock;
