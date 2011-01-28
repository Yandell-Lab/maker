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
use Perl::Unsafe::Signals;

BEGIN {
    $SIG{QUIT} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{KILL} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{TERM} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{STOP} = sub{$LOCK->unlock if($LOCK); exit(0)};
    $SIG{INT}  = sub{$LOCK->unlock if($LOCK); exit(0)};

    $SIG{'__WARN__'} = sub {
	warn $_[0] if ( $_[0] !~ /Not a CODE reference/ &&
			$_[0] !~ /Can\'t store item CODE/ &&
			$_[0] !~ /Find\:\:skip_pattern|File\/Find\.pm/ &&
			$_[0] !~ /Ran into unknown state/
			);
    };
}

my $pid = shift;
my $time = shift;
my $serial = shift;

UNSAFE_SIGNALS {
    select((select(STDOUT), $|=1)[0]);
    
    die "ERROR: Lock serialization does not exist\n\n"
	if(! defined($serial));
    die "ERROR:  Lacking input for lock maintainer\n\n"
	if(! defined($time) || ! defined($pid));
    
    $serial = uri_unescape($serial);
    $LOCK = Storable::thaw($serial);
    
    die "ERROR: Could not retrieve lock" if(! $LOCK);
    
    while(-f $LOCK->{lock_file}){
	if(! Proc::Signal::exists_proc_by_id($pid)){
	    $LOCK->unlock if($LOCK);
	    exit(0);
	}
	if($LOCK && ! $LOCK->still_mine){
	    exit(1);
	}
	exit (1) unless($LOCK->refresh);
	sleep $time;
    }
    
    $LOCK->unlock if($LOCK);

    exit(0);
};
    
DESTROY {
    $LOCK->unlock if($LOCK);
}
    
