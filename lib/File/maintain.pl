#!/usr/bin/env perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../../perl/lib";
use File::NFSLock;
use Storable;
use vars qw($LOCK);
use Proc::Signal;
use URI::Escape;
use Perl::Unsafe::Signals;

BEGIN {
    $SIG{QUIT} = sub{exit()};
    $SIG{ABRT} = sub{exit()};
    $SIG{KILL} = sub{exit()};
    $SIG{TERM} = sub{exit()};
    $SIG{STOP} = sub{exit()};
    $SIG{INT}  = sub{exit()};
    $SIG{USR1} = sub{exit()};

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

select((select(STDOUT), $|=1)[0]);

die "ERROR: Lock serialization does not exist\n\n"
    if(! defined($serial));
die "ERROR:  Lacking input for lock maintainer\n\n"
    if(! defined($time) || ! defined($pid));

$serial = uri_unescape($serial);
$LOCK = Storable::thaw($serial);

die "ERROR: Could not retrieve lock" if(! $LOCK);

$LOCK->_unlink_block(1);
my $step = time();
while(-f $LOCK->{lock_file}){
    if(getppid != $pid){
	$LOCK->unlock(1) if($LOCK);
	exit(0);
    }

    #these checks are actually kind of expensive (so don't do them often)
    if(abs(time() - $step) > 60){
	if(! Proc::Signal::exists_proc_by_id($pid)){
	    $LOCK->unlock(1) if($LOCK);
	    exit(0);
	}
	if(my $p = Proc::Signal::get_proc_by_id($pid)){
	    if($p->state eq 'zombie'){
		$LOCK->unlock(1) if($LOCK);
		exit(0);
	    }
	}
	$step = time();
    }

    exit (1) unless($LOCK->refresh);

    sleep $time;
}

exit(1);
