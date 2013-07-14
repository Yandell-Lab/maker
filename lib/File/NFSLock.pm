# -*- perl -*-#
# This module is a rewrite of the File::NFSLock utility from CPAN
# that was originally developed by Paul T Seamons and Rob B Brown 
################################################################

package File::NFSLock;

use strict;
use POSIX;
use POSIX qw(:sys_wait_h);
use File::Copy;
use Storable;
use IPC::Open3;
use Proc::Signal;
use URI::Escape;
use Sys::Hostname;
use Digest::MD5;
use Fcntl;
use vars qw(@ISA @EXPORT $VERSION %TYPES %LOCK_LIST $TYPE
            $LOCK_EXTENSION $HOSTNAME @CATCH_SIGS);

require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(&LOCK_SH &LOCK_EX &LOCK_NB);
$VERSION = '1.20';
$TYPE = 'NFS';

#Constants
sub LOCK_SH {Fcntl::LOCK_SH}
sub LOCK_EX {Fcntl::LOCK_EX}
sub LOCK_NB {Fcntl::LOCK_NB}

$LOCK_EXTENSION = '.NFSLock'; # customizable extension
$HOSTNAME = &Sys::Hostname::hostname();

### Number converion for lock_type
%TYPES = (LOCK_EX     => LOCK_EX,
	  BLOCKING    => LOCK_EX,
	  BL          => LOCK_EX,
	  EXCLUSIVE   => LOCK_EX,
	  EX          => LOCK_EX,
	  LOCK_NB     => LOCK_NB,
	  NONBLOCKING => LOCK_NB,
	  NB          => LOCK_NB,
	  LOCK_SH     => LOCK_SH,
	  SHARED      => LOCK_SH,
	  SH          => LOCK_SH,);

#signals to catch and replace
@CATCH_SIGS = qw(TERM INT QUIT ABRT KILL TERM);
my $graceful_sig = sub{ exit; };#ensures that DESTROY and END are called

#-------------------------------------------------------------------------------
#----------------------------------- METHODS -----------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $self = {};
    my $class = shift;

    bless($self, $class);

    $self->{file}               = shift || '';
    $self->{lock_type}          = shift || '';
    $self->{blocking_timeout}   = shift || 0;
    $self->{stale_lock_timeout} = shift || 0;
    $self->{lock_pid} = $$;
    $self->{hostname} = $HOSTNAME;
    $self->{unlocked} = 1;
    ($self->{id} = rand()) =~ s/\.//;
    $self->{quit_time} = ($self->{blocking_timeout}) ? time+$self->{blocking_timeout} : 0;

    if(ref($self->{file}) eq 'HASH'){
	my $h = $self->{file};
	$self->{file}               = $h->{file} || '';
        $self->{lock_type}          = $h->{lock_type} || '';
        $self->{stale_lock_timeout} = $h->{stale_lock_timeout} || 0;

	if($h->{blocking_timeout}){
	    $self->{blocking_timeout} = $h->{blocking_timeout};
	    $self->{quit_time} = time + $self->{blocking_timeout};
	}
	elsif($h->{quit_time}){
	    $self->{quit_time} = $h->{quit_time};
	    $self->{blocking_timeout} = $self->{quit_time} - time || 1;
	    return if($self->{blocking_timeout} < 0); #already finished
	}
	else{
	    $self->{blocking_timeout} = 0;
	    $self->{quit_time} = 0;
	}
    }

    $self->{lock_type} = $TYPES{$self->{lock_type}} if($self->{lock_type} !~ /^\d+/);

    die "ERROR: You must specify a file to lock.\n"unless(length($self->{file}));

    ### choose the lock filename
    my $lock_file = $self->{file}.$LOCK_EXTENSION;
    my ($name) = $lock_file =~ /([^\/]+)$/;
    if($name !~ /^\.NFSLock/){
	$lock_file =~ s/([^\/]+)$/\.NFSLock\.$1/;
    }
    $self->{lock_file} = $lock_file;
    $self->{rand_file} = $self->_rand_name($self->{file});

    #replace signals
    foreach my $signal (@CATCH_SIGS) {
	next if(!$TYPE); #locking is off
	$SIG{$signal} = $graceful_sig if(!$SIG{$signal} || $SIG{$signal} eq 'DEFAULT');
    }

    #make lock
    if($self->{lock_type} == LOCK_SH){
	$self->_new_SH;
    }
    elsif($self->{lock_type} == LOCK_EX){
	$self->_new_EX;
    }
    elsif($self->{lock_type} == LOCK_NB){
	delete($self->{blocking_timeout});
	delete($self->{quit_time});
	$self->_new_NB;
    }
    else{
	die "ERROR: Invalid lock type\n";
    }

    if($self->{unlocked}){
	# Revert handler
	foreach my $signal (@CATCH_SIGS) {
	    delete $SIG{$signal} if($SIG{$signal} && $SIG{$signal} eq $graceful_sig);
	}

	return undef;
    }

    return $self;
}

#-------------------------------------------------------------------------------
sub _new_EX {
    my $self       = shift;
    my $lock_file  = $self->{lock_file};
    my $rand_file  = $self->{rand_file};
    my $pid        = $self->{lock_pid};
    my $hostname   = $self->{hostname};
    my $quit_time  = $self->{quit_time};
    my $stale_lock_timeout = $self->{stale_lock_timeout};
 
    my $stack;
    my $gfail = 0; #counter to trigger ghost file check
    while (1) {
	#if locking has been turned off then just pretend success
	if(!$TYPE){
	    $self->_create_magic($lock_file);
	    last;
	}

	if($quit_time && time > $quit_time){
	    unlink($rand_file);
	    return undef;
	}

	if(-f $lock_file){
	    #create lock stack in IO
	    if($LOCK_LIST{"$lock_file\_EX"} && $self->still_mine){
		last if($self->_create_magic($lock_file));
		next;
	    }
	    elsif(! $stack){
		local $LOCK_EXTENSION = '.STACK';
		$stack = new File::NFSLock({file      => _fix_name($lock_file),
					    lock_type => LOCK_EX,
					    quit_time => $quit_time,
					    stale_lock_timeout => $stale_lock_timeout});
		next if(! $stack);
	    }
	    
	    ### remove an old lockfile if it is older than the stale_timeout
	    if($stale_lock_timeout && (stat($lock_file))[9] + $stale_lock_timeout < time){
		if(!$stack->still_mine){
		    undef $stack;
		    next;
		}
		
		#get the lock before releasing the stack
		last if($self->_create_magic($lock_file));
	    }

	    ### If lock exists and is readable, see who is mooching on the lock
	    if(open (my $FH,"+<$lock_file")){
		my @mine;
		my @them;
		my @dead;
		local $/ = "\n"; #just incase the user changed this
		while(defined(my $line=<$FH>)){
		    next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
		    
		    my $l_host = $1;
		    my $l_pid  = $2;
		    my $l_time = $3;
		    my $l_id   = $4;   
		    my $l_kind = $5;
		    
		    if($l_host eq $hostname){
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			elsif ($l_pid == $pid) { # This is me.
			    push(@mine, $line);
			    last;
			}
			elsif(kill(0, $l_pid)) { # Still running on this host.
			    push(@them, $line);
			    last;
			}
			else{ # Finished running on this host.
			    push(@dead, $line);
			    next;
			}
		    }
		    else{
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			else{
			    push(@them, $line);
			    last;
			}
		    }
		}
		close($FH);
		
		### If there was a stale lock discovered...
		if(!@them && !@mine){
		    last if($self->_create_magic($lock_file));
		}
		
		### If attempting to acquire the same type of lock
		### that it is already locked with, and I've already
		### locked it myself, then it is safe to lock again.
		### Just kick out successfully without really locking.
		### Assumes locks will be released in the reverse
		### order from how they were established.
		if (@mine && !@them && !$self->_check_shared){
		    last if($self->_create_magic($lock_file));
		}
	    }
	}
	else{
	    last if($self->_create_magic($rand_file) && $self->_do_lock);
	    $gfail++; #hmm, either race condidtion or ghost file

	    #check for ghost file error (corrupt filesystem)
	    die "FATAL: Potentially corrupt filesystem. The file $lock_file\n".
                "is listed within the directory, but cannot lstat. You may need\n".
                "to run fsck to repair it.\n"
		if($gfail % 10 == 0 && $self->ghost_check($lock_file));
	}

	### wait a moment or timeout
	my $time = time();
	my $dob = (stat($lock_file))[9];
	next if(!defined $dob); #file no longer exists

	my $age = $time - $dob;
	my $sleep_max = $age/100; #1% of current lock time
	$sleep_max = 1 if($sleep_max < 1);
	$sleep_max = 10 if($sleep_max > 10);
	$sleep_max = 1 if($quit_time && $quit_time - $time < $sleep_max);
	$sleep_max = 1 if($stale_lock_timeout && $stale_lock_timeout - $age < $sleep_max);
	
	if($quit_time && $time > $quit_time){
	    unlink($rand_file);
	    return undef;
	}
	else{
	    sleep $sleep_max;
	    next;
	}
    }

    $self->uncache; #trick for NFS to update cache   
    
    if($self->still_mine){
	unlink($rand_file) if($TYPE);
	$LOCK_LIST{"$lock_file\_EX"}++;
	delete($self->{unlocked});
	return $self;
    }
    else{
	return $self->_new_EX();
    }
}

#-------------------------------------------------------------------------------
sub _new_NB {
    my $self       = shift;
    my $lock_file  = $self->{lock_file};
    my $rand_file  = $self->{rand_file};
    my $pid        = $self->{lock_pid};
    my $hostname   = $self->{hostname};
    my $stale_lock_timeout = $self->{stale_lock_timeout};

    my $stack;
    { #single iteration block (so last will work))
	#if locking has been turned off then just pretend success
	if(!$TYPE){
	    $self->_create_magic($lock_file);
	    last;
	}
	if(-f $lock_file){
	    if($LOCK_LIST{"$lock_file\_EX"} && $self->still_mine){
		last if($self->_create_magic($lock_file));		
	    }
	    else{
		#create lock stack in IO
		local $LOCK_EXTENSION = '.STACK';
		$stack = new File::NFSLock(_fix_name($lock_file), LOCK_NB, 0, $stale_lock_timeout);
		if(! $stack){
		    unlink($rand_file);
		    return undef;
		}
	    }
	    
	    ### remove an old lockfile if it is older than the stale_timeout
	    if($stale_lock_timeout && (stat($lock_file))[9] + $stale_lock_timeout < time){
		if(!$stack->still_mine){
		    unlink($rand_file);
		    return undef;
		}
		
		#get the lock before releasing the stack
		last if($self->_create_magic($lock_file));
	    }

	    if(open (my $FH,"+<$lock_file")){
		my @mine;
		my @them;
		my @dead;		
		local $/ = "\n"; #just incase the user changed this		
		while(defined(my $line=<$FH>)){
		    next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
		    
		    my $l_host = $1;
		    my $l_pid  = $2;
		    my $l_time = $3;
		    my $l_id   = $4;   
		    my $l_kind = $5;
		    
		    if($l_host eq $hostname){
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			elsif ($l_pid == $pid) { # This is me.
			    push(@mine, $line);
			    last;
			}
			elsif(kill(0, $l_pid)) { # Still running on this host.
			    push(@them, $line);
			    last;
			}
			else{ # Finished running on this host.
			    push(@dead, $line);
			    next;
			}
		    }
		    else{
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			else{ # assume it is still running.
			    push(@them, $line);
			    last;
			}
		    }
		}
		close($FH);
		
		### If there was a stale lock discovered...
		last if(!@them && !@mine && $self->_create_magic($lock_file));
		
		### If attempting to acquire the same type of lock
		### that it is already locked with, and I've already
		### locked it myself, then it is safe to lock again.
		### Just kick out successfully without really locking.
		### Assumes locks will be released in the reverse
		### order from how they were established.
		last if (@mine && !@them && !$self->_check_shared);
	    }
	    
	    unlink($rand_file);
	    return undef;
	}
	else{
	    return undef unless($self->_create_magic($rand_file) && $self->_do_lock);
	}
    }
       
    $self->uncache; #trick for NFS to update cache
 
    if($self->still_mine){
	unlink($rand_file) if($TYPE);
	$LOCK_LIST{"$lock_file\_EX"}++;
	delete($self->{unlocked});
	return $self;
    }
    else{
	return $self->_new_NB();
    }
}
#-------------------------------------------------------------------------------
sub _new_SH {
    my $self       = shift;
    my $lock_file  = $self->{lock_file};
    my $rand_file  = $self->{rand_file};
    my $pid        = $self->{lock_pid};
    my $hostname   = $self->{hostname};
    my $quit_time  = $self->{quit_time};
    my $stale_lock_timeout = $self->{stale_lock_timeout};
 
    my $stack;
    my $gfail = 0; #counter to trigger ghost file check
    while (1) {
	#if locking has been turned off then just pretend success
	if(!$TYPE){
	    $self->_create_magic($lock_file);
	    last;
	}

	if($quit_time && time > $quit_time){
	    unlink($rand_file);
	    return undef
	}

	my $exists;
	my $shared;
	if(-f $lock_file){
	    $exists = 1;
	    $shared = $self->_check_shared;
	}

	if($exists && $shared){ #this is a shared lock so just join in
	    last if($self->_do_lock_shared);
	}
	elsif($exists && !$shared){ #exclusive lock, see if I can take it
	    #create lock stack in IO
	    local $LOCK_EXTENSION = '.STACK';
	    if(!$stack){
		$stack = new File::NFSLock({file      => _fix_name($lock_file),
					    lock_type => LOCK_EX,
					    quit_time => $quit_time,
					    stale_lock_timeout => $stale_lock_timeout});
		
		next if(! $stack);
	    }

	    ### remove an old lockfile if it is older than the stale_timeout
	    if($stale_lock_timeout && (stat($lock_file))[9] + $stale_lock_timeout < time){
		if(!$stack->still_mine){
		    unlink($rand_file);
		    return undef;
		}
		
		#get the lock before releasing the stack
		last if($self->_create_magic($lock_file));
	    }

	    if(open (my $FH,"+<$lock_file")){
		my @mine;
		my @them;
		my @dead;		
		local $/ = "\n"; #just incase the user changed this
		while(defined(my $line=<$FH>)){
		    next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
		    
		    my $l_host = $1;
		    my $l_pid  = $2;
		    my $l_time = $3;
		    my $l_id   = $4;   
		    my $l_kind = $5;
		    
		    if($l_host eq $hostname){
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			elsif ($l_pid == $pid) { # This is me.
			    push(@mine, $line);
			    last;
			}
			elsif(kill(0, $l_pid)) { # Still running on this host.
			    push(@them, $line);
			    last;
			}
			else{ # Finished running on this host.
			    push(@dead, $line);
			    next;
			}
		    }
		    else{
			# see if I can break it
			if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
			    push(@dead, $line);
			    next;
			}
			else{ # assume it is still running.
			    push(@them, $line);
			    last;
			}
		    }
		}
		close($FH);
		
		### If there was a stale lock discovered...
		last if(!@them && !@mine && $self->_create_magic($lock_file));
	    }
	}
	else{
	    last if($self->_create_magic($rand_file) && $self->_do_lock_shared);
	    $gfail++; #hmm, either race condidtion or ghost file

	    die "FATAL: Potentially corrupt filesystem. The file $lock_file\n".
		"is listed within the directory, but cannot lstat. You may need\n".
		"to run fsck to repair it.\n"
		if($gfail % 10 == 0 && $self->ghost_check($lock_file));

	    #semi-recoverable if rand_file is ghost
	    if($gfail % 10 == 0 && $self->ghost_check($rand_file)){
		warn "WARNING: Potentially corrupt filesystem. The file $rand_file\n".
		     "is listed within the directory, but cannot lstat. You may need\n".
		     "to run fsck to repair it.\n";
		$self->{rand_file} = $self->_rand_name($lock_file);
	    }
	}

	### wait a moment or timeout
	my $time = time();
	my $dob = (stat($lock_file))[9];
	next if(!defined $dob); #file no longer exists

	my $age = $time - $dob;
	my $sleep_max = $age/100; #1% of current lock time
	$sleep_max = 1 if($sleep_max < 1);
	$sleep_max = 10 if($sleep_max > 10);
	$sleep_max = 1 if($quit_time && $quit_time - $time < $sleep_max);
	$sleep_max = 1 if($stale_lock_timeout && $stale_lock_timeout - $age < $sleep_max);
	
	if($quit_time && $time > $quit_time){
	    unlink($rand_file);
	    return undef;
	}
	else{
	    sleep $sleep_max;
	    next;
	}
    }

    $self->uncache; #trick for NFS to update cache   
    
    if($self->still_mine){
	unlink($rand_file) if($TYPE);
	$LOCK_LIST{"$lock_file\_SH"}++;
	delete($self->{unlocked});
	return $self;
    }
    else{
	return $self->_new_SH();
    }
}

#-------------------------------------------------------------------------------
#these happen frequently on corrupt NFS filesystems
#can be called as function as well
sub ghost_check {
    my $self = shift;
    my $loc = shift;

    #for function call
    $loc = $self if($self && ! ref($self));
    return 0 if(! defined($loc));

    #check for ghost file (they can cause infinite loop)
    $loc =~ s/\/+$//;
    my ($dir, $name) = $loc =~ /(.*\/)([^\/]*)$/;
    $dir = './' if(! defined($dir) && defined($name));
    return 0 if(! defined($dir));

    opendir(DIR, $dir);
    my $ghost = grep {$_ eq $name && ! -e $loc} readdir(DIR);
    close(DIR);

    #I seem to have found one (make extra sure)
    if($ghost){
	sleep 5; #ghosts are persistent so long wait is sufficient
	opendir(DIR, $dir);
	$ghost = grep {$_ eq $name && ! -e $loc} readdir(DIR);
	close(DIR);
    }
    close(DIR);
    
    return $ghost;
}

#-------------------------------------------------------------------------------
sub unlock {
    my $self = shift;
    my $force = shift;

    #pretend success when locking is off
    if(!$TYPE){
	($self->{lock_type} == LOCK_SH) ? $self->_do_unlock_shared : $self->_do_unlock;
        $self->{unlocked} = 1;
        return 1;
    }

    $self->_unlink_block(0) if($force);

    local $?; #localize to clear killed locks
    my $stat;
    if($self->{unlocked}){
	#clean up any children that are floating around
	do {
	    $stat = waitpid(-1, WNOHANG);
	} while $stat > 0;    
	return 1;
    }

    #remove any temporary files
    unlink($self->{rand_file});

    #remove maintainer if running
    my $m_pid = $self->{_maintain};
    if(defined($m_pid) && Proc::Signal::id_matches_pattern($m_pid, 'maintain\.pl|defunct')){
	close($self->{_IN}) if(ref $self->{_IN} eq 'GLOB');

	#signal maintainer
	kill(SIGUSR1, $self->{_maintain});
	$stat = waitpid($self->{_maintain}, WNOHANG);

	#attempt kill multiple times if still running
	my $count = 0;
	while($stat == 0 && $count < 200){
	    kill(SIGUSR1, $self->{_maintain}); #try multiple signal ending in signal 9
	    usleep(0.1) if($stat == 0);
	    $stat = waitpid($self->{_maintain}, WNOHANG);
	    $count++;
	}
	
	#if still running, do this
	if($stat == 0){
	    kill(SIGKILL, $self->{_maintain});
	    waitpid($self->{_maintain}, 0);
	}
	
	$self->{_maintain} = undef;
    }

    #clean up any children that are floating around
    do {
	$stat = waitpid(-1, WNOHANG);
    } while $stat > 0;    

    my $ok = ($self->{lock_type} == LOCK_SH) ? $self->_do_unlock_shared : $self->_do_unlock;
    $self->{unlocked} = 1;
    
    # Revert handler once all locks are finished
    if(! grep {$_ > 0} values %LOCK_LIST){
	foreach my $signal (@CATCH_SIGS) {
	    delete $SIG{$signal} if($SIG{$signal} && $SIG{$signal} eq $graceful_sig);
	}
    }
    
    return $ok;
}
#-------------------------------------------------------------------------------
sub _check_shared {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    return ((stat($lock_file))[2] & 1); 
}

#-------------------------------------------------------------------------------
sub _rand_name {
    my $self = shift;
    my $file = shift;
    
    my $rand;
    while(!$rand || -f $rand){
	$rand = "$file.tmp.".(time % 10000).".$$.".(rand()*10000);
	$rand =~ s/([^\/]+)$/\.NFSLock\.$1/;
    }

    return $rand;
}
#-------------------------------------------------------------------------------
sub _create_magic {
  my $self = shift;
  my $append_file = shift;
  my $shared = ($self->{lock_type} == LOCK_SH) ? 1 : 0;
  my $pid = $self->{lock_pid};
  my $hostname   = $self->{hostname};
  my $id = $self->{id};

  $self->{lock_line} = "$hostname $pid ".time()." $id $shared\n";

  #pretend success when locking is off
  if(!$TYPE){
      if($shared){
	  $self->{id_line} ||= Digest::MD5::md5_hex($self->{lock_line})."\n";
      }
      return 1;
  }

  my $first  = ($append_file eq $self->{rand_file});

  if($shared){
      open(my $FH,"+>>$append_file") or return undef;
      local $/ = "\n"; #just incase the user changed this

      #reading
      if(! $first){
	  seek $FH, 0, 0;
	  my $line = <$FH>; #should always be first line
	  $self->{id_line} = $line if($line !~ / /);
      }

      $self->{id_line} ||= Digest::MD5::md5_hex($self->{lock_line})."\n";

      #writing
      print $FH $self->{id_line} if($first);
      print $FH $self->{lock_line};
      close($FH);

      chmod((0600 | 1), $append_file) or die "ERROR: I can't chmod files for locking";
  }
  else{
      open(my $FH,">$append_file") or return undef;
      print $FH $self->{lock_line};
      close($FH);
      chmod(0600, $append_file)
  }

  return 1;
}
#-------------------------------------------------------------------------------
sub _do_lock {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    my $rand_file = $self->{rand_file};

    #pretend success when locking is off
    return 1 if(!$TYPE);
    return if(! -f $rand_file);

    my $success;
    if($TYPE eq 'NFS'){
	$success = link( $rand_file, $lock_file ) && (stat($lock_file))[3] == 2;
	unlink($rand_file);
    }
    elsif($TYPE eq 'POSIX'){
	open(my $fh, ">>$lock_file");
	my $locked = flock($fh, LOCK_EX);
	$success = File::Copy::move($rand_file, $lock_file) if($locked);
	unlink($rand_file);
	close($fh);
    }

    return $success;
}
#-------------------------------------------------------------------------------
sub _do_lock_shared {
    my $self = shift;
    my $lock_file  = $self->{lock_file};
    my $rand_file  = $self->{rand_file};

    #pretend success when locking is off
    return 1 if(!$TYPE);

    ### Try to create $lock_file from the special
    ### file with the magic SHAREBIT set.
    my $success;
    if($TYPE eq 'NFS' && -f $rand_file){
	$success = link($rand_file, $lock_file) && (stat($lock_file))[3] == 2;
	unlink($rand_file);
    }
    elsif($TYPE eq 'POSIX' && -f $rand_file){
	open(my $fh, ">>$lock_file");
	my $locked = flock($fh, LOCK_EX);
	$success = File::Copy::move($rand_file, $lock_file) if($locked);
	unlink($rand_file);
	close($fh);
    }

    ### lock the locking process
    local $LOCK_EXTENSION = ".share";
    my $lock = new File::NFSLock(_fix_name($lock_file), LOCK_EX, 62, 60);    

    if(!$success && -f $lock_file && !$self->_check_shared){
	#lock exists as exclusive
	return undef;
    }
    elsif(!$success) {
	### it's already shared or is yet to exist so append my lock
	$self->_create_magic($lock_file);
    }
    elsif($success){
	return 1;
    }
}
#-------------------------------------------------------------------------------
sub _do_unlock {
  my $self = shift;

  return 1 if($self->{unlocked});

  $LOCK_LIST{$self->{lock_file}."_EX"}--;
  $self->{unlocked} = 1;

  #pretend success when locking is off
  return 1 if(!$TYPE);

  return 1 if(! -f $self->{lock_file});
  return 1 if($self->_check_shared);
  return 1 unless($LOCK_LIST{$self->{lock_file}."_EX"} <= 0);

  $LOCK_LIST{$self->{lock_file}."_EX"} = 0; # just in case it's less than 0
  my $stat = unlink($self->{lock_file}) if(!$self->_unlink_block && $self->still_mine);

  return $stat;
}
#-------------------------------------------------------------------------------
sub _do_unlock_shared {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    my $lock_line = $self->{lock_line};
    my $stale_lock_timeout = $self->{stale_lock_timeout};
    my $hostname   = $self->{hostname};
    my $pid = $self->{lock_pid};
    my $id = $self->{id};    

    return 1 if($self->{unlocked});

    $LOCK_LIST{$self->{lock_file}."_SH"}--;
    $self->{unlocked} = 1;

    #pretend success when locking is off
    return 1 if(!$TYPE);

    return 1 if(! -f $self->{lock_file});
    return 1 if(! $self->_check_shared);

    ### lock the unlocking process
    local $LOCK_EXTENSION = '.share';
    my $lock = new File::NFSLock (_fix_name($lock_file), LOCK_EX, 62, 60);

    ### get the handle on the lock file
    my $FH;
    if( ! open ($FH,"+<$lock_file") ){
	if( ! -f $lock_file || $! =~ /No such file or directory/){
	    return 1;
	}else{
	    die "Could not open for writing shared lock file $lock_file ($!)";
	}
    }
    
    ### read existing file
    my $id_line = '';
    my $content = '';
    my $time = time();
    local $/ = "\n"; #just incase the user changed this    
    while(defined(my $line=<$FH>)){
	if($line !~ / /){
	    $id_line .= $line;
	    next;
	}

	next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
	my $l_host = $1;
	my $l_pid  = $2;
	my $l_time = $3;
	my $l_id   = $4;
	my $l_kind = $5;

	if($l_host eq $hostname){
	    if($l_pid == $pid){
		next if($LOCK_LIST{$lock_file."_SH"} <= 0);
		next if($l_id == $id);
	    }
	    elsif(!kill(0, $l_pid)){
		next;
	    }
	    $content .= $line;
	}
	else{
	    #see if it's old
	    next if($stale_lock_timeout && $l_time + $stale_lock_timeout < $time);
	    $content .=$line;
	}
    }
    
    ### other shared locks exist
    if(length($content)){
	seek     $FH, 0, 0;
	print    $FH $id_line;
	print    $FH $content;
	truncate $FH, length($id_line.$content);
	close($FH);

	return 1;
    }
    else{  ### only I exist
	close($FH);
	$LOCK_LIST{$self->{lock_file}."_SH"} = 0; # just in case it's less than 0
	my $stat = unlink($self->{lock_file}) unless($self->_unlink_block);
	return $stat;
    }
}
#-------------------------------------------------------------------------------
sub _unlink_block {
    my $self = shift;
    my $arg = shift;

    if(defined $arg){
	$self->{_unlink_block} = $arg;
    }

    return $self->{_unlink_block};
}
#-------------------------------------------------------------------------------
sub uncache {
    my $self = shift;
    my $file = shift || $self->{lock_file};
    my $rand_file = $self->_rand_name($file);

    #pretend success when locking is off
    return 1 if(!$TYPE);

    ### hard link and rm trigger NFS to refresh
    my $stat;
    if($TYPE eq 'NFS'){
	$stat = link($file, $rand_file);
	unlink($rand_file);
    }
    elsif($TYPE eq 'POSIX'){
	open(my $fh, ">>$file");
	my $locked = flock($fh, LOCK_EX);
	$stat = File::Copy::move($rand_file, $file) if($locked);
	unlink($rand_file);
	close($fh);
    }

    return $stat;
}
#-------------------------------------------------------------------------------
sub maintain {
    my $self = shift;
    my $time = shift; #refresh interval
    my $pid = $self->{lock_pid};

    #pretend success when locking is off
    return 1 if(!$TYPE);

    die "ERROR: No time interval given to maintain lock\n\n"
	if(! $time || $time < 1);

    #clean up old maintainers
    local $?; #localize to clear locks

    my $stat;
    my $p = Proc::Signal::get_proc_by_id($self->{_maintain}) if(defined $self->{_maintain});

    if($p && $p->state ne 'zombie'){
	#clean up any children that are floating around
	do {
	    $stat = waitpid(-1, WNOHANG);
	} while $stat > 0;

	return 1;
    }
    elsif($p){
	close($self->{_IN}) if(ref $self->{_IN} eq 'GLOB');

	#signal maintainer
	kill(SIGUSR1, $self->{_maintain});
	$stat = waitpid($self->{_maintain}, WNOHANG);

	#attempt kill multiple times if still running
	my $count = 0;
	while($stat == 0 && $count < 200){
	    kill(SIGUSR1, $self->{_maintain}); #try multiple signal ending in signal 9
	    usleep(0.1) if($stat == 0);
	    $stat = waitpid($self->{_maintain}, WNOHANG);
	    $count++;
	}
	
	#if still running, do this
	if($stat == 0){
	    kill(SIGKILL, $self->{_maintain});
	    waitpid($self->{_maintain}, 0);
	}
	
	$self->{_IN} = undef;
	$self->{_maintain} = undef;
    }

    my $mine = $self->refresh; #refresh once on it's own

    #clean up any children that are floating around
    do {
	$stat = waitpid(-1, WNOHANG);
    } while $stat > 0;

    return 0 if(!$mine);

    #create lock serialization
    my $serial = Storable::freeze($self);
    $serial = uri_escape($serial, "\0-\377");

    #run maintainer executable in background
    my $exe = $INC{'File/NFSLock.pm'};
    die "ERROR: NFSLock does not appear to be loaded via use File::NFSLock\n\n"
	if(! $exe || ! -f $exe);

    $exe =~ s/NFSLock\.pm$/maintain\.pl/;
    $exe =~ s/ /\\ /g; #for paths with spaces
    my $m_pid = open3(my $IN, '>&STDOUT', '>&STDERR', "$^X $exe $pid $time $serial");

    #attach maintainer
    if($m_pid && $self->still_mine){
	$self->{_maintain} = $m_pid;
	$self->{_IN} = $IN;

	return 1;
    }
    else{
	close($IN);
	
	#signal maintainer
	kill(SIGUSR1, $self->{_maintain});
	$stat = waitpid($self->{_maintain}, WNOHANG);
	
	#attempt kill multiple times if still running
	my $count = 0;
	while($stat == 0 && $count < 200){
	    kill(SIGUSR1, $self->{_maintain}); #try multiple signal ending in signal 9
	    usleep(0.1) if($stat == 0);
	    $stat = waitpid($self->{_maintain}, WNOHANG);
	    $count++;
	}
	
	#if still running, do this
	if($stat == 0){
	    kill(SIGKILL, $self->{_maintain});
	    waitpid($self->{_maintain}, 0);
	}
	
	$self->{_IN} = undef;
	$self->{_maintain} = undef;

	return 0;
    }
}
#-------------------------------------------------------------------------------
sub refresh {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    my $old_lock_line = $self->{lock_line};
    my $id = $self->{id};
    my $pid = $self->{lock_pid};
    my $hostname   = $self->{hostname};
    my $shared = ($self->{lock_type} == LOCK_SH) ? 1 : 0;
    my $stale_lock_timeout = $self->{stale_lock_timeout};

    ### lock the locking process
    local $LOCK_EXTENSION = '.share';
    my $lock = new File::NFSLock (_fix_name($lock_file),LOCK_EX,62,60);
    
    #make sure the lock really is still mine (also refreshes NFS cache)
    return undef if(! $self->still_mine);

    $self->{lock_line} = "$hostname $pid ".time()." $id $shared\n";

    #pretend success when locking is off
    return 1 if(!$TYPE);    

    ### get the handle on the lock file
    my $FH;
    if(! open ($FH,"+<$lock_file") ){
	if( ! -f $lock_file || $! =~ /No such file or directory/){
	    return 0;
	}else{
	    die "Could not open for writing lock file $lock_file ($!)";
	}
    }
    
    my $id_line = '';
    my $content = '';
    ### read existing file
    if($shared){
	local $/ = "\n"; #just incase the user changed this
	my %seen;
	while(defined(my $line=<$FH>)){
	    next if($seen{$line}++);
	    if($line !~ / /){
		$id_line .= $line;
		next;
	    }

	    next if($line eq $old_lock_line);
	    next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
	    
	    my $l_host = $1;
	    my $l_pid  = $2;
	    my $l_time = $3;
	    my $l_id   = $4;   
	    my $l_kind = $5;
	    
	    if($l_host eq $hostname){
		if ($l_pid == $pid) { # This is me.
		    $content .= $line;
		}
		elsif(kill(0, $l_pid)) { # Still running on this host.
		    $content .=$line;
		}
		else{ # Finished running on this host.
		    next;
		}
	    }
	    else{
		# see if it's stale
		if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
		    next;
		}
		else{ # assume it is still running.
		    $content .= $line;
		}
	    }	    
	}
    }
    ### add new tag
    $content .= $self->{lock_line};
    
    ### fix file contents
    seek     $FH, 0, 0;
    print    $FH $id_line;
    print    $FH $content;
    truncate $FH, length($id_line.$content);
    close    $FH;
    
    return 1;
}
#-------------------------------------------------------------------------------
sub still_mine {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    my $lock_type = $self->{lock_type};
    my $lock_line = $self->{lock_line};
    my $pid = $self->{lock_pid};
    my $hostname   = $self->{hostname};
    my $id = $self->{id};

    #pretend success when locking is off
    return 1 if(!$TYPE);

    return 0 if(!-f $lock_file);

    my $shared = $self->_check_shared;
    return 0 if(!$shared && $lock_type == LOCK_SH);
    return 0 if($shared && $lock_type != LOCK_SH);

    #refresh NFS cache on the lock file
    $self->uncache($self->{lock_file});

    ### get the handle on the lock file
    my $FH;
    if(! open($FH,"+< $lock_file")){
	if(!-f $lock_file || $! =~ /No such file or directory/){
	    return 0;
	}else{
	    die "Could not open for reading the lock file $lock_file ($!)";
	}
    }

    ### read existing file
    my $mine = 0;
    local $/ = "\n"; #just incase the user changed this
    while(defined(my $line = <$FH>)){
	if ($line =~ /^$hostname $pid \d+/){
	    $mine = 1;
	    last;
	}
    }

    close($FH);

    return $mine;
}
#-------------------------------------------------------------------------------
#checks that lock is still yours when executing a command
sub do_safe {
    my $self = shift;
    my $subr = shift;
    
    if(ref $subr ne 'CODE'){
	die "ERROR: You must supply a CODE glob to File::NFSLock::do_safe\n"
    }

    if($self->still_mine){
	&{$subr}; #execute
    }
    else{
	die "ERROR: The lock is not yours to run the CODE glob\n";
    }
}
#-------------------------------------------------------------------------------
sub owners {
    my $self = shift;
    my $lock_file = $self->{lock_file};
    my $lock_line = $self->{lock_line};
    my $stale_lock_timeout = $self->{stale_lock_timeout};
    my $hostname   = $self->{hostname};
    my $pid = $self->{lock_pid};
    my $id = $self->{id};

    #pretend success when locking is off
    return 1 if(!$TYPE);

    return 1 unless ($self->{lock_type} eq LOCK_SH);

    ### lock the parsing process
    local $LOCK_EXTENSION = '.share';
    my $lock = new File::NFSLock (_fix_name($lock_file),LOCK_EX,62,60);

    #refresh NFS cache on the lock file
    $self->uncache( $self->{lock_file});

    #ignore other shared hosts where lock appears to be stale
    if ($stale_lock_timeout && time()-(stat($lock_file))[9] > $stale_lock_timeout){
	return $self->still_mine;
    }

    ### get the handle on the lock file
    my $FH;
    if( ! open ($FH,"+< $lock_file") ){
	if( ! -f $lock_file || $! =~ /No such file or directory/){
	    return 0;
	}else{
	    die "Could not open for reading the lock file $lock_file ($!)";
	}
    }

    my %seen;
    local $/ = "\n"; #just incase the user changed this
    while(defined(my $line=<$FH>)){
	next unless($line =~ /^([^\s]+) (\d+) (\d+) (\d+) (\d+)/);
	
	my $l_host = $1;
	my $l_pid  = $2;
	my $l_time = $3;
	my $l_id   = $4;   
	my $l_kind = $5;
	
	if($l_host eq $hostname){
	    if ($l_pid == $pid) { # This is me.
		$seen{"$l_host\_$l_pid"}++;
		next;
	    }
	    elsif($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
		next;
	    }
	    elsif(kill(0, $l_pid)) { # Still running on this host.
		$seen{"$l_host\_$l_pid"}++;
		next;
	    }
	    else{ # Finished running on this host.
		next;
	    }
	}
	else{
	    # see if it's old
	    if($stale_lock_timeout && $l_time + $stale_lock_timeout < time){
		next;
	    }
	    else{
		$seen{"$l_host\_$l_pid"}++;
		next;
	    }
	}
    }
    close($FH);

    return scalar(keys %seen);
}
#-------------------------------------------------------------------------------
sub shared_id {
    my $self = shift;
    my $line = $self->{id_line};
    return undef if(!$line);

    chomp $line;
    return $line
}
#-------------------------------------------------------------------------------
sub usleep {
    my $time = shift;
    select(undef,undef,undef,$time);
}
#-------------------------------------------------------------------------------
sub _fix_name {
    my $file = shift;
    my ($name) = $file =~ /([^\/]+)$/;
    if(length($name) > 150){	
	$name = Digest::MD5::md5_hex($name);
	$file =~ s/([^\/]+)$/$name/;
    }
    return($file);
}
#-------------------------------------------------------------------------------
sub DESTROY {
    shift()->unlock();
}
#-------------------------------------------------------------------------------
END {
    
}
#-------------------------------------------------------------------------------
1;
