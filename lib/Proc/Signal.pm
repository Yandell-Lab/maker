#! /usr/bin/perl -w
package Proc::Signal;

require Exporter;
use strict;
use vars qw(@EXPORT @EXPORT_OK @ISA $VERSION $PS);
use URI::Escape;
use File::Which;
use Proc::ProcessTable_simple;

@EXPORT_OK=qw(signal killall signalall exists_proc_by_id exists_proc_by_name);
@ISA=qw(Exporter);

$VERSION='1.0';
$PS = File::Which::which('ps');
#-----------------------------------------------------------------
sub VERSION {
    return $VERSION;
}
#-----------------------------------------------------------------
#sends a signal to proccess id
#basically an alias to kill (just sounds less dangerous)
sub signal {
    return kill(@_);
}
#-----------------------------------------------------------------
#sends a signal to all proccesses with a given name
#basically an alias to killall (just sounds less dangerous)
sub signalall {
    return killall(@_);
}
#-----------------------------------------------------------------
sub id_matches_pattern{
    my $pid = shift;
    my $pattern = shift;

    my $stat;   
    if($PS){
	$stat= grep {/$pattern/} `ps -o args -p $pid 2> /dev/null`;
    }
    else{
	my $obj = new Proc::ProcessTable_simple;	
	foreach my $p (@{$obj->table}) {
	    if ($p->pid == $pid){
		$stat = grep {/$pattern/} $p->cmndline();
		last;
	    }
	}
    }
    
    return $stat;
}
#-----------------------------------------------------------------
#checks to see if a process exists by id
sub exists_proc_by_id {
    my ($found, $signaled) = exists_kill(0, @_);

    return $found;
}
#-----------------------------------------------------------------
#checks to see if a process exists by name
#returns the number of procceses found
sub exists_proc_by_name {
    my ($found, $signaled) = exists_killall(0, @_);

    return $found;
}
#-----------------------------------------------------------------
#this function will perform killall similar to the linux command
#killall which kills all processes of a given name, but it will
#act like the perl command kill (i.e. take 0) and it returns the
#number of procceses signaled (not necessarily that died).
sub killall {
    my ($found, $signaled) = exists_killall(@_);

    return $signaled;
}
#-----------------------------------------------------------------
#checks to see if a process exists by id
#essentially kill with a signal of 0
sub exists_kill {
    my $signal = shift;
    my $id = shift;

    return (0,0) if(! $id || $id == $$);


    my $found = 0;
    my $signaled = 0;
    if($PS){
        if(grep {/$id/} `ps -o pid -p $id 2> /dev/null`){
	    $found++;
	    $signaled = kill($signal, $id);
	}
    }
    else{    
	my $obj = new Proc::ProcessTable_simple;
	foreach my $p (@{$obj->table}) {
	    #now check for the id
	    if ($p->pid == $id){
		$found++;
		$signaled = kill($signal, $id);
		last;
	    }
	}
    }

    return ($found, $signaled);
}
#-----------------------------------------------------------------
#returns process name for a given ID
sub get_name_by_id {
    my $id = shift;

    my $p = get_proc_by_id($id);
    my $name = '';
    if($p){
	$name   = $p->fname() || '';
    }

    return $name;
}
#-----------------------------------------------------------------
#returns parent process name for a given child ID
sub get_pname_by_id {
    my $id = shift;

    my $p = get_pproc_by_id($id);
    my $name = '';
    if($p){
	$name   = $p->fname() || '';
    }

    return $name;
}
#-----------------------------------------------------------------
#returns process table for a given process id
sub get_proc_by_id {
    my $id = shift;

    my $select;
    my $obj = new Proc::ProcessTable_simple;
    return $obj->get_proc_by_id($id);
}
#-----------------------------------------------------------------
#returns parent process table for a given child id
sub get_pproc_by_id {
    my $id = shift;
    my $p = get_proc_by_id($id);

    return get_proc_by_id($p->ppid) if($p);
    return undef;
}
#-----------------------------------------------------------------
#this function is used by killall and exists_proc_by_name. it
#returns processes that were found and how many can be signaled
sub exists_killall {
    my $signal = shift;
    my $name = shift || '';
    my $parent = shift || '';
    my $wait = shift || 0;

    return (0, 0) if($name eq '');

    my $err; #holds all $!
    my $found = 0; #holds count of found procceses
    my $signaled = 0; #holds count of killed procceses

    my $obj = new Proc::ProcessTable_simple;
    foreach my $p (@{$obj->table}) {
	#get all possible name variations for proccess
	my $cmdline  = $p->cmndline() || '';
	my $exe      = ($p->cmndline() =~ /^([^\s]+)/) ? $1 : '';
	my $f_name   = ($exe =~ /([^\/]+)$/) ? $1 : '';
	my $e_name   = $p->fname() || ''; #sometimes this is the value after perl or python etc.
	my $script   = '';
	my $s_name   = '';

	#split command line into arguments
	if(($e_name =~ /^perl$/ ||
	    $exe    =~ /^perl$/ ||
	    $f_name =~ /^perl$/) &&
	   $name ne 'perl'
	  ){
	    #must handles escaped characters
	    while($cmdline =~ /(\\.)/){
		my $replace = uri_escape($1, "\0-\377");
		$cmdline =~ s/\\./$replace/;
	    }

	    #must handle quotes
	    my @quotes = grep {$_ =~ s/(\'|\"|\`)/\\$1/} $cmdline =~ /(.)/g;
	    while(@quotes){
		my $q = shift @quotes;
		if($cmdline =~ /($q[^$q]+$q)/){
		    my $replace = uri_escape($1, "\0-\377");
		    $cmdline =~ s/$q[^$q]+$q/$replace/;
		    @quotes = grep {/\'|\"|\`/} $cmdline =~ /(.)/g;
		}
		@quotes = grep {$_ =~ s/(\'|\"|\`)/\\$1/} $cmdline =~ /(.)/g;
	    }

	    #now I can split on spaces and recover values
	    my @args = split(/\s+/, $cmdline);
	    shift @args; #remove first value because it is perl
	    foreach my $arg (@args){
		$arg = uri_unescape($arg);
		if($arg !~ /^-/){
		    $script = $arg;
		    ($s_name) = $arg =~ /([^\/]+)$/; 
		    last;
		}
	    }
	}

	#now check for the names
	if ($e_name =~ /^$name$/ ||
	    $exe    =~ /^$name$/ ||
	    $f_name =~ /^$name$/ ||
	    $s_name =~ /^$name$/ ||
	    $script =~ /^$name$/
	   ) {
	    next if($p->pid == $$ || $p->pid == 0);

	    next if($parent && $parent != $p->ppid);

	    my $num = kill ($signal, $p->pid);
	    waitpid($p->pid, 0) if($wait);

	    $found++;
	    $signaled += $num;
	    $err = $! if(! $num); 
	}
    }

    $! = $err if($err);

    return ($found, $signaled);
}
#-----------------------------------------------------------------
#reeps zombies belonging to me
sub reap_children_by_name {
    my $signal = shift;
    my $name = shift;

    return exists_killall($signal, $name, $$, 1);
}

#-----------------------------------------------------------------
1;
