#! /usr/bin/perl -w
package Proc::Signal;

require Exporter;
use strict;
use vars qw(@EXPORT @EXPORT_OK @ISA $VERSION);
use Proc::ProcessTable;
use URI::Escape;

@EXPORT_OK=qw(signal killall signalall exists_proc_by_id exists_proc_by_name);
@ISA=qw(Exporter);

$VERSION='1.0';
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
    my $id = shift;
    my $signal = shift;

    return (0,0) if(! $id || $id == $$);

    my $obj = new Proc::ProcessTable;
    my $found = 0;
    my $signaled = 0;

    foreach my $p (@{$obj->table}) {
	#now check for the id
	if ($p->pid == $id){
	    $found++;
	    $signaled = kill($signal, $id);
	    last;
	}
    }

    return ($found, $signaled);
}
#-----------------------------------------------------------------
#returns process name for a given ID
sub get_name_by_id {
    my $id = shift;

    my $obj = new Proc::ProcessTable;
    my $name = '';
    foreach my $p (@{$obj->table}) {
	#now check for the id
	if ($p->pid == $id){
	    $name   = $p->fname() || '';
	}
    }

    return $name;
}
#-----------------------------------------------------------------
#this function is used by killall and exists_proc_by_name. it
#returns processes that were found and how many can be signaled
sub exists_killall {
    my $signal = shift;
    my $name = shift || '';

    return (0, 0) if($name eq '');

    my $obj = new Proc::ProcessTable;

    my $err; #holds all $!
    my $found = 0; #holds count of found procceses
    my $signaled = 0; #holds count of killed procceses

    foreach my $p (@{$obj->table}) {
	#get all possible name variations for proccess
	my $cmdline  = $p->cmndline() || '';
	my ($front)  = $p->cmndline() =~ /^([^\s]+)/ || '';
	my ($f_name) = $front =~ /([^\/]+)$/ || '';
	my $exe      = ''; #$p->exec() || '';
	my $e_name   = $p->fname() || '';
	my $script   = '';
	my $s_name   = '';

	#split command line into arguments
	if(($e_name =~ /^perl$/ ||
	    $exe    =~ /^perl$/ ||
	    $f_name =~ /^perl$/ ||
	    $front  =~ /^perl$/) &&
	   $name ne 'perl'
	  ){
	    #must handles escaped characters
	    while($cmdline =~ /(\\.)/){
		my $replace = uri_escape($1, "\0-\377");
		$cmdline =~ s/\\./$replace/;
	    }

	    #must handle quotes
	    my @quotes = grep {$_ =~ s/(\'|\"|\`)/\\$1/} $cmdline =~ /(.)/g;
	    while(my $q = shift @quotes){
		if($cmdline =~ /($q[^$q]+$q)/){
		    my $replace = uri_escape($1, "\0-\377");
		    $cmdline =~ s/\\./$replace/;
		    @quotes = grep {/\'|\"|\`/} $cmdline =~ /(.)/g;
		}
	    }

	    #now I can split on spaces and recover values
	    my @args = split(/\s+/, $cmdline);
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
	    $front  =~ /^$name$/ ||
	    $s_name =~ /^$name$/ ||
	    $script =~ /^$name$/
	   ) {
	    next if($p->pid == $$ || $p->pid == 0);
	    my $num = kill ($signal, $p->pid);
	    $found++;
	    $signaled += $num;
	    $err = $! if(! $num); 
	}
    }

    $! = $err if($err);

    return ($found, $signaled);
}
#-----------------------------------------------------------------
1;
