#! /usr/bin/perl 
package Proc::ProcessTable_simple;

require Exporter;
use strict;
use warnings;
use File::Which;
use Carp;
use base qw();

our $VERSION='1.0';
our $PS = File::Which::which('ps');

sub new {
    if($PS){
	my $self = {};
	my $class = shift;
	bless($self, $class);
	return $self;
    }
    else{
	eval 'require Proc::ProcessTable';
	return Proc::ProcessTable->new(@_);
    }
}

sub get_proc_by_id {
    my $self = shift;
    my $id = shift;

    return unless($id);

    my $list = $self->_make_procs($id);
    return shift @$list;
}

sub table {
    my $self = shift;
    return $self->_make_procs();
}

sub _make_procs {
    my $self = shift;
    my $pid = shift || '';

    my $cmd = "$PS";
    $cmd .= ($pid) ? " p$pid" : " ax"; 
    $cmd .= " -o pid=pid".("_"x10).
	    " -o ppid=ppid".("_"x10).
	    " -o state=state".("_"x10).
	    " -o vsz=size".("_"x10).
	    " -o comm=fname".("_"x140).
	    " -o args=cmndline".
	    " -w -w 2> /dev/null |";

    my $stat = open(my $EX, $cmd);

    #what to do on failure
    if(!$stat){
	sleep 5;
	$stat = open($EX, $cmd);
	die "ERROR: Could not query process table: $!" if(!$stat);
    }

    my $i = 0;
    my @tags;
    my @lens;
    my @table;
    local $/ = "\n"; #just in case
    while(my $line = <$EX>){
	if(!$i++){ #first line
	    @tags = map {s/[\s\n]//g; $_} split(/_+/, $line);
	    while($line =~ s/^(\s*[^\s]+)//){
		push(@lens, length($1));
	    }
	    next;
	}

	my %proc;
	for(my $j = 0; $j < @lens; $j++){
	    my $l = $lens[$j];
	    my $val;
	    if($j == $#lens){
		$val = $line;
		chomp($val);
	    }
	    elsif($line =~ s/^(.{$l})//){
		$val = $1;
	    }
	    $proc{$tags[$j]} = $val;
	}
	next unless(defined $proc{state});

	map {s/^\s+|\s+$//g} values %proc;

	if($proc{state} =~ /^T/i){
	    $proc{state} = 'stopped';
	}
	elsif($proc{state} =~ /^Z/i){
	    $proc{state} = 'zombie';
	}
	else{
	    $proc{state} = '';
	}

	$proc{size} *= 1000; #convert to bytes

	#weird trailing character bug
	chomp($proc{cmndline});
	chomp($proc{fname});
	$proc{fname} =~ s/\s\s\s\s+.//;

	my $obj = Proc::Process_simple->new(\%proc);
	push(@table, $obj);
    }
    close($EX);
    return \@table;
}

sub Destroy {
}

package Proc::Process_simple;

require Exporter;
use strict;
use warnings;
use Carp;
use base qw();

sub new {
    my $self = {};
    my $class = shift;
    bless($self, $class);

    my $hash = shift;

    foreach my $key (keys %$hash){
	$self->{$key} = $hash->{$key};
    }
    
    return $self;
}

our $AUTOLOAD;
sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) or croak "$self is not an object";
    (my $name = $AUTOLOAD) =~ s/.*://; # strip fully-qualified portion
    return if($name eq 'DESTROY');

    die "ERROR: $name called from $type" if(! exists $self->{$name});
    if(@_){
	return $self->{$name} = shift;
    }
    else{
        return $self->{$name};
    }
}

sub DESTROY {
}

1;
