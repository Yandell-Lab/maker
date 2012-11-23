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
	require Proc::ProcessTable;
	return Proc::ProcessTable->new(@_);
    }
}

sub table {
    my $self = shift;

    my $cmd = "$PS -ax -o pid=pid".("_"x10).
	      " -o ppid=ppid".("_"x10).
	      " -o state=state".("_"x10).
	      " -o comm=fname".("_"x145).
	      " -o args=cmndline -w -w 2> /dev/null |";

    open(PS, $cmd);
    my $i = 0;
    my @cat;
    my @table;
    local $/ = "\n"; #just in case
    while(my $line = <PS>){
	if(!$i++){ #first line
	    @cat = $line =~ /^(\s*[^\s]+)(\s*[^\s]+)(\s*[^\s]+)(\s*[^\s]+)/;
	    @cat = map {length($_)} @cat;
	    next;
	}

	my @F = $line =~ /^(.{$cat[0]})(.{$cat[1]})(.{$cat[2]})(.{$cat[3]})(.*)$/;
	map {s/^\s+|\s+$//g} @F;

	next unless(defined $F[2]);

	if($F[2] =~ /^T/i){
	    $F[2] = 'stopped';
	}
	elsif($F[2] =~ /^Z/i){
	    $F[2] = 'zombie';
	}
	else{
	    $F[2] = '';
	}

	my $obj = Proc::Process_simple->new(@F);
	push(@table, $obj);
    }
    close(PS);
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

    $self->{pid}      = shift;
    $self->{ppid}     = shift;
    $self->{state}    = shift;
    $self->{fname}    = shift;
    $self->{cmndline} = shift;

    #weird trailing bug
    chomp($self->{fname});
    $self->{fname} =~ s/\s\s\s\s+.//;
    
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
