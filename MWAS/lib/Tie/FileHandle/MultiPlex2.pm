#! /usr/bin/perl

package Tie::FileHandle::MultiPlex2;

use base qw(Tie::FileHandle::Base);
use vars qw($VERSION);
$VERSION = 0.1;

# TIEHANDLE
# Usage: tie *HANDLE, 'Tie::FileHandle::Log', 'file_name', &GLOB, ...
sub TIEHANDLE {
    #test open
    my $self = [];
    foreach(@_[1..$#_]){
	if(ref(\$_) eq 'SCALAR'){
	    open(my $FH, "> $_") || die "ERROR: Could not open the log file:$file\n";
	    push(@$self, \*$FH);
	}
	elsif(ref(\$_) eq 'GLOB'){
	    push(@$self, \*$_);
	}
    }

    bless($self, $_[0]);
}

# PRINT
sub PRINT {
    print $_ $_[1] for @{ $_[0] };
}

1;
