#!/usr/local/bin/perl

=head1 NAME

Tie::FileHandle::Base - a base class to simplify filehandle tie module implementation

=head1 DESCRIPTION

By noting the redundancies inherent in the filehandle tie methods, this
module seeks to aid in implementation of new modules by reducing the number
of required functions.

Care should be taken by classes that use AUTOLOAD.  Make sure to predeclare
subroutines that will be autoloaded - as in:

    sub PRINT;

Otherwise this module will make incorrect presumptions and your module will
not function as you intend.

=head2 OUTPUT FUNCTIONS

Since PRINT, PRINTF, and WRITE are all quite similar in scope, any one of
these can be implemented from any of the others.  So, you only need implement
one of the above.

=head2 INPUT FUNCTIONS

By implementing READ or GETC, you can get the entire complement of READ,
READLINE, and GETC.   Note however that READ and GETC cannot be derived
nicely from READLINE.

=head2 OTHERS

EOF can be implemented crudely if given READ or GETC along
with a backwards supporting SEEK.

=head1 HISTORY

=over 2

=item *

03/09/02 - Robby Walker - did the output stuff - version 0.1

=item *

02/13/02 - Robby Walker - created the file - version 0.001

=back

=head1 METHODS

=over 4

=cut
#----------------------------------------------------------

package Tie::FileHandle::Base;

use vars qw($VERSION %loop);
use strict;

$VERSION = 0.1;
%loop = ();

# ------------------------------------------------------------------------
# METHOD: _fh_error
# ------------------------------------------------------------------------
# Our replacement for 'croak' or 'die' that errors to the proper level
sub _fh_error {
    my $error_string = shift;
    my $i = 1;
    while (my ($package, $filename, $line, $subroutine,
	               undef, undef, undef, undef, undef, undef) = caller($i++) )
    {
	if ( $package ne 'Tie::FileHandle::Base' ) {
	    $subroutine =~ s/.*\:([^:]+)/$1/;
	    die "Cannot execute filehandle method $subroutine : $error_string at $filename line $line\n";
	}
    }
}

# ------------------------------------------------------------------------
# METHOD: PRINT
# ------------------------------------------------------------------------

=item PRINT

Implements PRINT based on WRITE or PRINTF.

=cut
sub PRINT {
    my $self = shift;
    my $result = 0;

    # guard against loops
    _fh_error( "function not defined" ) if ( $loop{$self} );
    $loop{$self} = 1;

    if ( $self->can('WRITE') != \&WRITE ) {
	# loop over the strings
	$result = 1;
	foreach my $str ( @_ ) {
	    # print each string carefully
	    my $offset = 0;
	    my $ln = length( $str );
	    # loop until all characters are printed
	    while ( $offset != $ln ) {
		my $ret = $self->WRITE( $str, $ln - $offset, $offset );
		unless ( $ret ) {
		    $result = undef;
		    last;
		}
	    };
	    # see if we exited early
	    last if $offset != $ln;
	}

    } elsif ( $self->can('PRINTF') != \&PRINTF ) {
	$result = $self->PRINTF( '%s' x (@_+0), @_ );

    } else {
	_fh_error( "function not defined" );
    }

    $loop{$self} = 0;
    1;
}

# ------------------------------------------------------------------------
# METHOD: PRINTF
# ------------------------------------------------------------------------

=item PRINTF

Implements PRINTF based off of PRINT, which may in turn base itself off of WRITE.

=cut
sub PRINTF {
    ( shift )->PRINT( sprintf @_[1..$#_] );
}

# ------------------------------------------------------------------------
# METHOD: WRITE
# ------------------------------------------------------------------------

=item WRITE

Implements WRITE based off of PRINT, which may in turn base itself off of PRINTF.

=cut
sub WRITE {
    my ($self, $var, $length, $offset) = @_;
    my $ln = $length || length( $var );
    $self->PRINT( substr $var, $offset || 0, $ln ) && $ln;
}

# ------------------------------------------------------------------------
# METHOD: GETC
# ------------------------------------------------------------------------

=item GETC

=cut
sub GETC {

}

# ------------------------------------------------------------------------
# METHOD: READ
# ------------------------------------------------------------------------

=item READ

=cut
sub READ {

}

# ------------------------------------------------------------------------
# METHOD: READLINE
# ------------------------------------------------------------------------

=item READLINE

=cut
sub READLINE {

}

# ------------------------------------------------------------------------
# METHOD: EOF
# ------------------------------------------------------------------------

=item EOF

Crude EOF implemented using READ and SEEK.

=cut
sub EOF {
    my $self = shift;
    if ( $self->can('SEEK') ) {
	my $temp;
	# test EOF by reading
	return 1 unless $self->READ( $temp, 1 );
	$self->SEEK( 1, -1 ); # go back to where we were
    }
    return 0;
}



1;

__END__

=pod

=back

=head1 TODO

=over 4

=item *

Input stuff.

=item *

test.pl

=back

=head1 BUGS

This is a new module and has not been thoroughly tested.

=head1 AUTHORS AND COPYRIGHT

Written by Robby Walker ( robwalker@cpan.org ) for Point Writer ( http://www.pointwriter.com/ ).

You may redistribute/modify/etc. this module under the same terms as Perl itself.
