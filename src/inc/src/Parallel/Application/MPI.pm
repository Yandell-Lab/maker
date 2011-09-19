package Parallel::Application::MPI;

use strict;
use Carp;
use vars qw(@ISA $VERSION %EXPORT_TAGS @EXPORT_OK);
use Storable qw(nfreeze thaw); #for complex datastructures
use Perl::Unsafe::Signals; #stops zombie processes under hydra MPICH2
require Exporter;

@ISA = qw(Exporter);
$VERSION = '0.01';

#for Parallel::Application specific interface
my @app_functions = qw(App_Init
		       App_Finalize
		       App_Comm_rank
		       App_Comm_size
		       App_Send
		       App_Recv);

my @app_constants = qw(App_ANY_SOURCE
                       App_ANY_TAG);

#for independent interface
my @functions = qw(MPI_Init
		   MPI_Finalize
		   MPI_Comm_rank
		   MPI_Comm_size
		   MPI_Send
		   MPI_Recv);

my @constants = qw(MPI_ANY_SOURCE
                   MPI_ANY_TAG);

%EXPORT_TAGS = ( all => [@functions, @constants], app => [@app_functions, @app_constants] );
@EXPORT_OK = ( @functions, @constants, @app_functions, @app_constants );

#for Parallel::Application specific interface
sub App_Init {
    require Proc::Signal;
    if(Proc::Signal::get_pname_by_id($$) =~ /^(mpiexec|mpirun|mpdrun|mpdexec)$/){
	MPI_Init(@_);
	return 1;
    }

    return 0;
}
sub App_Finalize { return MPI_Finalize(@_); }
sub App_Comm_rank { return MPI_Comm_rank(@_); }
sub App_Comm_size { return MPI_Comm_size(@_); }
sub App_Send { return MPI_Send(@_); }
sub App_Recv { return MPI_Recv(@_); }
sub APP_ANY_TAG { return MPI_ANY_TAG(@_); }
sub APP_ANY_SOURCE { return MPI_ANY_SOURCE(@_); }

#for independent interface
use MAKER::ConfigData;
use Inline::Files;
use Inline (C => Config =>
	    CC => MAKER::ConfigData->config('MPICC'),
	    LD => MAKER::ConfigData->config('MPICC'),
	    INC => '-I'.MAKER::ConfigData->config('MPIDIR'));
use Inline C => 'BELOW';

sub MPI_ANY_SOURCE { return _MPI_ANY_SOURCE(); }
sub MPI_ANY_TAG { return _MPI_ANY_TAG(); }
sub MPI_Init { UNSAFE_SIGNALS { return _MPI_Init(); }; }
sub MPI_Finalize { UNSAFE_SIGNALS { return _MPI_Finalize(); }; }

sub MPI_Comm_rank {
    my $rank;
    UNSAFE_SIGNALS {
	$rank = _MPI_Comm_rank();
    };
    return $rank;
}

sub MPI_Comm_size {
    my $size;
    UNSAFE_SIGNALS {
	$size = _MPI_Comm_size();
    };
    return $size;
}

sub MPI_Send {
    my $buf = shift;
    my $dest = shift;
    my $tag = shift;

    #MPI_Send only works with references
    confess "ERROR: not a reference\n" if(ref($buf) eq '');
    confess "ERROR: must be a reference to a SCALAR or REF\n" if(ref($buf) !~ /^(REF|SCALAR)$/);

    #freeze anything that is not a scalar reference
    my $freeze = 0; #flag to let receiving node know whether to call thaw
    my $msg = $buf;

    if(ref($msg) ne 'SCALAR'){
	$freeze = 1;
	$msg = \ (nfreeze($$buf)); #MPI_Send only works with references
    }

    die if($msg =~ /\0/);
    my $len = length($$msg);
    UNSAFE_SIGNALS {
	_MPI_Send($msg, $len, $dest, $tag, $freeze);
    };
}

sub MPI_Recv {
    my $buf = shift;
    my $source = shift;
    my $tag = shift;
    my $freeze;
    
    confess "ERROR: Not a reference to a SCALAR\n" if (ref($buf) ne 'SCALAR' && ref($buf) ne 'REF');

    UNSAFE_SIGNALS {
	_MPI_Recv($buf, $source, $tag, \$freeze);
    };

    $$buf = thaw($$buf) if($freeze);
}

#binding to C MPI library
__C__
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <mpi.h>

double _MPI_ANY_SOURCE () {
    return (double)MPI_ANY_SOURCE;
}

double _MPI_ANY_TAG () {
    return (double)MPI_ANY_TAG;
}

int _MPI_Init () {
    return MPI_Init(&PL_origargc, &PL_origargv);
}

int _MPI_Finalize () {
    return MPI_Finalize();
}

int _MPI_Comm_rank () {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int _MPI_Comm_size () {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void _MPI_Send(SV *buf, int len, int dest, int tag, int freeze) {
    int flags[2] = { len, freeze };
    MPI_Send(flags, 2, MPI_INT, dest, tag, MPI_COMM_WORLD);
    MPI_Send(SvPVX(SvRV(buf)), len, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
}

void _MPI_Recv(SV* buf, int source, int tag, SV* freeze) {
    MPI_Status status;
    int flags[2];
     MPI_Recv(flags, 2, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
    int len = flags[0];
    sv_setiv(SvRV(freeze), flags[1]);
    char *msg = (char*)malloc((len+1)*sizeof(char));
    MPI_Recv(msg, len, MPI_CHAR, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
    sv_setpvn(SvRV(buf), msg, len);
    free(msg);
}
