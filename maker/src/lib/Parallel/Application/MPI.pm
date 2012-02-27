package Parallel::Application::MPI;

use strict;
use Carp;
use vars qw(@ISA $VERSION %EXPORT_TAGS @EXPORT_OK $CODE $LOADED $WARNED);
use Storable qw(nfreeze thaw); #for complex datastructures
use Perl::Unsafe::Signals; #stops zombie processes under hydra MPICH2
require Exporter;

@ISA = qw(Exporter);
$VERSION = '0.01';

use Inline;

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

%EXPORT_TAGS = ( all => [@functions, @constants],
		 app => [@app_functions, @app_constants] );
@EXPORT_OK = ( @functions, @constants, @app_functions, @app_constants );

#for Parallel::Application specific interface
sub App_Init { return MPI_Init(@_); }
sub App_Finalize { return MPI_Finalize(@_); }
sub App_Comm_rank { return MPI_Comm_rank(@_); }
sub App_Comm_size { return MPI_Comm_size(@_); }
sub App_Send { return MPI_Send(@_); }
sub App_Recv { return MPI_Recv(@_); }
sub APP_ANY_TAG { return MPI_ANY_TAG(@_); }
sub APP_ANY_SOURCE { return MPI_ANY_SOURCE(@_); }

#for independent interface
sub MPI_ANY_SOURCE {
    die "LOGIC ERROR: Cannot call MPI_ANY_SOURCE if not\n".
	"compiled and running under mpiexec" if(! _load());

    return _MPI_ANY_SOURCE();
}
sub MPI_ANY_TAG {
    die "LOGIC ERROR: Cannot call MPI_ANY_TAG if not\n".
	"compiled and running under mpiexec" if(! _load());

    return _MPI_ANY_TAG();
}
sub MPI_Init {
    if(_load()){
	UNSAFE_SIGNALS {
	    _MPI_Init();
	};

	return 1;
    }
    else{
	return 0;
    }
}
sub MPI_Finalize {
    if(_load()){
	UNSAFE_SIGNALS {
	    _MPI_Finalize();
	};
    }
}
sub MPI_Comm_rank {
    my $rank = 0;
    if(_load()){
	UNSAFE_SIGNALS {
	    $rank = _MPI_Comm_rank();
	};
    }
    return $rank;
}
sub MPI_Comm_size {
    my $size = 1;
    if(_load()){
	UNSAFE_SIGNALS {
	    $size = _MPI_Comm_size();
	};
    }
    return $size;
}
sub MPI_Send {
    my $buf = shift;
    my $dest = shift;
    my $tag = shift;

    die "LOGIC ERROR: Cannot run MPI_Send if not\n".
	"compiled and running under mpiexec" if(! _load());

    #MPI_Send only works with references
    confess "ERROR: not a reference\n" if(ref($buf) eq '');
    confess "ERROR: must be a reference to a SCALAR or REF\n"
	if(ref($buf) !~ /^(REF|SCALAR)$/);

    #freeze anything that is not a scalar reference
    my $freeze = 0; #flag to let receiving node know whether to call thaw
    my $msg = $buf;

    if(ref($msg) ne 'SCALAR'){
	$freeze = 1;
	$msg = \ (nfreeze($$buf)); #MPI_Send only works with references
    }

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

    die "LOGIC ERROR: Cannot run MPI_Recv if not\n".
	"compiled and running under mpiexec" if(! _load());    

    confess "ERROR: Not a reference to a SCALAR\n"
	if (ref($buf) ne 'SCALAR' && ref($buf) ne 'REF');

    UNSAFE_SIGNALS {
	_MPI_Recv($buf, $source, $tag, \$freeze);
    };

    $$buf = thaw($$buf) if($freeze);
}

#binding to C MPI library
sub _load {
    return 1 if($LOADED);

    require Proc::Signal;

    my $name = Proc::Signal::get_pname_by_id($$);
    if($name =~ /^(mpiexec|mpirun|mpdrun|mpdexec|mpd|smpd|orted|hydra_pmi_proxy)$/){
	require MAKER::ConfigData;
	my $mpi_support = MAKER::ConfigData->feature('mpi_support');
	if(! $mpi_support){
	    warn "** WARNING: You have not configured MAKER to run under MPI.\n".
		 "** Yet you are attempting to do so!!\n".
		 "**\n".
		 "** You need to configure MAKER by executing -->\n".
		 "**     perl $FindBin::Bin/../src/Build.PL\n".
		 "** Then run -->\n".
		 "**     $FindBin::Bin/../src/Build install\n\n" unless($WARNED);
	    $WARNED = 1; #turn this warning off now
	    return 0;
	}

        #Find self location so inline can use it
        my $loc = $INC{'Parallel/Application/MPI.pm'};
	$loc =~ s/\/*(lib\/)?Parallel\/Application\/MPI\.pm$//;

	#lock for first compilation only
	my $lock;
	if(! -f "$loc/lib/auto/Parallel/Application/MPI/MPI.bundle"){
	    require File::NFSLock;
	    $lock = new File::NFSLock("$loc/_MPI", 'EX', 300, 40) while(!$lock);
	}

	_bind(MAKER::ConfigData->config('MPICC'),
	      MAKER::ConfigData->config('MPIDIR'),
	      $loc);

	$LOADED = 1;
	$lock->unlock() if($lock);

	return 1;
    }

    return 0;
}

sub _bind {
    my $mpicc = shift;
    my $mpidir = shift;
    my $loc = shift;

    eval{
	#this comment is just a way to force Inline::C to recompile on changing MPICC and MPIDIR
	my $comment = "void _comment() {\nchar comment[] = \"MPICC=$mpicc, MPIDIR=$mpidir\";\n}\n"; 
	Inline->bind(C => $CODE . $comment,
		     NAME => 'Parallel::Application::MPI',
		     DIRECTORY => $loc,
		     CC => $mpicc,
		     LD => $mpicc,
		     INC => '-I'.$mpidir,);
    };
    
    my $err = $@;

    if($err =~ /cannot open shared object file|mca_base_param_reg_int/){	
	$err .= "\n\n".
	    "** If you are running using OpenMPI, you may have to preload object files\n".
	    "** for shared libraries to work. For bash, try executing a command\n".
	    "** similar to the following with the appropriate file location.\n".
	    "** Example --> export LD_PRELOAD=.../openmpi/lib/libmpi.so\n".
	    "** Please do this before trying to run MAKER again!!\n\n";
    }
    
    die $err if($err);
    
    return 1;
}

$CODE = <<END;

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
END

1;
