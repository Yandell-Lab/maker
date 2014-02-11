package Parallel::Application::MPI;

use strict;
use Carp;
use vars qw(@ISA $VERSION %EXPORT_TAGS @EXPORT_OK $CODE $LOADED $WARNED $INITIALIZED $FINALIZED);
use Storable qw(nfreeze thaw); #for complex datastructures
use Perl::Unsafe::Signals; #stops zombie processes under hydra MPICH2
use File::Temp qw(tempdir);
use IPC::Open3;
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
    confess "LOGIC ERROR: Cannot call MPI_ANY_SOURCE if not\n".
	"compiled and running under mpiexec" if(! _load());

    return C_MPI_ANY_SOURCE();
}
sub MPI_ANY_TAG {
    confess "LOGIC ERROR: Cannot call MPI_ANY_TAG if not\n".
	"compiled and running under mpiexec" if(! _load());

    return C_MPI_ANY_TAG();
}
sub MPI_SUCCESS {
    confess "LOGIC ERROR: Cannot call MPI_SUCCESS if not\n".
	"compiled and running under mpiexec" if(! _load());

    return C_MPI_SUCCESS();
}
sub MPI_Init {
    my $stat = 0;
    if($$ != 0 && !$INITIALIZED && _load()){
	# allow signals to interupt blocked MPI calls
	UNSAFE_SIGNALS {
	    $stat = C_MPI_Init();
	};
	$INITIALIZED = 1;
    }
    return $stat;
}
sub MPI_Finalize {
    my $stat = 1;
    if($$ != 0 && ! $FINALIZED && _load()){
	# allow signals to interupt blocked MPI calls
	UNSAFE_SIGNALS {
	    $stat = C_MPI_Finalize();
	};
	$FINALIZED = 1;
    }
    return $stat;
}
sub MPI_Comm_rank {
    my $rank = 0;
    # allow signals to interupt blocked MPI calls
    if(_load()){
	UNSAFE_SIGNALS {
	    $rank = C_MPI_Comm_rank();
	};
    }
    return $rank;
}
sub MPI_Comm_size {
    my $size = 1;
    # allow signals to interupt blocked MPI calls
    if(_load()){
	UNSAFE_SIGNALS {
	    $size = C_MPI_Comm_size();
	};
    }
    return $size;
}
sub MPI_Send {
    my $buf = shift;
    my $dest = shift;
    my $tag = shift;

    confess "LOGIC ERROR: Cannot run MPI_Send if not\n".
	"compiled and running under mpiexec" if(! _load());

    #MPI_Send only works with references
    confess "ERROR: not a reference\n" if(ref($buf) eq '');
    confess "ERROR: must be a reference to a SCALAR or REF\n"
	if(ref($buf) !~ /^(REF|SCALAR)$/);

    my $msg = nfreeze($buf); #always serialize the message
    my $len = length($msg);
    my $stat;
    # allow signals to interupt blocked MPI calls
    UNSAFE_SIGNALS {
	$stat = C_MPI_Send(\$msg, $len, $dest, $tag);
    };
    confess "ERROR: MPI_Send failed with status $stat" if($stat != MPI_SUCCESS);

    return $stat;
}

sub MPI_Recv {
    my $buf = shift;
    my $source = shift;
    my $tag = shift;

    confess "LOGIC ERROR: Cannot run MPI_Recv if not\n".
	    "compiled and running under mpiexec" if(! _load());    

    confess "ERROR: Not a reference to a SCALAR\n"
	if (ref($buf) ne 'SCALAR' && ref($buf) ne 'REF');
    
    my $stat;
    # allow signals to interupt blocked MPI calls
    UNSAFE_SIGNALS {
	$stat = C_MPI_Recv($buf, $source, $tag);
    };
    confess "ERROR: MPI_Recv failed with status $stat" if($stat != MPI_SUCCESS);
    $$buf = ${thaw($$buf)};

    return $stat;
}

#binding to C MPI library
sub _load {
    return 1 if($LOADED);

    require Proc::Signal;

    my $name = Proc::Signal::get_pname_by_id($$);
    if($name =~ /(mpiexec|mpirun|mpdrun|mpdexec|mpd|smpd|orted|orterun|pmi_proxy|hydra|mpispawn|exec\d+)$/){
	require MAKER::ConfigData;
	my $mpi_support = MAKER::ConfigData->feature('mpi_support');
	if(! $mpi_support){
	    warn "** WARNING: You have not configured MAKER to run under MPI.\n".
		 "** Yet you are attempting to do so!!\n".
		 "**\n".
		 "** You need to configure MAKER by executing -->\n".
		 "**     perl $FindBin::RealBin/../src/Build.PL\n".
		 "** Then run -->\n".
		 "**     $FindBin::RealBin/../src/Build install\n\n" unless($WARNED);
	    $WARNED = 1; #turn this warning off now
	    return 0;
	}

        #Find self location so inline can use it
        my $loc = $INC{'Parallel/Application/MPI.pm'};
	$loc =~ s/\/*(lib\/)?Parallel\/Application\/MPI\.pm$//;

	#if mpi_override it set, I'm using non-default options
	if(MAKER::ConfigData->feature('mpi_override')){
	    my $I = "$loc/lib"; #self location
	    $loc = tempdir("MPI_XXXXXX", CLEANUP => 1, TMPDIR => 1); #override location
	    my $M = 'Parallel::Application::MPI'; #self
	    my $mpicc  = MAKER::ConfigData->config('MPICC');
	    my $mpidir = MAKER::ConfigData->config('MPIDIR');
	    my $extra  = MAKER::ConfigData->config('CCFLAGSEX') || '';

	    #first call in separate executable to avoid setting $& (messes up regex)
	    my $cmd = "$^X -I'${I}' -M${M} -e '${M}::_bind(qw($mpicc $mpidir $loc $extra))'";
	    my $pid = open3('<&STDIN', '>&STDOUT', my $ERR, $cmd);
	    my $err = join('', <$ERR>);
	    waitpid($pid, 0);
	    croak $err if($?);
	    print STDERR $err if($err);
	}

	#lock for first compilation only
	my $lock;
	if(! -f "$loc/lib/auto/Parallel/Application/MPI/MPI.so" &&
	   ! -f "$loc/lib/auto/Parallel/Application/MPI/MPI.bundle"){
	    require File::NFSLock;
	    $lock = new File::NFSLock("$loc/_MPI", 'EX', 300, 40) while(!$lock);
	}

	_bind(MAKER::ConfigData->config('MPICC'),
	      MAKER::ConfigData->config('MPIDIR'),
	      $loc,
	      MAKER::ConfigData->config('CCFLAGSEX'));

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
    my $extra = shift || '';

    eval{
	#this comment is just a way to force Inline::C to recompile on changing MPICC and MPIDIR
	my $comment = "void _comment() {\nchar comment[] = \"MPICC=$mpicc, MPIDIR=$mpidir, CCFLAGSEX=$extra\";\n}\n"; 
	Inline->bind(C => $CODE . $comment,
		     NAME => 'Parallel::Application::MPI',
		     DIRECTORY => $loc,
		     CC => $mpicc,
		     LD => $mpicc,
		     CCFLAGSEX => $extra,
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
    
    confess $err if($err);
    
    return 1;
}

$CODE = <<END;

#include <mpi.h>

double C_MPI_ANY_SOURCE () {
    return (double)MPI_ANY_SOURCE;
}

double C_MPI_ANY_TAG () {
    return (double)MPI_ANY_TAG;
}

double C_MPI_SUCCESS () {
    return (double)MPI_SUCCESS;
}

int C_MPI_Init () {
    int stat;
    stat = MPI_Init(&PL_origargc, &PL_origargv);
    return stat;
}

int C_MPI_Finalize () {
    return MPI_Finalize();
}

int C_MPI_Comm_rank () {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int C_MPI_Comm_size () {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

int C_MPI_Send(SV *buf, int len, int dest, int tag) {
    STRLEN len2  = (STRLEN)len;
    SV* scalar   = SvRV(buf);
    char* string = SvPV(scalar, len2);

    int stat;
    stat = MPI_Send(string, len, MPI_CHAR, dest, tag, MPI_COMM_WORLD);

    return stat;
}

int C_MPI_Recv(SV* buf, int source, int tag) {
    int len;
    MPI_Status status;
    MPI_Probe(source, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &len);
    char *msg = (char*)malloc((len)*sizeof(char));

    int stat;
    stat = MPI_Recv(msg, len, MPI_CHAR, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);

    if(stat != MPI_SUCCESS){
	free(msg);
        return stat;
    }

    SV* scalar = SvRV(buf);
    sv_setpvn(scalar, msg, len);
    free(msg);

    return stat;
}
END

1;
