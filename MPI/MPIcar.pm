package Parallel::MPIcar;

use strict;
use Carp;
use vars qw($VERSION @ISA %EXPORT_TAGS @EXPORT_OK $AUTOLOAD 
	    $errno $errstr $exceptions);

# whether or not to throw exceptions in the MPI functions.
$exceptions = 1;

# these will be set (regardless of the above setting) if an error occurs in
# an MPI function.
$errno = 0;
$errstr = undef;

require Exporter;
require DynaLoader;
require AutoLoader;

@ISA = qw(Exporter DynaLoader);

my %constants = qw(MPI_2COMPLEX            MPI_Datatype
		   MPI_2DOUBLE_COMPLEX     MPI_Datatype
		   MPI_2DOUBLE_PRECISION   MPI_Datatype
		   MPI_2INT                MPI_Datatype
		   MPI_2INTEGER            MPI_Datatype
		   MPI_2REAL               MPI_Datatype
		   MPI_COMPLEX             MPI_Datatype
		   MPI_DATATYPE_NULL       MPI_Datatype
		   MPI_DOUBLE              MPI_Datatype
		   MPI_DOUBLE_COMPLEX      MPI_Datatype
		   MPI_DOUBLE_INT          MPI_Datatype
		   MPI_DOUBLE_PRECISION    MPI_Datatype
		   MPI_FLOAT               MPI_Datatype
		   MPI_FLOAT_INT           MPI_Datatype
		   MPI_INT                 MPI_Datatype
		   MPI_INTEGER             MPI_Datatype
		   MPI_BYTE                MPI_Datatype
		   MPI_CHAR                MPI_Datatype
		   MPI_CHARACTER           MPI_Datatype
		   MPI_LOGICAL             MPI_Datatype
		   MPI_LONG                MPI_Datatype
		   MPI_LONG_DOUBLE         MPI_Datatype
		   MPI_LONG_DOUBLE_INT     MPI_Datatype
		   MPI_LONG_INT            MPI_Datatype
		   MPI_LONG_LONG_INT       MPI_Datatype
		   MPI_REAL                MPI_Datatype
		   MPI_SHORT               MPI_Datatype
		   MPI_SHORT_INT           MPI_Datatype
                   MPI_STRING              MPI_Datatype
		   MPI_UNSIGNED            MPI_Datatype
		   MPI_UNSIGNED_CHAR       MPI_Datatype
		   MPI_UNSIGNED_LONG       MPI_Datatype
		   MPI_UNSIGNED_SHORT      MPI_Datatype
		   MPI_ANY_SOURCE          MPI_Status
		   MPI_ANY_TAG             MPI_Status
		   MPI_BAND                MPI_Op
		   MPI_BOR                 MPI_Op
		   MPI_BXOR                MPI_Op
		   MPI_LAND                MPI_Op
		   MPI_LOR                 MPI_Op
		   MPI_LXOR                MPI_Op
		   MPI_MAX                 MPI_Op
		   MPI_MAXLOC              MPI_Op
		   MPI_MIN                 MPI_Op
		   MPI_MINLOC              MPI_Op
		   MPI_OP_NULL             MPI_Op
		   MPI_PROD                MPI_Op
		   MPI_SUM                 MPI_Op
		   MPI_COMM_NULL           MPI_Comm
		   MPI_COMM_SELF           MPI_Comm
		   MPI_COMM_WORLD          MPI_Comm
		   MPI_CONGRUENT           undef
		   MPI_IDENT               undef
		   MPI_SIMILAR             undef
		   MPI_UNEQUAL             undef
		   MPI_VERSION             undef);

my @funcs =     qw(&MPI_Send
		   &MPI_Recv
		   &MPI_Barrier
		   &MPI_Bcast
		   &MPI_Comm_size
		   &MPI_Comm_rank
		   &MPI_Wtime
		   &MPI_Wtick
		   &MPI_Init
		   &MPI_Finalize
		   &MPI_Initialized
		   &MPI_Abort
		   &MPI_Reduce
		   &MPI_Allreduce
		   &MPI_Scatter
		   &MPI_Gather
		   &MPI_Sendrecv);

%EXPORT_TAGS = ( all => [ keys %constants, @funcs ] );
@EXPORT_OK = ( keys %constants, @funcs );

$VERSION = '0.04mod';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.  If a constant is not found then control is passed
    # to the AUTOLOAD in AutoLoader.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    my $val = constant($constname, @_ ? $_[0] : 0);

    if ($! != 0) {
	if ($! =~ /Invalid/) {
	    $AutoLoader::AUTOLOAD = $AUTOLOAD;
	    goto &AutoLoader::AUTOLOAD;
	}
	else {
	    croak "Your vendor has not defined MPI constant $constname";
	}
    }

    # some constants need to be blessed references to allow type checking.
    if ($constants{$constname} ne "undef") {
	eval "sub $constname {  my \$v = $val; my \$v2 = \\\$v; bless \$v2, \"$constants{$constname}\"; }";	
    } else {
	eval "sub $AUTOLOAD { $val }";
    }
    goto &$AUTOLOAD;
}

bootstrap Parallel::MPIcar $VERSION;

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__

=head1 NAME

Parallel::MPIcar - Perl interface to the MPI message passing system

=head1 SYNOPSIS

  use Parallel::MPIcar;
  MPI_Init();
  . . .
  MPI_Finalize();

=head1 DESCRIPTION

The following is a summary of the available constants and functions:

=head1 Error Handling

If an MPI error occurs, set:
$Parallel::MPIcar::errno
$Parallel::MPIcar::errstr

$Parallel::MPIcar::exceptions: if set, toss an exception when an error occurs.

=head1 Exported constants

   Datatypes (not all are supported!)

   MPI_2COMPLEX           
   MPI_2DOUBLE_COMPLEX    
   MPI_2DOUBLE_PRECISION  
   MPI_2INT               
   MPI_2INTEGER           
   MPI_2REAL              
   MPI_COMPLEX            
   MPI_DATATYPE_NULL      
   MPI_DOUBLE             
   MPI_DOUBLE_COMPLEX     
   MPI_DOUBLE_INT         
   MPI_DOUBLE_PRECISION   
   MPI_FLOAT              
   MPI_FLOAT_INT          
   MPI_INT                
   MPI_INTEGER            
   MPI_BYTE               
   MPI_CHAR               
   MPI_CHARACTER          
   MPI_LOGICAL            
   MPI_LONG               
   MPI_LONG_DOUBLE        
   MPI_LONG_DOUBLE_INT    
   MPI_LONG_INT           
   MPI_LONG_LONG_INT      
   MPI_REAL               
   MPI_SHORT              
   MPI_SHORT_INT          
   MPI_UNSIGNED           
   MPI_UNSIGNED_CHAR      
   MPI_UNSIGNED_LONG       
   MPI_UNSIGNED_SHORT

   New Datatypes

   MPI_STRING
   
   Status
   
   MPI_ANY_SOURCE    
   MPI_ANY_TAG       
   
   Operations
   
   MPI_BAND   
   MPI_BOR    
   MPI_BXOR   
   MPI_LAND   
   MPI_LOR    
   MPI_LXOR   
   MPI_MAX    
   MPI_MAXLOC 
   MPI_MIN    
   MPI_MINLOC 
   MPI_OP_NULL
   MPI_PROD   
   MPI_SUM    
   
   Communicators
   
   MPI_COMM_NULL
   MPI_COMM_SELF
   MPI_COMM_WORLD
   
   Communicator and Group Comparisons
   
   MPI_CONGRUENT 
   MPI_IDENT    
   MPI_SIMILAR  
   MPI_UNEQUAL  
   MPI_VERSION  

=head1 Exported functions

   MPI_Init()
   MPI_Finalize()
   MPI_Initialized()

   MPI_Comm_rank(communicator)
   MPI_Comm_size(communicator)
   
   MPI_Send(\$message, length, datatype, destination, tag, communicator)
   MPI_Recv(\$message, length, datatype, source, tag, communicator)
   MPI_Sendrecv(\$message, length, datatype, destination, tag, communicator)
   
   MPI_Barrier(comm)
   MPI_Bcast(\$from, count, datatype, root, communicator)
 
   MPI_Wtime()
   MPI_Wtick()

   MPI_Abort(communicator, errorcode)

   MPI_Reduce(\$from, \$to, count, datatype, operation, root, communicator)
   MPI_Allreduce(\$from, \$to, count, datatype, operation, communicator)
   MPI_Scatter(\$from, count, type, \$to, count, type, root, communicator)
   MPI_Gather(\$from, count, type, \$to, count, type, root, communicator)


=head1 AUTHORS

Josh Wilmes and Chris Stevens
Modified to current version by Carson Holt

=head1 SEE ALSO

MPI man pages.
The paper, "Parallel::MPI - An MPI Binding for Perl", included in the
  Parallel::MPI distribution

=cut
