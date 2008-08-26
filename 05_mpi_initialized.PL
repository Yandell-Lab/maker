#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();  

my $rank = MPI_Comm_rank(MPI_COMM_WORLD);

eval {
  $init = MPI_Initialized();
};

if ($@ || ! $init) {
   print "$@\n";
   print "not procok $rank\n" 
} else {
   print "procok $rank\n";
}

MPI_Finalize();
