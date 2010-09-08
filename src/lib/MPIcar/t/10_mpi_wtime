#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

my $rank = MPI_Comm_rank(MPI_COMM_WORLD);

eval {
  $start = MPI_Wtime();
};

if($rank == 0) {
    if ($@ || ! $start) {
       print "$@\n";
       print "not ok 1\n"
    } else {
       print "ok 1\n";
    }
}

sleep(1);


eval {
  $stop = MPI_Wtime();
};

print "Start=$start, Stop=$stop\n";

if($rank == 0) {
    if ($@ || ! $stop || $stop < $start) {
       print "$@\n";
       print "not ok 2\n" 
    } else {
       print "ok 2\n";
    }
}

MPI_Finalize();
