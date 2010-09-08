#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
my $from = $rank;

MPI_Reduce(\$from, \$to, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

if ($rank == 0) {
     print "Reduction produced $to\n";    
     if ($to == 1) {
        print "ok 1\n";
     } else {
        print "not ok 1\n";     
     }
}

MPI_Finalize();

