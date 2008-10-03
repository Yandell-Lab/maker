#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
my $x = $rank;
my $y;

MPI_Allreduce(\$x, \$y, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

print "rank=$rank, x=$x, y=$y\n";
if ($y == 6) {
  print "procok $rank\n";
} else {
  print "not procok $rank\n";
}

MPI_Finalize();
