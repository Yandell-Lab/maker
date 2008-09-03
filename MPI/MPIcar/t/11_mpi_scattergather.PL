#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();
$my_rank = MPI_Comm_rank(MPI_COMM_WORLD);

my $string = "12345678";

my $sendbuf = $string;
my $recvbuf;

print "RANK $my_rank: SCATTER\n";

MPI_Scatter(\$sendbuf, 4, MPI_CHAR,
	    \$recvbuf, 4, MPI_CHAR, 0, MPI_COMM_WORLD);

print "RANK $my_rank:  RECVBUF=\"$recvbuf\"\n";

$sendbuf = "scribble";

MPI_Gather(\$recvbuf, 4, MPI_CHAR,
	   \$sendbuf, 4, MPI_CHAR, 0, MPI_COMM_WORLD);

if ($my_rank == 0) {
    printf "FINAL BUF = $sendbuf\n";
    if ($sendbuf eq $string) {
    	print "ok 1\n";
    } else {
    	print "not ok 1\n";
    }
}

MPI_Finalize();

