#!/usr/bin/perl -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/arch -I/uufs/inscc.utah.edu/common/home/u0045039/Users/cholt/Desktop/Parallel-MPI-0.03/blib/lib
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

$my_rank = MPI_Comm_rank(MPI_COMM_WORLD);
$p = MPI_Comm_size(MPI_COMM_WORLD);

#print "pid = $$, rank = $my_rank\n";
$tag = 0;
if ($my_rank != 0) {
    @message = (31337, 666);
    $dest = 0;
    
    MPI_Send(\@message, 5, MPI_INT, $dest, $tag, MPI_COMM_WORLD);
} else { 
    # my_rank == 0
    $source = 1;

    @status = MPI_Recv(\@message, 5, MPI_INT, $source, $tag, MPI_COMM_WORLD);	
	
    print("Recieved: \"@message\" from $source\n");
    print("Status: (", (join ', ',@status), ")\n");
    if($message[0] == 31337 && $message[1] == 666) {
	print "ok 1\n";
    } else {
	print "MESSAGE: @message\nnot ok 1\n";
    }
    # (count,MPI_SOURCE,MPI_TAG,MPI_ERROR)
}

MPI_Finalize();

