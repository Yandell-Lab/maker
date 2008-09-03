#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/06_mpi_sendrecv.pl";
local(*OUTF);

open(OUTF, ">$file") or die "Cannot open $file for writing: $!\n";

print OUTF $Config{startperl} . " -I" . Cwd::cwd . "/blib/arch" .
				" -I" . Cwd::cwd . "/blib/lib\n";
print "Writing $file\n";
while(<DATA>) { print OUTF $_ }
close(OUTF);
chmod(0755, $file);
__END__
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

$my_rank = MPI_Comm_rank(MPI_COMM_WORLD);
$p = MPI_Comm_size(MPI_COMM_WORLD);

#print "pid = $$, rank = $my_rank\n";
$tag = 0;
if ($my_rank != 0) {
    $message = sprintf("Greetings from process %d!", $my_rank);
    $dest = 0;
    
    MPI_Send(\$message, length($message), MPI_CHAR, $dest, $tag,
	     MPI_COMM_WORLD);
} else { 
    # my_rank == 0
    for $source (1..$p-1) {
	$s = 0;
	@status = MPI_Recv(\$message, 100, MPI_CHAR, $source, $tag, 
			   MPI_COMM_WORLD);	
	
	printf("Recieved: \"%s\" from $source\n", $message);
	printf("Status: (" . (join ', ',@status) . ")\n");
	# (count,MPI_SOURCE,MPI_TAG,MPI_ERROR)
    }
}

MPI_Finalize();

