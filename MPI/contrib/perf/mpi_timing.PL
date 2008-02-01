#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "mpi_timing.pl";
local(*OUTF);

open(OUTF, ">$file") or die "Cannot open $file for writing: $!\n";

print OUTF $Config{startperl}, "\n\n";
print OUTF "use lib qw(", Cwd::cwd, "/../../blib/arch ",
	                  Cwd::cwd, "/../../blib/lib);\n\n";

print "Writing $file\n";
while(<DATA>) { print OUTF $_ }
close(OUTF);
chmod(0755, $file);
__END__

# /***************************************************************************
#  * MPI Timing Program - Perl Version
#  * CONVERTED TO Perl: Josh Wilmes
#  ***************************************************************************/

$|=1;
use Parallel::MPIcar qw(:all);

$NUMBER_REPS  =  200;
$MESSAGE_SIZE =  128;

MPI_Init();

$numtasks = MPI_Comm_size(MPI_COMM_WORLD);
$taskid = MPI_Comm_rank(MPI_COMM_WORLD);

die "Error; Run with num procs = 2\n" unless ($numtasks == 2);

$reps = $NUMBER_REPS;
$type = 1;
$outmsg = "x" x $MESSAGE_SIZE;
$inmsg = "x" x $MESSAGE_SIZE;

if ($taskid == 0) {
    # round-trip timing test 
    printf("Doing round trip test, minimal message size, %d reps.\n",$reps);
    $dest = 1;
    $source = 1;

    for $n (1..$reps) {
        MPI_Barrier(MPI_COMM_WORLD);
	$starttime = MPI_Wtime();
	
	# send message to worker - message type set to 1. 
	MPI_Send(\$outmsg, $MESSAGE_SIZE, MPI_CHAR, $dest, $type, MPI_COMM_WORLD);
	
	# Now wait to receive the echo reply from the worker 
	MPI_Recv(\$inmsg, $MESSAGE_SIZE, MPI_CHAR, $source, $type, MPI_COMM_WORLD);
	
	# calculate round trip time and print 
	$endtime = MPI_Wtime();
	$dt1 = ($endtime - $starttime) * 1000000;
	printf("round trip# %2d   uSec = %8d\n", $n, $dt1);
	
	$at1 += $dt1;
    }
    printf("\n*** Round Trip Avg uSec = %d\n", $at1 / $reps);
} elsif ($taskid == 1) {
    $dest = 0;
    $source = 0;
    for $n (1..$reps) {
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(\$inmsg, $MESSAGE_SIZE, MPI_CHAR, $source, $type, MPI_COMM_WORLD);
	MPI_Send(\$outmsg, $MESSAGE_SIZE, MPI_CHAR, $dest, $type, MPI_COMM_WORLD);
    }
}

MPI_Finalize();
exit(0);

