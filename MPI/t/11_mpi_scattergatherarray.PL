#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/11_mpi_scattergatherarray.pl";
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

my @message = (0.31337, 0.666, 0.42448, 0.777);

my @sendbuf = @message;
my @recvbuf;

print "RANK $my_rank: SCATTER\n";

MPI_Scatter(\@sendbuf, 2, MPI_FLOAT,
	    \@recvbuf, 2, MPI_FLOAT, 0, MPI_COMM_WORLD);

print "RANK $my_rank:  RECVBUF=\"@recvbuf\"\n";

@sendbuf = ();

MPI_Gather(\@recvbuf, 2, MPI_FLOAT,
	   \@sendbuf, 2, MPI_FLOAT, 0, MPI_COMM_WORLD);

if ($my_rank == 0) {
    my $flag;
    printf "FINAL BUF = @sendbuf\n";
    foreach (0..$#sendbuf) {
	if($sendbuf[$_] != $message[$_]) {
    	    $flag = 1; last;
	}
    }
    unless ($flag) {
    	print "ok 1\n";
    } else {
	print "not ok 1\n";
    }
}

MPI_Finalize();

