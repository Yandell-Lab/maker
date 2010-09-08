#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "cpi.pl";
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
$|=1;
use Parallel::MPIcar qw(:all);

sub f {
    my ($a) = @_;
    return (4.0 / (1.0 + $a * $a));
}

my $PI25DT = 3.141592653589793238462643;

MPI_Init();
$numprocs = MPI_Comm_size(MPI_COMM_WORLD);
$myid =     MPI_Comm_rank(MPI_COMM_WORLD);

#printf(STDERR "Process %d\n", $myid);

$n = 0;
while (1) {
    if ($myid == 0) {
	if ($n==0) { $n=100; } else { $n=0; }
	$startwtime = MPI_Wtime();
    }

    MPI_Bcast(\$n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    last if ($n == 0);

    $h   = 1.0 / $n;
    $sum = 0.0;
    
    for ($i = $myid + 1; $i <= $n; $i += $numprocs) {
	$x = $h * ($i - 0.5);
	$sum += f($x);
    }
    $mypi = $h * $sum;

    MPI_Reduce(\$mypi, \$pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
    if ($myid == 0) {
	printf("pi is approximately %.16f, Error is %.16f\n",
	       $pi, abs($pi - $PI25DT));
	$endwtime = MPI_Wtime();
	printf("wall clock time = %f\n", $endwtime - $startwtime);	       
    }
}

MPI_Finalize();

