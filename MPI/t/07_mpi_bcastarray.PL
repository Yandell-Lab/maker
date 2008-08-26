#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/07_mpi_bcastarray.pl";
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

my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
my @x = (($rank+0.31337), 0.666);

MPI_Bcast(\@x, 2, MPI_FLOAT, 1, MPI_COMM_WORLD);

if ($x[0] == 1.31337 and $x[1] == 0.666) {
  print "procok $rank\n";
} else {
  print "x=@x\n";
  print "not procok $rank\n";
}

MPI_Finalize();
