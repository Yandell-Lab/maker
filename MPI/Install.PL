#! /usr/bin/perl -w 

use strict;
use FindBin;
use lib"$FindBin::Bin/../lib";
use lib"$FindBin::Bin/../perl/lib";

chdir "$FindBin::Bin/MPIcar";
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("MPIcar/Makefile.old") if(-e "MPIcar/Makefile.old");

eval{
    require Parallel::MPIcar;
};

die "ERROR: MPI setup failed for maker\n" if $@;

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/");
system("cp $FindBin::Bin/mpi_evaluator $FindBin::Bin/../bin/");
system("cp $FindBin::Bin/mpi_iprscan $FindBin::Bin/../bin/");

print STDERR "\n\nMPI setup for maker was successful!!\n\n";
