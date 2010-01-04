#! /usr/bin/perl -w 

use strict;
use FindBin;
use lib"$FindBin::Bin/../lib";
use lib"$FindBin::Bin/../perl/lib";

chdir ("$FindBin::Bin/threads-1.71");
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("threads-1.71/Makefile.old") if(-e "threads-1.71/Makefile.old");

chdir "$FindBin::Bin/MPIcar";
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("MPIcar/Makefile.old") if(-e "MPIcar/Makefile.old");

eval{
    require threads;
};

die "ERROR: MPI setup failed for maker, could not install threads\n" if $@;

eval{
    require Parallel::MPIcar;
};

die "ERROR: MPI setup failed for maker, could not install Parallel::MPIcar\n" if $@;

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/");
system("cp $FindBin::Bin/mpi_evaluator $FindBin::Bin/../bin/");
system("cp $FindBin::Bin/mpi_iprscan $FindBin::Bin/../bin/");

print STDERR "\n\nMPI setup for maker was successful!!\n\n";
