#! /usr/bin/env perl -w 

use strict;
use FindBin;
use lib"$FindBin::Bin/../lib";
use lib"$FindBin::Bin/../perl/lib";

eval{
    require threads;
};

die "ERROR: MPI setup failed for maker, could not load threads\n".
    "Are you using perl compiled for threads?\n\n"if $@;

unless($threads::VERSION >= 1.67){
   die "mpi_maker requires threads version 1.67 or greater.\n",
   "You have version ". $threads::VERSION ."\n",
   "Please update threads to the newest version.\n",
   "Then try installing mpi_maker again\n";
}

chdir "$FindBin::Bin/MPIcar";
system("$^X Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("MPIcar/Makefile.old") if(-e "MPIcar/Makefile.old");

eval{
    require Parallel::MPIcar;
};

die "ERROR: MPI setup failed for maker, could not install Parallel::MPIcar\n" if $@;

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/") && die "ERROR: Could not copy MAKER MPI scripts";
system("cp $FindBin::Bin/mpi_evaluator $FindBin::Bin/../bin/") && die "ERROR: Could not copy MAKER MPI scripts";;
system("cp $FindBin::Bin/mpi_iprscan $FindBin::Bin/../bin/") && die "ERROR: Could not copy MAKER MPI scripts";;

print STDERR "\n\nMPI setup for maker was successful!!\n\n";
