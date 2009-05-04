#! /usr/bin/perl -w 

use strict;
use FindBin;

chdir "$FindBin::Bin/MPIcar";
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("MPIcar/Makefile.old") if(-e "MPIcar/Makefile.old");

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/");
system("cp $FindBin::Bin/mpi_evaluator $FindBin::Bin/../bin/");
