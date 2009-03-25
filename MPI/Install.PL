#! /usr/bin/perl -w 

use strict;
use FindBin;

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

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/");
