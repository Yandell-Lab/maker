#! /usr/bin/perl -w 

use strict;
use FindBin;

chdir ("$FindBin::Bin/threads-1.71");
system("perl Makefile.PL");
system("make install");
chdir("$FindBin::Bin");

chdir "$FindBin::Bin/MPIcar";
system("perl Makefile.PL");
system("make install");
chdir("$FindBin::Bin");

system("cp $FindBin::Bin/mpi_maker $FindBin::Bin/../bin/");
