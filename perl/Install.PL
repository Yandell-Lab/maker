#! /usr/bin/perl -w 

use strict;
use FindBin;

chdir ("$FindBin::Bin/threads-1.71");
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("threads-1.71/Makefile.old") if(-e "threads-1.71/Makefile.old");

chdir "$FindBin::Bin/DBI-1.607";
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("DBI-1.607/Makefile.old") if (-e "DBI-1.607/Makefile.old");

chdir ("$FindBin::Bin/DBD-SQLite-0.31");
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("DBD-SQLite-0.31/Makefile.old") if (-e "DBD-SQLite-0.31/Makefile.old");


