#! /usr/bin/perl -w 

use strict;
use FindBin;

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

chdir ("$FindBin::Bin/Proc-ProcessTable-0.41");
system("perl Makefile.PL");
system("make install");
system("make clean");
chdir("$FindBin::Bin");
unlink("Proc-ProcessTable-0.41/Makefile.old") if (-e "Proc-ProcessTable-0.41/Makefile.old");
unlink("Proc-ProcessTable-0.41/Process/Makefile.old") if (-e "Proc-ProcessTable-0.41/Process/Makefile.old");
