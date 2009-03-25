#! /usr/bin/perl -w 

use strict;
use FindBin;

chdir "$FindBin::Bin/DBI-1.607";
system("perl Makefile.PL");
system("make install");
chdir("$FindBin::Bin");

chdir ("$FindBin::Bin/DBD-SQLite-0.31");
system("perl Makefile.PL");
system("make install");
chdir("$FindBin::Bin");

