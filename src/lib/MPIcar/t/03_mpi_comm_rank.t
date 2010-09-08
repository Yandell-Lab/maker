#!/bin/bin/perl

$|=1;

print "1..2\n";

open (F, "mpirun -machinefile t/machinefile -np 2 t/03_mpi_comm_rank.pl|")
    or die $!;

while (<F>) {
        print;
        if (/^procok (\d+)/) {
            $procok{$1}++;
        }
}

for (0..1) {
   if ($procok{$_}) {
     print "ok " . ($_+1) . "\n";
   } else {
     print "not ok " . ($_+1) . "\n";   
   }
}
