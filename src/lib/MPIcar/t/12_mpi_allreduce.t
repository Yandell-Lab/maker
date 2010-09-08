#!/bin/bin/perl

$|=1;

print "1..4\n";

open (F, "mpirun -machinefile t/machinefile -np 4 t/12_mpi_allreduce.pl|") or die;

while (<F>) {
	print;
	if (/^procok (\d+)/) {
	    $procok{$1}++;
	}
}

for (0..3) {
   if ($procok{$_}) {
     print "ok " . ($_+1) . "\n";
   } else {
     print "not ok " . ($_+1) . "\n";   
   }
}
