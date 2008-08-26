#!/bin/bin/perl

$|=1;

print "1..3\n";

$ret = system("mpirun -machinefile t/machinefile -np 2 t/10_mpi_wtime.pl");

if ($ret) {
   print "$@\n";
   print "not ok 3\n"
} else {
   print "ok 3\n";
}
