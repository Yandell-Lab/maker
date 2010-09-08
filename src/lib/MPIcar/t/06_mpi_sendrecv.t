#!/bin/bin/perl

$|=1;

print "1..1\n";

$ret = system("mpirun -machinefile t/machinefile -np 2 t/06_mpi_sendrecv.pl");

if ($ret) {
   print "$@\n";
   print "not ok 1\n" 
} else {
   print "ok 1\n";
}
