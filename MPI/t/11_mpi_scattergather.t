#!/usr/bin/perl

$|=1;

print "1..2\n";

$ret = system("mpirun -machinefile t/machinefile -np 2 t/11_mpi_scattergather.pl");

if ($ret) {
   print "$@\n";
   print "not ok 2\n" 
} else {
   print "ok 2\n";
}
