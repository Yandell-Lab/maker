test: heat_plate_mpi.pl
	echo "10000 0.2" | mpirun -np 4 -machinefile ../../t/machinefile heat_plate_mpi.pl

heat_plate_mpi.pl: heat_plate_mpi.PL
	./heat_plate_mpi.PL

clean:
	rm -f *.o *~ heat_plate_mpi.pl
