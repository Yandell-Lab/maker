test: mpi_timing mpi_timing.pl
	./runtiming.pl

mpi_timing.pl: mpi_timing.PL
	./mpi_timing.PL

mpi_timing: mpi_timing.c
	mpicc mpi_timing.c -o mpi_timing

clean:
	rm -f *.o *~ mpi_timing.pl mpi_timing
