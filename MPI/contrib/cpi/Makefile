test: cpi cpi.pl
	mpirun -machinefile ../t/machinefile -np 2 ./cpi
	mpirun -machinefile ../t/machinefile -np 2 ./cpi.pl

cpi.pl: cpi.PL
	./cpi.PL

cpi: cpi.c
	mpicc cpi.c -o cpi

clean:
	rm -f *.o *~ cpi.pl cpi
