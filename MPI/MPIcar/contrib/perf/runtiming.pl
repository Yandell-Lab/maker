#!/usr/local/bin/perl


doit("mpi_timing");
doit("mpi_timing.pl");

sub doit {
    my ($prog) = @_;
    my @times;

    open(F,"mpirun -np 2 -machinefile ../t/machinefile $prog|")
      || die "Error: $prog, $!\n";
    
    while (<F>) {
	if (/uSec =\s+(\d+)/) {
	    push @times, $1;
	}
    }

    @times = reverse sort { $a <=> $b } @times;
    
    # skip some of the highest elements.
    shift @times;
    shift @times;
    shift @times;
    shift @times;
    shift @times;
    shift @times;
    
    my $sum;
    my $max = 0; 
    my $min = 99999999;
    foreach (@times) { 
	$sum += $_;
	$max = $_ if ($_ > $max);
	$min = $_ if ($_ < $min);
    }
    my $avg = int($sum/scalar(@times));

    printf("%-15.15s $avg uSec (avg) ($min - $max)\n","$prog:");
}


