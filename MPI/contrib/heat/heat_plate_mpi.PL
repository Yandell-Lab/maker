#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "heat_plate_mpi.pl";
local(*OUTF);

open(OUTF, ">$file") or die "Cannot open $file for writing: $!\n";

print OUTF $Config{startperl}, "\n\n";
print OUTF "use lib qw(", Cwd::cwd, "/../../blib/arch ",
                          Cwd::cwd, "/../../blib/lib);\n\n";

print "Writing $file\n";
while(<DATA>) { print OUTF $_ }
close(OUTF);
chmod(0755, $file);
__END__
$| = 1;

use Parallel::MPIcar qw(:all);

use vars qw($t $s $r $t1 $t2 $h $xl $yl $tol $h2 $aa $sum $global_sum
	    $nx $ny $i $j $ij $nproc $iter $last_iter $niter $n_left $npt $nxy
	    $nproc $my_rank $root @istatus @local_n @local_a);

MPI_Init();

$nproc   = MPI_Comm_size(MPI_COMM_WORLD);
$my_rank = MPI_Comm_rank(MPI_COMM_WORLD);

$xl = 10.0;
$yl = 10.0;
$tol = 0.001;

# initial set up done by my_rank = 0 (i.e. master) process only

if($my_rank eq 0) {
    print "there are $nproc processes\n";
    print "enter number of iterations and h:";
    ($niter, $h) = split(' ', scalar(<STDIN>));

    #  determine the workload for each processor
    $nx = $xl/$h + 1;
    $ny = $nx;
    print "nx=ny=$nx\n";

    foreach(0..($nproc-1)) {
	use integer;
	$local_n[$_] = $nx/$nproc;
    }

    # if workload cannot be divided evenly
    # to each process, assign to the
    # first n_left processes with an
    # additional piece of work

    if(($nx % $nproc) != 0) {
	use integer;
	$n_left = $nx - $nproc*($nx/$nproc);
	foreach(0..($n_left-1)) {
	    $local_n[$_]++;
	}

	$local_a[0] = 0.0;
	foreach(1..($nproc-1)) {
	      $local_a[$_] = $local_a[$_-1] + $h*$local_n[i-1]
	}
    }

    print "process #,  local workload,  start location\n";
    foreach(0..($nproc-1)) {
	print $_, "\t", $local_n[$_], "\t", $local_a[$_], "\n";
    }
}

#  now all the processes begin to work together

MPI_Bcast(\$h,    1,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(\$niter,1,MPI_INT,   0,MPI_COMM_WORLD);
MPI_Scatter(\@local_n,1,MPI_INT,   \$npt,1,MPI_INT,   0,MPI_COMM_WORLD);
MPI_Scatter(\@local_a,1,MPI_DOUBLE,\$aa, 1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#  add one additional point to the end strips and 
#  add additional 2 points to the interior strips 
#  to facilitate overlap      

$yl = 10.0;
$ny = $yl/$h + 1;

if($my_rank == 0 or $my_rank == ($nproc-1)) {
    $nx = $npt + 1;
} else {
    $nx = $npt + 2;
}

print "I am process $my_rank  nx,ny= $nx,$ny\n";
print "[$my_rank] aa=$aa\n";
print "[$my_rank] h=$h\n";
print "[$my_rank] niter=$niter\n";

$h2  = $h*$h;
$nxy = $nx*$ny;

#   source term

foreach $j (0..($ny-1)) {
    foreach $i (0..($nx-1)) {
	$t[$i][$j] = 20.0;
	$s[$i][$j] = 0.0;
	$r[$i][$j] = 0.0;
    }
}

#  fix the boundary conditions
#
#  domain is divided into horizontal stripes

#    along the west and east walls

#foreach $j (0..($ny-1)) {
#    $t[0][$j] = 20.0;
#    $t[$nx-1][$j] = 20.0;
#}
 
#   along the south and north walls

#  south boundary is owned by the first processor

if($my_rank == 0) {
    foreach $i (0..($nx-1)) {
	$t[$i][1] = 20.0;
    }
}

#  north boundary is owned by the last processor
 
if($my_rank == ($nproc-1)) {
    foreach $i (0..($ny-1)) {
	$xx = $i * $h;
	if($xx >= 3.0 and $xx <= 7.0) {
	    $t[$nx-1][$i] = 100.0;
        } else {
            $t[$nx-1][$i] = 20.0;
	}
    }
}

MPI_Barrier(MPI_COMM_WORLD);

if($my_rank == 0) {
    print "finish updating boundary conditions\n";
    $t1 = MPI_Wtime()
}


foreach $iter (1..$niter) {
    # send and receive interfacial values from nearest neighbours
    if($my_rank != 0 and $my_rank != ($nproc-1)) {
	# along the south side
        @status = MPI_Sendrecv($t[1],$ny,MPI_DOUBLE,$my_rank-1,10,
                               $t[0],$ny,MPI_DOUBLE,$my_rank-1,20,
                               MPI_COMM_WORLD);

	# along the north side
	@status = MPI_Sendrecv($t[$nx-2],$ny,MPI_DOUBLE,$my_rank+1,20,
                               $t[$nx-1],$ny,MPI_DOUBLE,$my_rank+1,10,
			       MPI_COMM_WORLD);

    }

    if($my_rank == 0) {
	# along the north side
        $root = 1;
        @status = MPI_Sendrecv($t[$nx-2],$ny,MPI_DOUBLE,$root,20,
                               $t[$nx-1],$ny,MPI_DOUBLE,$root,10,
			       MPI_COMM_WORLD);

    }


    if($my_rank == ($nproc-1)) {
	# along the south side
	$root = $nproc-2;
        @status = MPI_Sendrecv($t[1],$ny,MPI_DOUBLE,$root,10,
                               $t[0],$ny,MPI_DOUBLE,$root,20,
			       MPI_COMM_WORLD);
    }

    $sum = 0.0;

    foreach $j (1..($ny-2)) {
	foreach $i (1..($nx-2)) {
	    $r[$i][$j] = $s[$i][$j] * $h2 - (
                              $t[$i+1][$j] + $t[$i][$j+1]
                              - 4.0 * $t[$i][$j] +
                              $t[$i-1][$j] + $t[$i][$j-1]);
	    $sum = $sum + $r[$i][$j] * $r[$i][$j];
	}
    } 

    #  update temperatures
    foreach $j (0..($ny-1)) {
	foreach $i (0..($nx-1)) {
	    $t[$i][$j] = $t[$i][$j] - 0.25 * $r[$i][$j];
	}
    }

    # obtain global sum
    MPI_Allreduce(\$sum,\$global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    $global_sum = sqrt($global_sum);

    if($my_rank == 0) {
	if(($iter % 50) == 0) {
	    print "iter=$iter\nglobal_sum=$global_sum\n";
	}
    }

    $last_iter = $iter;
    last if $global_sum <= $tol;
}

if($my_rank == 0) {
    $t2 = MPI_Wtime();
    print "no. of iterations took = $last_iter\n";
    print "final L2 norm of residual: $global_sum\n";
    print "wall clock time in seconds= ",($t2-$t1), "\n";
}

#sleep($my_rank);
#if($my_rank == 0) {
#    foreach $i (0..$nx-2) {
#	print "[$my_rank]";
#	foreach $j (0..$ny-1) {
#	    printf "%3d ", $t[$i][$j];
#	}
#	print "\n";
#    }
#} elsif($my_rank == $nproc-1) {
#    foreach $i (1..$nx-1) {
#	print "[$my_rank]";
#	foreach $j (0..$ny-1) {
#	    printf "%3d ", $t[$i][$j];
#	}
#	print "\n";
#    }
#} else {
#    foreach $i (1..$nx-2) {
#	print "[$my_rank]";
#	foreach $j (0..$ny-1) {
#	    printf "%3d ", $t[$i][$j];
#	}
#	print "\n";
#    }
#}
#print "\n" if $my_rank == 1;

MPI_Finalize();
