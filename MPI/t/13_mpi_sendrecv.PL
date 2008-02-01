#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/13_mpi_sendrecv.pl";
local(*OUTF);

open(OUTF, ">$file") or die "Cannot open $file for writing: $!\n";

print OUTF $Config{startperl} . " -I" . Cwd::cwd . "/blib/arch" .
				" -I" . Cwd::cwd . "/blib/lib\n";
print "Writing $file\n";
while(<DATA>) { print OUTF $_ }
close(OUTF);
chmod(0755, $file);
__END__
$|=1;
use Parallel::MPIcar qw(:all);

MPI_Init();

$my_rank = MPI_Comm_rank(MPI_COMM_WORLD);
$p = MPI_Comm_size(MPI_COMM_WORLD);

#print "pid = $$, rank = $my_rank\n";
if ($my_rank != 0) {
    @message = qw(1 2);
    my @recv;
    
    MPI_Sendrecv(\@message, 2, MPI_INT, 0, 10,
                 \@recv,    2, MPI_INT, 0, 11,
	         MPI_COMM_WORLD);
    print "[$my_rank] @recv\n";
    if($recv[0] == 3 and $recv[1] == 4) {
	print "procok $my_rank\n";
    } else {
	print "not procok $my_rank\n";
    }
} else { 
    # $my_rank == 0
    @message = qw(3 4);
    my @recv;
    
    MPI_Sendrecv(\@message, 2, MPI_INT, 1, 11,
                 \@recv,    2, MPI_INT, 1, 10,
	         MPI_COMM_WORLD);
    print "[$my_rank] @recv\n";
    if($recv[0] == 1 and $recv[1] == 2) {
	print "procok $my_rank\n";
    } else {
	print "not procok $my_rank\n";
    }
}

MPI_Finalize();

