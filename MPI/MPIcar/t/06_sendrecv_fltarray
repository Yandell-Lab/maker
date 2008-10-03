#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/06_sendrecv_fltarray.pl";
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
$tag = 0;
if ($my_rank != 0) {
    @message = qw(0.31337 0.666);
    $dest = 0;
    
    MPI_Send(\@message, 2, MPI_FLOAT, $dest, $tag, MPI_COMM_WORLD);
} else { 
    # my_rank == 0
    my $source = 1;
    @status = MPI_Recv(\@message, 2, MPI_FLOAT, $source, $tag, MPI_COMM_WORLD);

    print("Recieved: \"@message\" from $source\n");
    print("Status: (", (join ', ',@status), ")\n");
    if($message[0] == 0.31337 && $message[1] == 0.666) {
	print "ok 1\n";
    } else {
	print "MESSAGE: @message\nnot ok 1\n";
    }
    # (count,MPI_SOURCE,MPI_TAG,MPI_ERROR)
}

MPI_Finalize();

