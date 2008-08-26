#!/usr/bin/perl -w

use strict;
use Config;
use Cwd;

my $file = "t/01_mpi_init.pl";
local(*OUTF);

open(OUTF, ">$file") or die "Cannot open $file for writing: $!\n";

print OUTF $Config{startperl}, "\n\n";
print OUTF "use lib qw(", Cwd::cwd, "/blib/arch ",
	                  Cwd::cwd, "/blib/lib);\n\n";

print "Writing $file\n";
while(<DATA>) { print OUTF $_ }
close(OUTF);
chmod(0755, $file);
__END__
$|=1;
use Parallel::MPIcar qw(:all);

eval {
  MPI_Init();  
};

if ($@) {
    my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
    print "$@\n";
    print "not procok $rank\n" 
} else {
    my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
    print "procok $rank\n";
}

MPI_Finalize();
