#! /usr/bin/perl

use warnings;
use strict;
use Cwd;
use lib Cwd::cwd()."/inc/lib"; #this PL file always run by Build, so cwd should be constant
use MAKER::Build;
use FindBin;
use File::Copy;

my $build = MAKER::Build->resume (
    properties => {
	config_dir => '_build',
    },
    );

if($build->feature('mpi_support')){
    File::Copy::copy("inc/bin/mpi_maker", "$FindBin::Bin/");
    File::Copy::copy("inc/bin/mpi_evaluator", "$FindBin::Bin/");
    File::Copy::copy("inc/bin/mpi_iprscan", "$FindBin::Bin/");
}
else{ #clean up if they already exist
    unlink("$FindBin::Bin/mpi_maker");
    unlink("$FindBin::Bin/mpi_evaluator");
    unlink("$FindBin::Bin/mpi_iprscan");
}
