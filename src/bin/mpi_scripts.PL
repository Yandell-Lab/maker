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
    File::Copy::copy("inc/bin/mpi_evaluator", "$FindBin::RealBin/");
    File::Copy::copy("inc/bin/mpi_iprscan", "$FindBin::RealBin/");
}
else{ #clean up if they already exist
    unlink("$FindBin::RealBin/mpi_maker");
    unlink("$FindBin::RealBin/mpi_evaluator");
    unlink("$FindBin::RealBin/mpi_iprscan");
}
