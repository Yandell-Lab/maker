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
    File::Copy::copy("inc/src/Parallel/MPIcar.pm", "$FindBin::Bin/");
    File::Copy::copy("inc/src/Parallel/MPIcar.xs", "$FindBin::Bin/");
    File::Copy::copy("inc/src/Parallel/typemap", "$FindBin::Bin/");
    File::Copy::copy("inc/src/Parallel/utils.c", "$FindBin::Bin/");
}
else{ #clean up if they already exist
    unlink("$FindBin::Bin/MPIcar.pm");
    unlink("$FindBin::Bin/MPIcar.xs");
    unlink("$FindBin::Bin/typemap");
    unlink("$FindBin::Bin/utils.c");
    my @rm = grep {!/\.PL$/} <$FindBin::Bin/MPIcar.*>;
    map {unlink($_)} @rm;
}
