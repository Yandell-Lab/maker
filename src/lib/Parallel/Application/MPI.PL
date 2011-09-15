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
    File::Copy::copy("inc/src/Parallel/Application/MPI.pm", "$FindBin::Bin/");
}
else{ #clean up if they already exist
    unlink("$FindBin::Bin/MPI.pm");
    my @rm = grep {!/\.PL$/} <$FindBin::Bin/MPI.*>;
    map {unlink($_)} @rm;
}
