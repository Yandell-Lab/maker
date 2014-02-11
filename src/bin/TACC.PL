#! /usr/bin/perl

use warnings;
use strict;
use Cwd;
use lib Cwd::cwd()."/inc/lib"; #this PL file always run by Build, so cwd should be constant
use MAKER::Build;
use FindBin;
use File::Copy;
use Sys::Hostname;

my $build = MAKER::Build->resume (
    properties => {
	config_dir => '_build',
    },
    );

#see if running on TACC (replaces default ibrun command)
if(hostname() =~ /tacc\.utexas\.edu/){
    File::Copy::copy("inc/bin/ibrun", "$FindBin::RealBin/");
}
else{ #clean up if they already exist
    unlink("$FindBin::RealBin/ibrun");
}
