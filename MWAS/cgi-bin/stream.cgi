#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::RealBin/../../../maker/lib";
use lib "$FindBin::RealBin/../../../maker/perl/lib";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/perl/lib";
         
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use stream;

my $webapp = stream->new(die_on_bad_params => 0,
			 cache             => 0,
			 );

$webapp->run();

