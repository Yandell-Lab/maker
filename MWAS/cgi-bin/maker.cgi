#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../../../maker/lib";
use lib "$FindBin::Bin/../../../maker/perl/lib";
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/perl/lib";
         
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use MWS;

my $webapp = MWS->new(die_on_bad_params => 0,
		      cache             => 0,
		     );

$webapp->run();

