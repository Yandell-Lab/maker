#!/usr/bin/perl
eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#Hack to get around findbin error where broken Carp is preloaded
BEGIN {
    if(@ARGV == 1 && $ARGV[0] eq 'findbin'){
        eval 'require FindBin';
        print $FindBin::RealBin;
        exit;
    }

    my $carp_version = `$^X -MCarp -e 'print \$Carp::VERSION'`;
    if($carp_version < 1.24){
        my $bin = `$0 findbin`;
        eval "use lib '$bin/../lib'";
        eval "use lib '$bin/lib'";
        eval "use lib '$bin/../perl/lib'";
        eval "use lib '$bin/inc/lib'";
        eval "use lib '$bin/inc/perl/lib'";
    }
}

use FindBin;
use lib "$FindBin::RealBin/../lib"; #libary just needed for configuration
use lib "$FindBin::RealBin/lib"; #libary just needed for configuration
use lib "$FindBin::RealBin/../perl/lib"; #libary just needed for configuration
use lib "$FindBin::RealBin/inc/lib"; #libary just needed for configuration
use lib "$FindBin::RealBin/inc/perl/lib"; #libary just needed for configuration

use strict;
use warnings;
use MAKER::Build;

chdir($FindBin::RealBin); #alway go here for safety

#configure Build
my $build = MAKER::Build->new(
     module_name => 'MAKER',
     dist_author => 'Yandell Lab <maker-devel@yandell-lab.org>',
     dist_version_from  => "$FindBin::RealBin/version",
     dist_abstract => ('MAKER is a portable and easily configurable genome annotation pipeline. '.
		       'It\'s purpose is to allow smaller eukaryotic and prokaryotic genome '.
		       'projects to independently annotate their genomes and to create genome '.
		       'databases. MAKER identifies repeats, aligns ESTs and proteins to a '.
		       'genome, produces ab-initio gene predictions and automatically '.
		       'synthesizes these data into gene annotations having evidence-based '.
		       'quality values.'
		       ),
    dynamic_config => 1,
    license  => 'artistic_2',
    auto_configure_requires => 1,
    requires => {
	'perl'               => '5.8.8',
	'File::Which'        => '0',
	'Bio::Root::Version' => '1.006901',
	'IO::All'            => '0',
	'IO::Prompt'         => '0',
	'DBD::SQLite'        => '0',
	'Perl::Unsafe::Signals' => '0',
	'Bit::Vector'        => '0',
	'Inline::C'          => '0.53',
	'forks'              => '0.34',
	'forks::shared'      => '0.34',
	'Want'               => '0',
	#'Proc::ProcessTable' => '0',
	#'Other::Module' => '>= 1.2, != 1.5, < 2.0',
    },
    exe_requires => {
	'BLAST' => 'xdformat, formatdb, makeblastdb',
	'RepeatMasker' => 'RepeatMasker',
	'Exonerate' => 'exonerate',
	'SNAP'  => 'snap',
	#'Augustus' => 'augustus',
    },
    lib_requires => {
    },
    build_requires => {
	'Test::More'         => '0',
	'ExtUtils::CBuilder' => '0',
	'DBI'                => '1.57',
	'File::Which'        => '0',
	'LWP::Simple'        => '0',
    },
    config_requires => {
	'File::Which'        => '0',
	'Term::ReadKey'      => '0',
    },
    recommends => {
	'DBD::Pg'            => '0',
	'Carp'               => '1.24',
	#'Other::Module' => '>= 1.2, != 1.5, < 2.0',
    },
    install_base => "$FindBin::RealBin/../",
    install_base_relpaths => {
	arch   => [qw(perl/lib)],
	lib    => [qw(perl/lib)],
	bindoc => [qw(perl/man)],
	libdoc => [qw(perl/man)],
    },
    PL_files => {
	'bin/TACC.PL' => ['bin/ibrun'],
    },
    get_options => {
	AnyDBM_ISA =>  { type  => '=s' },
	MWAS       =>  {},
    }
    );

#get command line AnyDBM_File override (note NDBM_File is not supported)
if(my $any = $build->args('AnyDBM_ISA')){
    $build->log_info("Configuring alternate AnyDBM_File order...");
    my @dbms = grep {$_} split(/[\,\s\t]/, $any);

    my $err;
    map {$err .= "ERROR: $_ is invalid\n" if(!/^(DB|GDBM|SDBM)_File$/)} @dbms;
    if($err){
	$build->log_info("FAILED\n");
	die $err."ERROR: Valid entries are DB_File, GDBM_File, and SDBM_File\n";
    }
    if(!grep {eval "require $_"} @dbms){
	$build->log_info("FAILED\n");
	die "ERROR: None of your entries exist on this system\n";
    }
    $build->config_data(AnyDBM_ISA => \@dbms);
    $build->log_info("DONE\n");
    $build->log_info("New order=".join(',', @dbms)."\n");
}

#ask for MPI
my $go = $build->y_n("\nMAKER supports distributed parallelization via MPI.\n".
		     "Would you like to configure MAKER for MPI (This\n".
		     "requires that you have an MPI client installed)?", 'N');
$build->feature(mpi_support => $go);
$build->config_mpi() if($go);


#configure for MWAS
if($build->args('MWAS')){
    die "ERROR: You must have apache installed and configured to use the MWAS interface\n"
	unless (File::Which::which('httpd') || -f '/etc/init.d/httpd' || -f '/etc/init.d/apache2');
    $build->feature(mwas_support => 1); #set mwas flag for Build
}
else{
    $build->feature(mwas_support => 0); #set mwas flag for Build
}

if($build->feature('mwas_support')){
    #add required CGI dependencies
    $build->add_build_requires('Test::WWW::Mechanize' => '0');
    $build->add_build_requires('Data::Random' => '0');
    $build->add_build_requires('Date::Parse' => '0');
    $build->add_build_requires('Net::OpenID::Consumer' => '0');
    $build->add_build_requires('Archive::Zip' => '0');
    $build->add_requires('CGI::Application' => '0');
    $build->add_requires('CGI::Application::Plugin::Session' => '0');
    $build->add_requires('CGI::Application::Plugin::Authentication' => '0');
    $build->add_requires('CGI::Application::Plugin::CAPTCHA' => '0');
    $build->add_requires('CGI::Application::Plugin::DBH' => '0');
    $build->add_requires('CGI::Application::Plugin::TT' => '0');
    $build->add_requires('CGI::Application::Plugin::Redirect' => '0');
    $build->add_requires('CGI::Application::Plugin::Stream' => '0');
    $build->add_requires('Mail::Sender' => '0');
    $build->add_requires('Bio::Graphics::Browser2' => '2.15');
    $build->add_requires('XML::Twig' => '0');
    $build->add_requires('GD' => '0');
    $build->add_requires('Devel::Size' => '0');
    $build->add_requires('Heap::Simple' => '0');
    $build->add_requires('Heap::Simple::XS' => '0');
    $build->add_exe_requires('Apollo' => 'apollo');
    $build->add_exe_requires('Apache Ant' => 'ant');
    $build->add_exe_requires('JBrowse' => 'flatfile-to-json.pl');
    $build->add_lib_requires('LibPNG' => 'png');
    $build->add_lib_requires('LibGD' => 'gd');

    $build->set_gmod_locs();
}

$build->config_data(PATH => $ENV{PATH});
$build->config_data(build_done => 0);
$build->create_build_script;

$build->log_info("\n\nThe file 'Build' has been created for you to finish installing MAKER.\n");

$build->maker_status;
