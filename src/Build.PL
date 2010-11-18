#!/usr/bin/perl
eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/inc/lib"; #libary just needed for configuration
use MAKER::Build;

#configure Build
my $build = MAKER::Build->new
    (
     module_name => 'MAKER',
     dist_author => 'Yandell Lab <maker-devel@yandell-lab.org>',
     dist_version_from  => 'version',
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
	 'perl'               => '5.8.0',
	 'Proc::ProcessTable' => '0',
	 'File::Which'        => '0',
	 'Bio::Root::Version' => '1.006',
	 'IO::All'            => '0',
	 'IO::Prompt'         => '0',
	 'DBD::SQLite'        => '0',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     exe_requires => {
	 'BLAST' => 'xdformat, formatdb, makeblastdb',
	 'RepeatMasker' => 'RepeatMasker',
	 'Exonerate' => 'exonerate',
	 'SNAP'  => 'snap',
	 'Augustus' => 'augustus',
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
     },
     recommends => {
	 'threads'  => '1.67',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     install_base => "$FindBin::Bin/../",
     install_base_relpaths => {
	 arch   => [qw(perl/lib)],
	 lib    => [qw(perl/lib)],
	 bindoc => [qw(perl/man)],
	 libdoc => [qw(perl/man)],
     },
     PL_files => {'lib/Parallel/MPIcar.PL' => ['lib/Parallel/MPIcar.pm',
						'lib/Parallel/MPIcar.xs',
						'lib/Parallel/utils.c',
						'lib/Parallel/typemap',],
		  'bin/mpi_scripts.PL' => ['bin/mpi_maker',
					   'bin/mpi_evaluator',
					   'bin/mpi_iprscan',],
     },
    );

#see if threaded and supports MPI
if(! $build->thread_support){
    print "\nWARNING:  Your version of perl does not seem to have been compiled\n".
	  "for threads; certain MAKER accessory scripts and features (such as\n".
	  "MPI support) will be not be available to you.\n\n";

    $build->feature(mpi_support => 0); #set mpi flag for Build
}
elsif(my @mpi = $build->mpi_support){
    my $go = $build->y_n("\nYour system appears to be configured for MPICH2. Would\n".
			 "you like to install the MPI version of MAKER?", 'Y');

    if($go){
	my $cc = shift @mpi;
	my $mpidir = shift @mpi;
	$build->feature(mpi_support => 1); #set mpi flag for Build
	$build->extra_compiler_flags("-I$mpidir -DFLOAT_HACK"); #add mpidir to includes
	$build->config('cc' => 'mpicc'); #mpicc now default compiler
	$build->config('ld' => 'mpicc');
	$build->add_requires('threads' => '1.67');	#add threads to require
	#$build->>add_to_cleanup();
    }
    else{
	$build->feature(mpi_support => 0); #set mpi flag for Build
    }
}
else{
    $build->feature(mpi_support => 0); #set mpi flag for Build
}

if(File::Which::which('httpd') || -f '/etc/init.d/httpd' || -f '/etc/init.d/apache2'){
    my $go = 0;
#    my $go = $build->y_n("\nWould you like to install the web interface to MAKER (MWAS) on this machine\n".
#			 "This requires that you be logged in as root or use sudo for \'./Build install\'?", 'Y');
    
    if($go){
	$build->feature(mwas_support => 1); #set mwas flag for Build

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
	$build->add_requires('Mail::Sender' => '0');
	$build->add_requires('Bio::Graphics::Browser2' => '2.15');
	$build->add_exe_requires('Apollo' => 'apollo');
	$build->add_exe_requires('JBrowse' => 'flatfile-to-json.pl');
    }
    else{
	$build->feature(mwas_support => 0); #set mpi flag for Build
    }
}

$build->create_build_script;

print "\n\nThe file 'Build' has been created for you to finish installing MAKER.\n";

$build->maker_status;
