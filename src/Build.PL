#! /usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/"; #libary just needed for configuration
use lib "$FindBin::Bin/../perl/lib/"; #libary just needed for configuration
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
	 'DBI'                => '0',
	 'DBD::SQLite'        => '0',
	 'Proc::ProcessTable' => '0',
	 'IO::All'            => '0',
	 'IO::Prompt'         => '0',
	 'File::Which'        => '0',
	 'Bio::Root::Version' => '1.006',
	 'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     exe_requires => {
	 'BLAST' => 'xdformat, formatdb',
	 'RepeatMasker' => 'RepeatMasker',
	 'Exonerate' => 'exonerate',
	 'SNAP'  => 'snap',
	 'Augustus' => 'augustus',
     },
     build_requires => {
	 'File::Which'        => '0',
	 'ExtUtils::CBuilder' => '0',
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
     PL_files => { 'bin/mpi_scripts.PL' => ['bin/mpi_maker',
					    'bin/mpi_evaluator',
					    'bin/mpi_iprscan',],
		   'lib/Parallel/MPIcar.PL' => ['lib/Parallel/MPIcar.pm',
						'lib/Parallel/MPIcar.xs',
						'lib/Parallel/utils.c',
						'lib/Parallel/typemap',]
     },
    );

#see if threaded
my $threaded = 1;
eval "require threads";
if($@){
    print "WARNING:  Your version of perl does not seem to have been compiled\n".
	  "for threads; certain MAKER accessory scripts and features (such as\n".
	  "MPI support) will be not be available to you\n";
    $threaded = 0;
}


#decide whether to configure maker MPI support
if($threaded && (my @mpi = $build->mpi_support)){
    my $go = $build->y_n("\nYour system appears to be configured for MPI. Would\n".
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

$build->create_build_script;

#term($build, 0);


sub term {
    my $build = shift;
    my $stat = shift;

    my @perl = map {keys %{$_->{requires}}} $build->prereq_failures();
    my @exes = $build->missing_exes();
    my $mpi = $build->check_installed_status('Parallel::MPIcar', '0')->{ok};

    $build->create_build_script;

    print "\n\nThe file 'Build' has been created for you to finish or rerun any\n".
	"steps that did not complete or failed during initial setup.\n\n";

    print "\n\n";
    print "================================================================\n";
    print "STATUS\n";
    print "================================================================\n";
    print "PERL Dependencies:\t";
    print ((@perl) ? 'MISSING' : 'INSTALLED');
    print"\n";
    print "\t!\t". join("\n\t!\t", @perl) ."\n\n" if(@perl);
    print "External Dependencies:\t";
    print ((@exes) ? 'MISSING' : 'INSTALLED');
    print "\n";
    print "\t!\t". join("\n\t!\t", @exes) ."\n\n" if(@exes);
    print "MPI MAKER (optional):\t";
    print (($mpi) ? 'INSTALLED' : 'NOT CONFIGURED');
    print "\n";
    print "MAKER:\t\t\t";
    print ((@perl || @exes) ? 'MISSING PREREQUISITES' : 'READY TO USE');
    print "\n";

    print "\n\nImportant Commands:\n".
	"\t./Build installdeps\t\#installs missing perl dependencies\n".
	"\t./Build installprereqs\t\#installs missing external program dependencies\n".
	"\t./Build mpi\t\t\#installs MPI MAKER and accessory scripts\n\n".
	"Other Commands:\n".
	"\t./Build repeatmasker\t\#installs just RepeatMasker\n".
	"\t./Build blast\t\t\#installs just BLAST\n".
	"\t./Build exonerate\t\#installs just Exonerate\n".
	"\t./Build snap\t\t\#installs just SNAP\n".
	"\t./Build augustus\t\#installs just Augustus\n";

    

    exit($stat);
}
