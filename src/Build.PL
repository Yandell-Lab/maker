#! /bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/inc/"; #libary just needed for configuration
use MAKER::Build;
eval 'require File::Which; 1;';

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
     requires => {
	 'perl'               => '5.8.0',
	 'DBI'                => '0',
	 'DBD::SQLite'        => '0',
	 'Proc::ProcessTable' => '0',
	 'IO::All'            => '0',
	 'IO::Prompt'         => '0',
	 'File::Which'        => '0',
	 'Bio::Root::Version' => '1.006',
	 'Module::Build::Functions'        => '0',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     build_requires => {
	 'File::Which'        => '0',
	 'ExtUtils::CBuilder' => '0',
     },
     recommends => {
	 'threads'  => '1.67',
     },
     install_base => "$FindBin::Bin/../perl",
     install_base_relpaths => {
	 arch   => [qw(lib)],
	 lib    => [qw(lib)],
	 bindoc => [qw(man)],
	 libdoc => [qw(man)],
     },
    );

#install missing perl dependencies
if($build->prereq_failures()){
    if($build->y_n("WARNING: You cannot continue without these dependencies\n".
		   "Do you want to let MAKER install these for you now?", 'Y')
       ){
	$build->dispatch('installdeps');

	print "\nRechecking dependencies to see if installation was successful\n";
	$build->check_prereq;
	if($build->prereq_failures()){
	    print "MAKER could not install these dependencies.  You will\n".
		"have to install them yourself befeore you can continue.\n";
	    exit(1);
	}
    }
    else{
	print "\nFAILURE: You will have to install these dependencies before you can continue.\n\n";
	exit(1);
    }
}

#install missing external program dependencies
$build->dispatch('checkprereqs');

#install maker MPI
if(my $dir = File::Which::which('mpicc')){
    $dir =~ s/[^\/]+\/mpicc$/include/;
    my ($MPIDIR) = grep {-f "$_/mpi.h"} (</usr/mpi*/include>,
					 </usr/local/mpi*/include>,
					 </usr/include/mpi*>,
					 </usr/local/include/mpi*>,
					 </usr/lib/mpi*/include>,
					 </usr/local/lib/mpi*/include>,
					 $dir);

    my $go = $build->y_n("Your system appears to be configured for MPI. Would\n".
			 "you like to install the MPI version of MAKER?", 'Y') if($MPIDIR);
    if($go){
	$build->config(cc => 'mpicc');
	$build->config(ld => 'mpicc');
	$build->extra_compiler_flags("-I$MPIDIR -DFLOAT_HACK");

	#Hack to force Build to use new config flags
	#This is needed because of a bug in Module::Build::Base::cbuilder
	$build->{stash}->{_cbuilder} = undef;
	
	my $ok = $build->check_installed_status('threads', '1.67')->{ok};
	
	if(! $ok){
	    $build->install_prereqs('threads');
	    $ok = $build->check_installed_status('threads', '1.67')->{ok};
	}
	
	if(! $ok){
	    print "\nCannot install MPI MAKER without threads\n\n";
	}
	else{
	    $build->dispatch('mpi');
	}
    }
}

$build->create_build_script;
#$build->dispatch('realclean');
