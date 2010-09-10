#! /bin/env perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib/"; #libary just needed for configuration
use lib "$FindBin::Bin/../perl/lib/"; #libary just needed for configuration
use lib "$FindBin::Bin/inc/"; #libary just needed for configuration
use MAKER::Build;
eval 'require File::Which; 1;';
eval 'require local::lib; 1;';

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
     config => {cc => 'mpicc',
		ld => 'mpicc'},
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
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     build_requires => {
	 'File::Which'        => '0',
	 'ExtUtils::CBuilder' => '0',
	 'LWP::Simple'        => '0',
     },
     recommends => {
	 'threads'  => '1.67',
     },
     install_base => "$FindBin::Bin/../",
     install_base_relpaths => {
	 arch   => [qw(perl/lib)],
	 lib    => [qw(perl/lib)],
	 bindoc => [qw(perl/man)],
	 libdoc => [qw(perl/man)],
     },
    );

#install missing perl dependencies
if($build->prereq_failures()){
    if($build->y_n("WARNING: You cannot continue without these dependencies\n".
		   "Do you want to let MAKER install these for you now?", 'Y')
       ){
	#try and use local lib
	my $ok = $build->check_installed_status('local::lib', '0')->{ok};
        if(! $ok){
            $build->install_prereq('local::lib', 1);
            $ok = $build->check_installed_status('local::lib', '0')->{ok};
        }
	require local::lib if($ok);

	#install dependencies
	$build->dispatch('installdeps');

	print "\nRechecking dependencies to see if installation was successful\n";
	$build->check_prereq;
	if($build->prereq_failures()){
	    print "\nMAKER could not install these dependencies.  You will\n".
		"have to install them yourself before you can continue.\n";
	    term($build, 1);
	}
    }
    else{
	print "\nFAILURE: You will have to install these dependencies before you can continue.\n\n";
	    term($build, 1);
    }
}

#install missing external program dependencies
print "Locating MAKER required external prerequisite programs:\n";
if(my @missing = $build->missing_exes){
    print "\t!\t".join("\n\t!\t", @missing)."\n";
    my $go = $build->y_n("\nDo you want MAKER to the missing programs for you?", 'Y');
    if($go){
	$build->dispatch('installprereqs');
    }
    if($build->missing_exes){
	print "\nYou will have to install these programs manually\n";
    }
}

#install maker MPI
my $mpi = $build->check_installed_status('Parallel::MPIcar', '0')->{ok};
if(my $cc = File::Which::which('mpicc') && ! $mpi){
    my $go = $build->y_n("\nYour system appears to be configured for MPI. Would\n".
			 "you like to install the MPI version of MAKER?", 'Y') if($cc);
    if($go){	
	my $ok = $build->check_installed_status('threads', '1.67')->{ok};
	
	if(! $ok){
	    $build->install_prereq('threads');
	    $ok = $build->check_installed_status('threads', '1.67')->{ok};
	}
	
	if(! $ok){
	    print "\nCannot install MPI MAKER without threads\n\n";
	}
	else{
	    $build->dispatch('mpi');
	    $mpi = $build->check_installed_status('Parallel::MPIcar', '0')->{ok};
	    if(! $mpi){
		print "\nERROR: Installation of MPI MAKER failed.\n";
		term($build, 1);
	    }
	}
    }
}

term($build, 0);
#$build->dispatch('realclean');

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
