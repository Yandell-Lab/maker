#! /bin/env perl
use Module::Build;
use FindBin;

#find MAKER version
my $VERSION = version();

my $build = Module::Build->new
    (
     module_name => 'MAKER',
     dist_author => 'Yandell Lab <maker-devel@yandell-lab.org>',
     dist_version  => $VERSION,
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
	 'perl'          => '5.8.0',
	 'DBI'  => '0',
	 'DBD::SQLite'  => '0',
	 'Proc::ProcessTable'  => '0',
	 'IO::All'  => '0',
	 'IO::Prompt'  => '0',
	 'Bio::Root::Version'  => '1.006',
	 'File::Which' => '0',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     recommends => {
	 'threads'  => '1.67',
     },
     auto_features => {
	 mpi_support => {
	     description => 'Optimize for parrallelization using MICH2',
	     requires => {
		 'threads' => '1.67',
	     },
	 },
     },
    );

#install missing perl dependencies
if($build->prereq_failures()){
    my $go = $build->y_n("Do you want to let MAKER install these dependecies for you now?", 'y');
    $build->dispatch('installdeps') if($go);
}


$build->create_build_script;


sub version {
    my $ver;

    if(open(IN, "< $FindBin::Bin/lib/version")){
	$ver = <IN>;
	close(IN);
    }
    else{
	$ver = '2'; # will be MAKER 2.X or greater
    }

    chomp($ver);
    ($ver) = $ver =~ /VERSION=(.*)/;
    $ver = '2' if(! $ver || ! $ver =~ /^\d+\.\d+$/ || $ver < 2); #must be at least 2

    return $ver;
}
