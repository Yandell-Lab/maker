#! /bin/env perl
use Module::Build;
use File::Which;
use POSIX;
use FindBin;
use LWP::Simple;
use Archive::Tar;

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
    $build->dispatch('installdeps')
	if($build->y_n("\nDo you want to let MAKER install these dependecies for you now?", 'y'));
}

#check external program dependencies
print STDERR "Checking external program dependencies.\n";
my @exes =qw(RepeatMasker
	     BLAST
	     exonerate
	     snap
	     augustus
	     );

foreach my $exe (@exes){
    print STDERR "Checking if you have $exe.\n";
    
    if($exe eq 'BLAST'){
	if(which('xdformat') || which('formatdb') || -f "$FindBin::Bin/prereqs/blast/bin/formatdb"){
	    print "\tSuccess: $exe found!\n\n";
	}
	else{
	    print "\tFailure: $exe does not appear to be installed!\n";
	    install($exe)if($build->y_n("Do you want MAKER to install $exe for you now?", 'y'));
	}
    }
    else{
	if(which($exe) || -f "$FindBin::Bin/$exe/bin/$exe" || -f "$FindBin::Bin/$exe/$exe"){
	    print "\tSuccess: $exe found!\n\n";
	}
	else{
	    print "\tFailure: $exe does not appear to be installed!\n";
	    install($exe)if($build->y_n("Do you want MAKER to install $exe for you now?", 'y'));
	}
    }
}

#set up MPI support
if(which('mpdboot')){
    $build->dispatch('mpi')
	if($build->y_n("\nYou appear to have MPICH2 installed on your system, do you".
		       " want to install MPI MAKER and accessory scripts?", 'y'));
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

sub install {
    my $exe  = shift; #prerequisite
    
    #get OS and architecture
    my ($OS, $ARC) = (POSIX::uname())[0,4];
    
    #get urls of exectuables to be installed
    my $locations;
    open(LOC, "< $FindBin::Bin/lib/location.conf") || die "ERROR: Could not open locations to download dependencies\n";
    my $line = <LOC>;
    if($line =~ /^\#\! \/bin\/env perl/){
	my $data = join('', <LOC>);
	eval $data;
    }
    close(LOC);

    chdir($FindBin::Bin);
    if($exe eq 'RepeatMasker'){
        my $path = "$FindBin::Bin/RepeatMasker"; #path for installation

	#RepeatMasker
	my $file = 'RepeatMasker.tar.gz'; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to RepeatMasker for OS
        getstore($url, $file) or fail($exe, $path);
        Archive::Tar->extract_archive($file) or fail($exe, $path);
        unlink($file);
        chdir($path);
	fail($exe, $path) if(! -f $exe);

	#RMBlast
        chdir($path);
	$file = 'rmblast.tar.gz'; #file to save to
        $url = $data->{rmblast}{"$OS\_$ARC"}; #url to rmblast for OS
        ($url, $file) = @{$src}{"$OS\_$ARC", 'file'} ;
        getstore($url, $file) or fail('RMBlast', $path);
        Archive::Tar->extract_archive($file) or fail('RMBlast', $path);
        unlink($file);
        my ($dir) = grep {-d $_} <rmblast-*-ncbi-blast*>;
        chdir($dir);
	fail('RMBlast', $path) if(! -f 'rmblastn');

	#TRF
        chdir($path);
	$file = 'trf'; #file to save to
        $url = $data->{trf}{"$OS\_$ARC"}; #url to rmblast for OS
        getstore($url, $file) or fail('TRF', $path);
        chmod(0755, $file) or fail('TRF', $path);
	fail('TRF', $path) if(! -f 'trf');
    }
    elsif($exe eq 'blast'){
        my $path = "$FindBin::Bin/blast"; #path for installation

	#BLAST+
	my $file = 'blast.tar.gz'; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
        getstore($url, $file) or fail($exe, $path);
        Archive::Tar->extract_archive($file) or fail($exe, $path);
        unlink($file);
        chdir($path);
	fail($exe, $path) if(! -f 'bin/blastn');
    }

    chdir($FindBin::Bin);
}

sub fail {
    chdir($FindBin::Bin);
}
