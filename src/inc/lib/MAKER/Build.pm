#------------------------------------------------------------------------
#----                          MAKER::Build                          ----
#------------------------------------------------------------------------
package MAKER::Build;
use strict;
use warnings;
use POSIX;
use Archive::Tar;
use File::Copy;
use File::Path;
use File::Which; #bundled with MAKER

BEGIN{
    #prepare correct version of Module Build for building everything
    my $Bundled_MB = 0.3607;  #version included in my distribution                                                                 
    # Find out what version of Module::Build is installed or fail quietly.
    # This should be cross-platform.
    my $Installed_MB =`$^X -e "eval q{require Module::Build; print Module::Build->VERSION} or exit 1"`;
    chomp $Installed_MB;
    $Installed_MB = 0 if $?;

    # Use the bundled copy of Module::Build if it's newer than the installed.
    unshift @INC, "inc/bundle/" if $Bundled_MB > $Installed_MB;

    require Module::Build;
}

use base qw(Module::Build);
__PACKAGE__->add_property( 'exe_requires' );

eval 'require LWP::Simple';

#------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------
#------------------------------------------------------------------------
sub new {
    my $class = shift @_;

    my $self = $class->SUPER::new(@_);
    $self->install_base_relpaths('exe' => 'exe');
    bless($self, $class);

    #performs a check for eternal algorithm dependencies
    $self->check_exes;

    return $self;
}

#returns MPI compiler and includes directory (undef when one not found)
sub mpi_support {
    my $self = shift;

    my $cc = File::Which::which('mpicc') || $self->config('cc');
    return unless($cc =~ /mpicc$/);
    my $ccdir = $cc;
    $cc =~ s/[^\/]+\/[^\/]+$/include/;

    my ($MPIDIR) = grep {-f "$_/mpi.h"} (</usr/include>,
                                         </usr/include/mpi*>,
					 </usr/mpi*/include>,
                                         </usr/local/include>,
                                         </usr/local/include/mpi*>,
                                         </usr/local/mpi*/include>,
                                         </usr/lib/>,
                                         </usr/lib/include/mpi*>,
                                         </usr/lib/mpi*/include>,
					 </usr/local/lib>,
                                         </usr/local/lib/include/mpi*>,
                                         </usr/local/lib/mpi*/include>,
					 $ccdir);

    return if(! $MPIDIR);

    return ($cc, $MPIDIR);
}

#add a Module to the requires list
sub add_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{requires}{$_} = 0} @_;
    }
    
}

#replaces Module::Build's config method
sub config {
    my $self = shift;
    
    #hack to get around bug in Module::Build 
    #otherwise it will ignore changes to config 'cc' adn 'ld'
    $self->{stash}->{_cbuilder} = undef;

    return $self->SUPER::config(@_);
}

#for building a release version of MAKER, updates MAKER version part of the commit
sub ACTION_commit {
    my $self = shift;
    
    #$self->depends_on("test");
    #$self->do_system(qw(svn commit));
}

#will install MAKER, but first installs the dependencies
#replaces Module::Build's ACTION_install
sub ACTION_install{
    my $self = shift;

    $self->depends_on('installdeps');
    $self->depends_on('installprereqs');
}

#replaces Module::Build's ACTION_installdeps
#so MAKER prereqs install locally inside of maker/perl/lib
sub ACTION_installdeps{
    my $self = shift;

    my $global = $self->y_n("By default MAKER will install these dependencies locally (i.e. only\n".
			    "for MAKER), would you rather install these globally (usually requires\n".
			    "that you be logged in as 'root' or run with sudo)?", 'N');
    
    my $prereq = $self->prereq_failures();
    
    if($prereq && $prereq->{requires}){
	foreach my $m (keys %{$prereq->{requires}}){    
	    $self->cpan_install($m, $global);
	}
    }

    print "\nRechecking dependencies to see if installation was successful\n";
    $self->check_prereq;

    if($self->prereq_failures()){
	my ($usr_id) = (getpwnam('root'))[2];
	print "WARNING: ";
	print "Try installing the missing packages as 'root' or using sudo.\n" if($< != $usr_id);
	print "You may need to configure and install these packages manually.\n";
	return 0;
    }
}

#these install individual external algorithms
sub ACTION_repeatmasker{ shift->_exe_action('RepeatMasker'); }
sub ACTION_blast{ shift->_exe_action('blast'); }
sub ACTION_exonerate{ shift->_exe_action('exonerate'); }
sub ACTION_snap{ shift->_exe_action('snap'); }
sub ACTION_augustus{ shift->_exe_action('augustus'); }

#runs all the algorithm installs that are missing
sub ACTION_installexes{
    my $self = shift;

    $self->dispatch("repeatmasker") if(! $self->_found_exe('RepeatMasker'));
    $self->dispatch("blast")  if(! $self->_found_exe('blast'));
    $self->dispatch("exonerate") if(! $self->_found_exe('exonerate'));
    $self->dispatch("snap") if(! $self->_found_exe('snap'));
    $self->dispatch("augustus") if(! $self->_found_exe('augustus'));
    
    $self->check_exes;
}

#prints out a simle configuration status message
sub ACTION_status {
    my $self = shift;
    
    $self->maker_status;
}

#checks external algorithms to see if they're present. Anologous to check_prereqs
sub check_exes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();
    return if(! $exe_failures);

    print "Checking external program dependencies...\n";
    while(my $cat = each %$exe_failures){
	my $s = ($cat eq 'exe_requires') ? '!' : '*';
	$cat =~ /exe_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$exe_failures->{$cat}}){
	    print "    $s  ".$exe_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the programs\n".
	"indicated above before proceeding with this installation\n".
	"\nRun 'Build installexes' to install missing prerequisites.\n";
}

#returns missing exes, anologous to prereq_failures
sub exe_failures {
    my $self = shift;
    my %exes = %{$self->exe_requires};
    
    my %exe_failures;
    while (my $name = each %exes){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $exes{$name});
	my $base = $self->install_destination('exe');
	my $dest = "$base/$name";

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {-f $_} ("$dest/$cmd", "$dest/bin/$cmd", File::Which::which($cmd)));
	}

        if(! $loc){
	    $exe_failures{'exe_requires'}{$name}{have} = '<none>';
	    $exe_failures{'exe_requires'}{$name}{message} = "$name is not installed";
	    $exe_failures{'exe_requires'}{$name}{need} = $exes{$name};
	}
    }

    return (keys %exe_failures) ? \%exe_failures : undef;
}

#hidden connection entry between ACTION_??? algorithm install
#checks to see if already installed and then run the install method
sub _exe_action{
    my $self = shift;
    my $exe = shift;

    if($self->_found_exe($exe)){
	my $go = $self->y_n("WARNING: $exe was already found on this system.\n".
			    "Do you still want MAKER to install $exe for you?", 'N');
	$self->_install_exe($exe) if($go);
    }
    else{
	$self->_install_exe($exe);
    }
}

#does actual installation of all external algorithms
sub _install_exe {
    my $self = shift;
    my $exe  = shift;
    my $base = $self->base_dir()."/../external/";
    my $path = "$base/$exe";

    #get OS and architecture
    my %os_ok = (Linux_x86_64 => 1,
		 Linux_i386 => 1,
		 Darwin_i386 => 1,
		 src => 1); 
    my ($OS, $ARC) = (POSIX::uname())[0,4];

    ($OS, $ARC) = ('src', '') if(! $os_ok{"$OS\_$ARC"}); #use source code for unknown architectures

    #get url for exectuable to be installed
    my $data;
    open(LOC, '<', $self->base_dir()."/locations")
	or die "ERROR: Could not open locations to download dependencies\n";
    my $line = <LOC>;
    if($line =~ /^\#\#MAKER/){
	$data = join('', <LOC>);
	eval $data;
    }
    close(LOC);

    #maker prerequisite installation directory
    if(! -d $base){
	mkdir($base) ||
	    die "ERROR could not create directory $base for installing external dependencies\n";
    }

    #install
    chdir($base);
    my @unlink;
    if($exe eq 'RepeatMasker'){
	#RepeatMasker
	my $file = "$base/$exe.tar.gz"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to RepeatMasker for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push(@unlink, $file);
	return $self->fail($exe, $path) if(! -f "$path/$exe");

	#RMBlast
	my $req = 'rmblast';
        chdir($path);
	File::Path::rmtree("$path/$req");
	$file = "$path/$req.tar.gz"; #file to save to
        $url = $data->{$req}{"$OS\_$ARC"}; #url to rmblast for OS
	print "Downloading $req...\n";
        $self->getstore($url, $file) or return $self->fail($req, $path);
	print "Unpacking $req tarball...\n";
        $self->extract_archive($file) or return $self->fail($req, $path);
        push(@unlink, $file);
	my ($dir) = grep {-d $_} <rmblast-*-ncbi-blast*>;
	if(-d "$dir/c++"){ #this is the source code and must be compiled
	    chdir("$dir/c++");
	    print "Configuring $req...\n";
	    $self->do_system("./configure --with-mt --prefix=$path/$dir --without-debug") or return $self->fail($req, $path);
	    $self->do_system("make") or return $self->fail($req, $path);
	    $self->do_system("make install") or return $self->fail($req, $path);
	}
        chdir($path);
	File::Copy::move($dir, $req);
	return $self->fail($req, $path) if(! -f "$path/$req/bin/rmblastn");

	#TRF
	$req = 'trf';
        chdir($path);
	unlink("$path/$req");
	$file = "$path/$req"; #file to save to
        $url = $data->{$req}{"$OS\_$ARC"}; #url to rmblast for OS
	print "Downloading $req...\n";
        $self->getstore($url, $file) or return $self->fail($req, $path);
        chmod(0755, $file) or return $self->fail($req, $path);
	return $self->fail($req, $path) if(! -f "$path/$req");

	#Configure RepeatMasker
	chdir($path);
	print "Configuring $exe...\n";
	my $tmp = "$path/.config.stdin";
	open(my $fh, "> $tmp");
	print $fh "\n";
	print $fh "$^X\n";
	print $fh "$path\n";
	print $fh "$path\n";
	print $fh "2\n";
	print $fh "$path/rmblast/bin\n";
	print $fh "Y\n";
	print $fh "5\n";
	close($fh);
	$self->do_system("$^X ./configure < $tmp > /dev/null") or return $self->fail($exe, $path);
        push(@unlink, $file);
        push(@unlink, $tmp);
    }
    elsif($exe eq 'blast'){
	#BLAST+
	File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <ncbi-blast-*>;
	if(-d "$dir/c++"){ #this is the source code and must be compiled
	    chdir("$dir/c++");
	    print "Configuring $exe...\n";
	    $self->do_system("./configure --prefix=$path") or return $self->fail($exe, $path);
	    $self->do_system("make") or return $self->fail($exe, $path);
	    $self->do_system("make install") or return $self->fail($exe, $path);
	}
        chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);;
	return $self->fail($exe, $path) if(! -f "$path/bin/blastn");
    }
    elsif($exe eq 'exonerate'){
	#exonerate
        File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <exonerate-*>;
	if(-d "$dir/src"){ #this is the source code and must be compiled
	    chdir($dir);
	    print "Configuring $exe...\n";
	    $self->do_system("./configure --prefix=$dir") or return $self->fail($exe, $path);
	    $self->do_system("make") or return $self->fail($exe, $path);
	    $self->do_system("make install") or return $self->fail($exe, $path);
	}
        chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);

	return $self->fail($exe, $path) if(! -f "$path/bin/$exe");
    }
    elsif($exe eq 'snap'){
	#snap
        File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	chdir($path);
	print "Configuring $exe...\n";
	$self->do_system("make") or return $self->fail($exe, $path);
	return $self->fail($exe, $path) if(! -f "$path/$exe");
    }
    elsif($exe eq 'augustus'){
	#augustus
        &File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	my ($dir) = grep {-d $_} <augustus-*>;
	chdir("$dir/src");
	print "Configuring $exe...\n";
	$self->do_system("make") or return $self->fail($exe, $path);
	chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
	return $self->fail($exe, $path) if(! -f "$path/bin/$exe");
    }
    else{
	die "ERROR: No install method defined for $exe in MAKER::Build::_install_exe.\n";
    }

    #report status
    print "Finished installing $exe.\n";

    #remove all tarballs/etc
    map{unlink($_)} @unlink;

    #change back to base
    chdir($self->base_dir());
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl, has flag for local install
sub cpan_install {
    my ($self, $desired, $global) = @_;
    
    # Here we use CPAN to actually install the desired module
    require Cwd;
    require CPAN;
    
    # Save this because CPAN will chdir all over the place.
    my $cwd = Cwd::cwd();
    
    if(! $global){
	my $base = $self->base_dir;
	CPAN::HandleConfig->load;
	$CPAN::Config->{make_install_arg} = "DESTDIR=$base/../perl/ INSTALLDIRS=perl INSTALLMAN1DIR=man".
	    " INSTALLMAN3DIR=man INSTALLARCHLIB=lib INSTALLPRIVLIB=lib INSTALLBIN=bin";
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }
    CPAN::Shell->install($desired);
    
    my $ok;
    my $expanded = CPAN::Shell->expand("Module", $desired);
    if ($expanded && $expanded->uptodate) {
	print "$desired installed successfully\n";
	$ok = 1;
    }
    else {
	print "$desired failed to install\n";
	$ok = 0;
    }
    
    chdir $cwd or die "Cannot chdir() back to $cwd: $!";
    return $ok;
}

#untars a package. Tries to use tar first then moves to the perl package untar Module.
sub extract_archive {
    my $self = shift;
    my $file = shift;

    return 0 if(! $file);
    
    if(File::Which::which('tar')){
	return $self->do_system("tar -xf $file"); #fast
    }
    else{
	return (Archive::Tar->extract_archive($file)) ? 1 : 0; #slow
    }
}

#downloads files from the internet.  Tries to use wget, then curl,
#and finally LWP::Simple
sub getstore {
    my $self = shift;
    my $url = shift;
    my $file = shift;
    
    if(File::Which::which('wget')){#Linux
	return $self->do_system("wget $url -c -O $file"); #gives status and can continue partial
    }
    elsif(File::Which::which('curl')){#Mac
	return $self->do_system("curl $url -C - -o $file"); #gives status and can continue partial
    }
    else{
	return LWP::Simple::getstore($url, $file); #just gets the file with no features
    }
}

#prints a nice status message for package configuration and install
sub maker_status {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my $mpi = $self->feature('mpi_support');

    $self->create_build_script;

    print "\n\nThe file 'Build' has been created for you to finish installing MAKER.\n";

    print "\n\n";
    print "================================================================\n";
    print "STATUS\n";
    print "================================================================\n";
    print "PERL Dependencies:\t";
    print ((@perl) ? 'MISSING' : 'INSTALLED');
    print"\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @perl) ."\n\n" if(@perl);
    print "External Dependencies:\t";
    print ((@exes) ? 'MISSING' : 'INSTALLED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @exes) ."\n\n" if(@exes);
    print "MPI SUPPORT:\t\t";
    print (($mpi) ? 'CONFIGURED' : 'NOT CONFIGURED');
    print "\n";
    print "MAKER:\t\t\t";
    print ((@perl || @exes) ? 'MISSING PREREQUISITES' : 'READY TO INSTALL');
    print "\n";

    print "\n\nImportant Commands:\n".
        "\t./Build installdeps\t\#installs missing perl dependencies\n".
        "\t./Build installexes\t\#installs missing external program dependencies\n".
        "\t./Build install\t\t\#installs MAKER\n".
        "\t./Build status\t\t\#Shows this status menu\n\n".
        "Other Commands:\n".
        "\t./Build repeatmasker\t\#installs just RepeatMasker\n".
        "\t./Build blast\t\t\#installs just BLAST\n".
        "\t./Build exonerate\t\#installs just Exonerate\n".
        "\t./Build snap\t\t\#installs just SNAP\n".
        "\t./Build augustus\t\#installs just Augustus\n";
}

1;
