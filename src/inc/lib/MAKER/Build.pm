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
use File::Which;

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

    my $ok = $self->check_installed_status('File::Which', '0')->{ok};
    if($ok){
	$self->check_exes;
    }
    else{
	die "ERROR: File::Which is required by MAKER::Build. Please install.\n";
    }

    return $self;
}

#returns MPI compiler and includes directroy if availblae
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

sub add_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	while(@_){$self->{properties}{requires}{shift @_} = shift @_}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{requires}{$_} = 0} @_;
    }
    
}

sub config {
    my $self = shift;
    
    #hack to get around bug in Module::Build 
    #otherwise it will ignore changes to config 'cc' adn 'ld'
    $self->{stash}->{_cbuilder} = undef;

    return $self->SUPER::config(@_);
}

sub ACTION_commit {
    my $self = shift;
    
    #$self->depends_on("test");
    #$self->do_system(qw(svn commit));
}

sub ACTION_repeatmasker{ shift->_exe_action('RepeatMasker'); }
sub ACTION_blast{ shift->_exe_action('blast'); }
sub ACTION_exonerate{ shift->_exe_action('exonerate'); }
sub ACTION_snap{ shift->_exe_action('snap'); }
sub ACTION_augustus{ shift->_exe_action('augustus'); }

sub ACTION_install{
    my $self = shift;

    $self->depends_on('installdeps');
    $self->depends_on('installprereqs');
}

sub ACTION_installdeps{
    my $self = shift;

    my $global = $self->y_n("By default MAKER will install these dependencies locally (i.e. only\n".
			     "for MAKER), would you rather install these globally (usually requires\n".
			     "that you be logged in as 'root' or run with sudo)?", 'N');
    
    my @prereq = $self->prereq_failures();

    foreach my $m (@prereq){
	
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

sub ACTION_mpi{
    my $self = shift;

    $self->depends_on('installdeps');

    my $ccdir = $self->config('cc') || '';
    $ccdir =~ s/[^\/]+\/mpicc$/include/;
    my $dir = File::Which::which('mpicc') || '';
    $dir =~ s/[^\/]+\/mpicc$/include/;

    my ($MPIDIR) = grep {-f "$_/mpi.h"} (</usr/mpi*/include>,
                                         </usr/local/mpi*/include>,
                                         </usr/include/mpi*>,
                                         </usr/local/include/mpi*>,
                                         </usr/lib/mpi*/include>,
                                         </usr/local/lib/mpi*/include>,
					 $ccdir,
                                         $dir);

    die "ERROR: Can't find mpi.h\n" if(! $MPIDIR);
    
    $self->extra_compiler_flags("-I$MPIDIR -DFLOAT_HACK");
    
    $self->SUPER::ACTION_install;
}

sub ACTION_installprereqs{
    my $self = shift;

    $self->depends_on("repeatmasker") if(! $self->_found_exe('RepeatMasker'));
    $self->depends_on("blast")  if(! $self->_found_exe('blast'));
    $self->depends_on("exonerate") if(! $self->_found_exe('exonerate'));
    $self->depends_on("snap") if(! $self->_found_exe('snap'));
    $self->depends_on("augustus") if(! $self->_found_exe('augustus'));
}

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
    if($line =~ /^\#\! \/bin\/env perl/){
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
# borrowed and modified from BioPerl
sub cpan_install {
    my ($self, $desired, $local) = @_;
    
    my $do_install = $self->y_n("The module $desired is must be installed to continue for\n".
				"the selected options. Shall I install via CPAN?", 'Y');
    
    if ($do_install) {
	# Here we use CPAN to actually install the desired module
	require Cwd;
	require CPAN;
	
	# Save this because CPAN will chdir all over the place.
	my $cwd = Cwd::cwd();
	
	if($local){
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
    else {
	print "You chose not to install $desired.\n";
	return 0;
    }
}

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

sub getstore {
    my $self = shift;
    my $url = shift;
    my $file = shift;
    
    if(File::Which::which('wget')){
	return $self->do_system("wget $url -c -O $file"); #gives status and can continue partial
    }
    else{
	return LWP::Simple::getstore($url, $file); #just gets the file with no features
    }
}

1;
