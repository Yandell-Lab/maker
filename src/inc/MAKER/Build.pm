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

BEGIN{
    #prepare correct version of Module Build for building everything
    my $Bundled_MB = 0.3607;  #version included in my distribution                                                                 
    # Find out what version of Module::Build is installed or fail quietly.
    # This should be cross-platform.
    my $Installed_MB =`$^X -e "eval q{require Module::Build; print Module::Build->VERSION} or exit 1"`;
    chomp $Installed_MB;
    $Installed_MB = 0 if $?;

    # Use the bundled copy of Module::Build if it's newer than the installed.
    unshift @INC, "inc/Module-Build/lib/perl5" if $Bundled_MB > $Installed_MB;

    require Module::Build;
}

use base qw(Module::Build);

eval 'require File::Which; 1;';
eval 'require LWP::Simple; 1;';

#------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------
#------------------------------------------------------------------------

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

sub missing_exes{
    my $self = shift;

    my @exes = qw(RepeatMasker
                  blast
                  exonerate
                  snap
                  augustus);

    my @missing;
    foreach my $exe (@exes){
        if(! $self->_found_exe($exe)){
            push(@missing, $exe);
        }
    }

    return @missing;
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

sub _found_exe {
    my $self = shift;
    my $exe = shift;
    my $base = $self->base_dir()."/../external";
    my $path = "$base/$exe";

    my $found;
    if($exe eq 'blast'){
	$found = File::Which::which('xdformat') || File::Which::which('formatdb') || -f "$path/bin/formatdb";
    }
    else{
	$found = File::Which::which($exe) || -f "$path/bin/$exe" || -f "$path/$exe";
    }

    return ($found) ? 1 : 0;
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

    #report status
    print "Finished installing $exe.\n";

    #remove all tarballs/etc
    map{unlink($_)} @unlink;

    #change back to base
    chdir($self->base_dir());
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl
sub install_prereq {
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
