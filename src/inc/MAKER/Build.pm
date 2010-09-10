#------------------------------------------------------------------------
#----                          MAKER::Build                          ----
#------------------------------------------------------------------------
package MAKER::Build;
use strict;
use warnings;

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
    $self->depends_on('checkprereqs');
}

sub ACTION_mpi{
    my $self = shift;

    $self->SUPER::ACTION_install;
}


sub ACTION_checkprereqs {
    my $self = shift;

    print "Locating MAKER required prerequisite programs:\n";

    my @exes = qw(RepeatMasker
		  blast
		  exonerate
		  snap
		  augustus);

    my $missing;
    foreach my $exe (@exes){
	if(! $self->_found_exe($exe)){
	    print "\t!\t$exe\n";
	    $missing = 1;
	}
    }

    print "\n";
    my $go = $self->y_n("Do you want MAKER to the missing programs for you?", 'Y') if($missing);
    $self->dispatch('installprereqs') if($go);
}

sub ACTION_installprereqs{
    my $self = shift;

    $self->depends_on("repeatmasker") if(! $self->_found_exe('RepeatMasker'));
    $self->depends_on("blast")  if(! $self->_found_exe('blast'));
    $self->depends_on("exonerate") if(! $self->_found_exe('exonerate'));
    $self->depends_on("snap") if(! $self->_found_exe('snap'));
    $self->depends_on("augustus") if(! $self->_found_exe('augustus'));
}

sub _get_missing_exes{
    my $self = shift;
}

sub _exe_action{
    my $self = shift;
    my $exe = shift;

    if($self->_found_exe($exe)){
	my $go = $self->y_n("WARNING: $exe was already found on this system.\n",
			    "Do you still want MAKER to install $exe for you?", 'y');
	$self->_install_exe($exe) if($go);
    }
    else{
	$self->_install_exe($exe);
    }
}

sub _found_exe {
    my $self = shift;
    my $exe = shift;
    my $base = $self->base_dir();
    my $path = "$base/external/$exe";

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
    my $base = $self->base_dir();
    my $path = "$base/external/$exe";

    #get OS and architecture
    my ($OS, $ARC) = (POSIX::uname())[0,4];
    
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

    chdir($base);
    if($exe eq 'RepeatMasker'){
	#RepeatMasker
	my $file = 'RepeatMasker.tar.gz'; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to RepeatMasker for OS
        getstore($url, $file) or $self->fail($exe, $path);
        Archive::Tar->extract_archive($file) or $self->fail($exe, $path);
        unlink($file);
        chdir($path);
	$self->fail($exe, $path) if(! -f $exe);

	#RMBlast
        chdir($path);
	$file = 'rmblast.tar.gz'; #file to save to
        $url = $data->{rmblast}{"$OS\_$ARC"}; #url to rmblast for OS
        getstore($url, $file) or $self->fail('RMBlast', $path);
        Archive::Tar->extract_archive($file) or $self->fail('RMBlast', $path);
        unlink($file);
        my ($dir) = grep {-d $_} <rmblast-*-ncbi-blast*>;
        chdir($dir);
	$self->fail('RMBlast', $path) if(! -f 'rmblastn');

	#TRF
        chdir($path);
	$file = 'trf'; #file to save to
        $url = $data->{trf}{"$OS\_$ARC"}; #url to rmblast for OS
        getstore($url, $file) or $self->fail('TRF', $path);
        chmod(0755, $file) or $self->fail('TRF', $path);
	$self->fail('TRF', $path) if(! -f 'trf');
    }
    elsif($exe eq 'blast'){
        my $path = $self->base_dir()."/../external/blast"; #path for installation

	#BLAST+
	my $file = 'blast.tar.gz'; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
        getstore($url, $file) or $self->fail($exe, $path);
        Archive::Tar->extract_archive($file) or $self->fail($exe, $path);
        unlink($file);
        chdir($path);
	$self->fail($exe, $path) if(! -f 'bin/blastn');
    }

    chdir($self->base_dir());
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl
sub install_prereq {
    my ($self, $desired) = @_;
    
    my $do_install = $self->y_n("The module $desired is must be installed to continue for\n".
				"the selected options. Shall I install via CPAN?", 'Y');
    
    if ($do_install) {
	# Here we use CPAN to actually install the desired module
	require Cwd;
	require CPAN;
	
	# Save this because CPAN will chdir all over the place.
	my $cwd = Cwd::cwd();
	
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

1;
