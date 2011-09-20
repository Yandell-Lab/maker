#------------------------------------------------------------------------
#----                          MAKER::Build                          ----
#------------------------------------------------------------------------
package MAKER::Build;
use strict;
use warnings;
use POSIX;
use Config;
use DynaLoader;
use File::Copy;
use File::Path;
use File::Which; #bundled with MAKER
use vars qw($BIN);
use Cwd ();

BEGIN{
    #prepare correct version of Module Build for building everything
    my $Bundled_MB = 0.3607;  #version included in my distribution

    # Find out what version of Module::Build is installed or fail quietly.
    # This should be cross-platform.
    my $Installed_MB =`$^X -e "eval q{require Module::Build; print Module::Build->VERSION} or exit 1"`;
    chomp $Installed_MB;
    $Installed_MB = 0 if $?;

    $BIN = Cwd::cwd();
    if (! $ENV{CGL_SO_SOURCE}) {
	$ENV{CGL_SO_SOURCE} = "$BIN/../lib/CGL/so.obo";
	$ENV{CGL_SO_SOURCE} = "$BIN/../../lib/CGL/so.obo"
	    if(! -f $ENV{CGL_SO_SOURCE});
    }
    if (! $ENV{CGL_GO_SOURCE}) {
	$ENV{CGL_GO_SOURCE} = "$BIN/../lib/CGL/gene_ontology.obo";
	$ENV{CGL_GO_SOURCE} = "$BIN/../../lib/CGL/gene_ontology.obo"
	    if(! -f $ENV{CGL_GO_SOURCE});
    }

    # Use the bundled copy of Module::Build if it's newer than the installed.
    if ($Bundled_MB > $Installed_MB){
	unshift @INC, "$BIN/inc/bundle" unless($INC[0] eq "$BIN/inc/bundle");
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = "$BIN/inc/bundle:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }

    require Module::Build;
}

use base qw(Module::Build);
__PACKAGE__->add_property( 'exe_requires' );
__PACKAGE__->add_property( 'lib_requires' );

eval 'require LWP::Simple';
eval 'require Archive::Tar';

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
    $self->check_libs;

    return $self;
}

#returns MPI compiler and includes directory
sub config_mpi {
    my $self = shift;

    my $ebase = $self->install_destination('exe');
    my @exes = grep {/(^|[\/])mpicc$/} (<$FindBin::Bin/../exe/*/*>, <$FindBin::Bin/../exe/*/bin/*>);
    my $mpicc = "$ebase/mpich2/bin/mpicc" if(-f "$ebase/mpich2/bin/mpicc");
    ($mpicc) = grep {/(^|[\/])mpicc$/} (<$FindBin::Bin/../exe/*/*>, <$FindBin::Bin/../exe/*/bin/*>) if(!$mpicc);
    $mpicc = $self->config('cc') if(! $mpicc && $self->config('cc') =~ /(^|[\/])mpicc$/);
    ($mpicc) = File::Which::where('mpicc') if(!$mpicc || ! -f $mpicc);

    $mpicc = $self->prompt("\nPlease specify the path to 'mpicc' on your system:", $mpicc);

    while($mpicc !~ /(^|[\/])mpicc$/ || ! -f $mpicc){
	$mpicc = $self->prompt("\nCannot find 'mpicc'.\n".
			    "Please specify the path (leave blank to cancel):", '');
	return if(! $mpicc);
    }

    my $ccdir = $mpicc;
    $ccdir =~ s/\/+[^\/]+\/[^\/]+$//;

    #directories to search for mpi.h
    my @includes = (<$ccdir/include>,
		    <$ebase/mpich2/include>,
		    <$ebase/*/include>,
		    </usr/include>,
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
		    </usr/local/lib/mpi*/include>);

    my ($MPIDIR) = grep {-f "$_/mpi.h"} @includes;

    $MPIDIR = $self->prompt("\nPlease specify the path to the directory containing 'mpi.h':", $MPIDIR);

    while(! -f "$MPIDIR/mpi.h"){
	$MPIDIR = $self->prompt("\nCannot find 'mpi.h'\n.".
				"Please specify the containing directory path (leave blank to cancel):", '');
	return if(! $MPIDIR);
    }

    $self->feature(mpi_support => 1);
    $self->config_data(MPIDIR => $MPIDIR);
    $self->config_data(MPICC => $mpicc);

    $self->add_exe_requires(mpicc => $mpicc);
    $self->add_lib_requires(MPI => "$MPIDIR/mpi.h");
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

#add a Module to the build_requires list
sub add_build_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{build_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{build_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{build_requires}{$_} = 0} @_;
    }
}

#add a Module to the exe_requires list
sub add_exe_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{exe_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{exe_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{exe_requires}{$_} = 0} @_;
    }
}

#add a Library to the lib_requires list
sub add_lib_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{lib_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{lib_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{lib_requires}{$_} = 0} @_;
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

#override default build
sub ACTION_build {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();
    
    if($self->feature('mpi_support')){
	$self->log_info("Configuring " . $self->dist_name . " with MPI support\n");
	die "\n* MISSING MAKER/MPI PREREQUISITES - CANNOT CONTINUE!!\n"
	    if(scalar(grep {!/^(CGI|Mail|Bio\:\:Graphics)/} @perl) ||
	       scalar(grep {/^MPI/} @exes, @libs));
    }
    else{
	$self->log_info("Building " . $self->dist_name . "\n");
	die "\n* MISSING MAKER PREREQUISITES - CANNOT CONTINUE!!\n"
	    if(scalar(grep {!/^(CGI|Mail|Bio\:\:Graphics)/} @perl));
    }
    $self->depends_on('code');
    $self->depends_on('docs');

    if($self->feature('mwas_support')){
	$self->log_info("Building MWAS\n");

        die "\n* MISSING MWAS PREREQUISITES - CANNOT CONTINUE!!\n"
            if(scalar(grep {!/^MPI/} @libs)||
               scalar(grep {/^(CGI|Mail|Bio\:\:Graphics)/} @perl));
    }    

    if($self->feature('mwas_support') && $self->invoked_action() eq 'build'){	
	require MWAS_util;

	while(! -d $self->config_data('APOLLO_ROOT')){
	    $self->set_gmod_locs('APOLLO_ROOT');
	}
	while(! -d $self->config_data('JBROWSE_ROOT')){
	    $self->set_gmod_locs('JBROWSE_ROOT');
	}
	while(! -f $self->config_data('GBROWSE_MASTER')){
	    $self->set_gmod_locs('GBROWSE_MASTER');
	}

	require GI;
	
	$main::server = 1;
	my $c_dir = MWAS_util::config_loc(); #configuration file directory
	my @files = ("$c_dir/maker_opts.ctl",
		     "$c_dir/maker_bopts.ctl",
		     "$c_dir/maker_exe.ctl",
		     "$c_dir/server.ctl",
		     "$c_dir/menus.ctl",
		     );

	my @exists = grep {-f $_} @files;
	my %O = GI::load_server_files(\@exists);

	$O{apache_user} = $self->prompt("\nWho is your apache user?", $O{apache_user});
	$O{web_address} = $self->prompt("\nWhat is the base URL to the server hosting MWAS?",
					 $O{web_address});
	$O{cgi_dir}     = $self->prompt("\nWeb accessible directory to install CGI files?",
					 $O{cgi_dir});
	$O{cgi_web}     = $self->prompt("\nURL to CGI directory (can be relative to base URL)?",
					 $O{cgi_web});
	$O{html_dir}    = $self->prompt("\nWeb accessible directory to install HTML files?",
					 $O{html_dir});
	$O{html_web}    = $self->prompt("\nURL to CGI directory (can be relative to base URL)?",
					 $O{html_web});
	$O{data_dir}    = $self->prompt("\nData directory for saving user files/results?",
					 $O{data_dir});
	$O{font_file}   = $self->prompt("\nFont file for webpage CAPTCHA?",
					 $O{font_file});
	$O{soba_url}    = $self->prompt("\nURL to Sequence Ontology SOBA CGI script?",
					 $O{soba_url});
	$O{DBI}         = $self->prompt("\nDBI interface type to use for database?",
					 $O{DBI});
	$O{dbname}      = $self->prompt("\nDatabase name?",
					 $O{dbname});
	$O{host}        = $self->prompt("\nDatabase host?",
					 $O{host}) unless($O{DBI} eq 'SQLite');
	$O{port}        = $self->prompt("\nDatabase port?",
					 $O{port}) unless($O{DBI} eq 'SQLite');
	$O{username}    = $self->prompt("\nDatabase username?",
					 $O{username});
	$O{password}    = $self->prompt("\nDatabase password?",
					 $O{password});
	$O{admin_email} = $self->prompt("\nAdministartive e-mail address for status/error information?",
					 $O{admin_email});
	$O{smtp_server} = $self->prompt("\nOutgoing e-mail SMTP server?",
					 $O{smtp_server});

	if(! -f '_mwas'){
	    mkdir('_mwas') or die "ERROR: Cannot save configuration for Build\n".
		"Please check you user prvleges and/or drive space.\n";
	}
	GI::generate_control_files('_mwas','server', \%O);
	GI::generate_control_files('_mwas','all', \%O);
	GI::generate_control_files('_mwas','menus', \%O);

	$self->log_info("\n==============================================================================\n".
			"Finished configuring MWAS settings. Following installation you can change these\n".
			"settings and others by editing the MWAS configuration files:\n".
			"     --> $c_dir/server.ctl\n".
			"     --> $c_dir/menus.ctl\n".
			"     --> $c_dir/maker_opts.ctl\n".
			"     --> $c_dir/maker_bopts.ctl\n".
			"     --> $c_dir/maker_exe.ctl\n".
			"\n==============================================================================\n\n");
    }
}

#commit current MAKER to subversion repository
sub ACTION_commit {
    my $self = shift;

    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update', '');
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('commit');

    #there were changes so re-run install
    if($s_svn != $f_svn){
	$self->dispatch('clean');
	$self->dispatch('install');
    }
}

#update current MAKER from the MAKER repository
sub ACTION_update {
    my $self = shift;

    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update');
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;

    #there were changes so re-run install
    if($s_svn != $f_svn){
	$self->dispatch('clean');
	$self->dispatch('install');
    }

    print "\nSVN STATUS:\n";
    $self->svn_w_args('status');
}

#syncronize the maker/src/bin and maker/src/inc/bin directories
#to maker/bin because of user edits to scripts
sub ACTION_sync {
    my $self = shift;

    $self->sync_bins();
}

#creates a versioned MAKER release
sub ACTION_release {
    my $self = shift;

    #update current MAKER
    print "\nUpdating to most current repository...\n";
    $self->sync_bins();
    my ($s_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    $self->svn_w_args('update', '');

    #doing
    print "\nPre-release commit of any user changes...\n";
    $self->svn_w_args('commit', '-m "pre-release commit"');
    $self->svn_w_args('update', '');

    #clean and check versions for release
    $self->dispatch('clean');
    my $ver = $self->check_update_version(); #returns MAKER version
    $self->{properties}->{dist_version} = $ver;

    File::Which::which('tar') || die "ERROR: Cannot find tar to build the release\n";
    File::Which::which('svn') || die "ERROR: Cannot find the executable svn\n";
    
    #build tarball for users to download
    my $cwd = $self->base_dir;
    my $tgz = "$cwd/maker-$ver.tgz";    
    if(! -f $tgz){
	my ($dir, $base) = $cwd =~ /^(.*\/)([^\/]+)\/src$/;
	
	my $exclude = `cd $dir; svn status $base`;
	$exclude = join("\n", ($exclude =~ /\?\s+([^\n]+)/g), "maker/src/maker-$ver.tgz") ."\n";
	open(OUT, "> .exclude~");
	print OUT $exclude;
	close(OUT);
	
	print "\nBuilding tarball for distribution...\n";
	my $command = "tar -C $dir -zcf $tgz $base --exclude \"*~\" --exclude \".svn\" --exclude \"maker-*.tgz\" --exclude-from .exclude~";
	system($command) && unlink($tgz);
	unlink(".exclude~");
	die "ERROR: tarball creation failed\n" if(! -f $tgz);
    }

    #there were changes so re-run install (updates version info in scripts)
    my ($f_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    if($s_svn != $f_svn){
	print "\nNow reinstalling MAKER scripts to reflect version changes...\n";
	sleep 1;
	$self->dispatch('realclean'); #clean up all old files
	$self->create_build_script; #update stored Build script
	$self->dispatch('install');
    }
}

#replacement for Module::Build's ACTION_install
sub ACTION_install {
    my $self = shift;
    
    $self->SUPER::ACTION_install();

    if($self->feature('mwas_support')){
	require GI;
	require MWAS_util;
	
	my $c_dir = MWAS_util::config_loc(); #configuration file directory
	my @files = ("_mwas/maker_opts.ctl",
		     "_mwas/maker_bopts.ctl",
		     "_mwas/maker_exe.ctl",
		     "_mwas/server.ctl",
		     "_mwas/menus.ctl",
		     );

	my @exists = grep {-f $_} @files;
	my %O = GI::load_server_files(\@exists);
	
	if(! -e $c_dir){
	    mkdir($c_dir) or die "ERROR: Could not create directory $c_dir\n".
		"Please check your user permisions and/or drive space.\n";
	}

	GI::generate_control_files($c_dir,'server', \%O);
	GI::generate_control_files($c_dir,'all', \%O);
	GI::generate_control_files($c_dir,'menus', \%O);    

	$self->log_info("Installing MWAS database and files\n");
	MWAS_util::mwas_setup(\%O);
	$self->log_info("Installing Apollo - Java Web Start\n");
	MWAS_util::apollo_setup(\%O);
	$self->log_info("Setting GBrowse Configuration\n");
	MWAS_util::gbrowse_setup(\%O);
	$self->log_info("Setting JBrowse Configuration\n");
	MWAS_util::jbrowse_setup(\%O);

	print "\n\nMWAS must be started using the following command\n".
	      "before jobs submitted to the server will run:\n".
	      "     sudo maker -MWAS START\n\n";
    }

    $self->maker_status();
}

#replaces Module::Build's ACTION_installdeps
#so MAKER prereqs install locally inside of maker/perl/lib
sub ACTION_installdeps{
    my $self = shift;

    my $prereq = $self->prereq_failures();
    if($prereq && $prereq->{requires}){
	my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
	my $access = 1 if(-w $Config{installsitelib} && -w $Config{installarchlib});
	if(!$access){
	    my ($root_id, $grp_id) = (getpwnam('root'))[2,3];
	    $access = 1 if($< == $root_id || $> == $root_id); 
	}

	my $local;
	if(! $access && ! grep {/Bio\:\:Graphics\:\:Browser2/} @perl){
	    $local = $self->y_n("You do not have write access to install missing Modules.\n".
				"I can try and install these locally (i.e. only for MAKER)\n".
				"in the .../maker/perl/lib directory, or you can run\n".
				"'./Build installdeps' as root or using sudo and try again.\n".
				"Do want MAKER to try and build a local installation?", 'N');
	}	

	die "\n\nERROR: You do not have write access to install missing Modules.\n".
	    "Please run './Build installdeps' as root or using sudo.\n\n" if(!$access && !$local);

	foreach my $m (keys %{$prereq->{build_requires}}){    
	    $self->cpan_install($m, $local);
	}
	foreach my $m (keys %{$prereq->{requires}}){    
	    $self->cpan_install($m, $local);
	}
	
	print "\nRechecking dependencies to see if installation was successful\n";
	$self->check_prereq;

	$self->maker_status;
	
	if($self->prereq_failures()){
	    my ($usr_id) = (getpwnam('root'))[2];
	    print "WARNING: Installation failed (please review any previous errors).\n";
	    print "Try installing the missing packages as 'root' or using sudo.\n" if($< != $usr_id);
	    print "You may need to configure and install these packages manually.\n";
	    return 0;
	}
    }
}

#these install individual external algorithms
sub ACTION_repeatmasker{ shift->_exe_action('RepeatMasker'); }
sub ACTION_blast{ shift->_exe_action('blast'); }
sub ACTION_exonerate{ shift->_exe_action('exonerate'); }
sub ACTION_snap{ shift->_exe_action('snap'); }
sub ACTION_augustus{ shift->_exe_action('augustus'); }
sub ACTION_jbrowse{ shift->_exe_action('jbrowse'); }
sub ACTION_apollo{ shift->_exe_action('apollo'); }
sub ACTION_mpich2{ shift->_exe_action('mpich2', 'mpich2version'); }
sub ACTION_gbrowse{
    my $self = shift;

    if($self->check_installed_status('Bio::Graphics::Browser2', '0')->{ok}){
        my $go = $self->y_n("WARNING: GBrowse was already found on this system.\n".
                            "Do you still want MAKER to install GBrowse for you?", 'N');

	$self->cpan_install('Bio::Graphics::Browser2', 0) if($go);
    }
    else{
	$self->cpan_install('Bio::Graphics::Browser2', 0);
    }
}


#runs all the algorithm installs that are missing
sub ACTION_installexes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();

    return if(! $exe_failures || ! $exe_failures->{exe_requires});
    foreach my $name (grep {!/jbrowse|apollo/i} keys %{$exe_failures->{exe_requires}}){
	$name = lc($name); #all dispatch parameters are lower case
	$self->dispatch($name);
    }

    foreach my $name (grep {/jbrowse|apollo/i} keys %{$exe_failures->{exe_requires}}){
	$name = lc($name); #all dispatch parameters are lower case
	$self->dispatch($name);
    }
    
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

#checks external libraries to see if they're present. Anologous to check_prereqs
sub check_libs{
    my $self = shift;

    my $lib_failures = $self->lib_failures();
    return if(! $lib_failures);

    print "Checking external library dependencies...\n";
    while(my $cat = each %$lib_failures){
	my $s = ($cat eq 'lib_requires') ? '!' : '*';
	$cat =~ /exe_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$lib_failures->{$cat}}){
	    print "    $s  ".$lib_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the libraries\n".
	"indicated above before proceeding with this installation\n".
	"\nMAKER cannot do this for you.  You will have to do it manually.\n";
}

#returns missing exes, anologous to prereq_failures
{
my $set_PATH = 0;
sub exe_failures {
    my $self = shift;
    my $other = shift || {};
    my %exes = (%{$self->exe_requires}, %{$other});

    if(! $set_PATH && $self->config_data('PATH')){
	$ENV{PATH} = ($ENV{PATH}) ?
	    $ENV{PATH}.":".$self->config_data('PATH') : $self->config_data('PATH');
	$set_PATH = 1;
    }

    my %exe_failures;
    while (my $name = each %exes){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $exes{$name});
	my $base = $self->install_destination('exe');
	my $dest = (-d "$base/$name") ? "$base/$name" : "$base/".lc($name);

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {-f $_} ("$dest/$cmd", "$dest/bin/$cmd", File::Which::which($cmd)));
	    last if(($loc) = grep {-f $_ && -x $_} ($cmd));
	}

        if(! $loc){
	    $exe_failures{'exe_requires'}{$name}{have} = '<none>';
	    $exe_failures{'exe_requires'}{$name}{message} = "$name is not installed";
	    $exe_failures{'exe_requires'}{$name}{need} = $exes{$name};
	}
    }

    return (keys %exe_failures) ? \%exe_failures : undef;
}
}

#returns missing libs, anologous to prereq_failures
sub lib_failures {
    my $self = shift;
    my $other = shift || {};
    my %libs = (%{$self->lib_requires}, %{$other});
    
    #for libpng
    push(@DynaLoader::dl_library_path, '/usr/X11/include', '/usr/X11/include', '/sw/lib', '/sw/include');

    my %lib_failures;
    while (my $name = each %libs){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $libs{$name});

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {-f $_} (DynaLoader::dl_findfile($cmd), $cmd));
	}

        if(! $loc){
	    $lib_failures{'lib_requires'}{$name}{have} = '<none>';
	    $lib_failures{'lib_requires'}{$name}{message} = "$name is not installed";
	    $lib_failures{'lib_requires'}{$name}{need} = $libs{$name};
	}
    }

    return (keys %lib_failures) ? \%lib_failures : undef;
}


#hidden connection entry between ACTION_??? algorithm install
#checks to see if already installed and then run the install method
sub _exe_action{
    my $self = shift;
    my $label = shift;
    my $script = shift;

    my $fail = ($script) ? $self->exe_failures({$label => $script}) : $self->exe_failures();
    my @list = keys %{$fail->{exe_requires}} if($fail->{exe_requires});
    if(grep {/^$label$/i} @list){
	$self->_install_exe($label);
    }
    else{
	my $go = $self->y_n("WARNING: $label was already found on this system.\n".
			    "Do you still want MAKER to install $label for you?", 'N');
	$self->_install_exe($label) if($go);
    }
}

#does actual installation of all external algorithms
sub _install_exe {
    my $self = shift;
    my $exe  = shift;
    my $base = $self->install_destination('exe');
    my $path = "$base/$exe";

    #get OS and architecture
    my %os_ok = (Linux_x86_64  => 1,
		 Linux_i386    => 1,
		 Darwin_i386   => 1,
		 Darwin_x86_64 => 1,
		 src           => 1); 
    my ($OS, $ARC) = (POSIX::uname())[0,4];

    ($OS, $ARC) = ('src', '') if(! $os_ok{"$OS\_$ARC"}); #use source code for unknown architectures

    #add fink paths just in case
    if($OS eq 'Darwin'){
	$ENV{C_INCLUDE_PATH} .= ":" if($ENV{C_INCLUDE_PATH});
	$ENV{LIBRARY_PATH} .= ":" if($ENV{LIBRARY_PATH});
	$ENV{C_INCLUDE_PATH} .= "/sw/include:/usr/X11/include";
	$ENV{LIBRARY_PATH} .= "/sw/lib:/usr/X11/lib";
    }

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
	    die "ERROR could not create directory $base for installing external program dependencies\n";
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

	#RepBase
	my $req = 'RepBase';
	my $go = $self->y_n("\nIf are a registered user of RepBase, then MAKER can\n".
			    "download and install RepBase for RepeatMasker for you.\n".
			    "Do you want to do this?", 'Y');

	if($go){
	    print "\n* NOTE: Register at http://www.girinst.org/\n\n";
	    my $user = $self->prompt("Please enter your username:", '');
	    my $pass;
	    if(!$self->check_installed_status('Term::ReadKey', '0')->{ok}){
		$pass = $self->safe_prompt("Please enter your Password:", '');
	    }
	    else{
		$pass = $self->prompt("Please enter your Password:", '');
	    }
	    chdir($path);
	    $file = "$path/$req.tar.gz"; #file to save to
	    $url = $data->{$req}{"$OS\_$ARC"}; #url to rmblast for OS
	    print "Downloading $req...\n";
	    $self->getstore($url, $file, $user, $pass) or return $self->fail($req, $path);
	    print "Unpacking $exe tarball...\n";
	    $self->extract_archive($file) or return $self->fail($req, $path);
	    push(@unlink, $file);
	}
	else{
	    print "\n\n*You will have to install RepBase manually.\n\n";
	}

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

	#RMBlast
	$req = 'rmblast';
        chdir($path);
	File::Path::rmtree("$path/$req");
	$file = "$path/$req.tar.gz"; #file to save to
        $url = $data->{$req}{"$OS\_$ARC"}; #url to rmblast for OS

	if($OS eq 'Linux'){
	    #glibc 2.5 or greater is required for Linux binary
	    my $ldd = File::Which::which('ldd');
	    my ($ver) = (`$ldd --version` =~ /(\d+\.\d+)/) if($ldd);
	    $url = $data->{$req}{src} if(!$ver || $ver < 2.5);
	}

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
	    $self->do_system("./configure --prefix=$path") or return $self->fail($exe, $path);
	    $self->do_system("make") or return $self->fail($exe, $path);
	    $self->do_system("make install") or return $self->fail($exe, $path);
	    chdir($base);
	}
	else{
	    chdir($base);
	    File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
	}

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
	my ($dir) = grep {-d $_} <augustus*>;
	chdir("$dir/src");
	print "Configuring $exe...\n";
	if($OS eq 'Darwin'){
	    File::Copy::move('Makefile', 'Makefile.bak') or return $self->fail($exe, $path);
	    open(IN, "< Makefile.bak");
	    open(OUT, "> Makefile");
	    while(my $line = <IN>){
		$line =~ s/\-static//g;
		print OUT $line;
	    }
	    close(OUT);
	    close(IN);
	}
	$self->do_system("make") or return $self->fail($exe, $path);
	chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
	return $self->fail($exe, $path) if(! -f "$path/bin/$exe");
    }
    elsif($exe eq 'mpich2'){
	#MPICH2
	&File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
	$self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
	$self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	my ($dir) = grep {-d $_} <mpich2*>;
	print "Configuring $exe...\n";
	chdir($dir);
	my %shared = (Linux  => '--enable-sharedlibs=gcc',
		      Darwin => '--enable-sharedlibs=osx-gcc --disable-f77 --disable-fc',
		      src    => '');
	$self->do_system("./configure --prefix=$path --enable-shared $shared{$OS}") or return $self->fail($exe, $path);
	$self->do_system("make") or return $self->fail($exe, $path);
	$self->do_system("make install") or return $self->fail($exe, $path);
	chdir($base);
	File::Path::rmtree($dir) or return $self->fail($exe, $path);
    }
    elsif($exe eq 'apache ant'){
	#do nothing
    }
    elsif($exe eq 'apollo'){
	#apollo
	if(grep {/Apache Ant/} map {keys %{$_->{exe_requires}}} $self->exe_failures()){
	    die "ERROR: Cannot install Apollo without Apache Ant\n".
		"You will need to install it first manually.\n";
	}

	&File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
	$self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
	$self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	my ($dir) = grep {-d $_} <trunk*>;
	print "Configuring $exe...\n";
	File::Copy::copy("$base/../GMOD/Apollo/gff3.tiers", "$dir/conf/gff3.tiers");
	chdir("$dir/src/java");
	$ENV{APOLLO_ROOT} = "$base/$dir";
	$self->do_system("ant jar") or return $self->fail($exe, $path);
	chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
    }
    elsif($exe eq 'jbrowse'){
	#jbrowse
	if(grep {/LibPNG/} map {keys %{$_->{exe_requires}}} $self->lib_failures()){
	    die "ERROR: Cannot install JBrowse without LibPNG\n".
		"You will need to install it first manually.\n";
	}

	&File::Path::rmtree($path);
	my $file = "$base/$exe.tar.gz"; #file to save to
	my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
	$self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
	$self->extract_archive($file) or return $self->fail($exe, $path);
	push (@unlink, $file);
	my ($dir) = grep {-d $_} <jbrowse*>;
	print "Configuring $exe...\n";
	#File::Copy::copy("$base/../GMOD/JBrowse/genome.css", "$dir/genome.css");
	chdir($dir);
	if($OS eq 'Darwin'){
	    $self->do_system('./configure CXXFLAGS=-I/usr/X11/include LDFLAGS=-L/usr/X11/lib')
		or return $self->fail($exe, $path);
	}
	else{
	    $self->do_system('./configure') or return $self->fail($exe, $path);
	}
	$self->do_system('make') or return $self->fail($exe, $path);
	chdir($base);
	File::Copy::move($dir, $exe) or return $self->fail($exe, $path);
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

#fail/cleanup method for installing exes
sub fail {
    my $self = shift;
    my $exe = shift;
    my $path = shift;

    print "\n\nERROR: Failed installing $exe, now cleaning installation path...\n".
	"You may need to install $exe manually.\n\n";

    File::Path::rmtree($path);
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl, has flag for local install
sub cpan_install {
    my ($self, $desired, $local) = @_;

    unless(File::Which::which('make')){
	die "\n\nERROR: Cannot find 'make' on your system. If this is a Mac you\n".
	    "will need to install xcode developer tools before you can\n".
	    "continue. This can be done from the OS X installation disk.\n\n";
    }

    if($desired eq 'Bio::Graphics::Browser2' && $local){
	print "\nWARNING: Bio::Graphics::Browser2 can only be installed globally so\n".
	    "you may need to be logged in as root or use sudo, otherwise this\n".
	    "installation will probably fail.\n\n";
	$local = 0;
    }

    if(! $local){
	my $loc = $self->module_overide($desired);

	if($loc && $desired eq 'Bio::Graphics::Browser2'){
	    print "\n\nWARNING: There is another version of the module Bio::Graphics::Browser2\n".
		"installed on this machine that will supercede a globally installed version.\n".
		"If you want to install the MWAS interface to MAKER, you will have to dlete the\n".
		"offending module before continuing.\n".
		"Location: $loc\n";
	    
	    die "\nWARNING: You will need to delete $loc before continuing.\n";
	}
	elsif($loc){
	    print "\n\nWARNING: There is another version of the module $desired\n".
		  "installed on this machine that will supercede a globally installed version.\n".
		  "I can only continue by installing a local MAKER-only version. If you want a\n".
		  "global installation of the module, you will have to quit and delete the\n".
		  "offending module.\n".
                  "Location: $loc\n";

	    $local = $self->y_n("Do you want to continue with a local installation?", 'N');

	    die "\nWARNING: You will need to delete $loc before continuing.\n"if(! $local);
	}
    }

    #set up PERL5LIB environmental varable since CPAN doesn't see my 'use lib'
    if($local){
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = $self->base_dir."/../perl/lib:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }

    # Here we use CPAN to actually install the desired module
    require CPAN;

    # Save this because CPAN will chdir all over the place.
    my $cwd = getcwd();

    #update CPAN if needed to avoid other installation issues with prereqs
    CPAN::Shell->install('Bundle::CPAN') if (! $self->check_installed_status('CPAN::HandleConfig', '5.5003')->{ok});
    die "\n\nERROR: Cannot load correct CPAN version or related modules\n".
	"Try running 'install Bundle::CPAN' from CPAN manually as either\n".
	"the root user or using sudo\n\n"
	if (! $self->check_installed_status('CPAN::HandleConfig', '0')->{ok});

    #fix weird Mac compilation errors
    my ($OS, $ARC) = (POSIX::uname())[0,4];
    if($OS eq 'Darwin'){
	#$ENV{CFLAGS} .= ($ENV{CFLAGS}) ? ' -m32' : '-m32';
	#$ENV{CCFLAGS} .= ($ENV{CCFLAGS}) ? ' -m32' : '-m32';
	#$ENV{CPPFLAGS} .= ($ENV{CPPFLAGS}) ? ' -m32' : '-m32';
	#$ENV{LDFLAGS} .= ($ENV{LDFLAGS}) ? ' -m32' : '-m32';
	#$ENV{LDDLFLAGS} .= ($ENV{LDDLFLAGS}) ? ' -m32' : '-m32';
    }

    #set up a non-global local module library for MAKER
    if($local){
	my $base = $self->base_dir;
	CPAN::HandleConfig->load;
	$CPAN::Config->{makepl_arg} = "DESTDIR=$base/../perl/ INSTALLDIRS=site INSTALLSITEMAN1DIR=man INSTALLSITEMAN3DIR=man".
	    " INSTALLSITEARCH=lib INSTALLSITELIB=lib INSTALLSITEBIN=lib/bin INSTALLSITESCRIPT=lib/bin";
	$CPAN::Config->{mbuildpl_arg} = "--install_base $base/../perl/ --installdirs site".
	    " --install_path libdoc=$base/../perl/man --install_path bindoc=$base/../perl/man".
	    " --install_path lib=$base/../perl/lib --install_path  arch=$base/../perl/lib".
	    " --install_path bin=$base/../perl/lib/bin --install_path script=$base/../perl/lib/bin ";
	$CPAN::Config->{prefs_dir} = "$ENV{HOME}/.cpan/prefs" if(! -w $CPAN::Config->{prefs_dir});
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }
    else{
	CPAN::HandleConfig->load;
	$CPAN::Config->{makepl_arg} = "INSTALLDIRS=site";
	$CPAN::Config->{mbuildpl_arg} = "--installdirs site";
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }

    #install YAML if needed to avoid other installation issues with prereqs
    CPAN::Shell->install('YAML') if (! $self->check_installed_status('YAML', '0')->{ok});    

    #CPAN::Shell->expand("Module", $desired)->cpan_version <= 2.16;
    #CPAN::Shell->install($desired);
    CPAN::Shell->force('install', $desired);

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
	my $command;
	my $u = scalar getpwuid($>);
	my $g = scalar getgrgid($));
	if($file =~ /\.gz$|\.tgz$/){
	    $command = "tar -zxm -f $file";
	}
	elsif($file =~ /\.bz2?$|\.tbz2?$/){
	    $command = "tar -jxm -f $file";
	}
	else{
	    $command = "tar -xm -f $file";
	}
	$command .= " --owner $u --group $g" unless((POSIX::uname())[0] =~ /darwin/i);

	return $self->do_system($command); #fast
    }
    else{
	die "ERROR: Archive::Tar required to unpack missing executables.\n".
	    "Try running ./Build installdeps first.\n"
	    if(!$self->check_installed_status('Archive::Tar', '0')->{ok});

	return (Archive::Tar->extract_archive($file)) ? 1 : 0; #slow
    }
}

#downloads files from the internet.  Tries to use wget, then curl,
#and finally LWP::Simple
sub getstore {
    my $self = shift;
    my $url = shift;
    my $file = shift;
    my $user = shift;
    my $pass = shift;

    if(File::Which::which('wget')){ #Linux
	my $command = "wget $url -c -O $file --no-check-certificate";
	$command .= " --user $user --password $pass" if(defined($user) && defined($pass));
	return $self->do_system($command); #gives status and can continue partial
    }
    elsif(File::Which::which('curl')){ #Mac
	my $command = "curl --connect-timeout 30 -f -L $url -o $file";
	$command .= " --user $user:$pass" if(defined($user) && defined($pass));
	my $continue = " -C -";

	#gives status and can continue partial
	my $stat = $self->do_system($command . $continue);
	#just redo if continue fails
	$stat = $self->do_system($command) if(! $stat);
	return $stat;
    }
    else{
	die "ERROR: LWP::Simple required to download missing executables\n".
	    "Try running ./Build installdeps first.\n"
	    if(!$self->check_installed_status('LWP::Simple', '0')->{ok});

	$url =~ s/^([^\:]\;\/\/)/$1\:\/\/$user\:$pass\@/ if(defined($user) && defined($pass));
	return LWP::Simple::getstore($url, $file); #just gets the file with no features
    }
}

#prints a nice status message for package configuration and install
sub maker_status {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();

    my $dist_name = $self->dist_name;
    my $dist_version = $self->dist_version;
    my $mpi = ($self->feature('mpi_support')) ? 'CONFIGURED' : 'DISABLED';
    $mpi = 'MISSING PREREQUISITES' if($self->feature('mpi_support') && scalar(grep {/^MPI/} @exes, @libs));
    my $mwas = ($self->feature('mwas_support')) ? 'READY TO INSTALL' : 'NOT CONFIGURED';
    $mwas = 'INSTALLED' if ($self->check_installed_status('Parallel::Application::MPI', '0')->{ok});
    $mwas = 'MISSING PREREQUISITES' if($self->feature('mwas_support') &&
				       (scalar(grep {!/^MPI/} @libs)||
					scalar(grep {/JBrowse|Apollo/} @exes)||
					scalar(grep {/^(CGI|Mail|Bio\:\:Graphics)/} @perl)));
    my $maker = (-f $self->base_dir."/../bin/maker") ? 'INSTALLED' : 'READY TO INSTALL';
    $maker =  'MISSING PREREQUISITES' if(scalar(grep {!/^(CGI|Mail|Bio\:\:Graphics)/} @perl) ||
					 scalar(grep {!/^MPI|JBrowse|Apollo/} @exes));

    print "\n\n";
    print "==============================================================================\n";
    print "STATUS $dist_name $dist_version\n";
    print "==============================================================================\n";
    print "PERL Dependencies:\t";
    print ((@perl) ? 'MISSING' : 'INSTALLED');
    print"\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @perl) ."\n\n" if(@perl);
    print "External Programs:\t";
    print ((@exes) ? 'MISSING' : 'INSTALLED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @exes) ."\n\n" if(@exes);
    print "External C Libraries:\t";
    print ((@libs) ? 'MISSING' : 'INSTALLED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @libs) ."\n\n" if(@libs);
    print "MPI SUPPORT:\t\t";
    print $mpi;
    print "\n";
    print "MWAS Web Interface:\t";
    print $mwas;
    print "\n";
    print "MAKER:\t\t\t";
    print $maker;
    print "\n";

    print "\n\nImportant Commands:\n".
        "\t./Build installdeps\t\#installs missing PERL dependencies\n".
        "\t./Build installexes\t\#installs all missing external programs\n".
        "\t./Build install\t\t\#installs MAKER\n".
        "\t./Build status\t\t\#Shows this status menu\n\n".
        "Other Commands:\n".
        "\t./Build repeatmasker\t\#installs RepeatMasker (asks for RepBase)\n".
        "\t./Build blast\t\t\#installs BLAST (NCBI BLAST+)\n".
        "\t./Build exonerate\t\#installs Exonerate (v2 on UNIX / v1 on Mac OSX)\n".
        "\t./Build snap\t\t\#installs SNAP\n".
        "\t./Build augustus\t\#installs Augustus\n".
        "\t./Build apollo\t\t\#installs Apollo\n".
        "\t./Build gbrowse\t\t\#installs GBrowse (must be root)\n".
        "\t./Build jbrowse\t\t\#installs JBrowse (MAKER copy, not web accecible)\n".
        "\t./Build mpich2\t\t\#installs MPICH2 (but manual install recommended)\n";
}

#test if there is another version of the module overriding the CPAN install
sub module_overide {
    my $self = shift;
    my $desired = shift;

    my $mod = $desired; #holds expected .pm file name
    $mod =~ s/\:\:/\//g;
    $mod .= '.pm';
    
    my $test=  qq(\@INC = qw($Config{installsitelib}
			     $Config{installsitearch}
			     $Config{installvendorlib}
			     $Config{installvendorarch}
			     $Config{installprivlib}
			     $Config{installarchlib});
		  require $desired;
		  print \$INC{\"$mod\"};
		  );
    
    my $ok = `$^X -e 'eval q{$test} or exit 1'`;
    my $loc = $self->module_loc($desired) if($ok);
    
    return ($loc && $loc ne $ok) ? $loc : undef;
}

#gets the location of a module
sub module_loc {
    my $self = shift;
    my $desired = shift;

    return if(! $desired);

    eval "require $desired"; #loads module into \%INC    

    $desired =~ s/\:\:/\//g;
    $desired .= ".pm";

    return $INC{$desired};
}

sub svn_w_args {
    my $self = shift;
    my $param = shift;
    my $o_args = shift;

    my $svn = File::Which::which("svn");
    if($svn){
	#get message off command line
	$svn .= " $param";
	if(defined($o_args)){
	    $svn .= " $o_args";
	}
	else{
	    foreach my $arg (@{$self->args->{ARGV}}){
		if($arg =~ /[\s\t]/){
		    $arg =~ s/\'/\\\'/g;
		    $arg = "'$arg'" 
		    }
		$svn .= " $arg";
	    }
	}
	$svn .= " ".$self->base_dir."/../";

	$self->do_system($svn);
    }
    else{
	die "ERROR: Cannot find the executable svn (subversion respository tool)\n";
    }
}

sub load_w_o_header {
    my $file = shift;

    my $data;
    open(IN, "< $file");
    while(my $line = <IN>){	
	#strip of perl shebang portiion 
	while($line =~ /^\#!.*perl/ ||
	      $line =~ /exec\s+.*perl/ ||
	      $line =~ /if\s*0\;/ ||
	      $line =~ /^[\s\t\n]+$/
	      ){
	    $line = <IN>;
	}
	$data = join('', $line, <IN>);
    }
    close(IN);

    return $data;
}

sub sync_bins {
    my $self = shift;
    my $cwd = $self->base_dir;

    my $BIN  = "$cwd/../bin";
    my $sbin = "$cwd/bin";
    my $ibin = "$cwd/inc/bin";

    my @ifiles = map {$_ =~ /([^\/]+)$/} grep {-f $_ && !/\~$/ && !/\.PL$/} <$ibin/*>;
    my @sfiles = map {$_ =~ /([^\/]+)$/} grep {-f $_ && !/\~$/ && !/\.PL$/} <$sbin/*>;

    foreach my $file (@ifiles, @sfiles){
	my $bfile = "$BIN/$file";
	my $rfile = (-e "$ibin/$file") ? "$ibin/$file" : "$sbin/$file";

	next if(! -e $bfile);
	#-w permission must be set before these files are edited by the user
	next unless(sprintf("%04o", (stat($bfile))[2] & 07777) =~ /[2367]/);

	my $bdata = load_w_o_header($bfile);
	my $rdata = load_w_o_header($rfile);

	#scripts have been altered by user
	if($bdata ne $rdata){
	    my $bmod = (stat($bfile))[9];
	    my $rmod = (stat($rfile))[9];
	    
	    if($bmod > $rmod){
		print "copying $bfile  -->  $rfile\n";

		#backup incase of failure
		File::Copy::move($rfile, "$rfile.bk~");

		if(open(IN, "> $rfile")){
		    print IN "\#!\\usr\\bin\\perl\n\n";
		    print IN $bdata;
		    close(IN);
		}

		#restore on failure
		if(! -f $rfile && -f "$rfile.bk~"){
		    File::Copy::move("$rfile.bk~", $rfile);
		}
		else{
		    unlink("$rfile.bk~");
		}
	    }
	}	
    }
}

sub check_update_version {
    my $self = shift;

    #get current subversion version
    my ($svn) = `svn info` =~ /Revision\:\s*(\d+)/;
    die "ERROR: Could not query subversion repository\n" if(!$svn);

    #get old version information for last stable release
    open(IN, "< version") or die "ERROR: Could not open MAKER version file\n";
    my $data = join("\n", <IN>);
    my ($old_svn) = $data =~ /\$SVN=(\d+)/;
    my ($old_version) = $data =~ /\$VERSION=([\d\.]+)/;
    close(IN);

    #check if update is really needed
    my $version = $old_version;
    if($old_svn == $svn){
	print "MAKER is already up to date as stable release $version\n";

	return $version;
    }
    else{
	#set new version
	$old_version =~ /(.*)\.(\d+)$/;
	$version = $1;
	my $s = $2; #sub version
	my $n = sprintf ('%02s', $s + 1); #new sub version
	$s = ".$s"; #add decimal
	$n = ".$n"; #add decimal
	
	#if version iteration results in lower value then make sub iterator
	#this means major version numbers can only be changed by the user
	if($n < $s){
	    $n = "$s.01";
	}
	$version .= $n;

	#output what will be next version to file
	#then another commit will be performed to
	#sync subverion with the release
	my $commit_svn = $svn;
	do{
	    $svn = $commit_svn;
	    $svn++;
	    open(OUT, "> version");
	    print OUT "\$VERSION=$version\n";
	    print OUT "\$SVN=$svn\n";
	    close(OUT);

	    #files to fix version for
	    my $cwd = $self->base_dir;
	    my @files = ("$cwd/bin/maker",
			 "$cwd/bin/evaluator",
			 "$cwd/bin/iprscan_wrap",
			 "$cwd/inc/bin/mpi_evaluator",
			 "$cwd/inc/bin/mpi_iprscan",
			 "$cwd/../lib/GI.pm"
			 );

	    #changing script version here
	    foreach my $file (@files){
		open(IN, "< $file");
		unlink($file);
		open(OUT, "> $file");
		while(my $line = <IN>){
		    $line =~ s/\$VERSION\s*\=\s*\'[\d\.]+\'/\$VERSION = \'$version\'/;
		    print OUT $line;
		}
		close(OUT);
		close(IN);
	    }

	    $self->svn_w_args('commit', "-m \"MAKER stable release version $version\"");
	    $self->svn_w_args('update', '');
	    ($commit_svn) = `svn info` =~ /Revision\:\s*(\d+)/;
	    die "ERROR: Could not query subversion repository\n" if(!$commit_svn);
	}while($svn != $commit_svn);

	print "MAKER has been updated to stable release $version\n";

	return $version;
    }
}

sub safe_prompt {
    require Term::ReadKey;

    my $self = shift;
    my $m = shift;
    my $d = shift || '';

    print "$m [$d ]";

    my $key = 0;
    my $r = "";
    #Start reading the keys
    ReadMode(4); #Disable the control keys (raw mode)

    while(ord($key = ReadKey(0)) != 10) { #This will continue until the Enter key is pressed (decimal value of 10)
	if(ord($key) == 127 || ord($key) == 8) { #DEL/Backspace was pressed
	    if(length($r) > 0){
		#1. Remove the last char from the password
		chop($r);
		#2 move the cursor back by one, print a blank character, move the cursor back by one
		print "\b \b";
	    }
	} elsif(ord($key) < 32) {
	    # Do nothing with these control characters
	} else {
	    $r = $r.$key;
	    print "*";
	}
    }
    print "\n"; #because the user pressed enter
    ReadMode(0); #Reset the terminal once we are done

    return $r; #Return the response
}

sub set_gmod_locs {
    my $self = shift;
    my $type = shift || 'ALL';

    if($type eq 'ALL' || $type eq 'APOLLO_ROOT'){
	my $APOLLO_ROOT = "$BIN/../exe/apollo";
	$APOLLO_ROOT = $ENV{APOLLO_ROOT} if(! -d $APOLLO_ROOT && $ENV{APOLLO_ROOT} && -d $ENV{APOLLO_ROOT});
	$APOLLO_ROOT = (File::Which::which('apollo') || '') =~ /^([^\n]+)\/bin\/apollo\n?$/ if(!-d $APOLLO_ROOT);
	$APOLLO_ROOT = '/usr/local/gmod/apollo' if(! -d $APOLLO_ROOT);
	$APOLLO_ROOT = '' if(! -d $APOLLO_ROOT);
	$APOLLO_ROOT = $self->prompt("\nPlease provide the base directory for your Apollo installation?",
				      $APOLLO_ROOT);
	$self->config_data(APOLLO_ROOT => $APOLLO_ROOT);
    }

    if($type eq 'ALL' || $type eq 'JBROWSE_ROOT'){
	my $JBROWSE_ROOT ="$BIN/../exe/jbrowse";
	$JBROWSE_ROOT = '/var/www/html/jbrowse' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = '/Library/WebServer/Documents/jbrowse' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = '/var/www/jbrowse' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = '/usr/local/gmod/jbrowse' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = '/data/var/www/jbrowse' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = '' if(! -d $JBROWSE_ROOT);
	$JBROWSE_ROOT = $self->prompt("\nPlease provide the base directory for your JBrowse installation?",
				       $JBROWSE_ROOT);
	$self->config_data(JBROWSE_ROOT => $JBROWSE_ROOT);
    }
    
    if($type eq 'ALL' || $type eq 'GBROWSE_MASTER'){
	my $GBROWSE_MASTER = '/etc/gbrowse/GBrowse.conf';
	$GBROWSE_MASTER = '/etc/gbrowse2/GBrowse.conf' if(! -f $GBROWSE_MASTER);
	$GBROWSE_MASTER = '' if(! -f $GBROWSE_MASTER);

	$GBROWSE_MASTER = $self->prompt("\nPlease provide the path to your GBrowse.conf file?",
					 $GBROWSE_MASTER);
	$self->config_data(GBROWSE_MASTER => $GBROWSE_MASTER);
    }
}

1;
