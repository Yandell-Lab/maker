package MWAS_util;

use strict;
use warnings;

use FindBin;
use GI;
use File::NFSLock;
use Data::Dumper;
use File::Copy;
use File::Glob;

=head1
#-----------------------------------------------------------------------------
#this method collects values pointed to by menu/file names for a type
#and returns an array reference
sub get_fileValue_from_fileName {
   my $dbh = shift;
   my $name = shift;
   my $type = shift;
   my $user_id = shift;

   die "ERROR: Invalid source value '$source'" if($source !~ /^all$|^server$|^user$/);

   #get value from user file options
   my ($val) = $self->dbh->selectrow_arrayref(qq{SELECT value FROM files WHERE name='$name'}.
					      qq{ AND type='$type'AND iduser='$id'}
					     ) unless($source eq 'server' || ! $id);

   #user values always override system values
   return $val if($val);

   #get value from global menu options
   ($val) = $self->dbh->selectrow_array(qq{SELECT value FROM menus WHERE name='$name'}.
					qq{ AND type='$type'}
				       ) unless($source eq 'user');

   return $val;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_server_default_options {
   my $dbh = shift;

   my $def_opt = $self->dbh->selectrow_hashref(qq{SELECT * FROM all_default_opt});

   return $def_opt;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_user_from_login {
   my $dbh = shift;

   my $username = $self->session->param('AUTH_USERNAME');
   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE login='$username'});

   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_user_id_from_login {
   my $self = shift;

   return $self->get_user_info()->{user_id};
}
#-----------------------------------------------------------------------------
sub get_files_from_user_id{
   my $self = shift;
   my $id = $self->get_user_id();

   my $files = $self->dbh->selectall_arrayref(qq{SELECT * FROM files WHERE user_id='$id'},
					      {Slice => {}}
					     );

   return $files;
}
=cut
#-----------------------------------------------------------------------------
sub config_loc {
    return '/etc/mwas';
}
#-----------------------------------------------------------------------------
#setup method for MWAS server and files
sub mwas_setup {
    my %CTL_OPT = %{shift @_};

    my $c_dir = config_loc();
    my $cgi_dir  = $CTL_OPT{cgi_dir};  #get directory to build CGI files
    my $html_dir = $CTL_OPT{html_dir}; #get directory to build HTML files
    my $data_dir = $CTL_OPT{data_dir}; #get directory for storing data

    #make directories if they don't exist
    foreach my $dir ($c_dir, $cgi_dir, $html_dir, $data_dir){
        my $base = $dir;
        $base =~ s/[^\/]+\/?$//;
        if(-d $base && ! -d $dir){
            mkdir($dir) or warn "ERROR: Could not create directory. Do you have\n".
                "sufficient privileges to do this?\n\n";
        }
    }

    die "ERROR: You must supply a directory for cgi_dir\n".
        "to install CGI files\n" if(! $cgi_dir);
    die "ERROR: You must supply a directory for html_dir\n".
        "to install HTML files\n" if(! $html_dir);
    die "ERROR: You must supply a directory for data_dir\n".
        "to store user files and jobs\n" if(! $data_dir);

    die "ERROR: The directory '$cgi_dir' does not exist\n" if(! -d $cgi_dir);
    die "ERROR: The directory '$html_dir' does not exist\n" if(! -d $html_dir);
    die "ERROR: The directory '$data_dir' does not exist\n" if(! -d $data_dir);

    #set up correct ownership
    system("chown -R :$CTL_OPT{apache_user} $cgi_dir") &&
        die "ERROR: Could not establish ownership for $cgi_dir\n".
        "You may need to log in as 'root' for setup\n\n";

    system("chown -R :$CTL_OPT{apache_user} $html_dir") &&
        die "ERROR: Could not establish ownership for $html_dir\n".
	"You may need to log in as 'root' for setup\n\n";

    system("chown :$CTL_OPT{apache_user} $data_dir") &&
        die "ERROR: Could not establish ownership for $data_dir\n".
        "You may need to log in as 'root' for setup\n\n";

    #set directory guid to preserve permissions on new files in directory
    system("chmod 2775 $cgi_dir") &&
        die "ERROR: Could not establish permissions for $cgi_dir\n".
        "You may need to log in as 'root' for setup\n\n";

    system("chmod 2775 $html_dir") &&
        die "ERROR: Could not establish permissions for $html_dir\n".
        "You may need to log in as 'root' for setup\n\n";

    system("chmod 2775 $data_dir") &&
        die "ERROR: Could not establish permissions for $data_dir\n".
        "You may need to log in as 'root' for setup\n\n";

    #copy all cgi, html, and support files to a web accessible directory
    (my $b_dir = $FindBin::RealBin) =~ s/\/(src|bin|MWAS\/bin)\/?$/\/MWAS\/bin/;
    my $m_dir = "$b_dir/../../"; #maker package base direcory
    my $co_dir = "$b_dir/../cgi-bin/"; #original direcory
    my $ho_dir = "$b_dir/../html/"; #original direcory
    my $m_lib = "$b_dir/../../lib/"; #/maker/lib
    my $p_lib = "$b_dir/../../perl/lib/"; #/maker/perl/lib
    system("cp", "-R", File::Glob::bsd_glob("$co_dir/*"), "$cgi_dir/") && die("ERROR: Copying files to $cgi_dir failed\n");
    system("cp", "-R", File::Glob::bsd_glob("$ho_dir/*"), "$html_dir/") && die("ERROR: Copying files to $html_dir failed\n");
    system("cp", "-R", File::Glob::bsd_glob("$m_lib/*"), "$cgi_dir/lib/") && die("ERROR: Copying files to $cgi_dir failed\n");
    if(-d $p_lib){ #only exists when MPI is installed
        system("cp", "-R", File::Glob::bsd_glob("$p_lib/*"), "$cgi_dir/lib/") && die("ERROR: Copying files to $cgi_dir failed\n");
    }

    #copy maker executables
    mkdir("$data_dir/maker") if(! -d "$data_dir/maker");
    system("cp", "-R", File::Glob::bsd_glob("$m_dir/*"), "$data_dir/maker/") && die("ERROR: Copying files to $data_dir/maker/ failed\n");

    #recursively set group permission to write for all files in the directory
    system("chmod -R g+w $cgi_dir") &&
        die "ERROR: Could not establish permissions for $cgi_dir\n".
	"You may need to log in as 'root' for setup\n\n";

    system("chmod -R g+w $html_dir") &&
        die "ERROR: Could not establish permissions for $html_dir\n".
	"You may need to log in as 'root' for setup\n\n";

    system("chmod -R g+w $data_dir/maker") &&
        die "ERROR: Could not establish permissions for $data_dir\n".
	"You may need to log in as 'root' for setup\n\n";

    #set correct exe location for moved files in maker/exe
    my %E = GI::parse_ctl_files([config_loc()."/maker_exe.ctl"]);
    while(my $key = each %E){
	$E{$key} = GI::s_abs_path($E{$key}) if($E{$key});
	if($E{$key} =~ /\/exe\/([^\/]+\/[^\/]+|[^\/]+\/bin\/[^\/]+)$/ &&
	   -e "$data_dir/maker/exe/$1"
	  ){
	    $E{$key} = "$data_dir/maker/exe/$1";

	    if($key eq 'RepeatMasker'){
		#Configure RepeatMasker
		my $cwd = Cwd::cwd;
		my $path = "$data_dir/maker/exe/RepeatMasker";
		chdir($path);
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
		system("$^X ./configure < $tmp &> /dev/null") &&
		    die "ERROR: Could not configure RepeatMasker for its new location\n";
		chdir($cwd);
	    }
	}
    }

    if($CTL_OPT{JBROWSE_ROOT} =~ /\/exe\/jbrowse$/ &&
           -e "$data_dir/maker/exe/jbrowse"
       ){
	$CTL_OPT{JBROWSE_ROOT} = "$data_dir/maker/exe/jbrowse";
    }

    #fixed values
    %CTL_OPT = (%CTL_OPT, %E);

    #generate the control files stored with website executables
    GI::generate_control_files(config_loc(),'all', \%CTL_OPT);
    GI::generate_control_files(config_loc(),'server', \%CTL_OPT);
    GI::generate_control_files(config_loc(),'menus', \%CTL_OPT);

    #build redirection page (so user can see directories if apache is not configured)
    open(OUT, "> $html_dir/index.html");
    my $url = ($CTL_OPT{cgi_web} =~ /\:\/\//) ?  $CTL_OPT{cgi_web} : "$CTL_OPT{web_address}/$CTL_OPT{cgi_web}";
    $url .= "/maker.cgi";
    $url =~ s/([^:])\/\//$1\//g;
    print OUT "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n".
	      "<HTML>\n".
	      "<meta http-equiv=\"REFRESH\" content=\"0;url=$url\">\n".
              "</HTML>";
    close(OUT);

    #DB error checking
    if($CTL_OPT{DBI} && $CTL_OPT{dbname}){
        #try to connect to the database as a test
        my $dsn = "DBI:$CTL_OPT{DBI}:dbname=";
        $dsn .= "$CTL_OPT{data_dir}/" if($CTL_OPT{DBI} eq 'SQLite');
        $dsn .= "$CTL_OPT{dbname};";
        $dsn .= "host=$CTL_OPT{host};" if($CTL_OPT{host});
        $dsn .= "port=$CTL_OPT{port};" if($CTL_OPT{host} && $CTL_OPT{port});
	if(!$CTL_OPT{username}){ #warn if no db_username specified
            warn "WARNING: No username for connecting to the database is specified.\n".
                "Without a username, online requests to MWAS will fail because of\n".
                "permission issues. Please make the necessary changes to server.ctl.\n\n";
        }

        my $dbh = DBI->connect($dsn, $CTL_OPT{username}, $CTL_OPT{password}, {RaiseError => 1})
            or die "<br>Got error $DBI::errstr when connecting to database\n";

        #fix permissions on SQLite database
        if($CTL_OPT{DBI} eq 'SQLite'){
            my $db = "$CTL_OPT{data_dir}/$CTL_OPT{dbname}";
            unless((sprintf "%04o", (stat $db)[2] & 07777)  == 0774){
                system("chmod 774 $db") &&
                    die "ERROR: Could not establish permissions for $db\n".
                    "You may need to log in as 'root' for setup\n\n";
            }
        }

	#setup the database
        check_table_structure($dbh, \%CTL_OPT);
        reset_queue($dbh, $CTL_OPT{data_dir});
        $dbh->disconnect;

        #set up correct ownership
        system("chown -R :$CTL_OPT{apache_user} $cgi_dir") &&
            die "ERROR: Could not establish ownership for $cgi_dir\n".
            "You may need to log in as 'root' for setup\n\n";

        system("chown -R :$CTL_OPT{apache_user} $html_dir") &&
            die "ERROR: Could not establish ownership for $html_dir\n".
            "You may need to log in as 'root' for setup\n\n";

        system("chown -R :$CTL_OPT{apache_user} $data_dir/maker") &&
            die "ERROR: Could not establish ownership for $data_dir/maker\n".
            "You may need to log in as 'root' for setup\n\n";

        system("chown -R :$CTL_OPT{apache_user} $data_dir/maker") &&
            die "ERROR: Could not establish ownership for $data_dir/maker\n".
            "You may need to log in as 'root' for setup\n\n";
    }
    else{
        die "ERROR: Could not setup the database structure because no database\n".
            "or database interface is specified in server.ctl.  Please finish\n".
            "editing the control files and run SETUP again.\n\n";
    }
}
#-----------------------------------------------------------------------------
sub apollo_setup {
    my %CTL_OPT = %{shift @_};

    #go the extra mile and try and setup apollo webstart
    if($CTL_OPT{APOLLO_ROOT}){
	my $cgi_dir  = $CTL_OPT{cgi_dir};
        my $html_dir = $CTL_OPT{html_dir};
	my $data_dir = $CTL_OPT{data_dir};

        #set up directory for web
        my (undef, $t_xml) = File::Temp::tempfile();
        mkdir("$html_dir/jars") if(! -d "$html_dir/jars");

        #create generate Apollo webstart
	(my $b_dir = $FindBin::RealBin) =~ s/\/(src|bin|MWAS\/bin)\/?$/\/MWAS\/bin/;
        my $tt_xml = "$b_dir/../../GMOD/Apollo/apollo_webstart.xml.tt";

	#webstart_generator.pl can't handle space in directory names
	if($tt_xml =~ / /){
	    unlink("./.apollo_tt_xml");
	    symlink($tt_xml, "./.apollo_tt_xml");
	    $tt_xml = "./.apollo_tt_xml";
	}

	my $jars = "$CTL_OPT{APOLLO_ROOT}/jars";
	if($jars =~ / /){
	    unlink("./.apollo_jars");
            symlink($jars, "./.apollo_jars");
            $jars = "./.apollo_jars";
	}

	my $c_dir = config_loc();
        $ENV{APOLLO_ROOT} = $CTL_OPT{APOLLO_ROOT};
	$ENV{PERL5LIB} = ($CTL_OPT{PERL5LIB}) ?
	    $CTL_OPT{PERL5LIB}.":$b_dir/../../perl/lib" : ":$b_dir/../../perl/lib";

        system("$CTL_OPT{APOLLO_ROOT}/bin/webstart_generator.pl","-i","$tt_xml","-d",
               $jars,"-o","$c_dir/apollo.jnlp","-D","$html_dir/jars") &&
               die "ERROR: Generating Apollo webstart jars and jnlp file failed\n";

	unlink("./.apollo_tt_xml");
	unlink("./.apollo_jars");
    }
    else{
        die "ERROR: You must suply a value for APOLLO_ROOT in server.ctl to\n".
            "setup Apollo for use with the MAKER web interface\n\n";
    }
}
#-----------------------------------------------------------------------------
sub gbrowse_setup {
    my %CTL_OPT= %{shift @_};

    #go the extra mile and try and configure GBROWSE
    if($CTL_OPT{GBROWSE_MASTER}){
        my $old = $CTL_OPT{GBROWSE_MASTER};
        my $new = "$CTL_OPT{GBROWSE_MASTER}.new";

        die "ERROR: The file specified in GBROWSE_MASTER '$old' does not exist\n"
            if(! -e $old);

        #make a new configuration file
        open(IN, "< $old") or die "ERROR: Could not open $old\n";
        open(OUT, "> $new") or die "ERROR: Could not open $new\n";

        my $section = ''; #keeps track of what section in the config file I am in
        while(my $line = <IN>){
            #identify section header
            if($line =~ /^[^\#\=]*\[([^\]]+)\][^\=]*$/){
                $section = $1;
            }

            print OUT $line unless($section =~ /^\=\~MWAS_/); #skip MAKER section
        }

        close(IN);

        #now add the MAKER section to GBrowse config file
        print OUT "\n[=~MWAS_(\\d+)_(\\d+)]\n";
        print OUT "description = MAKER Web Annotation Service\n";
        print OUT "path        = $CTL_OPT{cgi_dir}/stream.cgi type gbrowse user_id \$1 job_id \$2 |\n";
        close(OUT);

        #now replace old file
        system("mv $new $old") && die "ERROR: Could not replace existing GBrowse configuration file\n";

	#place configuration file
	my $cgi_dir  = $CTL_OPT{cgi_dir};
        my $html_dir = $CTL_OPT{html_dir};
	my $data_dir = $CTL_OPT{data_dir};

	(my $b_dir = $FindBin::RealBin) =~ s/\/(src|bin|MWAS\/bin)\/?$/\/MWAS\/bin/;
        my $tt_conf = "$b_dir/../../GMOD/GBrowse/gbrowse.conf.tt";
	my $c_dir = config_loc();

	if(! -f "$c_dir/gbrowse.conf.tt"){
	    system("cp", "$tt_conf", "$c_dir/")
		&& die "ERROR: Preparing GBrowse configuration file failed\n";
	}
    }
    else{
        die "ERROR: You must suply a value for GBROWSE_MASTER in server.ctl to\n".
            "setup GBROWSE for use with the MAKER web interface\n\n";
    }
}
#-----------------------------------------------------------------------------
sub jbrowse_setup {
    my %CTL_OPT= %{shift @_};

    #configure GBROWSE
    if($CTL_OPT{JBROWSE_ROOT}){
        #copy JBrowse configuration file
	my $cgi_dir  = $CTL_OPT{cgi_dir};
	my $data_dir = $CTL_OPT{data_dir};
        (my $b_dir = $FindBin::RealBin) =~ s/\/(src|bin|MWAS\/bin)\/?$/\/MWAS\/bin/;
        my $c_dir = config_loc();
	if(! -f "$c_dir/maker.css"){
	    system("cp", "$b_dir/../../GMOD/JBrowse/maker.css", "$c_dir/")
		&& die "ERROR: Preparing JBrowse configuration file failed\n";
	}
    }
    else{
        die "ERROR: You must suply a value for JBROWSE_ROOT in server.ctl to\n".
            "setup JBROWSE for use with the MAKER web interface\n\n";
    }
}
#-----------------------------------------------------------------------------
#check to see if a job with the exact same control file options exists
#and returns its id
sub package_already_exists{
    my $dbh = shift @_;
    my %CTL_OPT = %{shift @_};
    my $user_id = shift @_;

    my %def_opt = (GI::set_defaults('opts'), GI::set_defaults('bopts')); #get system produced CTL_OPT

    my @to_map =  grep {defined $CTL_OPT{$_} && !/gmhmm_e|gmhmm_p/i} keys %def_opt;
    my @set = map {"ctl_opt.".lc($_)." \= '".$CTL_OPT{$_}."'" } @to_map;

    my $dsn = "SELECT ctl_opt.job_id FROM ctl_opt JOIN jobs ON ctl_opt.job_id=jobs.job_id ".
	      "WHERE ".join(' AND ', @set)." AND jobs.is_packaged=1";

    my ($job_id) = $dbh->selectrow_array($dsn);

    return $job_id;
}
#-----------------------------------------------------------------------------
#copies a finhed job_package to another job_id
sub copy_package{
    my $dbh = shift;
    my $job_old = shift;
    my $job_new = shift;
    my $data_dir = get_data_dir($dbh);

    #get job run directory
    my $job_dir = "$data_dir/jobs/$job_old/";
    my $new_dir = "$data_dir/jobs/$job_new/";

    #get new result directory
    my $r_dir = "$job_new.maker.output";

    #copy and rename files
    File::Path::mkpath($new_dir) if(! -d $new_dir);

    my @files =  <$job_dir/*>;

    @files = grep {!/\.tar\.gz$/} @files;

    system("cp", "-R", @files, "$new_dir/");
    @files = (<$new_dir/*/$job_old*>,<$new_dir/$job_old*>);
    foreach my $f (@files){
	my $new = $f;
	$new =~ s/\/$job_old([\.\_][^\/]+)\/?$/\/$job_new$1/;
	move($f, $new);
    }

    #fix log contents on copy
    open(my $IN, "< $job_dir/$job_old.maker.output/$job_old\_master_datastore_index.log");
    open(my $OUT, "> $new_dir/$job_new.maker.output/$job_new\_master_datastore_index.log");
    while(my $line = <$IN>){
	my @F = split("\t", $line);
	$F[1] =~ s/^$job_old\_/$job_new\_/;
	print $OUT join("\t", @F);
    }
    close($OUT);
    close($IN);

    #re-tar everything
    system("cd $new_dir\n".
	   "tar --exclude \"run.log\" --exclude \"theVoid\*\" --exclude \"seen.dbm\"".
	   " --exclude \"mpi_blastdb\" -zcf $r_dir.tar.gz $r_dir") &&
	   die("ERROR: Building tarball for job '$job_new' failed\n");
}

#-----------------------------------------------------------------------------
sub get_data_dir {
    my $dbh = shift;

    my ($data_dir) = $dbh->selectrow_array("SELECT data_dir FROM all_default_opt");

    return $data_dir;
}
#-----------------------------------------------------------------------------
sub get_length_for_value {
    my $dbh = shift || die "ERROR: No dbh provided\n";
    my $value = shift;

    #get value from user file options
    my ($len) = $dbh->selectrow_array(qq{SELECT length FROM files WHERE value='$value'});

    #user values always override system values
    return $len if($len);

    #get value from global menu options
    ($len) = $dbh->selectrow_array(qq{SELECT length FROM menus WHERE value='$value'});

    return $len;
}
#-----------------------------------------------------------------------------
#get length of a fasta file in base pairs
sub fasta_length{
    my $file = shift;
    my $total = 0;

    open(my $IN, "< $file");
    while(my $line = <$IN>){
	next if($line =~ /^>/);
	chomp $line;
	$total += length($line);
    }
    close($IN);

    return $total;
}
#-----------------------------------------------------------------------------
#sends messages via e-mail i.e. errors or status reports
sub send_message{
    my $address = shift @_;
    my $smtp    = shift @_;
    my $subject = shift @_;
    my $msg     = shift @_;

    my $sender = Mail::Sender->new({smtp => $smtp});

    my $sq = $sender->MailMsg({to      => $address,
                               from    => "no-reply\@$smtp",
                               subject => $subject,
                               msg     => $msg
			       });
}
#-----------------------------------------------------------------------------
sub date_time {
    #get time values
    my ($sec, $min, $hour, $mday, $month,
        $year, $wday, $yday, $iddst) = localtime(time);

    #fix values
    $year += 1900;
    $month++;

    #fix digit spacing
    foreach my $v ($month, $mday, $hour, $min){
        if ($v < 10) {
            $v = "0$v";
        }
    }

    return "$month/$mday/$year $hour:$min";
}
#-----------------------------------------------------------------------------
#lock for updating the database
sub lockDB {
    my $dir = shift;
    return File::NFSLock->new("$dir/.dblock", 'EX', 300, 300);
}
#-----------------------------------------------------------------------------
#standard message for not getting the lock
sub lock_error {
    print STDERR "ERROR: Could not get lock on database:\t".date_time()."\n";
}
#-----------------------------------------------------------------------------
#standard message for not getting the lock
sub dbh_do_commit {
    my $dbh = shift || return;
    my $do_string = shift || return;
    my $dir = shift; #for lock

    #my $lock = lockDB($dir) || ((print STDERR lock_error()) && (return));
    $dbh->do($do_string);
    #$dbh->commit;
    #$lock->unlock;

    return 1;
}
#------------------------------------------------------------------------
#will verify the existance of all tables and columns and create them if necessary
sub check_table_structure {
    my $dbh = shift @_;
    my %O = %{shift @_};

    #checking if CTL_OPT is loaded
    my %def_opt = (GI::set_defaults('opts'), GI::set_defaults('bopts')); #get system produced CTL_OPT
    die "ERROR: CTL_OPT must be loaded before calling main::check_table_structure\n"
        if(@{[keys %def_opt]} > @{[keys %O]});

    my @tables = $dbh->tables;

    #list of expected tables
    my @expected = qw(users
                      menus
                      files
                      jobs
                      ctl_opt
                      all_default_opt
                      stat
                      id_store
		      );

    #if an expected table is not found, create it
    ####===lock the database so no changes are made until all tables are in place
    #my $lock = lockDB($O{data_dir}) || die lock_error();

    #use different datatypes depending on type of SQL used in DBI
    my $t_type = ($dbh->{Driver}->{Name} =~ /mysql/i) ? 'VARCHAR(100)' : 'TEXT';

    #weird extra tables in Pg
    @tables = grep {/^public\./} @tables if($dbh->{Driver}->{Name} =~ /^pg$/i);

    #create users table
    my $dsn = (qq{CREATE TABLE users (user_id INTEGER PRIMARY KEY, login $t_type UNIQUE, }.
		      qq{password $t_type, first $t_type, last $t_type, e_mail $t_type, }.
		      qq{institution $t_type, is_guest INTEGER, is_admin INTEGER, }.
		      qq{date_created $t_type, last_login $t_type)});
    $dbh->do($dsn) if(! grep {/users/} @tables);
    add_missing_columns($dbh, $dsn, 'users');

    #create files table
    $dsn = (qq{CREATE TABLE files (name $t_type, value $t_type UNIQUE, contig_count $t_type, }.
		   qq{type $t_type, length $t_type, user_id INTEGER, description $t_type)});
    $dbh->do($dsn) if(! grep {/files/} @tables);
    add_missing_columns($dbh, $dsn, 'files');

    #create jobs table
    $dsn = (qq{CREATE TABLE jobs (job_id INTEGER PRIMARY KEY, user_id INTEGER, submit_id INTEGER, length $t_type, }.
		   qq{type $t_type, old_id INTEGER, is_queued INTEGER, is_started INTEGER, is_running INTEGER, }.
		   qq{is_finished INTEGER, is_error INTEGER, is_packaged INTEGER, is_saved INTEGER, admin_block INTEGER, }.
		   qq{is_tutorial INTEGER, cpus INTEGER, start_time $t_type, finish_time $t_type, name $t_type, note $t_type)});
    $dbh->do($dsn) if(! grep {/jobs/} @tables);
    add_missing_columns($dbh, $dsn, 'jobs');

    #create ctl_opt table
    my @defaults = map {lc("\"$_\"")." $t_type"} keys %def_opt; #get table column names and datatype
    $dsn = "CREATE TABLE ctl_opt (job_id INTEGER UNIQUE, ".join(', ', @defaults).')';
    $dbh->do($dsn) if(! grep {/ctl_opt/} @tables);
    add_missing_columns($dbh, $dsn, 'ctl_opt');

    #create/re-create menus table (always)
    $dbh->do(qq{DROP TABLE menus}) if(grep {/menus/} @tables);
    $dbh->do(qq{CREATE TABLE menus (name $t_type, value $t_type, contig_count $t_type,}.
		    qq{type $t_type, length $t_type, is_tutorial INTEGER, description $t_type)});

    while(my $type = each %{$O{menus}}){
        while(my $name = each %{$O{menus}{$type}}){
            my $value = $O{menus}{$type}{$name};
            my $is_tutorial = 0;

            if($value =~ /\/dpp_contig.fasta$/ ||
               $value =~ /\/dpp_protein.fasta$/ ||
               $value =~ /\/dpp_est.fasta$/ ||
               $value =~ /\/pyu-contig.fasta$/ ||
               $value =~ /\/pyu-est.fasta$/ ||
               $value =~ /\/pyu-protein.fasta$/ ||
               $value =~ /\/pyu.hmm$/ ||
               $value =~ /\/pyu.mod$/ ||
               $value =~ /\/legacy-contig.fasta$/ ||
               $value =~ /\/legacy-set1.gff$/ ||
               $value =~ /\/legacy-set2.gff$/ ||
               $value =~ /\/pass-contig.fasta$/ ||
               $value =~ /\/pass-mRNAseq.gff$/ ||
               $value =~ /\/ecoli-contig.fasta$/ ||
               $value =~ /\/ecoli-est.fasta$/ ||
               $value =~ /\/ecoli-protein.fasta$/ ||
               $value =~ /\/ecoli.mod$/
               ){
		$is_tutorial++;
	    }

            if(! -e $value && $type !~ /^augustus_species$|^model_org$/){
                warn "WARNING: The menu option '$name => $value'\n".
		    "is not a valid file and will be ignored";
                next;
            }

            my $length = @{[stat($value)]}[7] || 0;
            my $c_count = 0;
            if($type =~ /^fasta$|^est$|^protein$|^repeat_protein$|^altest$|^genome$/){
                $length = MWAS_util::fasta_length($value);
                $c_count = `grep -c ">" $value`;
            }

            my ($check) = $dbh->selectrow_array(qq{SELECT value FROM menus WHERE type='$type' AND value='$value'});

            $dbh->do(qq{INSERT INTO menus (name, value, type, length, contig_count, is_tutorial) }.
                     qq{VALUES ('$name', '$value', '$type', '$length', '$c_count', $is_tutorial)}) unless($check);
        }
    }

    #log/re-log last unique key ids used (this is to help generate new keys)
    #Most versions of SQL can do this automatically but not all, so I will handle it in perl
    #please note that job_ids 1-5 aare reserved for tutorial data and user_id 1 for administrator
    #$dbh->do(qq{DROP TABLE id_store}) if(grep {/id_store/} @tables);
    $dsn = (qq{CREATE TABLE id_store (last_user_id INTEGER, last_job_id INTEGER, last_submit_id INTEGER)});
    if(! grep {/id_store/} @tables){
        $dbh->do($dsn);

        #fix value if null
        my ($last_user_id) = $dbh->selectrow_array(qq{SELECT user_id FROM users ORDER BY user_id DESC LIMIT 1}) || 1;
        my ($last_job_id)  = $dbh->selectrow_array(qq{SELECT job_id FROM jobs ORDER BY job_id DESC LIMIT 1}) || 1;
        my ($last_submit_id)  = $dbh->selectrow_array(qq{SELECT submit_id FROM jobs ORDER BY submit_id DESC LIMIT 1}) || 1;
        $last_job_id = 6 if($last_job_id < 6);
        $last_submit_id = $last_job_id if(! $last_submit_id);
        $dbh->do(qq{INSERT INTO id_store (last_user_id, last_job_id, last_submit_id) VALUES ($last_user_id, $last_job_id, $last_submit_id)});
    }
    add_missing_columns($dbh, $dsn, 'id_store');

    #add administrator to database
    my ($user_id) = $dbh->selectrow_array(qq{SELECT user_id FROM users WHERE user_id=1});
    $dbh->do(qq{INSERT INTO users (user_id, login, password, first, last, e_mail, is_guest, is_admin, date_created) }.
             qq{VALUES (1, 'admin', 'password', 'Server', 'Administrator', '', 0, 1, \'}.MWAS_util::date_time().qq{\')}) if(! $user_id);

    #update administrator e-mail
    if($dbh->selectrow_array(qq{SELECT e_mail FROM users WHERE e_mail='$O{admin_email}' and user_id != 1})){
        warn "WARNING:  The administrative e-mail address provided in the server.ctl file is already\n".
             "registered to another user.  The e-mail address will be set to empty and the control files\n".
             "will be updated accordingly\n";

        $O{admin_email} = '';
      GI::generate_control_files(config_loc(), 'server', \%O)
      }
    else{
        $dbh->do(qq{UPDATE users SET e_mail='$O{admin_email}' WHERE user_id=1}); #update admin e-mail address
    }

    #create/re-create all_default_opt table
    $dbh->do(qq{DROP TABLE all_default_opt}) if(grep {/all_default_opt/} @tables); #drop table
    @defaults = map {lc("\"$_\"")." $t_type"} grep {!/^menus$|^stat$/i} keys %O; #get table column names and datatype
    $dsn = 'CREATE TABLE all_default_opt ('.join(', ', @defaults).')';
    $dbh->do($dsn);

    @defaults = grep {!/^menus$|^stat$/i} keys %O; #keys to add
    my @lc_defaults = map {lc($_)}  @defaults;
    $dsn = "INSERT INTO all_default_opt (".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ('".join("', '", @O{@defaults})."')";
    $dbh->do($dsn); #add ctl_opt

    #create/re-create stat table
    my %stat = %{$O{STAT}};
    $dbh->do(qq{DROP TABLE stat}) if(grep {/stat/} @tables); #drop table
    @defaults = map {lc("\"$_\"")." $t_type"}  keys %stat; #get table column names and datatype
    $dsn = 'CREATE TABLE stat ('.join(', ', @defaults).')';
    $dbh->do($dsn);
    @defaults = keys %stat; #keys to add
    @lc_defaults = map {lc($_)}  @defaults;
    $dsn = "INSERT INTO stat (".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ('".join("', '", @stat{@defaults})."')";
    $dbh->do($dsn); #add ctl_opt

    #update values for possible null columns
    $dbh->do("UPDATE jobs SET type = 'maker' WHERE type IS NULL");
    $dbh->do("UPDATE jobs SET submit_id = job_id WHERE submit_id IS NULL");

    #--add tutorial jobs (always adds/re-adds)
    if($O{tutorials}){
        ### JOB 1 ###
        my $id = 1;
        my $t_message = "This example job will annotate the region of Drosophila melanogaster ".
            "chromosome 2L that encodes the gene decapentaplegic.  This is a gene ".
            "with multiple confirmed alternately spliced transcripts, and it ".
            "illustrates how EST evidence can suggest differtent splice forms ".
            "and UTR variation";

        #drop job from table
        $dbh->do(qq{DELETE FROM jobs WHERE job_id=$id}); #drop
        $dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$id}); #drop

        #add/re-add control file options
        my @files = @{$dbh->selectcol_arrayref(qq{SELECT value FROM menus WHERE is_tutorial=1 and }.
                                               qq{(value LIKE '%dpp_contig.fasta' or }.
                                                      qq{value LIKE '%dpp_est.fasta' or }.
                                                      qq{value LIKE '%dpp_protein.fasta')}
                                               )};
        if(@files == 3 && -f $files[0] && -f $files[1] && -f $files[2]){ #re-add here
            my %job_opt = %O; #make a copy of control file settings

            #add needed values
            ($job_opt{genome})  = grep {/dpp_contig.fasta/} @files;
            ($job_opt{est})     = grep {/dpp_est.fasta/} @files;
            ($job_opt{protein}) = grep {/dpp_protein.fasta/} @files;
            $job_opt{organism_type} = 'eukaryotic';
            $job_opt{est2genome} = 1;

            #add to ctl_opt table
            @defaults = (keys %def_opt); #keys to add
            my @lc_defaults = map {lc($_)}  @defaults;
            $dsn = "INSERT INTO ctl_opt (job_id, ".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ($id, '".join("', '", @job_opt{@defaults})."')";
            $dbh->do($dsn); #add ctl_opt
        }

        #add job but only if control file options were ok
        my ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$id});
        if($job_id){
            my $genome = $dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$id});
            ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM jobs WHERE job_id=$id});
            $dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, type, is_queued, is_started, }.
			    qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
			    qq{is_tutorial,  cpus, start_time, finish_time, name, note)}.
                     qq{VALUES($id, 1, \'}.MWAS_util::get_length_for_value($dbh, $genome).
			    qq{\', 'maker', 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, '}.MWAS_util::date_time.qq{', '}.
                     MWAS_util::date_time . qq{', 'D. melanogaster : Dpp','$t_message')}
                     ) if(! $job_id); #add job
        }

        ### JOB 2 ###
        $id = 2;
        $t_message = "This example job will annotate a region of the Pythium ultimum genome. ".
            "P. ultimum is an emerging model organism genome and this example illustrates ".
            "how to perform de novo annotaion on a new genome";

        #drop job from table
	$dbh->do(qq{DELETE FROM jobs WHERE job_id=$id}); #drop
        $dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$id}); #drop

        #add/re-add control file options
        @files = @{$dbh->selectcol_arrayref(qq{SELECT value FROM menus WHERE is_tutorial=1 and }.
                                            qq{(value LIKE '%pyu-contig.fasta' or }.
						   qq{value LIKE '%pyu-est.fasta' or }.
						   qq{value LIKE '%pyu-protein.fasta' or }.
						   qq{value LIKE '%pyu.hmm' or }.
						   qq{value LIKE '%pyu.mod')}
                                            )};
        if(@files >= 3 && -f $files[0] && -f $files[1] && -f $files[2]){ #re-add here
            my %job_opt = %O; #make a copy of control file settings

            #add needed values
            ($job_opt{genome})  = grep {/pyu-contig.fasta/} @files;
            ($job_opt{est})     = grep {/pyu-est.fasta/} @files;
            ($job_opt{protein}) = grep {/pyu-protein.fasta/} @files;
            ($job_opt{snaphmm}) = grep {/pyu.hmm/} @files;
            ($job_opt{gmhmm}) = grep {/pyu.mod/} @files;
            $job_opt{organism_type} = 'eukaryotic';
            $job_opt{est2genome} = 1;

            #add to ctl_opt table
            @defaults = (keys %def_opt); #keys to add
            my @lc_defaults = map {lc($_)}  @defaults;
            $dsn = "INSERT INTO ctl_opt (job_id, ".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ($id, '".join("', '", @job_opt{@defaults})."')";
            $dbh->do($dsn); #add ctl_opt
        }

        #add job but only if control file options were ok
        ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$id});
        if($job_id){
            my $genome = $dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$id});
            ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM jobs WHERE job_id=$id});
            $dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, type, is_queued, is_started, }.
			    qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
			    qq{is_tutorial,  cpus, start_time, finish_time, name, note)}.
                     qq{VALUES($id, 1, \'}.MWAS_util::get_length_for_value($dbh, $genome).
			    qq{\', 'maker', 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, '}.MWAS_util::date_time.qq{', '}.
                     MWAS_util::date_time . qq{', 'P. ultimum : De Novo Annotation','$t_message')}
                     ) if(! $job_id); #add job
        }

        ### JOB 3 ###
	$id = 3;
        $t_message = "This example job will annotate a region of the Pythium ultimum genome. ".
            "In the example mRNAseq data crossing spice site junctions is integrated ".
            "into a MAKER run. The data data was processed and aligned by an external ".
            "program and will be passed to MAKER in a GFF3 formated file. This example ".
            "demonstrates how you might use the GFF3 pass-through capability of MAKER.";

        #drop job from table
        $dbh->do(qq{DELETE FROM jobs WHERE job_id=$id}); #drop
        $dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$id}); #drop

        #add/re-add control file options
        @files = @{$dbh->selectcol_arrayref(qq{SELECT value FROM menus WHERE is_tutorial=1 and }.
                                            qq{(value LIKE '%pass-contig.fasta' or }.
						   qq{value LIKE '%pyu-est.fasta' or }.
						   qq{value LIKE '%pyu-protein.fasta' or }.
						   qq{value LIKE '%pass-mRNAseq.gff' or }.
						   qq{value LIKE '%pyu.hmm' or }.
						   qq{value LIKE '%pyu.mod')}
                                            )};
        if(@files >= 4 && -f $files[0] && -f $files[1] && -f $files[2] && -f $files[3]){ #re-add here
            my %job_opt = %O; #make a copy of control file settings

            #add needed values
            ($job_opt{genome})  = grep {/pass-contig.fasta/} @files;
            ($job_opt{est})     = grep {/pyu-est.fasta/} @files;
            ($job_opt{protein}) = grep {/pyu-protein.fasta/} @files;
            ($job_opt{est_gff}) = grep {/pass-mRNAseq.gff/} @files;
            ($job_opt{snaphmm}) = grep {/pyu.hmm/} @files;
            ($job_opt{gmhmm}) = grep {/pyu.mod/} @files;
            $job_opt{organism_type} = 'eukaryotic';
            $job_opt{est2genome} = 1;

            #add to ctl_opt table
            @defaults = (keys %def_opt); #keys to add
            my @lc_defaults = map {lc($_)}  @defaults;
            $dsn = "INSERT INTO ctl_opt (job_id, ".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ($id, '".join("', '", @job_opt{@defaults})."')";
            $dbh->do($dsn); #add ctl_opt
        }

        #add job but only if control file options were ok
        ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$id});
        if($job_id){
            my $genome = $dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$id});
            ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM jobs WHERE job_id=$id});
            $dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, type, is_queued, is_started, }.
			    qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
			    qq{is_tutorial,  cpus, start_time, finish_time, name, note)}.
                     qq{VALUES($id, 1, \'}.MWAS_util::get_length_for_value($dbh, $genome).
			    qq{\', 'maker', 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, '}.MWAS_util::date_time.qq{', '}.
                     MWAS_util::date_time . qq{', 'P. ultimum : Data Pass-through','$t_message')}
                     ) if(! $job_id); #add job
        }

        ### JOB 4 ###
        $id = 4;
        $t_message = "This example job will annotate a region of the Pythium ultimum genome. ".
            "In this example two pre-existing legacy annotation sets will be merged ".
            "and integrated with EST and protein homology data to produce a consensus ".
            "set of annotations.  This example demonstrates how legacy annotations can be ".
            "integrated into a MAKER run. Legacy annotation sets are passed to MAKER in ".
            "GFF3 format.";

        #drop job from table
        $dbh->do(qq{DELETE FROM jobs WHERE job_id=$id}); #drop
        $dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$id}); #drop

        #add/re-add control file options
        @files = @{$dbh->selectcol_arrayref(qq{SELECT value FROM menus WHERE is_tutorial=1 and }.
                                            qq{(value LIKE '%legacy-contig.fasta' or }.
						   qq{value LIKE '%pyu-est.fasta' or }.
						   qq{value LIKE '%pyu-protein.fasta' or }.
						   qq{value LIKE '%legacy-set1.gff' or }.
						   qq{value LIKE '%legacy-set2.gff' or }.
						   qq{value LIKE '%pyu.hmm' or }.
						   qq{value LIKE '%pyu.mod')}
                                            )};
        if(@files >= 5 && -f $files[0] && -f $files[1] && -f $files[2] && -f $files[3] && -f $files[4]){ #re-add here
            my %job_opt = %O; #make a copy of control file settings

            #add needed values
            ($job_opt{genome})  = grep {/legacy-contig.fasta/} @files;
            ($job_opt{est})     = grep {/pyu-est.fasta/} @files;
            ($job_opt{protein}) = grep {/pyu-protein.fasta/} @files;
            ($job_opt{model_gff}) = grep {/legacy-set1.gff/} @files;
            ($job_opt{model_gff}) .= ",".grep {/legacy-set2.gff/} @files;
            ($job_opt{snaphmm}) = grep {/pyu.hmm/} @files;
            ($job_opt{gmhmm}) = grep {/pyu.mod/} @files;
            $job_opt{organism_type} = 'eukaryotic';
            $job_opt{est2genome} = 1;

            #add to ctl_opt table
            @defaults = (keys %def_opt); #keys to add
            my @lc_defaults = map {lc($_)}  @defaults;
            $dsn = "INSERT INTO ctl_opt (job_id, ".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ($id, '".join("', '", @job_opt{@defaults})."')";
            $dbh->do($dsn); #add ctl_opt
	}

        #add job but only if control file options were ok
        ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$id});
        if($job_id){
            my $genome = $dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$id});
            ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM jobs WHERE job_id=$id});
            $dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, type, is_queued, is_started, }.
			    qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
			    qq{is_tutorial,  cpus, start_time, finish_time, name, note)}.
                     qq{VALUES($id, 1, \'}.MWAS_util::get_length_for_value($dbh, $genome).
			    qq{\', 'maker', 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, '}.MWAS_util::date_time.qq{', '}.
                     MWAS_util::date_time . qq{', 'P. ultimum : Legacy Annotations','$t_message')}
                     ) if(! $job_id); #add job
        }

        ### JOB 5 ###
        $id = 5;
        $t_message = "This example job will annotate a region of the Escherichia coli genome. ".
            "This is an example of how to use MAKER for annotating prokaryotic genomes";

        #drop job from table
        $dbh->do(qq{DELETE FROM jobs WHERE job_id=$id}); #drop
        $dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$id}); #drop

        #add/re-add control file options
        @files = @{$dbh->selectcol_arrayref(qq{SELECT value FROM menus WHERE is_tutorial=1 and }.
                                            qq{(value LIKE '%ecoli-contig.fasta' or }.
						   qq{value LIKE '%ecoli-est.fasta' or }.
						   qq{value LIKE '%ecoli-protein.fasta' or }.
						   qq{value LIKE '%ecoli.mod')}
                                            )};
        if(@files >= 4 && -f $files[0] && -f $files[1] && -f $files[2] && -f $files[3]){ #re-add here
            my %job_opt = %O; #make a copy of control file settings

            #add needed values
            ($job_opt{genome})  = grep {/ecoli-contig.fasta/} @files;
            ($job_opt{est})     = grep {/ecoli-est.fasta/} @files;
            ($job_opt{protein}) = grep {/ecoli-protein.fasta/} @files;
            ($job_opt{gmhmm}) = grep {/ecoli.mod/} @files;
            $job_opt{organism_type} = 'prokaryotic';
	    $job_opt{protein2genome} = 1;

            #add to ctl_opt table
            @defaults = (keys %def_opt); #keys to add
            my @lc_defaults = map {lc($_)}  @defaults;
            $dsn = "INSERT INTO ctl_opt (job_id, ".join(", ", map {"\"$_\""} @lc_defaults).") VALUES ($id, '".join("', '", @job_opt{@defaults})."')";
            $dbh->do($dsn); #add ctl_opt
        }

        #add job but only if control file options were ok
        ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$id});
        if($job_id){
            my $genome = $dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$id});
            ($job_id) = $dbh->selectrow_array(qq{SELECT job_id FROM jobs WHERE job_id=$id});
            $dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, type, is_queued, is_started, }.
			    qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
			    qq{is_tutorial,  cpus, start_time, finish_time, name, note)}.
                     qq{VALUES($id, 1, \'}.MWAS_util::get_length_for_value($dbh, $genome).
			    qq{\', 'maker', 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, '}.MWAS_util::date_time.qq{', '}.
                     MWAS_util::date_time . qq{', 'E. coli : De Novo Annotation','$t_message')}
                     ) if(! $job_id); #add job
        }
    }
    ### END ADD TUTORIAL JOBS ###

    #commit changes
    #$dbh->commit;
    #$lock->unlock;
    ####===unlock
}
#------------------------------------------------------------------------
sub add_missing_columns {
    my $dbh = shift;
    my $dsn = shift; # must be the statement used to build table
    my $table = shift;

    #check that all columns are present
    my $bak = $dsn;
    ($dsn) = $dsn =~ /\((.*)\)/;
    my %columns = map {/^\"?([^\s\"]+)\"?\s+(.*)/} split(/,\s*/, $dsn);

    my $sth = $dbh->prepare("SELECT * FROM $table LIMIT 1");
    $sth->execute;
    my @missing = ();
    my @found = map {lc($_)} @{$sth->{NAME}}; #make all lower case for any DBI
    my %found;
    @found{@found} = (); #build indexing hash

    #no columns were found to add, so delete table and rebuild
    if(! @found){
        $dbh->do("DROP TABLE $table");
        $dbh->do($bak); #create new table
    }

    foreach my $key (keys %columns){
        push(@missing, $key) if(! exists $found{$key});
    }

    #add missing columns
    foreach my $key (@missing){
        my $def = $columns{$key};
        $dbh->do("ALTER TABLE $table ADD $key $def");
    }
}
#------------------------------------------------------------------------
#puts jobs back into queue on server restart
sub reset_queue {
    my $dbh = shift;
    my $dir = shift;

    dbh_do_commit($dbh, qq{UPDATE jobs SET is_running=0, is_queued=1, is_started=1, is_finished=0, cpus=0 }.
                  qq{WHERE (is_started=1 OR is_running=1) and is_packaged=0 and is_error=0}, $dir);

    dbh_do_commit($dbh, qq{UPDATE jobs SET is_running=0, is_started=1, cpus=0 }.
                  qq{WHERE is_running=1}, $dir);
}
1;
