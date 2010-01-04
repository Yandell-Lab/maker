package MWS;

BEGIN{
   if (not ($ENV{CGL_SO_SOURCE})) {
      $ENV{CGL_SO_SOURCE} = "$FindBin::Bin/lib/CGL/so.obo";
   }
   if (not ($ENV{CGL_GO_SOURCE})) {
      $ENV{CGL_GO_SOURCE} = "$FindBin::Bin/lib/CGL/gene_ontology.obo"
   }
   #{ $ENV{'CAP_DEVPOPUP_EXEC'} = 1; }
}

use strict;
use warnings;

use base qw(CGI::Application);
use CGI::Application::Plugin::Session;
use CGI::Application::Plugin::Authentication;
use CGI::Application::Plugin::CAPTCHA;
use CGI::Application::Plugin::DBH (qw/dbh_config dbh/);
use CGI::Application::Plugin::TT;
use CGI::Application::Plugin::DevPopup;
use CGI::Application::Plugin::DevPopup::HTTPHeaders;
use CGI::Application::Plugin::Redirect;

use FindBin;
use GI;
use Data::Dumper;
use MWAS_util;
use URI::Escape;

#-----------------------------------------------------------------------------
sub cgiapp_init {
   my $self = shift;
   $self->SUPER::cgiapp_init;

   #load the server control files
   my %serv_opt = GI::set_defaults('server', {GI::parse_ctl_files(["$FindBin::Bin/config/server.ctl"])});

   #make sure required database values are setup
   if (! $serv_opt{DBI}) {
      die "ERROR: You must specify a DBI connection method in: $FindBin::Bin/config/server.ctl\n\n";
   }
   if (! $serv_opt{dbname}) {
      die "ERROR: You must specify a database to connect to in: $FindBin::Bin/config/server.ctl\n\n";
   }
    
   #connect to the database
   my $dsn = "DBI:$serv_opt{DBI}:dbname=$serv_opt{dbname};";
   $dsn .= "host=$serv_opt{host};" if($serv_opt{host});
   $dsn .= "port=$serv_opt{port};" if($serv_opt{host} && $serv_opt{port});

   $self->dbh_config($dsn, $serv_opt{username}, $serv_opt{password}, {AutoCommit => 0}) 
     or die "Got error $DBI::errstr when connecting to database\n";;

   #setup template params
   $self->tt_config(TEMPLATE_OPTIONS => {INCLUDE_PATH => "$FindBin::Bin/tt_templates/",
				         EVAL_PERL => 1});

   #setup authentication
   __PACKAGE__->authen->config(DRIVER => ['DBI',
					  DBH         => $self->dbh,
					  TABLE       => 'users',
					  CONSTRAINTS => {'users.login'    => '__CREDENTIAL_1__',
							  'users.password' => '__CREDENTIAL_2__',
							 }
					 ],
			       STORE => 'Session',
			       LOGIN_RUNMODE => 'home_login',
			       POST_LOGIN_RUNMODE => 'frontpage',
			       LOGOUT_RUNMODE => 'home_login',
			       LOGIN_SESSION_TIMEOUT => {IDLE_FOR => '30m',
							 EVERY => '1d'
							},
			      );
   
   $self->authen->protected_runmodes(qw(edit_account
					feedback
					filebox
					upload_file
					launch
					results
					submit_to_db
					delete_job
					job_create)
				     );

   #add default control options from server
   $self->param(server_opt => \%serv_opt);
}
#-----------------------------------------------------------------------------
sub cgiapp_prerun {
        my $self = shift;

        $self->tt_params({logged_in  => $self->authen->is_authenticated,
                          server_opt => $self->param('server_opt'),
                          session    => $self->session,
                          #query    => Dumper($self->query),
		      });
}
#-----------------------------------------------------------------------------
sub setup {
   my $self = shift;

   $self->start_mode('home_login');
   $self->run_modes([qw(login
			home_login
			guest_login
			authen_login
			forgot_login
			register
			launch
			help
			results
			frontpage
			upload_file
			job_create
			create_captcha
			filebox
			edit_account
			feedback
			submit_to_db
			queue
			delete_job
			job_status)]
		   );

   $self->captcha_config(IMAGE_OPTIONS => {width    => 150,
					   height   => 40,
					   lines    => 5,
					   font  => '/usr/share/fonts/bitstream-vera/VeraMono.ttf',
					   ptsize   => 18,
					   bgcolor  => "#FFFFFF"
					  },
			 CREATE_OPTIONS   => [ 'ttf', 'rect' ],
			 PARTICLE_OPTIONS => [ 300 ]
			);
}
#-----------------------------------------------------------------------------
sub teardown {
   my $self = shift;

   $self->dbh->disconnect if($self->dbh);
}
#-----------------------------------------------------------------------------
sub login {
   my $self = shift;
   my $add_text= shift;

   if($self->authen->is_authenticated){
       my $user_id = $self->get_user_id;
       return $self->redirect("maker.cgi?rm=frontpage&user_id=$user_id");
   }
   else{
       return $self->tt_process('maker_login.tt', {message => $add_text});
   }
}
#-----------------------------------------------------------------------------
sub delete_job {
   my $self = shift;

   my $q = $self->query;
   my $job_id = $q->param('job_id');
   my $user = $self->get_user_info();
   my $job = $self->get_job_info($job_id);

   if($job){
       die "ERROR: This is not your job to delete" if($user->{user_id} != $job->{user_id});
       
       $self->dbh->do(qq{DELETE FROM jobs WHERE job_id=$job_id});
       $self->dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$job_id});
       $self->dbh->commit();
       my %serv_opt = $self->param('server_opt');
       File::Path::rmtree("$serv_opt{data_dir}/job/$job_id");
   }

   my $user_id = $self->get_user_id();
   return $self->redirect("maker.cgi?rm=frontpage&user_id=$user_id");
}
#-----------------------------------------------------------------------------
sub home_login {
   my $self = shift;
    
   return $self->login("");
}
#-----------------------------------------------------------------------------  
sub authen_login {
   my $self = shift;

   return $self->login("Invalid username or password! Try again");
}
#-----------------------------------------------------------------------------
sub forgot_login {
   my $self = shift;

   my $q = $self->query();
   my $e_mail = $q->param('e_mail');
   my $serv_opt = $self->param('server_opt');

   my $user = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE e_mail='$e_mail'}) if($e_mail);

   my $mm;
   if(defined $user && defined $user->{login}){
       $mm = "MWAS: Forgotten login/password request\n\n".
	     "Username: ".$user->{login}."\n".
	     "Password: ".$user->{password}."\n";

       MWAS_util::send_message($user->{e_mail},
				  $serv_opt->{smtp},
				  'MWAS: Forgotten login/password request',
				  $mm
				 );

       return $self->login("Your username and password have been e-mailed to you");    
   }
   elsif(defined $e_mail){
       return $self->login("No account associated with e-mail address");    
   }    
   else{
       return $self->tt_process('forgot_login.tt');
   }
}
#-----------------------------------------------------------------------------
sub launch {
   my $self = shift;

   my $q = $self->query();
   my $job_id = $q->param('job_id');
   my $value = $q->param('contig');
   my $user_id = $self->get_user_id();
   $value =~ s/\/+$//;
   my ($name) = $value =~ /([^\/]+)$/;
   my $data_dir = $self->param('server_opt')->{data_dir};

   #make path
   File::Path::mkpath("/data/var/www/html/MWAS/users/$user_id/");
   my $gff = "http://malachite.genetics.utah.edu/MWAS/users/$user_id/$name.gff";
   my $xml = "http://malachite.genetics.utah.edu/MWAS/users/$user_id/$name.xml";

   open(OUT,  "> /data/var/www/html/MWAS/users/$user_id/$name.gff");
   open(IN, "< $data_dir/jobs/$job_id/$job_id.maker.output/$value/$name.gff");
   while(my $line = <IN>){
       print OUT $line;
   }
   close(OUT);
   close(IN);

   open(OUT,  "> /data/var/www/html/MWAS/users/$user_id/$name.xml");
   open(IN, "< $data_dir/jobs/$job_id/$job_id.maker.output/$value/$name.xml");
   while(my $line = <IN>){
       print OUT $line;
   }
   close(OUT);
   close(IN);

   if($q->param('apollo')){
       my $url_base_dir = "http://malachite.genetics.utah.edu/MWAS/";
       my $url_jnlp_base = "http://malachite.genetics.utah.edu/MWAS/users/$user_id/";
       my $url_gff3_file = $xml;;
       
       open(OUT,  "> /data/var/www/html/MWAS/users/$user_id/apollo.jnlp");
       open(IN, "< tt_templates/apollo.tt");
       while(my $line = <IN>){
	   $line =~ s/\[\% url_base_dir \%\]/$url_base_dir/g;
	   $line =~ s/\[\% url_jnlp_base \%\]/$url_jnlp_base/g;
	   $line =~ s/\[\% url_gff3_file \%\]/$url_gff3_file/g;
	   
	   print OUT $line;
       }
       close(OUT);
       close(IN);
       
       return $self->redirect("http://malachite.genetics.utah.edu/MWAS/users/$user_id/apollo.jnlp");
   }
   elsif($q->param('soba')){
      
       return $self->redirect("http://www.sequenceontology.org/cgi-bin/soba.cgi?rm=upload_urls&url=$gff");
   }
}
#-----------------------------------------------------------------------------
sub results {
   my $self = shift;
   my $serv_opt = $self->param('server_opt');
   my $q = $self->query();
   my $job_id = $q->param('job_id');
   my $user_id = $q->param('user_id');

   #raw output dir
   my $dir = $serv_opt->{data_dir}."/jobs/$job_id/$job_id.maker.output/";
   my $index = "$dir/$job_id\_master_datastore_index.log";

   die "ERROR: Result can not be found in $dir\n" if(! -d $dir);
   die "ERROR: Index file $index can not be found\n" if(! -f $index);

   my %values;
   my %status;
   open(my $IN, "< $index");
   while(my $line = <$IN>){
       chomp $line;
       my @data = split("\t", $line);

       $values{$data[0]} = $data[1] if($data[2] eq 'FINISHED');
       $status{$data[0]} = $data[2] unless($status{$data[0]} eq 'FINISHED');
   }
   close($IN);
   
   my @menus = map {{name => $_, value => $values{$_}}} keys %values;

   my %counts;
   @counts{qw(INCOMPLETE FINISHED FAILED SKIPPED TOTAL)} = map{0} qw(1 2 3 4 5); 
   foreach my $val (values %status){
       if($val eq 'RETRY' || $val eq 'STARTED'){
	   $counts{INCOMPLETE}++;
       }
       elsif($val eq 'FINISHED'){
	   $counts{FINISHED}++;
       }
       elsif($val =~ /DIED/){
	   $counts{FAILED}++;
       }
       elsif($val =~ /SKIPPED/){
	   $counts{SKIPPED}++;
       }
   }

   my ($genome) = $self->dbh->selectrow_array(qq{SELECT genome FROM ctl_opt WHERE job_id=$job_id});
   my ($c_count) = $self->dbh->selectrow_array(qq{SELECT contig_count FROM files WHERE value='$genome'});
   ($c_count) = $self->dbh->selectrow_array(qq{SELECT contig_count FROM menus WHERE value='$genome'}) if(! $c_count);
   
   $counts{TOTAL} = $c_count if($c_count);

   return $self->tt_process('results.tt', {menus => \@menus,
					   counts => \%counts,
					   user_id => $user_id,
					   job_id => $job_id});   
}
#-----------------------------------------------------------------------------
sub register {
   my $self = shift;

   my $serv_opt = $self->param('server_opt');
   my $q = $self->query();
   my $login = $q->param('login');
   my $first = $q->param('first');
   my $last = $q->param('last');
   my $institution = $q->param('institution');
   my $e_mail = $q->param('e_mail');
   my $password = $q->param('password');
   my $verify = $q->param('verify');
   my $captcha = $q->param("captcha");

   #no input this must be a new request
   if(!$login && !$first && !$last && !$institution &&
      !$e_mail && !$password && !$verify && !$captcha){
       return $self->tt_process('register.tt');
   }
   
   #holds any errors related to bad input
   my %errors;
   #captcha errors
   if(! $captcha){
       $errors{captcha} = 'Required';
   }
   elsif(!$self->captcha_verify($q->cookie("hash"), $captcha)){
       $errors{captcha} = 'Value did not match image';
   }
   #login errors
   if(! $login){
       $errors{login} = 'Required';
   }
   elsif($self->dbh->selectrow_array(qq{SELECT user_id FROM users WHERE login='$login'})){
       $errors{login} = 'ID already registered to another user';
   }
   #e_mail errors
   if(! $e_mail){
       $errors{e_mail} = 'Required';
   }
   elsif($self->dbh->selectrow_array(qq{SELECT e_mail FROM users WHERE e_mail='$e_mail'})){
       $errors{e_mail} = 'Address already used by another account';
   }
   #password errors
   if(! $password){
       $errors{password} = 'Required';
   }
   #password verify errors
   if($password && (!$verify || $password ne $verify)){
       $errors{verify} = 'Did not match password';
   }

   #depending on errors go back or create account
   if(keys %errors){ #there are errors
       return $self->tt_process('register.tt', {errors => \%errors,
						login => $login,
						first => $first,
						last  => $last,
						institution => $institution,
						e_mail => $e_mail,
						password => $password,
						verify => $verify
						}); #was general_add.tt
   }
   else{ #everything is fine, add user to database
#       my $lock = MWAS_util::lockDB($serv_opt->{data_dir});
       my ($new_id) = $self->dbh->selectrow_array(qq{SELECT last_user_id FROM id_store}); #get last new user_id
       $new_id++; #iterate the value
       $self->dbh->do(qq{UPDATE id_store SET last_user_id=$new_id}); #record new value
       $self->dbh->do(qq{INSERT INTO users (user_id, login, password, first, last, e_mail, }.
		      qq{institution, is_guest, is_admin, date_created, last_login) VALUES ($new_id, }.
		      qq{'$login', '$password', '$first', '$last', '$e_mail', '$institution', 0, 0, \'}.
		      MWAS_util::date_time().qq{\', '')});
       $self->dbh->commit();
#       $lock->unlock;

       #authenticate user
       $q->param(authen_username => $login);
       $q->param(authen_password => $password);
       $self->authen->{initialized} = 0; #force reauthentication
       $self->authen->initialize();

       my $user_id = $self->get_user_id();
       #forward user to account frontpage
       return $self->redirect("maker.cgi?rm=frontpage&user_id=$user_id");
   }
}
#-----------------------------------------------------------------------------
sub edit_account {
   my $self = shift;

   my $user = $self->get_user_info();
   my $serv_opt = $self->param('server_opt');
   my $q = $self->query();
   my $login = $q->param('login');
   my $first = $q->param('first');
   my $last = $q->param('last');
   my $institution = $q->param('institution');
   my $e_mail = $q->param('e_mail');
   my $old_password = $q->param('old_password');
   my $new_password = $q->param('new_password');
   my $verify = $q->param('verify');

   #no input this must be a new request
   if(!$login && !$first && !$last && !$institution &&
      !$e_mail && !$old_password && ! $new_password && !$verify
     ){
       return $self->tt_process('edit_account.tt', {user => $user});
   }
   
   #holds any errors related to bad input
   my %errors;
   #login errors
   if(! $login){
       $errors{login} = 'Empty value not permitted';
   }
   elsif($self->dbh->selectrow_array(qq{SELECT user_id FROM users WHERE login='$login' }.
				     qq{and user_id != }.$user->{user_id})
	){

       $errors{login} = 'Could not use new ID, already exists';
   }
   #e_mail errors
   if(! $e_mail){
       $errors{e_mail} = 'Empty value not permited';
   }
   elsif($self->dbh->selectrow_array(qq{SELECT e_mail FROM users WHERE e_mail='$e_mail'}.
				     qq{and user_id != }.$user->{user_id})
	){

       $errors{e_mail} = 'Could not use new address, already exists';
   }
   #old password errors
   if(!$user->{is_guest} && !$old_password){
       $errors{old_password} = 'You must enter your password to update the account';
   }
   elsif(!$user->{is_guest} && $old_password ne $user->{password}){
       $errors{old_password} = 'The password you entered was incorrect';
   }
   #new password errors
   if(! $new_password && $user->{is_guest}){
       $errors{old_password} = 'Empty value not permitted';
   }
   #password verify errors
   if($new_password && (!$verify || $new_password ne $verify)){
       $errors{verify} = 'Did not match password';
   }

   #depending on errors go back or create account
   if(keys %errors){ #there are errors
       return $self->tt_process('edit_account.tt', {errors => \%errors,
						    user => $user
						    }); #was general_add.tt
   }
   else{ #everything is fine, update user in database
#      my $lock = MWAS_util::lockDB($serv_opt->{data_dir});
       my $password = ($new_password) ? $new_password : $old_password; #select password to use
       $self->dbh->do(qq{UPDATE users SET login='$login', first='$first', last='$last', e_mail='$e_mail', }.
		      qq{password='$password', institution='$institution', is_guest=0 WHERE user_id=}.
		      $user->{user_id});
       $self->dbh->commit();
#       $lock->unlock;

       #update user info
       $user = $self->get_user_info();
       return $self->tt_process('edit_account.tt', {user => $user,
						    message => 'Account updated successfullyy'
						    }); #was general_add.tt
   }
}
#-----------------------------------------------------------------------------
sub guest_login {
    my $self = shift;

#   my $lock = MWAS_util::lockDB($serv_opt->{data_dir});
    my ($new_id) = $self->dbh->selectrow_array(qq{SELECT last_user_id FROM id_store}); #get last new user_id
    $new_id++; #iterate the value
    $self->dbh->do(qq{UPDATE id_store SET last_user_id=$new_id}); #record new value

    my $login = "guest_$new_id";
    my $password = '';

    $self->dbh->do(qq{INSERT INTO users (user_id, login, password, first, last, e_mail, }.
		   qq{institution, is_guest, is_admin, date_created, last_login) VALUES ($new_id, }.
		   qq{'$login', '$password', '', '', '', '', 1, 0, \'}.
		   MWAS_util::date_time().qq{\', '')});
    $self->dbh->commit();
#   $lock->unlock;

    #authenticate user
    $self->query->param(authen_username => $login);
    $self->query->param(authen_password => $password);
    $self->authen->{initialized} = 0; #force reauthentication
    $self->authen->initialize();
   my $user_id = $self->get_user_id();    
    #forward user to account frontpage
    return $self->redirect("maker.cgi?rm=frontpage&user_id=$user_id");
}
#-----------------------------------------------------------------------------
sub help {
   my $self = shift;
    
   return $self->tt_process('help.tt'); #was general_help.tt
}
#-----------------------------------------------------------------------------
sub frontpage {
   my $self = shift;
   my $message = shift;

   if(! $self->authen->is_authenticated){
       my $in_id = $self->query->param('user_id');
       my $q = $self->query;

       if($in_id){
	   my $user = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE user_id=$in_id});
	   return $self->home_login if(! $user || !$user->{is_guest});

	   #authenticate user
	   $q->param(authen_username => $user->{login});
	   $q->param(authen_password => $user->{password});
	   $self->authen->{initialized} = 0; #force reauthentication
	   $self->authen->initialize();
       }
       else{
	   return $self->home_login;
       }
   }

   my $user = $self->get_user_info();

   my $jobs = $self->dbh->selectall_arrayref(qq{SELECT * FROM jobs WHERE user_id=}.
					     $user->{user_id}.qq{ AND is_saved=1 ORDER BY job_id DESC},
					     {Slice => {}}
					    );

   #set job status
   foreach my $job (@$jobs){
       if($job->{is_packaged}){
	   $job->{status} = 'results ready';
       }
       elsif($job->{is_finished} && !$job->{is_running}){
	   $job->{status} = 'finishing';
       }
       elsif($job->{admin_block}){
	   $job->{status} = 'blocked';
       }
       elsif($job->{is_running}){
	   $job->{status} = 'running';
       }
       elsif($job->{is_started} && !$job->{is_queued}){
	   $job->{status} = 'idle';
       }
       elsif($job->{is_queued}){
	   $job->{status} = 'waiting in queue';
       }
       elsif($job->{is_saved}){
	   $job->{status} = 'edit';
       }
   }

   return $self->tt_process('frontpage.tt', {jobs => $jobs,
					     user => $user,
					     message => $message});
}
#-----------------------------------------------------------------------------
sub job_create {
   my $self = shift;
   my $message = shift;

   my $q = $self->query();
   my $job_id = $q->param('job_id');
   my $user = $self->get_user_info();
   my $job = $self->get_job_info($job_id);

   #load default values for all control files
   my %CTL_OPT = %{$self->get_server_default_options()};

   #load default values for stat
   my %STAT = %{$self->get_stat()};

   #get any logged control file options
   my %LOG_OPT = %CTL_OPT;
   if($job_id){
       my $ref = $self->dbh->selectrow_hashref(qq{SELECT * FROM ctl_opt WHERE job_id=$job_id});
       %LOG_OPT = (%CTL_OPT, %{$ref});
   }

   #build predictor hash
   my @predictors = split(",", $LOG_OPT{predictor});
   my %preds;
   @preds{@predictors} = map {1} @predictors;
   $LOG_OPT{predictor_hash} = \%preds;
   
   #reserve new job id value if needed (tutorials are started with new job_id)
   if(! $job_id || $job->{is_tutorial} || $job->{is_started}){
       ($job_id) = $self->dbh->selectrow_array(qq{SELECT last_job_id FROM id_store}); #get last job_id
       $job_id++; #iterate the value
       $self->dbh->do(qq{UPDATE id_store SET last_job_id=$job_id}); #record new value
       $self->dbh->commit();
   }

   #decide on organism type
   my $o_type = ($STAT{organism_type} !~ /STATIC|DISABLED/) ? $LOG_OPT{organism_type} : $CTL_OPT{organism_type};

   #collect options for drop down menus
   my %menus;
   
   #server menu options (menu options follow control files)
   $menus{alt_peptide}         = [qw(A C D E F G H I K L M N P Q R S T V W Y)];
   $menus{model_org}{server}   = $self->get_menus('model_org', 'server');
   $menus{snaphmm}{server}     = $self->get_menus('snaphmm', 'server');
   $menus{gmhmm}{server}     = $self->get_menus('gmhmm_E', 'server') if($o_type eq 'eukaryotic');
   $menus{gmhmm}{server}     = $self->get_menus('gmhmm_P', 'server') if($o_type eq 'prokaryotic');
   $menus{augustus_species}{server} = $self->get_menus('augustus_species', 'server');
   $menus{fgenesh_par_file}{server} = $self->get_menus('fgenesh_par_file', 'server');
   $menus{genome}{server}      = $self->get_menus('genome', 'server');
   $menus{est}{server}         = $self->get_menus('est', 'server');
   $menus{altest}{server}      = $self->get_menus('altest', 'server');
   $menus{protein}{server}     = $self->get_menus('protein', 'server');
   $menus{repeat_protein}{server} = $self->get_menus('repeat_protein', 'server');
   $menus{rmlib}{server}       = $self->get_menus('rmlib', 'server');
   $menus{model_gff}{server}   = $self->get_menus('model_gff', 'server');
   $menus{pred_gff}{server}    = $self->get_menus('pred_gff', 'server');
   $menus{est_gff}{server}     = $self->get_menus('est_gff', 'server');
   $menus{altest_gff}{server}  = $self->get_menus('altest_gff', 'server');
   $menus{protein_gff}{server} = $self->get_menus('protein_gff', 'server');
   $menus{rm_gff}{server}  = $self->get_menus('rm_gff', 'server');
   
   #user supplied menu options
   $menus{snaphmm}{user}  = $self->get_menus('snaphmm', 'user', $user->{user_id});
   $menus{augustus_species}{user}  = $self->get_menus('augustus_species', 'user', $user->{user_id});
   $menus{fgenesh_par_file}{user} = $self->get_menus('fgenesh_par_file', 'user', $user->{user_id});
   $menus{gmhmm}{user}  = $self->get_menus('gmhmm_E', 'user', $user->{user_id}) if($o_type eq 'eukaryotic');
   $menus{gmhmm}{user}  = $self->get_menus('gmhmm_P', 'user', $user->{user_id}) if($o_type eq 'prokaryotic');
   $menus{fastas}{user}   = $self->get_menus('fastas', 'user', $user->{user_id});
   $menus{gff3}{user}     = $self->get_menus('gff3', 'user', $user->{user_id});

   #example/tutorial menu options
   $menus{tutorials} = $self->get_tutorials();

   return $self->tt_process('job_create.tt',{menus   => \%menus,
					     ctl_opt => \%CTL_OPT, #default values
					     log_opt => \%LOG_OPT, #logged values
					     stat    => \%STAT, #opt status values
					     o_type  => $o_type,
					     user    => $user,
					     message    => $message,
					     job_id  => $job_id});
}
#------------------------------------------------------------------------------
sub submit_to_db {
   my $self = shift;

   my $q = $self->query;
   my $user_id = $self->get_user_id;
   my $is_queued = $q->param('add') ? 1 : 0;
   my $later = $q->param('later') ? 1 : 0; #save and come back later
   my $job_id = $q->param('job_id');

   #decide if this is a temporary ctl_opt cache update or not
   my ($is_saved) = $self->dbh->selectrow_array(qq{SELECT is_saved FROM jobs WHERE job_id=$job_id});
   if($q->param('just_update_opt') && ! $is_saved){
       $is_saved = 0; #only disable permanet save if user is just updating control file cache
   }
   else{
       $is_saved = 1;
   }

   #build new control options
   my %CTL_OPT = %{$self->get_server_default_options()};
   %CTL_OPT = (GI::set_defaults('opts', \%CTL_OPT),
	       GI::set_defaults('bopts', \%CTL_OPT)); #filter to needed sub-set
   while(my $key = each %CTL_OPT){
       my $value = $q->param($key);
       $CTL_OPT{$key} = $value if(defined $value);
   }

   #evaluate checkboxes
   $CTL_OPT{repeat_protein} = '' if(! $q->param('go_runner'));
   $CTL_OPT{model_org} = '' if(! $q->param('go_masker'));
   $CTL_OPT{snaphmm} = '' if(! $q->param('go_snap'));
   $CTL_OPT{gmhmm} = '' if(! $q->param('go_genemark'));
   $CTL_OPT{fgenesh_par_file} = '' if(! $q->param('go_fgenesh'));
   $CTL_OPT{augustus_species} = '' if(! $q->param('go_augustus'));

   #build predictors
   my @predictors = ($q->param('go_snap'),
		     $q->param('go_augustus'),
		     $q->param('go_fgenesh'),
		     $q->param('go_genemark'),
		     $q->param('est2genome'),
		     $q->param('protein2genome')
		     );
   @predictors = grep {$_} @predictors;
   push(@predictors, 'model_gff') if($CTL_OPT{model_gff});
   push(@predictors, 'pred_gff') if($CTL_OPT{pred_gff});
   $CTL_OPT{predictor} = join(",", @predictors);

   #get length of fasta file from db
   my $length = $self->get_length_for_value($CTL_OPT{genome});

   #get namen for job
   my $j_name = $self->get_name_for_value($CTL_OPT{genome});

   #if job exist update else make a new one
   if(my ($owner) = $self->dbh->selectrow_array(qq{SELECT user_id FROM jobs WHERE job_id=$job_id})){
       die "ERROR: This job does not belong to you\n" if($owner != $user_id);

       #update job
       $self->dbh->do(qq{UPDATE jobs SET length='$length', is_queued=$is_queued, name='$j_name', is_saved=$is_saved WHERE job_id=$job_id});
   }
   else{
       #add job
       $self->dbh->do(qq{INSERT INTO jobs (job_id, user_id, length, is_queued, is_started, }.
		      qq{is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
		      qq{is_tutorial, cpus, start_time, finish_time, name)}.
		      qq{VALUES($job_id, $user_id, '$length', $is_queued, 0, 0, 0, 0, 0, $is_saved, 0, 0, 0, '', '', '$j_name')}
		      );
   }

   #if ctl_opt exists update else make new entry
   if($self->dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$job_id})){
       #update control options for job
       my @defaults = (keys %CTL_OPT); #keys to add
       my @set = map {"$_ \= '$CTL_OPT{$_}'" } @defaults;
       $self->dbh->do("UPDATE ctl_opt SET ".join(", ", @set) . "WHERE job_id=$job_id");
   }
   else{
       #add control options for job
       my @defaults = (keys %CTL_OPT); #keys to add
       $self->dbh->do(qq{INSERT INTO ctl_opt (job_id, }.join(", ", @defaults).qq{) }.
		      qq{VALUES ($job_id, \'}.join("', '", @CTL_OPT{@defaults}).qq{\')}
		      );
   }

   $self->dbh->commit();

   if($later){ #save and move to frontpage
       return $self->frontpage("Job saved but not added to queue");
   }

   if($is_queued){ #save enqueue and move to front page
       return $self->frontpage("Job added to queue successfully");
   }

   #just go back to creating job
   return $self->job_create();
}
#---------------------------------------------------------------------------
#    my $sql_job= qq{INSERT INTO job (iduser,
# 				 description,
# 				 isstarted,
# 				 isended,
# 				 iscomplication,
# 				 iscomplicationNote,
# 				 note,
# 				 starttime,
# 				 endtime,
# 				 elapsetime,
# 				 totaltime,
# 				 isqued,
# 				 isflagged,
# 				 isevaluator,
# 				 jobstatus,
# 				 data_name,
# 				 data_format,
# 				 jobtype)
#                  VALUES($UID,
# 			'',
# 			0,
# 			0,
# 			0,
# 			0,
# 			'',
# 			NULL,
# 			NULL,
# 			0,
# 			0,
# 			0,
# 			0,
# 			0,
# 			'edit',
# 			'',
# 			'---',
# 			'$OT'
# 		       )
# 	      };
#---------------------------------------------------------------------------
sub feedback {
   my $self = shift;
   
   my $q = $self->query();
   my $feedback = $q->param('comment_text');
   
   use Mail::Sender;
   my $sender = Mail::Sender->new({smtp => 'm2.genetics.utah.edu'});
   
   my $mm = ("Maker Web Service\n".
	     "User feed back:\n".
	     "--------------------------------------\n".
	     "$feedback\n".
	     "--------------------------------------\n"
	    );
   
   my $sq = $sender->MailMsg({to      => "hadi\@genetics.utah.edu",
			      from    => "noreply\@genetics.utah.edu",
			      subject => "MWS:user Feed back",
			      msg     => $mm
			     });
      
   my $message="Your feed back has been sent.";
   
   return $self->tt_process('action_success.tt',{message=>$message});
}
#--------------------------------------------------------------------
sub update_user{
   my $self = shift;

=head1
   # Get CGI query object
   my $session = $self->session;
   my $q = $self->query();
	


   #	use Data::Dumper;
   #	my $st= "<br>";
   my $view = $q->param('pass');
   #	$st.=Dumper($session);
   #	return $st;

   my $UID= $self->get_user_id();

   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   #	my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore';

   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";








   if ($q->param('profile') eq "Change Profile") {

      my $student_id      = $q->param('student_id');
      my $first_name      = $q->param('first_name');
	    
	    
	    

	    
	    
      if ($first_name) {
      } else {
	 return $self->tt_process('general_add_error.tt',{message=>"First Name:$first_name is required"});
      }
	    
      my $last_name       = $q->param('last_name');
	    
      if ($last_name) {
      } else {
	 return $self->tt_process('general_add_error.tt',{message=>"Last Name: is required"});
      }
	    

      my $c_id            = $q->param('c_id');
	    
      if ($c_id) {
      } else {
	 return $self->tt_process('general_add_error.tt',{message=>"login: is required"});
      }
	



      my $e_mail          = $q->param('e_mail');
	    
      if ($e_mail) {
      } else {
	 return $self->tt_process('general_add_error.tt',{message=>"Email: is required"});
      }
	


	    
      my $same_login_sql="SELECT login                                                                                                                                                
           FROM user where login=\'".$c_id."\'";
      my @lg;

      my $sth = $dbh->prepare($same_login_sql);
      $sth->execute();

      my $ret;


      my @lg;


      while ( my $href = $sth->fetchrow_hashref ) {
	
	 return $self->tt_process('general_add_error.tt',{message=>"login: is already taken"});

		
	 last;
		
      }


      my $sql_i="update user                                                                                                                                                      
                               set first=\'$first_name\' ,                                                                                                                          
                               last=\'$last_name\',                                                                                                                                 
                               e_mail=\'$e_mail\',                                                                                                                                  
                               login= \'$c_id\'                                                                                                                                    
                                                                                                                                                                                    
                                                                                                                                                                                    
                      where iduser=$UID";
      #return $self->tt_process('action_success.tt',{message=>$sql_i});

      my $sth = $dbh->prepare($sql_i);
      my $rv = $sth->execute();
	    
      return $self->tt_process('action_success.tt',{message=>"Your profile has been updated"});

   }

	



   my $password        = $q->param('password');
	
   if ($password) {
   } else {
      return $self->tt_process('general_add_error.tt',{message=>"password: is required"});
   }
	    
	    
	    
   my $password_verify        = $q->param('password_verify');
	    
   if ($password_verify) {
   } else {
      return $self->tt_process('general_add_error.tt',{message=>"password verify: is required"});
   }
	    
	    
	    
   if ($password eq $password_verify) {
		
      #    return "<br>password matches";
   } else {
      return $self->tt_process('general_add_error.tt',{message=>"password: does not match"});
		
   }
	    
	    



   my $sql_up="update  user  set pass=\'$password\' where iduser=$UID";

   #	return $self->tt_process('action_success.tt',{message=>$sql_up});

   #print "<br>$sql_up";die; 
   my $sth = $dbh->prepare($sql_up);
   my $rv = $sth->execute();
   return $self->tt_process('action_success.tt',{message=>"Your password is updated"});

=cut
}
#-----------------------------------------------------------------------------
sub maker_add_one{

   my $self = shift;
   my $q = $self->query();

   my $noclue = $self->session->param('AUTH_USERNAME');

   my $session=$self->session;
   return $self->tt_process('maker_add_step1.tt',{
						  noclue=>$noclue,
						  session => $session});
}
#-----------------------------------------------------------------------------
sub maker_add_trained{

   my $self = shift;
   my $q = $self->query();
   my $noclue = $q->param('noclue');

   my $session=$self->session;
   return $self->tt_process('maker_add_trained.tt',{
						    noclue=>$noclue,

						    session => $session});


}
#------------------------------------------------------------------------------
sub maker_add_reannote{

   my $self = shift;
   my $q = $self->query();
   my $noclue = $q->param('noclue');

   my $session=$self->session;
   return $self->tt_process('maker_add_reannote.tt',{
						     noclue=>$noclue,
						     session => $session});
}
#-----------------------------------------------------------------------------
sub maker_add_notrained{

   my $self = shift;
   my $q = $self->query();
   my $noclue = $q->param('noclue');
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
 
 
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};
 


   #my @students = CourseDB::Student->retrieve_all();


   #    my @files=@{$self->get_file_box_entries()};

   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   #my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore';
   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";




   my $files_sql="SELECT distinct filename  FROM filebox  where iduser=$UID";



   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my @files;
   while ( my $href = $sth->fetchrow_hashref ) {

      push @files,$href->{filename};
   }

   my $count=$#files+1;






   my $session=$self->session;
   return $self->tt_process('maker_add_notrained.tt',{
						      noclue=>$noclue,
						      files=>\@files,
						      session => $session});
}
#-----------------------------------------------------------------------------
#returns the captcha image
sub create_captcha {
   my $self = shift;

   return $self->captcha_create;
}
#----------------------------------------------------------------------------
sub removeF{
   my $self = shift;
   my $errors = shift;

=head1
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
 
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};

   my $q = $self->query();
   my $fn = $q->param('fn');
   my $noclue=$q->param('noclue');

   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';

   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";

   my $del="delete from  filebox where iduser=\'$UID\' and filename=\'$fn\'";
   my $sth = $dbh->prepare($del);
   $sth->execute();

   my $ti="insert into  trash (iduser,filename)values ($UID,\'$fn\')";
   my $sth = $dbh->prepare($ti);
   $sth->execute();

   return $self->maker_file_view();    
=cut
}
#-----------------------------------------------------------------------------
sub removeJob{
   my $self = shift;
   my $errors = shift;

   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
 
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};

   my $q = $self->query();
   my $jn = $q->param('jn');
   my $noclue=$q->param('noclue');

   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   my $dbh = DBI->connect($connect_string)
     or die "<br>Got error $DBI::errstr when connecting to database\n";
    
   my $del="delete from  job  where iduser=$UID and idjob=$jn";
   my $sth = $dbh->prepare($del);
   $sth->execute();

   return $self->students_list();
}
#-----------------------------------------------------------------------------
sub maker_job_view{

=head1
   my $self = shift;
   my $errors = shift;

   my $q = $self->query();
   my $noclue = $q->param('noclue');
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';

   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";

   my $sqlhmm="select * from snaphmm";
   my $sth = $dbh->prepare($sqlhmm);
   $sth->execute();

   my @hmms;
   my $ret;
   while ( my $href = $sth->fetchrow_hashref ) {
      push @hmms,$href;
      $ret.= "\n".Dumper(\@hmms);
   }

   my $sqlaug="select * from augustus_species";
   my $sth = $dbh->prepare($sqlaug);
   $sth->execute();

   my @agu;
   my $ret;
   while ( my $href = $sth->fetchrow_hashref ) {
      push @agu,$href;
      $ret.= "\n".Dumper(\@agu);
   }

   my $sqlmodel="select * from model_org";
   my $sth = $dbh->prepare($sqlmodel);
   $sth->execute();

   my @morg;
   my $ret;
   while ( my $href = $sth->fetchrow_hashref ) {
      push @morg,$href;
      $ret.= "\n".Dumper(\@morg);
   }

   my $sqlpred="select * from predictor";
   my $sth = $dbh->prepare($sqlpred);
   $sth->execute();

   my @pred;

   while ( my $href = $sth->fetchrow_hashref ) {
      push @pred,$href;
      $ret.= "\n".Dumper(\@pred);
   }

   my $session=$self->session;

   return $self->tt_process('maker_add.tt',{snaphmms=>\@hmms,
					    predictor=>\@pred,
					    model_org=>\@morg,
					    agustus=>\@agu,
					    noclue=>$noclue,
					    session => $session,
					    errors  => $errors});
=cut
}
#------------------------------------------------------------------------------
sub maker_file_view_process{

=head1
   my $self = shift;
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
 
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};

   my $dir = "/home/apache/MWS/$UID/";
   unless(-d $dir){
      mkdir $dir or "Unable to create directory: ";
   }

   my $file_box_entry="default";
   my $file_box_desc="default";

   my $q = $self->query();

   my $jobid=$q->param('jid');

   my $FT=$q->param('FT');
   my $VL=$q->param('VL');

   #print 
   my $dumper=Dumper($self);
   #return "<br>chicka<hr>".$dumper;
   if ($q->param('gzip') eq 'Upload') {
      my $filename = $q->param('comp');
      my $output_file = $dir.$filename;
      #remove slashes
      my $file_desc=$q->param('file_desc');
    
      $file_box_entry = $filename;
      $file_box_desc =$file_desc;
    
      my $upload_filehandle = $q->upload('comp');
    
      open UPLOADFILE, ">>$output_file" ;
    
      binmode UPLOADFILE;
    
      while ( <$upload_filehandle> ) {
	 print UPLOADFILE;
      }
    
      close UPLOADFILE;
      if ($FT eq "gzip") {
	 use Archive::Tar;
	 my $tar = Archive::Tar->new;
	
	 $tar->read($output_file);
	 my $next = Archive::Tar->iter( $output_file, 1 );
	
	 my $nunu;
	
	 my $istar=0;
	 while ( my $f = $next->() ) {
	    $nunu .= $f->name;
	    $istar++;	    
	    my $buf=$tar->get_content($f->name);
	    open UPLOADFILE, ">>$output_file".$f->name ;
	    binmode UPLOADFILE;
	    print UPLOADFILE $buf;
	    close UPLOADFILE;
	 }
	
	 if ($istar < 1) {
	    return $self->tt_process('general_add_error.tt',
				     {message=>"Could not retrieve at least one file from compressed source. Must be a valid tar gzip archived file<hr><br><br>File:".$output_file."<hr><br>Please Examine the file and try again.Possible causes:<br>"
				     });	    
	 }

	 my $t=Dumper($tar);
	 return "<br>$nunu.yes";
      }
    
      if ($VL =~ /Fasta/) {
	 my $i=0;
	 my $output_file = $output_file;
	
	 open(my $OUT, $output_file) or die "Can't open  $output_file for writing\n$!";
	
	 my $isa=0;
	
	 while (<$OUT>) {
	    if (/^>\w/) {
	       $isa++;
	       if ($isa==1) {
		  print $OUT $_;
		  $i++;
	       } else {
		  return $self->tt_process('general_add_error.tt',{message=>"Could not validate as a fasta file<hr><br><br>File:".$output_file."<hr><br>Line".$i."+1: ".$_."<br>Please Examine the file and try again.Possible causes:<br>Your previous fasta header did not have any sequence. <br>Does not start with >seqid format<br>Must have neocleotides from the standard table"});
		    
	       }
		
	    } else {
	       if (/^[A-IK-NP-Z*]+$/i) {
		  print $OUT $_;
		  $i++;
		  $isa=0;
	       } else {
		  return $self->tt_process('general_add_error.tt',{message=>"Could not validate as a fasta file<hr><br><br>File:".$output_file."<hr><br>Line  ".$i."+1: ".$_."<br>Please Examine the file and try again.<br>You should not have any GAP(-)<br>Must have neocleotides from the standrd table"});
	       }
	    }
	 }
	 close $OUT;
      }
   }








   if ($q->param('copy') eq 'Save') {




      my $copy_name = $q->param('copy_name');
      my $copy_text=$q->param('copy_text');


      $file_box_entry=$copy_name;
      $file_box_desc="copy and paste";


    
      if ($copy_name) {
	 my $out_file = $dir.$copy_name;
	 open(my $OUT, '>', $out_file) or die "Can't open est $out_file for writing\n$!";
        


	 print $OUT $copy_text;
	 close $OUT;
      }

      #my $ss=Dumper($dir);
      #return "<br>$ss";
   }

   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   #my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore';

   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";

   my $sql_filebox="INSERT INTO filebox (iduser,filename,filedesc
				    )

				     VALUES
					  ($UID,\'$file_box_entry\',\'$file_box_desc\')";
					

   my $sth = $dbh->prepare($sql_filebox);
   my $rv = $sth->execute();
   my $file_id = $dbh->{insertid};
   $dbh->commit();






   #fileview after add
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   #my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore';
   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";




   my $files_sql="SELECT DISTINCT filename,filedesc  FROM filebox  where iduser=$UID";



   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my @files;
   while ( my $href = $sth->fetchrow_hashref ) {

      push @files,$href;
   }

   my $count=$#files+1;
   
   #return "<br>@files\n";
   #print add description with files
   return $self->tt_process('fileview.tt', {files => \@files,
					    UID      =>$UID,
					    count    =>$count,
					    first    =>$first_name.'   ',
					    lastt    =>$last_name,
					  
					    email    =>$e_mail,,
					    jobid   =>$jobid,
					    session  => $self->session});


   #get nd show done
   #x-
=cut
}
#-----------------------------------------------------------------------------
sub maker_user_view{
   my $self = shift;
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
    
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};
    
   my $q = $self->query();
    
   my $jid=$q->param('jid');
   my $jobid=$jid;
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
    
   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";
     
   my $files_sql="SELECT login  FROM user  where iduser=$UID";
    
   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my $cid = $sth->fetchrow_hashref;
   my $session = $self->session;
    
   #return "<br>$files_sql";
   return $self->tt_process('user_account.tt', {UID      =>$UID,
						first    =>$first_name,
						lastt    =>$last_name,
						email    =>$e_mail,
						jobid    =>$jobid,
						cid      =>$cid->{login},
						session  => $session});
}
#-----------------------------------------------------------------------------
sub maker_feedback_view{
   my $self = shift;
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
    
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};
    
   my $q = $self->query();
    
   my $jid=$q->param('jid');
   my $jobid=$jid;
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
    
   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";
    
   my $files_sql="SELECT login  FROM user  where iduser=$UID";
    
   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my $cid = $sth->fetchrow_hashref;
   my $session = $self->session;

   return $self->tt_process('user_feedback.tt', {UID      =>$UID,				  
						 first    =>$first_name,
						 lastt    =>$last_name,
						 email    =>$e_mail,
						 jobid    =>$jobid,
						 cid      =>$cid->{login},
						 session  => $session});
}
#-----------------------------------------------------------------------------
sub maker_help_view{
   my $self = shift;
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
    
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};
    
   my $q = $self->query();
    
   my $jid=$q->param('jid');
   my $jobid=$jid;
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';
   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";

   my $files_sql="SELECT login  FROM user  where iduser=$UID";

   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my $cid = $sth->fetchrow_hashref;
   my $session = $self->session;
    
   return $self->tt_process('help.tt', {UID      =>$UID,
					first    =>$first_name,
					lastt    =>$last_name,
					email    =>$e_mail,
					jobid    =>$jobid,
					cid      =>$cid->{login},
					session  => $session});
}

#-----------------------------------------------------------------------------
sub validate_file_view{
   my $self = shift;
   my $uref=$self->get_user_id("last");
   my $UID=$uref->{iduser};
    
   my $first_name=$uref->{first};
   my $last_name=$uref->{last};
   my $e_mail=$uref->{e_mail};
    
   my $q = $self->query();
    
   my $jid=$q->param('jid');
   my $jobid=$jid;
   my $connect_string = 'dbi:mysql:dbname=makerweb;user=bmoore;host=derringer.genetics.utah.edu';

   my $dbh = DBI->connect($connect_string)
     or die "Got error $DBI::errstr when connecting to database\n";
    
   my $files_sql="SELECT DISTINCT filename,filedesc  FROM filebox  where iduser=$UID";
    
   my $sth = $dbh->prepare($files_sql);
   $sth->execute();
    
   my @files;
   while ( my $href = $sth->fetchrow_hashref) {
	
      push @files,$href;
   }
    
   my $count=$#files+1;
    
   my $session = $self->session;
    
   return $self->tt_process('validate_file_view.tt', {files => \@files,
						      UID      =>$UID,
						      count    =>$count,
						      first    =>$first_name,
						      lastt    =>$last_name,
						      email    =>$e_mail,
						      jobid    =>$jobid,
						      session  => $session});
}
#-----------------------------------------------------------------------------
sub upload_file{
   my $self = shift;
   
   my $q = $self->query();
   my $from_job = $q->param('from_job'); #signal to use job create simple menu
   my ($filename) = $q->param("up_file") =~ /([^\/]+)$/;
   my $type = $q->param("type");
   my $name = $q->param("name") || $filename;
   my $fh = $q->upload('up_file');

   my $user_id = $self->get_user_id();
   my %serv_opt = %{$self->param("server_opt")};

   File::Path::mkpath("$serv_opt{data_dir}/users/$user_id");
   my $value = "$serv_opt{data_dir}/users/$user_id/$filename";

   open(my $OUT, "> $value");
   while(my $line = <$fh>){
       print $OUT $line;
   }
   close($OUT);
   close($fh);

   #should validate here
   #####################


   #get contig count and file length
   my $length = ($type eq 'fasta') ? MWAS_util::fasta_length($value) : @{[stat($value)]}[7];
   my $contig_count = ($type eq 'fasta') ? `grep -c ">" $value` : 0;
   chomp ($contig_count); #just in case


   $self->dbh->do(qq{INSERT INTO files (name, value, type, length, contig_count, user_id) }.
		  qq{VALUES ('$name', '$value', '$type', '$length', '$contig_count', $user_id)});

   $self->dbh->commit();

   if($from_job){
       $self->redirect("maker.cgi?rm=job_create&job_id=$from_job");
   }

   return $self->filebox;
}#-----------------------------------------------------------------------------
sub filebox{
   my $self = shift;
   
   my $q = $self->query();
   my $xfile = $q->param('xfile');
   my $from_job = $q->param('from_job'); #signal to use job create simple menu

   #delete the file indicated by user
   if($xfile){
       $self->dbh->do(qq{DELETE FROM files WHERE value='$xfile'});
       $self->dbh->commit();
       unlink($xfile) if(-e $xfile);
   }

   my $files = $self->get_all_files_info();

   return $self->tt_process('filebox.tt', {files => $files,
				           from_job => $from_job});
}


#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------
#this method collects all example/tutorial jobs for menu building
#and returns an array reference
sub get_tutorials {
   my $self  = shift;

   #get tutorial menu options
   my $tut = $self->dbh->selectall_arrayref(qq{SELECT name, job_id FROM jobs WHERE is_tutorial=1},
					    {Slice => {}}
					    );
   
   #make a value key and set it equal to job_id (for menu compatability)
   grep {$_->{value} = $_->{job_id}} @$tut;

   return [sort {$a->{name} cmp $b->{name}} (@$tut)];
}
#-----------------------------------------------------------------------------
#this method collects all menu names for a specific category
#and returns an array reference
sub get_menus {
   my $self  = shift;
   my $type  = shift;
   my $source = shift || 'all';
   my $uid   = shift;

   die "ERROR: Invalid source value '$source'" if($source !~ /^all$|^server$|^user$/);

   #get global menu options
   my $global = [];
   $global = $self->dbh->selectall_arrayref(qq{SELECT name, value, is_tutorial, length FROM menus WHERE type='$type'},
					    {Slice => {}}
					   ) unless($source eq 'user');
   
   #get user file based menu options
   my $user = [];
   $user = $self->dbh->selectall_arrayref(qq{SELECT name, value, length FROM files WHERE type='$type'}.
					  qq{ AND user_id='$uid'},
					  {Slice => {}}
					 ) unless($source eq 'server' || ! $uid);
    
   return [sort {$a->{name} cmp $b->{name}} (@$global, @$user)];
}
#-----------------------------------------------------------------------------
#this method collects values pointed to by menu/file names for a type
#and returns an array reference
sub get_value_for_name {
   my $self = shift;
   my $name = shift;
   my $type = shift;
   my $source = shift || 'all';
   my $id = $self->get_user_id();

   die "ERROR: Invalid source value '$source'" if($source !~ /^all$|^server$|^user$/);

   #get value from user file options
   my ($val) = $self->dbh->selectrow_arrayref(qq{SELECT value FROM files WHERE name='$name'}.
					      qq{ AND type='$type'AND user_id='$id'}
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
#this method collects values pointed to by menu/file names for a type
#and returns an array reference
sub get_name_for_value {
   my $self = shift;
   my $value = shift;

   #get value from user file options
   my ($nam) = $self->dbh->selectrow_arrayref(qq{SELECT name FROM files WHERE value='$value'});

   #user values always override system values
   return $nam if($nam);

   #get value from global menu options
   ($nam) = $self->dbh->selectrow_array(qq{SELECT name FROM menus WHERE value='$value'});
   
   return $nam;
}
#-----------------------------------------------------------------------------
#this method collects the fasta length for a file from the DB
sub get_length_for_value {
   my $self = shift;
   my $value = shift;

   #get value from user file options
   my ($len) = $self->dbh->selectrow_array(qq{SELECT length FROM files WHERE value='$value'});

   #user values always override system values
   return $len if($len);

   #get value from global menu options
   ($len) = $self->dbh->selectrow_array(qq{SELECT length FROM menus WHERE value='$value'});
   
   return $len;
}
#-----------------------------------------------------------------------------
#this method collects all information on the default ctl_opt from the database
#and returns a hash reference with all info
sub get_server_default_options {
   my $self = shift;

   my $def_opt = $self->dbh->selectrow_hashref(qq{SELECT * FROM all_default_opt});
    
   return $def_opt;
}
#-----------------------------------------------------------------------------
#this method collects all opt status
#and returns a hash reference with all info
sub get_stat {
   my $self = shift;

   my $stat = $self->dbh->selectrow_hashref(qq{SELECT * FROM stat});
    
   return $stat;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_user_info {
   my $self = shift;

   my $username = $self->session->param('AUTH_USERNAME');
   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE login='$username'});
    
   return $info;
}#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_user_info {
   my $self = shift;

   my $username = $self->session->param('AUTH_USERNAME');
   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE login='$username'});
    
   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user files
#and returns a hash reference with all infofile
sub get_file_info {
   my $self = shift;
   my $value = shift || return;

   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM files WHERE value=$value});
    
   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user files
#and returns a hash reference with all infofile
sub get_job_info {
   my $self = shift;
   my $job_id = shift || return;

   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM jobs WHERE job_id=$job_id});
    
   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user files
#and returns a hash reference with all infofile
sub get_all_files_info {
   my $self = shift;
   my $user_id = $self->get_user_id();

   my $info = $self->dbh->selectall_arrayref(qq{SELECT * FROM files WHERE user_id=$user_id},
					     {Slice => {}}
					    );

   foreach my $file (@$info){
       ($file->{filename}) = $file->{value} =~ /([^\/]+)$/;
   }
 
   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub get_user_id {
   my $self = shift;
    
   return $self->get_user_info()->{user_id};
}
#------------------------------------------------------------------------
sub get_file_box_entries{
   my $self = shift;
   my $id = $self->get_user_id();

   my $files = $self->dbh->selectall_arrayref(qq{SELECT * FROM files WHERE user_id='$id'},
					      {Slice => {}}
					     );
   
   return $files;    
}
#-----------------------------------------------------------------------------
1;
