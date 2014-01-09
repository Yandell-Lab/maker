package MWS;

BEGIN{
   if (not ($ENV{CGL_SO_SOURCE})) {
      $ENV{CGL_SO_SOURCE} = "$FindBin::Bin/lib/CGL/so.obo";
   }
   if (not ($ENV{CGL_GO_SOURCE})) {
      $ENV{CGL_GO_SOURCE} = "$FindBin::Bin/lib/CGL/gene_ontology.obo"
   }
   #$ENV{'CAP_DEVPOPUP_EXEC'} = 1;
   
   #must give group members access - i.e. apache group
   umask 0002;
}

use strict;
use warnings;

use base qw(CGI::Application);
use CGI::Application::Plugin::Session;
use CGI::Application::Plugin::Authentication;
use CGI::Application::Plugin::CAPTCHA;
use CGI::Application::Plugin::DBH (qw/dbh_config dbh/);
use CGI::Application::Plugin::TT;
#use CGI::Application::Plugin::DevPopup;
#use CGI::Application::Plugin::DevPopup::HTTPHeaders;
use CGI::Application::Plugin::Redirect;
use LWP::UserAgent;

use FindBin;
use GI;
use Data::Dumper;
use MWAS_util;
use URI::Escape;
use Mail::Sender;
use File::Basename;

#-----------------------------------------------------------------------------
sub cgiapp_init {
   my $self = shift;
   $self->SUPER::cgiapp_init;

   #load the server control files
   my %serv_opt;
   my $c_dir = MWAS_util::config_loc();

   if(-f "$c_dir/server.ctl"){
       %serv_opt = GI::set_defaults('server', {GI::parse_ctl_files(["$c_dir/server.ctl"])});
   }
   else{
       %serv_opt = GI::set_defaults('server', {GI::parse_ctl_files(["c_dir/server.ctl"])});
   }

   #make sure required database values are setup
   if (! $serv_opt{DBI}) {
      die "ERROR: You must specify a DBI connection method in: $c_dir/server.ctl\n\n";
   }
   if (! $serv_opt{dbname}) {
      die "ERROR: You must specify a database to connect to in: $c_dir/server.ctl\n\n";
   }
    
   #connect to the database
   my $dsn = "DBI:$serv_opt{DBI}:dbname=";
   $dsn .= "$serv_opt{data_dir}/" if($serv_opt{DBI} eq 'SQLite');
   $dsn .= "$serv_opt{dbname};";
   $dsn .= "host=$serv_opt{host};" if($serv_opt{host});
   $dsn .= "port=$serv_opt{port};" if($serv_opt{host} && $serv_opt{port});

   $self->dbh_config($dsn, $serv_opt{username}, $serv_opt{password}, {RaiseError => 1}) 
     or die "Got error $DBI::errstr when connecting to database\n";

   #reload default server options from server
   my %CTL_OPT = %{$self->get_server_default_options()}; #this is all control options
   @serv_opt{keys %serv_opt} = @CTL_OPT{keys %serv_opt}; #just get server options

   #setup template params
   $self->tt_config(TEMPLATE_OPTIONS => {INCLUDE_PATH => ["$FindBin::Bin/tt_templates/", "$c_dir/"],
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
					frontpage
					feedback
					filebox
					upload_file					
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
    my $q = $self->query;
    my $in_id = $q->param('guest_id');
    my $use_login = $self->param('server_opt')->{use_login};
    
    #authentication values
    my $is_authen = $self->authen->is_authenticated;
    my $user = $self->get_user_info;
    
    #false authentication, there is no user
    if($is_authen && ! $user){
	$self->authen->logout();
	$is_authen =$self->authen->is_authenticated;
    }

    #if no login is used then just use the administrator account
    if(! $use_login && (! $self->authen->is_authenticated || $self->get_user_id != 1)){
        my $user = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE user_id = 1});

        #authenticate user
        $q->param(authen_username => $user->{login});
        $q->param(authen_password => $user->{password});
        $self->authen->{initialized} = 0; #force reauthentication
        $self->authen->initialize();

        #return requested runmode as an authenticated user
        my $rm = $q->param('rm');
        $self->prerun_mode($rm) if($rm);
    }
    #if given guest_id for auto-login
    elsif($in_id && (! $self->authen->is_authenticated || $in_id != $self->get_user_id)){
        my $user = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE user_id=$in_id})
            if($in_id =~ /^\d+$/);
	
        #if not a valid guest make sure to reset authentication rather than returning active account
        if(! $user || !$user->{is_guest}){
            $self->authen->logout();
        }
        else{ #authenticate guest user
            $q->param(authen_username => $user->{login});
            $q->param(authen_password => $user->{password});
            $self->authen->{initialized} = 0; #force reauthentication
            $self->authen->initialize();

            #return requested runmode as an autheniticated user
            my $rm = $q->param('rm');
            $self->prerun_mode($rm) if($rm);
        }
    }

    $self->set_global_tt_params();
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
					   font  => $self->param('server_opt')->{font_file},
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
sub set_global_tt_params{
    my $self = shift;

    #set template values
    my $is_authen =$self->authen->is_authenticated;
    my $user = $self->get_user_info;

    #now get url information
    my $web_address = $self->param('server_opt')->{web_address};
    my $html_web    = $self->param('server_opt')->{html_web};
    my $cgi_web     = $self->param('server_opt')->{cgi_web};

    my $base_url_html = ($html_web =~ /http\:\/\//) ?
	"$html_web" : "$web_address/$html_web";
    $base_url_html =~ s/([^\:])\/+/$1\//g;
    
    my $base_url_cgi = ($cgi_web =~ /http\:\/\//) ?
	"$cgi_web" : "$web_address/$cgi_web";
    $base_url_cgi =~ s/([^\:])\/+/$1\//g;
    
    $self->param(base_url_html => $base_url_html);
    $self->param(base_url_cgi => $base_url_cgi);
    
    $self->tt_params({logged_in     => $is_authen,
		      server_opt    => $self->param('server_opt'), #server options file
		      session       => $self->session,
		      user          => $user,
		      base_url_html => $base_url_html,
		      base_url_cgi  => $base_url_cgi,
		      #authen       => Dumper($self->authen),
		  });
}

#-----------------------------------------------------------------------------
sub login {
   my $self = shift;
   my $add_text= shift;

   #already authenticated
   if($self->authen->is_authenticated){
       return $self->frontpage;
   }
   else{ #returns login screen
       return $self->tt_process('maker_login.tt', {message => $add_text});
   }
}
#-----------------------------------------------------------------------------
sub delete_job {
   my $self = shift;

   my $q = $self->query;
   my $job_id = $q->param('job_id');
   my $dest = $q->param('goto');
   my $user = $self->get_user_info();
   my $job = $self->get_job_info($job_id);

   if($job){
       die "ERROR: This is not your job to delete" if($user->{user_id} != $job->{user_id});
       
       $self->dbh->do(qq{DELETE FROM jobs WHERE job_id=$job_id});
       $self->dbh->do(qq{DELETE FROM ctl_opt WHERE job_id=$job_id});
       $self->dbh->commit();
       my %serv_opt = %{$self->param('server_opt')};
       File::Path::rmtree("$serv_opt{data_dir}/job/$job_id");
   }

   my $user_id = $self->get_user_id();
   return $self->queue if($dest && $dest eq 'queue');
   return $self->frontpage;
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

   #if no smtp server I can't send the username and password
   if(! $serv_opt->{smtp_server}){
       return $self->frontpage;
   }

   #evaluate captcha
   my $captcha = $q->param("captcha");
   if(! $self->captcha_verify($q->cookie("hash"), $captcha)){
       return $self->tt_process('forgot_login.tt', {message => "Value did not match image"});
   }

   #look up user
   my $user = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE e_mail='$e_mail'}) if($e_mail);

   #send mesage
   my $mm;
   if(defined $user && defined $user->{login}){
       $mm = "MWAS: Forgotten login/password request\n\n".
	     "Username: ".$user->{login}."\n".
	     "Password: ".$user->{password}."\n";

       MWAS_util::send_message($user->{e_mail},
				  $serv_opt->{smtp_server},
				  'MWAS: Forgotten login/password request',
				  $mm
				 );

       return $self->login("Your username and password have been e-mailed to you");    
   }
   elsif(defined $e_mail){
       return $self->tt_process('forgot_login.tt', {message => "No account associated with e-mail address"});
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
    my $user_id =  $q->param('user_id') || $self->get_user_id();
    my $value = $q->param('contig'); #datastore direcory for contig 
    $value =~ s/\/+$//;
    my ($name) = $value =~ /([^\/]+)$/;
    my %serv_opt = %{$self->param('server_opt')};
    my $data_dir = $serv_opt{data_dir};
    my $cgi_dir = $serv_opt{cgi_dir};
    my $c_dir = MWAS_util::config_loc();

    #GFF3 file
    my $gff_file = "$data_dir/jobs/$job_id/$job_id.maker.output/$value/$name.gff";
    if(! -f $gff_file){
	$name = uri_escape($name, '.');
	$gff_file = "$data_dir/jobs/$job_id/$job_id.maker.output/$value/$name.gff";
    }    

    #build a safename with '%' character escaped to avoid issues with browser interpretation of %
    my $safe_name = uri_escape($name, '_');
    $safe_name =~ s/\%/\_/g;
    
    #make path
    File::Path::mkpath("$serv_opt{html_dir}/users/$user_id/$job_id");
    my $gff_url = ($serv_opt{html_web} =~ /http\:\/\//) ?
	"$serv_opt{html_web}/users/$user_id/$job_id/$safe_name.gff" :
	"$serv_opt{web_address}/$serv_opt{html_web}/users/$user_id/$job_id/$safe_name.gff";
    
    #fix // in direcory structure
    $gff_url =~ s/([^\:])\/+/$1\//g;
    
    my $slink = "$serv_opt{html_dir}/users/$user_id/$job_id/$safe_name.gff";
    symlink($gff_file, $slink) if(! -e $slink);
    
    if($q->param('apollo')){
	#base URL for Apollo jars/images
	my $codebase = ($serv_opt{html_web} =~ /http\:\/\//) ? 
	    "$serv_opt{html_web}/" :
	    "$serv_opt{web_address}/$serv_opt{html_web}/";
	
	#fix // in direcory structure
	$codebase =~ s/([^\:])\/+/$1\//g;
	$codebase .= '/' if($codebase !~ /\/$/);
	
	#make jnlp file and URL
	my $jnlp_file = $slink;
	my $jnlp_url = $gff_url;
	$jnlp_file =~ s/gff$/jnlp/;
	$jnlp_url =~ s/gff$/jnlp/;
	
	#build jnlp content from template
	my $content = ${$self->tt_process('apollo.jnlp.tt', {codebase => $codebase,
							     gff_url => $gff_url})};
	
	open(JNLP, "> $jnlp_file");
	print JNLP $content;
	close(JNLP);
	
	return $self->redirect("$jnlp_url");       
    }
    elsif($q->param('soba')){
        my ($base, $name) = $gff_file =~ /(.*\/)([^\/]+)$/;
	if($name =~ /\%/){
	   $name =~ s/\%/\_/g;
	   my $link = $base . $name;
	   symlink($gff_file, $link);
	   $gff_file = $link;
	}
	
	#post the file to SOBA
	my $ua = LWP::UserAgent->new;
	my $response = $ua->post($serv_opt{soba_url},
				 Content_Type => 'form-data',
				 Content      => [ rm  => 'upload_files',
						   gff_file   => [$gff_file]]
	    );

	#fix the retunred content to remove relative URLs
	#my $content = $response->content;
	#my ($soba_base) = $serv_opt{soba_url} =~ /((http\:\/\/)?[^\/]+)/;
	#$content =~ s/\"\/([^\"])/\"$soba_base\/$1/g;
	
	#return the HTML content provided by SOBA
	#SOBA server must have 'Access-Control-Allow-Origin: *' in the
	#soba.cgi headers for cross site AJAX to work in FireFox
	#return $content;

	return $self->redirect("$serv_opt{soba_url}/?rm=reload_files&gff_file=$name");
    }
    elsif($q->param('gbrowse')){
	#show contig in GBrowse
	my $url = "$serv_opt{web_address}/cgi-bin/gb2/gbrowse/";
	$url =~ s/([^\:])\/+/$1\//g;
	$url .= "MWAS_$user_id\_$job_id/?name=$name";
	
	return $self->redirect($url);
    }
    elsif($q->param('jbrowse')){
	#show contig in JBrowse
	my $j_dir = $serv_opt{JBROWSE_ROOT};
	my $dir = "$serv_opt{html_dir}/users/$user_id/$job_id/";

	die "ERROR: JBROWSE_ROOT $j_dir does not exist\n" if(! -d $j_dir);
	
	#copy necessary JBrowse files if not yet copied
	File::Path::mkpath("$dir") if(! -d $dir);

	    
	#get all JBrowse contents
	my @files = map {File::Basename::basename($_)} <$j_dir/*>;
	my @to_copy = map {"$j_dir/$_" if(! -e "$dir/$_")} grep {$_ ne 'data'} @files;
	
	#get MAKER specific configuration file
	system("cp -R ".join(' ', @to_copy)." $dir");

	#add tracks if not currently added
	if(!-d "$dir/data"){
	    my $dstore = "$data_dir/jobs/$job_id/$job_id.maker.output/$job_id\_master_datastore_index.log";
	    system("cd $dir; ./bin/maker2jbrowse -d $dstore 1>&2");
	}

	my $url = ($serv_opt{html_web} =~ /http\:\/\//) ?
	    "$serv_opt{html_web}/users/$user_id/$job_id/" :
	    "$serv_opt{web_address}/$serv_opt{html_web}/users/$user_id/$job_id/";
	
	return $self->redirect($url);
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

       next unless(@data >= 3);

       $values{$data[0]} = $data[1] if($data[2] eq 'FINISHED');
       $status{$data[0]} = $data[2] unless($status{$data[0]} && $status{$data[0]} eq 'FINISHED');
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
   my ($type) = $self->dbh->selectrow_array(qq{SELECT type FROM jobs WHERE job_id=$job_id});
   my ($c_count) = $self->dbh->selectrow_array(qq{SELECT contig_count FROM files WHERE value='$genome'});
   ($c_count) = $self->dbh->selectrow_array(qq{SELECT contig_count FROM menus WHERE value='$genome'}) if(! $c_count);
   
   $counts{TOTAL} = $c_count if($c_count);

   return $self->tt_process('results.tt', {menus => \@menus,
					   counts => \%counts,
					   user_id => $user_id,
					   type => $type,
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

   #if no registration is allower redirect back to login
   if(! $serv_opt->{allow_register}){
       return $self->home_login;
   }

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

       my ($new_id) = $self->dbh->selectrow_array(qq{SELECT last_user_id FROM id_store}); #get last new user_id
       $new_id++; #iterate the value
       $self->dbh->do(qq{UPDATE id_store SET last_user_id=$new_id}); #record new value
       $self->dbh->do(qq{INSERT INTO users (user_id, login, password, first, last, e_mail, }.
		      qq{institution, is_guest, is_admin, date_created, last_login) VALUES ($new_id, }.
		      qq{'$login', '$password', '$first', '$last', '$e_mail', '$institution', 0, 0, \'}.
		      MWAS_util::date_time().qq{\', '')});
       $self->dbh->commit();


       #authenticate user for auto-login
       $q->param(authen_username => $login);
       $q->param(authen_password => $password);
       $self->authen->{initialized} = 0; #force reauthentication
       $self->authen->initialize();

       #forward user to account frontpage
       return $self->redirect("maker.cgi");
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

   #if not using login then you can't edit accounts
   if(! $serv_opt->{use_login}){
       return $self->frontpage;
   }

   #no input this must be a new request
   if(!$login && !$first && !$last && !$institution &&
      !$e_mail && !$old_password && ! $new_password && !$verify
     ){
       return $self->tt_process('edit_account.tt');
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
       return $self->tt_process('edit_account.tt', {errors => \%errors});
   }
   else{ #everything is fine, update user in database
       my $password = ($new_password) ? $new_password : $old_password; #select password to use
       $self->dbh->do(qq{UPDATE users SET login='$login', first='$first', last='$last', e_mail='$e_mail', }.
		      qq{password='$password', institution='$institution', is_guest=0 WHERE user_id=}.
		      $user->{user_id});
       $self->dbh->commit();

       #force reauthentication since login may have changed
       $q->param(authen_username => $login);
       $q->param(authen_password => $password);
       $self->authen->{initialized} = 0; #force reauthentication
       $self->authen->initialize();
       $self->set_global_tt_params; #reload tt_template options because of changes

       #update user info
       return $self->tt_process('edit_account.tt', {message => 'Account updated successfully',});
   }
}
#-----------------------------------------------------------------------------
sub guest_login {
    my $self = shift;

    my $serv_opt = $self->param('server_opt');

    return $self->home_login if(! $serv_opt->{allow_guest}); #redirect back to login if guest is not allowed

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

    #authenticate user, auto-login
    $self->query->param(authen_username => $login);
    $self->query->param(authen_password => $password);
    $self->authen->{initialized} = 0; #force reauthentication
    $self->authen->initialize();

    my $user_id = $self->get_user_id();    
    #forward user to account frontpage
    return $self->redirect("maker.cgi");
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
       return $self->home_login;
   }

   my $user = $self->get_user_info();

   #catch instance of false true on is_authenticated
   if(! $user){
       $self->authen->logout();
       return $self->frontpage();
   }

   my $dsn = "SELECT * FROM jobs WHERE user_id=".$user->{user_id}.
             " AND is_tutorial=0 AND is_saved=1 ORDER BY job_id DESC";
   my $jobs = $self->dbh->selectall_arrayref($dsn, {Slice => {}});

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
       elsif($job->{is_error}){
	   $job->{status} = 'idle';
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
					     message => $message});
}
#-----------------------------------------------------------------------------
sub queue {
   my $self = shift;
   my $message = shift;

   my $jobs = $self->dbh->selectall_arrayref(qq{SELECT * FROM jobs WHERE (is_queued=1 or is_running=1) }.
					     qq{and admin_block=0 and is_error=0 and }.
					     qq{is_finished=0 ORDER BY submit_id},
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
       elsif($job->{is_error}){
	   $job->{status} = 'idle';
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

   return $self->tt_process('queue.tt', {jobs => $jobs,
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
   my $func = $q->param('func');

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

   #the functional annotation option is only for finished jobs
   if($func && (! $job || ! $job->{is_packaged})){
       $func = 0;
   }
   elsif($job && $job->{type} eq 'functional'){
       $func = 1;
   }

   #log the old job_id for submiting new functional jobs
   my $old_id = '';
   if($func && $job->{type} ne 'functional'){
       $old_id = $job_id;
       $job_id = undef;
       $job = undef;
   }
   elsif($func){
       ($old_id) = $self->dbh->selectrow_array(qq{SELECT old_id FROM jobs WHERE job_id=$job_id});
   }

   #reserve new job id value if needed (tutorials are started with new job_id)
   if(! $job_id || $job->{is_tutorial} || $job->{is_started} || $job->{is_finished}){
       ($job_id) = $self->dbh->selectrow_array(qq{SELECT last_job_id FROM id_store}); #get last job_id
       $job_id++; #iterate the value
       $self->dbh->do(qq{UPDATE id_store SET last_job_id=$job_id}); #record new value
       $self->dbh->commit();
   }
   
   my %menus;   #collect options for drop down menus
   my $o_type; #organism type

   if(! $func){
       #decide on organism type
       $o_type = ($STAT{organism_type} !~ /STATIC|DISABLED/) ? $LOG_OPT{organism_type} : $CTL_OPT{organism_type};
       
       #server menu options (menu options follow control files)
       $menus{alt_peptide}         = [qw(A C D E F G H I K L M N P Q R S T V W Y)];
       $menus{model_org}{server}   = $self->get_menus('model_org', 'server');
       $menus{snaphmm}{server}     = $self->get_menus('snaphmm', 'server');
       $menus{gmhmm}{server}     = $self->get_menus('gmhmm_e', 'server') if($o_type eq 'eukaryotic');
       $menus{gmhmm}{server}     = $self->get_menus('gmhmm_p', 'server') if($o_type eq 'prokaryotic');
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
       $menus{gmhmm}{user}  = $self->get_menus('gmhmm', 'user', $user->{user_id});
       $menus{fastas}{user}   = $self->get_menus('fasta', 'user', $user->{user_id});
       $menus{gff3}{user}     = $self->get_menus('gff3', 'user', $user->{user_id});
       
       #example/tutorial menu options
       $menus{tutorials} = $self->get_tutorials();
   }

   return $self->tt_process('job_create.tt',{menus   => \%menus,
					     ctl_opt => \%CTL_OPT, #default values
					     log_opt => \%LOG_OPT, #logged values
					     stat    => \%STAT, #opt status values
					     o_type  => $o_type,
					     message => $message,
					     job_id  => $job_id,
					     old_id  => $old_id,
					     func    => $func,
					     job     => $job});
}
#------------------------------------------------------------------------------
sub submit_to_db {
   my $self = shift;

   my $q = $self->query;
   my $user_id = $self->get_user_id;
   my $is_queued = $q->param('add') ? 1 : 0;
   my $later = $q->param('later') ? 1 : 0; #save and come back later
   my $job_id = $q->param('job_id');
   my $serv_opt = $self->param('server_opt');
   my $func     = $q->param('func');
   my $old_id   = $q->param('old_id');

   #decide if this is a maker structural of functional job
   my $type = ($func) ? 'functional' : 'maker';

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

   if($func){
       my $job_ctl = $self->dbh->selectrow_hashref(qq{SELECT * FROM ctl_opt WHERE job_id=$old_id});
       %CTL_OPT = (%CTL_OPT, %$job_ctl);
       delete($CTL_OPT{tmp}) if(defined $CTL_OPT{tmp}); #temp
       delete($CTL_OPT{aed_threshold}) if(defined $CTL_OPT{aed_threshold}); #temp
       delete($CTL_OPT{job_id}) if(defined $CTL_OPT{job_id}); #temp
   }

   #join CTL_OPT in array to be seperated by comma
   while(my $key = each %CTL_OPT){
       my @values = $q->param($key);
       @values = grep {$_ ne ''} @values; #filter empty values
       $CTL_OPT{$key} = join(',', @values) if(@values);
       $CTL_OPT{$key} = '' if(!defined($CTL_OPT{$key}))
   }

   #evaluate checkboxes
   $CTL_OPT{est2genome} = ($q->param('est2genome')) ? 1 : 0;
   $CTL_OPT{protein2genome} = 1 if(! $q->param('est2genome'));
   $CTL_OPT{repeat_protein} = '' if(! $q->param('go_runner'));
   $CTL_OPT{model_org} = '' if(! $q->param('go_masker'));

   #get length of fasta file from db
   my $length = $self->get_length_for_value($CTL_OPT{genome});

   #get name for job
   my $j_name = $self->get_name_for_value($CTL_OPT{genome});

   #get new submit_id for job submission
   (my $submit_id) = $self->dbh->selectrow_array(qq{SELECT last_submit_id FROM id_store}); #get last submit_id
   $submit_id++; #iterate the value
   $self->dbh->do(qq{UPDATE id_store SET last_submit_id=$submit_id}); #record new value
   $self->dbh->commit();

   #if ctl_opt exists update else make new entry
   if($self->dbh->selectrow_array(qq{SELECT job_id FROM ctl_opt WHERE job_id=$job_id})){
       #update control options for job
       my @defaults = (keys %CTL_OPT); #keys to add
       my @set = map {"\"".lc($_)."\" = '$CTL_OPT{$_}'" } @defaults;
       $self->dbh->do("UPDATE ctl_opt SET ".join(", ", @set) . "WHERE job_id=$job_id");
   }
   else{
       #add control options for job
       my @defaults = (keys %CTL_OPT); #keys to add
       my @lc_defaults = map {lc($_)} @defaults;
       my $ver = GI::version();
       $self->dbh->do(qq{INSERT INTO ctl_opt (job_id, }.join(", ", map {"\"$_\""} @lc_defaults).qq{) }.
		      qq{VALUES ($job_id, '}.join("\', \'", @CTL_OPT{@defaults}).qq{\')}
		      );
   }

   #if job exist update else make a new one
   if(my ($owner) = $self->dbh->selectrow_array(qq{SELECT user_id FROM jobs WHERE job_id=$job_id})){
       die "ERROR: This job does not belong to you\n" if($owner != $user_id);

       #update job
       $self->dbh->do(qq{UPDATE jobs SET submit_id=$submit_id, length='$length', is_queued=$is_queued, }.
		      qq{name='$j_name', is_saved=$is_saved WHERE job_id=$job_id});
   }
   else{
       $j_name .= " - Post Processing" if ($j_name !~ /Post Processing$/ && $type eq 'functional');

       #add job
       $self->dbh->do(qq{INSERT INTO jobs (job_id, user_id, submit_id, length, type, is_queued, }.
	    qq{is_started, is_running, is_finished, is_error, is_packaged, is_saved, admin_block, }.
	    qq{is_tutorial, cpus, start_time, finish_time, name) }.
	    qq{VALUES ($job_id, $user_id, $submit_id, '$length', '$type', $is_queued, 0, }.
	    qq{0, 0, 0, 0, $is_saved, 0, 0, 0, '', '', '$j_name')}
	   );

       $self->dbh->do(qq{UPDATE jobs SET old_id = $old_id where job_id=$job_id}) if($type eq 'functional');
   }

   #save and move to frontpage
   if($later){
       $self->dbh->commit();
       return $self->frontpage("Job saved but not added to queue");
   }
   #save enqueue and move to front page
   elsif($is_queued){
       #get control parameters for job
       my $job_ctl = $self->dbh->selectrow_hashref(qq{SELECT * FROM ctl_opt WHERE job_id=$job_id});

       #check if a finished job used these exact same settings
       if(my $other_job_id = MWAS_util::package_already_exists($self->dbh, $job_ctl, $user_id)){
	   #make directories for running maker
	   my $user_dir = $serv_opt->{data_dir}."/users/$user_id";
	   my $job_dir = $serv_opt->{data_dir}."/jobs/$job_id";
	   File::Path::mkpath($user_dir) if(! -d $user_dir);
	   File::Path::mkpath($job_dir) if(! -d $job_dir);

	   MWAS_util::copy_package($self->dbh, $other_job_id, $job_id);
	   $self->dbh->do(qq{UPDATE jobs SET is_queued=0, is_started=1, is_finished=1, is_packaged=1,}.
			  qq{ start_time='}.MWAS_util::date_time . qq{', finish_time='}.
			  MWAS_util::date_time . qq{' WHERE job_id=$job_id});
       }
       elsif($type eq 'functional'){
	   MWAS_util::copy_package($self->dbh, $old_id, $job_id);
       }

       $self->dbh->commit;
       return $self->frontpage("Job added to queue successfully");
   }
   #just go back to creating job
   else{
       $self->dbh->commit;
       return $self->job_create();
   }
}
#---------------------------------------------------------------------------
sub feedback {
   my $self = shift;
   my $serv_opt = $self->param('server_opt');

   #if there is no admin e-mail and smtp server then I can't send feedback
   if(! $serv_opt->{smtp_server} && ! $serv_opt->{admin_email}){
       return $self->frontpage;
   }

   my $q = $self->query();
   my $comment = $q->param('comment_text');
   my $user_email = $q->param('e_mail') || "noreply@".$serv_opt->{smtp_server};

   if($comment){
       my $user = $self->get_user_info;
       my $sender = Mail::Sender->new({smtp => $serv_opt->{smtp_server}});

       my $mm = "Maker Web Annotation Service - User Feedback:\n".
	        "---------------------------------------------\n".
		"$comment\n".
		"---------------------------------------------\n";
       
       if($user && ! $user->{is_guest}){
	   $mm .= "User: ".$user->{user_id}."\n".
	          "First: ".$user->{first}."\n".
		  "Last: ".$user->{last}."\n".
		  "Institution: ".$user->{institution}."\n".
		  "E-mail: ".$user->{e_mail}."\n";
       }
       
       my $sq = $sender->MailMsg({to      => $serv_opt->{admin_email},
				  from    => $user_email,
				  subject => "MWAS: User Feedback",
				  msg     => $mm
				  });
       
       return $self->frontpage("Your message has been sent. Thank you for your feedback.");
   }
   else{
       return $self->tt_process('feedback.tt');
   }
}
#-----------------------------------------------------------------------------
#returns the captcha image
sub create_captcha {
   my $self = shift;

   return $self->captcha_create;
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
       $self->frontpage;
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
   my ($nam) = $self->dbh->selectrow_array(qq{SELECT name FROM files WHERE value='$value'});

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

   my $ctl_opt = $self->dbh->selectrow_hashref(qq{SELECT * FROM all_default_opt});
   my %def_opt = (GI::set_defaults('opts'), GI::set_defaults('bopts'),
		  GI::set_defaults('exe'), GI::set_defaults('server'));

   #fix control options to account for database lowercase restrictions
   while(my $key = each %def_opt){
       if(! exists $ctl_opt->{$key} && exists $ctl_opt->{lc($key)}){
	   $ctl_opt->{$key} = $ctl_opt->{lc($key)};
	   delete $ctl_opt->{lc($key)};
       }
   }   

   return $ctl_opt;
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

   my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM users WHERE login='$username'})
       if($username);
    
   return $info;
}
#-----------------------------------------------------------------------------
#this method collects all information on the user from the database
#and returns a hash reference with all info
sub is_guest {
    my $self = shift;

    return $self->get_user_info()->{is_guest};
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
