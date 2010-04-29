package stream;

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
use CGI::Application::Plugin::DBH (qw/dbh_config dbh/);
use CGI::Application::Plugin::TT;
use CGI::Application::Plugin::DevPopup;
use CGI::Application::Plugin::DevPopup::HTTPHeaders;
use CGI::Application::Plugin::Redirect;
use CGI::Application::Plugin::Stream (qw/stream_file/);

use FindBin;
use GI;
use Data::Dumper;
use MWAS_util;
use URI::Escape;
use Digest::MD5;

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
    
    #reload default server options from server
    my %CTL_OPT = %{$self->get_server_default_options()}; #this is all control options
    @serv_opt{keys %serv_opt} = @CTL_OPT{keys %serv_opt}; #just get server options
    
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
				LOGIN_URL => $serv_opt{html_web},
				POST_LOGIN_URL => $serv_opt{html_web},
				LOGOUT_URL => $serv_opt{html_web},
				LOGIN_SESSION_TIMEOUT => {IDLE_FOR => '30m',
							  EVERY => '1d'
							  },
				);

    $self->authen->protected_runmodes(qw());

    #add default control options from server
    $self->param(server_opt => \%serv_opt);
}
#-----------------------------------------------------------------------------
sub cgiapp_prerun {
        my $self = shift;

        my $is_authen = $self->authen->is_authenticated;
        my $user = $self->get_user_info;

        #false authentication / invalid user, logout
        if($is_authen && ! $user){
            $self->authen->logout();
            $is_authen =$self->authen->is_authenticated;
        }

        $self->tt_params({logged_in  => $is_authen,
                          server_opt => $self->param('server_opt'), #server options file                                                                                                                                                        
                          session    => $self->session,
                          user       => $user,
                          #authen     => Dumper($self->authen)
			 });
}
#-----------------------------------------------------------------------------
sub setup {
   my $self = shift;

   $self->start_mode('stream');
   $self->run_modes([qw(stream)]);
}
#-----------------------------------------------------------------------------
sub teardown {
   my $self = shift;

   $self->dbh->disconnect if($self->dbh);
}
#-----------------------------------------------------------------------------
sub stream {
    my $self = shift;
    my $q = $self->query();
    my $type = $q->param('type');

    my $job_id = $self->query->param('job_id') || return;
    my $user_id = $self->query->param('user_id') || return;
    my $md5 = $self->query->param('m');
    my $data_dir = $self->param('server_opt')->{data_dir};
    
    #for security remove non-digit characters
    $job_id =~ s/[^\d]+//g;
    $user_id = ~ s/[^\d]+//g;

    #validate parameters
    return if(! $job_id);
    return if(! $user_id);
    my $job_info = $self->get_job_info($job_id);
    return if(! $job_info || $user_id != $job_info->{user_id});
    
    #authenticate via login params or md5 digest
    if($m && $job_info){
	$key =  Digest::MD5::md5_hex(grep {defined($_)} values %$job_info);
	return if($key != $m);
    }
    else{
	my $user_info = $self->get_user_info();
	return if(! $user_info || $user_id != $user_info->{user_id});
    }

    #send file for requested file type
    if($type eq 'log'){
	my $file = "$data_dir/jobs/$job_id/job.log";

	return $self->_stream($file) if(-e $file);
    }
    elsif($type eq 'tarball'){
	my $job_id = $self->query->param('job_id') || return;
	my $user_id = $self->query->param('user_id') || return;
	my $data_dir = $self->param('server_opt')->{data_dir};
	my $file = "$data_dir/jobs/$job_id/$job_id.maker.output.tar.gz";

	return $self->_stream($file) if(-e $file);
    }

    return;
}
#-----------------------------------------------------------------------------
sub _stream{
    my $self = shift;
    my $file = shift;

    if($self->stream_file($file)){
	return;
    }
    else{
	return $self->error_mode();
    }
}
#-----------------------------------------------------------------------------
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
sub get_job_info {
    my $self = shift;
    my $job_id = shift || return;

    my $info = $self->dbh->selectrow_hashref(qq{SELECT * FROM jobs WHERE job_id=$job_id});

    return $info;
}
#-----------------------------------------------------------------------------
1;
