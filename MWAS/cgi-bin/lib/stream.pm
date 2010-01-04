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
use CGI::Application::Plugin::Session;
use CGI::Application::Plugin::Authentication;
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
 
  $self->authen->protected_runmodes(qw());

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

    if($type eq 'img'){
	my $src = $self->query->param('src') || return;
	my $file = "images/$src";

	return $self->_stream($file);
    }
    elsif($type eq 'log'){
	my $job_id = $self->query->param('job_id') || return;
	my $user_id = $self->query->param('user_id') || return;
	my $data_dir = $self->param('server_opt')->{data_dir};
	my $file = "$data_dir/jobs/$job_id/job.log";

	return $self->_stream($file);
    }
    elsif($type eq 'tarball'){
	my $job_id = $self->query->param('job_id') || return;
	my $user_id = $self->query->param('user_id') || return;
	my $data_dir = $self->param('server_opt')->{data_dir};
	my $file = "$data_dir/jobs/$job_id/$job_id.maker.output.tar.gz";

	return $self->_stream($file);
    }
    elsif($type eq 'jnlp'){
	my $job_id = $self->query->param('job_id') || return;
	my $user_id = $self->query->param('user_id') || return;
	my $value = $self->query->param('value') || return;
	my $data_dir = $self->param('server_opt')->{data_dir};
	my $file = "$data_dir/jobs/$job_id/$job_id.maker.output/$value/apollo.jnlp";

	return $self->_stream($file);
    }
    elsif($type eq 'gff3'){
	my $job_id = $self->query->param('job_id') || return;
	my $user_id = $self->query->param('user_id') || return;
	my $value = $self->query->param('value') || return;
	$value =~ s/\/+$//;
	my ($name) = $value =~ /([^\/]+)$/;
	my $data_dir = $self->param('server_opt')->{data_dir};
	my $file = "$data_dir/jobs/$job_id/$job_id.maker.output/$value/$name.gff";

	return $self->_stream($file);
    }
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
1;
