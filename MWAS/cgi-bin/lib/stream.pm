package stream;

BEGIN{
   if (not ($ENV{CGL_SO_SOURCE})) {
      $ENV{CGL_SO_SOURCE} = "$FindBin::RealBin/lib/CGL/so.obo";
   }
   if (not ($ENV{CGL_GO_SOURCE})) {
      $ENV{CGL_GO_SOURCE} = "$FindBin::RealBin/lib/CGL/gene_ontology.obo"
   }
   #{ $ENV{'CAP_DEVPOPUP_EXEC'} = 1; }
}

use strict;
use warnings;

use base qw(CGI::Application);
use CGI::Application::Plugin::DBH (qw/dbh_config dbh/);
use CGI::Application::Plugin::TT;
#use CGI::Application::Plugin::DevPopup;
#use CGI::Application::Plugin::DevPopup::HTTPHeaders;
use CGI::Application::Plugin::Redirect;
use CGI::Application::Plugin::Stream (qw/stream_file/);

use FindBin;
use GI;
use Data::Dumper;
use MWS;
use MWAS_util;
use URI::Escape;
use Digest::MD5;
use vars qw(@ISA);

@ISA = ('MWS');

#-----------------------------------------------------------------------------
sub cgiapp_init {
    my $self = shift;
    $self->SUPER::cgiapp_init;
   
#    my %serv_opt = %{$self->param('server_opt')};
#    
#    #setup authentication
#    __PACKAGE__->authen->config(DRIVER => ['DBI',
#					   DBH         => $self->dbh,
#					   TABLE       => 'users',
#					   CONSTRAINTS => {'users.login'    => '__CREDENTIAL_1__',
#							   'users.password' => '__CREDENTIAL_2__',
#						       }
#					   ],
#				STORE => 'Session',
#				LOGIN_URL => $serv_opt{html_web},
#				POST_LOGIN_URL => $serv_opt{html_web},
#				LOGOUT_URL => $serv_opt{html_web},
#				LOGIN_SESSION_TIMEOUT => {IDLE_FOR => '30m',
#							  EVERY => '1d'
#							  },
#				);
#
#    $self->authen->protected_runmodes(qw());
#
#    #add default control options from server
#    $self->param(server_opt => \%serv_opt);
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

    #hack - handle case for gbrowse calling stream and splitting on '='
    if(@ARGV && ! @{[grep {/\=/} @ARGV]}){
	my %params = @ARGV;
	map {$q->param($_ => $params{$_})} keys %params;
    }

    my $type = $q->param('type');
    my $job_id = $q->param('job_id');
    my $user_id = $q->param('user_id');
    my $md5 = $q->param('m');
    my $data_dir = $self->param('server_opt')->{data_dir};
    my $html_dir = $self->param('server_opt')->{html_dir};
    my $html_web = $self->param('server_opt')->{html_web};
    my $web_address = $self->param('server_opt')->{web_address};

    #for security remove non-digit characters
    $job_id =~ s/[^\d]+//g if($job_id);
    $user_id =~ s/[^\d]+//g if($user_id);

    my $job_info = $self->get_job_info($job_id) if($job_id);
    
    #authenticate via login params or md5 digest
#    if($md5 && $job_info){
#    	my $key =  Digest::MD5::md5_hex(grep {defined($_)} values %$job_info);
#	return if($key != $md5);
#    }
#    else{
#	my $user_info = $self->get_user_info();
#	return if(! $user_info || $user_id != $user_info->{user_id});
#    }

    #send file for requested file type
    if($type eq 'log'){
	#-return the contents of the MAKER produced STDERR log file
	return if(! $job_info || ! $user_id || $user_id != $job_info->{user_id});
	my $file = "$data_dir/jobs/$job_id/job.log";

	return '' if(! -e $file);

	my $name = "job_$job_id.log";
	my $new = "$html_dir/users/$user_id/$name";

	my $url = ($html_web =~ /http\:\/\//) ?
        "$html_web/users/$user_id/$name" :
        "$web_address/$html_web/users/$user_id/$name";

	File::Path::mkpath("$html_dir/users/$user_id/");

	symlink($file, $new);

	return $self->redirect($url);
    }
    elsif($type eq 'tarball'){
	#-returns the tarball of MAKER results

	return if(! $job_info || ! $user_id || $user_id != $job_info->{user_id});

	my $file = "$data_dir/jobs/$job_id/$job_id.maker.output.tar.gz";

	return $self->_stream($file) if(-e $file);
    }
    elsif($type eq 'gbrowse'){
	#-the gbrowse portion of stream is actually expected to
	#-be called outside of the webrowser so it will return
	#-STDOUT rather than a webpage and then exit

	return if(! $job_info || ! $user_id || $user_id != $job_info->{user_id});

	#create Gcd /vaBrowse configuration file for this Job
	my $content = ${$self->tt_process('gbrowse.conf.tt', {contig_dir => "$data_dir/jobs/$job_id/$job_id.maker.output/",
							      })};

	print $content;
	exit(0);
    }
    elsif($type eq 'file'){
	my $file = $q->param('value');

	return if(! $file || ! $user_id);

	my ($name) = $file =~ /([^\/]+)$/;
	my $new = "$html_dir/users/$user_id/$name";

	my $url = ($html_web =~ /http\:\/\//) ?
        "$html_web/users/$user_id/$name" :
        "$web_address/$html_web/users/$user_id/$name";

	File::Path::mkpath("$html_dir/users/$user_id/");

	symlink($file, $new);

	return $self->redirect($url);
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

    return if(! defined $username);

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
sub get_server_default_options {
    my $self = shift;

    my $def_opt = $self->dbh->selectrow_hashref(qq{SELECT * FROM all_default_opt});

    return $def_opt;
}
#-----------------------------------------------------------------------------
1;
