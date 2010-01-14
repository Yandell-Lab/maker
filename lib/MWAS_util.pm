package MWAS_util;

use strict;
use warnings;

use FindBin;
use GI;
use File::NFSLock;
use Data::Dumper;
use File::Copy;

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
#check to see if a job with the exact same control file options exists
#and returns its id
sub package_already_exists{
    my $dbh = shift @_;
    my %CTL_OPT = %{shift @_};
    my $user_id = shift @_;

    my %def_opt = (GI::set_defaults('opts'), GI::set_defaults('bopts')); #get system produced CTL_OPT
    my @set = map {"ctl_opt.".lc($_)." \= '".$CTL_OPT{lc($_)}."'" } grep {!/gmhmm_e|gmhmm_p/i}keys %def_opt;


    my $dsn = "SELECT ctl_opt.job_id FROM ctl_opt JOIN jobs ON ctl_opt.job_id=jobs.job_id WHERE ".join(' AND ', @set)." AND jobs.is_packaged=1";
#    $dsn .= ($user_id) ? " AND (user_id=$user_id OR is_tutorial=1)" : " AND is_tutorial=1";

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
    mkdir($new_dir) if(! -d $new_dir);

    my @files =  <$job_dir/*>;

    @files = grep {!/\.tar\.gz$/} @files;

    system("cp -R ".join(' ', @files)." $new_dir/");
    @files = (<$new_dir/*/$job_old*>,<$new_dir/$job_old*>);
    foreach my $f (@files){
	my $new = $f;
	$new =~ s/\/$job_old([\.\_][^\/]+)$/\/$job_new$1/;
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
	   "tar -zcf $r_dir.tar.gz $r_dir --exclude \"run.log\" --exclude ".
	   "\"theVoid\*\" --exclude \"seen.dbm\" --exclude \"mpi_blastdb\"") &&
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

    my $lock = lockDB() || ((print STDERR lock_error()) && (return));
    $dbh->do($do_string);
    $dbh->commit;

    return 1;
}
#-----------------------------------------------------------------------------
1;
