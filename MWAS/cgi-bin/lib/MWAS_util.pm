package MWAS_util;

use strict;
use warnings;

use FindBin;
use GI;
use File::NFSLock;

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
