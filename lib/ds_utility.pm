#------------------------------------------------------------------------
#----                       datastore_utility                        ---- 
#------------------------------------------------------------------------
package ds_utility;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Datastore::MD5;
use File::Path;
use Cwd;
use URI::Escape;
use File::NFSLock;

@ISA = qw(
       );

my $CWD = Cwd::getcwd();

#------------------------------------------------------------------------
#--------------------------- METHODS   ----------------------------------
#------------------------------------------------------------------------
sub new {
   my $class = shift;
   my @args = @_;
   my $self = {};

   bless ($self, $class);

   $self->_initialize(@args);

   return $self;
}
#------------------------------------------------------------------------
sub _initialize {
   my $self = shift @_;
   my %CTL_OPTIONS = %{shift @_};

   $CWD = $CTL_OPTIONS{CWD} || $CWD;

   my $out_base = $CTL_OPTIONS{out_base} || $CWD;
   my $out_name = $CTL_OPTIONS{out_name} || "output";
   my $ds_flag  = (exists($CTL_OPTIONS{datastore})) ? $CTL_OPTIONS{datastore} : 1;

   $self->{root} = "$out_base/$out_name\_datastore";
   $self->{log} = "$out_base/$out_name\_master_datastore_index.log";
   
   print STDERR "A data structure will be created for you at:\n".
   $self->{root}."\n\n".
   "To access files for individual sequences use the datastore index:\n".
   $self->{log}."\n\n" unless($main::qq);

   if($ds_flag){   
      $self->{ds_object} = new Datastore::MD5('root' => $self->{root},
					      'depth' => 2
					     );
   }
   else{
      $self->{ds_object} = undef;
   }

   #initialize a new blank log,
   #except when using the hidden chpc option
   #then just append to the existing log
   unless($CTL_OPTIONS{_chpc}){
       open(my $IN, ">", $self->{log});
       close($IN);
   }
}
#------------------------------------------------------------------------
sub get_index {
   my $self = shift;
   
   return $self->{log};
}
#------------------------------------------------------------------------
sub get_root {
   my $self = shift;
   
   return $self->{root};
}
#------------------------------------------------------------------------
sub mkdir {
   my $self = shift;
   my $id = shift;

   my $safe_id = uri_escape($id, 
			    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			   );

   my $dir = $self->{root}."/".$safe_id;

   if($self->{ds_object}){
      $dir = $self->{ds_object}->id_to_dir($safe_id);
      $self->{ds_object}->mkdir($safe_id) || die "ERROR: could not make datastore directory\n";
   }
   else{
      File::Path::mkpath($dir);
   }

   return $dir;
}
#------------------------------------------------------------------------
sub id_to_dir {
   my $self = shift;
   my $id = shift;

   my $safe_id = uri_escape($id, 
			    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			   );

   my $dir = $self->{root}."/".$safe_id;

   if($self->{ds_object}){
      $dir = $self->{ds_object}->id_to_dir($safe_id);
   }

   return $dir;
}
#------------------------------------------------------------------------
sub seq_dirs {
   my $self = shift;
   my $id = shift;

   my $out_dir = $self->mkdir($id);
   my $safe_id = uri_escape($id, 
			    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			   );
   my $the_void = "$out_dir/theVoid.$safe_id";
   File::Path::mkpath($the_void);

   return $out_dir, $the_void;
}
#------------------------------------------------------------------------
sub add_entry {
   my $self = shift;
   my @F = @_;

   #remove deep directory data so log is relative
   my $cwd = ($CWD) ? $CWD : Cwd::getcwd();
   my $entry = join("\t", @F);

   #maker/mpi_iprscan specific if statement
   if($entry =~ /\tFINISHED|\tSTARTED|\tDIED|\tSKIPPED|\tRETRY|\tFAILED/){
       $entry =~ s/$cwd\/.*\.maker\.output\/*|$cwd\/.*\.iprscan\.output\/*//;
   }

   #lock file so no one else writes to it (MPI safe)
   if(my $lock = new File::NFSLock($self->{log}, 'EX', 5, 5)){
       open(my $IN, ">>", $self->{log});
       print $IN $entry . "\n";
       close($IN);
       $lock->unlock;
   }
   else{
       die "ERROR: ds_utility::add_entry method timed out\n\n";
   }

}
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
1;
