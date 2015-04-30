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
use Carp;

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
   my $CTL_OPT = shift @_;

   $CWD = $CTL_OPT->{CWD} || $CWD;

   my $out_base = $CTL_OPT->{out_base} || $CWD;
   my $out_name = $CTL_OPT->{out_name} || "output";
   my $ds_flag  = (exists($CTL_OPT->{datastore})) ? $CTL_OPT->{datastore} : 1;

   $self->{out_base} = $out_base;
   $self->{root} = "$out_base/$out_name\_datastore";
   $self->{log} = "$out_base/$out_name\_master_datastore_index.log";
   $CTL_OPT->{SEEN_file} = $self->{SEEN_file} = "$out_base/seen.dbm";
   
   print STDERR "A data structure will be created for you at:\n".
   $self->{root}."\n\n".
   "To access files for individual sequences use the datastore index:\n".
   $self->{log}."\n\n" unless($main::qq);

   if($ds_flag){   
      carp "Calling Datastore::MD5::new" if($main::debug);
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
   if(!$CTL_OPT->{_multi_chpc} && !$CTL_OPT->{_resume}){
       open(IN, "> $self->{log}");
       open(IN, "> $self->{SEEN_file}");
       close(IN);
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
      carp "Calling Datastore::MD5::id_to_dir" if($main::debug);
      $dir = $self->{ds_object}->id_to_dir($safe_id);
      return $dir if(-d $dir); #avoid unecessary IO and system calls

      carp "Calling Datastore::MD5::mkdir" if($main::debug); 
     $self->{ds_object}->mkdir($safe_id) || die "ERROR: could not make datastore directory\n";
   }
   else{
      return $dir if(-d $dir); #avoid unecessary IO and system calls

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
      carp "Calling Datastore::MD5::id_to_dir" if($main::debug);
      $dir = $self->{ds_object}->id_to_dir($safe_id);
   }

   return $dir;
}
#------------------------------------------------------------------------
sub seq_dirs {
   my $self = shift;
   my $id = shift;

   carp "Calling Datastore::MD5::mkdir" if($main::debug);
   my $out_dir = $self->mkdir($id);
   carp "Calling uri_escape" if($main::debug);
   my $safe_id = uri_escape($id, 
			    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			   );
   my $the_void = "$out_dir/theVoid.$safe_id";
   carp "Calling File::Path::mkpath" if($main::debug);
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

   open(my $IN, ">>", $self->{log});
   print $IN $entry . "\n";
   close($IN);
}
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
1;
