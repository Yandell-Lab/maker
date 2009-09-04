#! /usr/bin/perl -w

package Process::MpiTiers;

use strict;
use Error qw(:try);
use Error::Simple;
use Process::MpiChunk;
use Storable;

#--set object variables for serialization of data
#this is needed when cloning an MPIChunk object
$Storable::forgive_me = 1; #allows serializaion of objects with code refs

#-----------------------------------------------------------------------------
#------------------------------GLOBAL VARIABLES-------------------------------
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------

#returns a new MpiTier object
sub new {
   my ($class, @args) = @_;

   my $self = {};
   
   bless ($self, $class);
   
   if (@args) {
      my $arg = shift @args;
      if (ref $arg eq 'Process::MpiTiers') {
	 $self = $arg->clone();
      }
      else {
	  my $VARS           = $arg; #this should be a hash ref
	  $self->{TIER_ID}   = shift @args;
	  $self->{TERMINATE} = 0;
	  $self->{FAILED}    = 0;

	  #optionaly override chunk type
	  $self->{CHUNK_REF} = shift @args || "Process::MpiChunk";
	  $self->{CHUNK_REF} = new $self->{CHUNK_REF}; #turn into object

	  #setup the tier
	  $self->_initialize_level(0);
	  $self->_prepare($VARS);
      }
   }

   return $self;
}
#--------------------------------------------------------------
#initializes variables in the MpiTier object
#called from $self->new
#this is a good place to set up any data needed prior to the chunk level

sub _prepare {
   my $self = shift;
   my $VARS = shift;

   try{
      $self->{CHUNK_REF}->prepare($VARS);
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, 'Failed in tier preparation');
   };

   return if($self->failed);

   if(! defined $VARS->{c_flag}){
       die "ERROR: VARS->{c_flag} is not set in MpiTier.\n".
	   "This value should be set by ".ref($self->{CHUNK_REF})."::prepare\n".
	   "durring initialization of the MpiTier. Please\n".
	   "edit the ".ref($self->{CHUNK_REF})."::prepare code to set this\n".
	   "value appropriately\n\n";
   }

   if($VARS->{c_flag} <= 0){
      $self->{TERMINATE} = 1;
      return;
   }

   $self->{VARS} = $VARS;
}
#--------------------------------------------------------------
#defines variables for new level and undefines last level
sub _initialize_level {
   my $self = shift;
   my $level = shift;

   #undef previous level
   if(exists ($self->{LEVEL}{CURRENT})){
       my $last = $self->{LEVEL}{CURRENT};
       $self->{LEVEL}{LAST} = $last;
       $self->{LEVEL}{$last} = undef;
   }

   #terminate
   if (! defined $level){ #level will be undefined when finished
      $self->{TERMINATE} = 1;
      $self->{DS} = "$self->{VARS}{seq_id}\t$self->{VARS}{out_dir}\tFINISHED";
      $level = $self->{LEVEL}{CURRENT};
   }

   #set current level
   $self->{LEVEL}{CURRENT} = $level;

   #--initiate level variables
   $self->{LEVEL}{$level} = {};
   $self->{LEVEL}{$level}{CHUNK_COUNT} = 0;
   $self->{LEVEL}{$level}{RESULT_COUNT} = 0;
   $self->{LEVEL}{$level}{CHUNKS} = [];
   $self->{LEVEL}{$level}{ERROR} = undef;
   $self->{LEVEL}{$level}{STARTED} = 0;
   $self->{LEVEL}{$level}{FINISHED} = 0;
}
#--------------------------------------------------------------
#jump to indicated level
#same as _initialize_level, defined for syntactic sugar
sub _go_to_level {
    my $self = shift;
    my $level = shift;

    return $self->_initialize_level($level);
}
#--------------------------------------------------------------
#itteratively gets the next chunk from the tier.  It returns
#undef if there is no chunk, or if the tier cannot yet advance

sub next_chunk {
   my $self = shift;

   return undef if ($self->terminated || $self->failed);

   #handle levels that have no chunks to run
   while($self->_level_finished){
       $self->_next_level;
       $self->_load_chunks;
   }

   #handle case where level needs to be initialized
   $self->_load_chunks if(! $self->_level_started);

   #get level after doing any necessary moves to next level
   my $level = $self->{LEVEL}{CURRENT};

   if (my $chunk = shift @{$self->{LEVEL}{$level}{CHUNKS}}) {
      return $chunk;
   }
   else{
      return undef;
   }
}

#--------------------------------------------------------------
#continues running until multiple chunks are encountered

sub run {
   my $self = shift;

   return if ($self->terminated && ! $self->failed);
   return if ($self->_level_started && ! $self->_level_finished);

   $self->_next_level if ($self->_level_finished);
   $self->_load_chunks if (!$self->_level_started);

   while ($self->num_chunks == 1){
       my $chunk = $self->next_chunk;
       $chunk->run($self->id);
       $self->update_chunk($chunk);

       return if ($self->terminated && ! $self->failed);
       return if ($self->_level_started && ! $self->_level_finished);

       $self->_next_level if ($self->_level_finished);
       $self->_load_chunks if (!$self->_level_started);
   }
}
#--------------------------------------------------------------
#runs the tier to completion without stopping

sub run_all {
   my $self = shift;

   return if ($self->terminated || $self->failed);
   return if ($self->_level_started && ! $self->_level_finished);

   $self->_next_level if ($self->_level_finished);
   $self->_load_chunks if (!$self->_level_started);

   while(my $chunk = $self->next_chunk){
       $chunk->run($self->id);
       $self->update_chunk($chunk);
   }
}

#--------------------------------------------------------------
#This method moves from one level to another

sub _next_level {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   return if ($self->terminated || $self->failed);

   if(! $self->_level_started){
       $self->_load_chunks;
       return;
   }

   return if (! $self->_level_finished);

   my $next_level;

   try {
      $next_level = $self->{CHUNK_REF}->flow($level, $self->{VARS});
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, 'Can not get next level');
   };

   return if($self->failed);

   $self->_go_to_level($next_level);
}

#--------------------------------------------------------------
#builds chunks into the current level

sub _load_chunks {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   return if ($self->terminated || $self->failed);
   return if ($self->_level_started);

   my $chunks;

   try {
      $chunks = $self->{CHUNK_REF}->loader($level, $self->{VARS}, $self->id);
   }
   catch Error::Simple with {
      my $E = shift;

      $self->{LEVEL}{$level}{CHUNKS} = [];
      $self->{LEVEL}{$level}{CHUNK_COUNT} = $self->num_chunks;
      $self->{LEVEL}{$level}{STARTED} = 1; #avoids infinite loop

      $self->_handler($E, 'Can not load chunks');
   };

   return if($self->failed);

   $self->{LEVEL}{$level}{CHUNKS} = $chunks;
   $self->{LEVEL}{$level}{CHUNK_COUNT} = $self->num_chunks;
   $self->{LEVEL}{$level}{STARTED} = 1;
}

#--------------------------------------------------------------
#takes a chunk and adds its results to the tier

sub update_chunk {
   my $self = shift;
   my $chunk = shift;

   #get chuck id
   my $id = $chunk->id();
   my ($tier_id, $level_num, $chunk_num) = split (":", $id);

   #check if chunk goes with tier/level
   if($tier_id != $self->id()){
      die "ERROR:  This MpiChunk is not part of this MpiTier\n";
   }
   if($level_num != $self->{LEVEL}{CURRENT}){
      die "ERROR:  This MpiChunk is not part of this level\n";
   }

   #check run status
   if($chunk->failed){
      my $E = $chunk->exception;
      $self->_handler($E, "Chunk failed at level $level_num\n");
   }
   else{
       #let the chunk add results to $self->{VARS}
       $chunk->result($self->{VARS});
   }

   $self->{LEVEL}{$level_num}{RESULT_COUNT}++;    
   $self->{ERROR} .= $chunk->error();

   print STDERR $chunk->error();
}
#--------------------------------------------------------------
#returns true if level is finished

sub _level_finished {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   return undef if (! $self->_level_started);

   if($self->{LEVEL}{$level}{RESULT_COUNT} == $self->{LEVEL}{$level}{CHUNK_COUNT}){
       $self->{LEVEL}{$level}{FINISHED} = 1;
   }
   
   return $self->{LEVEL}{$level}{FINISHED};
}

#--------------------------------------------------------------
#returns true if level has been been started, i.e. had chunks in it

sub _level_started {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   return $self->{LEVEL}{$level}{STARTED};
}
#--------------------------------------------------------------
#returns true if tier has terminated

sub terminated {
   my $self = shift;

   #_level_finished is required because there may still be chunks
   #that must be gathered and destroyed before terminating the tier
   if ($self->failed && $self->_level_finished){
       $self->{TERMINATE} = 1;
   }

   return $self->{TERMINATE};
}

#-------------------------------------------------------------
#returns true if tier failed

sub failed{
   my $self = shift;
   return $self->{FAILED};
}
#-------------------------------------------------------------
#returns whatevever is strored in $self->{ERROR}
#this is a good place for logging information

sub error{
   my $self = shift;
   return $self->{ERROR} || '';
}
#-------------------------------------------------------------
#returns the value of the current level

sub current_level{
   my $self = shift;
   return $self->{LEVEL}{CURRENT};
}
#--------------------------------------------------------------
#this reports the number of remaining chunks.

sub num_chunks {
   my $self = shift;

   my $level = $self->{LEVEL}{CURRENT};
   my $num = @{$self->{LEVEL}{$level}{CHUNKS}};

   return $num;
}

#-------------------------------------------------------------
#returns the id of the tier

sub id {
   my $self = shift @_;
   return $self->{TIER_ID};
}
#--------------------------------------------------------------
#returns a copy of this MpiTiers object using the Storable
#perl module.  As a result globs and code refs are ignored.

sub clone {
   my $self = shift;
    
   return Storable::dclone($self);
}
#--------------------------------------------------------------
#returns a line giving name of sequence and the location where
#analysis is stored.  Used for datastore implementation.

sub DS {
   my $self = shift;
   return $self->{DS} || $self->{VARS}{seq_id} . "\t" . $self->{VARS}{out_dir};
}
#-------------------------------------------------------------
sub fasta{
   my $self = shift;
   return $self->{VARS}{fasta};
}
#--------------------------------------------------------------
#exception handler

sub _handler {
   my $self = shift;
   my $E = shift;
   my $extra = shift || '';

   my $level = $self->current_level;
   my $seq_id = $self->{VARS}{seq_id};
   my $LOG = $self->{VARS}{LOG};

   print STDERR $E->{-text};
   print STDERR "ERROR: ".$extra."!!\n";
   print STDERR "FAILED CONTIG:$seq_id\n\n" if(defined $seq_id);

   if(defined $LOG){
      my $die_count = $LOG->get_die_count();
      $die_count++;
      $LOG->add_entry("DIED","RANK", $self->id);
      $LOG->add_entry("DIED","COUNT",$die_count);
   }

   #clear queue on error
   while(my $chunk = $self->next_chunk){
       $self->{LEVEL}{$level}{RESULT_COUNT}++;
   }

   $self->{DS} = "$self->{VARS}{seq_id}\t$self->{VARS}{out_dir}\tDIED";
   $self->{FAILED} = 1;
}

#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
1;
