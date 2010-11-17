#! /usr/bin/perl -w

package Process::MpiTiers;

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../perl/lib";
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
	  $self->{INTERRUPT} = 0;

	  #optionaly override chunk type
	  $self->{CHUNK_REF} = shift @args || "Process::MpiChunk";
	  eval "require ". $self->{CHUNK_REF};
	  $self->{CHUNK_REF} = $self->{CHUNK_REF}->new(); #turn into object

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

   $self->{VARS} = $VARS;

   try{
      $self->{CHUNK_REF}->_prepare($VARS);
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, 'Failed in tier preparation');
   };

   return if($self->failed);

   if($self->{CHUNK_REF}->_should_continue($self) <= 0){
      $self->_set_terminate(1);
      return;
   }
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
      $self->_set_terminate(1);
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
      $next_level = $self->{CHUNK_REF}->_flow($level, $self->{VARS});
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
      $chunks = $self->{CHUNK_REF}->_loader($level, $self->{VARS}, $self->id);
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
       $chunk->_result($self->{VARS});
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

   #set terminate if initialization step said not to continue
   #this is semi-redundant depending on how the MpiChunk::prepare
   #method was set up
   if($self->{CHUNK_REF}->_should_continue($self) <= 0){
       $self->_set_terminate(1);
   }

   #_level_finished is required because there may still be chunks
   #that must be gathered and destroyed before terminating the tier
   if ($self->failed && $self->_level_finished){
       $self->_set_terminate(1);
   }

   return $self->{TERMINATE};
}
#--------------------------------------------------------------
#sets $self->{TERMINATE}

sub _set_terminate {
   my $self = shift;
   my $arg = shift;

   if(! defined($arg) || $arg !~ /^0$|^1$/){
       die "FATAL:  Attempt to set TERMINATE to an illegal value\n".
	   "in Process::MpiTiers\n\n";
   }

   my $old = $self->{TERMINATE};
   $self->{TERMINATE} = $arg;

   if($arg == 1 &&(! defined($old) || $old != 1)){
       $self->_on_termination;
   }
}

#-------------------------------------------------------------
#returns true if tier failed

sub failed{
   my $self = shift;
   return $self->{FAILED};
}

#-------------------------------------------------------------
#sets $self->{FAILED}

sub _set_failed {
   my $self = shift;
   my $arg = shift;

   if(! defined($arg) || $arg !~ /^0$|^1$/){
       die "FATAL:  Attempt to set FAILED to an illegal value\n".
	   "in Process::MpiTiers\n\n";
   }

   my $old = $self->{FAILED};
   $self->{FAILED} = $arg;

   if($arg == 1 &&(! defined($old) || $old != 1)){
       $self->_on_failure;
   }
}
#-------------------------------------------------------------
#returns true if tier is interrupted

sub interrupt{
   my $self = shift;
   return $self->{INTERRUPT};
}
#-------------------------------------------------------------
#sets $self->{INTERRUPT}

sub _set_interrupt {
   my $self = shift;
   my $arg = shift;

   if(! defined($arg) || $arg !~ /^0$|^1$/){
       die "FATAL:  Attempt to set FAILED to an illegal value\n".
	   "in Process::MpiTiers\n\n";
   }

   my $old = $self->{INTERRUPT};
   $self->{INTERRUPT} = $arg;
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
#returns the MpiChunk::_on_failure method

sub _on_failure {
   my $self = shift;

   return $self->{CHUNK_REF}->_on_failure($self);
}
#--------------------------------------------------------------
#returns the MpiChunk::_on_termination method

sub _on_termination {
   my $self = shift;

   return $self->{CHUNK_REF}->_on_termination($self);
}
#-------------------------------------------------------------
#exception handler

sub _handler {
   my $self = shift;
   my $E = shift;
   my $extra = shift || '';

   my $level = $self->current_level;

   print STDERR $E->{-text};
   print STDERR "ERROR: ".$extra."!!\n" if($extra);

   #clear queue on error
   while(my $chunk = $self->next_chunk){
       $self->{LEVEL}{$level}{RESULT_COUNT}++;
   }

   $self->_set_failed(1);
}
#--------------------------------------------------------------

sub AUTOLOAD {
    my $self = shift;

    my $caller = caller();
    use vars qw($AUTOLOAD);
    my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
    $call =~/DESTROY/ && return;

    if (! exists $self->{VARS}{$call}) {
	die "Error: Invalid method \'$call\' in Process::MpiTiers\n".
	    "call to AutoLoader issued from: $caller\n\n";
    }

    if (@_) {
	$self->{VARS}{$call} = shift;
    }

    return $self->{VARS}{$call};
}
#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
1;
