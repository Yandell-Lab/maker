#! /usr/bin/perl -w

package Process::MpiTiers;

use strict;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../perl/lib";
use Error qw(:try);
use Error::Simple;
use Process::MpiChunk;
use Storable;
use Carp;

#--set object variables for serialization of data
#this is needed when cloning an MPIChunk object
$Storable::forgive_me = 1; #allows serializaion of objects with GLOBs
#$Storable::Deparse = 1; #now serializes CODE refs
#$Storable::Eval= 1;$ #now serializes CODE refs

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
	  $self->{VARS}      = {%$arg}; #forces copy of hash (1 level deep)
	  $self->{TIER_ID}   = shift @args || 0;
	  $self->{RANK}      = undef;
	  $self->{PARENT}    = $self->{TIER_ID};
	  $self->{TERMINATE} = 0;
	  $self->{FAILED}    = 0;
	  $self->{INTERRUPT} = 0;


	  #optionaly override chunk type
	  $self->{CHUNK_REF} = shift @args || "Process::MpiChunk";
	  eval "require ". $self->{CHUNK_REF};
	  $self->{CHUNK_REF} = $self->{CHUNK_REF}->new(); #turn into object

          #set up tier_type for tiers within tiers
          $self->{TIER_TYPE} = shift @args || 0;

	  #setup the tier
	  $self->_initialize_level(0);
	  $self->_prepare();
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
   my $VARS = $self->{VARS};
   my $tier_type = $self->{TIER_TYPE};

   try{
      $self->{CHUNK_REF}->_prepare($VARS, $tier_type);
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

   #undef forward levels when we move backwards
   if(exists ($self->{LEVEL}{CURRENT})){
       my $last = $self->{LEVEL}{CURRENT};
       $self->{LEVEL}{LAST} = $last;

       if(defined($level) && $last > $level){
	  for(my $i = $last; $i > $level; $i--){
	     $self->{LEVEL}{$i} = undef;
	  }
       }
   }

   #terminate
   if (! defined $level){ #level will be undefined when finished
      $self->_set_terminate(1);
      $level = $self->{LEVEL}{CURRENT};
      $self->{LEVEL} = undef;
      $self->{LEVEL}{CURRENT} = $level;
      $self->{LEVEL}{ALL}{CHUNK_COUNT} = 0;
      $self->{LEVEL}{ALL}{RESULT_COUNT} = 0;
      $self->{LEVEL}{ALL}{CHUNKS} = [];
      return;
   }

   #set current level
   $self->{LEVEL}{CURRENT} = $level;

   #--initiate level variables
   if(! $self->{LEVEL}{$level} || ! $self->_level_started($level) || $self->_level_finished($level)){
      $self->{LEVEL}{$level} = {};
      $self->{LEVEL}{$level}{CHUNK_COUNT} = 0;
      $self->{LEVEL}{$level}{RESULT_COUNT} = 0;
      $self->{LEVEL}{$level}{CHUNKS} = [];
      $self->{LEVEL}{$level}{ERROR} = undef;
      $self->{LEVEL}{$level}{STARTED} = 0;
      $self->{LEVEL}{$level}{FINISHED} = 0;
   }

   #--initiate all-level variables
   if(! $self->{LEVEL}{ALL}){
      $self->{LEVEL}{ALL}{CHUNK_COUNT} = 0;
      $self->{LEVEL}{ALL}{RESULT_COUNT} = 0;
      $self->{LEVEL}{ALL}{CHUNKS} = [];
   }

   #initiate earlier levels if I jump forward non-sequentially
   for(my $i = $level; $i >= 0; $i--){
      if(! $self->{LEVEL}{$i}){
	 $self->{LEVEL}{$i} = {};
	 $self->{LEVEL}{$i}{CHUNK_COUNT} = 0;
	 $self->{LEVEL}{$i}{RESULT_COUNT} = 0;
	 $self->{LEVEL}{$i}{CHUNKS} = [];
	 $self->{LEVEL}{$i}{ERROR} = undef;
	 $self->{LEVEL}{$i}{STARTED} = 1;
	 $self->{LEVEL}{$i}{FINISHED} = 0;
      }
   }
}
#--------------------------------------------------------------
#jump to indicated level
#same as _initialize_level, defined for syntactic sugar
sub go_to_level {
    my $self = shift;
    my $level = shift;

    return $self->_initialize_level($level);
}
#--------------------------------------------------------------
#itteratively gets the next chunk from the tier.  It returns
#undef if there is no chunk, or if the tier cannot yet advance

sub next_chunk {
   my $self = shift;
   my $what = shift; #tier or chunk

   confess "ERROR: Improper type specification in Process::MpiTiers:next_chunk\n"
       if($what && $what ne 'chunk' && $what ne 'tier');

   warn "WARNING: You must always set a rank before running MpiTiers\n"
       if(! defined $self->{RANK});

   return undef if ($self->terminated || $self->failed);

   #handle levels that have no chunks to run
   while($self->_level_finished){
       $self->_next_level;
       $self->_load_chunks;
       $self->_load_extra if($self->{LEVEL}{EXTRA});
       return undef if ($self->terminated || $self->failed); #check again
   }
   
   #handle case where level needs to be initialized
   $self->_load_chunks if(! $self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   #get level after doing any necessary moves to next level
   my $level = $self->{LEVEL}{CURRENT};

   #get chunk to return
   my $chunk;
   my @skip;
   while(($chunk = shift @{$self->{LEVEL}{ALL}{CHUNKS}})){
      last if(! $what);
      last if($what eq 'chunk' && ref($chunk) ne ref($self));
      last if($what eq 'tier' && ref($chunk) eq ref($self));
      
      push(@skip, $chunk);
      $chunk = undef;
   }

   unshift(@{$self->{LEVEL}{ALL}{CHUNKS}}, @skip);
   return $chunk;
}

#--------------------------------------------------------------
#fill in all chunks and levels to maximum without running

sub actualize {
   my $self = shift;
   my $rank = shift;

   $self->{RANK} = $rank if(defined($rank));
   warn "WARNING: You must always set a rank before running MpiTiers\n"
       if(! defined $self->{RANK});

   return if ($self->terminated || $self->failed);

   #handle levels that have no chunks to run
   while($self->_level_finished){
       $self->_next_level;
       $self->_load_chunks;
       $self->_load_extra if($self->{LEVEL}{EXTRA});
   }
   
   #handle case where level needs to be initialized
   $self->_load_chunks if(! $self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   #get level after doing any necessary moves to next level
   my $level = $self->{LEVEL}{CURRENT};
}

#--------------------------------------------------------------
#continues running until multiple chunks are encountered

sub run {
   my $self = shift;
   my $rank = shift;

   $self->{RANK} = $rank if(defined($rank));
   warn "WARNING: You must always supply a rank before running MpiTiers\n"
       if(! defined ($self->{RANK}));

   return if ($self->terminated && ! $self->failed);
   return if ($self->_level_started && ! $self->_level_finished);

   $self->_next_level if ($self->_level_finished);
   $self->_load_chunks if (!$self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   while ($self->num_chunks == 1 && ref($self->{LEVEL}{ALL}{CHUNKS}[0]) ne ref($self)){
       my $chunk = $self->next_chunk;
       $chunk->run($self->rank);
       $self->update_chunk($chunk);

       return if ($self->terminated && ! $self->failed);
       return if ($self->_level_started && ! $self->_level_finished);

       $self->_next_level if ($self->_level_finished);
       $self->_load_chunks if (!$self->_level_started);
       $self->_load_extra if($self->{LEVEL}{EXTRA});
   }
}
#--------------------------------------------------------------
#runs the tier to completion without stopping

sub run_all {
   my $self = shift;
   my $rank = shift;

   $self->{RANK} = $rank if(defined($rank));

   warn "WARNING: You must always supply a rank before running MpiTiers\n"
       if(! defined ($self->{RANK}));

   return if ($self->terminated || $self->failed);
   return if ($self->_level_started && ! $self->_level_finished);

   $self->_next_level if ($self->_level_finished);
   $self->_load_chunks if (!$self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   while(my $chunk = $self->next_chunk){
      if($chunk->can('run_all')){
	 $chunk->run_all($self->rank);
      }
      else{
	 $chunk->run($self->rank);
      }
      $self->update_chunk($chunk);
   }
}
#--------------------------------------------------------------
#This method moves from one level to another

sub _next_level {
   my $self = shift;
   my $VARS = $self->{VARS};
   my $level = $self->{LEVEL}{CURRENT};
   my $tier_type = $self->{TIER_TYPE};

   return if ($self->terminated || $self->failed);

   if(! $self->_level_started){
       $self->_load_chunks;
       $self->_load_extra if($self->{LEVEL}{EXTRA});
       return;
   }

   return if (! $self->_level_finished);

   my $next_level;

   try {
      $next_level = $self->{CHUNK_REF}->_flow($VARS, $level, $tier_type, $self);

      #---------------#temp
      #report the clock usage at each level in seconds
      #use Benchmark;
      #if($self->{BENCHMARK}){
      #    my $mark = new Benchmark;
      #    my $diff = timediff($mark, $self->{BENCHMARK});
      #    $self->{BENCHMARK} = $mark;
      #    print "LEVEL $level: Time taken was ", timestr($diff, 'all'), " seconds\n";
      #}
      #else{
      #    $self->{BENCHMARK} = new Benchmark;
      #}
      #---------------#temp
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, 'Can not get next level');
   };

   return if($self->failed);

   $self->go_to_level($next_level);
   
   #always reset global chunk variables when next level
   $self->{LEVEL}{ALL}{CHUNKS} = [];
   $self->{LEVEL}{ALL}{CHUNK_COUNT} = 0;
   $self->{LEVEL}{ALL}{RESULT_COUNT} = 0;
}

#--------------------------------------------------------------
#builds chunks into the current level

sub _load_chunks {
   my $self = shift;
   my $VARS = $self->{VARS};
   my $level = $self->{LEVEL}{CURRENT};
   my $tier_type = $self->{TIER_TYPE};

   return if ($self->terminated || $self->failed);
   return if ($self->_level_started);

   my $chunks;

   try {
      $chunks = $self->{CHUNK_REF}->_loader($VARS, $level, $tier_type, $self);
   }
   catch Error::Simple with {
      my $E = shift;

      $self->{LEVEL}{$level}{CHUNKS} = [];
      $self->{LEVEL}{$level}{CHUNK_COUNT} = $self->num_chunks($level);
      $self->{LEVEL}{$level}{STARTED} = 1; #avoids infinite loop
      $self->{LEVEL}{EXTRA} = 0;

      $self->_handler($E, 'Can not load chunks');
   };

   return if($self->failed);

   $self->{LEVEL}{$level}{CHUNKS} = $chunks;
   $self->{LEVEL}{$level}{CHUNK_COUNT} = $self->num_chunks($level);
   $self->{LEVEL}{$level}{STARTED} = 1;
   $self->{LEVEL}{EXTRA} = 1;
   $self->_load_extra;
}

sub _load_extra {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};
   my $tier_type = $self->{TIER_TYPE};

   return if ($self->terminated || $self->failed);

   #set up chunk pool
   if($self->{LEVEL}{EXTRA}){
      for(my $i = 0; $i <= $level; $i++){
	 $self->{LEVEL}{ALL}{CHUNK_COUNT} += @{$self->{LEVEL}{$i}{CHUNKS}};
	 push(@{$self->{LEVEL}{ALL}{CHUNKS}}, @{$self->{LEVEL}{$i}{CHUNKS}});
	 $self->{LEVEL}{$i}{CHUNKS} = [];
      }
      $self->{LEVEL}{EXTRA} = 0;

      #sort to specified order if any
      my $all = $self->{LEVEL}{ALL};
      
      try{
	  $all->{CHUNKS} = $self->{CHUNK_REF}->_sort_levels($all->{CHUNKS}, $level, $tier_type, $self)
	      if($self->{CHUNK_REF}->can('_sort_levels'));
      }
      catch Error::Simple with {
	  my $E = shift;
	  
	  $self->_handler($E, 'Failed sorting levels');
      };
   }
}
#--------------------------------------------------------------
#takes a chunk and adds its results to the tier

sub update_chunk {
   my $self = shift;
   my $chunk = shift;

   #get chuck id
   my $id = $chunk->id();
   my $parent = $chunk->parent();
   my ($tier_id, $level_num, $tier_type, $chunk_num) = split (":", $id);

   #handles older 3 level IDs
   if(! defined($chunk_num)){
      $chunk_num = $tier_type;
      $tier_type = 0;
   }

   #check if chunk goes with tier/level
   if($parent ne $self->id){
      confess "ERROR:  This MpiChunk does not belong to this MpiTier\n";
   }
   if(! $self->_level_started($level_num) || $self->_level_finished($level_num)){
      confess "ERROR:  This MpiChunk is not part of this level\n";
   }

   #check run status
   if($chunk->failed){
      my $E = $chunk->exception;
      $self->_handler($E, "Chunk failed at level:$level_num, tier_type:$tier_type");
   }
   elsif($chunk->terminated){
      #let the chunk add results to $self->{VARS}
      my $VARS = $self->{VARS};
      if(! $self->failed){
	  my $stat = $chunk->_result($VARS, $level_num, $tier_type, $self);
	  $chunk->_finalize($self) if($stat && $chunk->can('_finalize'));
      }
   }
   else{
       confess "ERROR: Logic error, attempt to update unfinished chunk\n";
   }

   $self->{LEVEL}{$level_num}{RESULT_COUNT}++;    
   $self->{LEVEL}{ALL}{RESULT_COUNT}++;    
   $self->{ERROR} .= $chunk->error();

   print STDERR $chunk->error();
}
#--------------------------------------------------------------
#method to allow Tiers to produce results like chunks
sub _result{
   my $self = shift;
   my $VARS = shift;            #this ia a hash reference
   my $level = shift;
   my $tier_type = shift;

   #only return results for finished/succesful chunks
   return if (! $self->terminated || $self->failed);

   #results set in $self->{CHUNK_REF}->_on_termination
   $self->{CHUNK_REF}->{RESULTS} = $self->{RESULTS};
   $self->{CHUNK_REF}->{VARS} = {};
   $self->{CHUNK_REF}->{FINISHED} = 1;

   #run result in CHUNK_REF
   my $stat;
   try {
       $stat = $self->{CHUNK_REF}->_result($VARS, $level, $tier_type, $self);
   }
   catch Error::Simple with {
       my $E = shift;
       
       $self->_handler($E, 'Failed gathering tier as chunk result');
   };

   #clear memory
   $self->{CHUNK_REF}->{RESULTS} = {};
   $self->{CHUNK_REF}->{VARS} = {};
   $self->{CHUNK_REF}->{FINISHED} = 0;

   return $stat;
}
#--------------------------------------------------------------
#returns true if level is finished

sub _level_finished {
   my $self = shift;
   my $level = shift;
   $level = $self->{LEVEL}{CURRENT} if(! defined($level));
   my $previous = $level - 1;

   return undef if (! $self->_level_started);
   return $self->{LEVEL}{$level}{FINISHED} if($self->{LEVEL}{$level}{FINISHED});

   if($self->{LEVEL}{$level}{RESULT_COUNT} == $self->{LEVEL}{$level}{CHUNK_COUNT} &&
      ($previous < 0 || $self->_level_finished($previous))
     ){
       $self->{LEVEL}{$level}{FINISHED} = 1;
   }
   
   return $self->{LEVEL}{$level}{FINISHED};
}

#--------------------------------------------------------------
#returns true if level has been been started, i.e. had chunks in it

sub _level_started {
   my $self = shift;
   my $level = shift;
   $level = $self->{LEVEL}{CURRENT} if(! defined($level));

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
   if ($self->failed &&
       ($self->_level_finished ||
	($self->chunk_total_count == $self->result_count))
      ){
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
       confess "FATAL:  Attempt to set TERMINATE to an illegal value\n".
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
       confess "FATAL:  Attempt to set FAILED to an illegal value\n".
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
       confess "FATAL:  Attempt to set FAILED to an illegal value\n".
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
#returns whatevever is strored in $self->{VARS}->{fasta}

sub fasta{
   my $self = shift;

   return $self->{VARS}->{fasta} || '';
}
#-------------------------------------------------------------
#returns whatevever is strored in $self->{VARS}->{q_def}

sub q_def{
   my $self = shift;

   return $self->{VARS}->{q_def} || '';
}
#-------------------------------------------------------------
#returns the value of the current level

sub current_level{
   my $self = shift;
   return $self->{LEVEL}{CURRENT};
}
#--------------------------------------------------------------
#this reports the number of chunks remaining to be processed

sub num_chunks {
   my $self = shift;
   my $level = shift;
   $level = 'ALL' if(! defined($level));

   #$self->_load_chunks if(! $self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   my $num = @{$self->{LEVEL}{$level}{CHUNKS}};

   return $num;
}
#--------------------------------------------------------------
#this reports the number of chunk results collected.

sub result_count {
   my $self = shift;
   my $level = shift;
   $level = 'ALL' if(! defined($level));

   #$self->_load_chunks if(! $self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   my $count = $self->{LEVEL}{$level}{RESULT_COUNT} || 0;

   return $count;
}
#--------------------------------------------------------------
#this reports the number of chunks that existed

sub chunk_total_count {
   my $self = shift;
   my $level = shift;
   $level = 'ALL' if(! defined($level));

   #$self->_load_chunks if(! $self->_level_started);
   $self->_load_extra if($self->{LEVEL}{EXTRA});

   my $count = $self->{LEVEL}{$level}{CHUNK_COUNT} || 0;

   return $count;
}
#-------------------------------------------------------------
#returns the id of the tier

sub id {
   my $self = shift @_;
   my $id = shift;

   if(defined($id)){
      $self->{TIER_ID} = $id;
   }

   return $self->{TIER_ID};
}
#--------------------------------------------------------------
#returns the parent id
#use this to identify parent of tiers in rtiers

sub parent {
   my $self = shift;
   my $arg = shift;

   if (defined($arg)) {
      $self->{PARENT} = $arg;
   }

   return $self->{PARENT};
}
#-------------------------------------------------------------
#returns the id of the tier (split on ':')

sub id_safe {
   my $self = shift @_;
   my ($id) = split(':', $self->{TIER_ID});

   return $id;
}

#-------------------------------------------------------------
#returns level the child tier belongs to in parent

sub level {
   my $self = shift @_;
   my ($id, $level) = split(':', $self->{TIER_ID});

   $level = $self->{LEVEL}{CURRENT} if(!$level);

   return $level;
}
#-------------------------------------------------------------
#returns tier_type value for the given tier

sub tier_type {
   my $self = shift @_;

   return $self->{TIER_TYPE};
}
#-------------------------------------------------------------
#returns the rank of the tier

sub rank {
   my $self = shift @_;
   my $arg = shift @_;

   if($arg){
       $self->{RANK} = $arg;
   }

   return $self->{RANK};
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
sub exception {
    my $self = shift;

    return $self->{EXCEPTION};
}
#-------------------------------------------------------------
#exception handler

sub _handler {
   my $self = shift;
   my $E = shift;
   my $extra = shift || '';

   my $level = $self->current_level;

   print STDERR $E->{-text};
   print STDERR "ERROR: ".$extra."\n" if($extra);

   #clear queue on error
   while(my $chunk = shift @{$self->{LEVEL}{ALL}{CHUNKS}}){
       my $clevel = $chunk->level;
       $self->{LEVEL}{$chunk->level}{RESULT_COUNT}++;
       $self->{LEVEL}{ALL}{RESULT_COUNT}++;
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
	confess "Error: Invalid method \'$call\' in Process::MpiTiers\n".
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
