#! /usr/bin/perl -w

package Process::IPRchunk;

use strict;

use Error qw(:try);
use Error::Simple;
use Storable;
use File::NFSLock;
use Fasta;
use FastaFile;
use File::Temp qw(tempfile);
use Widget::iprscan;
use URI::Escape;

#--set object variables for serialization of data
#this is needed when cloning an IPRchunk object
$Storable::forgive_me = 1; #allows serializaion of objects with code refs

#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------
sub new {
   my ($class, @args) = @_;

   my $self = {};

   bless ($self, $class);

   if (@args) {
      my $arg = shift @args;
      if (ref $arg eq 'Process::IPRchunk') {
	 $self = $arg->clone();
      }
      else {
	 $self->{LEVEL} = $arg;
	 my $VARS       = shift @args;
	 $self->{ID}    = shift @args || "0:".$self->{LEVEL}.":0";
	 $self->{RANK}  = shift @args || 0;
	 $self->{FINISHED} = 0;
	 $self->{FAILED}   = 0;
	 $self->{VARS} = {};
	 $self->{RESULTS} = {};

	 $self->_initialize($VARS);
      }
   }

   return $self;
}
#--------------------------------------------------------------
#this function/method is called by IPRtiers as part of IPRtiers
#initialization to prepare incoming data before building chunks.
#This method should not be called directly by the user or inside
#IPRchunks. Putting this preparation here as opposed to inside
#IPRtiers makes IPRchunks more portable and makes debugging
#easier.

sub prepare {
   my $self = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if (ref($self) ne "Process::IPRchunk") {
      $VARS = $self;
      $self = new Process::IPRchunk();
   }
   
   #==Prepare data here as part of initialization

   my $status = ''; #for error control
   #process fasta file
   try{
      #===
      $status = 'preparing iprscan job';
      #===
      
      my ($outname) = $VARS->{infile} =~ /([^\/]+)$/;
      $outname = $VARS->{outfile} if($VARS->{outfile});

      $VARS->{seq_id} = Fasta::getSeqID(\$VARS->{fasta});
      $VARS->{safe_id} = Fasta::seqID2SafeID($VARS->{seq_id});
      $VARS->{c_flag} = $VARS->{DS_CTL}->continue_flag($VARS->{seq_id});

      my $failed .= "$outname.failed/".$VARS->{safe_id}.".fasta";

      if(! $VARS->{c_flag}){
	  warn "WARNING: ".$VARS->{seq_id} ."failed 2 time and will not be tried again\n",
	       "The fasta sequence will be saved for debugging in $failed\n";

	  mkdir("$outname.failed") if(! -d "$outname.failed");

	  open(my $OUT, ">$failed");
	  print $OUT $VARS->{fasta};
	  close($OUT);
	  $VARS->{DS_CTL}->add_entry($VARS->{seq_id}, 'DIED_PERMANENT');
      }
      else{
	  unlink($failed) if(-e $failed); #remove old failed
      }

      $VARS->{DS_CTL}->add_entry($VARS->{seq_id}, 'STARTED');
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, $status, 'throw');
   };
   
   return 1;
}
#--------------------------------------------------------------
#This function is called by IPRtiers.  It returns all the
#chunks for a given level.  This allows for the number of
#chunks created per level to be controlled within IPRchunks.

sub loader {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;
   my $tID  = shift;

   #handle case of calling as function rather than method
   if (ref($self) ne "Process::IPRchunk") {
      $tID = $VARS;
      $VARS = $level;
      $level = $self;
      $self = new Process::IPRchunk();
   }

   my $chunks = $self->_go('load', $level, $VARS);
   for (my $i = 0; $i < @{$chunks}; $i++){
      $chunks->[$i]->id("$tID:$level:$i");
   }

   return $chunks;
}
#--------------------------------------------------------------
#Set up chunk specific variables.
#This method also reaches inside IPRtiers->{VARS} and pulls out
#all need variables.
#By Reaching into the IPRtiers object we keep control of the
#VARS datastructure within the IPRchunk (see method results).
#Another benefit is that variable order is no longer important,
#as it would be if variables were directly passed

sub _initialize{
   my $self       = shift;
   my $level      = $self->{LEVEL};
   my $VARS       = shift;	#this should be a hash reference

   my @args = @{$self->_go('init', $level, $VARS)};

   foreach my $key (@args) {
      if (! exists $VARS->{$key}) {
	 die "FATAL: argument \`$key\` does not exist in IPRtier object\n";
      }
      else {
	 $self->{VARS}{$key} = $VARS->{$key};
      }
   }
}
#--------------------------------------------------------------
#cals _go('run') to run code

sub run {
   my $self = shift;
   my $level = $self->{LEVEL};
   my $VARS = $self->{VARS};

   $self->{RANK} = shift || $self->{RANK};

   if ($self->{FINISHED} || $self->{FAILED}) {
      return undef;
   }

   my $results = $self->_go('run', $level, $VARS);

   $self->{VARS} = {};
   $self->{RESULTS} = $results;
   $self->{FINISHED} = 1;

   if(! $self->failed){
      return 1 ;
   }
   else{
      return undef;
   }
}
#--------------------------------------------------------------
#this funcion is called by MakerTiers.  It returns the flow of
#levels, i.e. order control, looping, etc.

sub flow {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if (ref($self) ne "Process::IPRchunk") {
      $VARS = $level;
      $level = $self;
      $self = new Process::IPRchunk();
   }

   return $self->_go('flow', $level, $VARS);
}
#--------------------------------------------------------------
#initializes chunk variables, runs code, or returns flow
#depending on flag.
#This method is the core of the IPRchunk object

sub _go {
   my $self = shift;
   my $flag = shift;
   my $level = shift @_;
   my $VARS = shift @_;

   my $next_level = $level + 1;
   my @chunks;
   my @args;
   my %results;

   my $level_status = '';

   try{
      if ($level == 0) {	#fixing fasta
	 $level_status = 'fixing fasta';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{fasta			
		        seq_id
		        safe_id}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $fasta = $VARS->{fasta};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_id = $VARS->{safe_id};
	    
	    #get seq
	    my $seq = Fasta::getSeqRef(\$fasta);
	    
	    #fix fasta seq
	    $$seq = uc($$seq);
	    $$seq =~ s/\*//g;
	    $$seq =~ s/[^ABCDEFGHIKLMNPQRSTVWXYZ]/C/g;

	    #make a safe fasta
	    my $safe_fasta = Fasta::toFasta('>'.$safe_id, $seq);
	    #-------------------------CODE
	 
	    #------------------------RESULTS
	    %results = (safe_fasta => $safe_fasta);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 1) {	#running iprscan
	 $level_status = 'running iprscan';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	     foreach my $app (@{$VARS->{appl}}) {
		 $VARS->{app} = $app;
		 my $chunk = new Process::IPRchunk($level, $VARS);
		 push(@chunks, $chunk);
	     }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{safe_fasta
			app
			params
			iprscan
			}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $safe_fasta = $VARS->{safe_fasta};
	    my $app        = $VARS->{app};
	    my $params     = $VARS->{params};
	    my $iprscan    = $VARS->{iprscan};

	    #make files
	    my (undef, $ifile) = tempfile();
	    my (undef, $ofile) = tempfile();
	    FastaFile::writeFile($safe_fasta, $ifile);

	    #build command
	    my $command = $iprscan . " " . join(' ', @$params);
	    $command .= " -appl $app -i $ifile -o $ofile";

	    my $w = new Widget::iprscan();
	    print STDERR "running iprscan $app...\n" unless($main::quiet);
	    $w->run($command);

	    unlink($ifile);
	    open(my $IN, "< $ofile");
	    my $result;
	    while(my $line = <$IN>){
		$result .= uri_unescape($line);
	    }
	    close($IN);
	    unlink($ofile);
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = ("$app" => $result
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 2) {	#blastx repeat mask
	 $level_status = 'collecting iprscan results';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{outfile},
		     @{$VARS->{appl}}
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $outfile = $VARS->{outfile};

	    my $lock = new File::NFSLock(".iprscan_lock", 'EX', 40, 40);

	    my $FH;
	    if($outfile){
		open($FH, ">> $outfile");
	    }
	    else{
		open($FH, ">&STDOUT");
	    }

	    while(my $key = each %{$VARS}){
		next if($key eq 'outfile');
		print $FH $VARS->{$key};
	    }
	    close($FH);

	    $lock->unlock;
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = ();
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    $next_level = undef;
	    #-------------------------NEXT_LEVEL
	 }
      }
      else {
	 warn "WARNING: Invalid level for method _go() in Process::IPRchunk\n";
	 return undef;
      }
   }
   catch Error::Simple with{
      my $E = shift;

      my $tag = ($flag eq 'run') ? 'handle' : 'throw';

      $self->_handler($E, $level_status, $tag);
   };

   #return args list for initializing
   return \@args if($flag eq 'init');
   #return results after running
   return \%results if($flag eq 'run');
   #return chunks for loader
   return \@chunks if($flag eq 'load');
   #return next_level for flow
   return $next_level if($flag eq 'flow');

   #should never reach this line
   die "FATAL: \'$flag\' is not a valid flag in IPRchunk _go!!\n";
}

#--------------------------------------------------------------
#reaches inside IPRtiers->{VARS} and places result in the 
#objects internal structure.
#By having this method work on IPRtiers->{VARS} rather than 
#simply returning a value for the results, we keep all contol
#of data structure inside of IPRchunk which makes for easier
#debugging

sub result {
   my $self = shift;
   my $VARS = shift;		#this ia a hash reference;
   my $RESULTS = $self->{RESULTS};	#this ia a hash reference

   #only return results for finished/succesful chunks
   return if (! $self->finished || $self->failed);

   while (my $key = each %{$RESULTS}) {
      $VARS->{$key} = $RESULTS->{$key};
   }

   return 1;
}
#--------------------------------------------------------------
#returns true if failed

sub failed {
   my $self = shift;
    
   return $self->{FAILED} || undef;
}
#--------------------------------------------------------------
#returns true if finished

sub finished {
   my $self = shift;

   return $self->{FINISHED} || undef;
}
#--------------------------------------------------------------
#returns Error::Simple exception object
#this is created after a failure

sub exception {
   my $self = shift;
    
   return $self->{EXCEPTION} || undef;
}
#--------------------------------------------------------------
#this is a good place to return STDERR or status messages
#for building log files
#you must fill $self->{ERROR} yourself

sub error {
   my $self = shift;
    
   return $self->{ERROR} || '';
}

#--------------------------------------------------------------
#returns the chunks id
#use this to keep track of chunks
#i.e. what tier they are from, what level they are,
#and what node they are running on

sub id {
   my $self = shift;
   my $arg = shift;

   if ($arg) {
      $self->{ID} = $arg;
   }

   return $self->{ID};
}

#--------------------------------------------------------------
#deep clone of self
#this uses the storable module, so globs and code refs are ignored

sub clone {
   my $self = shift;
    
   my $clone = Storable::dclone($self);

   return $clone;
}

#--------------------------------------------------------------
#exception handler

sub _handler {
   my $self = shift;
   my $E = shift;
   my $extra = shift;
   my $tag = shift || 'throw';

   $E->{-text} .= "ERROR: Failed while ".$extra."!!\n\n" if($extra);
   
   if($tag eq 'handle'){
      $self->{FAILED} = 1;
      $self->{EXCEPTION} = $E;
   }
   else{
      throw Error::Simple($E);
   }
}
#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
1;