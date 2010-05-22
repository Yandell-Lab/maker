#! /usr/bin/perl
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
use iprscan::runlog;

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
	 $self->{VARS}     = {};
	 $self->{RESULTS} = {};

	 $self->_initialize();
	 $self->_initialize_vars($VARS) if(ref($VARS) eq 'HASH');
      }
   }

   return $self;
}
#--------------------------------------------------------------
#Perform any user specified initialization steps here.
#This method always runs when a new chunk is created.
#This method does not take any arguements.
#$self->{VARS} has not been set yet, so if you need access
#to values in $self->{VARS}, add your code to the
#_initialize_vars method.

sub _initialize{
    my $self = shift;

    #do something here
}

#--------------------------------------------------------------
#Set up chunk specific variables.
#This method also reaches inside IPRtiers->{VARS} and pulls out
#all need variables.
#By Reaching into the IPRtiers object we keep control of the
#VARS datastructure within the IPRchunk (see method results).
#Another benefit is that variable order is no longer important,
#as it would be if variables were directly passed

sub _initialize_vars{
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
#gets called by MpiTiers object following a failure
sub _on_failure {
    my $self = shift;
    my $tier = shift;
    
    #handle case of calling as function rather than method
    if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
	$tier = $self;
	$self = new Process::IPRchunk();
    }
    
    my $level = $tier->current_level;
    my $seq_id = $tier->{VARS}{seq_id};
    my $out_dir = $tier->{VARS}{out_dir};
    my $LOG = $tier->{VARS}{LOG};
    my $DS_CTL = $tier->{VARS}{DS_CTL};
    
    print STDERR "FAILED CONTIG:$seq_id\n\n" if(defined $seq_id);
    
   if(defined $LOG){
       my $die_count = $LOG->get_die_count();
       $die_count++;
       $LOG->add_entry("DIED","RANK", $tier->id);
       $LOG->add_entry("DIED","COUNT",$die_count);
   }
    
    $DS_CTL->add_entry($seq_id, $out_dir, "FAILED");
    
    return;
}

#--------------------------------------------------------------
#gets called by MpiTiers object following termination
sub _on_termination {
    my $self = shift;
    my $tier = shift;

    #handle case of calling as function rather than method
    if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
	$tier = $self;
	$self = new Process::IPRchunk();
    }

    return if($tier->failed);
    return if($tier->{VARS}{c_flag} <= 0);

    #only reach this point if termination is due to success
    my $seq_id = $tier->{VARS}{seq_id};
    my $out_dir = $tier->{VARS}{out_dir};
    my $LOG = $tier->{VARS}{LOG};
    my $LOCK = $tier->{VARS}{LOCK};
    my $DS_CTL = $tier->{VARS}{DS_CTL};

    $DS_CTL->add_entry($seq_id, $out_dir, "FINISHED");
    $LOCK->unlock; #releases locks on the log file

    return;
}

#--------------------------------------------------------------
#gets called by MpiTiers object following termination
sub _should_continue {
    my $self = shift;
    my $tier = shift;

    #handle case of calling as function rather than method
    if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
	$tier = $self;
	$self = new Process::IPRchunk();
    }

    return $tier->{VARS}{c_flag};
}

#--------------------------------------------------------------
#this function/method is called by IPRtiers as part of IPRtiers
#initialization to prepare incoming data before building chunks.
#This method should not be called directly by the user or inside
#IPRchunks. Putting this preparation here as opposed to inside
#IPRtiers makes IPRchunks more portable and makes debugging
#easier. It will always run on the root node before distributing
#the tier.

sub _prepare {
   my $self = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
      $VARS = $self;
      $self = new Process::IPRchunk();
   }

   #instantiate empty LOG   
   $VARS->{LOG} = undef;

   #==Prepare data here as part of initialization

   #set up contig variables
   #===
   my $status = 'instantiating tier variables';
   #===

   #-set up variables that are the result of chunk accumulation
   #$VAR->{} = [];

   #--other variables
   $VARS->{c_flag} = 1; #always continue with this implementation and
                        #let child nodes decide on running chunk.
                        #Other implementations let the root node decide
   
   return 1;
}
#--------------------------------------------------------------
#This function is called by IPRtiers.  It returns all the
#chunks for a given level.  This allows for the number of
#chunks created per level to be controlled within IPRchunks.

sub _loader {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;
   my $tID  = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
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

sub _flow {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
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
      if ($level == 0) { #run log and initilaization
	 $level_status = 'initialization and checking run log';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{fasta
			CTL_OPT
			DS_CTL}
		     );
            #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $fasta   = Fasta::ucFasta(\$VARS->{fasta});
	    my $DS_CTL  =  $VARS->{DS_CTL};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	     
	    #get fasta parts
	    my $seq_id = Fasta::getSeqID(\$fasta);
	    my $safe_id = Fasta::seqID2SafeID($seq_id);
	    
	    #set up base and void directories for output
	    my ($out_dir, $the_void) = $DS_CTL->seq_dirs($seq_id);
	    
	    #contig combined output file
	    my $cfile = $out_dir."/".$safe_id.".ipr";
	    
	    my $LOG = iprscan::runlog->new(\%CTL_OPT,
					   {seq_id     => $seq_id,
					    out_dir    => $out_dir,
					    the_void   => $the_void,
					    fasta_ref  => \$fasta},
					   "$out_dir/run.log"
					   );
	    
	    my $LOCK = $LOG->strip_off_lock();
	     
	    my ($c_flag, $message) = $LOG->get_continue_flag();
	    $DS_CTL->add_entry($seq_id, $out_dir, $message) if($message);
	     
	    #process existing iprscan output
	    if($c_flag == 0){		   
		die "ERROR: Can't find $cfile yet iprscan::runlog says the contig is finished\n"
		    if(! -e $cfile);
		 
		my $lock = new File::NFSLock(".iprscan_lock", 'EX', 60, 60);
		my $outfile = $CTL_OPT{outfile};
		 
		my $FH;
		if($outfile){
		    open($FH, ">> $outfile");
		}
		else{
		    open($FH, ">&STDOUT");
		}
		 
		#open for reading
		open(my $CFH, "< $cfile");
		 
		while(my $line = <$CFH>){
		    print $FH $line;
		}
		 
		close($FH);
		close($CFH);
		 
		$lock->unlock;
	    }
	    #-------------------------CODE
	     
	    #------------------------RESULTS
	    %results = (seq_id => $seq_id,
			safe_id => $safe_id,
			the_void => $the_void,
			cfile => $cfile,
			c_flag => $c_flag,
			out_dir => $out_dir,
			LOG => $LOG,
			LOCK => $LOCK);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 1) {	#fixing fasta
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
		        safe_id
			the_void}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $fasta = $VARS->{fasta};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_id = $VARS->{safe_id};
	    my $the_void = $VARS->{the_void};

	    #safely escape %'s
	    $safe_id =~ s/x/\%78/g;
	    $safe_id =~ s/\%/x/g;

	    #get seq
	    my $seq = Fasta::getSeqRef(\$fasta);
	    
	    #fix fasta seq
	    $$seq = uc($$seq);
	    $$seq =~ s/\*//g;
	    $$seq =~ s/[^ABCDEFGHIKLMNPQRSTVWXYZ]/C/g;

	    #make a safe fasta
	    my $safe_fasta = Fasta::toFasta('>'.$safe_id, $seq);
	    my $fasta_file = "$the_void/query.fasta";
	    FastaFile::writeFile($safe_fasta, $fasta_file);
	    #-------------------------CODE
	 
	    #------------------------RESULTS
	    %results = (safe_fasta => $safe_fasta,
			fasta_file => $fasta_file);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 2) {	#running iprscan
	 $level_status = 'running iprscan';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	     foreach my $app (@{$VARS->{CTL_OPT}{appl}}) {
		 $VARS->{app} = $app;
		 my $chunk = new Process::IPRchunk($level, $VARS);
		 push(@chunks, $chunk);
	     }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{CTL_OPT
			safe_id
			fasta_file
			app
			params
			iprscan
			the_void
			LOG
			}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT    = %{$VARS->{CTL_OPT}};
	    my $fasta_file = $VARS->{fasta_file};
	    my $safe_id    = $VARS->{safe_id};
	    my $app        = $VARS->{app};
	    my $params     = $VARS->{params};
	    my $iprscan    = $VARS->{iprscan};
	    my $the_void   = $VARS->{the_void};
	    my $LOG        = $VARS->{LOG};
	    
	    #specify output files
	    my $ofile = "$the_void/$safe_id.$app";

	    #build command
	    my $command = $iprscan;
	    $command .= " -nocrc" if($CTL_OPT{nocrc});
	    $command .= " -seqtype $CTL_OPT{seqtype}" if(defined $CTL_OPT{seqtype});
	    $command .= " -trtable $CTL_OPT{trtable}" if(defined $CTL_OPT{trtable});
	    $command .= " -goterms" if($CTL_OPT{goterms});
	    $command .= " -iprlookup" if($CTL_OPT{iprlookup});
	    $command .= " -format $CTL_OPT{format}" if(defined $CTL_OPT{format});
	    $command .= " -verbose" if($CTL_OPT{verbose});
	    $command .= " " . join(' ', @$params);
	    $command .= " -appl $app -i $fasta_file -o $ofile";

	    $LOG->add_entry("STARTED", $ofile, "");

	    my $w = new Widget::iprscan();
	    if(-e $ofile){
		print STDERR "rereading $ofile...\n" unless($main::quiet);
	    }
	    else{
		print STDERR "running iprscan $app...\n" unless($main::quiet);
		$w->run($command);
	    }

	    $LOG->add_entry("FINISHED", $ofile, "");

	    open(my $IN, "< $ofile");
	    my $result;
	    while(my $line = <$IN>){
		#extra steps because InterProScan can't handle %
		$line =~ s/^([^\s\t]+)//;
		my $id = $1;
		$id =~ s/x/\%/g;
		$id =~ uri_unescape($id);
		$result .= $id . $line;
	    }
	    close($IN);
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
      elsif ($level == 3) {	#blastx repeat mask
	 $level_status = 'collecting iprscan results';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{CTL_OPT
			cfile
		       },		     
		     @{$VARS->{CTL_OPT}{appl}}
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $cfile = $VARS->{cfile}; #combined contig output file
	    my $outfile = $CTL_OPT{outfile}; #combined every contig output file

	    my $lock = new File::NFSLock(".iprscan_lock", 'EX', 40, 40);

	    my $FH;
	    if($outfile){
		open($FH, ">> $outfile");
	    }
	    else{
		open($FH, ">&STDOUT");
	    }

	    my $CFH;
	    open($CFH, "> $cfile");

	    foreach my $key (@{$CTL_OPT{appl}}){
		next if(! defined ($VARS->{$key}));
		print $FH $VARS->{$key};
		print $CFH $VARS->{$key};
	    }

	    close($FH);
	    close($CFH);

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

sub _result {
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
