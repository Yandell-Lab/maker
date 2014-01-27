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
use Carp;

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
         my $VARS           = $arg;
         $self->{LEVEL}     = shift @args;
         $self->{TIER_TYPE} = shift @args || 0;
         $self->{ID}        = shift @args || "0:".$self->{LEVEL}.":".$self->{TIER_TYPE}.":0";
         $self->{RANK}      = 0;
         $self->{FINISHED}  = 0;
         $self->{FAILED}    = 0;
         $self->{VARS}      = {};
         $self->{RESULTS}   = {};

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
   my $VARS       = shift;	#this should be a hash reference
   my $level      = $self->{LEVEL};
   my $tier_type  = $self->{TIER_TYPE};

   my @args = @{$self->_go('init', $VARS, $level, $tier_type)};

   foreach my $key (@args) {
      if (! exists $VARS->{$key}) {
	  confess "FATAL: argument \`$key\` does not exist in IPRtier object\n";
      }
      else {
	 $self->{VARS}{$key} = $VARS->{$key};
      }
   }
}
#--------------------------------------------------------------
#gets called by MpiTiers after result is processed

sub _finalize {
    my $self = shift;
    my $tier = shift;
    
    return if(! defined $tier->{VARS}->{extra_chunks});
    
    my $chunks = $tier->{VARS}->{extra_chunks};
    my $level = $tier->{LEVEL}{CURRENT};
    my $new_level = $tier->{VARS}->{extra_level};
    
    #initialize up to newest level
    $tier->go_to_level($new_level) if($new_level > $level);
    $tier->{LEVEL}{$new_level}{STARTED} = 1;
    $tier->{LEVEL}{EXTRA} = 1;
    
    #add chunks to appropriate level
    foreach my $c (@$chunks){
	push(@{$tier->{LEVEL}{$c->level}{CHUNKS}}, $c);
	$tier->{LEVEL}{$c->level}{CHUNK_COUNT}++;
    }

    #clear temporary variables
    delete($tier->{VARS}->{extra_chunks});
    delete($tier->{VARS}->{extra_level});
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
    my $q_def = $tier->{VARS}{q_def};
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

    if($tier->tier_type == 0){
	my $def = $tier->{VARS}->{q_def};

	$tier->{VARS} = {};
	$tier->{RESULTS} = {};

	$tier->{VARS}->{q_def} = $q_def;
	$DS_CTL->add_entry($seq_id, $out_dir, "FAILED");
    }
    else{
	$tier->{VARS} = {};
	$tier->{RESULTS} = {};
    }

   #restore logs
    $tier->{VARS}{LOG} = $LOG;
    $tier->{VARS}{DS_CTL} = $DS_CTL;
    
    return;
}

#--------------------------------------------------------------
#gets called by MpiTiers object following termination
sub _on_termination {
    my $self = shift;
    my $tier = shift;
    my $tier_type = $tier->{TIER_TYPE};

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

    #interrupted because another process is working on contig
    $tier->_set_interrupt(1) if($tier->{VARS}{c_flag} == -3);

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
   my $tier_type = shift || 0;

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
   my $VARS = shift;
   my $level = shift;
   my $tier_type = shift;
   my $tier  = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
       $tier = $tier_type;
       $tier_type = $level;
       $level = $VARS;
       $VARS = $self;
       $self = new Process::IPRchunk();
   }

   my $rank = $tier->rank;

   my $chunks = $self->_go('load', $VARS, $level, $tier_type);
   for (my $i = 0; $i < @{$chunks}; $i++){
       $chunks->[$i]->id("$rank:$level:$tier_type:$i");
       $chunks->[$i]->parent($tier->id);
   }

   return $chunks;
}
#--------------------------------------------------------------
#calls _go('run') to run code

sub run {
   my $self = shift;
   my $rank = shift;

   my $level = $self->{LEVEL};
   my $VARS = $self->{VARS};
   my $tier_type = $self->{TIER_TYPE};

   $self->{RANK} = shift || $self->{RANK};

   if ($self->{FINISHED} || $self->{FAILED}) {
      return undef;
   }

   my $results = $self->_go('run', $VARS, $level, $tier_type);

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
#for syntactic sugar, same as run

sub run_all {shift->run(@_);}
#--------------------------------------------------------------
#this funcion is called by MakerTiers.  It returns the flow of
#levels, i.e. order control, looping, etc.

sub _flow {
    my $self = shift;
    my $VARS = shift;
    my $level = shift;
    my $tier_type = shift;
    my $tier = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::IPRchunk" && ref($self) ne "Process::IPRchunk") {
       $tier = $tier_type;
       $tier_type = $level;
       $level = $VARS;
       $VARS = $self;
      $self = new Process::IPRchunk();
   }

   return $self->_go('flow', $VARS, $level, $tier_type);
}
#--------------------------------------------------------------
#initializes chunk variables, runs code, or returns flow
#depending on flag.
#This method is the core of the IPRchunk object

sub _go {
   my $self = shift;
   my $flag = shift;
   my $VARS = shift @_;
   my $level = shift @_;
   my $tier_type = shift;

   my $next_level = $level + 1;
   my $result_stat = 1;
   my @chunks;
   my @args;
   my %results;

   my $level_status = '';

   try{
      if ($level == 0) { #run log and initilaization
	 $level_status = 'initialization and checking run log';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($VARS, $level, $tier_type);
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
	    my $fasta   = Fasta::ucFastaRef($VARS->{fasta});
	    my $DS_CTL  =  $VARS->{DS_CTL};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	     
	    #get fasta parts
	    my $seq_id = Fasta::getSeqID($fasta);
	    my $seq_ref = Fasta::getSeqRef($fasta);
	    my $seq_length = length($$seq_ref);
	    my $safe_id = Fasta::seqID2SafeID($seq_id);
	    
	    #set up base and void directories for output
	    my ($out_dir, $the_void) = $DS_CTL->seq_dirs($seq_id);
	    
	    #contig combined output file
	    my $cfile = $out_dir."/".$safe_id.".ipr";
	    
	    my $LOG = iprscan::runlog->new(\%CTL_OPT,
					   {seq_id     => $seq_id,
					    seq_length => $seq_length,
					    out_dir    => $out_dir,
					    the_void   => $the_void,
					    fasta_ref  => $fasta},
					   "$out_dir/run.log"
					   );
	    
	    my $LOCK = $LOG->strip_off_lock();
	    
	    my ($c_flag, $message) = $LOG->get_continue_flag();

	    $DS_CTL->add_entry($seq_id, $out_dir, $message) if($message);
	     
	    #process existing iprscan output
	    if($c_flag == 0){
		die "ERROR: Can't find $cfile yet iprscan::runlog says the contig is finished\n"
		    if(! -e $cfile);

		#open for reading
		open(my $CFH, "< $cfile");		 
		my $data = '';
		$data .= join('', <$CFH>);
		close($CFH);

		#open for writing
		my $FH;
		my $outfile = $CTL_OPT{outfile};		 
		if($outfile){
		    open($FH, ">> $outfile");
		}
		else{
		    open($FH, ">&STDOUT");
		}		 

		#lock and write
		my $lock = new File::NFSLock("$CTL_OPT{out_base}/.output_lock", 'EX', 60, 60);
		print $FH $data;
		close($FH);
		$lock->unlock;
	    }
	    #-------------------------CODE
	     
	    #------------------------RETURN
	    %results = (seq_id => $seq_id,
			safe_id => $safe_id,
			seq_ref => $seq_ref,
			the_void => $the_void,
			cfile => $cfile,
			c_flag => $c_flag,
			out_dir => $out_dir,
			LOG => $LOG,
			LOCK => $LOCK);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
            #-------------------------RESULT
	     while (my $key = each %{$self->{RESULTS}}) {
		 $VARS->{$key} = $self->{RESULTS}->{$key};
	     }
            #-------------------------RESULT
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
	    my $chunk = new Process::IPRchunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{fasta			
		        seq_id
			seq_ref
		        safe_id
			the_void}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $fasta = $VARS->{fasta};
	    my $seq_id = $VARS->{seq_id};
	    my $seq_ref = $VARS->{seq_ref};
	    my $safe_id = $VARS->{safe_id};
	    my $the_void = $VARS->{the_void};

	    #safely escape %'s
	    $safe_id =~ s/x/\%78/g;
	    $safe_id =~ s/\%/x/g;

	    #fix fasta seq
	    $$seq_ref = uc($$seq_ref);
	    $$seq_ref =~ s/\*//g;
	    $$seq_ref =~ s/[^ABCDEFGHIKLMNPQRSTVWXYZ]/C/g;

	    #make a safe fasta
	    my $safe_fasta = Fasta::toFasta('>'.$safe_id, $seq_ref);
	    my $fasta_file = "$the_void/query.fasta";
	    FastaFile::writeFile($safe_fasta, $fasta_file);
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (safe_fasta => $safe_fasta,
			fasta_file => $fasta_file);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
            #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		$VARS->{$key} = $self->{RESULTS}->{$key};
	    }
            #-------------------------RESULT
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
	     foreach my $app (@{$VARS->{CTL_OPT}{_appl}}) {
		 $VARS->{app} = $app;
		 my $chunk = new Process::IPRchunk($VARS, $level, $tier_type);
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
		$id = uri_unescape($id);
		$result .= $id . $line;
	    }
	    close($IN);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ("$app" => $result
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
            #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		$VARS->{$key} = $self->{RESULTS}->{$key};
	    }
            #-------------------------RESULT
         }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 3) {	#collecting iprscan results
	 $level_status = 'collecting iprscan results';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::IPRchunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{CTL_OPT
			cfile
		       },		     
		     @{$VARS->{CTL_OPT}{_appl}}
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $cfile = $VARS->{cfile}; #combined contig output file
	    my $outfile = $CTL_OPT{outfile}; #combined every contig output file

	    #gather data
	    my $data = '';
	    foreach my $key (@{$CTL_OPT{_appl}}){
		next if(! defined ($VARS->{$key}));
		$data .= $VARS->{$key};
	    }

	    #write to the combined file
	    my $CFH;
	    open($CFH, "> $cfile");	    
	    print $CFH $data;
	    close($CFH);

	    #open global output file for witing
	    my $FH;
	    if($outfile){
		open($FH, ">> $outfile");
	    }
	    else{
		open($FH, ">&STDOUT");
	    }

	    #lock and write
	    my $lock = new File::NFSLock("$CTL_OPT{out_base}.output_lock", 'EX', 60, 60);
	    print $FH $data;
	    close($FH);
	    $lock->unlock;
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ();
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
            #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		$VARS->{$key} = $self->{RESULTS}->{$key};
	    }
            #-------------------------RESULT
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
   #return result_stat for result
   return $result_stat if($flag eq 'result');
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
   my $level = $self->{LEVEL};
   my $tier_type = shift;
   my $tier = shift;

   $level = $self->{LEVEL} if(!defined($level));
   $tier_type = $self->{TIER_TYPE} if(!defined($tier_type));

   #only return results for finished/succesful chunks
   return if (! $self->finished || $self->failed);

   #do level specific result processing
   return $self->_go('result', $VARS, $level, $tier_type);
}
#--------------------------------------------------------------
#sorts chunks for the tier into a more convenient order (optional)
sub _sort_levels{
    my $self = shift;
    my $chunks = shift;
    my $level = shift;
    my $tier_type = shift;
    my $tier = shift;
    
    #sorts by fasta chunk order then by level
    #so all levels on fasta chunk 1 get processed before chunk 2
    #but they can also be processed simultaneously
    @$chunks = sort {crit1($a) <=> crit1($b) || $a->level <=> $b->level || $a->id cmp $b->id} @$chunks;
    
    return $chunks;
}
#--------------------------------------------------------------
#first sort criteria
sub crit1 {
    my $chunk = shift;

    if(defined $chunk->{VARS}->{order}){
	return $chunk->{VARS}->{order};
    }
    else{
	return -1;
    }
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
#alias to finished (syntactic sugar)
sub terminated {return shift->finished;}
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
#returns the chunks number if part of a group

sub number {
    my $self = shift;
    my $arg = shift;

    my ($num) = $self->{ID} =~ /\:(\d+)$/;

    return $num;
}
#--------------------------------------------------------------
#returns the chunks rank
#i.e. on what MPI node the chunk ran

sub rank {
    my $self = shift;
    my $arg = shift;

    if (defined($arg)) {
	$self->{RANK} = $arg;
    }

    return $self->{RANK};
}
#--------------------------------------------------------------
#returns the parent id
#use this to identify correct tier to add to

sub parent {
    my $self = shift;
    my $arg = shift;

    if (defined($arg)) {
	$self->{PARENT} = $arg;
    }

    return $self->{PARENT};
}
#--------------------------------------------------------------
#return level the chunk is working on
sub level {
    my $self = shift;
    return $self->{LEVEL};
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
      $self->{RESULTS} = {};
      $self->{VARS} = {};
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
