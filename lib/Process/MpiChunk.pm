#! /usr/bin/perl -w

package Process::MpiChunk;

use strict;

use Error qw(:try);
use Error::Simple;
use Storable;

#--set object variables for serialization of data
#this is needed when cloning an MpiChunk object
$Storable::forgive_me = 1; #allows serializaion of objects with code refs

#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------
sub new { #do not change or edit this
   my ($class, @args) = @_;

   my $self = {};

   bless ($self, $class);

   if (@args) {
      my $arg = shift @args;
      if (ref $arg eq 'Process::MpiChunk') {
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
	 $self->{RESULTS}  = {};

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
#Set up chunk specific variables. Only runs when $VARS exist
#This method also reaches inside MpiTiers->{VARS} and pulls out
#all need variables.
#By Reaching into the MpiTiers object we keep control of the
#VARS datastructure within the MpiChunk (see method results).
#Another benefit is that variable order is no longer important,
#as it would be if variables were directly passed

sub _initialize_vars{
   my $self       = shift;
   my $level      = $self->{LEVEL};
   my $VARS       = shift;	#this should be a hash reference

   my @args = @{$self->_go('init', $level, $VARS)};

   foreach my $key (@args) {
      if (! exists $VARS->{$key}) {
	 die "FATAL: argument \`$key\` does not exist in MpiTier object\n";
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
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $self;
      $self = new Process::MpiChunk();
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
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $self;
      $self = new Process::MpiChunk();
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
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $self;
      $self = new Process::MpiChunk();
   }

   #interrupted because another process is working on contig
   $tier->_set_interrupt(1) if($tier->{VARS}{c_flag} == -3);

   return $tier->{VARS}{c_flag};
}
#--------------------------------------------------------------
#this function/method is called by MpiTiers as part of MpiTiers
#initialization to prepare incoming data before building chunks.
#This method should not be called directly by the user or inside
#MpiChunks. Putting this preparation here as opposed to inside
#MpiTiers makes MpiChunks more portable and makes debugging
#easier.  It will always run on the root node before distributing
#the tier.

sub _prepare {
   my $self = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $VARS = $self;
      $self = new Process::MpiChunk();
   }

   #instantiate empty LOG
   $VARS->{LOG} = undef;
   
   #==Prepare data here as part of initialization

   #set up contig variables
   #===
   my $status = 'instantiating tier variables';
   #===
   
   #-set up variables that are heldover from last chunk
   $VARS->{holdover_blastn}     = [];
   $VARS->{holdover_blastx}     = [];
   $VARS->{holdover_tblastx}    = [];
   $VARS->{holdover_pred}       = [];
   $VARS->{holdover_est_gff}    = [];
   $VARS->{holdover_altest_gff} = [];
   $VARS->{holdover_prot_gff}   = [];
   $VARS->{holdover_pred_gff}   = [];
   $VARS->{holdover_model_gff}  = [];
   
   #-set up variables that are the result of chunk accumulation
   my $empty;
   $VARS->{masked_total_seq} = \$empty; #empty scalar ref
   $VARS->{p_fastas} = {};
   $VARS->{t_fastas} = {};
   
   #--other variables
   $VARS->{res_dir} = undef;
   $VARS->{c_flag} = 1; #always continue with this implementation and
                        #let child nodes decide on running chunk.
                        #Other implementations let the root node decide

   return 1;
}
#--------------------------------------------------------------
#This function is called by MpiTiers.  It returns all the
#chunks for a given level.  This allows for the number of
#chunks created per level to be controlled within MpiChunks.

sub _loader {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;
   my $tID  = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tID = $VARS;
      $VARS = $level;
      $level = $self;
      $self = new Process::MpiChunk();
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
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $VARS = $level;
      $level = $self;
      $self = new Process::MpiChunk();
   }

   return $self->_go('flow', $level, $VARS);
}
#--------------------------------------------------------------
#initializes chunk variables, runs code, or returns flow
#depending on flag.
#This method is the core of the MpiChunk object

sub _go {
   my $self = shift;
   my $flag = shift;
   my $level = shift @_;
   my $VARS = shift @_;

   my $next_level = $level + 1;
   my $result_stat = 1;
   my @chunks;
   my @args;
   my %results;

   my $level_status = '';

   try{
      if ($level == 0) {	#examining contents of the fasta file and run log
	 $level_status = 'examining contents of the fasta file and run log';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
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
	    my $fasta = Fasta::ucFastaRef($VARS->{fasta}); #build uppercase fasta
	    my $DS_CTL = $VARS->{DS_CTL};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};


	    #get fasta parts
	    my $q_def = Fasta::getDef($fasta); #Get fasta header
	    my $q_seq_ref = Fasta::getSeqRef($fasta); #Get reference to fasta sequence
	    my $seq_id = Fasta::def2SeqID($q_def); #Get sequence identifier
	    my $safe_seq_id = Fasta::seqID2SafeID($seq_id); #Get safe named identifier	    

	    #set up base and void directories for output
	    my ($out_dir, $the_void) = $DS_CTL->seq_dirs($seq_id);
	    
	    #-build and proccess the run log
	    my $LOG = runlog->new(\%CTL_OPT,
				  {seq_id     => $seq_id,
				   seq_length => length(${$q_seq_ref}),
				   out_dir    => $out_dir,
				   the_void   => $the_void,
				   fasta_ref  => $fasta},
				  "$out_dir/run.log"
				  );

	    my $LOCK = $LOG->strip_off_lock();
	    
	    my ($c_flag, $message) = $LOG->get_continue_flag();
	    $DS_CTL->add_entry($seq_id, $out_dir, $message) if($message);
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (fasta => $fasta,
			q_def => $q_def,
			q_seq_ref => $q_seq_ref,
			seq_id => $seq_id,
			safe_seq_id => $safe_seq_id,
			out_dir => $out_dir,
			the_void => $the_void,
			c_flag => $c_flag,
			LOG => $LOG,
			LOCK => $LOCK
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
     elsif ($level == 1) {	#set up GFF3 output and fasta chunks
	 $level_status = 'setting up GFF3 output and fasta chunks';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{fasta
			CTL_OPT
			out_dir
			build
			seq_id
			safe_seq_id
			the_void
			q_seq_ref}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $out_dir = $VARS->{out_dir};
	    my $build = $VARS->{build};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $fasta = $VARS->{fasta};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    
	    
	    my $GFF3 = Dumper::GFF::GFFV3->new("$out_dir/$safe_seq_id.gff",
					       $build,
					       $the_void
					      );

	    $GFF3->set_current_contig($seq_id, $q_seq_ref);     
	    
	    my $fasta_chunker = new FastaChunker();
	    $fasta_chunker->parent_fasta($fasta);
	    $fasta_chunker->chunk_size($CTL_OPT{max_dna_len});
	    $fasta_chunker->min_size($CTL_OPT{split_hit});
	    $fasta_chunker->load_chunks();

	    #--build an index of the databases
	    my $proteins = $VARS->{CTL_OPT}{_p_db};
	    my $trans    = $VARS->{CTL_OPT}{_e_db};
	    my $altests  = $VARS->{CTL_OPT}{_a_db};
	    my $fasta_t_index = GI::build_fasta_index($trans) if($trans); 
	    my $fasta_p_index = GI::build_fasta_index($proteins) if($proteins);
	    my $fasta_a_index = GI::build_fasta_index($altests) if($altests); 	 

	    my $chunk = $fasta_chunker->next_chunk();
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (GFF3 => $GFF3,
			fasta_chunker => $fasta_chunker,
			chunk => $chunk,
			fasta_t_index => $fasta_t_index,
			fasta_p_index => $fasta_p_index,
			fasta_a_index => $fasta_a_index
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
	     #if no masking skip those steps
	     if(!$VARS->{CTL_OPT}{model_org} && !$VARS->{CTL_OPT}{rmlib} &&
		!$VARS->{CTL_OPT}{rm_gff} && !$VARS->{CTL_OPT}{repeat_protein}
	       ){
		 $VARS->{masked_total_seq} = $VARS->{q_seq_ref};
		 $next_level = 6;
	     }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 2) {	#do repeat masking
	 $level_status = 'doing repeat masking';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			the_void
			safe_seq_id
			q_seq_ref
			GFF_DB
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $LOG = $VARS->{LOG};
	    my $GFF_DB = $VARS->{GFF_DB};
	    my $chunk = $VARS->{chunk};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};

	    #-- repeatmask with gff3 input
	    my $rm_gff_keepers = [];
	    if ($CTL_OPT{go_gffdb}) {
	       $rm_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							    $q_seq_ref,
							    'repeat'
							   );
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_gff_keepers) if($CTL_OPT{organism_type} ne 'prokaryotic');
	    }
      
	    #-- repeatmask with RepeatMasker	 
	    my $rm_rb_keepers = []; #repeat masker RepBase
	    if ($CTL_OPT{model_org}) { #model organism repeats
	       $rm_rb_keepers = GI::repeatmask($chunk,
					       $the_void,
					       $safe_seq_id,
					       $CTL_OPT{model_org},
					       $CTL_OPT{RepeatMasker},
					       '',
					       $CTL_OPT{cpus},
					       $LOG
					      );
	 
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_rb_keepers) if($CTL_OPT{organism_type} ne 'prokaryotic');
	    }
	    my $rm_sp_keepers = []; #repeat masker species
	    if ($CTL_OPT{rmlib}) {  #species specific repeats;
	       $rm_sp_keepers = GI::repeatmask($chunk,
					       $the_void,
					       $safe_seq_id,
					       $CTL_OPT{model_org},
					       $CTL_OPT{RepeatMasker},
					       $CTL_OPT{rmlib},
					       $CTL_OPT{cpus},
					       $LOG
					      );
	 
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_sp_keepers) if($CTL_OPT{organism_type} ne 'prokaryotic');
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (rm_gff_keepers => $rm_gff_keepers,
			rm_rb_keepers => $rm_rb_keepers,
			rm_sp_keepers => $rm_sp_keepers,
			chunk => $chunk
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
      elsif ($level == 3) {	#blastx repeat mask
	 $level_status = 'doing blastx repeats';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{rm_blastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_r_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk},
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'repeatrunner');
		
		next if($fin{$blast_finished});

		$db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
		$VARS->{db} = $db;
		$VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
		$fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
		my $chunk = new Process::MpiChunk($level, $VARS);
		push(@chunks, $chunk);
	    }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{db
			chunk
			the_void
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    my $res_dir;
	    my $rm_blastx_keepers = [];
	    if ($db) {
	       GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
	       ($rm_blastx_keepers, $res_dir) = GI::repeatrunner_as_chunks($chunk,
									   $db,
									   $the_void,
									   $safe_seq_id,
									   \%CTL_OPT,
									   $self->{RANK},
									   $LOG,
									   $LOG_FLAG
									   );
		 
	       $rm_blastx_keepers = GI::clean_blast_hits($rm_blastx_keepers,
							 $CTL_OPT{pcov_rm_blastx},
							 $CTL_OPT{pid_rm_blastx},
							 $CTL_OPT{eval_rm_blastx},
							 0 #contiguity flag
							 );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( rm_blastx_keepers => $rm_blastx_keepers,
			 res_dir => $res_dir
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		push(@{$VARS->{$key}}, $self->{RESULTS}->{$key});
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 4) {	#collect blastx repeatmask
	 $level_status = 'collecting blastx repeatmasking';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			res_dir
		        rm_blastx_keepers
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $res_dir = $VARS->{res_dir};
	    my $rm_blastx_keepers = $VARS->{rm_blastx_keepers};
	    my $LOG = $VARS->{LOG};

	    $rm_blastx_keepers = GI::collect_blast($chunk,
						   $rm_blastx_keepers,
						   $res_dir,
						   $LOG
						   );
	    
	    #mask the chunk
	    $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_blastx_keepers) if($CTL_OPT{organism_type} ne 'prokaryotic');

	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( chunk => $chunk,
			 res_dir => $res_dir,
			 rm_blastx_keepers => $rm_blastx_keepers
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
      elsif ($level == 5) {	#process all repeats
	 $level_status = 'processing all repeats';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			rm_gff_keepers
			rm_rb_keepers
			rm_sp_keepers
			rm_blastx_keepers
			q_seq_ref
			GFF3
			masked_total_seq}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $chunk = $VARS->{chunk};
	    my $rm_gff_keepers = $VARS->{rm_gff_keepers};
	    my $rm_rb_keepers = $VARS->{rm_rb_keepers};
	    my $rm_sp_keepers = $VARS->{rm_sp_keepers};
	    my $rm_blastx_keepers = $VARS->{rm_blastx_keepers};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $GFF3 = $VARS->{GFF3};
	    my $masked_total_seq = $VARS->{masked_total_seq};


	    #-combine and cluster repeat hits for consensus
	    my $rm_keepers = repeat_mask_seq::process($rm_gff_keepers, 
						      $rm_rb_keepers, 
						      $rm_sp_keepers, 
						      $rm_blastx_keepers,
						      $q_seq_ref
						     );
      
	    #-add repeats to GFF3
	    $GFF3->add_repeat_hits($rm_keepers);
      
	    #-build big masked sequence
	    $$masked_total_seq .= $chunk->seq();
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (rm_keepers => $rm_keepers,
			masked_total_seq => $masked_total_seq,
			rm_gff_keepers => [], #clear memory
			rm_rb_keepers => [],  #clear memory
			rm_sp_keepers => [],  #clear memory
			rm_blastx_keepers => [] #clear memory
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
	    if ($VARS->{chunk} = $VARS->{fasta_chunker}->next_chunk) {
	       $next_level = 2;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 6) {	#prep masked sequence and abinits
	 $level_status = 'preparing masked sequence and ab-inits';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{the_void
			safe_seq_id
			q_def
			masked_total_seq
			fasta
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $q_def = $VARS->{q_def};
	    my $masked_total_seq = $VARS->{masked_total_seq};
	    my $fasta = $VARS->{fasta};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};	 
	 
	    my $masked_fasta = Fasta::toFastaRef($q_def.' masked', $masked_total_seq);
	    my $masked_file = $the_void."/query.masked.fasta";
	    FastaFile::writeFile($masked_fasta, $masked_file) unless($CTL_OPT{_no_mask});
	 
	    my $unmasked_file = $the_void."/query.fasta";
	    FastaFile::writeFile($fasta, $unmasked_file);

	    #==ab initio predictions here
	    #do masked predictions first
	    my $preds = [];
	    if(! $CTL_OPT{_no_mask}){ 
		$preds = GI::abinits($masked_file,
				     $the_void,
				     $safe_seq_id.".masked",
				     \%CTL_OPT,
				     $LOG,
				     $unmasked_file #for genemark
				     );
	    
		#add tag that says these are masked
		foreach my $p (@$preds){
		    my $alg = $p->algorithm();
		    next if ($alg eq 'genemark'); #genemark is never masked because it splits genes
		    $p->algorithm("$alg\_masked");
		    foreach my $h($p->hsps){
			$h->algorithm("$alg\_masked");
		    }
		}
	    }

	    #now do unmasked predictions
	    if($CTL_OPT{unmask} || $CTL_OPT{_no_mask}){
	       my $unmasked_preds = GI::abinits($unmasked_file,
						$the_void,
						$safe_seq_id,
						\%CTL_OPT,
						$LOG
						);

	       push(@$preds, @$unmasked_preds);
	    }

	    #==QRNA noncoding RNA prediction here
	    my $qra_preds = [];
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (masked_fasta => $masked_fasta,
			preds => $preds,
			qra_preds => $qra_preds
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
      elsif ($level == 7) {	#prep new fasta chunks
	 $level_status = 'preparing new fasta chunks';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{masked_fasta
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $masked_fasta = $VARS->{masked_fasta};


	    #--reset fastachunker for masked chunks
	    my $fasta_chunker = new FastaChunker();
	    $fasta_chunker = new FastaChunker();
	    $fasta_chunker->parent_fasta($masked_fasta);
	    $fasta_chunker->chunk_size($CTL_OPT{max_dna_len});
	    $fasta_chunker->min_size($CTL_OPT{split_hit});
	    $fasta_chunker->load_chunks();

	    my $chunk = $fasta_chunker->next_chunk();
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (fasta_chunker => $fasta_chunker,
			masked_fasta => $masked_fasta,
			chunk => $chunk
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
      elsif ($level == 8) {	#blastn
	 $level_status = 'doing blastn of ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{blastn_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_e_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk},
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'blastn');

		next if($fin{$blast_finished});

		$db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
		$VARS->{db} = $db;
		$VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
		$fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
		my $chunk = new Process::MpiChunk($level, $VARS);
		push(@chunks, $chunk);
	    }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{db
			chunk
			the_void
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    #==BLAST ANALYSIS HERE
	    #-blastn search the file against ESTs
	    my $res_dir;
	    my $blastn_keepers = [];
	    if ($db) {
	       GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
	       ($blastn_keepers, $res_dir) = GI::blastn_as_chunks($chunk,
								  $db,
								  $the_void,
								  $safe_seq_id,
								  \%CTL_OPT,
								  $self->{RANK},
								  $LOG,
								  $LOG_FLAG
								  );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (res_dir => $res_dir,
			blastn_keepers => $blastn_keepers
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		push(@{$VARS->{$key}}, $self->{RESULTS}->{$key});
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 9) {	#collect blastn
	 $level_status = 'collecting blastn reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			res_dir
			blastn_keepers
			masked_fasta
			fasta_t_index
			holdover_blastn
			the_void
			q_seq_ref
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $masked_fasta = $VARS->{masked_fasta};
	    my $fasta_t_index = $VARS->{fasta_t_index};
	    my $holdover_blastn = $VARS->{holdover_blastn};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};

	    $blastn_keepers = GI::collect_blast($chunk,
						$blastn_keepers,
						$res_dir,
						$LOG
						);

	    #==merge in heldover Phathits from last round
	    if ($chunk->number != 0) { #if not first chunk
		#reviews heldover blast hits,
		#then merges and reblasts them if they cross the divide
		$blastn_keepers = GI::merge_resolve_hits($masked_fasta,
							 $fasta_t_index,
							 $blastn_keepers,
							 $holdover_blastn,
							 $the_void,
							 \%CTL_OPT,
							 'blastn',
							 $LOG
							 );
	    }
	    
	    #separate out hits too close to chunk divide to be run with exonerate
	    $holdover_blastn = [];
	    if (not $chunk->is_last) {
		($blastn_keepers, $holdover_blastn) = GI::process_the_chunk_divide($chunk,
										   $CTL_OPT{'split_hit'},
										   $CTL_OPT{'pred_flank'},
										   1, #treat strands independently
										   [$blastn_keepers]
										   );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (holdover_blastn => $holdover_blastn,
			blastn_keepers => $blastn_keepers,
			res_dir => undef
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
      elsif ($level == 10) {	#tblastx
	 $level_status = 'doing tblastx of alt-ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{tblastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
            foreach my $db (@{$VARS->{CTL_OPT}{_a_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk},
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'tblastx');

                next if($fin{$blast_finished});

                $db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
                $VARS->{db} = $db;
                $VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
                $fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
                my $chunk = new Process::MpiChunk($level, $VARS);
                push(@chunks, $chunk);
            }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{db
			chunk
			the_void
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    my $res_dir;
	    my $tblastx_keepers = [];
	    if ($db) {
	       GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
	       ($tblastx_keepers, $res_dir) = GI::tblastx_as_chunks($chunk,
								    $db,
								    $the_void,
								    $safe_seq_id,
								    \%CTL_OPT,
								    $self->{RANK},
								    $LOG,
								    $LOG_FLAG
								    );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (res_dir => $res_dir,
			tblastx_keepers => $tblastx_keepers
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		push(@{$VARS->{$key}}, $self->{RESULTS}->{$key});
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 11) {	#collect tblastx
	 $level_status = 'collecting tblastx reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk			    
			res_dir
			tblastx_keepers
			masked_fasta
			fasta_a_index
			holdover_tblastx
			the_void
			q_seq_ref
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $masked_fasta = $VARS->{masked_fasta};
	    my $fasta_a_index = $VARS->{fasta_a_index};
	    my $holdover_tblastx = $VARS->{holdover_tblastx};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};

	    $tblastx_keepers = GI::collect_blast($chunk,
						 $tblastx_keepers,
						 $res_dir,
						 $LOG
						 );
	    
	    #==merge in heldover Phathits from last round
	    if ($chunk->number != 0) { #if not first chunk
		#reviews heldover blast hits,
		#then merges and reblasts them if they cross the divide
		$tblastx_keepers = GI::merge_resolve_hits($masked_fasta,
							  $fasta_a_index,
							  $tblastx_keepers,
							  $holdover_tblastx,
							  $the_void,
							  \%CTL_OPT,
							  'tblastx',
							  $LOG
							  );
	    }
	    
	    #separate out hits too close to chunk divide to be run with exonerate
	    $holdover_tblastx = [];
	    if (not $chunk->is_last) {
		($tblastx_keepers, $holdover_tblastx) = GI::process_the_chunk_divide($chunk,
										     $CTL_OPT{'split_hit'},
										     $CTL_OPT{'pred_flank'},
										     1, #treat strands independently
										     [$tblastx_keepers]
										     );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (holdover_tblastx => $holdover_tblastx,
			tblastx_keepers => $tblastx_keepers,
			res_dir => undef
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
      elsif ($level == 12) {	#blastx
	 $level_status = 'doing blastx of proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{blastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_p_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk},
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'blastx');
		
                next if($fin{$blast_finished});

                $db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
                $VARS->{db} = $db;
                $VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
                $fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
                my $chunk = new Process::MpiChunk($level, $VARS);
                push(@chunks, $chunk);
            }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{db
			chunk
			the_void
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    #==BLAST ANALYSIS HERE
	    #-blastx search the file against proteins
	    my $res_dir;
	    my $blastx_keepers = [];
	    if ($db) {
		GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
		($blastx_keepers, $res_dir) = GI::blastx_as_chunks($chunk,
								   $db,
								   $the_void,
								   $safe_seq_id,
								   \%CTL_OPT,
								   $self->{RANK},
								   $LOG,
								   $LOG_FLAG
								   );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (res_dir => $res_dir,
			blastx_keepers => $blastx_keepers
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		push(@{$VARS->{$key}}, $self->{RESULTS}->{$key});
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 13) {	#collect blastx
	 $level_status = 'collecting blastx reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			res_dir
			blastx_keepers
			masked_fasta
			fasta_p_index
			holdover_blastx
			the_void
			q_seq_ref
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $masked_fasta = $VARS->{masked_fasta};
	    my $fasta_p_index = $VARS->{fasta_p_index};
	    my $holdover_blastx = $VARS->{holdover_blastx};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};

	    $blastx_keepers = GI::collect_blast($chunk,
						$blastx_keepers,
						$res_dir,
						$LOG
						);

	    #==merge in heldover Phathits from last round
	    if ($chunk->number != 0) { #if not first chunk
		#reviews heldover blast hits,
		#then merges and reblasts them if they cross the divide
		$blastx_keepers = GI::merge_resolve_hits($masked_fasta,
							 $fasta_p_index,
							 $blastx_keepers,
							 $holdover_blastx,
							 $the_void,
							 \%CTL_OPT,
							 'blastx',
							 $LOG
							 );
	    }
	    
	    #separate out hits too close to chunk divide to be run with exonerate
	    $holdover_blastx = [];
	    if (not $chunk->is_last) {
		($blastx_keepers, $holdover_blastx) = GI::process_the_chunk_divide($chunk,
										   $CTL_OPT{'split_hit'},
										   $CTL_OPT{'pred_flank'},
										   1, #treat strands independently
										   [$blastx_keepers]
										   );
	    }
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = (holdover_blastx => $holdover_blastx,
			blastx_keepers => $blastx_keepers,
			res_dir => undef
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
      elsif ($level == 14) {	#exonerate ESTs
	 $level_status = 'polishig ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
		        blastn_keepers
			holdover_blastn
			the_void
			q_seq_ref
			fasta
			fasta_t_index
			seq_id
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $holdover_blastn = $VARS->{holdover_blastn};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $fasta = $VARS->{fasta};
	    my $fasta_t_index = $VARS->{fasta_t_index};
	    my $seq_id = $VARS->{seq_id};
	    my $LOG = $VARS->{LOG};

	    #-polish blastn hits with exonerate
	    my $exonerate_e_data = [];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
		$exonerate_e_data = GI::polish_exonerate($chunk,
							 $fasta,
							 $blastn_keepers,
							 $fasta_t_index,
							 $the_void,
							 'e',
							 $CTL_OPT{exonerate},
							 $CTL_OPT{pcov_blastn},
							 $CTL_OPT{pid_blastn},
							 $CTL_OPT{en_score_limit},
							 $CTL_OPT{en_matrix},
							 $CTL_OPT{pred_flank},
							 $CTL_OPT{est_forward},
							 $CTL_OPT{_to_exonerate}{$seq_id},
							 $LOG
							 );
	    }

	    #-clean the blastn hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning blastn...\n" unless $main::quiet;
	    $blastn_keepers = GI::clean_blast_hits($blastn_keepers,
						   $CTL_OPT{pcov_blastn},
						   $CTL_OPT{pid_blastn},
						   $CTL_OPT{eval_blastn},
						   1 #contiguity flag
						   );

	    $blastn_keepers = GI::combine($blastn_keepers, $holdover_blastn);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (exonerate_e_data => $exonerate_e_data,
			blastn_keepers => $blastn_keepers,
			holdover_blastn => $holdover_blastn,
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
      elsif ($level == 15) {	#exonerate alt-ESTs
	 $level_status = 'polishing alt-ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{tblastx_keepers
			holdover_tblastx
			the_void
			q_seq_ref
			fasta
			fasta_a_index
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $holdover_tblastx = $VARS->{holdover_tblastx};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $fasta = $VARS->{fasta};
	    my $fasta_t_index = $VARS->{fasta_t_index};
	    my $LOG = $VARS->{LOG};

	    #-polish tblastx hits with exonerate
	    my $exonerate_a_data = [];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
		#$exonerate_a_data = GI::polish_exonerate('',
		#                                        $fasta,
		#					 $tblastx_keepers,
		#					 $fasta_t_index,
		#					 $the_void,
		#					 'a',
		#					 $CTL_OPT{exonerate},
		#					 $CTL_OPT{pcov_tblastx},
		#					 $CTL_OPT{pid_tblastx},
		#					 $CTL_OPT{en_score_limit},
		#					 $CTL_OPT{en_matrix},
		#					 $CTL_OPT{pred_flank},
		#					 $CTL_OPT{est_forward},
		#                                        [], #empty
		#					 $LOG
		#					 );
	    }
	    
	    #-clean the tblastx hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning tblastx...\n" unless $main::quiet;
	    $tblastx_keepers = GI::clean_blast_hits($tblastx_keepers,
						    $CTL_OPT{pcov_tblastx},
						    $CTL_OPT{pid_tblastx},
						    $CTL_OPT{eval_tblastx},
						    1 #contiguity flag
						    );

	    $tblastx_keepers = GI::combine($tblastx_keepers, $holdover_tblastx);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (exonerate_a_data => $exonerate_a_data,
			tblastx_keepers => $tblastx_keepers,
			holdover_tblastx => $holdover_tblastx,
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
      elsif ($level == 16) {	#exonerate proteins
	 $level_status = 'polishing proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{blastx_keepers
			holdover_blastx
			the_void
			q_seq_ref
			fasta
			fasta_p_index
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};;
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $holdover_blastx = $VARS->{holdover_blastx};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $fasta = $VARS->{fasta};
	    my $fasta_p_index = $VARS->{fasta_p_index};
	    my $LOG = $VARS->{LOG};
      
	    #-make a multi-fasta of the seqs in the blastx_clusters 
	    #-polish the blastx hits with exonerate
	    my $exonerate_p_data =[];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
		$exonerate_p_data = GI::polish_exonerate('',
							 $fasta,
							 $blastx_keepers,
							 $fasta_p_index,
							 $the_void,
							 'p',
							 $CTL_OPT{exonerate},
							 $CTL_OPT{pcov_blastx},
							 $CTL_OPT{pid_blastx},
							 $CTL_OPT{ep_score_limit},
							 $CTL_OPT{ep_matrix},
							 $CTL_OPT{pred_flank},
							 $CTL_OPT{est_forward},
							 [], #empty
							 $LOG
							);
	    }

	    #-clean the blastx hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning blastx...\n" unless $main::quiet;
	    $blastx_keepers = GI::clean_blast_hits($blastx_keepers,
						   $CTL_OPT{pcov_blastx},
						   $CTL_OPT{pid_blastx},
						   $CTL_OPT{eval_blastx},
						   0 #contiguity flag
						   );

	    $blastx_keepers = GI::combine($blastx_keepers, $holdover_blastx);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (exonerate_p_data => $exonerate_p_data,
			blastx_keepers => $blastx_keepers,
			holdover_blastx => [],
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
      elsif ($level == 17) {	#process chunk divide
	 $level_status = 'processing the chunk divide';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			the_void
			q_seq_ref
			masked_fasta
			preds
			blastn_keepers
			blastx_keepers
			tblastx_keepers
			exonerate_e_data
			exonerate_a_data
			exonerate_p_data
			holdover_blastn
			holdover_blastx
			holdover_tblastx
			holdover_pred
			holdover_est_gff
			holdover_altest_gff
			holdover_prot_gff
			holdover_pred_gff
			holdover_model_gff
			fasta_p_index
			fasta_t_index
			fasta_a_index
			GFF_DB
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $masked_fasta = $VARS->{masked_fasta};
	    my $preds = $VARS->{preds};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $exonerate_a_data = $VARS->{exonerate_a_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $holdover_blastn = $VARS->{holdover_blastn};
	    my $holdover_blastx = $VARS->{holdover_blastx};
	    my $holdover_tblastx = $VARS->{holdover_tblastx};
	    my $holdover_pred = $VARS->{holdover_pred};
	    my $holdover_est_gff = $VARS->{holdover_est_gff};
	    my $holdover_altest_gff = $VARS->{holdover_altest_gff};
	    my $holdover_prot_gff = $VARS->{holdover_prot_gff};
	    my $holdover_pred_gff = $VARS->{holdover_pred_gff};
	    my $holdover_model_gff = $VARS->{holdover_model_gff};
	    my $fasta_p_index = $VARS->{fasta_p_index};
	    my $fasta_t_index = $VARS->{fasta_t_index};
	    my $fasta_a_index = $VARS->{fasta_a_index};
	    my $GFF_DB = $VARS->{GFF_DB};
	    my $LOG = $VARS->{LOG};
   
	    #-get only those predictions on the chunk
	    my $preds_on_chunk = GI::get_preds_on_chunk($preds,
							$chunk
						       );

	    #==GFF3 passthrough of evidence
	    my $prot_gff_keepers = [];
	    my $est_gff_keepers = [];
	    my $altest_gff_keepers = [];
	    my $model_gff_keepers = [];
	    my $pred_gff_keepers = [];
	    if ($CTL_OPT{go_gffdb}) {
	       #-protein evidence passthraough
	       $prot_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							      $q_seq_ref,
							      'protein'
							      );
	       #-est evidence passthrough
	       $est_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							     $q_seq_ref,
							     'est'
							     );
	       #-altest evidence passthrough
	       $altest_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
								$q_seq_ref,
								'altest'
								);
	       #-gff gene annotation passthrough here
	       $model_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							       $q_seq_ref,
							       'model'
							       );
	       #-pred passthrough
	       $pred_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							      $q_seq_ref,
							      'pred'
							      );
	   }
	    
	    #combine remaining holdover types
	    push(@{$preds_on_chunk}, @{$holdover_pred});
	    push(@{$pred_gff_keepers}, @{$holdover_pred_gff});
	    push(@{$est_gff_keepers}, @{$holdover_est_gff});
	    push(@{$altest_gff_keepers}, @{$holdover_altest_gff});
	    push(@{$prot_gff_keepers}, @{$holdover_prot_gff});
	    push(@{$model_gff_keepers}, @{$holdover_model_gff});
	    
	    #clear holdovers
	    @{$holdover_pred} = ();
	    @{$holdover_est_gff} = ();
	    @{$holdover_altest_gff} = ();
	    @{$holdover_prot_gff} = ();
	    @{$holdover_pred_gff} = ();
	    @{$holdover_model_gff} = ();
	 
	    #==PROCESS HITS CLOSE TO CODE DIVISIONS
	    #holdover hits that are too close to the divide for review with next chunk
	    if (not $chunk->is_last) { #if not last chunk
		($blastn_keepers,
		 $blastx_keepers,
		 $tblastx_keepers,
		 $preds_on_chunk,
		 $est_gff_keepers,
		 $altest_gff_keepers,
		 $prot_gff_keepers,
		 $pred_gff_keepers,
		 $model_gff_keepers,
		 $exonerate_e_data,
		 $exonerate_a_data,
		 $exonerate_p_data,
		 $holdover_blastn,
		 $holdover_blastx,
		 $holdover_tblastx,
		 $holdover_pred,
		 $holdover_est_gff,
		 $holdover_altest_gff,
		 $holdover_prot_gff,
		 $holdover_pred_gff,
		 $holdover_model_gff
		 ) = GI::process_the_chunk_divide($chunk,
						  $CTL_OPT{'split_hit'},
						  $CTL_OPT{'pred_flank'},
						  1,
						  [$blastn_keepers,
						   $blastx_keepers,
						   $tblastx_keepers,
						   $preds_on_chunk,
						   $est_gff_keepers,
						   $altest_gff_keepers,
						   $prot_gff_keepers,
						   $pred_gff_keepers,
						   $model_gff_keepers,
						   $exonerate_e_data,
						   $exonerate_a_data,
						   $exonerate_p_data]
						  );
	    }

	    #Shatter hits. This is only for prokaryotic organisms.
	    #Flip strand of blastn where appropriate.
	    #This is done on blastn hits because exonerate is skipped.
	    #I shatter after processing the chunk divide to avoid weird
	    #complications from flipping on only one side of a split HSP
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$blastn_keepers  = PhatHit_utils::shatter_all_hits($blastn_keepers);
		$blastx_keepers  = PhatHit_utils::shatter_all_hits($blastx_keepers);
		$tblastx_keepers = PhatHit_utils::shatter_all_hits($tblastx_keepers);
		#this checks the open reading frame and can flip the hit
		foreach my $phat_hit (@$blastn_keepers){
		    $phat_hit = PhatHit_utils::copy($phat_hit, 'both')
			if exonerate::splice_info::needs_to_be_revcomped($phat_hit);
		}
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (prot_gff_keepers => $prot_gff_keepers,
			est_gff_keepers => $est_gff_keepers,
			altest_gff_keepers => $altest_gff_keepers,
			model_gff_keepers => $model_gff_keepers,
			pred_gff_keepers => $pred_gff_keepers,
			preds_on_chunk => $preds_on_chunk,
			blastn_keepers => $blastn_keepers,
			blastx_keepers => $blastx_keepers,
			tblastx_keepers => $tblastx_keepers,
			exonerate_e_data => $exonerate_e_data,
			exonerate_a_data => $exonerate_a_data,
			exonerate_p_data => $exonerate_p_data,
			holdover_est_gff => $holdover_est_gff,
			holdover_altest_gff => $holdover_altest_gff,
			holdover_prot_gff => $holdover_prot_gff,
			holdover_pred_gff => $holdover_pred_gff,
			holdover_model_gff => $holdover_model_gff,
			holdover_pred => $holdover_pred,
			holdover_blastn => $holdover_blastn,
			holdover_blastx => $holdover_blastx,
			holdover_tblastx => $holdover_tblastx
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
      elsif ($level == 18) {	#prep hint clusters
	 $level_status = 'preparing evidence clusters for annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw{q_seq_ref
			tblastx_keepers
			blastx_keepers
			blastn_keepers
			exonerate_e_data
			exonerate_a_data
			exonerate_p_data
			preds_on_chunk
			est_gff_keepers
			altest_gff_keepers
			prot_gff_keepers
			pred_gff_keepers
			model_gff_keepers
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $exonerate_a_data = $VARS->{exonerate_a_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $preds_on_chunk = $VARS->{preds_on_chunk};
	    my $est_gff_keepers = $VARS->{est_gff_keepers};
	    my $altest_gff_keepers = $VARS->{altest_gff_keepers};
	    my $prot_gff_keepers = $VARS->{prot_gff_keepers};
	    my $pred_gff_keepers = $VARS->{pred_gff_keepers};
	    my $model_gff_keepers = $VARS->{model_gff_keepers};

	    #trim evidence down to size if specified
	    if($CTL_OPT{fast}){
		my @sets = ($blastn_keepers,
			    $tblastx_keepers,
			    $blastx_keepers,
			    $exonerate_e_data,
			    $exonerate_a_data,
			    $exonerate_p_data,
			    $est_gff_keepers,
			    $altest_gff_keepers,
			    $prot_gff_keepers);

		#replace actual values
		foreach my $set (@sets) {
		    @$set = map {@$_} @{clean_and_cluster($set, ,$q_seq_ref, 10)};
		}
	    }

	    #combine final data sets
	    print STDERR "Preparing evidence for hint based annotation\n" unless($main::quiet);
	    my $final_est = GI::combine($exonerate_e_data,
					$est_gff_keepers
					);

	    $final_est = GI::combine($final_est,
				     $blastn_keepers
				     ) if($CTL_OPT{organism_type} eq 'prokaryotic');


	    my $final_altest = GI::combine($tblastx_keepers,
					   $exonerate_a_data,
					   $altest_gff_keepers
					  );

	    my $final_prot = GI::combine($blastx_keepers,
					 $exonerate_p_data,
					 $prot_gff_keepers
					);

	    my $final_pred = GI::combine($preds_on_chunk,
					 $pred_gff_keepers
					);

	    #group evidence for annotation
	    my $all_data = maker::auto_annotator::prep_hits($final_prot,
							    $final_est,
							    $final_altest,
							    $final_pred,
							    $model_gff_keepers,
							    $q_seq_ref,
							    $CTL_OPT{single_exon},
							    $CTL_OPT{single_length},
							    $CTL_OPT{pred_flank},
							    $CTL_OPT{organism_type},
							    $CTL_OPT{est_forward}
							    );
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (all_data => $all_data);
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
      elsif ($level == 19) {	#annotate transcripts
	 $level_status = 'annotating transcripts';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{trans} = {}; #reset
	    my @data_sets;
	    for(my $i = 0; $i < @{$VARS->{all_data}}; $i++){
		my $j = $i % $VARS->{CTL_OPT}->{_mpi_size};
		push(@{$data_sets[$j]}, $VARS->{all_data}->[$i]);
	    }

	    foreach my $dc (@data_sets){
                $VARS->{dc} = $dc;
		$VARS->{LOG_FLAG} = (!@chunks) ? 1 : 0;
                my $chunk = new Process::MpiChunk($level, $VARS);
                push(@chunks, $chunk);
            }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{the_void
			q_def
			seq_id
			q_seq_ref
			masked_total_seq
			dc
			LOG
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $masked_total_seq = $VARS->{masked_total_seq};
	    my $q_def = $VARS->{q_def};
	    my $seq_id = $VARS->{seq_id};
	    my $dc = $VARS->{dc};
	    my $the_void = $VARS->{the_void};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};
	    
	    #####working here###########
	    #==MAKER hint based predictions and annotations built here
	    #process transcripts
	    print STDERR "Making transcripts\n" unless($main::quiet || !$LOG_FLAG);
	    my $trans = maker::auto_annotator::annotate_trans($q_seq_ref,
							      $masked_total_seq,
							      $q_def,
							      $seq_id,
							      $dc,
							      $the_void,
							      \%CTL_OPT,
							      $LOG
							      );
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (trans => $trans);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		while (my $src = each %{$self->{RESULTS}->{$key}}){
		    push(@{$VARS->{$key}->{$src}}, @{$self->{RESULTS}->{$key}->{$src}});
		}
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 20) {	#grouping transcripts into genes
	 $level_status = 'clustering transcripts into genes for annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw{trans
			all_data
			q_seq_ref
			seq_id
			chunk
			build
			the_void
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $trans = $VARS->{trans};
	    my $all_data = $VARS->{all_data};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $seq_id = $VARS->{seq_id};
	    my $chunk = $VARS->{chunk};
	    my $build = $VARS->{build};
	    my $the_void = $VARS->{the_void};

	    #group transcriipts into genes
	    print STDERR "Processing transcripts into genes\n" unless($main::quiet);
	    my $annotations = maker::auto_annotator::annotate_genes($trans,
								    $all_data,
								    $q_seq_ref,
								    $seq_id,
								    $chunk->number(),
								    $build,
								    $the_void,
								    \%CTL_OPT
								    );
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (annotations   => $annotations);
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
      elsif ($level == 21) {	#adding quality control statistics
	 $level_status = 'adding statistics to annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my @data_sets;
	    my $i = 0;
	    while (my $key = each %{$VARS->{annotations}}){ 
		foreach my $an (@{$VARS->{annotations}->{$key}}){
		    my $j = $i % $VARS->{CTL_OPT}->{_mpi_size};
		    push(@{$data_sets[$j]{$key}}, $an);
		    $i++;
		}
	    }
	     
	    foreach my $an (@data_sets){
		$VARS->{an} = $an;
		$VARS->{LOG_FLAG} = (!@chunks) ? 1 : 0;
		my $chunk = new Process::MpiChunk($level, $VARS);
		push(@chunks, $chunk);
	    }
	    $VARS->{annotations} = {}; #reset
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw{q_seq_ref
			an
			LOG_FLAG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $an = $VARS->{an};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    #adds AED and other quality control statistics (can change names)
	    print STDERR "Calculating annotation quality statistics\n" unless($main::quiet || !$LOG_FLAG);
	    my $annotations = maker::auto_annotator::annotate_stats($an,
								    $q_seq_ref,
								    \%CTL_OPT);
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (annotations => $annotations);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		while (my $src = each %{$self->{RESULTS}->{$key}}){
		    push(@{$VARS->{$key}->{$src}}, @{$self->{RESULTS}->{$key}->{$src}});
		}
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 22) {	#deciding on final annotations
	 $level_status = 'choosing best annotation set';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw{annotations
			out_dir
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $annotations = $VARS->{annotations};
	    my $out_dir = $VARS->{out_dir};

	    #get best annotations
	    print STDERR "Choosing best annotations\n" unless($main::quiet);
	    my $maker_anno = maker::auto_annotator::best_annotations($annotations,
								     $out_dir,
								     \%CTL_OPT
								    );
	    
	    
	    #get best non-overlapping ab-inits
	    my $non_over = maker::auto_annotator::get_non_overlaping_abinits($maker_anno,
									     $annotations->{abinit},
									     \%CTL_OPT
									     );
	    
	    #add non-overlapping to final set if specified
	    if($CTL_OPT{keep_preds}){
		push(@$maker_anno, @$non_over);
		$non_over = [];
	    }

	    #map old IDs forward if specified
	    if($CTL_OPT{map_forward}){
		$maker_anno = maker::auto_annotator::map_forward($maker_anno,
								 $annotations->{model_gff}
								 );
	    }

	    #run evaluator if specified
	    if($CTL_OPT{evaluate}){
		#evaluator::evaluate::evaluate_maker_annotations($maker_anno,
		#						$q_seq_ref,
		#						$out_dir,
		#						$the_void,
		#						\%CTL_OPT
		#						);
	    }

	    #get AED scored preds for GFF3
	    my @scored_preds;
	    while(my $key = each %$annotations){
		next unless($key =~ /_abinit$|^pred_gff$/);
		
		foreach my $g (@{$annotations->{$key}}){
		    my @p_bases = map {$_->{p_base}} @{$g->{t_structs}};
		    push(@scored_preds, @p_bases);
		}
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (maker_anno => $maker_anno,
			non_over => $non_over,
			scored_preds => \@scored_preds);
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
      elsif ($level == 23) {	#local output
	 $level_status = 'processing chunk output';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{chunk
			maker_anno
			non_over
			annotations
			blastx_keepers
			blastn_keepers
			tblastx_keepers
			exonerate_e_data
			exonerate_a_data
			exonerate_p_data
			est_gff_keepers
			altest_gff_keepers
			prot_gff_keepers
			scored_preds
			p_fastas
			t_fastas
			GFF3}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my $chunk       = $VARS->{chunk};
	    my $maker_anno  = $VARS->{maker_anno};
	    my $non_over    = $VARS->{non_over};
	    my $annotations = $VARS->{annotations};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $exonerate_a_data = $VARS->{exonerate_a_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $est_gff_keepers = $VARS->{est_gff_keepers};
	    my $altest_gff_keepers = $VARS->{altest_gff_keepers};
	    my $prot_gff_keepers = $VARS->{prot_gff_keepers};
	    my $scored_preds = $VARS->{scored_preds};
	    my $p_fastas = $VARS->{p_fastas};
	    my $t_fastas = $VARS->{t_fastas};
	    my $GFF3 = $VARS->{GFF3};


	    #==OUTPUT DATA HERE      
	    #--- GFF3
	    $GFF3->add_genes($maker_anno);
	    $GFF3->add_phathits($blastn_keepers);
	    $GFF3->add_phathits($tblastx_keepers);
	    $GFF3->add_phathits($blastx_keepers);
	    $GFF3->add_phathits($exonerate_e_data);
	    $GFF3->add_phathits($exonerate_a_data);
	    $GFF3->add_phathits($exonerate_p_data);
	    $GFF3->add_phathits($est_gff_keepers);
	    $GFF3->add_phathits($altest_gff_keepers);
	    $GFF3->add_phathits($prot_gff_keepers);
	    $GFF3->add_phathits($scored_preds);
	    $GFF3->resolved_flag if (not $chunk->is_last); #adds ### between contigs
            
	    #--- building fastas for annotations (grows with iteration)
	    GI::maker_p_and_t_fastas($maker_anno,
				     $non_over,
				     $annotations->{abinit},
				     $p_fastas,
				     $t_fastas,
				   );
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (p_fastas => $p_fastas,
			t_fastas => $t_fastas
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
	    if ($VARS->{chunk} = $VARS->{fasta_chunker}->next_chunk) {
	       $next_level = 8;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 24) {	#global output
	 $level_status = 'processing contig output';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{the_void
			out_dir
			seq_id
			safe_seq_id
			q_seq_ref
			preds		  
			p_fastas
			t_fastas
			GFF3
			DS_CTL
			CTL_OPT
			LOG
			LOCK}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $out_dir = $VARS->{out_dir};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $preds = $VARS->{preds};
	    my $p_fastas = $VARS->{p_fastas};
	    my $t_fastas = $VARS->{t_fastas};
	    my $GFF3 = $VARS->{GFF3};
	    my $DS_CTL = $VARS->{DS_CTL};
	    my $LOG = $VARS->{LOG};
	    my $LOCK = $VARS->{LOCK};

	    #--- write fastas for ab-initio predictions
	    GI::write_p_and_t_fastas($p_fastas, $t_fastas, $safe_seq_id, $out_dir);
	    
	    #--- write GFF3 file
	    $GFF3->finalize();
	    
	    #--cleanup maker files created with each fasta sequence
	    File::Path::rmtree ($the_void) if $CTL_OPT{clean_up}; #rm temp directory
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
	 warn "WARNING: Invalid level for method _go() in Process::MpiChunk\n";
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
   die "FATAL: \'$flag\' is not a valid flag in MpiChunk _go!!\n";
}

#--------------------------------------------------------------
#reaches inside MpiTiers->{VARS} and places result in the 
#objects internal structure.
#By having this method work on MpiTiers->{VARS} rather than 
#simply returning a value for the results, we keep all contol
#of data structure inside of MpiChunk which makes for easier
#debugging

sub _result {
   my $self = shift;
   my $VARS = shift;		#this ia a hash reference;
   my $level = $self->{LEVEL};

   #only return results for finished/succesful chunks
   return if (! $self->finished || $self->failed);

   #do level specific result processing
   return $self->_go('result', $level, $VARS);
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
