#! /usr/bin/perl -w

package Process::mapChunk;

use strict;

use Error qw(:try);
use Error::Simple;
use Storable;

#--set object variables for serialization of data
#this is needed when cloning an mapChunk object
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
      if (ref $arg eq 'Process::mapChunk') {
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
#VARS datastructure within the mapChunk (see method results).
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
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $tier = $self;
      $self = new Process::mapChunk();
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
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $tier = $self;
      $self = new Process::mapChunk();
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
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $tier = $self;
      $self = new Process::mapChunk();
   }

   return $tier->{VARS}{c_flag};
}
#--------------------------------------------------------------
#this function/method is called by MpiTiers as part of MpiTiers
#initialization to prepare incoming data before building chunks.
#This method should not be called directly by the user or inside
#mapChunks. Putting this preparation here as opposed to inside
#MpiTiers makes mapChunks more portable and makes debugging
#easier.  It will always run on the root node before distributing
#the tier.

sub _prepare {
   my $self = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $VARS = $self;
      $self = new Process::mapChunk();
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
   $VARS->{holdover_est_gff}    = [];
   $VARS->{holdover_prot_gff}   = [];
   
   #-set up variables that are the result of chunk accumulation
   $VARS->{masked_total_seq} = '';
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
#chunks created per level to be controlled within mapChunks.

sub _loader {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;
   my $tID  = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $tID = $VARS;
      $VARS = $level;
      $level = $self;
      $self = new Process::mapChunk();
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
   if ($self ne "Process::mapChunk" && ref($self) ne "Process::mapChunk") {
      $VARS = $level;
      $level = $self;
      $self = new Process::mapChunk();
   }

   return $self->_go('flow', $level, $VARS);
}
#--------------------------------------------------------------
#initializes chunk variables, runs code, or returns flow
#depending on flag.
#This method is the core of the mapChunk object

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
	if ($level == 0) {	#examining contents of the fasta file and run log
	    $level_status = 'examining contents of the fasta file and run log';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
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
		my $fasta = Fasta::ucFasta(\$VARS->{fasta}); #build uppercase fasta
		my $DS_CTL = $VARS->{DS_CTL};
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		
		#get fasta parts
		my $q_def = Fasta::getDef(\$fasta); #Get fasta header
		my $q_seq_ref = Fasta::getSeqRef(\$fasta); #Get reference to fasta sequence
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
				       fasta_ref  => \$fasta},
				      "$out_dir/run.log"
				      );
		
		my $LOCK = $LOG->strip_off_lock();
		
		my ($c_flag, $message) = $LOG->get_continue_flag();
		$DS_CTL->add_entry($seq_id, $out_dir, $message) if($message);

				
		my $GFF3 = Dumper::GFF::GFFV3->new("$out_dir/$safe_seq_id.gff",
						   1,
						   $the_void
						   );
		
		$GFF3->set_current_contig($seq_id, $q_seq_ref);     

		#--build an index of the databases
		my $prot_dir = $VARS->{CTL_OPT}{_protein};
		my $tran_dir    = $VARS->{CTL_OPT}{_est};
		my $fasta_t_index = GI::build_fasta_index($tran_dir) if($tran_dir); 
		my $fasta_p_index = GI::build_fasta_index($prot_dir) if($prot_dir);

		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (fasta => $fasta,
			    q_def => $q_def,
			    q_seq_ref => $q_seq_ref,
			    seq_id => $seq_id,
			    safe_seq_id => $safe_seq_id,
			    out_dir => $out_dir,
			    the_void => $the_void,
			    c_flag => $c_flag,
			    LOG => $LOG,
			    LOCK => $LOCK,
			    GFF3 => $GFF3,
			    fasta_t_index => $fasta_t_index,
			    fasta_p_index => $fasta_p_index,
			    masked_fasta => $fasta
			    );
		#------------------------RESULTS
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
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{fasta
			    CTL_OPT
			    out_dir
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
		my $seq_id = $VARS->{seq_id};
		my $safe_seq_id = $VARS->{safe_seq_id};
		my $the_void = $VARS->{the_void};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $fasta = $VARS->{fasta};
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		
		my $fasta_chunker = new FastaChunker();
		$fasta_chunker->parent_fasta($fasta);
		$fasta_chunker->chunk_size($CTL_OPT{max_dna_len});
		$fasta_chunker->min_size($CTL_OPT{split_hit});
		$fasta_chunker->load_chunks();
				
		my $chunk = $fasta_chunker->next_chunk();
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (fasta_chunker => $fasta_chunker,
			    chunk => $chunk,
			    );
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 2) {	#blastn
	    $level_status = 'doing blastn of ESTs';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{chunk
			    the_void
			    safe_seq_id
			    fasta
			    fasta_t_index
			    holdover_blastn
			    q_seq_ref
			    LOG
			    CTL_OPT}
			 );
		#------------------------ARGS_IN
	    }
	    elsif ($flag eq 'run') {
		#-------------------------CODE
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		my $chunk = $VARS->{chunk};
		my $db = $CTL_OPT{est};
		my $the_void = $VARS->{the_void};
		my $safe_seq_id = $VARS->{safe_seq_id};
		my $masked_fasta = $VARS->{fasta};
		my $fasta_t_index = $VARS->{fasta_t_index};
		my $holdover_blastn = $VARS->{holdover_blastn};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $LOG = $VARS->{LOG};
		my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;
		
		#==BLAST ANALYSIS HERE
		GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
		
		#-blastn search the file against ESTs
		my $blastn_keepers = [];
		if ($CTL_OPT{_est}) {
		    my $res_dir = GI::blastn_as_chunks($chunk,
						       $db,
						       $CTL_OPT{est},					       
						       $the_void,
						       $safe_seq_id,
						       \%CTL_OPT,
						       $self->{RANK},
						       $LOG,
						       $LOG_FLAG
						       );
		    
		    $blastn_keepers = GI::collect_blastn($chunk,
							 $res_dir,
							 \%CTL_OPT,
							 $LOG
							 );
		    
		    #==merge in heldover Phathits from last round
		    if ($chunk->number != 0) { #if not first chunk
			#reviews heldover blast hits,
			#then merges and reblasts them if they cross the divide
			$blastn_keepers = GI::merge_resolve_hits(\$masked_fasta,
								 $fasta_t_index,
								 $blastn_keepers,
								 $holdover_blastn,
								 $the_void,
								 \%CTL_OPT,
								 'blastn',
								 $LOG
								 );
		    }		    
		}

		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (blastn_keepers => $blastn_keepers,
			    );
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 3) {	#blastx
	    $level_status = 'doing blastx of proteins';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{chunk
			    the_void
			    safe_seq_id
			    fasta
			    fasta_p_index
			    holdover_blastx
			    q_seq_ref
			    LOG
			    CTL_OPT}
			 );
		#------------------------ARGS_IN
	    }
	    elsif ($flag eq 'run') {
		#-------------------------CODE
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		my $chunk = $VARS->{chunk};
		my $db = $CTL_OPT{protein};
		my $the_void = $VARS->{the_void};
		my $safe_seq_id = $VARS->{safe_seq_id};
		my $masked_fasta = $VARS->{fasta};
		my $fasta_p_index = $VARS->{fasta_p_index};
		my $holdover_blastx = $VARS->{holdover_blastx};
		my $q_seq_ref = $VARS->{q_seq_ref};
		
		my $LOG = $VARS->{LOG};
		my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;
		
		#-blastx search the file against ESTs
		GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
		
		my $blastx_keepers = [];
		if ($CTL_OPT{_protein}) {
		    my$res_dir = GI::blastx_as_chunks($chunk,
						      $db,
						      $CTL_OPT{protein},
						      $the_void,
						      $safe_seq_id,
						      \%CTL_OPT,
						      $self->{RANK},
						      $LOG,
						      $LOG_FLAG
						      );
		    
		    $blastx_keepers = GI::collect_blastx($chunk,
							 $res_dir,
							 \%CTL_OPT,
							 $LOG
							 );
		    
		    
		    #==merge in heldover Phathits from last round
		    if ($chunk->number != 0) { #if not first chunk
			#reviews heldover blast hits,
			#then merges and reblasts them if they cross the divide
			$blastx_keepers = GI::merge_resolve_hits(\$masked_fasta,
								 $fasta_p_index,
								 $blastx_keepers,
								 $holdover_blastx,
								 $the_void,
								 \%CTL_OPT,
								 'blastx',
								 $LOG
								 );
		    }
		}
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (blastx_keepers => $blastx_keepers,
			    );
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 4) {	#exonerate proteins
	    $level_status = 'doing exonerate of proteins';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{blastx_keepers
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
		my $the_void = $VARS->{the_void};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $fasta = $VARS->{fasta};
		my $fasta_p_index = $VARS->{fasta_p_index};
		my $LOG = $VARS->{LOG};
		
		#-make a multi-fasta of the seqs in the blastx_clusters 
		#-polish the blastx hits with exonerate
		my $exonerate_p_data =[];
		
		if($CTL_OPT{organism_type} eq 'eukaryotic'){
		    $exonerate_p_data = GI::polish_exonerate($fasta,
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
							     $LOG
							     );
		}
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (exonerate_p_data => $exonerate_p_data);
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 5) {	#exonerate ESTs
	    $level_status = 'doing exonerate of ESTs';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{blastn_keepers
			    the_void
			    q_seq_ref
			    fasta
			    fasta_t_index
			    LOG
			    CTL_OPT}
			 );
		#------------------------ARGS_IN
	    }
	    elsif ($flag eq 'run') {
		#-------------------------CODE
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		my $blastn_keepers = $VARS->{blastn_keepers};
		my $the_void = $VARS->{the_void};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $fasta = $VARS->{fasta};
		my $fasta_t_index = $VARS->{fasta_t_index};
		my $LOG = $VARS->{LOG};
		
		#-polish blastn hits with exonerate
		my $exonerate_e_data = [];
		if($CTL_OPT{organism_type} eq 'eukaryotic'){
		    $exonerate_e_data = GI::polish_exonerate($fasta,
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
							     $LOG
							     );
		}
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (exonerate_e_data => $exonerate_e_data);
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 6) {	#process chunk divide
	    $level_status = 'processing the chunk divide';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{chunk
			    the_void
			    q_seq_ref
			    fasta
			    blastn_keepers
			    blastx_keepers
			    exonerate_e_data
			    exonerate_p_data
			    holdover_blastn
			    holdover_blastx
			    fasta_p_index
			    fasta_t_index
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
		my $masked_fasta = $VARS->{fasta};
		my $blastn_keepers = $VARS->{blastn_keepers};
		my $blastx_keepers = $VARS->{blastx_keepers};
		my $exonerate_e_data = $VARS->{exonerate_e_data};
		my $exonerate_p_data = $VARS->{exonerate_p_data};
		my $holdover_blastn = $VARS->{holdover_blastn};
		my $holdover_blastx = $VARS->{holdover_blastx};
		my $fasta_p_index = $VARS->{fasta_p_index};
		my $fasta_t_index = $VARS->{fasta_t_index};
		my $LOG = $VARS->{LOG};
		
		#==PROCESS HITS CLOSE TO CODE DIVISIONS
		#holdover hits that are too close to the divide for review with next chunk
		if (not $chunk->is_last) { #if not last chunk
		    ($blastn_keepers,
		     $blastx_keepers,
		     $exonerate_e_data,
		     $exonerate_p_data,
		     $holdover_blastn,
		     $holdover_blastx,
		     ) = GI::process_the_chunk_divide($chunk,
						      $CTL_OPT{'split_hit'},
						      $CTL_OPT{'pred_flank'},
						      1,
						      [$blastn_keepers,
						       $blastx_keepers,
						       $exonerate_e_data,
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
		    #this checks the open reading frame and can flip the hit
		    foreach my $phat_hit (@$blastn_keepers){
			$phat_hit = PhatHit_utils::copy($phat_hit, 'both')
			    if exonerate::splice_info::needs_to_be_revcomped($phat_hit);
		    }
		}
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (blastn_keepers => $blastn_keepers,
			    blastx_keepers => $blastx_keepers,
			    exonerate_e_data => $exonerate_e_data,
			    exonerate_p_data => $exonerate_p_data,
			    holdover_blastn => $holdover_blastn,
			    holdover_blastx => $holdover_blastx
			    );
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 7) {	#annotations
	    $level_status = 'calculating annotations';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{chunk
			    the_void
			    out_dir
			    q_seq_ref
			    build
			    fasta
			    blastx_keepers
			    blastn_keepers
			    exonerate_e_data
			    exonerate_p_data
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
		my $out_dir = $VARS->{out_dir};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $fasta = $VARS->{fasta};
		my $masked_fasta = $VARS->{fasta};
		my $blastx_keepers = $VARS->{blastx_keepers};
		my $blastn_keepers = $VARS->{blastn_keepers};
		my $exonerate_e_data = $VARS->{exonerate_e_data};
		my $exonerate_p_data = $VARS->{exonerate_p_data};
		my $LOG = $VARS->{LOG};
		
		#combine final data sets
		my $final_est = $exonerate_e_data;
		
		if($CTL_OPT{organism_type} eq 'prokaryotic'){
		    $final_est = GI::combine($blastn_keepers,
					     $final_est
					     );
		}
		
		my $final_prot = GI::combine($blastx_keepers,
					     $exonerate_p_data
					     );
		
		#####working here###########
		#==MAKER annotations built here
		#-auto-annotate the input file
		my $annotations = maker::auto_annotator::annotate($fasta,
								  $masked_fasta,
								  $chunk->number(),
								  $final_prot,
								  $final_est,
								  [],
								  [],
								  [],
								  $the_void,
								  1,
								  \%CTL_OPT,
								  $LOG
								  );
		
		#get best annotations
		my $maker_anno = maker::auto_annotator::best_annotations($annotations,
									 $out_dir,
									 \%CTL_OPT
									 );
		
		#map old IDs forward if specified
		if($CTL_OPT{map_forward}){
		    $maker_anno = maker::auto_annotator::map_forward($maker_anno,
								     $annotations->{model_gff}
								     );
		}
		
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (maker_anno => $maker_anno);
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 8) {	#local output
	    $level_status = 'processing chunk output';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
		push(@chunks, $chunk);
		#-------------------------CHUNKER
	    }
	    elsif ($flag eq 'init') {
		#------------------------ARGS_IN
		@args = (qw{chunk
			    maker_anno
			    blastx_keepers
			    blastn_keepers
			    exonerate_p_data
			    exonerate_e_data
			    p_fastas
			    t_fastas
			    GFF3
			    q_seq_ref
			    q_def
			    CTL_OPT}
			 );
		#------------------------ARGS_IN
	    }
	    elsif ($flag eq 'run') {
		#-------------------------CODE
		my $chunk       = $VARS->{chunk};
		my $maker_anno  = $VARS->{maker_anno};
		my $blastx_keepers = $VARS->{blastx_keepers};
		my $blastn_keepers = $VARS->{blastn_keepers};
		my $exonerate_p_data = $VARS->{exonerate_p_data};
		my $exonerate_e_data = $VARS->{exonerate_e_data};
		my $p_fastas = $VARS->{p_fastas};
		my $t_fastas = $VARS->{t_fastas};
		my $q_seq_ref = $VARS->{q_seq_ref};
		my $q_def = $VARS->{q_def};
		my $GFF3 = $VARS->{GFF3};
		my %CTL_OPT = %{$VARS->{CTL_OPT}};
		
		#==OUTPUT DATA HERE      
		#--- GFF3
		$GFF3->add_genes($maker_anno);
		$GFF3->add_phathits($blastx_keepers);
		$GFF3->add_phathits($blastn_keepers);
		$GFF3->add_phathits($exonerate_p_data);
		$GFF3->add_phathits($exonerate_e_data);
		$GFF3->resolved_flag if (not $chunk->is_last); #adds ### between contigs
		
		#--- building fastas for annotations (grows with iteration)
		GI::maker_p_and_t_fastas($maker_anno,
					 [],
					 [],
					 $p_fastas,
					 $t_fastas,
					 );

		#--mask annotations in current contig
		my @phathits;
		foreach my $a (@$maker_anno){
		    foreach my $s (@{$a->{t_structs}}){
			push(@phathits, $s->{hit});
		    }
		}

		my $masked_total_seq = repeat_mask_seq::mask_seq($$q_seq_ref, \@phathits);
		my $fasta = Fasta::toFasta($q_def, \$masked_total_seq);
		my $q_seq_ref = Fasta::getSeqRef(\$fasta); #Get reference to fasta sequence
		#-------------------------CODE
		
		#------------------------RESULTS
		%results = (p_fastas     => $p_fastas,
			    t_fastas     => $t_fastas,
			    fasta        => $fasta,
			    masked_fasta => $fasta,
			    q_seq_ref    => $q_seq_ref,
			    CTL_OPT      => \%CTL_OPT,
			    );
		#------------------------RESULTS
	    }
	    elsif ($flag eq 'flow') {
		#-------------------------NEXT_LEVEL
		if ($VARS->{chunk} = $VARS->{fasta_chunker}->next_chunk) {
		    $next_level = 3;
		}
		#-------------------------NEXT_LEVEL
	    }
	}
	elsif ($level == 9) {	#global output
	    $level_status = 'processing contig output';
	    if ($flag eq 'load') {
		#-------------------------CHUNKER
		my $chunk = new Process::mapChunk($level, $VARS);
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
	    warn "WARNING: Invalid level for method _go() in Process::mapChunk\n";
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
    die "FATAL: \'$flag\' is not a valid flag in mapChunk _go!!\n";
}

#--------------------------------------------------------------
#reaches inside MpiTiers->{VARS} and places result in the 
#objects internal structure.
#By having this method work on MpiTiers->{VARS} rather than 
#simply returning a value for the results, we keep all contol
#of data structure inside of mapChunk which makes for easier
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
