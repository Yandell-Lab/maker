#! /usr/bin/perl -w

package Process::MpiChunk;

use strict;

use Error qw(:try :warndie);
use Error::Simple;
use Storable;
use runlog;
use GI;
use maker::auto_annotator;

#--set object variables for serialization of data
#this is needed when cloning an MpiChunk object
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
	 $self->{VARS} = {};
	 $self->{RESULTS} = {};

	 $self->_initialize($VARS);
      }
   }

   return $self;
}
#--------------------------------------------------------------
#this function/method is called by MpiTiers as part of MpiTiers
#initialization to prepare incoming data before building chunks.
#This method should not be called directly by the user or inside
#MpiChunks. Putting this preparation here as opposed to inside
#MpiTiers makes MpiChunks more portable and makes debugging
#easier.

sub prepare {
   my $self = shift;
   my $VARS = shift;

   #handle case of calling as function rather than method
   if (ref($self) ne "Process::MpiChunk") {
      $VARS = $self;
      $self = new Process::MpiChunk();
   }
   
   #==Prepare data here as part of initialization

   my $status = ''; #for error control
   #process fasta file
   try{
      #===
      $status = 'examining contents of the fasta file';
      #===

      #build uppercase fasta
      $VARS->{fasta} = Fasta::ucFasta(\$VARS->{fasta});

      #get fasta parts
      $VARS->{q_def} = Fasta::getDef(\$VARS->{fasta}); #Get fasta header
      $VARS->{q_seq_ref} = Fasta::getSeqRef(\$VARS->{fasta}); #Get reference to fasta sequence
      $VARS->{seq_id} = Fasta::def2SeqID($VARS->{q_def}); #Get sequence identifier
      $VARS->{safe_seq_id} = Fasta::seqID2SafeID($VARS->{seq_id}); #Get safe named identifier

      #set up base and void directories for output
      ($VARS->{out_dir}, $VARS->{the_void}) = $VARS->{DS_CTL}->seq_dirs($VARS->{seq_id});      

      #===
      $status = 'trying to build/process the run.log file';
      #===

      #-build and proccess the run log
      $VARS->{LOG} = runlog->new($VARS->{CTL_OPT},
				 {seq_id     => $VARS->{seq_id},
				  seq_length => length(${$VARS->{q_seq_ref}}),
				  out_dir    => $VARS->{out_dir},
				  the_void   => $VARS->{the_void},
				  fasta_ref  => \$VARS->{fasta}},
				 $VARS->{the_void}."/run.log"
				);
      
      ($VARS->{c_flag}, my $message) = $VARS->{LOG}->get_continue_flag();
      $VARS->{DS_CTL}->add_entry($VARS->{seq_id}, $VARS->{out_dir}, $message);
   }
   catch Error::Simple with {
      my $E = shift;

      $self->_handler($E, $status, 'throw');
   };

   return if($VARS->{c_flag} <= 0);

   #set up contig variables and GFF3 output object
   try{
      #===
      $status = 'setting ip contig variables and GFF3 output object';
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
      $VARS->{masked_total_seq} = '';
      $VARS->{p_fastas} = {};
      $VARS->{t_fastas} = {};

      #--prebuild an index of the databases
      my $proteins = $VARS->{CTL_OPT}{_protein};
      my $trans    = $VARS->{CTL_OPT}{_est};
      my $altests  = $VARS->{CTL_OPT}{_altest};
      GI::build_fasta_index($trans) if($trans); 
      GI::build_fasta_index($proteins) if($proteins);
      GI::build_fasta_index($altests) if($altests);

      #--other variables
      $VARS->{res_dir} = undef;
   }
   catch Error::Simple with{
      my $E = shift;

      $self->_handler($E, $status, 'throw');
   };
   
   return 1;
}
#--------------------------------------------------------------
#This function is called by MpiTiers.  It returns all the
#chunks for a given level.  This allows for the number of
#chunks created per level to be controlled within MpiChunks.

sub loader {
   my $self = shift;
   my $level = shift;
   my $VARS = shift;
   my $tID  = shift;

   #handle case of calling as function rather than method
   if (ref($self) ne "Process::MpiChunk") {
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
#Set up chunk specific variables.
#This method also reaches inside MpiTiers->{VARS} and pulls out
#all need variables.
#By Reaching into the MpiTiers object we keep control of the
#VARS datastructure within the MpiChunk (see method results).
#Another benefit is that variable order is no longer important,
#as it would be if variables were directly passed

sub _initialize{
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
   if (ref($self) ne "Process::MpiChunk") {
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
   my @chunks;
   my @args;
   my %results;

   my $level_status = '';

   try{
      if ($level == 0) {	#set up GFF3 output and fasta chunks
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
	    my $proteins = $VARS->{CTL_OPT}{_protein};
	    my $trans    = $VARS->{CTL_OPT}{_est};
	    my $altests  = $VARS->{CTL_OPT}{_altest};
	    my $fasta_t_index = GI::build_fasta_index($trans) if($trans); 
	    my $fasta_p_index = GI::build_fasta_index($proteins) if($proteins);
	    my $fasta_a_index = GI::build_fasta_index($altests) if($altests); 	 

	    my $chunk = $fasta_chunker->next_chunk();
	    #-------------------------CODE
	 
	    #------------------------RESULTS
	    %results = (GFF3 => $GFF3,
			fasta_chunker => $fasta_chunker,
			chunk => $chunk,
			fasta_t_index => $fasta_t_index,
			fasta_p_index => $fasta_p_index,
			fasta_a_index => $fasta_a_index
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 1) {	#do repeat masking
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
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_gff_keepers);
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
					       $CTL_OPT{force},
					       $LOG
					      );
	 
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_rb_keepers);
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
					       $CTL_OPT{force},
					       $LOG
					      );
	 
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_sp_keepers);
	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (rm_gff_keepers => $rm_gff_keepers,
			rm_rb_keepers => $rm_rb_keepers,
			rm_sp_keepers => $rm_sp_keepers,
			chunk => $chunk
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 2) {	#blastx repeat mask
	 $level_status = 'doing blastx repeats';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    foreach my $db (@{$VARS->{CTL_OPT}{r_db}}) {
	       $VARS->{db} = $db;
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
	    my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;

	    my $res_dir;
	    if ($CTL_OPT{_repeat_protein}) {
	       $res_dir = GI::blastx_as_chunks($chunk,
					       $db,
					       $the_void,
					       $safe_seq_id,
					       $CTL_OPT{_blastx},
					       $CTL_OPT{eval_rm_blastx},
					       $CTL_OPT{split_hit},
					       $CTL_OPT{cpus},
					       $CTL_OPT{repeat_protein},
					       $CTL_OPT{_formater},
					       $self->{RANK},
					       $CTL_OPT{force},
					       $LOG,
					       $LOG_FLAG
					      );
	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (res_dir => $res_dir);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 3) {	#collect blastx repeatmask
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
	    my $LOG = $VARS->{LOG};

	    my $rm_blastx_keepers = [];
	    if ($res_dir) {
	       $rm_blastx_keepers = GI::collect_blastx($chunk,
						       $res_dir,
						       $CTL_OPT{eval_rm_blastx},
						       $CTL_OPT{bit_rm_blastx},
						       $CTL_OPT{pcov_rm_blastx},
						       $CTL_OPT{pid_rm_blastx},
						       $CTL_OPT{split_hit},
						       $CTL_OPT{force},
						       $LOG
						      );
	 
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_blastx_keepers);
	    }
	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = ( rm_blastx_keepers=> $rm_blastx_keepers,
			 chunk => $chunk,
			 res_dir => $res_dir
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 4) {	#process all repeats
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
	    $masked_total_seq .= $chunk->seq();
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (rm_keepers => $rm_keepers,
			masked_total_seq => $masked_total_seq,
			rm_gff_keepers => [], #clear memory
			rm_rb_keepers => [],  #clear memory
			rm_sp_keepers => [],  #clear memory
			rm_blastx_keepers => [] #clear memory
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    if ($VARS->{chunk} = $VARS->{fasta_chunker}->next_chunk) {
	       $next_level = 1;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 5) {	#prep masked sequence and abinits
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
	 
	 
	    my $masked_fasta = Fasta::toFasta($q_def.' masked', \$masked_total_seq);
	    my $masked_file = $the_void."/query.masked.fasta";
	    FastaFile::writeFile(\$masked_fasta, $masked_file);
	 
	    my $unmasked_file = $the_void."/query.fasta";
	    FastaFile::writeFile(\$fasta, $unmasked_file);

	    #==ab initio predictions here
	    my $masked_preds = GI::abinits($masked_file,
					   $the_void,
					   $safe_seq_id.".masked",
					   \%CTL_OPT,
					   $LOG
					   );

            my $preds = GI::abinits($unmasked_file,
                                    $the_void,
                                    $safe_seq_id,
                                    \%CTL_OPT,
                                    $LOG
				    );

	    push(@$preds, @$masked_preds);

	    #==QRNA noncoding RNA prediction here
	    my $qra_preds = [];
	    #-------------------------CODE
	 
	    #------------------------RESULTS
	    %results = (masked_fasta => $masked_fasta,
			preds => $preds,
			qra_preds => $qra_preds
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 6) {	#prep new fasta chunks
	 $level_status = 'preparing new fasa chunks';
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

	    #------------------------RESULTS
	    %results = (fasta_chunker => $fasta_chunker,
			masked_fasta => $masked_fasta,
			chunk => $chunk
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 7) {	#blastn
	 $level_status = 'doing blastn of ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    foreach my $db (@{$VARS->{CTL_OPT}{e_db}}) {
	       $VARS->{db} = $db;
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
	    my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;

	    #==BLAST ANALYSIS HERE
	    #-blastn search the file against ESTs
	    my $res_dir;
	    if ($CTL_OPT{_est}) {
	       $res_dir = GI::blastn_as_chunks($chunk,
					       $db,
					       $the_void,
					       $safe_seq_id,
					       $CTL_OPT{_blastn},
					       $CTL_OPT{eval_blastn},
					       $CTL_OPT{split_hit},
					       $CTL_OPT{cpus},
					       $CTL_OPT{est},
					       $CTL_OPT{_formater},
					       $self->{RANK},
					       $CTL_OPT{force},
					       $LOG,
					       $LOG_FLAG
					      );

	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (res_dir => $res_dir,
			chunk => $chunk
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 8) {	#collect blastn
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
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};


	    my $blastn_keepers = [];
	    if ($res_dir) {
	       $blastn_keepers = GI::collect_blastn($chunk,
						    $res_dir,
						    $CTL_OPT{eval_blastn},
						    $CTL_OPT{bit_blastn},
						    $CTL_OPT{pcov_blastn},
						    $CTL_OPT{pid_blastn},
						    $CTL_OPT{split_hit},
						    $CTL_OPT{force},
						    $LOG
						   );
	    }
	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (blastn_keepers => $blastn_keepers,
			res_dir => $res_dir
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 9) {	#blastx
	 $level_status = 'doing blastx of proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    foreach my $db (@{$VARS->{CTL_OPT}{p_db}}) {
	       $VARS->{db} = $db;
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
	    my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;


	    #-blastx search the file against ESTs
	    my $res_dir;
	    if ($CTL_OPT{_protein}) {
	       $res_dir = GI::blastx_as_chunks($chunk,
					       $db,
					       $the_void,
					       $safe_seq_id,
					       $CTL_OPT{_blastx},
					       $CTL_OPT{eval_blastx},
					       $CTL_OPT{split_hit},
					       $CTL_OPT{cpus},
					       $CTL_OPT{protein},
					       $CTL_OPT{_formater},
					       $self->{RANK},
					       $CTL_OPT{force},
					       $LOG,
					       $LOG_FLAG
					      );

	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (res_dir => $res_dir);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 10) {	#collect blastx
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
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};


	    my $blastx_keepers = [];
	    if ($res_dir) {
	       $blastx_keepers = GI::collect_blastx($chunk,
						    $res_dir,
						    $CTL_OPT{eval_blastx},
						    $CTL_OPT{bit_blastx},
						    $CTL_OPT{pcov_blastx},
						    $CTL_OPT{pid_blastx},
						    $CTL_OPT{split_hit},
						    $CTL_OPT{force},
						    $LOG
						   );
	    }
	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (blastx_keepers => $blastx_keepers,
			res_dir => $res_dir
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 11) {	#tblastx
	 $level_status = 'doing tblastx of altESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    foreach my $db (@{$VARS->{CTL_OPT}{a_db}}) {
	       $VARS->{db} = $db;
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
	    my $LOG_FLAG = ($self->id =~ /^\d+\:\d+\:0$/) ? 1 : 0;


	    my $res_dir;
	    if ($CTL_OPT{_altest}) {
	       $res_dir = GI::tblastx_as_chunks($chunk,
						$db,
						$the_void,
						$safe_seq_id,
						$CTL_OPT{_tblastx},
						$CTL_OPT{eval_tblastx},
						$CTL_OPT{split_hit},
						$CTL_OPT{cpus},
						$CTL_OPT{altest},
						$CTL_OPT{_formater},
						$self->{RANK},
						$CTL_OPT{force},
						$LOG,
						$LOG_FLAG
					       );
	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (res_dir => $res_dir);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 12) {	#collect tblastx
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
			LOG
			CTL_OPT}
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};


	    my $tblastx_keepers = [];
	    if ($res_dir) {
	       $tblastx_keepers = GI::collect_tblastx($chunk,
						      $res_dir,
						      $CTL_OPT{eval_tblastx},
						      $CTL_OPT{bit_tblastx},
						      $CTL_OPT{pcov_tblastx},
						      $CTL_OPT{pid_tblastx},
						      $CTL_OPT{split_hit},
						      $CTL_OPT{force},
						      $LOG
						     );
	    }
	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (tblastx_keepers => $tblastx_keepers,
			res_dir => $res_dir
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 13) {	#process chunk divide
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
      
	    #==merge heldover Phathits from last round
	    if ($chunk->number != 0) { #if not first chunk
	       #reviews heldover blast hits,
	       #then merges and reblasts them if they cross the divide
	       ($blastn_keepers,
		$blastx_keepers,
		$tblastx_keepers) = GI::merge_resolve_hits(\$masked_fasta,
							   $fasta_t_index,
							   $fasta_p_index,
							   $fasta_a_index,
							   $blastn_keepers,
							   $blastx_keepers,
							   $tblastx_keepers,
							   $holdover_blastn,
							   $holdover_blastx,
							   $holdover_tblastx,
							   $the_void,
							   \%CTL_OPT,
							   $LOG
							  );
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
	    }
	 
	    #==PROCESS HITS CLOSE TO CODE DIVISIONS
	    #holdover hits that are too close to the divide for review with next chunk
	    if (not $chunk->is_last) { #if not last chunk
	       ($holdover_blastn,
		$holdover_blastx,
		$holdover_tblastx,
		$holdover_pred,
		$holdover_est_gff,
		$holdover_altest_gff,
		$holdover_prot_gff,
		$holdover_pred_gff,
		$holdover_model_gff,
		$blastn_keepers,
		$blastx_keepers,
		$tblastx_keepers,
		$preds_on_chunk,
		$est_gff_keepers,
		$altest_gff_keepers,
		$prot_gff_keepers,
		$pred_gff_keepers,
		$model_gff_keepers
	       ) = GI::process_the_chunk_divide($chunk,
						$CTL_OPT{'split_hit'},
						$blastn_keepers,
						$blastx_keepers,
						$tblastx_keepers,
						$preds_on_chunk,
						$est_gff_keepers,
						$altest_gff_keepers,
						$prot_gff_keepers,
						$pred_gff_keepers,
						$model_gff_keepers
					       );
	    }
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (prot_gff_keepers => $prot_gff_keepers,
			est_gff_keepers => $est_gff_keepers,
			altest_gff_keepers => $altest_gff_keepers,
			model_gff_keepers => $model_gff_keepers,
			pred_gff_keepers => $pred_gff_keepers,
			preds_on_chunk => $preds_on_chunk,
			blastn_keepers => $blastn_keepers,
			blastx_keepers => $blastx_keepers,
			tblastx_keepers => $tblastx_keepers,
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
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 14) {	#exonerate proteins
	 $level_status = 'doing exonerate of proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
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

	    #variables that are persistent outside of try block
	    my $blastx_data;
	    my $exonerate_p_data;
      
	    #-cluster the blastx hits
	    print STDERR "cleaning blastx...\n" unless $main::quiet;
      
	    my $blastx_clusters = cluster::clean_and_cluster($blastx_keepers,
							     $q_seq_ref,
							     10
							    );
	    undef $blastx_keepers; #free up memory
      
	    #-make a multi-fasta of the seqs in the blastx_clusters 
	    #-polish the blastx hits with exonerate
      
	    my $exoner_p_clust = GI::polish_exonerate($fasta,
						      $blastx_clusters,
						      $fasta_p_index,
						      $the_void,
						      5,
						      'p',
						      $CTL_OPT{exonerate},
						      $CTL_OPT{pcov_blastx},
						      $CTL_OPT{pid_blastx},
						      $CTL_OPT{ep_score_limit},
						      $CTL_OPT{ep_matrix},
						      $CTL_OPT{force},
						      $LOG
						     );
      
	    #flatten clusters
	    $blastx_data      = GI::flatten($blastx_clusters);
	    $exonerate_p_data = GI::flatten($exoner_p_clust, 'exonerate:p');
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (blastx_data => $blastx_data,
			exonerate_p_data => $exonerate_p_data
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 15) {	#exonerate ESTs
	 $level_status = 'doing exonerate of ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($level, $VARS);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw{tblastx_keepers
			blastn_keepers
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
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_ref = $VARS->{q_seq_ref};
	    my $fasta = $VARS->{fasta};
	    my $fasta_t_index = $VARS->{fasta_t_index};
	    my $LOG = $VARS->{LOG};


	    #-cluster the tblastx hits
	    print STDERR "cleaning tblastx...\n" unless $main::quiet;
	    my $tblastx_clusters = cluster::clean_and_cluster($tblastx_keepers,
							      $q_seq_ref,
							      10
							     );
	    undef $tblastx_keepers; #free up memory
	    
	    #flatten the clusters
	    my $tblastx_data = GI::flatten($tblastx_clusters);
	    
	    
	    #-cluster the blastn hits
	    print STDERR "cleaning blastn...\n" unless $main::quiet;
	    my $blastn_clusters = cluster::clean_and_cluster($blastn_keepers,
							     $q_seq_ref,
							     10
							    );
	    undef $blastn_keepers; #free up memory
	    
	    #-polish blastn hits with exonerate
	    my $exoner_e_clust = GI::polish_exonerate($fasta,
						      $blastn_clusters,
						      $fasta_t_index,
						      $the_void,
						      5,
						      'e',
						      $CTL_OPT{exonerate},
						      $CTL_OPT{pcov_blastn},
						      $CTL_OPT{pid_blastn},
						      $CTL_OPT{en_score_limit},
						      $CTL_OPT{en_matrix},
						      $CTL_OPT{force},
						      $LOG
						     ); 
	    
	    #flatten clusters
	    my $blastn_data      = GI::flatten($blastn_clusters);
	    my $exonerate_e_data = GI::flatten($exoner_e_clust,
					       'exonerate:e'
					      );
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (blastn_data => $blastn_data,
			exonerate_e_data => $exonerate_e_data,
			tblastx_data => $tblastx_data
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 16) {	#annotations
	 $level_status = 'calculating annotations';
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
			out_dir
			build
			fasta
			masked_fasta
			tblastx_data
			blastx_data
			exonerate_e_data
			exonerate_p_data
			preds_on_chunk
			est_gff_keepers
			altest_gff_keepers
			prot_gff_keepers
			pred_gff_keepers
			model_gff_keepers
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
	    my $build = $VARS->{build};
	    my $fasta = $VARS->{fasta};
	    my $masked_fasta = $VARS->{masked_fasta};
	    my $tblastx_data = $VARS->{tblastx_data};
	    my $blastx_data = $VARS->{blastx_data};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $preds_on_chunk = $VARS->{preds_on_chunk};
	    my $est_gff_keepers = $VARS->{est_gff_keepers};
	    my $altest_gff_keepers = $VARS->{altest_gff_keepers};
	    my $prot_gff_keepers = $VARS->{prot_gff_keepers};
	    my $pred_gff_keepers = $VARS->{pred_gff_keepers};
	    my $model_gff_keepers = $VARS->{model_gff_keepers};
	    my $LOG = $VARS->{LOG};

	    #combine final data sets
	    my $final_est = GI::combine($exonerate_e_data,
					$est_gff_keepers
				       );
	    my $final_altest = GI::combine($tblastx_data,
					   $altest_gff_keepers
					  );
	    my $final_prot = GI::combine($blastx_data,
					 $exonerate_p_data,
					 $prot_gff_keepers
					);
	    my $final_pred = GI::combine($preds_on_chunk,
					 $pred_gff_keepers
					);

	    #####working here###########
	    #==MAKER annotations built here
	    #-auto-annotate the input file
	    my $annotations = maker::auto_annotator::annotate($fasta,
							      $masked_fasta,
							      $chunk->number(),
							      $final_prot,
							      $final_est,
							      $final_altest,
							      $final_pred,
							      $model_gff_keepers,
							      $the_void,
							      $build,
							      \%CTL_OPT,
							      $LOG
							     );
	    
	    my $maker_anno = maker::auto_annotator::best_annotations($annotations,
								     $out_dir,
								     \%CTL_OPT
								    );

	    my $non_over = maker::auto_annotator::get_non_overlaping_abinits($maker_anno,
									     $annotations->{abinit}
									     );	    
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (maker_anno => $maker_anno, annotations => $annotations, non_over => $non_over);
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 17) {	#local output
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
			blastx_data
			blastn_data
			tblastx_data
			exonerate_p_data
			exonerate_e_data
			est_gff_keepers
			altest_gff_keepers
			prot_gff_keepers
			pred_gff_keepers
			preds_on_chunk
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
	    my $blastx_data = $VARS->{blastx_data};
	    my $blastn_data = $VARS->{blastn_data};
	    my $tblastx_data = $VARS->{tblastx_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $est_gff_keepers = $VARS->{est_gff_keepers};
	    my $altest_gff_keepers = $VARS->{altest_gff_keepers};
	    my $prot_gff_keepers = $VARS->{prot_gff_keepers};
	    my $pred_gff_keepers = $VARS->{pred_gff_keepers};
	    my $preds_on_chunk = $VARS->{preds_on_chunk};
	    my $p_fastas = $VARS->{p_fastas};
	    my $t_fastas = $VARS->{t_fastas};
	    my $GFF3 = $VARS->{GFF3};


	    #==OUTPUT DATA HERE      
	    #--- GFF3
	    $GFF3->add_genes($maker_anno);
	    $GFF3->add_phathits($blastx_data);
	    $GFF3->add_phathits($blastn_data);
	    $GFF3->add_phathits($tblastx_data);
	    $GFF3->add_phathits($exonerate_p_data);
	    $GFF3->add_phathits($exonerate_e_data);
	    $GFF3->add_phathits($est_gff_keepers);
	    $GFF3->add_phathits($altest_gff_keepers);
	    $GFF3->add_phathits($prot_gff_keepers);
	    $GFF3->add_phathits($preds_on_chunk);
	    $GFF3->add_phathits($pred_gff_keepers);
	    $GFF3->resolved_flag if (not $chunk->is_last); #adds ### between contigs
            
	    #--- building fastas for annotations (grows with iteration)
	    GI::maker_p_and_t_fastas($maker_anno,
				     $non_over,
				     $annotations->{abinit},
				     $p_fastas,
				     $t_fastas,
				   );
	    #-------------------------CODE

	    #------------------------RESULTS
	    %results = (p_fastas => $p_fastas,
			t_fastas => $t_fastas
		       );
	    #------------------------RESULTS
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    if ($VARS->{chunk} = $VARS->{fasta_chunker}->next_chunk) {
	       $next_level = 7;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($level == 18) {	#global output
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
			CTL_OPT}
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


	    #--- write fastas for ab-initio predictions
	    GI::write_p_and_t_fastas($p_fastas, $t_fastas, $safe_seq_id, $out_dir);
	    
	    #--- write GFF3 file
	    $GFF3->finalize();
	    
	    #--cleanup maker files created with each fasta sequence
	    File::Path::rmtree ($the_void) if $CTL_OPT{clean_up}; #rm temp directory
	    
	    #-- write to DS log the finished files
	    #$DS_CTL->add_entry($seq_id, $out_dir, 'FINISHED');
	    
	    #--- clear the log variable
	    $VARS->{LOG} = undef;
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