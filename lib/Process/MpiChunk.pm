package Process::MpiChunk;

use strict;

use Error qw(:try);
use Error::Simple;
use Storable qw(store);
use Process::MpiTiers;
use File::Copy;
use Carp;

#--set object variables for serialization of data
#this is needed when cloning an MpiChunk object
$Storable::forgive_me = 1; #allows serializaion of objects with GLOBs
$Storable::Deparse = 1; #now serializes CODE refs
$Storable::Eval= 1; #now serializes CODE refs

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
	 my $VARS           = {%$arg}; #forces copy of hash (1 level deep)	 
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
#Set up chunk specific variables. Only runs when $VARS exist
#This method also reaches inside MpiTiers->{VARS} and pulls out
#all need variables.
#By Reaching into the MpiTiers object we keep control of the
#VARS datastructure within the MpiChunk (see method results).
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
	 confess "FATAL: argument \`$key\` does not exist in MpiTier object\n";
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
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $self;
      $self = new Process::MpiChunk();
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
   my $tier_type = shift || 0;

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $VARS = $self;
      $self = new Process::MpiChunk();
   }

   #==Prepare data here as part of initialization

   $VARS->{c_flag} = 1; #always continue with this implementation and
   #let child nodes decide on running chunk.
   #Other implementations let the root node decide

   if($tier_type == 0){
      #instantiate empty LOG
      $VARS->{LOG} = undef;

      #get genome index
      my $g_index = $VARS->{g_index};
      delete($VARS->{g_index});

      #get fasta parts
      my $q_def = $VARS->{q_def};
      my $seq_id = Fasta::def2SeqID($q_def); #get sequence identifier
      my $safe_seq_id = Fasta::seqID2SafeID($seq_id); #Get safe named identifier
      my $q_seq_obj = $g_index->get_Seq_by_id($seq_id);

      #still no object? try rebuilding the index and try again
      if(!$q_seq_obj) {
	  for(my $i = 0; $i < 2 && !$q_seq_obj; $i++){
	      sleep 5;
	      print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
	      $g_index->reindex($i);
	      $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
	  }

	  if (not $q_seq_obj) {
	      print STDERR "stop here: $seq_id\n";
	      confess "ERROR: Fasta index error\n";
	  }
      }

      my $q_seq_length = $q_seq_obj->length(); #seq length

      #set values inside of tier
      $VARS->{seq_id} = $seq_id;
      $VARS->{safe_seq_id} = $safe_seq_id;
      $VARS->{q_seq_length} = $q_seq_length;

      #-set up variables that are the result of chunk accumulation
      $VARS->{p_fastas} = {};
      $VARS->{t_fastas} = {};
      $VARS->{n_fastas} = {};
      $VARS->{gff3_files} = [];
   }
   elsif($tier_type == 1){

   }
   elsif($tier_type == 2){

   }
   elsif($tier_type == 3){
      #-set up variables that are heldover from last chunk
      $VARS->{holdover_blastn}          = [];
      $VARS->{holdover_blastx}          = [];
      $VARS->{holdover_tblastx}         = [];
      $VARS->{holdover_exonerate_e}     = [];
      $VARS->{holdover_exonerate_p}     = [];
      $VARS->{holdover_exonerate_a}     = [];
      $VARS->{holdover_pred}            = [];
      $VARS->{holdover_est_gff}         = [];
      $VARS->{holdover_altest_gff}      = [];
      $VARS->{holdover_prot_gff}        = [];
      $VARS->{holdover_pred_gff}        = [];
      $VARS->{holdover_model_gff}       = [];

      #--other variables
      $VARS->{exonerate_e_data} = [];
      $VARS->{exonerate_p_data} = [];
      $VARS->{exonerate_a_data} = [];
      $VARS->{blastn_keepers}   = [];
      $VARS->{blastx_keepers}   = [];
      $VARS->{tblastx_keepers}  = [];
      $VARS->{res_dir} = [];
      $VARS->{edge_status} = {};
   }
   elsif($tier_type == 4){

   }

   return 1;
}
#--------------------------------------------------------------
#This function is called by MpiTiers.  It returns all the
#chunks for a given level.  This allows for the number of
#chunks created per level to be controlled within MpiChunks.
#called as {CHUNK_REF}

sub _loader {
   my $self = shift;
   my $VARS = shift;
   my $level = shift;
   my $tier_type = shift;
   my $tier  = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $tier_type;
      $tier_type = $level;
      $level = $VARS;
      $VARS = $self;
      $self = new Process::MpiChunk();
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
#cals _go('run') to run code

sub run {
   my $self = shift;
   my $rank = shift;

   my $VARS = $self->{VARS};
   my $level = $self->{LEVEL};
   my $tier_type = $self->{TIER_TYPE};

   $self->{RANK} = $rank if(defined($rank));

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
#called frm {CHUNK_REF}

sub _flow {
   my $self = shift;
   my $VARS = shift; 
   my $level = shift;
   my $tier_type = shift;
   my $tier = shift;

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $tier_type;
      $tier_type = $level;
      $level = $VARS;
      $VARS = $self;
      $self = new Process::MpiChunk();
   }

   return $self->_go('flow', $VARS, $level, $tier_type);
}
#--------------------------------------------------------------
#initializes chunk variables, runs code, or returns flow
#depending on flag.
#This method is the core of the MpiChunk object

sub _go {
   my $self = shift;
   my $flag = shift;
   my $VARS = shift;
   my $level = shift;
   my $tier_type = shift;

   $VARS = $self->{VARS} if(! defined($VARS));
   $level = $self->{LEVEL} if(! defined($level));
   $tier_type = $self->{TIER_TYPE} if(! defined($tier_type));

   my $next_level = $level + 1;
   my $result_stat = 1;
   my @chunks;
   my @args;
   my %results;

   my $level_status = '';

   try{
      ##
      ##TO TIER_TYPE 0
      ##
      if ($tier_type == 0 && $level == 0) { #examining contents of the fasta file and run log
	 $level_status = 'examining contents of the fasta file and run log';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(q_def
			seq_id
			safe_seq_id
			q_seq_length
			CTL_OPT
			DS_CTL)
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $q_def = $VARS->{q_def};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $DS_CTL = $VARS->{DS_CTL};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};

	    #set up base and void directories for output
	    my ($out_dir, $the_void) = $DS_CTL->seq_dirs($seq_id);
	    
	    #-build and proccess the run log
	    my $LOG = runlog->new(\%CTL_OPT,
				  {seq_id     => $seq_id,
				   seq_length => $q_seq_length,
				   out_dir    => $out_dir,
				   the_void   => $the_void},
				  "$out_dir/run.log"
				  );


	    my $LOCK = $LOG->strip_off_lock();
	    my ($c_flag, $message) = $LOG->get_continue_flag();
	    $DS_CTL->add_entry($seq_id, $out_dir, $message) if($message);

	    #write local contig fasta
	    my $g_index;
	    my $q_seq_obj;
	    my $unmasked_file = "$the_void/query.fasta";
	    if($c_flag > 0){

		#get genome index and object
		$g_index = GI::build_fasta_index($CTL_OPT{_g_db});
		$q_seq_obj = $g_index->get_Seq_by_id($seq_id);

		#no sequence, try again (probably NFS)
		for(my $i = 0; $i < 2 && !$q_seq_obj; $i++){
		    sleep 5;
		    print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
		    $g_index->reindex($i);
		    $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		}
		if (! $q_seq_obj) {
		    print STDERR "stop here: $seq_id\n";
		    confess "ERROR: Fasta index error\n";
		}

		open (my $FAS, "> $the_void/query.fasta");
		my $seq = $q_seq_obj->seq;
		print $FAS ${&Fasta::seq2fastaRef($seq_id, \$seq)};
		close ($FAS);
	    }
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = (unmasked_file => $unmasked_file,
			q_seq_obj => $q_seq_obj,
			out_dir => $out_dir,
			the_void => $the_void,
			c_flag => $c_flag,
			LOG => $LOG,
			LOCK => $LOCK,
			g_index => $g_index
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
      elsif ($tier_type == 0 && $level == 1) {	#set up GFF3 output and fasta chunks
	 $level_status = 'setting up GFF3 output and fasta chunks';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(CTL_OPT
			out_dir
			build
			seq_id
			q_seq_obj
			safe_seq_id
			the_void
			q_def)
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $out_dir = $VARS->{out_dir};
	    my $build = $VARS->{build};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $the_void = $VARS->{the_void};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    
	    my $GFF3 = Dumper::GFF::GFFV3->new("$out_dir/$safe_seq_id.gff",
					       $build,
					       $the_void
					       );
	    
	    $GFF3->set_current_contig($seq_id, $q_seq_obj);
	    
	    #build chunks
	    my $fasta_chunker = new FastaChunker();
	    $fasta_chunker->parent_def($q_def);
	    $fasta_chunker->parent_seq($q_seq_obj);
	    $fasta_chunker->chunk_size($CTL_OPT{max_dna_len});
	    $fasta_chunker->min_size($CTL_OPT{max_dna_len}-1);
	    $fasta_chunker->load_chunks();
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = (GFF3 => $GFF3,
			fasta_chunker => $fasta_chunker
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
	    if($VARS->{CTL_OPT}{_no_mask}){
		$VARS->{fasta_chunker}->reset;
		$VARS->{m_index} = $VARS->{g_index};
		$VARS->{m_seq_obj} = $VARS->{q_seq_obj};
		$VARS->{masked_file} = $VARS->{unmasked_file};

		$next_level = 4;
	    }
	    elsif(-f $VARS->{the_void}."/query.masked.fasta" &&
		  -f $VARS->{the_void}."/query.masked.gff"){
		#make masked index
		my $seq_id = $VARS->{seq_id};
		my $masked_file = $VARS->{the_void}."/query.masked.fasta";
		my $m_index = GI::build_fasta_index($masked_file);
		my $m_seq_obj = $m_index->get_Seq_by_id($seq_id);

		#still no sequence? try rebuilding the index and try again
		if(!$m_seq_obj) {
		    sleep 5;
		    $m_index->reindex(0);
		    $m_seq_obj = $m_index->get_Seq_by_id($seq_id);
		    
		    if (not $m_seq_obj) {
			unlink $VARS->{the_void}."/query.masked.fasta";
			unlink $VARS->{the_void}."/query.masked.fasta.index";
			unlink $VARS->{the_void}."/query.masked.fasta.index.dir";
			unlink $VARS->{the_void}."/query.masked.fasta.index.pag";
			unlink $VARS->{the_void}."/query.masked.gff";
		    }
		}

		if($m_seq_obj){
		    #build masked chunks    
		    my $fasta_chunker = new FastaChunker();
		    $fasta_chunker->parent_def($VARS->{q_def}." masked");
		    $fasta_chunker->parent_seq($m_seq_obj);
		    $fasta_chunker->chunk_size($VARS->{CTL_OPT}{max_dna_len});
		    $fasta_chunker->min_size($VARS->{CTL_OPT}{max_dna_len}-1);
		    $fasta_chunker->load_chunks();

		    $VARS->{masked_file} = $masked_file;
		    $VARS->{m_index} = $m_index;
		    $VARS->{m_seq_obj} = $m_seq_obj;
		    $VARS->{fasta_chunker} = $fasta_chunker;

		    #GFF3 file with masking data
		    $VARS->{gff3_files} = [$VARS->{the_void}."/query.masked.gff"];
		    
		    $next_level = 4;
		}
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 0 && $level == 2) {     #build masking tiers
	 $level_status = 'builing masking tiers';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    #local masked GFF3
	    my $gff3_file = $VARS->{the_void}."/query.masked.gff";
            my $GFF3_m = Dumper::GFF::GFFV3->new($gff3_file,
						 '',
						 $VARS->{the_void}
						 );

            $GFF3_m->set_current_contig($VARS->{seq_id});

	    while(my $fchunk = $VARS->{fasta_chunker}->next_chunk){
	       my $order = $fchunk->number;
	       my $the_void = $VARS->{the_void};
               my $subvoid = ($main::old_struct) ? $the_void : "$the_void/".int($fchunk->start/1000000);
               mkdir($subvoid) unless(-d $subvoid);
	       my %args = (chunk        => $fchunk,
			   order        => $order,
			   the_void     => $VARS->{the_void},
			   subvoid      => $subvoid,
			   seq_id       => $VARS->{seq_id},
			   safe_seq_id  => $VARS->{safe_seq_id},
			   q_seq_obj    => $VARS->{q_seq_obj},
			   q_seq_length => $VARS->{q_seq_length},
			   dbfile       => $VARS->{dbfile},
			   GFF3_m       => $GFF3_m,
			   gff3_file    => $gff3_file,
			   LOG          => $VARS->{LOG},
			   DS_CTL       => $VARS->{DS_CTL},
			   CTL_OPT      => $VARS->{CTL_OPT}
			   );

	       my $tier_type = 1;
	       my $tier = new Process::MpiTiers(\%args, $self->rank, $self->{CHUNK_REF}, $tier_type);
	       push(@chunks, $tier); #really a tier
	    }

	    #empty chunker to save memory
            my $c_array = $VARS->{fasta_chunker}->{chunks};
            @$c_array = map {undef} @$c_array;
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";

	    #-------------------------CODE
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #-------------------------CODE

	    #------------------------RETURN
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'chunk'){
		  $VARS->{fasta_chunker}->replace($self->{RESULTS}->{chunk});
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    if(defined($VARS->{GFF3_m})){
		$VARS->{GFF3_m}->finalize;
		delete($VARS->{GFF3_m});
	    }
	    #-------------------------NEXT_LEVEL
	 }	    
      }	  
      ##
      ##TO TIER_TYPE 1
      ##
      elsif ($tier_type == 1 && $level == 0) {	#do repeat masking
	 $level_status = 'doing repeat masking';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{LOG} = Storable::dclone($VARS->{LOG}); #ensures independent log for each chunk
	    $VARS->{LOG}->set_child($VARS->{order});
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			q_seq_obj
			the_void
			subvoid
			safe_seq_id
			q_seq_length
			dbfile
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $LOG = $VARS->{LOG};
	    my $dbfile = $VARS->{dbfile};
	    my $chunk = $VARS->{chunk};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_seq_length = $VARS->{q_seq_length};

	    my $TMP = GI::get_global_temp;
	    if($CTL_OPT{go_gffdb} && GI::is_NFS_mount($dbfile) && !GI::is_NFS_mount($TMP)){
		$dbfile = GI::localize_file($dbfile);
	    }
	    my $GFF_DB = new GFFDB($dbfile) if($CTL_OPT{go_gffdb});

	    #-- repeatmask with gff3 input
	    my $rm_gff_keepers = [];
	    if ($CTL_OPT{go_gffdb}) {
	       $rm_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							    $q_seq_obj,
							    'repeat',
							    $q_seq_length
							   );
	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_gff_keepers)
		 if($CTL_OPT{organism_type} ne 'prokaryotic');
	    }

	    #-- repeatmask with RepeatMasker
	    my $rm_rb_keepers = []; #repeat masker RepBase
	    if ($CTL_OPT{model_org}) { #model organism repeats
	       my @models = split(/\,/, $CTL_OPT{model_org});
	       foreach my $mod (@models){
		   my $keepers = GI::repeatmask($chunk,
						$subvoid,
						$safe_seq_id,
						$mod,
						$CTL_OPT{RepeatMasker},
						'',
						$CTL_OPT{cpus},
						$LOG);
		   push(@$rm_rb_keepers, @$keepers);
	       }

	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_rb_keepers)
		   if($CTL_OPT{organism_type} ne 'prokaryotic');
	    }
	    my $rm_sp_keepers = []; #repeat masker species
	    if ($CTL_OPT{rmlib}) {  #species specific repeats;
	       foreach my $db (@{$VARS->{CTL_OPT}{_m_db}}){
		   my $keepers = GI::repeatmask($chunk,
						$subvoid,
						$safe_seq_id,
						$CTL_OPT{model_org},
						$CTL_OPT{RepeatMasker},
						$db,
						$CTL_OPT{cpus},
						$LOG);
		   push(@$rm_sp_keepers, @$keepers);
	       }

	       #mask the chunk
	       $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_sp_keepers)
		   if($CTL_OPT{organism_type} ne 'prokaryotic');

	    }
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = ( rm_gff_keepers => $rm_gff_keepers,
			 rm_rb_keepers => $rm_rb_keepers,
			 rm_sp_keepers => $rm_sp_keepers,
			 chunk => $chunk,
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
      elsif ($tier_type == 1 && $level == 1) {	#blastx repeat mask
	 $level_status = 'doing blastx repeats';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{rm_blastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_r_db}}){
	       my $blast_finished = GI::get_blast_finished_name($VARS->{chunk}->number,
								$db,
								$VARS->{the_void},
								$VARS->{safe_seq_id},
								'repeatrunner');
	       
	       next if($fin{$blast_finished});
	       
	       $db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
	       $VARS->{db} = $db;
	       $VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
	       $fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
	       my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	       push(@chunks, $chunk);
	    }
	    delete($VARS->{db});
	    delete($VARS->{LOG_FLAG});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(db
			chunk
			the_void
			subvoid
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT)
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};
	    
	    my $res_dir;
	    my $rm_blastx_keepers = [];
	    if ($db) {
	       GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
	       ($rm_blastx_keepers, $res_dir) = GI::repeatrunner_as_chunks($chunk,
									   $db,
									   $subvoid,
									   $safe_seq_id,
									   \%CTL_OPT,
									   $LOG,
									   $LOG_FLAG
									   );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( res_dir => $res_dir,
			 rm_blastx_keepers => $rm_blastx_keepers,
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
      elsif ($tier_type == 1 && $level == 2) {	#collect blastx repeatmask
	 $level_status = 'collecting blastx repeatmasking';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    #build chunks
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			res_dir
		        rm_blastx_keepers
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $res_dir = $VARS->{res_dir};
	    my $rm_blastx_keepers = $VARS->{rm_blastx_keepers};
	    my $LOG = $VARS->{LOG};

	    $rm_blastx_keepers = GI::combine_blast_report($chunk,
							  $rm_blastx_keepers,
							  $res_dir,
							  $LOG
							  );
	    
	    #mask the chunk
	    $chunk = repeat_mask_seq::mask_chunk($chunk, $rm_blastx_keepers)
	      if($CTL_OPT{organism_type} ne 'prokaryotic');

	    $res_dir = undef;
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( chunk => $chunk,
			 res_dir => $res_dir,
			 rm_blastx_keepers => $rm_blastx_keepers,
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
      elsif ($tier_type == 1 && $level == 3) {	#process all repeats
	 $level_status = 'processing all repeats';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			rm_gff_keepers
			rm_rb_keepers
			rm_sp_keepers
			rm_blastx_keepers
			GFF3_m
		       )
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $chunk = $VARS->{chunk};
	    my $rm_gff_keepers = $VARS->{rm_gff_keepers};
	    my $rm_rb_keepers = $VARS->{rm_rb_keepers};
	    my $rm_sp_keepers = $VARS->{rm_sp_keepers};
	    my $rm_blastx_keepers = $VARS->{rm_blastx_keepers};
	    my $GFF3_m = $VARS->{GFF3_m};

	    #-combine and cluster repeat hits for consensus
	    my $rm_keepers = repeat_mask_seq::process($rm_gff_keepers, 
						      $rm_rb_keepers, 
						      $rm_sp_keepers, 
						      $rm_blastx_keepers
						     );
      
	    #-add repeats to GFF3
	    my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);
	    $GFF3_m->add_repeat_hits($rm_keepers, $uid);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = (rm_keepers => $rm_keepers,
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
	    $next_level = undef;
	    #-------------------------NEXT_LEVEL
	 }
      }
      ##
      ##RETURN TO TIER_TYPE 0
      ##
      elsif ($tier_type == 0 && $level == 3) {	#prep masked sequence
	 $level_status = 'preparing masked sequence';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(the_void
			safe_seq_id
			seq_id
			q_def
                        fasta_chunker
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $fasta_chunker = $VARS->{fasta_chunker};
	    my $q_def = $VARS->{q_def};
	    my $seq_id = $VARS->{seq_id};

	    #write local contig fasta
	    my $seq;
	    my $masked_file = "$the_void/query.masked.fasta";

	    $VARS->{fasta_chunker}->reset;
	    while(my $fchunk = $VARS->{fasta_chunker}->next_chunk){
		$seq .= ${$fchunk->seq};
		$fchunk->seq(''); #clear memory
	    }

	    open (my $FAS, "> $masked_file");
	    print $FAS ${&Fasta::seq2fastaRef($seq_id, \$seq)};
	    close ($FAS);

	    $seq = undef; #clear memory

	    #make masked index and seq object
	    my $m_index = GI::build_fasta_index($masked_file);
	    my $m_seq_obj = $m_index->get_Seq_by_id($seq_id);

            #still no sequence? try rebuilding the index and try again
	    if(!$m_seq_obj) {
		for(my $i = 0; $i < 2 && !$m_seq_obj; $i++){
		    sleep 5;
		    print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
		    $m_index->reindex($i);
		    $m_seq_obj = $m_index->get_Seq_by_id($seq_id);
		}
		
		if (not $m_seq_obj) {
		    print STDERR "stop here: $seq_id\n";
		    confess "ERROR: Fasta index error\n";
		}
	    }

	    #build masked chunks
	    $fasta_chunker = new FastaChunker();
	    $fasta_chunker->parent_def($q_def." masked");
	    $fasta_chunker->parent_seq($m_seq_obj);
	    $fasta_chunker->chunk_size($CTL_OPT{max_dna_len});
	    $fasta_chunker->min_size($CTL_OPT{max_dna_len}-1);
	    $fasta_chunker->load_chunks();
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (masked_file => $masked_file,
			m_index => $m_index,
			m_seq_obj => $m_seq_obj,
			fasta_chunker => $fasta_chunker
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
      elsif ($tier_type == 0 && $level == 4){ #preparing evidence tiers
	 $level_status = 'preparing new fasta chunk tiers';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};

	    #make sure that the abinits have not already been processed
	    my $missing = 0;
	    my @files;
	    my $total_chunks = $VARS->{fasta_chunker}->total_chunks;
	    for(my $i = 0; $i < $total_chunks; $i++){
	       my $the_void = $VARS->{the_void};
	       my $safe_seq_id = $VARS->{safe_seq_id};
	       my $section_file = "$the_void/$safe_seq_id.$i.pred.raw.section";
	       if(! -f $section_file){
		   $missing = 1;
		   undef @files;
		   last;
	       }
	       else{
		   push(@files, $section_file);
	       }
	    }
	    
	    if(!$missing){
		push(@{$VARS->{section_files}}, @files);
	    }
	    else{
		my $chunk_count = int($VARS->{q_seq_length}/10000000 + 0.5) || 1;
		my $chunk_size  = int($VARS->{q_seq_length}/$chunk_count);
		
		my $mfasta_chunker = new FastaChunker();
		$mfasta_chunker->parent_def($VARS->{q_def}." masked");
		$mfasta_chunker->parent_seq($VARS->{m_seq_obj});
		$mfasta_chunker->chunk_size($chunk_size);
		$mfasta_chunker->min_size($chunk_size-1);
		$mfasta_chunker->flank(1000000);
		$mfasta_chunker->load_chunks();
		
		my $qfasta_chunker = new FastaChunker();
		$qfasta_chunker->parent_def($VARS->{q_def});
		$qfasta_chunker->parent_seq($VARS->{q_seq_obj});
		$qfasta_chunker->chunk_size($chunk_size);
		$qfasta_chunker->min_size($chunk_size-1);
		$qfasta_chunker->flank(1000000);
		$qfasta_chunker->load_chunks();

		#first tier is a pred only tier
		my %args = (the_void       => $VARS->{the_void},
			    fasta_chunker  => $VARS->{fasta_chunker},
			    q_def          => $VARS->{q_def},
			    q_seq_obj      => $VARS->{q_seq_obj},
			    m_seq_obj      => $VARS->{m_seq_obj},
			    q_seq_length   => $VARS->{q_seq_length},
			    seq_id         => $VARS->{seq_id},
			    safe_seq_id    => $VARS->{safe_seq_id},
			    LOG            => $VARS->{LOG},
			    CTL_OPT        => $VARS->{CTL_OPT},
			    qfasta_chunker => $qfasta_chunker,
			    mfasta_chunker => $mfasta_chunker,);
		
		my $tier_type = 2;
		my $tier = new Process::MpiTiers(\%args, $self->rank, $self->{CHUNK_REF}, $tier_type);
		push(@chunks, $tier); #really a tier
	    }

	    #all other tiers are alignment evidence
	    $VARS->{fasta_chunker}->reset;
	    while(my $fchunk = $VARS->{fasta_chunker}->next_chunk){
	       my $order = $fchunk->number;
	       my $the_void = $VARS->{the_void};
	       my $safe_seq_id = $VARS->{safe_seq_id};
	       my $section_file = "$the_void/$safe_seq_id.$order.raw.section";
	       my $junction_start_file = "$the_void/$safe_seq_id.".($order-1)."-$order.raw.section";
	       my $junction_end_file = "$the_void/$safe_seq_id.$order-".($order+1).".raw.section";
	       my $gff3_file = "$the_void/evidence_$order.gff";

	       if(-f $section_file &&
		  -f $gff3_file &&
		  (-f $junction_start_file || $fchunk->is_first) &&
		  (-f $junction_end_file || $fchunk->is_last)
		 ){
		  push(@{$VARS->{section_files}}, $section_file);
		  push(@{$VARS->{section_files}}, $junction_start_file) if(!$fchunk->is_first);
		  push(@{$VARS->{section_files}}, $junction_end_file) if(!$fchunk->is_last);
		  push(@{$VARS->{gff3_files}}, $gff3_file);
		  next;
	       }

	       my $GFF3_e = Dumper::GFF::GFFV3->new($gff3_file,
						    '',
						    $the_void
						    );
	       $GFF3_e->set_current_contig($VARS->{seq_id});

	       my $subvoid = ($main::old_struct) ? $the_void : "$the_void/".int($fchunk->start/1000000);
               mkdir($subvoid) unless(-d $subvoid);

	       my %args = (chunk        => $fchunk,
                           order        => $order,
			   the_void     => $VARS->{the_void},
			   subvoid      => $subvoid,
			   q_def        => $VARS->{q_def},
			   seq_id       => $VARS->{seq_id},
			   safe_seq_id  => $VARS->{safe_seq_id},
			   q_seq_obj    => $VARS->{q_seq_obj},
			   m_seq_obj    => $VARS->{m_seq_obj},
			   q_seq_length => $VARS->{q_seq_length},
			   dbfile       => $VARS->{dbfile},
			   GFF3_e       => $GFF3_e,
			   gff3_file    => $gff3_file,
			   LOG          => $VARS->{LOG},
			   DS_CTL       => $VARS->{DS_CTL},
			   CTL_OPT      => $VARS->{CTL_OPT}
			   );

	       my $tier_type = 3;
               my $tier = new Process::MpiTiers(\%args, $self->rank, $self->{CHUNK_REF}, $tier_type);
               push(@chunks, $tier); #really a tier
	    }
	    delete($VARS->{preds});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #-------------------------CODE

	    #------------------------RETURN
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'section_files' || $key eq 'holdover_files' || $key eq 'gff3_files'){
		  push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 2 && $level == 0) {	#abinits
	 $level_status = 'preparing ab-inits';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	     while(1){
		 $VARS->{mchunk} = $VARS->{mfasta_chunker}->next_chunk;
		 $VARS->{qchunk} = $VARS->{qfasta_chunker}->next_chunk;
		 last if(!$VARS->{mchunk});
		 my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
		 push(@chunks, $chunk);		 
	     }
	     delete($VARS->{qchunk});
	     delete($VARS->{mchunk});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(the_void
			safe_seq_id
			seq_id
			q_def
			mchunk
                        qchunk
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $q_def = $VARS->{q_def};
	    my $mchunk = $VARS->{mchunk};
	    my $qchunk = $VARS->{qchunk};
	    my $the_void = $VARS->{the_void};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
            my $LOG = $VARS->{LOG};

	    #==ab initio predictions here
	    #do masked predictions first
	    my $t_dir = GI::get_global_temp();
	    
	    my @abinits;
	    if(!$CTL_OPT{_no_mask} && @{$CTL_OPT{_run}}){
		my $sid = $safe_seq_id.".abinit_masked.".$mchunk->number();
		my $t_file = "$t_dir/$sid";
		$mchunk->write_file_w_flank($t_file) if(! -f $t_file);
		my @args = ($t_file, $the_void, \%CTL_OPT, $LOG);
		push(@abinits, @{GI::snap(@args)})     if(grep {/snap/} @{$CTL_OPT{_run}});
		push(@abinits, @{GI::augustus(@args)}) if(grep {/augustus/} @{$CTL_OPT{_run}});
		push(@abinits, @{GI::fgenesh(@args)})  if(grep {/fgenesh/} @{$CTL_OPT{_run}});
		unlink($t_file);
	    }

	    #genemark is never masked
	    if(grep {/genemark/} @{$CTL_OPT{_run}}){
		my $sid = $safe_seq_id.".abinit_nomask.".$mchunk->number();
		my $t_file = "$t_dir/$sid";
		$qchunk->write_file_w_flank($t_file) if(! -f $t_file);
		my @args = ($t_file, $the_void, \%CTL_OPT, $LOG);
		push(@abinits, @{GI::genemark(@args)});
		unlink($t_file) unless($CTL_OPT{unmask} || $CTL_OPT{_no_mask} || $CTL_OPT{trna});
	    }

	    #now do other unmasked predictions if requested
	    if($CTL_OPT{unmask} || $CTL_OPT{_no_mask}){
		my $sid = $safe_seq_id.".abinit_nomask.".$qchunk->number();
		my $t_file = "$t_dir/$sid";
		my @args = ($t_file, $the_void, \%CTL_OPT, $LOG);
		$qchunk->write_file_w_flank($t_file) if(! -f $t_file);
		push(@abinits, @{GI::snap(@args)})     if(grep {/snap/} @{$CTL_OPT{_run}});
		push(@abinits, @{GI::augustus(@args)}) if(grep {/augustus/} @{$CTL_OPT{_run}});
		push(@abinits, @{GI::fgenesh(@args)})  if(grep {/fgenesh/} @{$CTL_OPT{_run}});
		unlink($t_file) unless(grep{/trnascan/} @{$CTL_OPT{_run}});
	    }

	    #==ncRNA prediction here
	    #tRNAscan
            if(grep{/trnascan/} @{$CTL_OPT{_run}}){
                my $sid = $safe_seq_id.".abinit_nomask.".$mchunk->number();
		my $t_file = "$t_dir/$sid";
                $qchunk->write_file_w_flank($t_file) if(! -f $t_file);
		my @args = ($t_file, $the_void, \%CTL_OPT, $LOG);
		push(@abinits, @{GI::trnascan(@args)});
		unlink($t_file) unless(grep{/snoscan/} @{$CTL_OPT{_run}});
            }

	    #snoscan
            if(grep{/snoscan/} @{$CTL_OPT{_run}}){
                my $sid = $safe_seq_id.".abinit_nomask.".$mchunk->number();
                my $t_file = "$t_dir/$sid";
                $qchunk->write_file_w_flank($t_file) if(! -f $t_file);
                my @args = ($t_file, $the_void, \%CTL_OPT, $LOG);
                push(@abinits, @{GI::snoscan(@args)});
                unlink($t_file);
            }

	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = ( abinits => \@abinits
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
		if($key eq 'abinits'){
		    push(@{$VARS->{$key}}, @{$self->{RESULTS}{$key}});
		}
		else{
		    $VARS->{$key} = $self->{RESULTS}{$key};
		}
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 2 && $level == 1) {	#merge abinits
	 $level_status = 'gathering ab-init output files';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	     my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	     push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(the_void
			safe_seq_id
			abinits
			fasta_chunker
			qfasta_chunker
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $abinits = $VARS->{abinits};
	    my $fasta_chunker = $VARS->{fasta_chunker};	 
	    my $qfasta_chunker = $VARS->{qfasta_chunker};
            my $LOG = $VARS->{LOG};

	    #divide by type and hmm
	    my %data;
	    foreach my $d (@$abinits){
		$d->[0] =~ /(masked|nomask)\.\d+\.[^\.]+\.([^\.]+)$/;
		my $base = "$2\_$1|".$d->[1];
		push(@{$data{$base}}, $d);
	    }

	    #sort each set of files
	    my $crit = sub {
		my $file = shift;
		my ($n) = $file =~ /\.(\d+)\.[^\.]+\.[^\.]+$/;
		return $n;
	    };
	    @$_ = sort {&$crit($a->[0]) <=> &$crit($b->[0])} @$_ foreach (values %data);

	    #parse, merge, and section files
	    my @buf;
	    my %pkeepers;
	    my %mkeepers;
	    my @section_files;
	    $fasta_chunker->reset; #for sectioning files
	    my $c_count = $qfasta_chunker->total_chunks;
	    for(my $i = 0; $i < $c_count; $i++){
		my $qchunk0 = $qfasta_chunker->get_chunk($i);
		my $qchunk1 = $qfasta_chunker->get_chunk($i+1) if($i+1 < $c_count);

		#walk down all predictions (one section at a time)
		foreach my $k (keys %data){
		    #==parse the hits
		    #only first chunk
		    if($i == 0){
			my $preds = GI::parse_abinit_file(@{$data{$k}[$i]}, $qchunk0);
			my ($pp, $pm) = PhatHit_utils::separate_by_strand('query', $preds);
			@$pp = sort {$a->start <=> $b->start} @$pp;
			@$pm = sort {$a->start <=> $b->start} @$pm;
			$data{$k}[$i] = [$pp, $pm];
		    }
		    my $c0 = $data{$k}[$i];
		    
		    #==resolve chunk divide
		    if($qchunk1){
			#other chunks
			{ #block forces lexicals
			    my $preds = GI::parse_abinit_file(@{$data{$k}[$i+1]}, $qchunk1);
			    my ($pp, $pm) = PhatHit_utils::separate_by_strand('query', $preds);
			    @$pp = sort {$a->start <=> $b->start} @$pp;
			    @$pm = sort {$a->start <=> $b->start} @$pm;
			    $data{$k}[$i+1] = [$pp, $pm];
			}
			my $c1 = $data{$k}[$i+1];

			#go up to junction on plus strand
			my $plimit = $qchunk0->end + 1;
			while(my $p = shift @{$c0->[0]}){
			    if($p->end < $plimit){
				push(@{$pkeepers{$k}}, $p);
			    }
			    else{
				if($p->start < $plimit){
				    push(@{$pkeepers{$k}}, $p);
				    $plimit = $p->end+1;
				}
				undef(@{$c0->[0]}); #drop the rest
				last;
			    }
			}
			
			#go up to junction on minus
			my $mlimit = $qchunk0->end + 1;
			while(my $p = shift @{$c0->[1]}){
			    if($p->end < $mlimit){
				push(@{$mkeepers{$k}}, $p);
			    }
			    else{
				if($p->start < $mlimit){
				    push(@{$mkeepers{$k}}, $p);
				    $mlimit = $p->end+1;
				}
				undef @{$c0->[1]}; #drop the rest
				last;
			    }
			}
		    
			#trim off leading extra in neighbor plus
			while(my $p = shift @{$c1->[0]}){
			    if($p->start >= $plimit){
				unshift(@{$c1->[0]}, $p);
				last;
			    }
			}
			
			#trim off leading extra in neighbor minus
			while(my $p = shift @{$c1->[1]}){
			    if($p->start >= $mlimit){
				unshift(@{$c1->[1]}, $p);
				last;
			    }
			}
		    }
		    else{ #for last chunk
			push(@{$pkeepers{$k}}, @{$c0->[0]});
			push(@{$mkeepers{$k}}, @{$c0->[1]});
		    }
		}

		#==section the results overlapping first big chunk
		while(my $fchunk = shift @buf || $fasta_chunker->next_chunk){
		    my $E = $fchunk->end;
		    if($E <= $qchunk0->end){ #must not cross into neighboring chunk
			my @on_chunk;
			my @ncrna;
			foreach my $k (keys %pkeepers){
			    if($k =~ /^(trnascan|snoscan)/){ #separate ncRNA
				push(@ncrna, shift @{$pkeepers{$k}}) while(@{$pkeepers{$k}} &&
									   $pkeepers{$k}[0]->start <= $E);
				next;
			    }

			    push(@on_chunk, shift @{$pkeepers{$k}}) while(@{$pkeepers{$k}} &&
									  $pkeepers{$k}[0]->start <= $E);
			}
			foreach my $k (keys %mkeepers){
			    if($k =~ /^(trnascan|snoscan)/){ #separate ncRNA
				push(@ncrna, shift @{$mkeepers{$k}}) while(@{$mkeepers{$k}} &&
									   $mkeepers{$k}[0]->start <= $E);
				next;
			    }
			    push(@on_chunk, shift @{$mkeepers{$k}}) while(@{$mkeepers{$k}} &&
									  $mkeepers{$k}[0]->start <= $E);
			}
			my $order = $fchunk->number;
			my $section_file = "$the_void/$safe_seq_id.$order.pred.raw.section";
			if(! -f $section_file){
			    my %section = (preds_on_chunk => \@on_chunk,
					   ncrna_on_chunk => \@ncrna);
			    $LOG->add_entry("STARTED", $section_file, "");
			    store (\%section, $section_file);
			    $LOG->add_entry("FINISHED", $section_file, "");
			}
			push(@section_files, $section_file);
		    }
		    else{
			push(@buf, $fchunk);
			last;
		    }
		}
	    }
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = ( section_files => \@section_files
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
	     $next_level = undef;
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 0) {	#blastn
	 $level_status = 'doing blastn of ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{LOG} = Storable::dclone($VARS->{LOG}); #ensures independent log for each chunk
	    $VARS->{LOG}->set_child($VARS->{order});
	    $VARS->{blastn_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset

	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_e_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk}->number,
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'blastn');

		next if($fin{$blast_finished});

		$db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
		$VARS->{db} = $db;
		$VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
		$fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
		my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
		push(@chunks, $chunk);
	    }
	    delete($VARS->{db});
	    delete($VARS->{LOG_FLAG});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(db
			chunk
			the_void
			subvoid
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
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
								  $subvoid,
								  $safe_seq_id,
								  \%CTL_OPT,
								  $LOG,
								  $LOG_FLAG
								  );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( res_dir => $res_dir,
			 blastn_keepers => $blastn_keepers,
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
	    if(!$VARS->{res_dir} || !@{$VARS->{res_dir}}){
		$next_level = 4;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 1) {	#collect blastn
	 $level_status = 'collecting blastn reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			res_dir
			blastn_keepers
			m_seq_obj
			the_void
			subvoid
			safe_seq_id
			q_def
			q_seq_length
			LOG
			CTL_OPT
			edge_status)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $m_seq_obj = $VARS->{m_seq_obj};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};
	    my $edge_status = $VARS->{edge_status};

	    $blastn_keepers = GI::combine_blast_report($chunk,
						       $blastn_keepers,
						       $res_dir,
						       $LOG
						       );

	    #separate out hits too close to chunk divide to be run with exonerate
	    my $holdover_blastn = [];
	    ($blastn_keepers, $holdover_blastn) = GI::process_the_chunk_divide($chunk,
									       $CTL_OPT{'split_hit'},
									       $CTL_OPT{'pred_flank'},
									       0,
									       0,
									       0,
									       [$blastn_keepers]
									       );

	    #get start and end holdovers files
	    my @start;
	    my @end;
	    my @span;
	    foreach my $h (@$holdover_blastn){
	       if($h->{_holdover} == 1){
		   push(@start, $h);
	       }
	       elsif($h->{_holdover} == 2){
		   push(@end, $h);
	       }
	       elsif($h->{_holdover} == 3){
		   push(@span, $h);
	       }
	       else{
		   confess "ERROR:Holdover hits not labeled. This shouldn't happen?";
	       }
	    }

	    my $order = $chunk->number;
	    my $start_file = "$the_void/$safe_seq_id.$order.start.blastn.holdover";
	    my $end_file   = "$the_void/$safe_seq_id.$order.end.blastn.holdover";
	    my $start_neighbor = "$the_void/$safe_seq_id.".($order+1).".start.blastn.holdover";
	    my $end_neighbor   = "$the_void/$safe_seq_id.".($order-1).".end.blastn.holdover";

	    $holdover_blastn = [];
	    if(! $chunk->is_first){
	       my $lock1;
	       $lock1 = new File::NFSLock($start_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
	       $LOG->add_entry("STARTED", $start_file, "");
	       if(! -f $start_file){
		   store (\@start, $start_file);
	       }
	       $LOG->add_entry("FINISHED", $start_file, "");

	       my $start_junction = "$the_void/$safe_seq_id.".($order-1).".$order.junction.blastn.holdover";
	       my $lock2;
	       if(($lock2 = new File::NFSLock($start_junction, 'EX', 300, 40)) && (-f $end_neighbor) && $lock2->maintain(30)){
		   my $lock3;
		   $lock3 = new File::NFSLock($end_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		   my $neighbor = retrieve($end_neighbor);
		   unlink($start_file, $end_neighbor);

		   #merge and reblast
		   my $keepers = GI::merge_resolve_hits($m_seq_obj,
							$chunk->start,
							$q_def,
							$q_seq_length,
							$VARS->{CTL_OPT}{_e_db},
							\@start,
							$neighbor,
							$subvoid,
							\%CTL_OPT,
							'blastn',
							$LOG);
		   push(@$blastn_keepers, @$keepers);

		   $edge_status->{blastn_keepers}{start}++;
		   $edge_status->{exonerate_e_data}{start}++;
		   $lock1->unlock;
		   $lock3->unlock;
		   $lock2->unlock;
	       }
	       else{
		   $lock1->unlock;
	       }
	    }
	    if(! $chunk->is_last){
	       my $lock1;
	       $lock1 = new File::NFSLock($end_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
	       $LOG->add_entry("STARTED", $end_file, "");
	       if(! -f $end_file){
		  store (\@end, $end_file);
	       }
	       $LOG->add_entry("FINISHED", $end_file, "");

	       my $end_junction = "$the_void/$safe_seq_id.$order.".($order+1).".junction.blastn.holdover";
	       my $lock2;
	       if(($lock2 = new File::NFSLock($end_junction, 'EX', 300, 40)) && (-f $start_neighbor) && $lock2->maintain(30)){
		  my $lock3;
		  $lock3 = new File::NFSLock($start_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		  my $neighbor = retrieve($start_neighbor);
		  unlink($end_file, $start_neighbor);

		  #merge and reblast
		  my $keepers = GI::merge_resolve_hits($m_seq_obj,
						       $chunk->end,
						       $q_def,
						       $q_seq_length,
						       $VARS->{CTL_OPT}{_e_db},
						       \@end,
						       $neighbor,
						       $subvoid,
						       \%CTL_OPT,
						       'blastn',
						       $LOG);
		  push(@$blastn_keepers, @$keepers);
		  
		  $edge_status->{blastn_keepers}{end}++;
		  $edge_status->{exonerate_e_data}{end}++;
		  $lock1->unlock;
		  $lock3->unlock;
		  $lock2->unlock;
	       }
	       else{
		   $lock1->unlock;
	       }
	    }

	    #trim combined clusters
	    $blastn_keepers = cluster::shadow_cluster($CTL_OPT{depth_blastn}, $blastn_keepers);
	    $blastn_keepers = GI::flatten($blastn_keepers);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( blastn_keepers => $blastn_keepers,
			 res_dir => undef,
			 holdover_files => [$start_file, $end_file],
			 edge_status => $edge_status
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'holdover_files'){
		  push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 2) {	#exonerate ESTs
	 $level_status = 'polishig ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $seq_id = $VARS->{seq_id};
	    my $ex = $VARS->{CTL_OPT}{_to_exonerate}{$seq_id} || [];
	    my $ql = $VARS->{q_seq_length};
	    my $l = $VARS->{chunk}->length;

	    #set max chunks
	    my $altsize = int((@{$VARS->{blastn_keepers}} + @{$ex} * $l/$ql)/20);
	    my $size = ($altsize < $VARS->{CTL_OPT}->{_mpi_size}) ? $altsize : $VARS->{CTL_OPT}->{_mpi_size};
	    $size = 1 if(! $size);
	    $size = 10 if($size > 10);

	    my @data_sets;
	    for(my $i = 0; $i < @{$VARS->{blastn_keepers}}; $i++){
		my $j = $i % $size;
		push(@{$data_sets[$j]}, $VARS->{blastn_keepers}->[$i]);
	    }
	    for(my $i = 0; $i < @{$ex}; $i++){
		my $j = $i % $size;
		push(@{$data_sets[$j]}, $ex->[$i]);
	    }
	    for(my $i = 0; $i < @data_sets; $i++){
                $VARS->{dc}  = $data_sets[$i];
                $VARS->{id} = $i;
                my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
                push(@chunks, $chunk);
            }
	    delete($VARS->{dc});
	    delete($VARS->{id});

	    $VARS->{exonerate_e_data} = []; #reset
	    $VARS->{blastn_keepers} = []; #reset
	    $VARS->{exonerate_e_clusters} = []; #reset
	    $VARS->{blastn_clusters} = []; #reset
	    $VARS->{_clust_flag} = (@chunks > 10) ? 1 : 0;
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
		        dc
			the_void
			subvoid
			q_seq_length
			q_def
			q_seq_obj
			seq_id
			GFF3_e
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $dc = $VARS->{dc};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $seq_id = $VARS->{seq_id};
	    my $LOG = $VARS->{LOG};
	    my $GFF3_e = $VARS->{GFF3_e};

	    #-polish blastn hits with exonerate
	    my $exonerate_e_data = [];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
		#handle mising db in seq (IO error)
		if(! $q_seq_obj->{db}){
		    my $g_index = GI::build_fasta_index($CTL_OPT{_g_db});
		    $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		    for(my $i = 0; $i < 2 && !$q_seq_obj; $i++){
			sleep 5;
			print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
			$g_index->reindex($i);
			$q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		    }
		    if (! $q_seq_obj) {
			print STDERR "stop here: $seq_id\n";
			confess "ERROR: Fasta index error\n";
		    }
		}

		$exonerate_e_data = GI::polish_exonerate($chunk,
							 $q_seq_obj,
							 $q_seq_length,
							 $q_def,
							 $dc,
							 $CTL_OPT{_e_db},
							 $subvoid,
							 'e',
							 $CTL_OPT{exonerate},
							 $CTL_OPT{pcov_blastn},
							 $CTL_OPT{pid_blastn},
							 $CTL_OPT{en_score_limit},
							 $CTL_OPT{split_hit},
							 $CTL_OPT{min_intron},
							 $CTL_OPT{en_matrix},
							 $CTL_OPT{pred_flank},
							 $CTL_OPT{est_forward},
							 $LOG);
	    }

	    #-clean the blastn hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning blastn...\n" unless $main::quiet;
	    my $blastn_keepers = [grep {ref($_)} @{$dc}];

	    #Shatter hits. This is only for prokaryotic organisms.
	    #Flip strand of blastn where appropriate.
	    #This is done on blastn hits because exonerate is skipped.
	    #I shatter after processing the chunk divide to avoid weird
	    #complications from flipping on only one side of a split HSP
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$blastn_keepers  = PhatHit_utils::shatter_all_hits($blastn_keepers);

		#this checks the open reading frame and can flip the hit
		foreach my $phat_hit (@$blastn_keepers){
		    $phat_hit = PhatHit_utils::copy($phat_hit, 'both')
			if exonerate::splice_info::needs_to_be_revcomped($phat_hit, $q_seq_obj);
		}
	    }

	    $blastn_keepers = GI::clean_blast_hits($blastn_keepers,
						   $CTL_OPT{pcov_blastn},
						   $CTL_OPT{pid_blastn},
						   $CTL_OPT{eval_blastn},
						   1 #contiguity flag
						   );

	    my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);
	    $GFF3_e->add_phathits($blastn_keepers, $uid);
	    $GFF3_e->add_phathits($exonerate_e_data, $uid);

	    #blastn will be empty from this point on in the script if eukaryotic
	    my $blastn_clusters = []; #will stay empty for eukaryotes
	    my $exonerate_e_clusters = []; #will stay empty for prokaryotes
	    my $depth = ($CTL_OPT{organism_type} eq 'eukaryotic') ? 20 : 0;
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
	       $exonerate_e_clusters = cluster::clean_and_cluster($depth, $exonerate_e_data)
	    }
	    else{
	       $blastn_clusters = cluster::clean_and_cluster($depth, $blastn_keepers);
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_e_clusters => $exonerate_e_clusters,
			 blastn_clusters => $blastn_clusters,
			);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 3) {	#further cluster
	 $level_status = 'flattening EST clusters';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    delete($VARS->{_clust_flag});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(exonerate_e_clusters
		        blastn_clusters
			_clust_flag)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $exonerate_e_clusters = $VARS->{exonerate_e_clusters}; #array of overlapping clusters
	    my $blastn_clusters = $VARS->{blastn_clusters}; #array of overlapping clusters
	    my $flag = $VARS->{_clust_flag};	    

	    #further combine and cluster
	    my $depth = ($VARS->{CTL_OPT}->{organism_type} eq 'eukaryotic') ? 20 : 0;
	    if($flag){
	       $blastn_clusters = cluster::clean_and_cluster($depth, $blastn_clusters);
	       $exonerate_e_clusters = cluster::clean_and_cluster($depth, $exonerate_e_clusters);
	    }

	    my $blastn_keepers = GI::flatten($blastn_clusters);
	    my $exonerate_e_data = GI::flatten($exonerate_e_clusters);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_e_clusters => [], #clear memory
			 blastn_clusters => [], #clear memory
			 exonerate_e_data => $exonerate_e_data,
			 blastn_keepers => $blastn_keepers,
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
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }  
      elsif ($tier_type == 3 && $level == 4) {	#tblastx
	 $level_status = 'doing tblastx of alt-ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{tblastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
            foreach my $db (@{$VARS->{CTL_OPT}{_a_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk}->number,
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'tblastx');

                next if($fin{$blast_finished});

                $db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
                $VARS->{db} = $db;
                $VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
                $fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
                my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
                push(@chunks, $chunk);
            }
	    delete($VARS->{db});
	    delete($VARS->{LOG_FLAG});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(db
			chunk
			the_void
			subvoid
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $LOG = $VARS->{LOG};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    my $res_dir;
	    my $tblastx_keepers = [];
	    if ($db) {
	       GI::set_global_temp($CTL_OPT{_TMP}) if($CTL_OPT{_TMP});
	       ($tblastx_keepers, $res_dir) = GI::tblastx_as_chunks($chunk,
								    $db,
								    $subvoid,
								    $safe_seq_id,
								    \%CTL_OPT,
								    $LOG,
								    $LOG_FLAG
								    );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( res_dir => $res_dir,
			 tblastx_keepers => $tblastx_keepers,
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
	    if(!$VARS->{res_dir} || !@{$VARS->{res_dir}}){
		$next_level = 8;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 5) {	#collect tblastx
	 $level_status = 'collecting tblastx reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			res_dir
			tblastx_keepers
			m_seq_obj
			the_void
			subvoid
			safe_seq_id
			q_def
                        q_seq_length
			edge_status
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $m_seq_obj = $VARS->{m_seq_obj};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};
	    my $edge_status = $VARS->{edge_status};

	    $tblastx_keepers = GI::combine_blast_report($chunk,
							$tblastx_keepers,
							$res_dir,
							$LOG
							);

	    #separate out hits too close to chunk divide to be run with exonerate
	    my $holdover_tblastx = [];
	    ($tblastx_keepers, $holdover_tblastx) = GI::process_the_chunk_divide($chunk,
										 $CTL_OPT{'split_hit'},
										 $CTL_OPT{'pred_flank'},
										 0,
										 0,
										 0,
										 [$tblastx_keepers]
										 );

	    #get start and end holdovers files
	    my @start;
	    my @end;
	    my @span;
	    foreach my $h (@$holdover_tblastx){
	       if($h->{_holdover} == 1){
		   push(@start, $h);
	       }
	       elsif($h->{_holdover} == 2){
		   push(@end, $h);
	       }
	       elsif($h->{_holdover} == 3){
		   push(@span, $h);
	       }
	       else{
		   confess "ERROR:Holdover hits not labeled. This shouldn't happen?";
	       }
	    }

	    my $order = $chunk->number;
	    my $start_file = "$the_void/$safe_seq_id.$order.start.tblastx.holdover";
	    my $end_file   = "$the_void/$safe_seq_id.$order.end.tblastx.holdover";
	    my $start_neighbor = "$the_void/$safe_seq_id.".($order+1).".start.tblastx.holdover";
	    my $end_neighbor   = "$the_void/$safe_seq_id.".($order-1).".end.tblastx.holdover";

	    $holdover_tblastx = [];	    
	    if(! $chunk->is_first){
	       my $lock1;
	       $lock1 = new File::NFSLock($start_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
	       $LOG->add_entry("STARTED", $start_file, "");
	       if(! -f $start_file){
		  store (\@start, $start_file);
	       }
	       $LOG->add_entry("FINISHED", $start_file, "");

	       my $start_junction = "$the_void/$safe_seq_id.".($order-1).".$order.junction.tblastx.holdover";
	       my $lock2;
	       if(($lock2 = new File::NFSLock($start_junction, 'EX', 300, 40)) && (-f $end_neighbor) && $lock2->maintain(30)){
		  my $lock3;
		  $lock3 = new File::NFSLock($end_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		  my $neighbor = retrieve($end_neighbor);
		  unlink($start_file, $end_neighbor);

		  #merge and reblast
		  my $keepers = GI::merge_resolve_hits($m_seq_obj,
						       $chunk->start,
						       $q_def,
						       $q_seq_length,
						       $VARS->{CTL_OPT}{_a_db},
						       \@start,
						       $neighbor,
						       $subvoid,
						       \%CTL_OPT,
						       'tblastx',
						       $LOG);
		  push(@$tblastx_keepers, @$keepers);

		  $edge_status->{tblastx_keepers}{start}++;
		  $edge_status->{exonerate_a_data}{start}++;
		  $lock1->unlock;
		  $lock3->unlock;
		  $lock2->unlock;
	       }
	       else{
		   $lock1->unlock;
	       }
	    }
	    if(! $chunk->is_last){
	       my $lock1;
	       $lock1 = new File::NFSLock($end_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
	       $LOG->add_entry("STARTED", $end_file, "");
	       if(! -f $end_file){
		  store (\@end, $end_file);
	       }
	       $LOG->add_entry("FINISHED", $end_file, "");

	       my $end_junction = "$the_void/$safe_seq_id.$order.".($order+1).".junction.tblastx.holdover";
	       my $lock2;
	       if(($lock2 = new File::NFSLock($end_junction, 'EX', 300, 40)) && (-f $start_neighbor) && $lock2->maintain(30)){
		  my $lock3;
		  $lock3 = new File::NFSLock($start_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		  my $neighbor = retrieve($start_neighbor);
		  unlink($end_file, $start_neighbor);

		  #merge and reblast
		  my $keepers = GI::merge_resolve_hits($m_seq_obj,
						       $chunk->end,
						       $q_def,
						       $q_seq_length,
						       $VARS->{CTL_OPT}{_a_db},
						       \@end,
						       $neighbor,
						       $subvoid,
						       \%CTL_OPT,
						       'tblastx',
						       $LOG);
		  push(@$tblastx_keepers, @$keepers);

		  $edge_status->{tblastx_keepers}{end}++;
		  $edge_status->{exonerate_a_data}{end}++;
		  $lock1->unlock;
		  $lock3->unlock;
		  $lock2->unlock;
	       }
	       else{
		   $lock1->unlock;
	       }
	    }

	    #trim combined clusters
	    $tblastx_keepers = cluster::shadow_cluster($CTL_OPT{depth_tblastx}, $tblastx_keepers);
	    $tblastx_keepers = GI::flatten($tblastx_keepers);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( tblastx_keepers => $tblastx_keepers,
			 res_dir => undef,
                         holdover_files => [$start_file, $end_file],
			 edge_status => $edge_status
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'holdover_files'){
		  push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 6) {	#exonerate alt-ESTs
	 $level_status = 'polishing alt-ESTs';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    #set max chunks
	    my $altsize = int(@{$VARS->{tblastx_keepers}}/20);
	    my $size = ($altsize < $VARS->{CTL_OPT}->{_mpi_size}) ? $altsize : $VARS->{CTL_OPT}->{_mpi_size};
	    $size = 1 if(! $size);
	    $size = 10 if($size > 10);

	    my @data_sets;
	    for(my $i = 0; $i < @{$VARS->{tblastx_keepers}}; $i++){
	       my $j = $i % $size;
	       push(@{$data_sets[$j]}, $VARS->{tblastx_keepers}->[$i]);
	    }
	    for(my $i = 0; $i < @data_sets; $i++){
	       $VARS->{dc}  = $data_sets[$i];
	       $VARS->{id} = $i;
	       my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	       push(@chunks, $chunk);
	    }
	    delete($VARS->{dc});
	    delete($VARS->{id});

	    $VARS->{exonerate_a_data} = []; #reset
	    $VARS->{tblastx_keepers} = []; #reset
	    $VARS->{exonerate_a_clusters} = []; #reset
	    $VARS->{tblastx_clusters} = []; #reset
	    $VARS->{_clust_flag} = (@chunks > 10) ? 1 : 0;
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			dc
			the_void
			subvoid
			q_seq_length
			q_def
			q_seq_obj
			seq_id
			GFF3_e
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $dc = $VARS->{dc};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $seq_id = $VARS->{seq_id};
	    my $LOG = $VARS->{LOG};
	    my $GFF3_e = $VARS->{GFF3_e};

	    #-polish tblastx hits with exonerate
	    my $exonerate_a_data = [];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
	       #handle mising db in seq (IO error)
	       if(! $q_seq_obj->{db}){
		   my $g_index = GI::build_fasta_index($CTL_OPT{_g_db});
		   $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		   for(my $i = 0; $i < 2 && !$q_seq_obj; $i++){
		       sleep 5;
		       print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
		       $g_index->reindex($i);
		       $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		   }
		   if (! $q_seq_obj) {
		       print STDERR "stop here: $seq_id\n";
		       confess "ERROR: Fasta index error\n";
		   }
	       }

	       $exonerate_a_data = GI::polish_exonerate($chunk,
	                                                $q_seq_obj,
							$q_seq_length,
							$q_def,
							$dc,
							$CTL_OPT{_a_db},
							$subvoid,
							'a',
							$CTL_OPT{exonerate},
							$CTL_OPT{pcov_tblastx},
							$CTL_OPT{pid_tblastx},
							$CTL_OPT{en_score_limit},
	                                                $CTL_OPT{split_hit},
	                                                $CTL_OPT{min_intron},
							$CTL_OPT{en_matrix},
							$CTL_OPT{pred_flank},
							$CTL_OPT{est_forward},
							$LOG);
	    }
	    
	    #-clean the tblastx hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning tblastx...\n" unless $main::quiet;
	    my $tblastx_keepers = [grep {ref($_)} @{$dc}];

	    #Shatter hits. This is only for prokaryotic organisms.
	    #Flip strand of blastn where appropriate.
	    #This is done on blastn hits because exonerate is skipped.
	    #I shatter after processing the chunk divide to avoid weird
	    #complications from flipping on only one side of a split HSP
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$tblastx_keepers  = PhatHit_utils::shatter_all_hits($tblastx_keepers);

		#this checks the open reading frame and can flip the hit
		foreach my $phat_hit (@$tblastx_keepers){
		    $phat_hit = PhatHit_utils::copy($phat_hit, 'both')
			if exonerate::splice_info::needs_to_be_revcomped($phat_hit, $q_seq_obj);
		}
	    }

	    $tblastx_keepers = GI::clean_blast_hits($tblastx_keepers,
						    $CTL_OPT{pcov_tblastx},
						    $CTL_OPT{pid_tblastx},
						    $CTL_OPT{eval_tblastx},
						    1 #contiguity flag
						    );

	    my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);
	    $GFF3_e->add_phathits($tblastx_keepers, $uid);
	    $GFF3_e->add_phathits($exonerate_a_data, $uid);

	    my $tblastx_clusters = []; #will stay empty for eukaryotes
	    my $exonerate_a_clusters = []; #will stay empty for prokaryotes
	    my $depth = ($CTL_OPT{organism_type} eq 'eukaryotic') ? 20 : 0;
            if($CTL_OPT{organism_type} eq 'eukaryotic'){
               $exonerate_a_clusters = cluster::clean_and_cluster($depth, $exonerate_a_data)
            }
            else{
		$tblastx_clusters = cluster::clean_and_cluster($depth, $tblastx_keepers);
            }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_a_clusters => $exonerate_a_clusters,
			 tblastx_clusters => $tblastx_clusters
			);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
               push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 7) {	#further cluster and flatten
	 $level_status = 'flattening altEST clusters';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    delete($VARS->{_clust_flag});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(exonerate_a_clusters
		        tblastx_clusters
			_clust_flag)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $exonerate_a_clusters = $VARS->{exonerate_a_clusters}; #array of overlapping clusters
	    my $tblastx_clusters = $VARS->{tblastx_clusters}; #array of overlapping clusters
	    my $flag = $VARS->{_clust_flag};

	    #further combine and cluster
            my $depth = ($VARS->{CTL_OPT}->{organism_type} eq 'eukaryotic') ? 20 : 0;
	    if($flag){
	       $tblastx_clusters = cluster::clean_and_cluster($depth, $tblastx_clusters);
	       $exonerate_a_clusters = cluster::clean_and_cluster($depth, $exonerate_a_clusters);
	    }

	    my $tblastx_keepers = GI::flatten($tblastx_clusters);
	    my $exonerate_a_data = GI::flatten($exonerate_a_clusters);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_a_clusters => [], #clear memory
			 tblastx_clusters => [], #clear memory
			 exonerate_a_data => $exonerate_a_data,
			 tblastx_keepers => $tblastx_keepers,
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
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }  
      elsif ($tier_type == 3 && $level == 8) {	#blastx
	 $level_status = 'doing blastx of proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{blastx_keepers} = []; #reset
	    $VARS->{res_dir} = []; #reset
	    
	    #only create all chunks if not already finished
	    my %fin;
	    foreach my $db (@{$VARS->{CTL_OPT}{_p_db}}){
		my $blast_finished = GI::get_blast_finished_name($VARS->{chunk}->number,
								 $db,
								 $VARS->{the_void},
								 $VARS->{safe_seq_id},
								 'blastx');
		
                next if($fin{$blast_finished});

                $db =~ /\.mpi\.(\d+)\.(\d+)(\:.*)?$/;
                $VARS->{db} = $db;
                $VARS->{LOG_FLAG} = (!$2) ? 1 : 0;
                $fin{$blast_finished} = -e $blast_finished if($VARS->{LOG_FLAG});
                my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
                push(@chunks, $chunk);
            }
	    delete($VARS->{db});
	    delete($VARS->{LOG_FLAG});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(db
			chunk
			the_void
			subvoid
			safe_seq_id
			LOG
			LOG_FLAG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $db = $VARS->{db};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
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
								   $subvoid,
								   $safe_seq_id,
								   \%CTL_OPT,
								   $LOG,
								   $LOG_FLAG
								   );
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( res_dir => $res_dir,
			 blastx_keepers => $blastx_keepers,
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
	    if(!$VARS->{res_dir} || !@{$VARS->{res_dir}}){
		$next_level = 12;
	    }
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 9) {	#collect blastx
	 $level_status = 'collecting blastx reports';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			res_dir
			blastx_keepers
			m_seq_obj
			the_void
			subvoid
			safe_seq_id
			q_def
			q_seq_length
			edge_status
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $res_dir = $VARS->{res_dir};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $m_seq_obj = $VARS->{m_seq_obj};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $LOG = $VARS->{LOG};
	    my $chunk = $VARS->{chunk};
	    my $edge_status = $VARS->{edge_status};

	    $blastx_keepers = GI::combine_blast_report($chunk,
						       $blastx_keepers,
						       $res_dir,
						       $LOG
						       );
	    
	    #separate out hits too close to chunk divide to be run with exonerate
	    my $holdover_blastx = [];
	    ($blastx_keepers, $holdover_blastx) = GI::process_the_chunk_divide($chunk,
									       $CTL_OPT{'split_hit'},
									       $CTL_OPT{'pred_flank'},
									       0,
									       0,
									       0,
									       [$blastx_keepers]
									       );

	    #get start and end holdovers files
	    my @start;
	    my @end;
	    my @span;
	    foreach my $h (@$holdover_blastx){
	       if($h->{_holdover} == 1){
		   push(@start, $h);
	       }
	       elsif($h->{_holdover} == 2){
                   push(@end, $h);
               }
               elsif($h->{_holdover} == 3){
                   push(@span, $h);
               }
               else{
                   confess "ERROR:Holdover hits not labeled. This shouldn't happen?";
               }
	    }

	    my $order = $chunk->number;
	    my $start_file = "$the_void/$safe_seq_id.$order.start.blastx.holdover";
	    my $end_file   = "$the_void/$safe_seq_id.$order.end.blastx.holdover";
	    my $start_neighbor = "$the_void/$safe_seq_id.".($order+1).".start.blastx.holdover";
	    my $end_neighbor   = "$the_void/$safe_seq_id.".($order-1).".end.blastx.holdover";

	    $holdover_blastx = [];	    
	    if(! $chunk->is_first){
	       my $lock1;
	       $lock1 = new File::NFSLock($start_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
	       $LOG->add_entry("STARTED", $start_file, "");
	       if(! -f $start_file){
		   store (\@start, $start_file);
	       }
	       $LOG->add_entry("FINISHED", $start_file, "");
	       
	       my $start_junction = "$the_void/$safe_seq_id.".($order-1).".$order.junction.blastx.holdover";
	       my $lock2;
	       if(($lock2 = new File::NFSLock($start_junction, 'EX', 300, 40)) && (-f $end_neighbor) && $lock2->maintain(30)){
		   my $lock3;
		   $lock3 = new File::NFSLock($end_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		   my $neighbor = retrieve($end_neighbor);
		   unlink($start_file, $end_neighbor);

		   #merge and reblast
		   my $keepers = GI::merge_resolve_hits($m_seq_obj,
							$chunk->start,
							$q_def,
							$q_seq_length,
							$VARS->{CTL_OPT}{_p_db},
							\@start,
							$neighbor,
							$subvoid,
							\%CTL_OPT,
							'blastx',
							$LOG);		   
		   push(@$blastx_keepers, @$keepers);

		   $edge_status->{blastx_keepers}{start}++;
		   $edge_status->{exonerate_p_data}{start}++;
		   $lock1->unlock;
		   $lock3->unlock;
		   $lock2->unlock;
	       }
	       else{
		   $lock1->unlock;
	       }
	    }
	    if(! $chunk->is_last){
		my $lock1;
		$lock1 = new File::NFSLock($end_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
		$LOG->add_entry("STARTED", $end_file, "");
		if(! -f $end_file){
		    store (\@end, $end_file);
		}
		$LOG->add_entry("FINISHED", $end_file, "");

		my $end_junction = "$the_void/$safe_seq_id.$order.".($order+1).".junction.blastx.holdover";
		my $lock2;
		if(($lock2 = new File::NFSLock($end_junction, 'EX', 300, 40)) && (-f $start_neighbor) && $lock2->maintain(30)){
		    my $lock3;
		    $lock3 = new File::NFSLock($start_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		    my $neighbor = retrieve($start_neighbor);
		    unlink($end_file, $start_neighbor);

		    #merge and reblast
		    my $keepers = GI::merge_resolve_hits($m_seq_obj,
							 $chunk->end,
							 $q_def,
							 $q_seq_length,
							 $VARS->{CTL_OPT}{_p_db},
							 \@end,
							 $neighbor,
							 $subvoid,
							 \%CTL_OPT,
							 'blastx',
							 $LOG);
		    push(@$blastx_keepers, @$keepers);

		    $edge_status->{blastx_keepers}{end}++;
		    $edge_status->{exonerate_p_data}{end}++;
		    $lock1->unlock;
		    $lock3->unlock;
		    $lock2->unlock;
		}
		else{
		    $lock1->unlock;
		}
	    }

	    #trim combined clusters
	    $blastx_keepers = cluster::shadow_cluster($CTL_OPT{depth_blastx}, $blastx_keepers);
	    $blastx_keepers = GI::flatten($blastx_keepers);
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = ( blastx_keepers => $blastx_keepers,
			 res_dir => undef,
                         holdover_files => [$start_file, $end_file],
			 edge_status =>  $edge_status
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'holdover_files'){
		  push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 10) {	#exonerate proteins
	 $level_status = 'polishing proteins';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    #set max chunks
	    my $altsize = int(@{$VARS->{blastx_keepers}}/20);
	    my $size = ($altsize < $VARS->{CTL_OPT}->{_mpi_size}) ? $altsize : $VARS->{CTL_OPT}->{_mpi_size};
	    $size = 1 if(! $size);
	    $size = 10 if($size > 10);

	    my @data_sets;
	    for(my $i = 0; $i < @{$VARS->{blastx_keepers}}; $i++){
		my $j = $i % $size;
		push(@{$data_sets[$j]}, $VARS->{blastx_keepers}->[$i]);
	    }
	    for(my $i = 0; $i < @data_sets; $i++){
                $VARS->{dc}  = $data_sets[$i];
                $VARS->{id} = $i;
		my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
		push(@chunks, $chunk);
	    }
	    delete($VARS->{dc});
	    delete($VARS->{id});

	    $VARS->{exonerate_p_data} = []; #reset
	    $VARS->{blastx_keepers} = []; #reset
	    $VARS->{exonerate_p_clusters} = []; #reset
	    $VARS->{blastx_clusters} = []; #reset
	    $VARS->{_clust_flag} = (@chunks > 10) ? 1 : 0;
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			dc
			the_void
			subvoid
			q_seq_length
			q_def
			q_seq_obj
			seq_id
			GFF3_e
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $chunk = $VARS->{chunk};
	    my $dc = $VARS->{dc};
	    my $the_void = $VARS->{the_void};
	    my $subvoid = $VARS->{subvoid};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $q_def = $VARS->{q_def};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $seq_id = $VARS->{seq_id};
	    my $LOG = $VARS->{LOG};
	    my $GFF3_e = $VARS->{GFF3_e};

	    #-polish blastx hits with exonerate
	    my $exonerate_p_data = [];
	    if($CTL_OPT{organism_type} eq 'eukaryotic'){
		#handle mising db in seq (IO error)
		if(! $q_seq_obj->{db}){
		    my $g_index = GI::build_fasta_index($CTL_OPT{_g_db});
		    $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		    for(my $i = 0; $i < 2 && !$q_seq_obj; $i++){
			sleep 5;
			print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n" if($i);
			$g_index->reindex($i);
			$q_seq_obj = $g_index->get_Seq_by_id($seq_id);
		    }
		    if (! $q_seq_obj) {
			print STDERR "stop here: $seq_id\n";
			confess "ERROR: Fasta index error\n";
		    }
		}

		$exonerate_p_data = GI::polish_exonerate($chunk,
							 $q_seq_obj,
							 $q_seq_length,
							 $q_def,
							 $dc,
							 $CTL_OPT{_p_db},
							 $subvoid,
							 'p',
							 $CTL_OPT{exonerate},
							 $CTL_OPT{pcov_blastx},
							 $CTL_OPT{pid_blastx},
							 $CTL_OPT{ep_score_limit},
							 $CTL_OPT{split_hit},
							 $CTL_OPT{min_intron},
							 $CTL_OPT{ep_matrix},
							 $CTL_OPT{pred_flank},
							 $CTL_OPT{est_forward},
							 $LOG);
	    }

	    #-clean the blastx hits
	    #this must happen after exonerate (otherwise I filter out good hits)
	    print STDERR "cleaning blastx...\n" unless $main::quiet;
	    my $blastx_keepers = [grep {ref($_)} @{$dc}];

	    #Shatter hits. This is only for prokaryotic organisms.
	    #Flip strand of blastn where appropriate.
	    #This is done on blastn hits because exonerate is skipped.
	    #I shatter after processing the chunk divide to avoid weird
	    #complications from flipping on only one side of a split HSP
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$blastx_keepers  = PhatHit_utils::shatter_all_hits($blastx_keepers);
	    }

	    $blastx_keepers = GI::clean_blast_hits($blastx_keepers,
						   $CTL_OPT{pcov_blastx},
						   $CTL_OPT{pid_blastx},
						   $CTL_OPT{eval_blastx},
						   0 #contiguity flag
						   );

	    my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);
	    $GFF3_e->add_phathits($blastx_keepers, $uid);
	    $GFF3_e->add_phathits($exonerate_p_data, $uid);

	    #blastx will be empty from this point on in the script if eukaryotic
	    my $depth = ($CTL_OPT{organism_type} eq 'eukaryotic') ? 20 : 0;
	    my $blastx_clusters = cluster::clean_and_cluster($depth, $blastx_keepers);
	    my $exonerate_p_clusters = cluster::clean_and_cluster($depth, $exonerate_p_data);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_p_clusters => $exonerate_p_clusters,
			 blastx_clusters => $blastx_clusters,
			);
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
            while (my $key = each %{$self->{RESULTS}}) {
	       push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
            }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 11) {	#further cluster
	 $level_status = 'flattening protein clusters';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    delete($VARS->{_clust_flag});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(exonerate_p_clusters
		        blastx_clusters
			_clust_flag)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $exonerate_p_clusters = $VARS->{exonerate_p_clusters}; #array of overlapping clusters
	    my $blastx_clusters = $VARS->{blastx_clusters}; #array of overlapping clusters
	    my $flag = $VARS->{_clust_flag};

	    #further combine and cluster
            my $depth = ($VARS->{CTL_OPT}->{organism_type} eq 'eukaryotic') ? 20 : 0;
	    if($flag){
	       $blastx_clusters = cluster::clean_and_cluster($depth, $blastx_clusters);
	       $exonerate_p_clusters = cluster::clean_and_cluster($depth, $exonerate_p_clusters);
	    }

	    my $blastx_keepers = GI::flatten($blastx_clusters);
	    my $exonerate_p_data = GI::flatten($exonerate_p_clusters);
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( exonerate_p_clusters => [], #clear memory
			 blastx_clusters => [], #clear memory
			 exonerate_p_data => $exonerate_p_data,
			 blastx_keepers => $blastx_keepers,
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
	     $VARS->{LOG}->add_entry("###"); #indicate progress checkpoint
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 3 && $level == 12) {	#prepare section files
	 $level_status = 'prepare section files';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
                        the_void
			safe_seq_id
                        q_seq_obj
                        q_seq_length
			blastn_keepers
			blastx_keepers
			tblastx_keepers
			exonerate_e_data
			exonerate_a_data
			exonerate_p_data
			dbfile
			GFF3_e
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $chunk = $VARS->{chunk};
	    my $the_void = $VARS->{the_void};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $q_seq_length = $VARS->{q_seq_length};
	    my $blastn_keepers = $VARS->{blastn_keepers};
	    my $blastx_keepers = $VARS->{blastx_keepers};
	    my $tblastx_keepers = $VARS->{tblastx_keepers};
	    my $exonerate_e_data = $VARS->{exonerate_e_data};
	    my $exonerate_a_data = $VARS->{exonerate_a_data};
	    my $exonerate_p_data = $VARS->{exonerate_p_data};
	    my $dbfile = $VARS->{dbfile};
	    my $GFF3_e = $VARS->{GFF3_e};
	    my $LOG = $VARS->{LOG};
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};

	    my $TMP = GI::get_global_temp;
	    if($CTL_OPT{go_gffdb} && GI::is_NFS_mount($dbfile) && !GI::is_NFS_mount($TMP)){
		$dbfile = GI::localize_file($dbfile);
	    }
	    my $GFF_DB = new GFFDB($dbfile) if($CTL_OPT{go_gffdb});

	    #==GFF3 passthrough of evidence
	    my $prot_gff_keepers = [];
	    my $est_gff_keepers = [];
	    my $altest_gff_keepers = [];
	    my $model_gff_keepers = [];
	    my $pred_gff_keepers = [];
	    if ($CTL_OPT{go_gffdb}) {
	       print STDERR "Gathering GFF3 input into hits - chunk:".$chunk->number."\n"
		   unless($main::quiet);

	       my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);

	       #-protein evidence passthraough
	       $prot_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							      $q_seq_obj,
							      'protein',
							      $q_seq_length
							      );
	       $GFF3_e->add_phathits($prot_gff_keepers, $uid);

	       #-est evidence passthrough
	       $est_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							     $q_seq_obj,
							     'est',
                                                              $q_seq_length
							     );
	       $GFF3_e->add_phathits($est_gff_keepers, $uid);

	       #-altest evidence passthrough
	       $altest_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
								$q_seq_obj,
								'altest',
								$q_seq_length
								);
	       $GFF3_e->add_phathits($altest_gff_keepers, $uid);

	       #-gff gene annotation passthrough here
	       $model_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							       $q_seq_obj,
							       'model',
							       $q_seq_length
							       );
	       #$GFF3_e->add_phathits($model_gff_keepers, $uid);

	       #-pred passthrough
	       $pred_gff_keepers = $GFF_DB->phathits_on_chunk($chunk,
							      $q_seq_obj,
							      'pred',
                                                              $q_seq_length
							      );
	       #$GFF3_e->add_phathits($pred_gff_keepers, $uid);

	       #-other passthrough
	       my $lines = $GFF_DB->lines_for_chunk($chunk, 'other');
	       $GFF3_e->add_lines($lines);
	    }
	    $GFF3_e->finalize; #write_file

	    #trim evidence down to size if specified
	    my @sets = ($est_gff_keepers,
			$altest_gff_keepers,
			$prot_gff_keepers);
	    
	    #replaces actual values
	    my $depth = ($CTL_OPT{organism_type} eq 'eukaryotic') ? 20 : 0;
	    foreach my $set (@sets) {
		@$set = map {@$_} @{cluster::clean_and_cluster($depth, $set, 0, 1)};
	    }

	    my %section = (est_gff_keepers => $est_gff_keepers,
			   altest_gff_keepers => $altest_gff_keepers,
			   prot_gff_keepers => $prot_gff_keepers,
			   pred_gff_keepers => $pred_gff_keepers,
			   model_gff_keepers => $model_gff_keepers,
			   blastn_keepers => $blastn_keepers,
			   blastx_keepers => $blastx_keepers,
			   tblastx_keepers => $tblastx_keepers,
			   exonerate_e_data => $exonerate_e_data,
			   exonerate_a_data => $exonerate_a_data,
			   exonerate_p_data => $exonerate_p_data
			   );
	    
	    #separate out junction crossing hits
	    my %start_junction;
	    my %end_junction;
	    my ($aB, $aE) = ($chunk->offset +1, $chunk->offset + $chunk->length);
	    while(my $key = each %section){
	       my @keepers;
	       foreach my $h (@{$section{$key}}){
		  my ($bB, $bE) = PhatHit_utils::get_span_of_hit($h,'query');
		  ($bB, $bE) = ($bE, $bB) if($bB > $bE);
		  if($bB < $aB && ! $chunk->is_first) {
		     push(@{$start_junction{$key}}, $h);
		  }
		  elsif($bE > $aE && ! $chunk->is_last) {
		     push(@{$end_junction{$key}}, $h);
		  }
		  else{
		     push(@keepers, $h);
		  }
	       }
	       $section{$key} = \@keepers;
	    }

	    #make section files
	    my @all_files;
	    my $order = $chunk->number;

	    my $section_file = "$the_void/$safe_seq_id.$order.raw.section";
	    my $junction_start_file = "$the_void/$safe_seq_id.".($order-1)."-$order.raw.section";
	    my $junction_end_file = "$the_void/$safe_seq_id.$order-".($order+1).".raw.section";
	    $LOG->add_entry("STARTED", $section_file, "");
	    store (\%section, $section_file) unless(-f $section_file);
	    $LOG->add_entry("FINISHED", $section_file, "");
	    push(@all_files, $section_file);

	    #make combined junction files
	    my $start_file = "$the_void/$safe_seq_id.$order.start.section.holdover";
	    my $end_file = "$the_void/$safe_seq_id.$order.end.section.holdover";
	    my $start_neighbor = "$the_void/$safe_seq_id.".($order+1).".start.section.holdover";
	    my $end_neighbor   = "$the_void/$safe_seq_id.".($order-1).".end.section.holdover";

	    if(! $chunk->is_first && ! -f $junction_start_file){
		my $lock1;
		$lock1 = new File::NFSLock($start_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
		$LOG->add_entry("STARTED", $start_file, "");
		if(! -f $start_file){
		    store (\%start_junction, $start_file);
		}
		$LOG->add_entry("FINISHED", $start_file, "");

		my $lock2;
		if(($lock2 = new File::NFSLock($junction_start_file, 'EX', 300, 40)) && (-f $end_neighbor) && $lock2->maintain(30)){
		    my $lock3;
		    $lock3 = new File::NFSLock($end_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		    my $neighbor = retrieve($end_neighbor);
		    unlink($start_file, $end_neighbor);
		    while(my $key = each %$neighbor){
			push(@{$start_junction{$key}}, @{$neighbor->{$key}});
		    }
		    $LOG->add_entry("STARTED", $junction_start_file, "");
		    store(\%start_junction, $junction_start_file);
		    $LOG->add_entry("FINISHED", $junction_start_file, "");

		    push(@all_files, $junction_start_file) if(-f $junction_start_file);

		    $lock1->unlock;
		    $lock3->unlock;
		    $lock2->unlock;
		}
		else{
		    $lock1->unlock;
		}
	    }

	    if(! $chunk->is_last && ! -f $junction_end_file){
		my $lock1;
		$lock1 = new File::NFSLock($end_file, 'EX', 300, 40) while(! $lock1 || ! $lock1->maintain(30));
		$LOG->add_entry("STARTED", $end_file, "");
		if(! -f $end_file){
		    store (\%end_junction, $end_file);
		}
		$LOG->add_entry("FINISHED", $end_file, "");
		
		my $lock2;
		if(($lock2 = new File::NFSLock($junction_end_file, 'EX', 300, 40)) && (-f $start_neighbor) && $lock2->maintain(30)){
		    my $lock3;
		    $lock3 = new File::NFSLock($start_neighbor, 'EX', 300, 40) while(! $lock3 || ! $lock3->maintain(30));
		    my $neighbor = retrieve($start_neighbor);
		    unlink($end_file, $start_neighbor);
		    while(my $key = each %$neighbor){
			push(@{$end_junction{$key}}, @{$neighbor->{$key}});
		    }
		    $LOG->add_entry("STARTED", $junction_end_file, "");
		    store(\%end_junction, $junction_end_file);
		    $LOG->add_entry("FINISHED", $junction_end_file, "");

		    push(@all_files, $junction_end_file) if(-f $junction_end_file);

		    $lock1->unlock;
		    $lock3->unlock;
		    $lock2->unlock;
		}
		else{
		    $lock1->unlock;
		}
	    }
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = ( holdover_files => [$start_file, $end_file],
			 section_files => \@all_files
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
               if($key eq 'holdover_files'){
                  push(@{$VARS->{$key}}, @{$self->{RESULTS}->{$key}});
               }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    $next_level = undef;
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 0 && $level == 5) {	#process section files
	 $level_status = 'processing the chunk divide';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(fasta_chunker
			the_void
                        safe_seq_id
			m_seq_obj
			section_files
			LOG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $fasta_chunker = $VARS->{fasta_chunker};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $m_seq_obj = $VARS->{m_seq_obj};
	    my $section_files = $VARS->{section_files};
	    my $LOG = $VARS->{LOG};

	    #check if sections are already processed
	    $fasta_chunker->reset;
	    my @new_files;
	    while(my $chunk = $fasta_chunker->next_chunk){
		my $order = $chunk->number;
		my $section_file = "$the_void/$safe_seq_id.$order.final.section";
		if(-f $section_file){
		    push(@new_files, $section_file);
		}
		else{
		    undef @new_files;
		    last;
		}
	    }

	    #process data from each chunk linearly
	    $fasta_chunker->reset;
	    my $holdovers = {}; #holds heldover data on each chunk
	    while(my $chunk = $fasta_chunker->next_chunk){
	       last if(@new_files == $fasta_chunker->total_chunks); #already finished

	       #get evidence alignments from chunk section files
	       my $order = $chunk->number;
               my $s_file = "$the_void/$safe_seq_id.$order.raw.section";
               my $j_file = "$the_void/$safe_seq_id.$order-".($order+1).".raw.section";
	       ($s_file) = grep{$_ eq $s_file} @$section_files;
	       ($j_file) = grep{$_ eq $j_file} @$section_files;

	       #missing junction
	       if(!$chunk->is_last && (! $j_file || ! -f $j_file)){
		   unlink($j_file) if($j_file);
		   unlink($s_file) if($s_file);
		   confess "ERROR: Missing junction file for $order\-".($order+1)."\n";
	       }

	       my $section = retrieve($s_file);
	       my $junction = ($j_file) ? retrieve($j_file) : {};

	       #fix weird storable is ARRAY not HASH error
	       if (ref($junction) ne 'HASH'){
		   unlink($j_file);
		   confess STDERR "ERROR: Storable retrieval error.\n".
		                  "Must remove file before trying again.\n";
	       }

	       #merge the junction data onto the rest of the chunk section
	       while(my $key = each %$junction){
		  push(@{$section->{$key}}, @{$junction->{$key}});
	       }

	       #merge the junction data onto the rest of the chunk section
	       while(my $key = each %$holdovers){
		  push(@{$section->{$key}}, @{$holdovers->{$key}});
	       }

               my $sp_file = "$the_void/$safe_seq_id.$order.pred.raw.section";
	       ($sp_file) = grep{$_ eq $sp_file} @$section_files;

	       my $psection = retrieve($sp_file);

	       #merge the prediction data onto the rest of the chunk section
	       push(@{$section->{preds_on_chunk}}, @{$psection->{preds_on_chunk}})
		   if($psection->{preds_on_chunk});
	       push(@{$section->{ncrna_on_chunk}}, @{$psection->{ncrna_on_chunk}})
		   if($psection->{ncrna_on_chunk});

	       #keys to grab out of $section hash
	       my @keys = qw(blastn_keepers
			     blastx_keepers
			     tblastx_keepers
			     preds_on_chunk
			     ncrna_on_chunk
			     est_gff_keepers
			     altest_gff_keepers
			     prot_gff_keepers
			     pred_gff_keepers
			     model_gff_keepers
			     exonerate_e_data
			     exonerate_a_data
			     exonerate_p_data);

	       #==PROCESS HITS CLOSE TO CODE DIVISIONS
	       #holdover hits that are too close to the divide for review with next chunk
	       if (not $chunk->is_last) { #if not last chunk
		  (@{$section}{@keys},
		   @{$holdovers}{@keys}) = GI::process_the_chunk_divide($chunk,
									$CTL_OPT{'split_hit'},
									$CTL_OPT{'pred_flank'},
									1,
									1,
									1,
									[@{$section}{@keys}]
									);
	       }

	       #write final section file
	       my $section_file = "$the_void/$safe_seq_id.$order.final.section";
	       $LOG->add_entry("STARTED", $section_file, "");
	       store($section, $section_file) unless(-f $section_file);
	       $LOG->add_entry("FINISHED", $section_file, "");
	       push (@new_files, $section_file);
	    }
	    $section_files = \@new_files;
	    #-------------------------CODE
	    
	    #------------------------RETURN
	    %results = ( section_files => $section_files,
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
      elsif ($tier_type == 0 && $level == 6) {     #build annotation tiers
	 $level_status = 'builing annotation tiers';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    unlink(@{$VARS->{holdover_files}}) if($VARS->{holdover_files});
	    delete($VARS->{holdover_files});

	    return [] if(!@{$VARS->{CTL_OPT}{_run}} && 
			 !@{$VARS->{CTL_OPT}{_predictor}});

	    my $section_files = $VARS->{section_files};
	    $VARS->{fasta_chunker}->reset;
	    while(my $fchunk = $VARS->{fasta_chunker}->next_chunk){
	       my $order = $fchunk->number;
               my $the_void = $VARS->{the_void};
	       my $subvoid = ($main::old_struct) ? $the_void : "$the_void/".int($fchunk->start/1000000);
               mkdir($subvoid) unless(-d $subvoid);
	       my ($file) = grep {/\.$order\.final\.section$/} @$section_files;
	       my %args = (section_file       => $file,
			   chunk              => $fchunk,
			   order              => $order,
			   the_void           => $VARS->{the_void},
			   subvoid            => $subvoid,
			   safe_seq_id        => $VARS->{safe_seq_id},
			   seq_id             => $VARS->{seq_id},
			   q_def              => $VARS->{q_def},
			   out_dir            => $VARS->{out_dir},
			   q_seq_obj          => $VARS->{q_seq_obj},
			   m_seq_obj          => $VARS->{m_seq_obj},
			   GFF3               => $VARS->{GFF3},
			   LOG                => $VARS->{LOG},
			   DS_CTL             => $VARS->{DS_CTL},
			   CTL_OPT            => $VARS->{CTL_OPT});

	       my $tier_type = 4;
	       my $tier = new Process::MpiTiers(\%args, $self->rank, $self->{CHUNK_REF}, $tier_type);
	       push(@chunks, $tier); #really a tier
	    }
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    confess "ERROR: Logic error in tier_type:$tier_type, level:$level, flag:$flag.\n";
	    #-------------------------CODE

	    #------------------------RETURN
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'p_fastas' || $key eq 't_fastas' || $key eq 'n_fastas'){
		  while(my $key2 = each %{$self->{RESULTS}->{$key}}){
		     $VARS->{$key}->{$key2} .= $self->{RESULTS}->{$key}->{$key2};
		  }
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 4 && $level == 0) {	#prep hint clusters
	 $level_status = 'preparing evidence clusters for annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{LOG} = Storable::dclone($VARS->{LOG}); #ensures independent log for each chunk
	    $VARS->{LOG}->set_child($VARS->{order});
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw(section_file
			seq_id
			safe_seq_id
			q_seq_obj
                        chunk
                        the_void
			CTL_OPT
			LOG)
		     );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $chunk = $VARS->{chunk};
	    my $section_file = $VARS->{section_file};
	    my $LOG = $VARS->{LOG};

	    my $section = Storable::retrieve($section_file);
	    my $tblastx_keepers    = $section->{tblastx_keepers};
	    my $blastx_keepers     = $section->{blastx_keepers};
	    my $blastn_keepers     = $section->{blastn_keepers};
	    my $exonerate_e_data   = $section->{exonerate_e_data};
	    my $exonerate_a_data   = $section->{exonerate_a_data};
	    my $exonerate_p_data   = $section->{exonerate_p_data};
	    my $preds_on_chunk     = $section->{preds_on_chunk};
	    my $ncrna_on_chunk     = $section->{ncrna_on_chunk};
	    my $est_gff_keepers    = $section->{est_gff_keepers};
	    my $altest_gff_keepers = $section->{altest_gff_keepers};
	    my $prot_gff_keepers   = $section->{prot_gff_keepers};
	    my $pred_gff_keepers   = $section->{pred_gff_keepers};
	    my $ncrna_gff_keepers  = $section->{ncrna_gff_keepers};
	    my $model_gff_keepers  = $section->{model_gff_keepers};

	    #combine final data sets
	    print STDERR "Preparing evidence for hint based annotation\n" unless($main::quiet);
	    my $final_est = GI::combine($exonerate_e_data,
					$est_gff_keepers);

	    #add unpolished alignments (for prokaryotes)
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$final_est = GI::combine($final_est, $blastn_keepers);
	    }
	    else{ #remove unpolisjed alignments (from gff3)
		@$final_est = grep {$_->algorithm !~ /^blastn/} @$final_est;
	    }

	    my $final_altest = GI::combine($exonerate_a_data,
					   $altest_gff_keepers);

	    #add unpolished alignments (for prokaryotes)
	    if($CTL_OPT{organism_type} eq 'prokaryotic'){
		$final_altest = GI::combine($final_altest, $tblastx_keepers);
	    }
	    else{ #remove unpolished alignments (from gff3)
		@$final_altest = grep {$_->algorithm !~ /^tblastx/} @$final_altest;
	    }

	    my $final_prot = GI::combine($blastx_keepers,
					 $exonerate_p_data,
					 $prot_gff_keepers);

	    my $final_pred = GI::combine($preds_on_chunk,
					 $pred_gff_keepers);

	    my $final_ncrna = GI::combine($ncrna_on_chunk,
					  $ncrna_gff_keepers);

	    #run evm now that the evidence is aligned and the predictors have run                                                       
            my $evm_preds = [];#makes an empty array ref                                                                                
	    if(grep{/evm/} @{$CTL_OPT{_run}}){
		my $t_dir = GI::get_global_temp();
		my $sid = $safe_seq_id.".abinit_nomask.".$chunk->number();                                                                
		my $t_file = "$t_dir/$sid";
		my $seq = $q_seq_obj->seq();
		$seq = Fasta::toFastaRef('>'.$seq_id, \$seq); #over writes $seq to save memory                                       
		FastaFile::writeFile($seq, $t_file); #takes the fasta ref and writes it out as a wraped fasta file                      
		
		$evm_preds = GI::evm($t_file,
				     $the_void,
				     \%CTL_OPT,
				     $LOG,
				     $final_prot,
				     $final_est,
				     $final_altest,
				     $final_pred,
				     $q_seq_obj,
				     $seq_id);
		
		push(@$final_pred, @$evm_preds);# $evm_preds onto $final_pred dereference both of them                                 
	    }
	    
	    #group evidence for annotation
	    my $all_data = maker::auto_annotator::prep_hits($final_prot,
							    $final_est,
							    $final_altest,
							    $final_pred,
							    $final_ncrna,
							    $model_gff_keepers,
							    $q_seq_obj,
							    $CTL_OPT{single_exon},
							    $CTL_OPT{single_length},
							    $CTL_OPT{pred_flank},
							    $CTL_OPT{organism_type},
							    $CTL_OPT{est_forward},
							    $CTL_OPT{correct_est_fusion});
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = (all_data => $all_data,
			evm_preds => $evm_preds #check
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
      elsif ($tier_type == 4 && $level == 1) {	#annotate transcripts
	 $level_status = 'annotating transcripts';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    $VARS->{trans} = {}; #reset

            my $altsize = int(@{$VARS->{all_data}}/20);
	    my $size = ($altsize < $VARS->{CTL_OPT}->{_mpi_size}) ? $altsize : $VARS->{CTL_OPT}->{_mpi_size};
            $size = 1 if(! $size);
            $size = 10 if($size > 10);
	    $size = 1;

	    my @data_sets;
	    for(my $i = 0; $i < @{$VARS->{all_data}}; $i++){
		my $j = $i % $size;
		push(@{$data_sets[$j]}, $VARS->{all_data}->[$i]);
	    }
	    foreach my $dc (@data_sets){
		$VARS->{dc} = $dc;
		$VARS->{LOG_FLAG} = (!@chunks) ? 1 : 0;
		my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
		push(@chunks, $chunk);
            }
	    delete($VARS->{dc});
	    delete($VARS->{LOG_FLAG});
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(the_void
			q_def
			seq_id
			q_seq_obj
			m_seq_obj
			dc
			LOG
			LOG_FLAG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $m_seq_obj = $VARS->{m_seq_obj};
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
	    my $trans = maker::auto_annotator::annotate_trans($q_seq_obj,
							      $m_seq_obj,
							      $q_def,
							      $seq_id,
							      $dc,
							      $the_void,
							      \%CTL_OPT,
							      $LOG
							      );
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( trans => $trans
		       );
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
      elsif ($tier_type == 4 && $level == 2) {	#grouping transcripts into genes
	 $level_status = 'clustering transcripts into genes for annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw(trans
			all_data
			q_seq_obj
			seq_id
			chunk
			the_void
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $trans = $VARS->{trans};
	    my $all_data = $VARS->{all_data};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $seq_id = $VARS->{seq_id};
	    my $chunk = $VARS->{chunk};
	    my $the_void = $VARS->{the_void};

	    #group transcriipts into genes
	    print STDERR "Processing transcripts into genes\n" unless($main::quiet);
	    my $annotations = maker::auto_annotator::annotate_genes($trans,
								    $all_data,
								    $q_seq_obj,
								    $seq_id,
								    $chunk->number(),
								    $the_void,
								    \%CTL_OPT
								    );
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = ( annotations   => $annotations
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
      elsif ($tier_type == 4 && $level == 3) {	#adding quality control statistics
	 $level_status = 'adding statistics to annotations';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $altsize = 0;
	    while (my $key = each %{$VARS->{annotations}}){
		$altsize += @{$VARS->{annotations}->{$key}};
	    }
	    
	    $altsize = int($altsize/20);
	    my $size = ($altsize < $VARS->{CTL_OPT}->{_mpi_size}) ? $altsize : $VARS->{CTL_OPT}->{_mpi_size};
	    $size = 1 if(! $size);
	    $size = 10 if($size > 10);
	    $size = 1;

	    my @data_sets;
	    my $i = 0;
	    while (my $key = each %{$VARS->{annotations}}){ 
		foreach my $an (@{$VARS->{annotations}->{$key}}){
		    my $j = $i % $size;
		    push(@{$data_sets[$j]{$key}}, $an);
		    $i++;
		}
	    }
	     
	    foreach my $an (@data_sets){
		$VARS->{an} = $an;
		$VARS->{LOG_FLAG} = (!@chunks) ? 1 : 0;
		my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
		push(@chunks, $chunk);
	    }

	    delete($VARS->{an});
	    delete($VARS->{LOG_FLAG});
	    $VARS->{annotations} = {}; #reset
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw(q_seq_obj
			an
			LOG_FLAG
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $an = $VARS->{an};
	    my $q_seq_obj = $VARS->{q_seq_obj};
	    my $LOG_FLAG = $VARS->{LOG_FLAG};

	    #adds AED and other quality control statistics (can change names)
	    print STDERR "Calculating annotation quality statistics\n" unless($main::quiet || !$LOG_FLAG);
	    my $annotations = maker::auto_annotator::annotate_stats($an,
								    $q_seq_obj,
								    \%CTL_OPT);
	    #-------------------------CODE
	 
	    #------------------------RETURN
	    %results = ( annotations => $annotations
		       );
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
      elsif ($tier_type == 4 && $level == 4) {	#deciding on final annotations
	 $level_status = 'choosing best annotation set';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
      	 elsif ($flag eq 'init') {
            #------------------------ARGS_IN
	    @args = (qw(annotations
			out_dir
			CTL_OPT)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $annotations = $VARS->{annotations};
	    my $out_dir = $VARS->{out_dir};

	    #get best annotations
	    print STDERR "Choosing best annotations\n" unless($main::quiet);
	    my $maker_anno = maker::auto_annotator::best_annotations($annotations,
								     \%CTL_OPT
								    );
	    
	    #get best non-overlapping ab-inits
	    my $non_over = maker::auto_annotator::get_non_overlaping_abinits($maker_anno,
									     $annotations,
									     \%CTL_OPT
									     );

	    #get non-coding annotations
	    my $non_coding = [];
	    my @nc_keys = grep {/^ncrna_|_ncrna$/} keys %$annotations;
	    foreach my $k (@nc_keys){
		push(@$non_coding, @{$annotations->{$k}});
	    }
	    
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
		#						$q_seq_obj,
		#						$out_dir,
		#						$the_void,
		#						\%CTL_OPT
		#						);
	    }

	    #get AED scored preds for GFF3
	    my @scored_preds;
	    while(my $key = each %$annotations){
		next unless($key =~ /^(pred|model)_gff/ || $key =~ /_abinit$/);
		
		foreach my $g (@{$annotations->{$key}}){
		    if($key =~ /^model_gff/){
			my @models = map {$_->{hit}} @{$g->{t_structs}};
			push(@scored_preds, @models);
		    }
		    else{
			my @p_bases = map {$_->{p_base}} @{$g->{t_structs}};
			push(@scored_preds, @p_bases);
		    }
		}
	    }
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( maker_anno => $maker_anno,
			 non_over => $non_over,
			 non_coding => $non_coding,
			 scored_preds => \@scored_preds
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
      elsif ($tier_type == 4 && $level == 5) {	#local output
	 $level_status = 'processing chunk output';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------Chunker
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(chunk
			maker_anno
			non_over
                        non_coding
			annotations
			scored_preds
			GFF3)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my $chunk        = $VARS->{chunk};
	    my $maker_anno   = $VARS->{maker_anno};
	    my $non_over     = $VARS->{non_over};
	    my $non_coding   = $VARS->{non_coding};
	    my $annotations  = $VARS->{annotations};
	    my $scored_preds = $VARS->{scored_preds};
	    my $GFF3         = $VARS->{GFF3};

	    #==OUTPUT DATA HERE      
	    #--- GFF3
	    my $uid = join('.', $tier_type, $level, $self->number, $chunk->number);
	    $GFF3->add_genes($maker_anno);
	    $GFF3->add_ncgenes($non_coding);
	    $GFF3->add_phathits($scored_preds, $uid);
	    $GFF3->resolved_flag if (not $chunk->is_last); #adds ### between contigs
            
	    #--- building fastas for annotations (grows with iteration)
	    my $p_fastas = {};
	    my $t_fastas = {};
	    my $n_fastas = {};
	    GI::maker_p_and_t_fastas($maker_anno,
				     $non_coding,
				     $non_over,
				     $annotations,
				     $p_fastas,
				     $t_fastas,
				     $n_fastas,
				   );
	    #-------------------------CODE

	    #------------------------RETURN
	    %results = ( p_fastas => $p_fastas,
			 t_fastas => $t_fastas,
			 n_fastas => $n_fastas,
			 maker_anno         => [], #clear memory
			 blastn_keepers     => [], #clear memory
			 tblastx_keepers    => [], #clear memory
			 blastx_keepers     => [], #clear memory
			 exonerate_e_data   => [], #clear memory
			 exonerate_a_data   => [], #clear memory
			 exonerate_p_data   => [], #clear memory
			 est_gff_keepers    => [], #clear memory
			 altest_gff_keepers => [], #clear memory
			 prot_gff_keepers   => [], #clear memory
			 scored_preds       => [], #clear memory
		       );
	    #------------------------RETURN
	 }
	 elsif ($flag eq 'result') {
	    #-------------------------RESULT
	    while (my $key = each %{$self->{RESULTS}}) {
	       if($key eq 'p_fastas' || $key eq 't_fastas' || $key eq 'n_fastas'){
		  while(my $key2 = each %{$self->{RESULTS}->{$key}}){
		     $VARS->{$key}->{$key2} .= $self->{RESULTS}->{$key}->{$key2};
		  }
	       }
	       else{
		  $VARS->{$key} = $self->{RESULTS}->{$key};
	       }
	    }
	    #-------------------------RESULT
	 }
	 elsif ($flag eq 'flow') {
	    #-------------------------NEXT_LEVEL
	    $next_level = undef;
	    #-------------------------NEXT_LEVEL
	 }
      }
      elsif ($tier_type == 0 && $level == 7) {	#global output
	 $level_status = 'processing contig output';
	 if ($flag eq 'load') {
	    #-------------------------CHUNKER
	    my $chunk = new Process::MpiChunk($VARS, $level, $tier_type);
	    push(@chunks, $chunk);
	    #-------------------------CHUNKER
	 }
	 elsif ($flag eq 'init') {
	    #------------------------ARGS_IN
	    @args = (qw(the_void
			out_dir
			seq_id
			safe_seq_id
			p_fastas
			t_fastas
			n_fastas
			GFF3
			gff3_files
			DS_CTL
			CTL_OPT
			LOG
			LOCK)
		    );
	    #------------------------ARGS_IN
	 }
	 elsif ($flag eq 'run') {
	    print STDERR "$level_status\n";
	    #-------------------------CODE
	    my %CTL_OPT = %{$VARS->{CTL_OPT}};
	    my $the_void = $VARS->{the_void};
	    my $out_dir = $VARS->{out_dir};
	    my $seq_id = $VARS->{seq_id};
	    my $safe_seq_id = $VARS->{safe_seq_id};
	    my $p_fastas = $VARS->{p_fastas};
	    my $t_fastas = $VARS->{t_fastas};
	    my $n_fastas = $VARS->{n_fastas};
	    my $GFF3 = $VARS->{GFF3};
	    my $gff3_files = $VARS->{gff3_files};
	    my $DS_CTL = $VARS->{DS_CTL};
	    my $LOG = $VARS->{LOG};
	    my $LOCK = $VARS->{LOCK};

	    #--- write fastas for ab-initio predictions
	    GI::write_p_and_t_fastas($p_fastas, $t_fastas, $n_fastas, $safe_seq_id, $out_dir);
	    
	    #--- write GFF3 file
	    $GFF3->merge($gff3_files);
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
   confess "FATAL: \'$flag\' is not a valid flag in MpiChunk _go!!\n";
}
#--------------------------------------------------------------
#gets called by MpiTiers object following termination
#called from {CHUNK_REF}
sub _on_termination {
   my $self = shift;
   my $tier = shift;
   my $tier_type = $tier->{TIER_TYPE};

   #handle case of calling as function rather than method
   if ($self ne "Process::MpiChunk" && ref($self) ne "Process::MpiChunk") {
      $tier = $self;
      $self = new Process::MpiChunk();
   }

   return if($tier->failed);
   return if($tier->{VARS}{c_flag} <= 0);

   #only reach this point if termination is due to success
   $self->{RESULTS} = {};
   my $LOG = $tier->{VARS}{LOG};
   my $DS_CTL = $tier->{VARS}{DS_CTL};

   if($tier_type == 0){
      $tier->{VARS}{DS_CTL}->add_entry($tier->{VARS}{seq_id},
				       $tier->{VARS}{out_dir},
				       "FINISHED");
      $tier->{VARS}{LOCK}->unlock; #releases locks on the log file
      
      $tier->{RESULTS} = {};
   }
   elsif($tier_type == 1){
      $tier->{RESULTS}->{chunk} = $tier->{VARS}{chunk};
      $tier->{RESULTS}->{GFF3_m} = $tier->{VARS}{GFF3_m};
      $tier->{RESULTS}->{gff3_files} = [$tier->{VARS}{gff3_file}];
   }
   elsif($tier_type == 2){
       $tier->{RESULTS}->{section_files} = $tier->{VARS}{section_files};
   }
   elsif($tier_type == 3){
      $tier->{RESULTS}->{holdover_files} = $tier->{VARS}{holdover_files};
      $tier->{RESULTS}->{section_files} = $tier->{VARS}{section_files};
      $tier->{RESULTS}->{gff3_files} = [$tier->{VARS}{gff3_file}];
   }
   elsif($tier_type == 4){
      $tier->{RESULTS}->{p_fastas} = $tier->{VARS}->{p_fastas};
      $tier->{RESULTS}->{t_fastas} = $tier->{VARS}->{t_fastas};
      $tier->{RESULTS}->{n_fastas} = $tier->{VARS}->{n_fastas};
   }

   delete($tier->{VARS});

   #restore
   $tier->{VARS}{LOG} = $LOG;
   $tier->{VARS}{DS_CTL} = $DS_CTL;

   return;
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
   my $level = shift;
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

   if (defined($arg)) {
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

   $E->{-text} .= "ERROR: Failed while ".$extra."\n" if($extra);
   
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
sub retrieve {
    my @args = @_;
    
    try {
	sleep 5 if(! -e $args[0]);
	confess "ERROR: No such file or directory at $args[0]\n" if(! -e $args[0]);
	return Storable::retrieve(@args);
    }
    catch Error::Simple with {
        my $E = shift;
	
	print STDERR "Removing file: $args[0]\n";
	unlink($args[0]);
	throw $E;
    };
}
#-----------------------------------------------------------------------------
sub store {
    my @args = @_;

    try {
	my $file = $args[1];
	$args[1] =~ s/([^\/]+)$/.storetmp.$1.storetmp/;
	Storable::store(@args);
	move($args[1], $file)
	    or die "ERROR: Storable::store failed on $file\n";
    }
    catch Error::Simple with {
        my $E = shift;
	
	unlink($args[1]);
	throw $E;
    };
}
#-----------------------------------------------------------------------------
1;
