#! /usr/bin/perl -w

package Process::MakerTiers;

use FindBin;
use lib "$FindBin::Bin/../..";

use strict;
use Process::MakerChunk;
use File::Path;
use URI::Escape;

#-----------------------------------------------------------------------------
#------------------------------GLOBAL VARIABLES-------------------------------
#-----------------------------------------------------------------------------
my %SEEN;

#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------

#returns a new MakerTier object
sub new {
   my ($class, @args) = @_;

   my $self = {};

   bless ($self, $class);

   if (@args) {
      my $arg = shift @args;
      if (ref $arg eq 'Process::MakerChunk') {
	 $self = $arg->clone();
      }
      else {
	 $self->{VARS}{fasta} = $arg;
	 $self->{VARS}{CTL_OPTIONS} = \%{shift @args};
	 $self->{VARS}{OPT} = shift @args;
	 $self->{TIER_ID} = shift @args;
	 $self->{TERMINATE} = 0;
	 $self->{LEVEL}{CURRENT} = -1;
	 $self->_initialize();
      }
   }

   return $self;
}

#--------------------------------------------------------------
#initializes variables in the MakerTier object
#called from $self->new

sub _initialize {
   my $self        = shift;
   my $fasta       = $self->{VARS}{fasta};
   my %CTL_OPTIONS = %{$self->{VARS}{CTL_OPTIONS}};
   my %OPT         = %{$self->{VARS}{OPT}};

   $self->{VARS}{query_def} = Fasta::getDef($fasta); #Get fasta header
   ($self->{VARS}{seq_id})  = $self->{VARS}{query_def} =~ /^>(\S+)/; #Get identifier
   $self->{VARS}{query_seq} = Fasta::getSeq($fasta); #Get fasta Sequence
    
   #--error checking for non-unique ids in multi-fasta
   if (exists $SEEN{$self->{VARS}{seq_id}}) {
      warn "ERROR:  The multi-fasta file contains non-unique sequence ids.\n",
      "The id " . $self->{VARS}{seq_id} . " occurs more than once.\n";

      $SEEN{$self->{VARS}{seq_id}}++;
      my $new_id = $self->{VARS}{seq_id}.".V".$SEEN{$self->{VARS}{seq_id}};

      warn "\n\n\nThe id ".$self->{VARS}{seq_id}." will now be changed to $new_id \n";

      $self->{VARS}{seq_id} = $new_id;
      $SEEN{$new_id}++;
   }
   else {
      $SEEN{$self->{VARS}{seq_id}}++;
   }

   #--build a safe name for file names from the sequence identifier  
   $self->{VARS}{seq_out_name} = uri_escape($self->{VARS}{seq_id},
					    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
					   );

   #set up base directory for output
   $self->{VARS}{out_dir} = $CTL_OPTIONS{'out_base'};

   #use datastore as base if datastore flag is set
   if (exists $CTL_OPTIONS{'datastore'}) {
      $self->{VARS}{out_dir} = $CTL_OPTIONS{'datastore'}->id_to_dir($self->{VARS}{seq_out_name});
      $CTL_OPTIONS{'datastore'}->mkdir($self->{VARS}{seq_out_name}) ||
          die "ERROR: could not make directory $self->{VARS}{out_dir}\n";
   }

   #skip this fasta if gff3 output already exists
   if (-e $self->{VARS}{out_dir}."/".$self->{VARS}{seq_out_name}.".gff" && ! $OPT{f}){
      print STDERR "#----------------------------------------------------------------------\n",
                   "The contig:".$self->{VARS}{seq_id}." has already been processed!!\n",
                   "Maker will now skip to the next contig.\n",
                   "Run maker with the -f flag to force Maker to recompute all contig data.\n",
		   "#----------------------------------------------------------------------\n\n\n";
      $self->{TERMINATE} = 1;
      return;
   }

   #--set up void directory where analysis is stored
   $self->{VARS}{the_void}  = build_the_void($self->{VARS}{seq_out_name},
					     $self->{VARS}{out_dir}
					    );

   #-set up variables that are heldover from last chunk
   $self->{VARS}{holdover_chunk} = undef;

   #-set up variables that are the result of chunk accumulation
   $self->{VARS}{masked_total_seq} = '';
   $self->{VARS}{p_fastas} = '';
   $self->{VARS}{t_fastas} = '';
   $self->{VARS}{p_snap_fastas} = '';
   $self->{VARS}{t_snap_fastas} = '';

   $self->{VARS}{GFF3} = new Dumper::GFF::GFFV3();    
   $self->{VARS}{GFF3}->seq($self->{VARS}{query_seq});
   $self->{VARS}{GFF3}->seq_id($self->{VARS}{seq_id});

   return;
}
#--------------------------------------------------------------
#runs the MakerTier until multiple chunks are available to distribute

sub run {
   my $self = shift;
   my $current_level = $self->{LEVEL}{CURRENT};

   #---debug
   print STDERR "\n\n\nNow in LEVEL: $current_level\n\n\n";
   #---debug

   return undef if ($self->terminated);
   return undef if ($self->_level_started && ! $self->_level_finished);
   return $self->run if($self->_next_level);

   my %OPT = %{$self->{VARS}{OPT}};
   my %CTL_OPTIONS = %{$self->{VARS}{CTL_OPTIONS}};

   if ($current_level == -1){
      #==DECIDE REPEAT MASKING HERE
      if($OPT{R}){
	 print STDERR "Repeatmasking skipped!!\n";
	 $self->{VARS}{masked_total_seq} = ${$self->{VARS}{query_seq}};

	 $self->{LEVEL}{CURRENT} = 3;
	 $self->_initiate_level(3);
	 return $self->run;
      }
      elsif($OPT{GFF}){
	 $self->{VARS}{masked_total_seq} = repeat_mask_seq::gff(uc(${$self->{VARS}{query_seq}}), 
								$self->{VARS}{seq_id},
								$CTL_OPTIONS{'rm_gff'}
							       );

	 $self->{LEVEL}{CURRENT} = 3;
	 $self->_initiate_level(3);
	 return $self->run;
      }
      else{
	 #---
	 my $fasta_chunker = new FastaChunker();
	 $fasta_chunker->parent_fasta($self->{VARS}{fasta});
	 $fasta_chunker->chunk_size($CTL_OPTIONS{'max_dna_len'});
	 $fasta_chunker->min_size($CTL_OPTIONS{'split_hit'});
	 $fasta_chunker->load_chunks();

	 my $chunk_count = 0;

	 $self->{VARS}{fasta_chunker} = $fasta_chunker;
	 $self->{VARS}{chunk_count} = $chunk_count;
	 $self->{VARS}{f_chunk} = $fasta_chunker->get_chunk($chunk_count);
	 #---

	 $self->{LEVEL}{CURRENT} = 0;
	 $self->_initiate_level(0);

	 return $self->run;
      }
   }
   elsif($current_level == 0){#repeat masking individual chunk
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 1){#repeatmask blastx multiple makerChunks
      $self->_load_chunks_for_level($current_level);

      return; #pause here;
   }
   elsif($current_level == 2){#collecting repeatmask blastx results
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 3){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      my $chunk_count = ++$self->{VARS}{chunk_count};
      
      if ($self->{VARS}{f_chunk} = $self->{VARS}{fasta_chunker}->get_chunk($chunk_count)){
	 $self->{LEVEL}{CURRENT} = 0;
	 $self->_initiate_level(0);
      }

      return $self->run;
   }
   elsif($current_level == 4){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 5){
      #---get curent fasta chunk
      my $chunk_count = $self->{VARS}{chunk_count};
      $self->{VARS}{f_chunk} = $self->{VARS}{fasta_chunker}->get_chunk($chunk_count);
      #---

      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 6){
      $self->_load_chunks_for_level($current_level);

      return undef; #pause here
   }
   elsif($current_level == 7){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 8){
      $self->_load_chunks_for_level($current_level);

      return undef; #pause here
   }
   elsif($current_level == 9){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 10){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 11){
      $self->_load_chunks_for_level($current_level);
      my $chunk = $self->next_chunk;
      $chunk->run($self->id);
      $self->update_chunk($chunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 12){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 13){
      $self->_load_chunks_for_level($current_level);
      my $chunk = $self->next_chunk;
      $chunk->run($self->id);
      $self->update_chunk($chunk);
      $self->_polish_results();

      return $self->run;
   }
   elsif($current_level == 14){
      $self->_load_chunks_for_level($current_level);
      my $chunk = $self->next_chunk;
      $chunk->run($self->id);
      $self->update_chunk($chunk);
      $self->_polish_results();

      my $chunk_count = ++$self->{VARS}{chunk_count};
      
      if ($self->{VARS}{f_chunk} = $self->{VARS}{fasta_chunker}->get_chunk($chunk_count)){
	 $self->{LEVEL}{CURRENT} = 4;
	 $self->_initiate_level(4);
      }

      return $self->run;
   }
   elsif($current_level == 15){
      $self->_load_chunks_for_level($current_level);
      my $mChunk = $self->next_chunk;
      $mChunk->run($self->id);
      $self->update_chunk($mChunk);
      $self->_polish_results();

      #output error log for this contig
      my $err_file = $self->{VARS}{out_dir}."/".$self->{VARS}{seq_out_name}.".error";
      open(OUT, "> $err_file");
      print OUT $self->error;
      close(OUT);

      return $self->run;
   }
   else{
      $self->{TERMINATE} = 1;

      return undef; #stop here
   }
}

#--------------------------------------------------------------
#itteratively gets the next chunk from the tier
#returns undef if there is no chunk or the tier can not yet advance
#calls $self->run to try and advance the MakerTier

sub next_chunk {
   my $self = shift;
   my $current_level = $self->{LEVEL}{CURRENT};

   if ($current_level == -1 || ! $self->_level_started){
      $self->run;
      $current_level = $self->{LEVEL}{CURRENT};
   }
   elsif($self->_level_finished && ! $self->terminated){
      $self->run;
      $current_level = $self->{LEVEL}{CURRENT};
   }

   return undef if ($self->terminated);


   if (my $chunk = shift @{$self->{LEVEL}{$current_level}{CHUNKS}}) {
      return $chunk;
   }
   else{
      return undef;
   }
}

#--------------------------------------------------------------
#moves tier up one level is current level is finished

sub _next_level {
   my $self  = shift;
   my $level = $self->{LEVEL}{CURRENT};
    
   if ($level  == -1) {
      return undef;
   }
   elsif($self->terminated){
      return undef;
   }
   elsif (! $self->_level_started){
      return undef;
   }
   elsif (! $self->_level_finished) {
      return undef;
   }

   #--get results for current level
   $self->_polish_results();

   #--now go up one level
   $level++;
   $self->{LEVEL}{CURRENT} = $level;
   $self->_initiate_level($level);

   return 1;
}
#--------------------------------------------------------------
#defines variables for new level and undefines last level

sub _initiate_level {
   my $self = shift;
   my $level = shift;

   #undef previous level and current level
   $self->{LEVEL}{$level - 1} = undef;
   $self->{LEVEL}{$level} = undef;

   #--initiate level variables
   $self->{LEVEL}{$level}{CHUNK_COUNT} = 0;
   $self->{LEVEL}{$level}{RESULT_COUNT} = 0;
   $self->{LEVEL}{$level}{CHUNKS} = undef;
   $self->{LEVEL}{$level}{RESULTS} = undef;
   $self->{LEVEL}{$level}{ERROR} = undef;
   $self->{LEVEL}{$level}{STARTED} = undef;
}
#--------------------------------------------------------------
#gets variables needed for chunks from the current level
#also constructs chunks using _build_chunk 

sub _load_chunks_for_level {
   my $self = shift;
   my $level = shift;

   $self->{LEVEL}{$level}{STARTED} = 1;

   #--select variables to send to Process::MakerChunk object
   if ($level == 0) {
      #----------------------CLEAR_MEMORY
      $self->{VARS}{rma_keepers} = [];
      $self->{VARS}{repeat_blastx_keepers} = [];
      #----------------------CLEAR_MEMORY

      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{the_void},
		  $self->{VARS}{seq_out_name},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 1) {
      foreach my $repeat_protein (@{$self->{VARS}{CTL_OPTIONS}{repeat_protein}}) {
	 #------------------------ARGS_IN
	 my @args =( $self->{VARS}{f_chunk},
		     $repeat_protein,
		     $self->{VARS}{the_void},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	 #------------------------ARGS_IN
		
	 #-------------------------CHUNK
	 $self->_build_chunk($level,\@args);
	 #-------------------------CHUNK
      }
	
      return 1;
   }
   elsif ($level == 2) {
      	 #------------------------ARGS_IN
	 my @args =( $self->{VARS}{f_chunk},
		     $self->{VARS}{rep_blastx_res_dir},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	 #------------------------ARGS_IN
		
	 #-------------------------CHUNK
	 $self->_build_chunk($level,\@args);
	 #-------------------------CHUNK
   }
   elsif ($level == 3) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{rma_keepers},
		  $self->{VARS}{repeat_blastx_keepers},
		  $self->{VARS}{GFF3},
		  $self->{VARS}{query_def},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{masked_total_seq},
		  $self->{VARS}{the_void},
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 4) {
      #----------------------CLEAR_MEMORY
      $self->{VARS}{rma_keepers} = [];
      $self->{VARS}{repeat_blastx_keepers} = [];
      #----------------------CLEAR_MEMORY

      #------------------------ARGS_IN
      my @args =( $self->{VARS}{masked_total_seq},
		  $self->{VARS}{the_void},
		  $self->{VARS}{seq_out_name},
		  $self->{VARS}{query_def},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 5) {
      #----------------------CLEAR_MEMORY
      $self->{VARS}{blastn_keepers} = [];
      $self->{VARS}{blastx_keepers} = [];
      $self->{VARS}{snaps_on_chunk} = [];
      $self->{VARS}{augus_on_chunk} = [];
      $self->{VARS}{blastn_data} = [];
      $self->{VARS}{blastx_data} = [];
      $self->{VARS}{exonerate_e_data} = []; 
      $self->{VARS}{exonerate_p_data} = []; 
      $self->{VARS}{annotations} = []; 
      #----------------------CLEAR_MEMORY

      #------------------------ARGS_IN
      my @args =( $self->{VARS}{holdover_chunk},
		  $self->{VARS}{f_chunk}
		);
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 6) {
      foreach my $transcripts (@{$self->{VARS}{CTL_OPTIONS}{est}}) {
	 #------------------------ARGS_IN
	 my @args =( $self->{VARS}{f_chunk},
		     $transcripts,
		     $self->{VARS}{the_void},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	 #------------------------ARGS_IN

	 #-------------------------CHUNK
	 $self->_build_chunk($level,\@args);
	 #-------------------------CHUNK
      }

      return 1;
   }
   elsif ($level == 7) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{blastn_res_dir},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 8) {
      foreach my $proteins (@{$self->{VARS}{CTL_OPTIONS}{protein}}) {
	 #------------------------ARGS_IN
	 my @args =( $self->{VARS}{f_chunk},
		     $proteins,
		     $self->{VARS}{the_void},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	 #------------------------ARGS_IN

	 #-------------------------CHUNK
	 $self->_build_chunk($level,\@args);
	 #-------------------------CHUNK
      }

      return 1;
   }
   elsif ($level == 9) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{blastx_res_dir},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 10) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{holdover_chunk},
		  $self->{VARS}{snaps},
		  $self->{VARS}{augus},
		  $self->{VARS}{blastx_keepers},
		  $self->{VARS}{blastn_keepers},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{CTL_OPTIONS}{split_hit},
		  $self->{VARS}{the_void}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 11) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{fasta},
		  $self->{VARS}{blastx_keepers},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 12) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{fasta},
		  $self->{VARS}{blastn_keepers},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 13) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{fasta},
		  $self->{VARS}{masked_fasta},
		  $self->{VARS}{f_chunk}->number,
		  $self->{VARS}{exonerate_p_data},
		  $self->{VARS}{exonerate_e_data},
		  $self->{VARS}{blastx_data},
		  $self->{VARS}{snaps_on_chunk},
		  $self->{VARS}{augus_on_chunk},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f},
		  $self->{VARS}{OPT}{SNAPS}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 14) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{blastx_data},
		  $self->{VARS}{blastn_data},
		  $self->{VARS}{exonerate_p_data},
		  $self->{VARS}{exonerate_e_data},
		  $self->{VARS}{annotations},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{snaps_on_chunk},
		  $self->{VARS}{augus_on_chunk},
		  $self->{VARS}{p_fastas},
		  $self->{VARS}{t_fastas},
		  $self->{VARS}{p_snap_fastas},
		  $self->{VARS}{t_snap_fastas},
		  $self->{VARS}{p_augus_fastas},
		  $self->{VARS}{t_augus_fastas},
		  $self->{VARS}{GFF3}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 15) {
      #----------------------CLEAR_MEMORY
      $self->{VARS}{blastn_keepers} = [];
      $self->{VARS}{blastx_keepers} = [];
      $self->{VARS}{snaps_on_chunk} = [];
      $self->{VARS}{augus_on_chunk} = [];
      $self->{VARS}{blastn_data} = [];
      $self->{VARS}{blastx_data} = [];
      $self->{VARS}{exonerate_e_data} = []; 
      $self->{VARS}{exonerate_p_data} = []; 
      $self->{VARS}{annotations} = []; 
      #----------------------CLEAR_MEMORY

      #------------------------ARGS_IN
      my @args =( $self->{VARS}{p_fastas},
		  $self->{VARS}{t_fastas},
		  $self->{VARS}{p_snap_fastas},
		  $self->{VARS}{t_snap_fastas},
		  $self->{VARS}{p_augus_fastas},
		  $self->{VARS}{t_augus_fastas},
		  $self->{VARS}{GFF3},
		  $self->{VARS}{seq_out_name},
		  $self->{VARS}{out_dir},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   else {
      return undef;
   }

   return 1;
}
#--------------------------------------------------------------
#processes results from each level once all chunks are in

sub _polish_results{
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   if ($level < 0 || ! $self->_level_finished()){
      return undef;
   }

   if (not $self->{LEVEL}{$level}{RESULTS}) {
      return undef;
   }

   #---clear value in specific variable before redefining
   if ($level == 1) {
       $self->{VARS}{repeat_blastx_keepers} = [];
   }
   elsif ($level == 7) {
       $self->{VARS}{blastn_keepers} = [];
   }
   elsif ($level == 9) {
       $self->{VARS}{blastx_keepers} = [];
   }

   #collect values from result
   foreach my $result (@{$self->{LEVEL}{$level}{RESULTS}}) {
      my @results = @{$result};

      if ($level == 0) {
	 #------------------------RESULTS
	 $self->{VARS}{f_chunk}     = shift @results;
	 $self->{VARS}{rma_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 1) {
	 #------------------------RESULTS
	 $self->{VARS}{rep_blastx_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 2) {
	 #------------------------RESULTS
	 $self->{VARS}{repeat_blastx_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 3) {
	 #------------------------RESULTS
	 $self->{VARS}{f_chunk}          = shift @results;
	 $self->{VARS}{masked_total_seq} = shift @results;
	 $self->{VARS}{GFF3}             = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 4) {
	 #------------------------RESULTS
	 $self->{VARS}{masked_fasta}  = shift @results;
	 $self->{VARS}{snaps}         = shift @results;
	 $self->{VARS}{augus}         = shift @results;
	 $self->{VARS}{fasta_chunker} = shift @results;
	 $self->{VARS}{chunk_count}   = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 5) {
	 #------------------------RESULTS
	 $self->{VARS}{f_chunk} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 6) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 7) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 8) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 9) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 10) {
	 #------------------------RESULTS
	 $self->{VARS}{holdover_chunk} = shift @results;
	 $self->{VARS}{blastx_keepers} = shift @results;
	 $self->{VARS}{blastn_keepers} = shift @results;
	 $self->{VARS}{snaps_on_chunk} = shift @results;
	 $self->{VARS}{augus_on_chunk} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 11) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_data}      = shift @results;
	 $self->{VARS}{exonerate_p_data} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 12) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_data}      = shift @results;
	 $self->{VARS}{exonerate_e_data} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 13) {
	 #------------------------RESULTS
	 $self->{VARS}{annotations} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 14) {
	 #------------------------RESULTS
	 $self->{VARS}{GFF3}           = shift @results;
	 $self->{VARS}{p_fastas}       = shift @results;
	 $self->{VARS}{t_fastas}       = shift @results;
	 $self->{VARS}{p_snap_fastas}  = shift @results;
	 $self->{VARS}{t_snap_fastas}  = shift @results;
	 $self->{VARS}{p_augus_fastas} = shift @results;
	 $self->{VARS}{t_augus_fastas} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 15) {
	 #------------------------RESULTS
	 #------------------------RESULTS
      }
   }

   return 1;
}
#--------------------------------------------------------------
#builds each individual chunk for processing
#called from $self->_load_chunks_for_level

sub _build_chunk {
   my $self = shift;
   my $level = shift;
   my $args = shift;

   my $chunk_id = $self->id().":".$level.":". $self->{LEVEL}{$level}{CHUNK_COUNT};
    
   my $chunk = Process::MakerChunk->new($level, $args, $chunk_id);
   push (@{$self->{LEVEL}{$level}{CHUNKS}}, $chunk);
   $self->{LEVEL}{$level}{CHUNK_COUNT}++;
}
#--------------------------------------------------------------
#returns a line giving name of sequence and the location where
#analysis is stored.  Used for datastore implementation.

sub DS {
   my $self = shift;
   return $self->{VARS}{seq_id} . "\t" . $self->{VARS}{out_dir};
}

#--------------------------------------------------------------
#returns the number of chunks currently available for processing

sub num_chunks {
   my $self = shift;

   $self->_next_level();
   my $current = $self->{LEVEL}{CURRENT};

   my $num = @{$self->{LEVEL}{$current}{CHUNKS}};

   return $num;
}
#--------------------------------------------------------------
#returns a copy of this MakerTiers object

sub clone {
   my $self = shift;
    
   my $clone = dclone($self);

   return $clone;
}
#--------------------------------------------------------------
#returns true if level is finished

sub _level_finished {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   if ($self->{LEVEL}{$level}{CHUNK_COUNT} == $self->{LEVEL}{$level}{RESULT_COUNT}) {
      return 1;
   }

   return 0;
}

#--------------------------------------------------------------
#returns true if level has been been started, i.e. had chunks in it

sub _level_started {
   my $self = shift;
   my $level = $self->{LEVEL}{CURRENT};

   return undef if (! exists $self->{LEVEL}{$level});
   return $self->{LEVEL}{$level}{STARTED};
}
#--------------------------------------------------------------
#returns true if tier has terminated

sub terminated {
   my $self = shift;

   return $self->{TERMINATE};
}

#--------------------------------------------------------------
#takes a chunk and adds its results to the tier

sub update_chunk {
   my $self = shift;
   my $chunk = shift;

   my $id = $chunk->id();
   my $result = $chunk->result();

   my ($tier_id, $level_num, $chunk_num) = split (":", $id);
    
   push (@{$self->{LEVEL}{$level_num}{RESULTS}}, $result);
   $self->{LEVEL}{$level_num}{RESULT_COUNT}++;
    
   $self->{LEVEL}{$level_num}{ERROR} .= $chunk->error();
   $self->{ERROR} .= $chunk->error();
   print STDERR $chunk->error();
}
#-------------------------------------------------------------
sub error{
   my $self = shift;
   return $self->{ERROR};
}
#-------------------------------------------------------------
#returns the id of the tier

sub id {
   my $self = shift @_;
   return $self->{TIER_ID};
}

#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
sub build_the_void {
   my $seq_id  = shift;
   my $out_dir = shift;

   $out_dir =~ s/\/$//;

   my $vid = "theVoid\.$seq_id";
   my $the_void = "$out_dir/$vid";
   mkpath ($the_void);

   return $the_void;
}

#-----------------------------------------------------------------------------
1;
