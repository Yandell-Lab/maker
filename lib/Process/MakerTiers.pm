#! /usr/bin/perl -w

package Process::MakerTiers;

use FindBin;
use lib "$FindBin::Bin/../..";

use strict;
use Error qw(:try);
use Error::Simple;
use Process::MakerChunk;
use File::Path;
use URI::Escape;
use runlog;
use Shared_Functions;
use File::Temp qw(tempfile);
use File::Copy;

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
	 $self->{FAILED} = 0;
	 $self->{LEVEL}{CURRENT} = -1;
	 $self->_initialize();
	 $self->_continue();
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

   try{
      $self->{VARS}{query_def} = Fasta::getDef($fasta); #Get fasta header
      ($self->{VARS}{seq_id})  = $self->{VARS}{query_def} =~ /^>(\S+)/; #Get identifier
      $self->{VARS}{query_seq} = Fasta::getSeq($fasta); #Get fasta Sequence
      
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
	    die "ERROR: could not make directory ".$self->{VARS}{out_dir}."\n";
      }
      
      
      #--set up void directory where analysis is stored
      $self->{VARS}{the_void}  = Shared_Functions::build_the_void($self->{VARS}{seq_out_name},
								  $self->{VARS}{out_dir}
								 );

   }
   catch Error::Simple with{
      my $E = shift;
      
      print STDERR $E->{-text};
      print STDERR "\n\nMaker failed while examining contents of the fasta file!!\n\n";
      my $code = 2;
      $code = $E->{-value} if (defined($E->{-value}));
      
      exit($code);
   };

   try{
      #-build and proccess the run log
      $self->{LOG} = runlog->new(\%CTL_OPTIONS, \%OPT, $self->{VARS}{the_void}, "run.log");
   }
   catch Error::Simple with{
      my $E = shift;
      
      print STDERR $E->{-text};
      print STDERR "\n\nMaker failed while trying to building/processing the run.log file!!\n\n";
      my $code = 2;
      $code = $E->{-value} if (defined($E->{-value}));
      
      exit($code);
   };
   
   return;
}
#--------------------------------------------------------------
#decide whether to continue with contig based on length and LOG continue flag
sub _continue {
   my $self = shift;

   if (defined($self->{CONTINUE})){
      $self->{TERMINATE} = 1 if ($self->{CONTINUE} <= 0);
      return $self->{CONTINUE};
   }

   my $LOG = $self->{LOG};

   #get variables
   my %OPT = %{$self->{VARS}{OPT}};
   my %CTL_OPTIONS = %{$self->{VARS}{CTL_OPTIONS}};
   my $seq_id = $self->{VARS}{seq_id};
   my $out_dir = $self->{VARS}{out_dir};
   my $seq_out_name = $self->{VARS}{seq_out_name};

   #skip contig if too short
   if (length(${$self->{VARS}{query_seq}}) < $CTL_OPTIONS{min_contig}){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Skipping the contig \'$seq_id\' because it is too short\n",
                   "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tSKIPPED_SMALL\n";
      }
      
      $self->{CONTINUE} = -2; #skipped signal is -2
      $self->{TERMINATE} = 1;
      return $self->{CONTINUE};
   }

   #==Decide whether to skip the current contig based on log
   my $continue_flag = $LOG->get_continue_flag();
   if ($continue_flag == 1) {
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now starting the contig:$seq_id!!\n",
                   "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tSTARTED";
      }
   }
   elsif ($continue_flag == 0) {
      print STDERR "#---------------------------------------------------------------------\n",
      "The contig:$seq_id has already been processed!!\n",
      "Maker will now skip to the next contig.\n",
      "Run maker with the -f flag to force Maker to recompute all contig data.\n",
      "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tFINISHED";
      }

      $self->{TERMINATE} = 1;
   }
   elsif ($continue_flag == -1) {
      print STDERR "#---------------------------------------------------------------------\n",
      "The contig:$seq_id failed on the last run!!\n",
      "Maker will now skip to the next contig rather than try again.\n",
      "Run maker with the -f flag to force Maker to recompute all contig data.\n",
      "Run maker with the -died flag to have Maker retry data that failed.\n",
      "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tDIED_SKIPPED";
      }

      $self->{TERMINATE} = 1;
   }
   elsif ($continue_flag == 2) {
      print STDERR "#---------------------------------------------------------------------\n",
      "The failed contig:$seq_id will now run again!!\n",
      "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tRETRY";
      }
   }
   elsif ($continue_flag == -2) {
      print STDERR "#---------------------------------------------------------------------\n",
      "Skipping the contig:$seq_id!!\n",
      "However this contig is still not finished!!\n",
      "#---------------------------------------------------------------------\n\n\n";

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tSKIPPED";
      }

      $self->{TERMINATE} = 1;
   }
   elsif ($continue_flag == -3) {
      my $die_count = $LOG->get_die_count();
      print STDERR "#---------------------------------------------------------------------\n",
      "The contig:$seq_id failed $die_count time!!\n",
      "Maker will not try again!!\n",
      "The contig will be stored in $out_dir/$seq_out_name.died.fasta\n",
      "You can use this fasta file to debug and re-run this sequence\n",
      "#---------------------------------------------------------------------\n\n\n";

      open (DFAS, "> $out_dir/$seq_out_name.died.fasta");
      print DFAS $self->{VARS}{fasta};
      close (DFAS);

      if ($OPT{d}) {
	 $self->{DS} = "$seq_id\t$out_dir\tDIED_SKIPPED_PERMANENT";
      }

      $self->{TERMINATE} = 1;
   }

   $self->{CONTINUE} = $continue_flag;
   
   return $self->{CONTINUE};
}
#--------------------------------------------------------------
#exception handler
sub _handler {
   my $self = shift;
   my $E = shift;
   my $extra = shift;

   my $seq_id = $self->{VARS}{seq_id};
   my $LOG = $self->{LOG};

   print STDERR $E->{-text};
   print STDERR $extra;
   print STDERR "FAILED CONTIG:$seq_id\n\n";
   
   my $die_count = $LOG->get_die_count();
   $die_count++;
   
   $LOG->add_entry("DIED","RANK", $self->id);
   $LOG->add_entry("DIED","COUNT",$die_count);

   $self->{DS} = "$self->{VARS}{seq_id}\t$self->{VARS}{out_dir}\tDIED";
   $self->{FAILED} = 1;
}
#--------------------------------------------------------------
sub run {
   my $self = shift;

   my $ret = $self->_run();

   return $ret;
}
#--------------------------------------------------------------
#runs the MakerTier until multiple chunks are available to distribute

sub _run {
   my $self = shift;
   my $current_level = $self->{LEVEL}{CURRENT};

   #---debug
   #print STDERR "\n\n\nNow in LEVEL: $current_level\n\n\n";
   #---debug

   return undef if ($self->terminated || $self->failed);
   return undef if ($self->_level_started && ! $self->_level_finished);
   return $self->run if($self->_next_level);

   my %OPT = %{$self->{VARS}{OPT}};
   my %CTL_OPTIONS = %{$self->{VARS}{CTL_OPTIONS}};

   if ($current_level == -1){#initiation level no chunk associated
      #-set up variables that are heldover from last chunk
      $self->{VARS}{holdover_blastn} = [];
      $self->{VARS}{holdover_blastx} = [];
      $self->{VARS}{holdover_tblastx} = [];
      $self->{VARS}{holdover_preds} = [];

      #-set up variables that are the result of chunk accumulation
      $self->{VARS}{masked_total_seq} = '';
      $self->{VARS}{p_fastas} = '';
      $self->{VARS}{t_fastas} = '';
      
      $self->{VARS}{GFF3} = new Dumper::GFF::GFFV3();    
      $self->{VARS}{GFF3}->seq($self->{VARS}{query_seq});
      $self->{VARS}{GFF3}->seq_id($self->{VARS}{seq_id});

      #==DECIDE REPEAT MASKING HERE
      if($OPT{R}){
	 try{
	    print STDERR "Repeatmasking skipped!!\n";
	    $self->{VARS}{masked_total_seq} = ${$self->{VARS}{query_seq}};
	    
	    $self->{LEVEL}{CURRENT} = 4;
	    $self->_initiate_level(4);
	 }
	 catch Error::Simple with{
	    my $E = shift;
	    $self->_handler($E, "\n\nMaker failed after skip repeat masking prep!!\n");
	 };

	 return $self->run;
      }
      elsif($CTL_OPTIONS{rm_gff}){
	 try{
	    $self->{VARS}{masked_total_seq} = repeat_mask_seq::gff(uc(${$self->{VARS}{query_seq}}), 
								   $self->{VARS}{seq_id},
								   $CTL_OPTIONS{'rm_gff'}
								  );
	    
	    $self->{LEVEL}{CURRENT} = 4;
	    $self->_initiate_level(4);
	 }
	 catch Error::Simple with{
	    my $E = shift;
	    $self->_handler($E, "\n\nMaker failed at repeat GFF proccessing!!\n");
	 };

	 return $self->run;
      }
      else{
	 try{
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
	 }
	 catch Error::Simple with{
	    my $E = shift;
	    $self->_handler($E, "\n\nMaker failed at rpeat masking preperation!!\n");
	 };

	 return $self->run;
      }
   }
   elsif($current_level == 0){#repeat masking individual chunk
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at RepeatMasker!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 1){#repeatmask blastx multiple makerChunks
      try{
	 $self->_load_chunks_for_level($current_level);
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at transposable element masking!!\n");
      };

      return; #pause here;
   }
   elsif($current_level == 2){#collecting repeatmask blastx results
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed while collecting transposable element alignments!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 3){
      try{
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
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at masking chunk with repeats!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 4){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at ab-initio gene predictions!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 5){
      try{
	 #---get curent fasta chunk
	 my $chunk_count = $self->{VARS}{chunk_count};
	 $self->{VARS}{f_chunk} = $self->{VARS}{fasta_chunker}->get_chunk($chunk_count);
	 #---
	 
	 $self->_load_chunks_for_level($current_level);
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at performing blastn of ESTs!!\n");
      };

      return undef; #pause here
   }
   elsif($current_level == 6){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at collecting blastn of ESTs!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 7){
      try{
	 $self->_load_chunks_for_level($current_level);
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at performing blastx of proteins!!\n");
      };

      return undef; #pause here
   }
   elsif($current_level == 8){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at collecting blastx of proteins!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 9){
       try{
	   $self->_load_chunks_for_level($current_level);
       }
       catch Error::Simple with{
	   my $E = shift;
	   $self->_handler($E, "\n\nMaker failed at performing tblastx of alt-ESTs!!\n");
       };

       return undef; #pause here                                                                                                                                                                                                                                                     
   }
   elsif($current_level == 10){
       try{
	   $self->_load_chunks_for_level($current_level);
	   my $mChunk = $self->next_chunk;
	   $mChunk->run($self->id);
	   $self->update_chunk($mChunk);
	   $self->_polish_results();
       }
       catch Error::Simple with{
	   my $E = shift;
	   $self->_handler($E, "\n\nMaker failed at collecting tblastx of alt-ESTs!!\n");
       };

       return $self->run;
   }
   elsif($current_level == 11){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed while proccessing the chunk divide!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 12){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $chunk = $self->next_chunk;
	 $chunk->run($self->id);
	 $self->update_chunk($chunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at exonerate against proteins!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 13){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at exonerate against transcripts!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 14){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $chunk = $self->next_chunk;
	 $chunk->run($self->id);
	 $self->update_chunk($chunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at annotation generation!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 15){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $chunk = $self->next_chunk;
	 $chunk->run($self->id);
	 $self->update_chunk($chunk);
	 $self->_polish_results();
	 
	 my $chunk_count = ++$self->{VARS}{chunk_count};
	 
	 if ($self->{VARS}{f_chunk} = $self->{VARS}{fasta_chunker}->get_chunk($chunk_count)){
	    $self->{LEVEL}{CURRENT} = 5;
	    $self->_initiate_level(5);
	 }
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed at processing annotations on chunk!!\n");
      };

      return $self->run;
   }
   elsif($current_level == 16){
      try{
	 $self->_load_chunks_for_level($current_level);
	 my $mChunk = $self->next_chunk;
	 $mChunk->run($self->id);
	 $self->update_chunk($mChunk);
	 $self->_polish_results();
      }
      catch Error::Simple with{
	 my $E = shift;
	 $self->_handler($E, "\n\nMaker failed while writing final data!!\n");
      };

      return $self->run;
   }
   else{
      $self->{TERMINATE} = 1;
      $self->{DS} = "$self->{VARS}{seq_id}\t$self->{VARS}{out_dir}\tFINISHED";

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

      if($self->failed){ #clear the queue
	 $self->update_chunk($chunk);
	 while (my $chunk = shift @{$self->{LEVEL}{$current_level}{CHUNKS}}) {
	    $self->update_chunk($chunk);
	 }
	 return undef;
      }

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
   elsif($self->terminated || $self->failed){
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
   $self->{LEVEL}{$level}{CHUNKS} = [];
   $self->{LEVEL}{$level}{RESULTS} = [];
   $self->{LEVEL}{$level}{ERROR} = undef;
   $self->{LEVEL}{$level}{STARTED} = undef;
}
#--------------------------------------------------------------
#gets variables needed for chunks from the current level
#also constructs chunks using _build_chunk 

sub _load_chunks_for_level {
   my $self = shift;
   my $level = shift;

   return undef if ($self->terminated || $self->failed);

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
      $self->{VARS}{snaps_on_chunk} = [];
      $self->{VARS}{augus_on_chunk} = [];
      $self->{VARS}{blastn_keepers} = [];
      $self->{VARS}{blastx_keepers} = [];
      $self->{VARS}{tblastx_keepers} = [];
      $self->{VARS}{blastn_data} = [];
      $self->{VARS}{blastx_data} = [];
      $self->{VARS}{tblastx_data} = [];
      $self->{VARS}{exonerate_e_data} = []; 
      $self->{VARS}{exonerate_p_data} = []; 
      $self->{VARS}{annotations} = [];
      #----------------------CLEAR_MEMORY

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
   elsif ($level == 6) {
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
   elsif ($level == 7) {
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
   elsif ($level == 8) {
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
   elsif ($level == 9) {
       foreach my $alt_ests (@{$self->{VARS}{CTL_OPTIONS}{alt_est}}) {
         #------------------------ARGS_IN
	   my @args =( $self->{VARS}{f_chunk},
		       $alt_ests,
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
   elsif ($level == 10) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{tblastx_res_dir},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 11) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{f_chunk},
		  $self->{VARS}{masked_fasta},
		  $self->{VARS}{snaps},
		  $self->{VARS}{augus},
		  $self->{VARS}{blastn_keepers},
		  $self->{VARS}{blastx_keepers},
		  $self->{VARS}{tblastx_keepers},
		  $self->{VARS}{fasta_t_index},
		  $self->{VARS}{fasta_p_index},
		  $self->{VARS}{fasta_a_index},
		  $self->{VARS}{holdover_blastn},
		  $self->{VARS}{holdover_blastx},
		  $self->{VARS}{holdover_tblastx},
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
		  $self->{VARS}{blastx_keepers},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{fasta_p_index},
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
		  $self->{VARS}{blastn_keepers},
		  $self->{VARS}{tblastx_keepers},
		  $self->{VARS}{query_seq},
		  $self->{VARS}{fasta_t_index},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 14) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{fasta},
		  $self->{VARS}{masked_fasta},
		  $self->{VARS}{f_chunk}->number,
		  $self->{VARS}{exonerate_p_data},
		  $self->{VARS}{exonerate_e_data},
		  $self->{VARS}{blastx_data},
		  $self->{VARS}{preds_on_chunk},
		  $self->{VARS}{the_void},
		  $self->{VARS}{CTL_OPTIONS},
		  $self->{VARS}{OPT}{f},
		  $self->{VARS}{OPT}{PREDS}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 15) {
      #------------------------ARGS_IN
      my @args =( $self->{VARS}{blastx_data},
		  $self->{VARS}{blastn_data},
		  $self->{VARS}{tblastx_data},
		  $self->{VARS}{exonerate_p_data},
		  $self->{VARS}{exonerate_e_data},
		  $self->{VARS}{annotations},
		  $self->{VARS}{p_fastas},
		  $self->{VARS}{t_fastas},
		  $self->{VARS}{GFF3}
		);
      #------------------------ARGS_IN

      #-------------------------CHUNK
      $self->_build_chunk($level,\@args);
      #-------------------------CHUNK
   }
   elsif ($level == 16) {
      #----------------------CLEAR_MEMORY
      $self->{VARS}{blastn_keepers} = [];
      $self->{VARS}{blastx_keepers} = [];
      $self->{VARS}{tblastx_keepers} = [];
      $self->{VARS}{snaps_on_chunk} = [];
      $self->{VARS}{augus_on_chunk} = [];
      $self->{VARS}{blastn_data} = [];
      $self->{VARS}{blastx_data} = [];
      $self->{VARS}{tblastx_data} = [];
      $self->{VARS}{exonerate_e_data} = []; 
      $self->{VARS}{exonerate_p_data} = []; 
      $self->{VARS}{annotations} = []; 
      #----------------------CLEAR_MEMORY

      #------------------------ARGS_IN
      my @args =( $self->{VARS}{snaps},
		  $self->{VARS}{augus},
		  $self->{VARS}{p_fastas},
		  $self->{VARS}{t_fastas},
		  $self->{VARS}{GFF3},
		  $self->{VARS}{seq_out_name},
		  $self->{VARS}{out_dir},
		  $self->{VARS}{the_void},
		  $self->{VARS}{query_seq},
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

   return undef if($self->terminated || $self->failed);

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
   elsif ($level == 6) {
       $self->{VARS}{blastn_keepers} = [];
   }
   elsif ($level == 8) {
       $self->{VARS}{blastx_keepers} = [];
   }
   elsif ($level == 10) {
       $self->{VARS}{tblastx_keepers} = [];
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
	 $self->{VARS}{fasta_t_index}   = shift @results;
	 $self->{VARS}{fasta_p_index}   = shift @results;
	 $self->{VARS}{fasta_a_index}   = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 5) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 6) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 7) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 8) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 9) {
	 #------------------------RESULTS
	 $self->{VARS}{tblastx_res_dir} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 10) {
	 #------------------------RESULTS
	 $self->{VARS}{tblastx_keepers} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 11) {
	 #------------------------RESULTS
	 $self->{VARS}{holdover_blastn}  = shift @results;
	 $self->{VARS}{holdover_blastx}  = shift @results;
	 $self->{VARS}{holdover_tblastx} = shift @results;
	 $self->{VARS}{holdover_preds}   = shift @results;
	 $self->{VARS}{blastn_keepers}   = shift @results;
	 $self->{VARS}{blastx_keepers}   = shift @results;
	 $self->{VARS}{tblastx_keepers}  = shift @results;
	 $self->{VARS}{preds_on_chunk}   = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 12) {
	 #------------------------RESULTS
	 $self->{VARS}{blastx_data}      = shift @results;
	 $self->{VARS}{exonerate_p_data} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 13) {
	 #------------------------RESULTS
	 $self->{VARS}{blastn_data}      = shift @results;
	 $self->{VARS}{tblastx_data}      = shift @results;
	 $self->{VARS}{exonerate_e_data} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 14) {
	 #------------------------RESULTS
	 $self->{VARS}{annotations} = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 15) {
	 #------------------------RESULTS
	 $self->{VARS}{GFF3}           = shift @results;
	 $self->{VARS}{p_fastas}       = shift @results;
	 $self->{VARS}{t_fastas}       = shift @results;
	 #------------------------RESULTS
      }
      elsif ($level == 16) {
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
    
   my $chunk = Process::MakerChunk->new($level, $args, $chunk_id,  $self->{LOG});
   push (@{$self->{LEVEL}{$level}{CHUNKS}}, $chunk);
   $self->{LEVEL}{$level}{CHUNK_COUNT}++;
}
#--------------------------------------------------------------
#returns a line giving name of sequence and the location where
#analysis is stored.  Used for datastore implementation.

sub DS {
   my $self = shift;
   return $self->{DS} || $self->{VARS}{seq_id} . "\t" . $self->{VARS}{out_dir};
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

   $self->{TERMINATE} = 1 if ($self->failed && $self->_level_finished);

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
    
   #$self->{LEVEL}{$level_num}{ERROR} .= $chunk->error();
   $self->{ERROR} .= $chunk->error();

   if($chunk->failed){
      my $E = $chunk->exception;
      $self->_handler($E, "\n\nMaker failed at chunk update in level $level_num!!\n");
   }
   print STDERR $chunk->error();
}
#-------------------------------------------------------------
sub failed{
   my $self = shift;
   return $self->{FAILED};
}
#-------------------------------------------------------------
sub fasta{
   my $self = shift;
   return $self->{VARS}{fasta};
}
#-------------------------------------------------------------
sub error{
   my $self = shift;
   return $self->{ERROR} || '';
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
1;
