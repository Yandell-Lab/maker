#! /usr/bin/perl -w

package Process::MakerChunk;

use strict;
use vars qw($TMP);

use Error qw(:try);
use Error::Simple;
use Storable qw (freeze thaw dclone);

use FindBin;
use lib "$FindBin::Bin/../..";
use File::Temp qw(tempfile tempdir);
use Dumper::GFF::GFFV3;
use Dumper::XML::Game;
use Datastore::MD5;
use URI::Escape;
use File::Path;
use Data::Dumper;
use Getopt::Long;
use FileHandle;
use PostData;
use Cwd qw(cwd abs_path);
use Fasta;
use Iterator::Fasta;
use FastaChunker;
use Widget::RepeatMasker;
use Widget::blastx;
use Widget::tblastx;
use Widget::blastn;
use Widget::snap; 
use PhatHit_utils;
use Shadower;
use Bio::DB::Fasta;
use polisher::exonerate::protein;
use polisher::exonerate::est;
use maker::auto_annotator;
use cluster;
use repeat_mask_seq;
use maker::sens_spec;
use runlog;
use Shared_Functions;

$TMP = tempdir("maker_XXXXXX", CLEANUP => 1, TMPDIR => 1);

#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------
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
	 $self->_initialize($arg, @args);
      }
   }

   return $self;
}

#--------------------------------------------------------------
sub _initialize{
   my $self       = shift;
   $self->{LEVEL} = shift;
   $self->{VARS}  = shift; #this should be an array reference
   $self->{ID}    = shift || "0:".$self->{LEVEL}.":0";
   $self->{LOG}   = shift;
   $self->{RANK}  = shift || 0;

   $self->{LOG_FLAG} = ($self->{ID} =~ /\:0$/) ? 1 : 0;

}

#--------------------------------------------------------------
sub run {
   my $self = shift;
   $self->{RANK} = shift || $self->{RANK} || 0;

   my $ret;			#return value

   try{
      $ret = $self->_run();
   }
   catch Error::Simple with{
      my $E = shift;
      $self->{FAILED} = 1;
      $self->{EXCEPTION} = $E;
   };

   return $ret;
}

#--------------------------------------------------------------
sub _run {
   my $self = shift;

   if (exists $self->{RESULT}) {
      return;
   }

   my $level = $self->{LEVEL};
   my $vars = $self->{VARS};
   my @results;

   if ($level == 0) {
      #------------------------ARGS_IN
      my $chunk        = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      my $rma_keepers = [];

      #-- repeatmask the input file
      if(! $CTL_OPTIONS{rmlib_only}){
	 $chunk->seq(uc($chunk->seq())); #must be upper case before soft masking
	 
	 my $rma_keepers1 = Shared_Functions::repeatmask($chunk,
							 $the_void,
							 $seq_out_name,
							 $CTL_OPTIONS{'model_org'},
							 $CTL_OPTIONS{'RepeatMasker'},
							 '',
							 $CTL_OPTIONS{'cpus'},
							 $opt_f,
							 $self->{LOG}
							);
	 push(@{$rma_keepers}, @{$rma_keepers1});

	 #-mask the chunk using repeatmasker hits
	 $chunk = repeat_mask_seq::mask_chunk($chunk, $rma_keepers1);
      }

      #-mask species specific repeats;
      if($CTL_OPTIONS{rmlib}){
	 my $rma_keepers2 = Shared_Functions::repeatmask($chunk,
							 $the_void,
							 $seq_out_name,
							 $CTL_OPTIONS{'model_org'},
							 $CTL_OPTIONS{'RepeatMasker'},
							 $CTL_OPTIONS{'rmlib'},
							 $CTL_OPTIONS{'cpus'},
							 $opt_f,
							 $self->{LOG}
							);
	 push(@{$rma_keepers}, @{$rma_keepers2});

	 #-mask the chunk using repeatmasker hits
	 $chunk = repeat_mask_seq::mask_chunk($chunk, $rma_keepers2);
      }
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($chunk, $rma_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 1) {
      #------------------------ARGS_IN
      my $chunk             = shift @{$vars};
      my $repeat_protein    = shift @{$vars};
      my $the_void          = shift @{$vars};
      my $seq_out_name      = shift @{$vars};
      my %CTL_OPTIONS       = %{shift @{$vars}};
      my $opt_f             = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- blastx against a repeat library (for better masking)
      my $rep_blastx_res_dir = '';
       $rep_blastx_res_dir = Shared_Functions::blastx_as_chunks($chunk,
								$repeat_protein,
								$the_void,
								$seq_out_name,
								$CTL_OPTIONS{blastx},
								$CTL_OPTIONS{eval_blastx},
								$CTL_OPTIONS{split_hit},
								$CTL_OPTIONS{cpus},
								$CTL_OPTIONS{old_repeat_protein},
								$CTL_OPTIONS{xdformat},
								$CTL_OPTIONS{alt_peptide},
								$self->{RANK},
								$opt_f,
								$self->{LOG},
								$self->{LOG_FLAG}
							       ) if ($CTL_OPTIONS{old_repeat_protein});
      #-------------------------CHUNK
	    
      #------------------------RESULTS
      @results = ($rep_blastx_res_dir);
      #------------------------RESULTS
   }
   elsif ($level == 2) {
      #------------------------ARGS_IN
      my $chunk              = shift @{$vars};
      my $rep_blastx_res_dir = shift @{$vars};
      my %CTL_OPTIONS        = %{shift @{$vars}};
      my $opt_f              = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- merge and collect blastx results
      my $repeat_blastx_keepers = [];
      $repeat_blastx_keepers = Shared_Functions::collect_blastx($chunk,
								$rep_blastx_res_dir,
								$CTL_OPTIONS{eval_blastx},
								$CTL_OPTIONS{bit_blastx},
								$CTL_OPTIONS{percov_blastx},
								$CTL_OPTIONS{percid_blastx},
								$CTL_OPTIONS{split_hit},
								$opt_f,
								$self->{LOG}
							       ) if($CTL_OPTIONS{old_repeat_protein});
      #-------------------------CHUNK
	
      #------------------------RESULTS
      @results = ($repeat_blastx_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 3) {
      #------------------------ARGS_IN
      my $chunk                 = shift @{$vars};
      my $rma_keepers           = shift @{$vars};
      my $repeat_blastx_keepers = shift @{$vars};
      my $GFF3                  = shift @{$vars};
      my $query_def             = shift @{$vars};
      my $query_seq             = shift @{$vars};
      my $masked_total_seq      = shift @{$vars};
      my $the_void              = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-mask the chunk using blastx hits
      $chunk = repeat_mask_seq::mask_chunk($chunk, $repeat_blastx_keepers);
      
      #-combine and cluster blastx and repeatmasker hits
      #-to get consensus repeat hits for gff3 and XML annotations
      my $rm_keepers = repeat_mask_seq::process($rma_keepers, 
						$repeat_blastx_keepers,
						$query_seq
					       );
	 
      #-add repeats to GFF3
      $GFF3->repeat_hits($rm_keepers);

      #-build/fill big masked sequence
      $masked_total_seq .= $chunk->seq();
      #-------------------------CHUNK
	    
      #------------------------RESULTS
      @results = ($chunk, $masked_total_seq, $GFF3);
      #------------------------RESULTS
   }
   elsif ($level == 4) {
      #------------------------ARGS_IN
      my $masked_total_seq = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my $query_def    = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};      
      #------------------------ARGS_IN

      #-------------------------CHUNK
      my $masked_fasta = Fasta::toFasta($query_def.' masked', \$masked_total_seq);
      FastaFile::writeFile($masked_fasta ,$the_void."/query.masked.fasta");
      
      #==SNAP ab initio here
      my $snaps = [];
      $snaps = Shared_Functions::snap($masked_fasta,
				      $the_void,
				      $seq_out_name,
				      $CTL_OPTIONS{snap},
				      $CTL_OPTIONS{'snaphmm'},
				      $opt_f,
				      $self->{LOG}
				      ) if ($CTL_OPTIONS{'snap'});

      #==AUGUSTUS ab initio here
      my $augus = [];
      $augus = Shared_Functions::augustus($masked_fasta,
					  $the_void,
					  $seq_out_name,
					  $CTL_OPTIONS{'augustus'},
					  $CTL_OPTIONS{'augustus_species'},
					  $opt_f,
					  $self->{LOG}
					 ) if ($CTL_OPTIONS{'augustus'});

      #-- build an index of the databases
      my $t_dir = $TMP."/rank".$self->{RANK};
      File::Path::mkpath($t_dir);

      $CTL_OPTIONS{old_est} =~ /([^\/]+)$/;
      my $t_name = $1;

      $CTL_OPTIONS{old_protein} =~ /([^\/]+)$/;
      my $p_name = $1;

      my $a_name;
      if($CTL_OPTIONS{old_alt_est}){
	  $CTL_OPTIONS{old_alt_est} =~ /([^\/]+)$/;  
	  $a_name = $1;;
      }

      my $trans_file = $t_dir."/".$t_name;
      my $prot_file = $t_dir."/".$p_name;
      my $alt_est_file = $t_dir."/".$a_name if($CTL_OPTIONS{old_alt_est});

      if (! -e $trans_file) {
	  system("cp $CTL_OPTIONS{old_est} $trans_file");
      }
      if (! -e $prot_file) {
          system("cp $CTL_OPTIONS{old_protein} $prot_file");
      }
      if ($CTL_OPTIONS{old_alt_est} && ! -e $alt_est_file) {
          system("cp $CTL_OPTIONS{old_alt_est} $alt_est_file");
      }

      my $fasta_t_index     = Shared_Functions::build_fasta_index($trans_file);
      my $fasta_p_index     = Shared_Functions::build_fasta_index($prot_file);
      my $fasta_a_index;
      $fasta_a_index = Shared_Functions::build_fasta_index($alt_est_file) if ($CTL_OPTIONS{old_alt_est});

      #--set up new chunks for remaining levels
      my $fasta_chunker = new FastaChunker();
      $fasta_chunker->parent_fasta($$masked_fasta);
      $fasta_chunker->chunk_size($CTL_OPTIONS{'max_dna_len'});
      $fasta_chunker->min_size($CTL_OPTIONS{'split_hit'});
      $fasta_chunker->load_chunks();
      
      my $chunk_count = 0;
      
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($masked_fasta, $snaps, $augus, $fasta_chunker, $chunk_count, $fasta_t_index, $fasta_p_index, $fasta_a_index);
      #------------------------RESULTS
   }
   elsif ($level == 5) {
      #------------------------ARGS_IN
      my $chunk        = shift @{$vars};
      my $transcripts  = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- blastn search the file against ESTs
      my $blastn_res_dir = Shared_Functions::blastn_as_chunks($chunk,
							      $transcripts,
							      $the_void,
							      $seq_out_name,
							      $CTL_OPTIONS{blastn},
							      $CTL_OPTIONS{eval_blastn},
							      $CTL_OPTIONS{split_hit},
							      $CTL_OPTIONS{cpus},
							      $CTL_OPTIONS{old_est},
							      $CTL_OPTIONS{xdformat},
							      $self->{RANK},
							      $opt_f,
							      $self->{LOG},
							      $self->{LOG_FLAG}
							     );
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastn_res_dir);
      #------------------------RESULTS
   }
   elsif ($level == 6) {
      #------------------------ARGS_IN
      my $chunk          = shift @{$vars};
      my $blastn_res_dir = shift @{$vars};
      my %CTL_OPTIONS    = %{shift @{$vars}};
      my $opt_f          = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- merge and collect blastn results
      my $blastn_keepers = Shared_Functions::collect_blastn($chunk, 
							    $blastn_res_dir,
							    $CTL_OPTIONS{eval_blastn},
							    $CTL_OPTIONS{bit_blastn},
							    $CTL_OPTIONS{percov_blastn},
							    $CTL_OPTIONS{percid_blastn},
							    $CTL_OPTIONS{split_hit},
							    $opt_f,
							    $self->{LOG}
							   );      
      #-------------------------CHUNK
	
      #------------------------RESULTS
      @results = ($blastn_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 7) {
      #------------------------ARGS_IN
      my $chunk        = shift @{$vars};
      my $proteins     = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- blastx search  the masked input file
      my $blastx_res_dir = Shared_Functions::blastx_as_chunks($chunk,
							      $proteins,
							      $the_void,
							      $seq_out_name,
							      $CTL_OPTIONS{blastx},
							      $CTL_OPTIONS{eval_blastx},
							      $CTL_OPTIONS{split_hit},
							      $CTL_OPTIONS{cpus},
							      $CTL_OPTIONS{old_protein},
							      $CTL_OPTIONS{xdformat},
							      $CTL_OPTIONS{alt_peptide},
							      $self->{RANK},
							      $opt_f,
							      $self->{LOG},
							      $self->{LOG_FLAG}
							     );
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastx_res_dir);
      #------------------------RESULTS
   }
   elsif ($level == 8) {
      #------------------------ARGS_IN
      my $chunk          = shift @{$vars};
      my $blastx_res_dir = shift @{$vars};
      my %CTL_OPTIONS    = %{shift @{$vars}};
      my $opt_f          = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- merge and collect blastx results
      my $blastx_keepers = Shared_Functions::collect_blastx($chunk,
							    $blastx_res_dir,
							    $CTL_OPTIONS{eval_blastx},
							    $CTL_OPTIONS{bit_blastx},
							    $CTL_OPTIONS{percov_blastx},
							    $CTL_OPTIONS{percid_blastx},
							    $CTL_OPTIONS{split_hit},
							    $opt_f,
							    $self->{LOG}
							   );
      #-------------------------CHUNK
	
      #------------------------RESULTS
      @results = ($blastx_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 9) {
      #------------------------ARGS_IN
      my $chunk        = shift @{$vars};
      my $alt_ests     = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- blastx search  the masked input file
      my $tblastx_res_dir = '';
      $tblastx_res_dir = Shared_Functions::tblastx_as_chunks($chunk,
							     $alt_ests,
							     $the_void,
							     $seq_out_name,
							     $CTL_OPTIONS{tblastx},
							     $CTL_OPTIONS{eval_tblastx},
							     $CTL_OPTIONS{split_hit},
							     $CTL_OPTIONS{cpus},
							     $CTL_OPTIONS{old_alt_est},
							     $CTL_OPTIONS{xdformat},
							     $self->{RANK},
							     $opt_f,
							     $self->{LOG},
							     $self->{LOG_FLAG}
							    ) if($alt_ests);
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($tblastx_res_dir);
      #------------------------RESULTS
   }
   elsif ($level == 10) {
      #------------------------ARGS_IN
      my $chunk           = shift @{$vars};
      my $tblastx_res_dir = shift @{$vars};
      my %CTL_OPTIONS     = %{shift @{$vars}};
      my $opt_f           = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- merge and collect blastx results
      my $tblastx_keepers = [];
      $tblastx_keepers = Shared_Functions::collect_tblastx($chunk,
							   $tblastx_res_dir,
							   $CTL_OPTIONS{eval_tblastx},
							   $CTL_OPTIONS{bit_tblastx},
							   $CTL_OPTIONS{percov_tblastx},
							   $CTL_OPTIONS{percid_tblastx},
							   $CTL_OPTIONS{split_hit},
							   $opt_f,
							   $self->{LOG}
							  ) if($CTL_OPTIONS{old_alt_est});
      #-------------------------CHUNK
	
      #------------------------RESULTS
      @results = ($tblastx_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 11) {
      #------------------------ARGS_IN
      my $chunk            = shift @{$vars};
      my $masked_fasta     = shift @{$vars};
      my $snaps            = shift @{$vars};
      my $augus            = shift @{$vars};
      my $blastn_keepers   = shift @{$vars};
      my $blastx_keepers   = shift @{$vars};
      my $tblastx_keepers  = shift @{$vars};
      my $fasta_t_index    = shift @{$vars};
      my $fasta_p_index    = shift @{$vars};
      my $fasta_a_index    = shift @{$vars};
      my $holdover_blastn  = shift @{$vars};
      my $holdover_blastx  = shift @{$vars};
      my $holdover_tblastx = shift @{$vars};
      my $the_void         = shift @{$vars};
      my %CTL_OPTIONS      = %{shift @{$vars}};
      my $opt_f            = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK

      #-- decide which gene finder to use to build annotations 
      my $preds_on_chunk;

      if ($CTL_OPTIONS{predictor} eq 'augustus') {
	 $preds_on_chunk = Shared_Functions::get_preds_on_chunk($augus,
								$chunk
							       );
      }
      elsif ($CTL_OPTIONS{predictor} eq 'snap') {
	 $preds_on_chunk = Shared_Functions::get_preds_on_chunk($snaps,
								$chunk
							       );
      }
      elsif ($CTL_OPTIONS{predictor} eq 'est2genome') {
	 $preds_on_chunk = [];
      }
      else {
	 die "ERROR: invalid predictor type: $CTL_OPTIONS{predictor}\n";
      }
      
      #==merge heldover Phathits from last round
      if ($chunk->number != 0) { #if not first chunk
	 ($blastn_keepers,
	  $blastx_keepers,
	  $tblastx_keepers) = Shared_Functions::merge_and_resolve_hits($masked_fasta,
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
								       \%CTL_OPTIONS,
								       $opt_f,
								       $self->{LOG}
								      );
      }
      
      #==PROCESS HITS CLOSE TOO CHUNK DIVISIONS 
      my $holdover_preds = [];
      $holdover_blastn = [];
      $holdover_blastx = [];
      $holdover_tblastx = [];
  
      if (not $chunk->is_last) { #if not last chunk
	 ($holdover_blastn,
	  $holdover_blastx,
	  $holdover_tblastx,
	  $holdover_preds,
	  $blastn_keepers,
	  $blastx_keepers,
	  $tblastx_keepers,
	  $preds_on_chunk) = Shared_Functions::process_the_chunk_divide($chunk,
									$CTL_OPTIONS{'split_hit'},
									$blastn_keepers,
									$blastx_keepers,
									$tblastx_keepers,
									$preds_on_chunk
								       );
      }
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($holdover_blastn, $holdover_blastx, $holdover_tblastx, $holdover_preds,
		  $blastn_keepers, $blastx_keepers, $tblastx_keepers, $preds_on_chunk);
      #------------------------RESULTS
   }
   elsif ($level == 12) {
      #------------------------ARGS_IN
      my $fasta           = shift @{$vars};
      my $blastx_keepers  = shift @{$vars};
      my $query_seq       = shift @{$vars};
      my $fasta_p_index   = shift @{$vars};
      my $the_void        = shift @{$vars};
      my %CTL_OPTIONS     = %{shift @{$vars}};
      my $opt_f           = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-cluster the blastx hits
      print STDERR "cleaning blastx...\n" unless($main::quiet);
      my $blastx_clusters = cluster::clean_and_cluster($blastx_keepers,
						       $query_seq,
						       10);

      #-- make a multi-fasta of the seqs in the blastx_clusters 
      #-- polish the blastx hits with exonerate
      my $exonerate_p_clusters = Shared_Functions::polish_exonerate($fasta,
								    $blastx_clusters,
								    $fasta_p_index,
								    $the_void,
								    5,
								    'p',
								    $CTL_OPTIONS{exonerate},
								    $CTL_OPTIONS{percov_blastx},
								    $CTL_OPTIONS{percid_blastx},
								    $CTL_OPTIONS{ep_score_limit},
								    $CTL_OPTIONS{ep_matrix},
								    $opt_f,
								    $self->{LOG}
								   );
      
      my $blastx_data      = Shared_Functions::flatten($blastx_clusters);
      my $exonerate_p_data = Shared_Functions::flatten($exonerate_p_clusters, 'exonerate:p');
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastx_data, $exonerate_p_data);
      #------------------------RESULTS
   }
   elsif ($level == 13) {
      #------------------------ARGS_IN
      my $fasta           = shift @{$vars};
      my $blastn_keepers  = shift @{$vars};
      my $tblastx_keepers = shift @{$vars};
      my $query_seq       = shift @{$vars};
      my $fasta_t_index   = shift @{$vars};
      my $the_void        = shift @{$vars};
      my %CTL_OPTIONS     = %{shift @{$vars}};
      my $opt_f           = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-cluster the tblastx hits
      print STDERR "cleaning tblastx...\n" unless $main::quiet;
      my $tblastx_clusters = cluster::clean_and_cluster($tblastx_keepers,
							$query_seq,
							10
						       );

      undef $tblastx_keepers; #free up memory
      my $tblastx_data      = Shared_Functions::flatten($tblastx_clusters);

      #-- Cluster the blastn hits
      print STDERR "cleaning blastn...\n" unless($main::quiet);
      my $blastn_clusters = cluster::clean_and_cluster($blastn_keepers,
						       $query_seq,
						       10
						      );

      #-- polish blastn hits with exonerate
      my $exonerate_e_clusters = Shared_Functions::polish_exonerate($fasta,
								    $blastn_clusters,
								    $fasta_t_index,
								    $the_void,
								    5,
								    'e',
								    $CTL_OPTIONS{exonerate},
								    $CTL_OPTIONS{percov_blastn},
								    $CTL_OPTIONS{percid_blastn},
								    $CTL_OPTIONS{en_score_limit},
								    $CTL_OPTIONS{en_matrix},
								    $opt_f,
								    $self->{LOG}
								   );

      my $blastn_data      = Shared_Functions::flatten($blastn_clusters);
      my $exonerate_e_data = Shared_Functions::flatten($exonerate_e_clusters, 'exonerate:e');
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastn_data, $tblastx_data, $exonerate_e_data);
      #------------------------RESULTS
   }
   elsif ($level == 14) {
      #------------------------ARGS_IN
      my $fasta            = shift @{$vars};
      my $masked_fasta     = shift @{$vars};
      my $c_number         = shift @{$vars};
      my $exonerate_p_data = shift @{$vars};
      my $exonerate_e_data = shift @{$vars};
      my $blastx_data      = shift @{$vars};
      my $preds_on_chunk   = shift @{$vars};
      my $the_void         = shift @{$vars};
      my %CTL_OPTIONS      = %{shift @{$vars}};
      my $opt_f            = shift @{$vars};
      my $opt_preds        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #==MAKER annotations built here

      #-- decide which gene finder to use to build annotations 
      my $pred_command;

      if ($CTL_OPTIONS{predictor} eq 'augustus') {
	 $pred_command = $CTL_OPTIONS{augustus} .' --species='.$CTL_OPTIONS{augustus_species};
      }
      elsif ($CTL_OPTIONS{predictor} eq 'snap') {
	 $pred_command = $CTL_OPTIONS{snap}.' '.$CTL_OPTIONS{snaphmm};
      }
      elsif ($CTL_OPTIONS{predictor} eq 'est2genome') {
	 $pred_command = '';
      }
      else {
	 die "ERROR: invalid predictor type: $CTL_OPTIONS{predictor}\n";
      }
      
      #-auto-annotate the input file

      my $annotations = maker::auto_annotator::annotate($fasta,
							$$masked_fasta,
							$c_number,
							$exonerate_p_data,
							$exonerate_e_data,
							$blastx_data,
							$preds_on_chunk,
							$the_void,
							$pred_command,
							$CTL_OPTIONS{'snap_flank'},
							$CTL_OPTIONS{'single_exon'},
							$opt_f,
							$opt_preds,
							$CTL_OPTIONS{predictor},
							$self->{LOG},
							\%CTL_OPTIONS,
						       );
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($annotations);
      #------------------------RESULTS
   }
   elsif ($level == 15) {
      #------------------------ARGS_IN
      my $blastx_data      = shift @{$vars};
      my $blastn_data      = shift @{$vars};
      my $tblastx_data     = shift @{$vars};
      my $exonerate_p_data = shift @{$vars};
      my $exonerate_e_data = shift @{$vars};
      my $annotations      = shift @{$vars};
      my $p_fastas         = shift @{$vars};
      my $t_fastas         = shift @{$vars};
      my $GFF3             = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #--- GFF3      
      $GFF3->genes($annotations);
      $GFF3->phat_hits($blastx_data);
      $GFF3->phat_hits($blastn_data);
      $GFF3->phat_hits($tblastx_data);
      $GFF3->phat_hits($exonerate_p_data);
      $GFF3->phat_hits($exonerate_e_data);
	
      #--- building fastas for annotations
      my ($p_fasta, $t_fasta) = Shared_Functions::get_maker_p_and_t_fastas($annotations);
      $p_fastas .= $p_fasta;
      $t_fastas .= $t_fasta;
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($GFF3, $p_fastas, $t_fastas);
      #------------------------RESULTS
   }
   elsif ($level == 16) {
      #------------------------ARGS_IN
      my $snaps          = shift @{$vars};
      my $augus          = shift @{$vars};
      my $p_fastas       = shift @{$vars};
      my $t_fastas       = shift @{$vars};
      my $GFF3           = shift @{$vars};
      my $seq_out_name   = shift @{$vars};
      my $out_dir        = shift @{$vars};
      my $the_void       = shift @{$vars};
      my $query_seq      = shift @{$vars};
      my %CTL_OPTIONS    = %{shift @{$vars}};
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      #--- building fastas of predictions
      my ($p_snap_fastas, $t_snap_fastas) = Shared_Functions::get_snap_p_and_t_fastas($query_seq, $snaps);
      my ($p_augus_fastas, $t_augus_fastas) = Shared_Functions::get_snap_p_and_t_fastas($query_seq, $augus);

      #--Write fasta files and gff3 files now that all chunks are finished
      FastaFile::writeFile(\$p_fastas ,"$out_dir\/$seq_out_name\.maker.proteins.fasta");
      FastaFile::writeFile(\$t_fastas ,"$out_dir\/$seq_out_name\.maker.transcripts.fasta");
      if ($CTL_OPTIONS{'snap'}) {
	  FastaFile::writeFile(\$p_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.proteins.fasta");
	  FastaFile::writeFile(\$t_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.transcript.fasta");
      }
      if ($CTL_OPTIONS{'augustus'}) {
	 FastaFile::writeFile(\$p_augus_fastas ,"$out_dir\/$seq_out_name\.maker.augus.proteins.fasta");
	 FastaFile::writeFile(\$t_augus_fastas ,"$out_dir\/$seq_out_name\.maker.augus.transcript.fasta");
      }
      $GFF3->predictions($snaps);
      $GFF3->predictions($augus);
      $GFF3->print($out_dir."/".$seq_out_name.".gff");
      
      #--cleanup maker files created with each fasta sequence
      File::Path::rmtree ($the_void) if $CTL_OPTIONS{clean_up};	#rm temp directory
      #-------------------------CHUNK
      
      #------------------------RESULTS
      @results = ();
      #------------------------RESULTS
   }
   else {
      warn "Error: Invalid argument for method run() in Process::MakerChunk\n";
      return undef;
   }
   
   $self->{RESULT} = \@results;
}

#--------------------------------------------------------------
sub result {
   my $self = shift;
    
   return $self->{RESULT} || [];
}
#--------------------------------------------------------------
sub failed {
   my $self = shift;
    
   return $self->{FAILED} || undef;
}
#--------------------------------------------------------------
sub exception {
   my $self = shift;
    
   return $self->{EXCEPTION} || undef;
}
#--------------------------------------------------------------
sub error {
   my $self = shift;
    
   return $self->{ERROR} || '';
}

#--------------------------------------------------------------
sub id {
   my $self = shift;
   my $arg = shift;

   return $self->{ID};
}

#--------------------------------------------------------------
sub clone {
   my $self = shift;
    
   my $clone = dclone($self);

   return $clone;
}

#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
1;
