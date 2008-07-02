#! /usr/bin/perl -w

package Process::MakerChunk;

use strict;
use Storable qw (freeze thaw dclone);

use FindBin;
use lib "$FindBin::Bin/../..";

use File::Util;
use File::Temp qw(tempfile);
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
   $self->{VARS}  = shift;	#this should be an array reference
   $self->{ID}    = shift || 0;
   $self->{RANK}  = shift || 0;
}

#--------------------------------------------------------------
sub run {
   my $self = shift;
   $self->{RANK} = shift || $self->{RANK} || 0;

   if (exists $self->{RESULT}) {
      return;
   }

   my $level = $self->{LEVEL};
   my $vars = $self->{VARS};
   my @results;

   #--redirect STDERR to a log file
   open (OLDERR, ">&STDERR");
   my (undef, $t_name) = tempfile();
   close(STDERR);
   open(STDERR, "| tee $t_name >&2");
   select((select(STDERR), $|=1)[0]);

   if ($level == 0) {
      #------------------------ARGS_IN
      my $chunk        = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      my $opt_f        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- repeatmask the input file
      $chunk->seq(uc($chunk->seq())); #must be upper case before soft masking

      my $rma_keepers = repeatmask($chunk, 
				   $the_void,
				   $seq_out_name,
				   $CTL_OPTIONS{'model_org'},
				   $CTL_OPTIONS{'RepeatMasker'},
				   $CTL_OPTIONS{'rmlib'},
				   $CTL_OPTIONS{'cpus'},
				   $opt_f
				  );

      #-mask the chunk using repeatmasker hits
      $chunk = repeat_mask_seq::mask_chunk($chunk, $rma_keepers);
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
      my $repeat_blastx_keepers = blastx($chunk,
					 $repeat_protein,
					 $the_void,
					 $seq_out_name,
					 $CTL_OPTIONS{blastx},
					 $CTL_OPTIONS{eval_blastx},
					 $CTL_OPTIONS{bit_blastx},
					 $CTL_OPTIONS{percov_blastx},
					 $CTL_OPTIONS{percid_blastx},
					 $CTL_OPTIONS{split_hit},
					 $CTL_OPTIONS{cpus},
					 $CTL_OPTIONS{old_repeat_protein},
					 $CTL_OPTIONS{xdformat},
					 $self->id(),
					 $self->{RANK},
					 $opt_f
					);
      #-------------------------CHUNK
	    
      #------------------------RESULTS
      @results = ($repeat_blastx_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 2) {
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
      #merge blast reports
      my($f_util) = File::Util->new();
      my(@dirs) = $f_util->list_dir($the_void, '--dirs-only');
      @dirs = grep (/\.blastx\.temp_dir$/, @dirs);
      
      foreach my $dir (@dirs) {
	 my $blast_finished = $dir;
	 $blast_finished =~ s/\.temp_dir$//;
	 system ("cat $the_void/$dir/*.blastx > $the_void/$blast_finished");
	 rmtree ("$the_void/$dir");
      }

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
   elsif ($level == 3) {
      #------------------------ARGS_IN
      my $masked_total_seq = shift @{$vars};
      my $the_void     = shift @{$vars};
      my $seq_out_name = shift @{$vars};
      my $query_def    = shift @{$vars};
      my %CTL_OPTIONS  = %{shift @{$vars}};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      my $masked_fasta = Fasta::toFasta($query_def.' masked', \$masked_total_seq);
      FastaFile::writeFile($masked_fasta ,$the_void."/query.masked.fasta");
      
      #==SNAP ab initio here
      my $snaps = snap($masked_fasta,
		       $the_void,
		       $seq_out_name,
		       $CTL_OPTIONS{snap},
		       $CTL_OPTIONS{'snaphmm'}
		      );

      #--set up new chunks for remaining levels
      my $fasta_chunker = new FastaChunker();
      $fasta_chunker->parent_fasta($$masked_fasta);
      $fasta_chunker->chunk_size($CTL_OPTIONS{'max_dna_len'});
      $fasta_chunker->min_size($CTL_OPTIONS{'split_hit'});
      $fasta_chunker->load_chunks();
      
      my $chunk_count = 0;
      
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($masked_fasta, $snaps, $fasta_chunker, $chunk_count);
      #------------------------RESULTS
   }
      elsif ($level == 4) {
      #------------------------ARGS_IN
      my $holdover_chunk = shift @{$vars};
      my $chunk          = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-merge heldover chunk sequence from last round
      unless($chunk->number == 0){
	 $chunk = merge_chunks($holdover_chunk, $chunk);
      } 
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($chunk);
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
      my $blastn_keepers = blastn($chunk,
				  $transcripts,
				  $the_void,
				  $seq_out_name,
				  $CTL_OPTIONS{blastn},
				  $CTL_OPTIONS{eval_blastn},
				  $CTL_OPTIONS{bit_blastn},
				  $CTL_OPTIONS{percov_blastn},
				  $CTL_OPTIONS{percid_blastn},
				  $CTL_OPTIONS{split_hit},
				  $CTL_OPTIONS{cpus},
				  $CTL_OPTIONS{old_est},
				  $CTL_OPTIONS{xdformat},
				  $self->id(),
				  $self->{RANK},
				  $opt_f
				 );
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastn_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 6) {
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
      my $blastx_keepers = blastx($chunk, 
				  $proteins,
				  $the_void,
				  $seq_out_name,
				  $CTL_OPTIONS{blastx},
				  $CTL_OPTIONS{eval_blastx},
				  $CTL_OPTIONS{bit_blastx},
				  $CTL_OPTIONS{percov_blastx},
				  $CTL_OPTIONS{percid_blastx},
				  $CTL_OPTIONS{split_hit},
				  $CTL_OPTIONS{cpus},
				  $CTL_OPTIONS{old_protein},
				  $CTL_OPTIONS{xdformat},
				  $self->id(),
				  $self->{RANK},
				  $opt_f
				 );
      #-------------------------CHUNK
	
      #------------------------RESULTS
      @results = ($blastx_keepers);
      #------------------------RESULTS
   }
   elsif ($level == 7) {
      #------------------------ARGS_IN
      my $chunk          = shift @{$vars};
      my $holdover_chunk = shift @{$vars};
      my $snaps          = shift @{$vars};
      my $blastx_keepers = shift @{$vars};
      my $blastn_keepers = shift @{$vars};
      my $query_seq      = shift @{$vars};
      my $split_hit      = shift @{$vars};
      my $the_void       = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #merge blast reports
      my($f_util) = File::Util->new();
      my(@dirs) = $f_util->list_dir($the_void, '--dirs-only');
      @dirs = grep (/blast.\.temp_dir$/, @dirs);

      foreach my $dir (@dirs) {
	 my $blast_finished = $dir;
	 $blast_finished =~ s/\.temp_dir$//;
	 system ("cat $the_void/$dir/*blast* > $the_void/$blast_finished");
	 rmtree ("$the_void/$dir");
      }
    
      #==get just the snaps that overlap this chunk
      my $snaps_on_chunk = get_snaps_on_chunk($snaps,
					      $chunk
					     );

      #==PROCESS HITS CLOSE TOO CHUNK DIVISIONS
      if(not $chunk->is_last){ #if not last chunk
	 ($holdover_chunk,
	  $blastx_keepers,
	  $blastn_keepers,
	  $snaps_on_chunk) = process_the_chunk_divide($chunk,
						      $split_hit,
						      $blastx_keepers,
						      $blastn_keepers,
						      $snaps_on_chunk
						     );
      }
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($holdover_chunk, $blastx_keepers, $blastn_keepers, $snaps_on_chunk);
      #------------------------RESULTS
   }
   elsif ($level == 8) {
      #------------------------ARGS_IN
      my $fasta           = shift @{$vars};
      my $blastx_keepers  = shift @{$vars};
      my $query_seq       = shift @{$vars};
      my $the_void        = shift @{$vars};
      my %CTL_OPTIONS     = %{shift @{$vars}};
      my $opt_f           = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- build an index of the databases
      my $fasta_p_index     = build_fasta_index($CTL_OPTIONS{old_protein});

      #-cluster the blastx hits
      print STDERR "cleaning blastx...\n";
      my $blastx_clusters = cluster::clean_and_cluster($blastx_keepers,
						       $query_seq,
						       10);

      #-- make a multi-fasta of the seqs in the blastx_clusters 
      #-- polish the blastx hits with exonerate
      my $exonerate_p_clusters = polish_exonerate($fasta,
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
						  $opt_f
						 );

      my $blastx_data      = flatten($blastx_clusters);
      my $exonerate_p_data = flatten($exonerate_p_clusters, 'exonerate:p');
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastx_data, $exonerate_p_data);
      #------------------------RESULTS
   }
   elsif ($level == 9) {
      #------------------------ARGS_IN
      my $fasta           = shift @{$vars};
      my $blastn_keepers = shift @{$vars};
      my $query_seq       = shift @{$vars};
      my $the_void        = shift @{$vars};
      my %CTL_OPTIONS     = %{shift @{$vars}};
      my $opt_f           = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #-- build an index of the databases
      my $fasta_t_index     = build_fasta_index($CTL_OPTIONS{old_est});

      #-- Cluster the blastn hits
      print STDERR "cleaning blastn...\n";
      my $blastn_clusters = cluster::clean_and_cluster($blastn_keepers,
						       $query_seq,
						       10);

      #-- polish blastn hits with exonerate
      my $exonerate_e_clusters = polish_exonerate($fasta,
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
						  $opt_f
						 );

      my $blastn_data      = flatten($blastn_clusters);
      my $exonerate_e_data = flatten($exonerate_e_clusters, 'exonerate:e');
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($blastn_data, $exonerate_e_data);
      #------------------------RESULTS
   }
   elsif ($level == 10) {
      #------------------------ARGS_IN
      my $fasta            = shift @{$vars};
      my $masked_fasta     = shift @{$vars};
      my $c_number         = shift @{$vars};
      my $exonerate_p_data = shift @{$vars};
      my $exonerate_e_data = shift @{$vars};
      my $blastx_data      = shift @{$vars};
      my $snaps_on_chunk   = shift @{$vars};
      my $the_void         = shift @{$vars};
      my %CTL_OPTIONS      = %{shift @{$vars}};
      my $opt_f            = shift @{$vars};
      my $opt_snaps        = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #==MAKER annotations built here
      
      #-auto-annotate the input file
      my $snap_command = $CTL_OPTIONS{snap}.' '.$CTL_OPTIONS{snaphmm};

      my $annotations = maker::auto_annotator::annotate($fasta,
							$$masked_fasta,
							$c_number,
							$exonerate_p_data,
							$exonerate_e_data,
							$blastx_data,
							$snaps_on_chunk,
							$the_void,
							$snap_command,
							$CTL_OPTIONS{'snap_flank'},
							$CTL_OPTIONS{'single_exon'},
							$opt_f,
							$opt_snaps,
							'snap'#temporary fix
						       );
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($annotations);
      #------------------------RESULTS
   }
   elsif ($level == 11) {
      #------------------------ARGS_IN
      my $blastx_data      = shift @{$vars};
      my $blastn_data      = shift @{$vars};
      my $exonerate_p_data = shift @{$vars};
      my $exonerate_e_data = shift @{$vars};
      my $annotations      = shift @{$vars};
      my $snaps            = shift @{$vars};
      my $query_seq        = shift @{$vars};
      my $snaps_on_chunk   = shift @{$vars};
      my $p_fastas         = shift @{$vars};
      my $t_fastas         = shift @{$vars};
      my $p_snap_fastas    = shift @{$vars};
      my $t_snap_fastas    = shift @{$vars};
      my $GFF3             = shift @{$vars};
      #------------------------ARGS_IN

      #-------------------------CHUNK
      #--- GFF3      
      $GFF3->genes($annotations);
      $GFF3->predictions($snaps_on_chunk);
      $GFF3->phat_hits($blastx_data);
      $GFF3->phat_hits($blastn_data);
      $GFF3->phat_hits($exonerate_p_data);
      $GFF3->phat_hits($exonerate_e_data);
	
      #--- building fastas
      my ($p_fasta, $t_fasta) = get_maker_p_and_t_fastas($annotations);
      $p_fastas .= $p_fasta;
      $t_fastas .= $t_fasta;
   
      my ($p_snap_fasta, $t_snap_fasta) = get_snap_p_and_t_fastas($query_seq, $snaps_on_chunk);
      $p_snap_fastas .= $p_snap_fasta;
      $t_snap_fastas .= $t_snap_fasta;
      #-------------------------CHUNK

      #------------------------RESULTS
      @results = ($GFF3, $p_fastas, $t_fastas, $p_snap_fastas, $t_snap_fastas);
      #------------------------RESULTS
   }
   elsif ($level == 12) {
      #------------------------ARGS_IN
      my $p_fastas      = shift @{$vars};
      my $t_fastas      = shift @{$vars};
      my $p_snap_fastas = shift @{$vars};
      my $t_snap_fastas = shift @{$vars};
      my $GFF3          = shift @{$vars};
      my $seq_out_name  = shift @{$vars};
      my $out_dir       = shift @{$vars};
      my $the_void      = shift @{$vars};
      my %CTL_OPTIONS   = %{shift @{$vars}};
      #------------------------ARGS_IN
      
      #-------------------------CHUNK
      #--Write fasta files and gff3 files now that all chunks are finished
      FastaFile::writeFile(\$p_fastas ,"$out_dir\/$seq_out_name\.maker.proteins.fasta");
      FastaFile::writeFile(\$t_fastas ,"$out_dir\/$seq_out_name\.maker.transcripts.fasta");
      FastaFile::writeFile(\$p_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.proteins.fasta");
      FastaFile::writeFile(\$t_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.transcript.fasta");
      $GFF3->print($out_dir."/".$seq_out_name.".gff");
      
      #--cleanup maker files created with each fasta sequence
      rmtree ($the_void) if $CTL_OPTIONS{clean_up}; #rm temp directory
      #-------------------------CHUNK
      
      #------------------------RESULTS
      @results = ();
      #------------------------RESULTS
   }
   else {
      warn "Error: Invalid argument for method run() in Process::MakerChunk\n";
      return undef;
   }
   
   #--redirect STDERR back to STDERR
   close(STDERR);
   open (STDERR, ">&OLDERR");
   close(OLDERR);

   #--collect STDERR log file data
   open (IN, "< $t_name");
   $self->{ERROR} = join('', <IN>);
   close(IN);
   $self->{RESULT} = \@results;
}

#--------------------------------------------------------------
sub result {
   my $self = shift;
    
   return $self->{RESULT};
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
sub get_snaps_on_chunk {
   my $snaps = shift;
   my $chunk = shift;

   my $c_pstart = $chunk->offset + 1;
   my $c_pend = $chunk->offset + $chunk->length;
   
   my $c_mstart = $chunk->offset + 1;
   my $c_mend = $chunk->offset + $chunk->length;

   if($chunk->p_cutoff || $chunk->m_cutoff){
      $c_pstart = $chunk->p_cutoff;
      $c_pend = $chunk->p_cutoff + $chunk->length - 1;
      
      $c_mstart = $chunk->m_cutoff;
      $c_mend = $chunk->m_cutoff + $chunk->length - 1;
   }

   my @keepers;
   foreach my $snap (@{$snaps}){
      my $s_start = $snap->start('query');

      if ($snap->strand('query') eq '1' && $c_pstart <= $s_start && $s_start <= $c_pend){
	 push (@keepers, $snap);
      }
      elsif ($snap->strand('query') eq '-1' && $c_mstart <= $s_start && $s_start <= $c_mend){
	 push (@keepers, $snap);
      }
   }

   return \@keepers;
}
#-----------------------------------------------------------------------------
sub get_maker_p_and_t_fastas {
   my $annotations = shift @_;
   
   my $p_fastas = '';
   my $t_fastas = '';
   
   foreach my $an (@$annotations) {
      foreach my $a (@{$an->{t_structs}}) {
	 my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	 $p_fastas .= $$p_fasta;
	 $t_fastas .= $$t_fasta;
      }
   }
   
   return $p_fastas, $t_fastas;
}
#-----------------------------------------------------------------------------
sub merge_chunks{
   my $chunk1 = shift @_;
   my $chunk2 = shift @_;

   if (ref($chunk1) eq 'FastaChunk' && ref($chunk2) ne 'FastaChunk'){
      return $chunk1;
   }
   elsif(ref($chunk2) eq 'FastaChunk' && ref($chunk1) ne 'FastaChunk'){
      die "ERROR: The second argument to main::merge_chunks\n",
          "must be a FastaChunk object";
   }
   elsif ($chunk1->length == 0 ){
      return $chunk2;
   }

   #chunks must be flush to each other and given in order
   die "ERROR: Can not merge chunks at main::mergechunks\n"
      if ($chunk1->offset > $chunk2->offset);
   
   $chunk2->seq($chunk1->seq . $chunk2->seq);
   $chunk2->length($chunk1->length + $chunk2->length);
   $chunk2->offset($chunk1->offset);
   $chunk2->p_cutoff($chunk1->p_cutoff);
   $chunk2->m_cutoff($chunk1->m_cutoff);

   return $chunk2;
}

#-----------------------------------------------------------------------------
sub process_the_chunk_divide{
   my $chunk = shift @_;
   my $split_hit = shift @_;
   my $hit_groups = \@_;

   my $phat_hits;

   foreach my $group (@{$hit_groups}){
      push(@{$phat_hits}, @{$group});
   }

   my ($p_hits, $m_hits) = PhatHit_utils::seperate_by_strand('query', $phat_hits);
   my $p_coors  = PhatHit_utils::to_begin_and_end_coors($p_hits, 'query');
   my $m_coors  = PhatHit_utils::to_begin_and_end_coors($m_hits, 'query');

   foreach my $p_coor (@{$p_coors}){
      $p_coor->[0] -= $chunk->offset();
      $p_coor->[1] -= $chunk->offset();
      $p_coor->[0] = $chunk->length if($p_coor->[0] > $chunk->length);
      $p_coor->[1] = $chunk->length if($p_coor->[1] > $chunk->length);
   }
   foreach my $m_coor (@{$m_coors}){
      $m_coor->[0] -= $chunk->offset();
      $m_coor->[1] -= $chunk->offset();
      $m_coor->[0] = $chunk->length if($m_coor->[0] > $chunk->length);
      $m_coor->[1] = $chunk->length if($m_coor->[1] > $chunk->length);
   }

   my $p_pieces = Shadower::getPieces(\($chunk->seq), $p_coors, 10);
      $p_pieces = [sort {$b->{e} <=> $a->{e}} @{$p_pieces}];
   my $m_pieces = Shadower::getPieces(\($chunk->seq), $m_coors, 10);
      $m_pieces = [sort {$b->{e} <=> $a->{e}} @{$m_pieces}];

   my @keepers;
   my $cutoff = $chunk->length + $chunk->offset - $split_hit;
   my $p_cutoff = $chunk->length + $chunk->offset + 1;
   my $m_cutoff = $chunk->length + $chunk->offset + 1;

   foreach my $p_piece (@{$p_pieces}){
      if ($p_piece->{e} + $chunk->offset >= $cutoff){
	 $p_cutoff = $p_piece->{b} + $chunk->offset;
      }
   }
   foreach my $m_piece (@{$m_pieces}){
      if ($m_piece->{e} + $chunk->offset >= $cutoff){
	 $m_cutoff = $m_piece->{b} + $chunk->offset;
      }
   }

   if ($p_cutoff <= 1 && $m_cutoff <= 1){  #too small, all are heldover for next round
      return $chunk, @keepers;
   }

   foreach my $group (@{$hit_groups}){
      my $group_keepers = [];

      foreach my $hit (@{$group}){
	 my $b = $hit->nB('query');
	 my $e = $hit->nE('query');
	 my $strand = $hit->strand;

	 ($b, $e) = ($e, $b) if $b > $e;
	 
	 if (($e < $p_cutoff && $strand eq '1') ||
	     ($e < $m_cutoff && $strand eq '-1')
	    ){
	       push(@{$group_keepers}, $hit);
	 }
      }

      push(@keepers, $group_keepers);
   }

   my $abs_cutoff = ($p_cutoff < $m_cutoff) ? $p_cutoff -200 : $m_cutoff -200;
   my $sub_strt = $abs_cutoff - $chunk->offset - 1;
      $sub_strt = 0 if($sub_strt < 0);

   my $new_chunk = Storable::dclone($chunk); 
      $new_chunk->p_cutoff($p_cutoff);
      $new_chunk->m_cutoff($m_cutoff);
      $new_chunk->length($chunk->length - ($abs_cutoff - $chunk->offset - 1));
      $new_chunk->offset($abs_cutoff - 1);
      $new_chunk->seq(substr($new_chunk->seq,
			     $sub_strt,
			     $new_chunk->length
			    )
		     );

   return $new_chunk, @keepers;
}
#-----------------------------------------------------------------------------
sub write_quality_data {
   my $quality_indices = shift;
   my $seq_id          = shift;

   my $out_file = $seq_id.'.maker.transcripts.qi';
   my $fh = new FileHandle();
   $fh->open(">$out_file");

   print $fh "genomic_seq\ttranscript\tquality_index\n";

   while (my $d = shift(@{$quality_indices})) {
      my $t_name = $d->[0];
      my $t_qi   = $d->[1];
	
      print $fh "$seq_id\t$t_name\t$t_qi\n";
   }
   $fh->close();
}
#-----------------------------------------------------------------------------
sub get_snap_p_and_t_fastas {
   my $seq   = shift;
   my $snaps = shift;
	
   my $p_fastas = '';
   my $t_fastas = '';
   foreach my $hit (@{$snaps}) {
      my $t_name = $hit->name(); # note this is being set in GFFV3::pred_data
      my $t_seq  = maker::auto_annotator::get_transcript_seq($hit, $seq);	
		
      my ($p_seq, $offset, $end) = 
      maker::auto_annotator::get_translation_seq($t_seq);
		
      my $score = 0;
      foreach my $hsp ($hit->hsps) {
	 $score += $hsp->score();
      }
		
      my $p_def = '>'.$t_name.' protein score:'.$score;
      my $t_def = '>'.$t_name.' snap.transcript offset:'.$offset;
      $t_def.= ' score:'.$score; 
		
      my $p_fasta = Fasta::toFasta($p_def, \$p_seq);
      my $t_fasta = Fasta::toFasta($t_def, \$t_seq);
		
      $p_fastas .= $$p_fasta;
      $t_fastas .= $$t_fasta;
		
   }
   return ($p_fastas, $t_fastas);
}

#-----------------------------------------------------------------------------
sub get_p_and_t_fastas {
   my $t_struct = shift;
	
   my $t_seq  = $t_struct->{t_seq};
   my $p_seq  = $t_struct->{p_seq};
   my $t_off  = $t_struct->{t_offset};
   my $t_name = $t_struct->{t_name};
	
   my $p_def = '>'.$t_name.' protein'; 
   my $t_def = '>'.$t_name.' transcript offset:'.$t_off;
	
   my $p_fasta = Fasta::toFasta($p_def, \$p_seq);
   my $t_fasta = Fasta::toFasta($t_def, \$t_seq);
	
   return($p_fasta, $t_fasta);
}

#----------------------------------------------------------------------------
sub load_anno_hsps {
   my $annotations = shift;
   my @coors;
   my $i = @{$annotations};
   foreach my $an (@$annotations) {
      foreach my $a (@{$an->[0]}) {
	 my $hit = $a->{hit};
	 foreach my $hsp ($hit->hsps()) {
	    push(@coors, [$hsp->nB('query'),
			  $hsp->nE('query'),
			 ]);
	 }
      }
   }
   return (\@coors, $i);;
}

#-----------------------------------------------------------------------------
sub load_clust_hsps {
   my $clusters = shift;
   my @coors;
   my $i = @{$clusters};
   foreach my $c (@$clusters) {
      foreach my $hit (@{$c}) {
	 foreach my $hsp ($hit->hsps()) {
	    push(@coors, [$hsp->nB('query'),
			  $hsp->nE('query'),
			 ]);
	 }
      }
   }
   return (\@coors, $i);
}

#-----------------------------------------------------------------------------
sub load_snap_hsps {
   my $snaps = shift;
   my @coors;
   my $i = @{$snaps};
   foreach my $hit (@{$snaps}) {
      foreach my $hsp ($hit->hsps()) {
	 push(@coors, [$hsp->nB('query'),
		       $hsp->nE('query'),
		      ]);
      }
   }
   return (\@coors, $i);
}

#-----------------------------------------------------------------------------
sub flatten {
   my $clusters = shift;
   my $type     = shift;
   my @hits;
   foreach my $c (@{$clusters}) {
      foreach my $hit (@{$c}) {
	 $hit->type($type) if defined($type);
	 push(@hits, $hit);
      }
   }
   return \@hits;
}

#-----------------------------------------------------------------------------
sub snap {
   my $fasta      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $snap = shift;
   my $snaphmm = shift;
   my $opt_f = shift;

   my %params;
   my $file_name = "$the_void/$seq_id.all";
   my $o_file    = "$the_void/$seq_id\.all\.snap";
   
   FastaFile::writeFile($fasta , $file_name);
		
   runSnap($file_name,
	   $o_file,
	   $snap,
	   $snaphmm,
	   $opt_f
	  );
		
   $params{min_exon_score}  = -100000;	    #-10000;
   $params{min_gene_score}  = -100000;	    #0;
		
   my $chunk_keepers = Widget::snap::parse($o_file,
					   \%params,
					   $file_name,
					  );

   unlink($file_name);

   return $chunk_keepers;
}

#-----------------------------------------------------------------------------
sub runSnap {
   my $q_file   = shift;
   my $out_file = shift;
   my $snap = shift;
   my $snaphmm = shift;
   my $opt_f = shift;

   my $command  = $snap;
   $command .= " $snaphmm";
   $command .= " $q_file";
   $command .= " > $out_file";
	
   my $w = new Widget::snap();
	
   if (-e $out_file && ! $opt_f) {
      print STDERR "re reading snap report.\n";
      print STDERR "$out_file\n";
   }
   else {
      print STDERR "running  snap.\n";
      $w->run($command);
   }
}

#-----------------------------------------------------------------------------
sub polish_exonerate {
   my $g_fasta           = shift;
   my $phat_hit_clusters = shift;
   my $db_index          = shift;
   my $the_void          = shift;
   my $depth             = shift;
   my $type              = shift;
   my $exonerate         = shift;
   my $percov            = shift;
   my $percid            = shift;
   my $score_limit       = shift;
   my $matrix            = shift;
   my $opt_f             = shift;

   my $def = Fasta::getDef($g_fasta);
   my $seq = Fasta::getSeq($g_fasta);
	
   my $exe = $exonerate;
	
	
   my @exonerate_clusters;
   my $i = 0;
   foreach my $c (@{$phat_hit_clusters}) {
      my $n = 0;
      my $got_some = 0;

      foreach my $hit (@{$c}) {
	 last if $n == $depth;

	 next if $hit->pAh < $percov;
	 next if $hit->hsp('best')->frac_identical < $percid;
	 
	 my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit,'query');

	 my @coors = [$nB, $nE];
	 my $p = Shadower::getPieces($seq, \@coors, 50);
	 my $p_def = $def." ".$p->[0]->{b}." ".$p->[0]->{e};
	 my $p_fasta = Fasta::toFasta($p_def, \$p->[0]->{piece});
	 my ($name) = $p_def =~ />([^\s\t\n]+)/;

	 #build a safe name for file names from the sequence identifier
	 my $safe_name = uri_escape($name,
				    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				   );
	 $safe_name .= '.fasta';
	 my $d_file = $the_void."/".$safe_name.'.'.$i.'.'.$n;
	 FastaFile::writeFile($p_fasta, $d_file);
	 my $offset = $p->[0]->{b} - 1;
	 my $id  = $hit->name();
	 $id =~ s/\s+/_/g;
	 $id =~ s/\|/_/g;
	 my $fastaObj = $db_index->get_Seq_by_id($hit->name);
	 if (not $fastaObj) {
	    print "stop here:".$hit->name."\n";
	    die;
	 }
	 my $seq      = $fastaObj->seq();
	 my $def      = $db_index->header($hit->name);
	 $def =~ s/\|/_/g;
	 my $fasta    = Fasta::toFasta('>'.$def, \$seq);

	 #build a safe name for file names from the sequence identifier
	 my $safe_id = uri_escape($id, 
				  '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				 );

	 my $t_file    = $the_void."/".$safe_id.'.'.$i.'.'.$n.'.fasta';
	 my $ext = "$i\.$n";
	 FastaFile::writeFile($fasta, $t_file);

	 my $exonerate_hits = to_polisher($d_file,
					  $t_file,
					  $the_void,
					  $offset,
					  $type,
					  $ext,
					  $exe,
					  $score_limit,
					  $matrix,
					  $opt_f
					 );

	 foreach my $exonerate_hit (@{$exonerate_hits}) {
	    if (defined($exonerate_hit) && exonerate_okay($exonerate_hit)) {
	       $n++;
	       push(@{$exonerate_clusters[$i]}, $exonerate_hit);
	       $got_some = 1;
	    }
	 }
      }
      $i++ if $got_some;
   }
   return \@exonerate_clusters;
}

#-----------------------------------------------------------------------------
sub exonerate_okay {
   my $hit  = shift;

   my $i = 0;
   foreach my $hsp ($hit->hsps()) {
      return 0 unless defined($hsp->nB('query'));
      return 0 unless defined($hsp->nE('query'));
      return 0 unless defined($hsp->nB('hit'));
      return 0 unless defined($hsp->nE('hit'));
      return 0 unless defined($hsp->strand('query'));
      return 0 unless defined($hsp->strand('query'));
      return 0 unless defined($hsp->strand('hit'));
      return 0 unless defined($hsp->strand('hit'));

      my $q_str = $hsp->query_string();
      my $h_str = $hsp->hit_string();
		
      if ($h_str =~ /Target Intron/) {
	 print STDERR "BADDD EXONERATE!\n";
	 sleep 4;
	 return 0;
      } elsif ($q_str =~ /Target Intron/) {
	 print STDERR "BADDD EXONERATE!\n";
	 sleep 4;
	 return 0;
      }
      $i++;
   }

   return 1 
}

#-----------------------------------------------------------------------------
sub to_polisher {
   my $d_file   = shift;
   my $t_file   = shift;
   my $the_void = shift;
   my $offset   = shift;
   my $type     = shift;
   my $ext      = shift;
   my $exe      = shift;
   my $score_limit = shift;
   my $matrix = shift;
   my $opt_f = shift;

   if ($type eq 'p') {
      return polisher::exonerate::protein::polish($d_file,
						  $t_file,
						  $the_void,
						  $offset,
						  $ext,
						  $exe,
						  $score_limit,
						  $matrix,
						  $opt_f
						 );
   } elsif ($type eq 'e') {
      return polisher::exonerate::est::polish($d_file,
					      $t_file,
					      $the_void,
					      $offset,
					      $ext,
					      $exe,
					      $score_limit,
					      $matrix,
					      $opt_f
					     );
   } else {
      die "unknown type:$type in sub to_polisher.\n";
   }
}

#-----------------------------------------------------------------------------
sub make_multi_fasta {
   my $index    = shift;
   my $clusters = shift;;
   my $fastas = '';
   foreach my $c (@{$clusters}) {
      foreach my $hit (@{$c}) {
	 my $id = $hit->name();
	 my $fastaObj = $index->get_Seq_by_id($id);
	 my $seq      = $fastaObj->seq(); 
	 my $def      = $index->header($id);
	 my $fasta    = Fasta::toFasta('>'.$def, \$seq);
	 $fastas     .= $$fasta; 
      }
   }
   return \$fastas;
}

#-----------------------------------------------------------------------------
sub build_fasta_index {
   my $db = shift;
   my $index = new Bio::DB::Fasta($db);
   return $index;
}

#-----------------------------------------------------------------------------
sub repeatmask {
   my $chunk        = shift;
   my $the_void     = shift;
   my $seq_id       = shift;
   my $model_org    = shift;
   my $RepeatMasker = shift;
   my $rmlib        = shift;
   my $cpus         = shift;
   my $opt_f        = shift;

   my $chunk_number = $chunk->number();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.out";
   my $q_length = $chunk->parent_seq_length();
   my $query_def = $chunk->parent_def();
   my $query_seq = $chunk->seq();

   $chunk->write_file($file_name);
		
   runRepeatMasker($file_name, 
		   $model_org, 
		   $the_void, 
		   $o_file,
		   $RepeatMasker,
		   $rmlib,
		   $cpus,
		   $opt_f
		  );		# -no_low
		
   my $rm_chunk_keepers = Widget::RepeatMasker::parse($o_file, 
						      $seq_id, 
						      $q_length,
						      $opt_f
						     );
  
   PhatHit_utils::add_offset($rm_chunk_keepers, 
			     $chunk->offset(),
			    );
   #     PhatHit_utils::merge_hits($rm_keepers,  
   # 			      $rm_chunk_keepers, 
   # 			      20,
   # 			     );

   $chunk->erase_fasta_file();
	
   return ($rm_chunk_keepers);
}

#-----------------------------------------------------------------------------
sub blastn {
   my $chunk      = shift;
   my $db         = shift;
   my $the_void    = shift;
   my $seq_id     = shift;
   my $blastn = shift;
   my $eval_blastn = shift;
   my $bit_blastn = shift,
   my $percov_blastn = shift;
   my $percid_blastn = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $old_db = shift;
   my $xdformat = shift;
   my $id = shift;
   my $rank = shift;
   my $opt_f = shift;

   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;
	
   my $chunk_number = $chunk->number();
   my $q_length = $chunk->parent_seq_length();
 
   my ($db_old_n) = $old_db =~ /([^\/]+)$/;
   $db_old_n  =~ s/\.fasta$//;
   my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.blastn";

   my $t_dir = "/tmp/rank".$rank;
   mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $o_file    = "$blast_finished\.temp_dir/$db_n\.blastn";

   $db =~ /([^\/]+$)/;
   my $tmp_db = "$t_dir/$1";

   if (-e $blast_finished && ! $opt_f) {
      $o_file = $blast_finished;
	    
	 return [] if ($id !~ /\:0$/);
      }
      elsif (! -e $blast_finished && ! @{[<$tmp_db.n??*>]}) {
	 system("cp $db $tmp_db");
	 system ("$xdformat -n $tmp_db");
      }
	
   $chunk->write_file($t_file_name);
   runBlastn($t_file_name,
	     $tmp_db,
	     $o_file,
	     $blastn,
	     $eval_blastn,
	     $cpus,
	     $opt_f
	    );

   my %params;
   $params{significance}  = $eval_blastn;
   $params{hsp_bit_min}   = $bit_blastn;
   $params{percov}        = $percov_blastn;
   $params{percid}        = $percid_blastn;
   $params{split_hit}     = $split_hit;

   my $chunk_keepers = Widget::blastn::parse($o_file,
					     \%params,
					    );
   
   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset(),
			    );

   $chunk->erase_fasta_file();
    
   if($chunk->p_cutoff || $chunk->m_cutoff){
      my @keepers;
      
      foreach my $hit (@{$chunk_keepers}){
	 if ($hit->strand('query') eq '1' && $hit->start('query') >= $chunk->p_cutoff){
	    push (@keepers, $hit)
	 }
	 elsif ($hit->strand('query') eq '-1' && $hit->start('query') >= $chunk->m_cutoff){
	    push (@keepers, $hit)
	 }
      }
      
      return \@keepers;
   }
   else{
      return $chunk_keepers
   }
}

#-----------------------------------------------------------------------------
sub runBlastn {
   my $q_file   = shift;
   my $db       = shift;
   my $out_file = shift;
   my $blastn = shift;
   my $eval_blastn = shift;
   my $cpus = shift;
   my $opt_f = shift;

   my $command  = $blastn;
   $command .= " $db $q_file B=10000 V=10000 E=$eval_blastn";
   $command .= " wordmask=seg";
   $command .= " R=3";
   $command .= " W=15";
   $command .= " M=1";
   $command .= " N=-3";
   $command .= " Q=3";
   $command .= " Z=128000000";
   $command .= " cpus=$cpus";	
   $command .= " topcomboN=1";
   $command .= " hspmax=100";
   $command .= " gspmax=100";
   $command .= " hspsepqmax=10000";
   $command .= " lcmask";
   $command .= " filter=seg";
   $command .= " gi";
   #$command .= " mformat=2"; # remove for full report
   $command .= " > $out_file";
	
   my $w = new Widget::blastn();
   if (-e $out_file && ! $opt_f) {
      print STDERR "re reading blast report.\n";
      print STDERR "$out_file\n";
   }
   else {
      print STDERR "running  blast search.\n";
      my $dir = $out_file;
      $dir =~ s/[^\/]+$//;
      mkpath ($dir);
      $w->run($command);
   }
}

#-----------------------------------------------------------------------------
sub blastx {
   my $chunk         = shift;
   my $db            = shift;
   my $the_void      = shift;
   my $seq_id        = shift;
   my $blastx        = shift;
   my $eval_blastx   = shift;
   my $bit_blastx    = shift;
   my $percov_blastx = shift;
   my $percid_blastx = shift;
   my $split_hit     = shift;
   my $cpus          = shift;
   my $old_db        = shift;
   my $xdformat      = shift;
   my $id            = shift;
   my $rank          = shift;
   my $opt_f         = shift;
	
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;

   my $q_length = $chunk->parent_seq_length();
   my $chunk_number = $chunk->number();
    
   my ($db_old_n) = $old_db =~ /([^\/]+)$/;
   $db_old_n  =~ s/\.fasta$//;
   my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.blastx";
    
   my $t_dir = "/tmp/rank".$rank;
   mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $o_file    = "$blast_finished\.temp_dir/$db_n\.blastx";
    
   $db =~ /([^\/]+$)/;
   my $tmp_db = "$t_dir/$1";
	     
   if (-e $blast_finished && ! $opt_f) {
      $o_file = $blast_finished ;
	    
      return [] if ($id !~ /\:0$/);
   }
   elsif (! @{[<$tmp_db.n??*>]}) {
      system("cp $db $tmp_db");
      system ("$xdformat -p $tmp_db");
   }
	     
   $chunk->write_file($t_file_name);
   runBlastx($t_file_name,
	     $tmp_db,
	     $o_file,
	     $blastx,
	     $eval_blastx,
	     $cpus,
	     $opt_f
	    );

   my %params;
   $params{significance}  = $eval_blastx;
   $params{hsp_bit_min}   = $bit_blastx;
   $params{percov}        = $percov_blastx;
   $params{percid}        = $percid_blastx;
   $params{split_hit}     = $split_hit;
   
   my $chunk_keepers = Widget::blastx::parse($o_file,
					     \%params,
					    );

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset(),
			    );

   $chunk->erase_fasta_file();

   if($chunk->p_cutoff || $chunk->m_cutoff){
      my @keepers;
      
      foreach my $hit (@{$chunk_keepers}){
	 if ($hit->strand('query') eq '1' && $hit->start('query') >= $chunk->p_cutoff){
	    push (@keepers, $hit)
	 }
	 elsif ($hit->strand('query') eq '-1' && $hit->start('query') >= $chunk->m_cutoff){
	    push (@keepers, $hit)
	 }
      }
      
      return \@keepers;
   }
   else{
      return $chunk_keepers
   }
}

#-----------------------------------------------------------------------------
sub runBlastx {
   my $q_file   = shift;
   my $db       = shift;
   my $out_file = shift;
   my $blastx = shift;
   my $eval_blastx = shift;
   my $cpus = shift;
   my $opt_f = shift;

   my $command  = $blastx;
   $command .= " $db $q_file B=10000 V=10000 E=$eval_blastx";
   $command .= " wordmask=seg";
   #$command .= " T=20";
   #$command .= " W=5";
   #$command .= " wink=5";
   $command .= " Z=300";
   $command .= " Y=500000000";
   $command .= " hspmax=100";
   $command .= " cpus=$cpus";
   $command .= " gspmax=100";
   $command .= " hspsepqmax=10000";
   $command .= " lcfilter";
   $command .= " filter=seg";
   $command .= " gi";
   #$command .= " mformat=2"; # remove for full report
   $command .= " > $out_file";
   my $w = new Widget::blastx();

   if (-e $out_file  && ! $opt_f) {
      print STDERR "re reading blast report.\n";
      print STDERR "$out_file\n";
   }
   else {
      print STDERR "running  blast search.\n";
      my $dir = $out_file;
      $dir =~ s/[^\/]+$//;
      mkpath ($dir);
      $w->run($command);
   }
}

#-----------------------------------------------------------------------------
sub runRepeatMasker {
   my $q_file   = shift;
   my $species  = shift;
   my $dir      = shift;
   my $o_file   = shift;
   my $RepeatMasker = shift;
   my $rmlib = shift;
   my $cpus = shift;
   my $opt_f = shift;
   my $no_low   = shift;
	
   my $command  = $RepeatMasker;
    
   if ($rmlib) {
      $command .= " $q_file -lib $rmlib -dir $dir ";    
   } else {
      $command .= " $q_file -species $species -dir $dir ";
   }
   $command .= " -nolow" if defined($no_low);
	
   my $w = new Widget::RepeatMasker();
   if (-e $o_file && ! $opt_f) {
      print STDERR "re reading repeat masker report.\n";
      print STDERR "$o_file\n";
   }
   else {
      print STDERR "running  repeat masker.\n";
      $w->run($command);
   }
}

#-----------------------------------------------------------------------------
sub build_the_void {
   my $seq_id  = shift;
   my $out_dir = shift;

   my $vid = "theVoid\.$seq_id";   
   my $the_void = "$out_dir"."$vid";
   mkpath ($the_void);

   return $the_void;
}

#-----------------------------------------------------------------------------
1;
