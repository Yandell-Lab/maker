#------------------------------------------------------------------------
#----                         Shared_Functios                        ---- 
#------------------------------------------------------------------------
package Shared_Functions;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
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
use Widget::augustus;
use Widget::xdformat;
use PhatHit_utils;
use Shadower;
use Bio::DB::Fasta;
use polisher::exonerate::protein;
use polisher::exonerate::est;
use maker::auto_annotator;
use cluster;
use repeat_mask_seq;
use maker::sens_spec;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- CLASS FUNCTIONS ----------------------------
#------------------------------------------------------------------------
sub get_preds_on_chunk {
   my $preds = shift;
   my $chunk = shift;

   my $c_start = $chunk->offset + 1;
   my $c_end = $chunk->offset + $chunk->length;
   
   my @keepers;
   foreach my $pred (@{$preds}) {
      my $s_start = $pred->start('query');

      if ($pred->strand('query') eq '1' && $c_start <= $s_start && $s_start <= $c_end) {
	 push (@keepers, $pred);
      }
   }

   return \@keepers;
}

#-----------------------------------------------------------------------------
sub merge_and_resolve_hits{
   my $fasta = shift @_;
   my $fasta_t_index = shift @_;
   my $fasta_p_index = shift @_;
   my $blastn_keepers = shift @_;
   my $blastx_keepers = shift @_;
   my $blastn_holdovers = shift @_;
   my $blastx_holdovers = shift @_;
   my $the_void = shift @_;
   my %CTL_OPTIONS = %{shift @_};
   my $OPT_f = shift @_;
   my $LOG = shift @_;

   PhatHit_utils::merge_hits($blastn_keepers,  
			     $blastn_holdovers, 
			     $CTL_OPTIONS{split_hit},
			    );
   @{$blastn_holdovers} = ();

   PhatHit_utils::merge_hits($blastx_keepers,  
			     $blastx_holdovers, 
			     $CTL_OPTIONS{split_hit},
			    );
   @{$blastx_holdovers} = ();

   $blastn_keepers = reblast_merged_hits($fasta,
					 $blastn_keepers,
					 $fasta_t_index,
					 $the_void,
					 'blastn',
					 \%CTL_OPTIONS,
					 $OPT_f,
					 $LOG
					);

   $blastx_keepers = reblast_merged_hits($fasta,
					 $blastx_keepers,
					 $fasta_p_index,
					 $the_void,
					 'blastx',
					 \%CTL_OPTIONS,
					 $OPT_f,
					 $LOG
					);

   return ($blastn_keepers, $blastx_keepers);
}
#-----------------------------------------------------------------------------
sub reblast_merged_hits {
   my $g_fasta     = shift @_;
   my $hits         = shift @_;
   my $db_index    = shift @_;
   my $the_void    = shift @_;
   my $type        = shift @_;
   my %CTL_OPTIONS = %{shift @_};
   my $OPT_f       = shift @_;
   my $LOG         = shift @_;

   #==get data from parent fasta

   #parent fasta get def and seq
   my $par_def = Fasta::getDef($$g_fasta);
   my $par_seq = Fasta::getSeq($$g_fasta);

   #get seq id off def line
   my ($p_id)  = $par_def =~ /^([^\s\t\n]+)/;
   $p_id =~ s/\|/_/g;

   #build a safe name for file names from the sequence identifier
   my $p_safe_id = uri_escape($p_id, 
			      '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			     );

   my @blast_keepers;

   #==check whether to re-blast each hit
   foreach my $hit (@{$hits}) {
      #if not a merged hit take as is
      if (not $hit->{'_sequences_was_merged'}) {
	 push (@blast_keepers, $hit);
	 next;
      }

      #== excise region of query fasta and build new chunk object for blast
      
      #find region of hit to get piece
      my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit,'query');   
      my @coors = [$nB, $nE];
      my $piece = Shadower::getPieces($par_seq, \@coors, $CTL_OPTIONS{split_hit});
      
      #get piece fasta def, seq, and offset
      my $piece_def = $par_def." ".$piece->[0]->{b}." ".$piece->[0]->{e};
      my $piece_seq = $piece->[0]->{piece};
      my $offset = $piece->[0]->{b} - 1;
      
      #make the piece into a fasta chunk
      my $chunk = new FastaChunk();
      $chunk->seq($piece_seq);
      $chunk->def($piece_def);
      $chunk->parent_def($par_def);
      $chunk->size(length($$par_seq));	     #the max size of a chunk
      $chunk->length(length($chunk->seq())); #the actual size of a chunk                                                
      $chunk->offset($offset);
      $chunk->number(0);
      $chunk->is_last(1);
      $chunk->parent_seq_length(length($$par_seq));

      #==build new fasta and db for blast search from hit name and db index
      
      #get name to search db index
      my $t_id  = $hit->name();
      $t_id =~ s/\s+/_/g;
      $t_id =~ s/\|/_/g;
      
      #build a safe name for file names from the sequence identifier
      my $t_safe_id = uri_escape($t_id, 
				 '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				);
      #search db index
      my $fastaObj = $db_index->get_Seq_by_id($hit->name);
      if (not $fastaObj) {
	 print STDERR "stop here:".$hit->name."\n";
	 die "ERROR: Fasta index error\n";
      }
      
      #get fasta def and seq
      my $t_seq      = $fastaObj->seq();
      my $t_def      = $db_index->header($hit->name);
      $t_def =~ s/\|/_/g;
      
      #write fasta file
      my $fasta = Fasta::toFasta('>'.$t_def, \$t_seq);
      my $t_file = $the_void."/".$t_safe_id.'.for_'.$type.'.fasta';
      FastaFile::writeFile($fasta, $t_file);
      
      #build db for blast using xdformat
      xdformat($CTL_OPTIONS{xdformat}, $t_file, $type, $CTL_OPTIONS{alt_peptide});
      
      #==run the blast search
      if ($type eq 'blastx') {
	  
	  print STDERR "re-running blast against $t_id...\n" unless $main::quiet;
	  my $keepers = blastx($chunk, 
			       $t_file,
			       $the_void,
			       $p_safe_id."-2-".$t_safe_id,
			       $CTL_OPTIONS{blastx},
			       $CTL_OPTIONS{eval_blastx},
			       $CTL_OPTIONS{bit_blastx},
			       $CTL_OPTIONS{percov_blastx},
			       $CTL_OPTIONS{percid_blastx},
			       $CTL_OPTIONS{split_hit},
			       $CTL_OPTIONS{cpus},
			       $OPT_f,
			       $LOG
			      );
	  
	  push(@blast_keepers, @{$keepers});
	  print STDERR "...finished\n" unless $main::quiet;
      }
      elsif ($type eq 'blastn') {
	  
	  print STDERR "re-running blast against $t_id...\n" unless $main::quiet;
	  my $keepers = blastn($chunk, 
			       $t_file,
			       $the_void,
			       $p_safe_id."-2-".$t_safe_id,,
			       $CTL_OPTIONS{blastn},
			       $CTL_OPTIONS{eval_blastn},
			       $CTL_OPTIONS{bit_blastn},
			       $CTL_OPTIONS{percov_blastn},
			       $CTL_OPTIONS{percid_blastn},
			       $CTL_OPTIONS{split_hit},
			       $CTL_OPTIONS{cpus},
			       $OPT_f,
			       $LOG
			      );
	  
	  push(@blast_keepers, @{$keepers});
	  print STDERR "...finished\n" unless $main::quiet;
      }
      else {
	  die "ERROR: Invaliv type \'$type\' in maker::reblast_merged_hit\n";
      }

      unlink <$t_file.x??>;
  }
   
   #==return hits
   return (\@blast_keepers);
}
#-----------------------------------------------------------------------------
sub process_the_chunk_divide{
   my $chunk = shift @_;
   my $split_hit = shift @_;
   my $hit_groups = \@_; #processed and returned in order given by user

   my $phat_hits;

   foreach my $group (@{$hit_groups}) {
      push(@{$phat_hits}, @{$group});
   }

   my ($p_hits, $m_hits) = PhatHit_utils::seperate_by_strand('query', $phat_hits);
   my $p_coors  = PhatHit_utils::to_begin_and_end_coors($p_hits, 'query');
   my $m_coors  = PhatHit_utils::to_begin_and_end_coors($m_hits, 'query');

   foreach my $p_coor (@{$p_coors}) {
      $p_coor->[0] -= $chunk->offset();
      $p_coor->[1] -= $chunk->offset();
      $p_coor->[0] = $chunk->length if($p_coor->[0] > $chunk->length);
      $p_coor->[1] = $chunk->length if($p_coor->[1] > $chunk->length);
   }
   foreach my $m_coor (@{$m_coors}) {
      $m_coor->[0] -= $chunk->offset();
      $m_coor->[1] -= $chunk->offset();
      #fix coordinates for hits outside of chunk end
      $m_coor->[0] = $chunk->length if($m_coor->[0] > $chunk->length);
      $m_coor->[1] = $chunk->length if($m_coor->[1] > $chunk->length);
      #fix coordinates for hits outside of chunk begin
      $m_coor->[0] = 0 if($m_coor->[0] < 0);
      $m_coor->[1] = 0 if($m_coor->[1] < 0);
   }

   my $p_pieces = Shadower::getPieces(\($chunk->seq), $p_coors, 10);
   $p_pieces = [sort {$b->{e} <=> $a->{e}} @{$p_pieces}];
   my $m_pieces = Shadower::getPieces(\($chunk->seq), $m_coors, 10);
   $m_pieces = [sort {$b->{e} <=> $a->{e}} @{$m_pieces}];

   my @keepers;
   my @holdovers;

   my $cutoff = $chunk->length + $chunk->offset - $split_hit;
   my $p_cutoff = $chunk->length + $chunk->offset + 1;
   my $m_cutoff = $chunk->length + $chunk->offset + 1;

   foreach my $p_piece (@{$p_pieces}) {
      if ($p_piece->{e} + $chunk->offset >= $cutoff) {
	 $p_cutoff = $p_piece->{b} + $chunk->offset;
      }
   }
   foreach my $m_piece (@{$m_pieces}) {
      if ($m_piece->{e} + $chunk->offset >= $cutoff) {
	 $m_cutoff = $m_piece->{b} + $chunk->offset;
      }
   }

   if ($p_cutoff <= 1 && $m_cutoff <= 1) { #too small, all are heldover for next round
      push (@holdovers, @{$hit_groups});
      return @holdovers, @keepers;
   }

   foreach my $group (@{$hit_groups}) {
      my $group_keepers = [];
      my $group_holdovers = [];

      foreach my $hit (@{$group}) {
	 my $b = $hit->nB('query');
	 my $e = $hit->nE('query');
	 my $strand = $hit->strand;

	 ($b, $e) = ($e, $b) if $b > $e;
	 
	 if (($e < $p_cutoff && $strand eq '1') ||
	     ($e < $m_cutoff && $strand eq '-1')
	    ) {
	    push(@{$group_keepers}, $hit);
	 }
	 else {
	    push(@{$group_holdovers}, $hit);
	 }
      }

      push(@keepers, $group_keepers);
      push(@holdovers, $group_holdovers);
   }

   my $abs_cutoff = ($p_cutoff < $m_cutoff) ? $p_cutoff -200 : $m_cutoff -200;
   my $sub_strt = $abs_cutoff - $chunk->offset - 1;
   $sub_strt = 0 if($sub_strt < 0);

   #hit holdovers and keepers are returned in same order given by user
   return @holdovers, @keepers; 
}
#-----------------------------------------------------------------------------
sub build_datastore {
   my %CTL_OPTIONS = %{shift @_};

   $CTL_OPTIONS{'dsroot'} = "$CTL_OPTIONS{'out_base'}/$CTL_OPTIONS{'out_name'}_datastore";
   $CTL_OPTIONS{'dsindex'} = "$CTL_OPTIONS{'out_base'}/$CTL_OPTIONS{'out_name'}_master_datastore.index";

   print STDERR "A data structure will be created for you at:\n".
                "$CTL_OPTIONS{'dsroot'}\n\n".
                "To access files for individual sequences use the datastore index:\n".
                "$CTL_OPTIONS{'dsindex'}\n\n";
    
   $CTL_OPTIONS{'datastore'} = new Datastore::MD5('root' => $CTL_OPTIONS{'dsroot'}, 'depth' => 2);

   return %CTL_OPTIONS;
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
sub create_blastdb {
   my %CTL_OPTIONS = %{shift @_};
   my %OPT =  %{shift @_};

   xdformat($CTL_OPTIONS{xdformat}, $CTL_OPTIONS{'protein'}, 'blastx', $CTL_OPTIONS{'alt_peptide'});
   xdformat($CTL_OPTIONS{xdformat}, $CTL_OPTIONS{'est'}, 'blastn');
   xdformat($CTL_OPTIONS{xdformat}, $CTL_OPTIONS{'repeat_protein'}, 'blastx', $CTL_OPTIONS{'alt_peptide'})
       if (! $OPT{R} && ! $OPT{GFF});
}
#----------------------------------------------------------------------------
sub create_mpi_blastdb {
   my $CTL_OPTIONS = shift @_;
   my $OPT = shift @_;
   my $mpi_size = shift @_;

   $CTL_OPTIONS->{'old_protein'}        = $CTL_OPTIONS->{'protein'};
   $CTL_OPTIONS->{'old_est'}            = $CTL_OPTIONS->{'est'};
   $CTL_OPTIONS->{'old_repeat_protein'} = $CTL_OPTIONS->{'repeat_protein'};
    
   $CTL_OPTIONS->{'protein'} = split_db($CTL_OPTIONS->{'protein'}, $mpi_size);
   $CTL_OPTIONS->{'est'} = split_db($CTL_OPTIONS->{'est'}, $mpi_size);
   
   $CTL_OPTIONS->{'repeat_protein'} = split_db($CTL_OPTIONS->{'repeat_protein'}, $mpi_size)
       unless($OPT->{R});
}
#----------------------------------------------------------------------------
sub split_db {
   my $file = shift @_;
   my $mpi_size = shift @_;

   my $fasta_iterator = new Iterator::Fasta($file);
   my $db_size = $fasta_iterator->number_of_entries();
   my $bins = $mpi_size;
   $bins = $db_size if ($db_size < $bins);

   my @fhs;
   my @db_files;

   if ($bins == 1) {
      push (@db_files, $file);
	
      return \@db_files;
   }

   my ($f_name) = $file =~ /([^\/]+$)/;
   $f_name =~ s/\.fasta$//;
    
   my $d_name = "$f_name\.mpi\.$mpi_size";
   my $b_dir = cwd(). "/mpi_blastdb";
   my $f_dir = "$b_dir/$d_name";
   my $t_dir = "/tmp/$d_name";

   if(-e "$f_dir"){
      @db_files = File::Find::Rule->file->name("*$d_name\.*")->in( "$f_dir");
      return \@db_files;
   }

   mkdir($t_dir);
   mkdir($b_dir) unless (-e $b_dir);

   for (my $i = 0; $i < $bins; $i++) {
      my $name = "$t_dir/$d_name\.$i\.fasta";
      my $fh;
      open ($fh, "> $name");

      push (@fhs, $fh);
   }

   while (my $fasta = $fasta_iterator->nextEntry()) {
      my $fh = shift @fhs;
      print $fh "$fasta\n";
      push (@fhs, $fh);
   }

   foreach my $fh (@fhs) {
      close ($fh);
   }

   print `mv $t_dir $f_dir`;

   if(-e "$f_dir"){
       @db_files = File::Find::Rule->file->name("*$d_name\.*")->in( "$f_dir");
       return \@db_files;
   }
   else{
       die "ERROR: Could not split db\n";
   }
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
   my $snap       = shift;
   my $snaphmm    = shift;
   my $opt_f = shift;
   my $LOG = shift;
	
   my %params;
   my $file_name = "$the_void/$seq_id.all";
   my $o_file    = "$the_void/$seq_id\.all\.snap";


   $LOG->add_entry("STARTED", $o_file, "");   

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

   $LOG->add_entry("FINISHED", $o_file, "");

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
      print STDERR "re reading snap report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  snap.\n" unless $main::quiet;
      $w->run($command);
   }
}
#-----------------------------------------------------------------------------
sub augustus {
   my $fasta      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $exe        = shift;
   my $org        = shift;
   my $opt_f = shift;
   my $LOG = shift;

   my %params;
   my $file_name = "$the_void/$seq_id.all";
   my $o_file    = "$the_void/$seq_id\.all\.augustus";

   $LOG->add_entry("STARTED", $o_file, ""); 

   FastaFile::writeFile($fasta , $file_name) unless -e $file_name;

   runAugustus($file_name,
               $o_file,
               $exe,
               $org,
	       $opt_f
	      );

   $params{min_exon_score}  = -100000;      #-10000;
   $params{min_gene_score}  = -100000;      #0;

   my $chunk_keepers = Widget::augustus::parse($o_file,
					       \%params,
					       $file_name
					      );

   $LOG->add_entry("FINISHED", $o_file, "");

   unlink($file_name);

   return $chunk_keepers;
}

#-----------------------------------------------------------------------------
sub runAugustus {
   my $q_file   = shift;
   my $out_file = shift;
   my $exe      = shift;
   my $org      = shift;
   my $opt_f = shift;

   my $command  = $exe;
   $command .= ' --species='."$org";
   $command .= " $q_file";
   $command .= " > $out_file";

   my $w = new Widget::augustus();

   if (-e $out_file && ! $opt_f) {
      print STDERR "re reading augustus report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  augustus.\n" unless $main::quiet;
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
   my $LOG               = shift;

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
	    die "ERROR: Fasta index error\n";
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
					  $opt_f,
					  $LOG
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
   my $LOG   = shift;

   if ($type eq 'p') {
      return polisher::exonerate::protein::polish($d_file,
						  $t_file,
						  $the_void,
						  $offset,
						  $ext,
						  $exe,
						  $score_limit,
						  $matrix,
						  $opt_f,
						  $LOG
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
					      $opt_f,
					      $LOG
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
sub build_all_indexes {
   my @dbs = @_;

   foreach my $db (@dbs) {
      my $index = new Bio::DB::Fasta($db);
   }
}
#-----------------------------------------------------------------------------
sub xdformat {
    my $xdformat = shift;
    my $file = shift;
    my $type = shift;
    my $alt_peptide = shift || 'c';

    die "ERROR: Can not find xdformat executable\n" if(! -e $xdformat);
    die "ERROR: Can not find the db file $file\n" if(! -e $file);
    die "ERROR: You must define a type (blastn|blastx)\n" if(! $type);

    my $command = $xdformat;
    if ( $type eq 'blastx') {
	$command .= " -p -C $alt_peptide $file";
    }
    elsif ($type eq 'blastn'){
	$command .= " -n $file";
    }
    else{
	die "ERROR: Unsuported type $type\n";
    }

    if (($type eq 'blastn' && ! -e $file.'.xnd') ||
	($type eq 'blastx' && ! -e $file.'.xpd')
       ) {
       my $w = new Widget::xdformat();
       print STDERR "running  xdformat.\n" unless $main::quiet;
       $w->run($command);
    }
}
#-----------------------------------------------------------------------------
sub blastn_as_chunks {
   my $chunk      = shift;
   my $db         = shift;
   my $the_void    = shift;
   my $seq_id     = shift;
   my $blastn = shift;
   my $eval_blastn = shift;
   my $cpus = shift;
   my $old_db = shift;
   my $xdformat = shift;
   my $rank = shift;
   my $opt_f = shift;
   my $LOG = shift;
   my $LOG_FLAG = shift;

   #build names for files to use and copy
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;
	
   my $chunk_number = $chunk->number();
 
   my ($db_old_n) = $old_db =~ /([^\/]+)$/;
   $db_old_n  =~ s/\.fasta$//;
   my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.blastn";

   my $t_dir = "/tmp/rank".$rank;
   File::Path::mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.blastn";

   $db =~ /([^\/]+)$/;
   my $tmp_db = "$t_dir/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #copy db to local tmp dir and run xdformat 
   if (! @{[<$tmp_db.xn?*>]} && (! -e $blast_finished || $opt_f) ) {
      system("cp $db $tmp_db");
      xdformat($xdformat, $tmp_db, 'blastn');
   }
   elsif (-e $blast_finished && ! $opt_f) {
      print STDERR "re reading blast report.\n" unless $main::quiet;
      print STDERR "$blast_finished\n" unless $main::quiet;
      return $blast_dir;
   }
	
   #call blast executable
   $chunk->write_file($t_file_name);  

   runBlastn($t_file_name,
	     $tmp_db,
	     $o_file,
	     $blastn,
	     $eval_blastn,
	     $cpus,
	     $opt_f
	    );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_blastn{
   my $chunk      = shift;
   my $blast_dir    = shift;
   my $eval_blastn = shift;
   my $bit_blastn = shift,
   my $percov_blastn = shift;
   my $percid_blastn = shift;
   my $split_hit = shift;
   my $opt_f = shift;
   my $LOG = shift;

   my $blast_finished = $blast_dir;
   $blast_finished =~ s/\.temp_dir$//;

   #merge blast reports
   if (! -e $blast_finished || $opt_f) {
      system ("cat $blast_dir/*blast* > $blast_finished");
      File::Path::rmtree ("$blast_dir");
   }

   my %params;
   $params{significance}  = $eval_blastn;
   $params{hsp_bit_min}   = $bit_blastn;
   $params{percov}        = $percov_blastn;
   $params{percid}        = $percid_blastn;
   $params{split_hit}     = $split_hit;

   my $chunk_keepers = Widget::blastn::parse($blast_finished,
					     \%params,
					    );

   $LOG->add_entry("FINISHED", $blast_finished, "");
   
   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset(),
			    );

   if ($chunk->p_cutoff || $chunk->m_cutoff) {
      my @keepers;
      
      foreach my $hit (@{$chunk_keepers}) {
	 if ($hit->strand('query') eq '1' && $hit->start('query') >= $chunk->p_cutoff) {
	    push (@keepers, $hit)
	 }
	 elsif ($hit->strand('query') eq '-1' && $hit->start('query') >= $chunk->m_cutoff) {
	    push (@keepers, $hit)
	 }
      }
      
      return \@keepers;
   }
   else {
      return $chunk_keepers
   }
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
   my $opt_f = shift;
   my $LOG = shift;

   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;
	
   my $chunk_number = $chunk->number();
   my $q_length = $chunk->parent_seq_length();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.blastn";

   $LOG->add_entry("STARTED", $o_file, ""); 

   $chunk->write_file($file_name);
   runBlastn($file_name,
	     $db,
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

   $LOG->add_entry("FINISHED", $o_file, "");

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset(),
			    );

   $chunk->erase_fasta_file();

   return $chunk_keepers
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
      print STDERR "re reading blast report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  blast search.\n" unless $main::quiet;
      my $dir = $out_file;
      $dir =~ s/[^\/]+$//;
      File::Path::mkpath($dir);
      $w->run($command);
   }

}
#-----------------------------------------------------------------------------
sub blastx_as_chunks {
   my $chunk         = shift;
   my $db            = shift;
   my $the_void      = shift;
   my $seq_id        = shift;
   my $blastx        = shift;
   my $eval_blastx   = shift;
   my $cpus          = shift;
   my $old_db        = shift;
   my $xdformat      = shift;
   my $alt_peptide    = shift;
   my $rank          = shift;
   my $opt_f         = shift;
   my $LOG = shift;
   my $LOG_FLAG = shift;

   #build names for files to use and copy	
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;

   my $chunk_number = $chunk->number();
    
   my ($db_old_n) = $old_db =~ /([^\/]+)$/;
   $db_old_n  =~ s/\.fasta$//;
   my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.blastx";
    
   my $t_dir = "/tmp/rank".$rank;
   File::Path::mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.blastx";
    
   $db =~ /([^\/]+)$/;
   my $tmp_db = "$t_dir/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG);

   #copy db to local tmp dir and run xdformat 
   if (! @{[<$tmp_db.xp?*>]} && (! -e $blast_finished || $opt_f) ) {
      system("cp $db $tmp_db");
      xdformat($xdformat, $tmp_db, 'blastx', $alt_peptide);
  }
   elsif (-e $blast_finished && ! $opt_f) {
      print STDERR "re reading blast report.\n" unless $main::quiet;
      print STDERR "$blast_finished\n" unless $main::quiet;
      return $blast_dir;
   }

   #call blast executable	     
   $chunk->write_file($t_file_name);

   runBlastx($t_file_name,
	     $tmp_db,
	     $o_file,
	     $blastx,
	     $eval_blastx,
	     $cpus,
	     $opt_f
	    );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_blastx{
   my $chunk      = shift;
   my $blast_dir    = shift;
   my $eval_blastx = shift;
   my $bit_blastx = shift,
   my $percov_blastx = shift;
   my $percid_blastx = shift;
   my $split_hit = shift;
   my $opt_f = shift;
   my $LOG = shift;

   my $blast_finished = $blast_dir;
   $blast_finished =~ s/\.temp_dir$//;

   #merge blast reports
   if (! -e $blast_finished || $opt_f) {
      system ("cat $blast_dir/*blast* > $blast_finished");
      rmtree ("$blast_dir");
   }

   my %params;
   $params{significance}  = $eval_blastx;
   $params{hsp_bit_min}   = $bit_blastx;
   $params{percov}        = $percov_blastx;
   $params{percid}        = $percid_blastx;
   $params{split_hit}     = $split_hit;
   
   my $chunk_keepers = Widget::blastx::parse($blast_finished,
					     \%params,
					    );

   $LOG->add_entry("FINISHED", $blast_finished, "");

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset()
			    );

   if ($chunk->p_cutoff || $chunk->m_cutoff) {
      my @keepers;
      
      foreach my $hit (@{$chunk_keepers}) {
	 if ($hit->strand('query') eq '1' && $hit->start('query') >= $chunk->p_cutoff) {
	    push (@keepers, $hit)
	 }
	 elsif ($hit->strand('query') eq '-1' && $hit->start('query') >= $chunk->m_cutoff) {
	    push (@keepers, $hit)
	 }
      }
      
      return \@keepers;
   }
   else {
      return $chunk_keepers;
   }
}
#-----------------------------------------------------------------------------

sub blastx {
   my $chunk      = shift;
   my $db         = shift;
   my $the_void    = shift;
   my $seq_id     = shift;
   my $blastx = shift;
   my $eval_blastx = shift;
   my $bit_blastx = shift;
   my $percov_blastx = shift;
   my $percid_blastx = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $opt_f = shift;
   my $LOG = shift;

   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;

   my $q_length = $chunk->parent_seq_length();
   my $chunk_number = $chunk->number();
		
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.blastx";

   $LOG->add_entry("STARTED", $o_file, ""); 

   $chunk->write_file($file_name);
   runBlastx($file_name,
	     $db,
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

   $LOG->add_entry("FINISHED", $o_file, "");

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset(),
			    );
   
   $chunk->erase_fasta_file();
   
   if ($chunk->p_cutoff || $chunk->m_cutoff) {
      my @keepers;
      
      foreach my $hit (@{$chunk_keepers}) {
	 if ($hit->strand('query') eq '1' && $hit->start('query') >= $chunk->p_cutoff) {
	    push (@keepers, $hit)
	 }
	 elsif ($hit->strand('query') eq '-1' && $hit->start('query') >= $chunk->m_cutoff) {
	    push (@keepers, $hit)
	 }
      }
      
      return \@keepers;
   }
   else {
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

   if (-e $out_file  && !  $opt_f) {
      print STDERR "re reading blast report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  blast search.\n" unless $main::quiet;
      my $dir = $out_file;
      $dir =~ s/[^\/]+$//;
      File::Path::mkpath($dir);
      $w->run($command);
   }
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
   my $LOG = shift;

   my $chunk_number = $chunk->number();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.out";
   my $q_length = $chunk->parent_seq_length();
   my $query_def = $chunk->parent_def();
   my $query_seq = $chunk->seq();

   $LOG->add_entry("STARTED", $o_file, ""); 

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
						      $q_length
						     );

   $LOG->add_entry("FINISHED", $o_file, ""); 
  
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
      $command .= " $q_file -lib $rmlib -dir $dir -pa $cpus";    
   } else {
      $command .= " $q_file -species $species -dir $dir -pa $cpus";
   }
   $command .= " -nolow" if defined($no_low);
	
   my $w = new Widget::RepeatMasker();
   if (-e $o_file && ! $opt_f) {
      print STDERR "re reading repeat masker report.\n" unless $main::quiet;
      print STDERR "$o_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  repeat masker.\n" unless $main::quiet;
      $w->run($command);
   }
}

#-----------------------------------------------------------------------------
sub build_the_void {
   my $seq_id  = shift;
   my $out_dir = shift;

   $out_dir =~ s/\/$//;

   my $vid = "theVoid\.$seq_id";   
   my $the_void = "$out_dir/$vid";
   File::Path::mkpath ($the_void);

   return $the_void;
}

#-----------------------------------------------------------------------------
sub load_control_files {
   my @ctlfiles = @{shift @_};
   my %OPT = %{shift @_};

   my %CTL_OPTIONS;
   my %OK_FIELDS;

   my @MAKER_OPTS_PARAMS = ('genome',
			    'est',
			    'protein',
			    'repeat_protein',
			    'rmlib',
			    'rm_gff',
			    'predictor',
			    'snaphmm',
			    'augustus_species',
			    'model_org',
			    'max_dna_len',
			    'min_contig',
			    'split_hit',
			    'snap_flank',
			    'te_remove',
			    'single_exon',
			    'use_seq_dir',
			    'clean_up',
			    'cpus',
			    'alt_peptide'
			   );
    
   my @MAKER_BOPTS_PARAMS = ('percov_blastn',
			     'percid_blastn',
			     'eval_blastn',
			     'bit_blastn',
			     'percov_blastx',
			     'percid_blastx',
			     'eval_blastx',
			     'bit_blastx',
			     'e_perc_cov',
			     'ep_score_limit',
			     'en_score_limit'
			    );

   my @MAKER_EXE_PARAMS = ('xdformat',
			   'blastn',
			   'blastx',
			   'snap',
			   'augustus',
			   'RepeatMasker',
			   'exonerate',
			  );


   foreach my $attr (@MAKER_OPTS_PARAMS, @MAKER_BOPTS_PARAMS, @MAKER_EXE_PARAMS) {
      $OK_FIELDS{$attr}++;
   }


   #set default values for certain control options
   $CTL_OPTIONS{'clean_up'} = 0;
   $CTL_OPTIONS{'max_dna_len'} = 100000;
   $CTL_OPTIONS{'min_contig'} = 10000;
   $CTL_OPTIONS{'percov_blastn'} = 0.80;
   $CTL_OPTIONS{'percid_blastn'} = 0.85;
   $CTL_OPTIONS{'eval_blastn'} = 1e-10;
   $CTL_OPTIONS{'bit_blastn'} = 40;
   $CTL_OPTIONS{'percov_blastx'} = 0.50;
   $CTL_OPTIONS{'percid_blastx'} = 0.40;
   $CTL_OPTIONS{'eval_blastx'} = 1e-6;
   $CTL_OPTIONS{'bit_blastx'} = 30;
   $CTL_OPTIONS{'e_perc_cov'} = 50;
   $CTL_OPTIONS{'alt_peptide'} = 'c';
   $CTL_OPTIONS{'en_score_limit'} = 20;
   $CTL_OPTIONS{'ep_score_limit'} = 20;
   $CTL_OPTIONS{'predictor'} = 'snap';
   
   #load values from control files
   foreach my $ctlfile (@ctlfiles) {
      open (CTL, "< $ctlfile") or die"ERROR: Could not open control file \"$ctlfile\".\n";
	
      while (my $line = <CTL>) {
	 chomp($line);
	    
	 if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:]+)\:([^\s\t\n]+)/) {
	    my $key = $1;
	    my $value = $2;
	    if (exists $OK_FIELDS{$key}) {
	       if ($value =~ /\$/) {
		  $value = `echo $value`;
		  chomp($value);
	       }
	       $CTL_OPTIONS{$key} = $value unless (not defined $value);
	    }
	    else {
	       warn "ERROR: Invalid option \"$key\" in control file $ctlfile\n";
	    }
	 }
      }
   }

   #use command line defined genome
   $CTL_OPTIONS{'genome'} = $OPT{g} if (defined $OPT{g});
    
   #use command line defined predictor
   $CTL_OPTIONS{'predictor'} = $OPT{predictor} if (defined $OPT{predictor});
   die "ERROR: Invalid predictor defined: $CTL_OPTIONS{'predictor'}\n".
       "Must be set to 'snap', 'augustus, or 'es2genome'\n"
       unless ($CTL_OPTIONS{'predictor'} =~ /snap|augustus|est2genome/);

   #validate required values from control files
   my @infiles = ('genome', 'protein', 'est', 'xdformat', 'blastn',
		  'blastx', 'exonerate', 'snap'
		 );

   #sometimes required
   push (@infiles, 'repeat_protein') if ($CTL_OPTIONS{te_remove});
   push (@infiles, 'RepeatMasker') unless($OPT{R} || $OPT{GFF});
   push (@infiles, 'rm_gff') if ($OPT{GFF});
   push (@infiles, 'augustus') if ($CTL_OPTIONS{predictor} eq 'augustus' || $CTL_OPTIONS{'augustus'});

   my $error;

   foreach my $in (@infiles) {
      if (not $CTL_OPTIONS{$in}) {
	 $error .= "You have failed to provide a value for \'$in\' in the control files\n";
	 next;
      }

      if (not -e $CTL_OPTIONS{$in}) {
	 $error .= "The \'$in\' file $CTL_OPTIONS{$in} does not exist.\n".
	 "Please check your control files: maker_opts.ctl, maker_bopts, or maker_exe.ctl\n";
	 next;
      }

      #set the absolute path to the file to reduce ambiguity
      $CTL_OPTIONS{$in} = Cwd::abs_path($CTL_OPTIONS{$in}) unless ($in =~ /^blastn$|^blastx$/);
   }

   die $error if (defined $error);

   if (! $OPT{R} && ! $CTL_OPTIONS{'model_org'}) {
      warn "There is no model specified for RepeatMasker in maker_opts.ctl : model_org.\n".
           "As a result the default (drosophila) will be used.\n";
      $CTL_OPTIONS{'model_org'} = "drosophila";
   }
   if ( ($CTL_OPTIONS{'predictor'} eq 'augustus' || $CTL_OPTIONS{'augustus'}) &&
	not $CTL_OPTIONS{'augustus_species'}
      ) {
      warn "There is no species specified for Augustus in maker_opts.ctl : augustus_species.\n".
           "As a result the default (fly) will be used.\n";
      $CTL_OPTIONS{'augustus_species'} = "fly";
   }
   if ( ($CTL_OPTIONS{'predictor'} eq 'augustus' || $CTL_OPTIONS{'augustus'}) &&
        (! $ENV{AUGUSTUS_CONFIG_PATH} || ! -e "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg")
	) {
       die "ERROR: The environmental variable AUGUSTUS_CONFIG_PATH has not been set or is not set correctly\n",
           "Please set this in your profile per Augustus installation instructions\n";
   }
   if (not $CTL_OPTIONS{'snaphmm'}) {
      warn "There is no model specified for for Snap in maker_opts.ctl : snaphmm.\n".
           "As a result, the default (fly) will be used.\n";
      $CTL_OPTIONS{'snaphmm'} = "fly";
   }
   if (! -e $CTL_OPTIONS{'snaphmm'} &&
       (! exists $ENV{'ZOE'} || ! -e $ENV{'ZOE'}."/HMM/".$CTL_OPTIONS{'snaphmm'})
      ) {
      
      die "ERROR: The snaphmm specified for Snap in maker_opts.ctl does not exist.\n";
   }
   if ($CTL_OPTIONS{'max_dna_len'} < 50000) {
      warn "ERROR: max_dna_len is set too low.  The minimum value permited is 50,000\n".
           "max_dna_len wil be reset to 50,000\n\n";
      $CTL_OPTIONS{'max_dna_len'} = 50000;
   }

   #set values for datastructure    
   $CTL_OPTIONS{'genome'} =~ /([^\/]+)$/;
   $CTL_OPTIONS{'out_name'} = $1;
   $CTL_OPTIONS{'out_name'} =~ s/\.[^\.]+$//;
   $CTL_OPTIONS{'out_base'} = Cwd::cwd();

   if ($CTL_OPTIONS{'use_seq_dir'}) {
      my @file_struct = split(/\//, $CTL_OPTIONS{'genome'});
      pop @file_struct;
      $CTL_OPTIONS{'out_base'} = join("/", @file_struct);
   }
    
   if (not $CTL_OPTIONS{'out_base'}) {
      die "No working directory, check your use_seq_dir option\n";
   }

   return %CTL_OPTIONS;
}

#-----------------------------------------------------------------------------
sub generate_control_files {
   #--build maker_opts.ctl file
   my $repeat_protein = Cwd::abs_path("$FindBin::Bin/../data/te_proteins.fasta");

   open (OUT, "> maker_opts.ctl");
   print OUT "#-----sequence and library files\n";
   print OUT "genome: #genome sequence file (required)\n";
   print OUT "est: #EST sequence file (required)\n";
   print OUT "protein:  #protein sequence file (required)\n";
   print OUT "repeat_protein:$repeat_protein #a database of transposable element proteins\n";
   print OUT "rmlib: #an organism specific repeat library (optional)\n";
   print OUT "rm_gff: #a gff3 format file of repeat elements (only used with -GFF flag)\n";
   print OUT "\n";
   print OUT "#-----external application specific options\n";
   print OUT "snaphmm:fly #SNAP HMM model\n";
   print OUT "augustus_species:fly #Augustus gene prediction model\n";
   print OUT "model_org:all #RepeatMasker model organism\n";
   print OUT "alt_peptide:c #amino acid used to replace non standard amino acids in xdformat\n";
   print OUT "cpus:1 #max number of cpus to use in BLAST and RepeatMasker\n";
   print OUT "\n";
   print OUT "#-----Maker specific options\n";
   print OUT "predictor:snap #identifies which gene prediction program to use for annotations\n";
   print OUT "te_remove:1 #mask regions with excess similarity to transposable element proteins\n";
   print OUT "max_dna_len:100000 #length for dividing up contigs into chunks (larger values increase memory usage)\n";
   print OUT "min_contig:10000 #all contigs from the input genome file below this size are skipped\n";
   print OUT "split_hit:10000 #length of the splitting of hits (max intron size for EST and protein alignments)\n";
   print OUT "snap_flank:200 #length of sequence surrounding EST and protein evidence used to extend gene predictions\n";
   print OUT "single_exon:0 #consider EST hits aligning to single exons when generating annotations, 1 = yes, 0 = no\n";
   print OUT "use_seq_dir:1 #place output files in same directory as sequence file: 1 = yes, 0 = no\n";
   print OUT "clean_up:0 #remove theVoid directory: 1 = yes, 0 = no\n";
   close (OUT);

   #--build maker_bopts.ctl file
   open (OUT, "> maker_bopts.ctl");
   print OUT "#-----BLAST and Exonerate statistics thresholds\n";
   print OUT "percov_blastn:0.80 #Blastn Percent Coverage Threhold EST-Genome Alignments\n";
   print OUT "percid_blastn:0.85 #Blastn Percent Identity Threshold EST-Genome Aligments\n";
   print OUT "eval_blastn:1e-10 #Blastn eval cutoff\n";
   print OUT "bit_blastn:40 #Blastn bit cutoff\n";
   print OUT "percov_blastx:0.50 #Blastx Percent Coverage Threhold Protein-Genome Alignments\n";
   print OUT "percid_blastx:0.40 #Blastx Percent Identity Threshold Protein-Genome Aligments\n";
   print OUT "eval_blastx:1e-6 #Blastx eval cutoff\n";
   print OUT "bit_blastx:30 #Blastx bit cutoff\n";
   print OUT "e_perc_cov:50 #Exonerate Percent Coverage Thresshold EST_Genome Alignments\n";
   print OUT "ep_score_limit:20 #Report  alignments scoring at least this percentage of the maximal score exonerate nucleotide\n";
   print OUT "en_score_limit:20 #Report  alignments scoring at least this percentage of the maximal score exonerate protein\n";
   close(OUT);

   #--build maker_exe.ctl file
   my %executables = ( xdformat => '',
		       blastn => '',
		       blastx => '',
		       snap => '',
		       augustus => '',
		       exonerate => '',
		       RepeatMasker => '',
		       exonerate => ''
		     );

   while (my $exe = each %executables) {
      my $loc = `which $exe`;
      chomp $loc;
      $executables{$exe} = $loc unless ($loc =~ /^no $exe/);
   }

   open (OUT, "> maker_exe.ctl");
   print OUT "#-----Location of executables required by Maker\n";
   print OUT "xdformat:".$executables{xdformat}." #location of xdformat executable\n";
   print OUT "blastn:".$executables{blastn}." #location of blastn executable\n";
   print OUT "blastx:".$executables{blastx}." #location of blastn executable\n";
   print OUT "snap:".$executables{snap}." #location of snap executable\n";
   print OUT "augustus:".$executables{augustus}." #location of augustus executable (optional)\n";
   print OUT "RepeatMasker:".$executables{RepeatMasker}." #location of RepeatMasker executable\n";
   print OUT "exonerate:".$executables{exonerate}." #location of exonerate executable\n";
   close(OUT);
}

1;
