#------------------------------------------------------------------------
#----                                GI                              ----
#------------------------------------------------------------------------
package GI;

use strict;
use vars qw(@ISA @EXPORT $VERSION $TMP);
use Exporter;
use FileHandle;
use File::Temp qw(tempfile tempdir);
use Dumper::GFF::GFFV3;
use Dumper::XML::Game;
use URI::Escape;
use File::Path;
use File::Copy;
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
use Widget::fgenesh;
use Widget::xdformat;
use Widget::formatdb;
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

$TMP = tempdir("maker_XXXXXX", CLEANUP => 1, TMPDIR => 1);
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
      my $s_end = $pred->end('query');

      ($s_start, $s_end) = ($s_end, $s_start) if ($s_end < $s_start);

      if ($c_start <= $s_start && $s_start <= $c_end) {
         push (@keepers, $pred);
      }
   }

   return \@keepers;
}
#-----------------------------------------------------------------------------
sub merge_resolve_hits{
   my $fasta = shift @_;
   my $fasta_t_index = shift @_;
   my $fasta_p_index = shift @_;
   my $fasta_a_index = shift @_;
   my $blastn_keepers = shift @_;
   my $blastx_keepers = shift @_;
   my $tblastx_keepers = shift @_;
   my $blastn_holdovers = shift @_;
   my $blastx_holdovers = shift @_;
   my $tblastx_holdovers = shift @_;
   my $the_void = shift @_;
   my %CTL_OPT = %{shift @_};
   my $LOG = shift @_;

   PhatHit_utils::merge_hits($blastn_keepers,  
			     $blastn_holdovers, 
			     $CTL_OPT{split_hit},
			    );
   @{$blastn_holdovers} = ();

   PhatHit_utils::merge_hits($blastx_keepers,  
			     $blastx_holdovers, 
			     $CTL_OPT{split_hit},
			    );
   @{$blastx_holdovers} = ();

   PhatHit_utils::merge_hits($tblastx_keepers,
                             $tblastx_holdovers,
			     $CTL_OPT{split_hit},
			    );
   @{$tblastx_holdovers} = ();

   $blastn_keepers = reblast_merged_hits($fasta,
					 $blastn_keepers,
					 $fasta_t_index,
					 $the_void,
					 'blastn',
					 \%CTL_OPT,
					 $LOG
					);

   $blastx_keepers = reblast_merged_hits($fasta,
					 $blastx_keepers,
					 $fasta_p_index,
					 $the_void,
					 'blastx',
					 \%CTL_OPT,
					 $LOG
					);

   $tblastx_keepers = reblast_merged_hits($fasta,
					  $tblastx_keepers,
					  $fasta_a_index,
					  $the_void,
					  'tblastx',
					  \%CTL_OPT,
					  $LOG
					 );

   return ($blastn_keepers, $blastx_keepers, $tblastx_keepers);
}
#-----------------------------------------------------------------------------
sub reblast_merged_hits {
   my $g_fasta     = shift @_;
   my $hits        = shift @_;
   my $db_index    = shift @_;
   my $the_void    = shift @_;
   my $type        = shift @_;
   my %CTL_OPT = %{shift @_};
   my $LOG         = shift @_;

   #==get data from parent fasta

   #parent fasta get def and seq
   my $par_def = Fasta::getDef($$g_fasta);
   my $par_seq = Fasta::getSeqRef($$g_fasta);

   #get seq id off def line
   my ($p_id)  = $par_def =~ /^([^\s\t\n]+)/;
   $p_id =~ s/^\>//g;		#just in case
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
      my $piece = Shadower::getPieces($par_seq, \@coors, $CTL_OPT{split_hit});
      
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
      
      #get name
      my $t_id  = $hit->name();
      $t_id =~ s/\|/_/g;

      #build a safe name for file names from the sequence identifier
      my $t_safe_id = uri_escape($t_id, 
				 '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				);
      #search db index
      my $fastaObj = $db_index->get_Seq_by_id($hit->name);
      if (not $fastaObj) {
	 #rebuild index and try again
	 my $db_file = $db_index->{dirname} ."/". $db_index->{offsets}->{__file_0};
	 $db_index = Bio::DB::Fasta->new($db_file, '-reindex' => 1);
	 $fastaObj = $db_index->get_Seq_by_id($hit->name);
	 if (not $fastaObj) {
	    print STDERR "stop here:".$hit->name."\n";
	    die "ERROR: Fasta index error\n";
	 }
      }
      
      #get fasta def and seq
      my $t_seq      = $fastaObj->seq();
      my $t_def      = $db_index->header($hit->name);

      #write fasta file
      my $fasta = Fasta::toFasta('>'.$t_def, \$t_seq);
      my $t_file = $the_void."/".$t_safe_id.'.for_'.$type.'.fasta';
      FastaFile::writeFile($fasta, $t_file);
      
      #build db for blast using xdformat or formatdb
      dbformat($CTL_OPT{_formater}, $t_file, $type);
      
      #==run the blast search
      if ($type eq 'blastx') {
	  
	 print STDERR "re-running blast against ".$hit->name."...\n" unless $main::quiet;
	 my $keepers = blastx($chunk, 
			      $t_file,
			      $the_void,
			      $p_safe_id.".".$piece->[0]->{b}.".".$piece->[0]->{e},
			      $CTL_OPT{_blastx},
			      $CTL_OPT{eval_blastx},
			      $CTL_OPT{bit_blastx},
			      $CTL_OPT{pcov_blastx},
			      $CTL_OPT{pid_blastx},
			      $CTL_OPT{split_hit},
			      $CTL_OPT{cpus},
			      $CTL_OPT{force},
			      $LOG
			     );
	  
	 push(@blast_keepers, @{$keepers});
	 print STDERR "...finished\n" unless $main::quiet;
      }
      elsif ($type eq 'blastn') {
	  
	 print STDERR "re-running blast against ".$hit->name."...\n" unless $main::quiet;
	 my $keepers = blastn($chunk, 
			      $t_file,
			      $the_void,
			      $p_safe_id.".".$piece->[0]->{b}.".".$piece->[0]->{e},
			      $CTL_OPT{_blastn},
			      $CTL_OPT{eval_blastn},
			      $CTL_OPT{bit_blastn},
			      $CTL_OPT{pcov_blastn},
			      $CTL_OPT{pid_blastn},
			      $CTL_OPT{split_hit},
			      $CTL_OPT{cpus},
			      $CTL_OPT{force},
			      $LOG
			     );
	  
	 push(@blast_keepers, @{$keepers});
	 print STDERR "...finished\n" unless $main::quiet;
      }
      elsif ($type eq 'tblastx') {

	 print STDERR "re-running blast against ".$hit->name."...\n" unless $main::quiet;
	 my $keepers = tblastx($chunk,
			       $t_file,
			       $the_void,
			       $p_safe_id.".".$piece->[0]->{b}.".".$piece->[0]->{e},
			       $CTL_OPT{_tblastx},
			       $CTL_OPT{eval_tblastx},
			       $CTL_OPT{bit_tblastx},
			       $CTL_OPT{pcov_tblastx},
			       $CTL_OPT{pid_tblastx},
			       $CTL_OPT{split_hit},
			       $CTL_OPT{cpus},
			       $CTL_OPT{force},
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
      #fix coordinates for hits outside of chunk end   
      $p_coor->[0] = $chunk->length if($p_coor->[0] > $chunk->length);
      $p_coor->[1] = $chunk->length if($p_coor->[1] > $chunk->length);
      #fix coordinates for hits outside of chunk begin
      $p_coor->[0] = 0 if($p_coor->[0] < 0);
      $p_coor->[1] = 0 if($p_coor->[1] < 0);
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

   #hit holdovers and keepers are returned in same order given by user
   return @holdovers, @keepers;
}

#-----------------------------------------------------------------------------
# sub write_quality_data {
#    my $quality_indices = shift;
#    my $seq_id          = shift;

#    my $out_file = $seq_id.'.maker.transcripts.qi';
#    my $fh = new FileHandle();
#    $fh->open(">$out_file");

#    print $fh "genomic_seq\ttranscript\tquality_index\n";

#    while (my $d = shift(@{$quality_indices})) {
#       my $t_name = $d->[0];
#       my $t_qi   = $d->[1];
	
#       print $fh "$seq_id\t$t_name\t$t_qi\n";
#    }
#    $fh->close();
# }
#-----------------------------------------------------------------------------
sub abinit_p_and_t_fastas {
   my $preds = shift;
   my $id = shift;
   my $seq_ref = shift;
   my $out_dir = shift;

   my %fhs;

   foreach my $hit (@{$preds}) {
      my $source = $hit->algorithm;
      $source = uri_escape($source,
			   '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
			  );
      my $t_name = $hit->name(); # note this is being set in GFFV3::pred_data
      my $t_seq  = maker::auto_annotator::get_transcript_seq($hit, $seq_ref);	
		
      my ($p_seq, $offset, $end) = 
      maker::auto_annotator::get_translation_seq($t_seq);
		
      my $score = 0;
      foreach my $hsp ($hit->hsps) {
	 $score += $hsp->score();
      }
      my $t_def = ">$t_name transcript offset:$offset score:$score";
      my $p_def = ">$t_name protein score:$score";
      
      if (! exists $fhs{$source}) {
	 open(my $t_fh, "> $out_dir/$id.maker.$source.transcripts.fasta");
	 open(my $p_fh, "> $out_dir/$id.maker.$source.proteins.fasta");
	 $fhs{$source}{t} = $t_fh;
	 $fhs{$source}{p} = $p_fh;
      }

      my $t_fh = $fhs{$source}{t};
      my $p_fh = $fhs{$source}{p};

      print $t_fh Fasta::toFasta($t_def, \$t_seq);
      print $p_fh Fasta::toFasta($p_def, \$p_seq);
   }
}
#-----------------------------------------------------------------------------
sub maker_p_and_t_fastas {
   my $annotations = shift @_;
   
   my $p_fastas = '';
   my $t_fastas = '';
   
   foreach my $an (@$annotations) {
      foreach my $a (@{$an->{t_structs}}) {
	 my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	 $p_fastas .= $p_fasta;
	 $t_fastas .= $t_fasta;
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
   my $CTL_OPT = shift @_;
   my $mpi_size = shift@_ || 1;

   ($CTL_OPT->{_protein}, $CTL_OPT->{p_db}) = split_db($CTL_OPT->{protein}, $mpi_size, $CTL_OPT->{alt_peptide});
   ($CTL_OPT->{_est}, $CTL_OPT->{e_db}) = split_db($CTL_OPT->{est}, $mpi_size);
   ($CTL_OPT->{_est_reads},  $CTL_OPT->{d_db}) = split_db($CTL_OPT->{est_reads}, $mpi_size);
   ($CTL_OPT->{_altest},  $CTL_OPT->{a_db}) = split_db($CTL_OPT->{altest}, $mpi_size);
   ($CTL_OPT->{_repeat_protein}, $CTL_OPT->{r_db}) = split_db($CTL_OPT->{repeat_protein}, $mpi_size, $CTL_OPT->{alt_peptide});
}
#----------------------------------------------------------------------------
sub split_db {
   my $file = shift @_;
   my $mpi_size = shift @_ || 1;
   my $alt = shift @_;

   return ('', []) if (not $file);

   my $fasta_iterator = new Iterator::Fasta($file);
   my $db_size = $fasta_iterator->number_of_entries();
   my $bins = $mpi_size;
   $bins = $db_size if ($db_size < $bins);

   my @fhs;
   my @db_files;

   my ($f_name) = $file =~ /([^\/]+$)/;
   $f_name =~ s/\.fasta$//;
    
   my $d_name = "$f_name\.mpi\.$mpi_size";
   my $b_dir = cwd(). "/mpi_blastdb";
   my $f_dir = "$b_dir/$d_name";
   my $t_dir = $TMP."/$d_name";
   my $t_full = $t_dir.".fasta";
   my $f_full = $f_dir.".fasta";

   if (-e "$f_dir") {
      my @t_db = <$f_dir/*$d_name\.*>;

      foreach my $f (@t_db) {
	 push (@db_files, $f) if (! -d $f);
      }
      
      return $f_full, \@db_files;
   }

   mkdir($t_dir);
   mkdir($b_dir) unless (-e $b_dir);

   for (my $i = 0; $i < $bins; $i++) {
      my $name = "$t_dir/$d_name\.$i\.fasta";
      my $fh;
      open ($fh, "> $name");

      push (@fhs, $fh);
   }
   
   my $FA;
   open($FA, "> $t_full");
   while (my $fasta = $fasta_iterator->nextEntry()) {
      my $def = Fasta::getDef(\$fasta);
      my $seq_ref = Fasta::getSeqRef(\$fasta);

      #fix non standard peptides
      if (defined $alt) {
	 $$seq_ref =~ s/\*//g;
	 $$seq_ref =~ s/[^abcdefghiklmnpqrstvwyzxABCDEFGHIKLMNPQRSTVWYZX\-\n]/$alt/g;
      }

      #Skip empty fasta entries
      next if($$seq_ref eq '');

      #reformat fasta, just incase
      my $fasta_ref = Fasta::toFastaRef($def, $seq_ref);

      my $fh = shift @fhs;
      print $fh $$fasta_ref;
      print $FA $$fasta_ref;
      push (@fhs, $fh);
   }
   close($FA);

   foreach my $fh (@fhs) {
      close ($fh);
   }

   system("mv $t_full $f_full");
   system("mv $t_dir $f_dir");

   if (-e "$f_dir") {
      my @t_db = <$f_dir/*$d_name\.*>;

      foreach my $f (@t_db) {
	 push (@db_files, $f) if (! -d $f);
      }

      return $f_full, \@db_files;
   }
   else {
      die "ERROR: Could not split db\n";
   }
}
#----------------------------------------------------------------------------
# sub load_anno_hsps {
#    my $annotations = shift;
#    my @coors;
#    my $i = @{$annotations};
#    foreach my $an (@$annotations) {
#       foreach my $a (@{$an->[0]}) {
# 	 my $hit = $a->{hit};
# 	 foreach my $hsp ($hit->hsps()) {
# 	    push(@coors, [$hsp->nB('query'),
# 			  $hsp->nE('query'),
# 			 ]);
# 	 }
#       }
#    }
#    return` (\@coors, $i);;
# }
#-----------------------------------------------------------------------------
# sub load_clust_hsps {
#    my $clusters = shift;
#    my @coors;
#    my $i = @{$clusters};
#    foreach my $c (@$clusters) {
#       foreach my $hit (@{$c}) {
# 	 foreach my $hsp ($hit->hsps()) {
# 	    push(@coors, [$hsp->nB('query'),
# 			  $hsp->nE('query'),
# 			 ]);
# 	 }
#       }
#    }
#    return (\@coors, $i);
# }
#-----------------------------------------------------------------------------
# sub load_snap_hsps {
#    my $snaps = shift;
#    my @coors;
#    my $i = @{$snaps};
#    foreach my $hit (@{$snaps}) {
#       foreach my $hsp ($hit->hsps()) {
# 	 push(@coors, [$hsp->nB('query'),
# 		       $hsp->nE('query'),
# 		      ]);
#       }
#    }
#    return (\@coors, $i);
# }
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
#------------------------------------------------------------------------
sub combine {
   my @bag;
   while (my $hits = shift(@_)) {
      foreach my $hit (@{$hits}) {
	 push(@bag, $hit);
      }
   }
   return \@bag;
}
#-----------------------------------------------------------------------------
sub abinits {
   my @preds;

   push(@preds, @{snap(@_)});
   push(@preds, @{augustus(@_)});
   push(@preds, @{fgenesh(@_)});
   push(@preds, @{twinscan(@_)});

   return \@preds;
}
#-----------------------------------------------------------------------------
sub snap {
   my $in_file     = shift;
   my $the_void    = shift;
   my $seq_id      = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;
   
   my $exe    = $CTL_OPT->{snap};
   my $snaphmm = $CTL_OPT->{snaphmm};
   
   return [] if(not $exe);
	
   my %params;
   my $out_file = "$the_void/$seq_id\.all\.snap";

   $LOG->add_entry("STARTED", $out_file, "");   

   my $command  = $exe;
   $command .= " $snaphmm";
   $command .= " $in_file";
   $command .= " > $out_file";
	
   my $w = new Widget::snap();
	
   if (-e $out_file && ! $CTL_OPT->{force}) {
      print STDERR "re reading snap report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  snap.\n" unless $main::quiet;
      $w->run($command);
   }
	
   $params{min_exon_score}  = -100000; #-10000;
   $params{min_gene_score}  = -100000; #0;
		
   my $keepers = Widget::snap::parse($out_file,
				     \%params,
				     $in_file,
				    );

   $LOG->add_entry("FINISHED", $out_file, "");

   return $keepers;
}
#-----------------------------------------------------------------------------
sub augustus {
   my $in_file     = shift;
   my $the_void    = shift;
   my $seq_id      = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;

   my $exe = $CTL_OPT->{augustus};
   my $org = $CTL_OPT->{augustus_species};

   return [] if(not $exe);

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.augustus";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = $exe;
   $command .= ' --species='."$org";
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::augustus();

   if (-e $out_file && ! $CTL_OPT->{force}) {
      print STDERR "re reading augustus report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  augustus.\n" unless $main::quiet;
      $w->run($command);
   }

   $params{min_exon_score}  = -100000; #-10000;
   $params{min_gene_score}  = -100000; #0;

   my $chunk_keepers = Widget::augustus::parse($out_file,
					       \%params,
					       $in_file
					      );

   $LOG->add_entry("FINISHED", $out_file, "");

   return $chunk_keepers;
}
#-----------------------------------------------------------------------------
sub fgenesh {
   my $in_file     = shift;
   my $the_void    = shift;
   my $seq_id      = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;

   my $exe = $CTL_OPT->{fgenesh};
   my $org = $CTL_OPT->{fgenesh_par_file};

   return [] if(not $exe);

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.fgenesh";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = $exe;
   $command .= " $org";
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::fgenesh();

   if (-e $out_file && ! $CTL_OPT->{force}) {
      print STDERR "re reading fgenesh report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  fgenesh.\n" unless $main::quiet;
      $w->run($command);
   }

   $params{min_exon_score}  = -100000; #-10000;
   $params{min_gene_score}  = -100000; #0;

   my $chunk_keepers = Widget::fgenesh::parse($out_file,
					      \%params,
					      $in_file
					     );

   $LOG->add_entry("FINISHED", $out_file, "");

   return $chunk_keepers;
}
#-----------------------------------------------------------------------------
sub twinscan {
   my $in_file     = shift;
   my $the_void    = shift;
   my $seq_id      = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;

   my $exe = $CTL_OPT->{twinscan};
   my $org = $CTL_OPT->{twinscan_species};

   return [] if(not $exe);

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.twinscan";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = $exe;
   $command .= ' --species='."$org";
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::twinscan();

   if (-e $out_file && ! $CTL_OPT->{force}) {
      print STDERR "re reading twinscan report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  twinscan.\n" unless $main::quiet;
      $w->run($command);
   }

   $params{min_exon_score}  = -100000; #-10000;
   $params{min_gene_score}  = -100000; #0;

   my $chunk_keepers = Widget::twinscan::parse($out_file,
					       \%params,
					       $in_file
					      );

   $LOG->add_entry("FINISHED", $out_file, "");

   return $chunk_keepers;
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
   my $pcov              = shift;
   my $pid               = shift;
   my $score_limit       = shift;
   my $matrix            = shift;
   my $opt_f             = shift;
   my $LOG               = shift;

   my $def = Fasta::getDef($g_fasta);
   my $seq = Fasta::getSeqRef($g_fasta);
	
   my $exe = $exonerate;
	
   my @exonerate_clusters;
   my $i = 0;
   foreach my $c (@{$phat_hit_clusters}) {
      my $n = 0;
      my $got_some = 0;

      foreach my $hit (@{$c}) {
	 last if $n == $depth;

	 next if $hit->pAh < $pcov;
	 next if $hit->hsp('best')->frac_identical < $pid;
	 
	 my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit,'query');

	 my @coors = [$nB, $nE];
	 my $p = Shadower::getPieces($seq, \@coors, 50);
	 my $p_def = $def." ".$p->[0]->{b}." ".$p->[0]->{e};
	 my $p_fasta = Fasta::toFasta($p_def, \$p->[0]->{piece});
	 my $name =  Fasta::def2SeqID($p_def);
	 my $safe_name = Fasta::seqID2SafeID($name);

	 my $d_file = $the_void."/".$safe_name.'.'.$p->[0]->{b}.'.'.$p->[0]->{e}.".fasta";

	 FastaFile::writeFile($p_fasta, $d_file);

	 my $offset = $p->[0]->{b} - 1;
	 my $id  = $hit->name();

	 my $fastaObj = $db_index->get_Seq_by_id($hit->name);
	 if (not $fastaObj) {
	    #rebuild index and try again
	    my $db_file = $db_index->{dirname} ."/".$db_index->{offsets}->{__file_0};
	    $db_index = Bio::DB::Fasta->new($db_file, '-reindex' => 1);
	    $fastaObj = $db_index->get_Seq_by_id($hit->name);
	    if (not $fastaObj) {
	       print STDERR "stop here:".$hit->name."\n";
	       die "ERROR: Fasta index error\n";
	    }
	 }
	 my $seq      = $fastaObj->seq();
	 my $def      = $db_index->header($hit->name);

	 my $fasta    = Fasta::toFasta('>'.$def, \$seq);

	 #build a safe name for file names from the sequence identifier
	 my $safe_id = Fasta::seqID2SafeID($id);

	 my $t_file    = $the_void."/".$safe_id.'.fasta';
	 FastaFile::writeFile($fasta, $t_file);

	 my $exonerate_hits = to_polisher($d_file,
					  $t_file,
					  $the_void,
					  $offset,
					  $type,
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
      }
      elsif ($q_str =~ /Target Intron/) {
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
						  $exe,
						  $score_limit,
						  $matrix,
						  $opt_f,
						  $LOG
						 );
   }
   elsif ($type eq 'e') {
      return polisher::exonerate::est::polish($d_file,
					      $t_file,
					      $the_void,
					      $offset,
					      $exe,
					      $score_limit,
					      $matrix,
					      $opt_f,
					      $LOG
					     );
   }
   else {
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
sub dbformat {
   my $command = shift;
   my $file = shift;
   my $type = shift;

   die "ERROR: Can not find xdformat or formatdb executable\n" if(! -e $command);
   die "ERROR: Can not find the db file $file\n" if(! -e $file);
   die "ERROR: You must define a type (blastn|blastx|tblastx)\n" if(! $type);


   if ($command =~ /xdformat/) {
      if (($type eq 'blastn' && ! -e $file.'.xnd') ||
	  ($type eq 'blastx' && ! -e $file.'.xpd') ||
	  ($type eq 'tblastx' && ! -e $file.'.xnd')
	 ) {
	 $command .= " -p" if($type eq 'blastx');
	 $command .= " -n" if($type eq 'blastn' || $type eq 'tblastx');
	 $command .= " $file";

	 my $w = new Widget::xdformat();
	 print STDERR "formating database...\n" unless $main::quiet;
	 $w->run($command);
      }
   }
   elsif ($command =~ /formatdb/) {
      if (($type eq 'blastn' && ! -e $file.'.nsq') ||
	  ($type eq 'blastx' && ! -e $file.'.psq') ||
	  ($type eq 'tblastx' && ! -e $file.'.nsq')
	 ) {
	 $command .= " -p T" if($type eq 'blastx');
	 $command .= " -p F" if($type eq 'blastn' || $type eq 'tblastx');
	 $command .= " -i $file";

	 my $w = new Widget::formatdb();
	 print STDERR "formating database...\n" unless $main::quiet;
	 $w->run($command);
      }
   }
   else {
      die "ERROR: databases can only be formated by xdformat or formatdb not \'$command\'\n";
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
   my $split_hit = shift;
   my $cpus = shift;
   my $old_db = shift;
   my $formater = shift;
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

   my $t_dir = $TMP."/rank".$rank;
   File::Path::mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.blastn";

   $db =~ /([^\/]+)$/;
   my $tmp_db = "$t_dir/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #copy db to local tmp dir and run xdformat or formatdb 
   if (! @{[<$tmp_db.xn?*>]} && (! -e $blast_finished || $opt_f) ) {
      copy($db, $tmp_db);
      dbformat($formater, $tmp_db, 'blastn');
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
	     $split_hit,
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
   my $pcov_blastn = shift;
   my $pid_blastn = shift;
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
   $params{percov}        = $pcov_blastn;
   $params{percid}        = $pid_blastn;
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
   my $pcov_blastn = shift;
   my $pid_blastn = shift;
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
	     $split_hit,
	     $cpus,
	     $opt_f
	    );

   my %params;
   $params{significance}  = $eval_blastn;
   $params{hsp_bit_min}   = $bit_blastn;
   $params{percov}        = $pcov_blastn;
   $params{percid}        = $pid_blastn;
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
   my $blast = shift;
   my $eval_blast = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $opt_f = shift;

   my $command  = $blast;
   if ($command =~ /blastn$/) {
      $command .= " $db $q_file B=100000 V=100000 E=$eval_blast";
      $command .= " wordmask=seg";
      $command .= " R=3";
      $command .= " W=15";
      $command .= " M=1";
      $command .= " N=-3";
      $command .= " Q=3";
      $command .= " Z=1000";
      $command .= " Y=500000000";
      $command .= " cpus=$cpus";	
      $command .= " topcomboN=1";
      $command .= " hspmax=100";
      $command .= " gspmax=100";
      $command .= " hspsepqmax=$split_hit";
      $command .= " lcmask";
      $command .= " wordmask=seg";
      $command .= " gi";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p blastn";
      $command .= " -d $db -i $q_file -b 100000 -v 100000 -e $eval_blast";
      $command .= " -E 3";
      $command .= " -W 15";
      $command .= " -r 1";
      $command .= " -q -3";
      $command .= " -G 3";
      $command .= " -z 1000";
      $command .= " -Y 500000000";
      $command .= " -a $cpus";	
      $command .= " -K 100";
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   else{
      die "ERROR: Must be a blastn executable";  
   }

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
   my $split_hit     = shift;
   my $cpus          = shift;
   my $old_db        = shift;
   my $formater      = shift;
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
    
   my $t_dir = $TMP."/rank".$rank;
   File::Path::mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.blastx";
    
   $db =~ /([^\/]+)$/;
   my $tmp_db = "$t_dir/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG);

   #copy db to local tmp dir and run xdformat or format db 
   if (! @{[<$tmp_db.xp?*>]} && (! -e $blast_finished || $opt_f) ) {
       copy($db, $tmp_db);
       dbformat($formater, $tmp_db, 'blastx');
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
	     $split_hit,
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
   my $pcov_blastx = shift;
   my $pid_blastx = shift;
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
   $params{percov}        = $pcov_blastx;
   $params{percid}        = $pid_blastx;
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
   my $pcov_blastx = shift;
   my $pid_blastx = shift;
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
	     $split_hit,
	     $cpus,
	     $opt_f
	    );

   my %params;
   $params{significance}  = $eval_blastx;
   $params{hsp_bit_min}   = $bit_blastx;
   $params{percov}        = $pcov_blastx;
   $params{percid}        = $pid_blastx;
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
   my $blast = shift;
   my $eval_blast = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $opt_f = shift;

   my $command  = $blast;
   if ($command =~ /blastx$/) {
      $command .= " $db $q_file B=10000 V=10000 E=$eval_blast";
      $command .= " wordmask=seg";
      #$command .= " T=20";
      #$command .= " W=5";
      #$command .= " wink=5";
      $command .= " Z=300";
      $command .= " Y=500000000";
      $command .= " hspmax=100";
      $command .= " cpus=$cpus";
      $command .= " gspmax=100";
      #$command .= " hspsepqmax=10000";
      $command .= " lcmask";
      $command .= " kap";
      $command .= " wordmask=seg";
      $command .= " gi";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p blastx";
      $command .= " -d $db -i $q_file -b 100000 -v 100000 -e $eval_blast";
      $command .= " -z 300";
      $command .= " -Y 500000000";
      $command .= " -a $cpus";	
      $command .= " -K 100";
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   else{
      die "ERROR: Must be a blastx executable";  
   }

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
sub tblastx_as_chunks {
   my $chunk      = shift;
   my $db         = shift;
   my $the_void    = shift;
   my $seq_id     = shift;
   my $tblastx = shift;
   my $eval_tblastx = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $old_db = shift;
   my $formater = shift;
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
   my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.tblastx";

   my $t_dir = $TMP."/rank".$rank;
   File::Path::mkpath($t_dir);

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.tblastx";

   $db =~ /([^\/]+)$/;
   my $tmp_db = "$t_dir/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #copy db to local tmp dir and run xdformat or formatdb
   if (! @{[<$tmp_db.xn?*>]} && (! -e $blast_finished || $opt_f) ) {
      copy($db, $tmp_db);
      dbformat($formater, $tmp_db, 'tblastx');
   }
   elsif (-e $blast_finished && ! $opt_f) {
      print STDERR "re reading blast report.\n" unless $main::quiet;
      print STDERR "$blast_finished\n" unless $main::quiet;
      return $blast_dir;
   }
	
   #call blast executable
   $chunk->write_file($t_file_name);  

   runtBlastx($t_file_name,
	      $tmp_db,
	      $o_file,
	      $tblastx,
	      $eval_tblastx,
	      $split_hit,
	      $cpus,
	      $opt_f
	     );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_tblastx{
   my $chunk      = shift;
   my $blast_dir    = shift;
   my $eval_tblastx = shift;
   my $bit_tblastx = shift,
   my $pcov_tblastx = shift;
   my $pid_tblastx = shift;
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
   $params{significance}  = $eval_tblastx;
   $params{hsp_bit_min}   = $bit_tblastx;
   $params{percov}        = $pcov_tblastx;
   $params{percid}        = $pid_tblastx;
   $params{split_hit}     = $split_hit;

   my $chunk_keepers = Widget::tblastx::parse($blast_finished,
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
sub tblastx {
   my $chunk      = shift;
   my $db         = shift;
   my $the_void    = shift;
   my $seq_id     = shift;
   my $tblastx = shift;
   my $eval_tblastx = shift;
   my $bit_tblastx = shift,
   my $pcov_tblastx = shift;
   my $pid_tblastx = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $opt_f = shift;
   my $LOG = shift;

   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n  =~ s/\.fasta$//;
	
   my $chunk_number = $chunk->number();
   my $q_length = $chunk->parent_seq_length();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.tblastx";

   $LOG->add_entry("STARTED", $o_file, ""); 

   $chunk->write_file($file_name);
   runtBlastx($file_name,
	      $db,
	      $o_file,
	      $tblastx,
	      $eval_tblastx,
	      $split_hit,
	      $cpus,
	      $opt_f
	     );

   my %params;
   $params{significance}  = $eval_tblastx;
   $params{hsp_bit_min}   = $bit_tblastx;
   $params{percov}        = $pcov_tblastx;
   $params{percid}        = $pid_tblastx;
   $params{split_hit}     = $split_hit;

   my $chunk_keepers = Widget::tblastx::parse($o_file,
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
sub runtBlastx {
   my $q_file   = shift;
   my $db       = shift;
   my $out_file = shift;
   my $blast = shift;
   my $eval_blast = shift;
   my $split_hit = shift;
   my $cpus = shift;
   my $opt_f = shift;

   my $command  = $blast;
   if ($command =~ /tblastx$/) {
      $command .= " $db $q_file B=100000 V=100000 E=$eval_blast";
      $command .= " wordmask=seg";
      #$command .= " W=15";
      $command .= " Z=1000";
      $command .= " Y=500000000";
      $command .= " cpus=$cpus";	
      $command .= " topcomboN=1";
      $command .= " hspmax=100";
      $command .= " gspmax=100";
      $command .= " hspsepqmax=$split_hit";
      $command .= " lcmask";
      $command .= " wordmask=seg";
      $command .= " gi";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p tblastx";
      $command .= " -d $db -i $q_file -b 100000 -v 100000 -e $eval_blast";
      #$command .= " -W 15";
      $command .= " -z 1000";
      $command .= " -Y 500000000";
      $command .= " -a $cpus";	
      $command .= " -K 100";
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   else{
      die "ERROR: Must be a tblastx executable";  
   }

   my $w = new Widget::tblastx();
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
   $file_name .= ".specific" if($rmlib);
   my $o_file    = "$file_name\.out";
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
   }
   else {
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
#returns a directory name to write analysis output to
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
#this function sets the defualt values for control options
#the function also determines what control options are valid, so
#to add control options edit this function
sub set_defaults {
   my $type = shift || 'all';

   if ($type !~ /^all$|^opts$|^bopt$|^exe$/) {
      warn "WARNING: Invalid type \'$type\' in S_Func ::set_defaults";
      $type = 'all';
   }

   my %CTL_OPT;

   #maker_opts
   if ($type eq 'all' || $type eq 'opts') {
      $CTL_OPT{'genome'} = '';
      $CTL_OPT{'genome_gff'} = '';
      $CTL_OPT{'est_pass'} = 0;
      $CTL_OPT{'altest_pass'} = 0;
      $CTL_OPT{'protein_pass'} = 0;
      $CTL_OPT{'rm_pass'} = 0;
      $CTL_OPT{'model_pass'} = 1;
      $CTL_OPT{'pred_pass'} = 0;
      $CTL_OPT{'other_pass'} = 1;
      $CTL_OPT{'est'} = '';
      $CTL_OPT{'est_reads'} = '';
      $CTL_OPT{'altest'} = '';
      $CTL_OPT{'est_gff'} = '';
      $CTL_OPT{'altest_gff'} = '';
      $CTL_OPT{'protein'} = '';
      $CTL_OPT{'protein_gff'} = '';
      $CTL_OPT{'model_org'} = 'all';
      $CTL_OPT{'repeat_protein'} = Cwd::abs_path("$FindBin::Bin/../data/te_proteins.fasta");
      $CTL_OPT{'rmlib'} = '';
      $CTL_OPT{'rm_gff'} = '';
      $CTL_OPT{'predictor'} = 'est2genome';
      $CTL_OPT{'snaphmm'} = 'fly';
      $CTL_OPT{'augustus_species'} = 'fly';
      $CTL_OPT{'fgenesh_par_file'} = 'Dicots';
      $CTL_OPT{'model_gff'} = '';
      $CTL_OPT{'pred_gff'} = '';
      $CTL_OPT{'other_gff'} = '';
      $CTL_OPT{'alt_peptide'} = 'c';
      $CTL_OPT{'cpus'} = 1;
      $CTL_OPT{'max_dna_len'} = 100000;
      $CTL_OPT{'min_contig'} = 1;
      $CTL_OPT{'split_hit'} = 10000;
      $CTL_OPT{'pred_flank'} = 200;
      $CTL_OPT{'single_exon'} = 0;
      $CTL_OPT{'keep_preds'} = 0;
      $CTL_OPT{'retry'} = 1;
      $CTL_OPT{'datastore'} = 0;
      $CTL_OPT{'clean_up'} = 0;
   }

   #maker_bopts
   if ($type eq 'all' || $type eq 'bopts') {
      $CTL_OPT{'blast_type'} = 'wublast';
      $CTL_OPT{'pcov_blastn'} = 0.80;
      $CTL_OPT{'pid_blastn'} = 0.85;
      $CTL_OPT{'eval_blastn'} = 1e-10;
      $CTL_OPT{'bit_blastn'} = 40;
      $CTL_OPT{'pcov_blastx'} = 0.50;
      $CTL_OPT{'pid_blastx'} = 0.40;
      $CTL_OPT{'eval_blastx'} = 1e-6;
      $CTL_OPT{'bit_blastx'} = 30;
      $CTL_OPT{'pcov_rm_blastx'} = 0.50;
      $CTL_OPT{'pid_rm_blastx'} = 0.40;
      $CTL_OPT{'eval_rm_blastx'} = 1e-6;
      $CTL_OPT{'bit_rm_blastx'} = 30;
      $CTL_OPT{'pcov_tblastx'} = 0.80;
      $CTL_OPT{'pid_tblastx'} = 0.85;
      $CTL_OPT{'eval_tblastx'} = 1e-10;
      $CTL_OPT{'bit_tblastx'} = 40;
      $CTL_OPT{'en_score_limit'} = 20;
      $CTL_OPT{'ep_score_limit'} = 20;
   }

   #maker_exe
   if ($type eq 'all' || $type eq 'exe') {
      my @exes = ('xdformat',
		  'formatdb',
		  'blastall',
		  'blastn',
		  'blastx',
		  'tblastx',
		  'RepeatMasker',
		  'exonerate',
		  'snap',
		  'augustus',
		  'fgenesh',
		  'twinscan',
		  'jigsaw',
		  'qrna',
		  'fathom',
		 );

      foreach my $exe (@exes) {
	 my $loc = `which $exe 2> /dev/null`;
	 chomp $loc;
	 if ($loc =~ /^no $exe/) {
	    $CTL_OPT{$exe} = '';
	 }
	 else {
	    $CTL_OPT{$exe} = $loc;
	 }
      }
   }

   #evaluator
   if ($type eq 'all' || $type eq 'eva') {
      $CTL_OPT{'eva_pcov_blastn'} = 0.80;
      $CTL_OPT{'eva_pid_blastn'} = 0.85;
      $CTL_OPT{'eva_eval_blastn'} = 1e-10;
      $CTL_OPT{'eva_bit_blastn'} = 40;
      $CTL_OPT{'side_thre'} = 5;
      $CTL_OPT{'eva_window_size'} = 70;
      $CTL_OPT{'eva_split_hit'} = 1;
      $CTL_OPT{'eva_hspmax'} = 100;
      $CTL_OPT{'eva_gspmax'} = 100;
      $CTL_OPT{'enable_fathom'} = 1;

   }
   
   return %CTL_OPT;
}
#-----------------------------------------------------------------------------
#this function parses the control files and sets up options for each maker run
#error checking for starup occurs here
sub load_control_files {
   my @ctlfiles = @{shift @_};
   my %OPT = %{shift @_};
   my $mpi_size = shift @_ || 1; 

   #--set default values and control structure
   my %CTL_OPT = set_defaults();

   my $error;	      #hold all fatal errors from control file parsing

   #--load values from control files
   foreach my $ctlfile (@ctlfiles) {
      open (CTL, "< $ctlfile") or die"ERROR: Could not open control file \"$ctlfile\".\n";
	
      while (my $line = <CTL>) {
	 chomp($line);
	    
	 if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:]+)\:([^\n\#]*)/) {
	    my $key = $1;
	    my $value = $2;

	    #remove preceding and trailing whitespace
	    $value =~ s/^[\s\t]+|[\s\t]+$//g;
	    
	    if (exists $CTL_OPT{$key}) { #should already exist or is a bad value
	       #resolve environmental variables
	       if ($value =~ /\$/) {
		  $value = `echo \"$value\"`;
		  chomp($value);
	       }

	       #require numerical values for certain options
	       if ($CTL_OPT{$key} =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/ &&
		   $value !~  /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/
		  ) {
		  $error .= "ERROR: valid number required for option \'$key\'\n\n"
	       }

	       #set value
	       $CTL_OPT{$key} = defined($value) ? $value : '';
	    }
	    else {
	       warn "WARNING: Invalid option \'$key\' in control file $ctlfile\n\n";
	    }
	 }
      }
   }

   #--load command line options
   $CTL_OPT{genome} = $OPT{genome} if (defined $OPT{genome});
   $CTL_OPT{force} = $OPT{force} if (defined $OPT{force});
   $CTL_OPT{predictor} = $OPT{predictor} if (defined $OPT{predictor});
   $CTL_OPT{retry} = $OPT{retry} if (defined $OPT{retry});
   $CTL_OPT{cpus} = $OPT{cpus} if (defined $OPT{cpus});

   #skip repeat masking command line option
   if ($OPT{R}) {
      $CTL_OPT{model_org} = '';
      $CTL_OPT{repeat_protein} = '';
      $CTL_OPT{rmlib} = '';
      $CTL_OPT{rm_gff} = '';
      $CTL_OPT{rm_pass} = 0;
   }

   #parse predictor and error check
   $CTL_OPT{predictor} =~ s/\s+//g;
   my @predictors = split(',', $CTL_OPT{predictor});
   $CTL_OPT{_predictor} = \ @predictors;

   foreach my $p (@predictors) {
      if ($p !~ /^snap$|^augustus$|^est2genome$|^fgenesh$/ &&
	  $p !~ /^twinscan$|^jigsaw$|^gff$|^abinit$/
	 ) {
	 $error .= "ERROR: Invalid predictor defined: $p\n".
	 "Valid entries are: est2genome, abinit, gff, snap, augustus,\n".
	 "or fgenesh\n\n";
      }
   }

   #check blast type validity and related values (NCBI vs. WUBLAST)
   if ($CTL_OPT{blast_type} !~ /^wublast$|^ncbi$/) {
      warn "WARNING: blast_type must be set to \'wublast\' or \'ncbi\'.\n",
      "The value $CTL_OPT{blast_type} is invalid.\n",
      "This will now be reset to the default 'wublast'";
      
      $CTL_OPT{blast_type} = 'wublast';
   }
   
   if ($CTL_OPT{blast_type} =~ /^wublast$/ &&
       ! -e $CTL_OPT{blastn} &&
       ! -e $CTL_OPT{blastx} &&
       ! -e $CTL_OPT{tblastx} &&
       -e $CTL_OPT{blastall}
      ) {
      warn "WARNING: blast_type is set to \'wublast\' but wublast executables\n",
      "can not be located.  NCBI blast will be used instead.\n";

      $CTL_OPT{blast_type} = 'ncbi';
   }

   if ($CTL_OPT{blast_type} =~ /^ncbi$/ &&
       ! -e $CTL_OPT{blastall} &&
       -e $CTL_OPT{blastn} &&
       -e $CTL_OPT{blastx} &&
       -e $CTL_OPT{tblastx}
      ) {
      warn "WARNING: blast_type is set to \'ncbi\' but ncbi executables\n",
      "can not be located.  WUBLAST blast will be used instead.\n";

      $CTL_OPT{blast_type} = 'wublast';
   }
   
   if ($CTL_OPT{blast_type} =~ /^wublast$/) {
      $CTL_OPT{_formater} = $CTL_OPT{xdformat};
      $CTL_OPT{_blastn} = $CTL_OPT{blastn};
      $CTL_OPT{_blastx} = $CTL_OPT{blastx};
      $CTL_OPT{_tblastx} = $CTL_OPT{tblastx};
   }
   elsif ($CTL_OPT{blast_type} =~ /^ncbi$/) {
      $CTL_OPT{_formater} = $CTL_OPT{formatdb};
      $CTL_OPT{_blastn} = $CTL_OPT{blastall};
      $CTL_OPT{_blastx} = $CTL_OPT{blastall};
      $CTL_OPT{_tblastx} = $CTL_OPT{blastall};
   }
   
   #--validate existence of required values from control files

   #required files in here
   my @infiles;

   #decide if required
   push (@infiles, '_blastn', '_formater') if($CTL_OPT{est});
   push (@infiles, '_blastn', '_formater') if($CTL_OPT{est_reads}); 
   push (@infiles, '_blastx', '_formater') if($CTL_OPT{protein}); 
   push (@infiles, '_blastx', '_formater') if($CTL_OPT{repeat_protein}); 
   push (@infiles, '_tblastx', '_formater') if($CTL_OPT{altest});
   push (@infiles, 'genome') if($CTL_OPT{genome});
   push (@infiles, 'genome') if(!$CTL_OPT{genome_gff});
   push (@infiles, 'exonerate') if($CTL_OPT{est}); 
   push (@infiles, 'exonerate') if($CTL_OPT{protein}); 
   push (@infiles, 'repeat_protein') if ($CTL_OPT{repeat_protein});
   push (@infiles, 'est') if($CTL_OPT{est}); 
   push (@infiles, 'protein') if($CTL_OPT{protein}); 
   push (@infiles, 'altest') if($CTL_OPT{altest}); 
   push (@infiles, 'est_reads') if($CTL_OPT{est_reads}); 
   push (@infiles, 'RepeatMasker') if($CTL_OPT{rmlib});
   push (@infiles, 'RepeatMasker') if($CTL_OPT{model_org});
   push (@infiles, 'rmlib') if ($CTL_OPT{rmlib});
   push (@infiles, 'snap') if (grep (/snap/, $CTL_OPT{predictor}));
   push (@infiles, 'snap') if ($CTL_OPT{snap});
   push (@infiles, 'augustus') if (grep (/augustus/, $CTL_OPT{predictor})); 
   push (@infiles, 'augustus') if ($CTL_OPT{augustus});
   push (@infiles, 'fgenesh') if (grep (/fgenesh/, $CTL_OPT{predictor}));
   push (@infiles, 'fgenesh') if ($CTL_OPT{fgenesh});
   push (@infiles, 'twinscan') if (grep (/twinscan/, $CTL_OPT{predictor}));
   push (@infiles, 'twinscan') if ($CTL_OPT{twinscan});
   push (@infiles, 'jigsaw') if (grep (/jigsaw/, $CTL_OPT{predictor}));
   push (@infiles, 'jigsaw') if ($CTL_OPT{jigsaw});
   push (@infiles, 'qrna') if ($CTL_OPT{qrna});
   push (@infiles, 'rm_gff') if($CTL_OPT{rm_gff});
   push (@infiles, 'est_gff') if($CTL_OPT{est_gff});
   push (@infiles, 'protein_gff') if($CTL_OPT{protein_gff});
   push (@infiles, 'genome_gff') if($CTL_OPT{genome_gff});
   push (@infiles, 'pred_gff') if($CTL_OPT{pred_gff});
   push (@infiles, 'model_gff') if ($CTL_OPT{model_gff});
   push (@infiles, 'model_gff') if (grep (/gff/, $CTL_OPT{predictor}) &&
				    (!$CTL_OPT{genome_gff} ||
				     !$CTL_OPT{model_pass})
				   );

   #uniq the array
   my %uniq;
   foreach my $in (@infiles) {
      $uniq{$in}++;
   }
   @infiles = keys %uniq;

   #verify existence of required values
   foreach my $in (@infiles) {
      if (not $CTL_OPT{$in}) {
	 $error .= "ERROR: You have failed to provide a value for \'$in\' in the control files\n\n";
	 next;
      }

      if (not -e $CTL_OPT{$in}) {
	 $error .= "ERROR: The \'$in\' file $CTL_OPT{$in} does not exist.\n".
	 "Please check your control files: maker_opts.ctl, maker_bopts, or maker_exe.ctl\n\n";
	 next;
      }
      
      #set the absolute path to the file to reduce ambiguity
      #this breaks blast which requires symbolic links
      $CTL_OPT{$in} = Cwd::abs_path($CTL_OPT{$in}) unless ($in =~ /^_blastn$|^_blastx$|^_tblastx$/);
   }
			       
   #--error check that values are meaningful
   if ($CTL_OPT{augustus} && not $CTL_OPT{augustus_species}) {
      warn "WARNING: There is no species specified for Augustus in maker_opts.ctl augustus_species.\n".
      "As a result the default (fly) will be used.\n\n";
      $CTL_OPT{augustus_species} = "fly";
   }
   if ($CTL_OPT{augustus} &&
       (! $ENV{AUGUSTUS_CONFIG_PATH} || ! -e "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg")
      ) {
      $error .= "ERROR: The environmental variable AUGUSTUS_CONFIG_PATH has not been set\n".
      "or is not set correctly Please set this in your profile per Augustus\n".
      "installation instructions\n\n";
   }
   if (($CTL_OPT{snap}||$CTL_OPT{enable_fathom}) && not $CTL_OPT{snaphmm}) {
      warn "WARNING: There is no model specified for for Snap in maker_opts.ctl snaphmm.\n".
      "As a result, the default (fly) will be used.\n\n";
      $CTL_OPT{snaphmm} = "fly";
   }

   if (($CTL_OPT{snap} || $CTL_OPT{enable_fathom}) &&
       ! -e $CTL_OPT{snaphmm} &&
       (! exists $ENV{ZOE} || ! -e $ENV{ZOE}."/HMM/".$CTL_OPT{snaphmm})
      ) {
      $error .= "ERROR: The snaphmm specified for Snap in maker_opts.ctl does not exist.\n\n";
   }
   if ($CTL_OPT{fgenesh}) {
      if (! $CTL_OPT{fgenesh_par_file}) {
	 $error .= "ERROR: There is no parameter file secified for for FgenesH in\n".
	 "maker_opts.ctl fgenesh_par_file\n\n";
      }
      elsif (! -e $CTL_OPT{fgenesh_par_file}) {
	 $error .= "ERROR: The parameter file specified for fgenesh in maker_opts.ctl does not exist.\n\n";
      }
   }
   if ($CTL_OPT{max_dna_len} < 50000) {
      warn "WARNING: \'max_dna_len\' is set too low.  The minimum value permited is 50,000\n".
      "max_dna_len will be reset to 50,000\n\n";
      $CTL_OPT{max_dna_len} = 50000;
   }
   if ($CTL_OPT{min_contig} <= 0) {
      warn "WARNING: \'min_contig\' must be set to 1 or higher.\n".
      "min_contig will be reset to 1\n\n";
      $CTL_OPT{min_contig} = 1;
   }
   if ($CTL_OPT{retry} < 0) {
      warn "WARNING: \'retry\' must be set to 0 or greater.\n\n";
      $CTL_OPT{retry} = 0;
   }

   die $error if ($error);   

   #--check genome fasta file
   my $iterator = new Iterator::Any( -fasta => $CTL_OPT{genome},
				     -gff => $CTL_OPT{genome_gff}
				   );

   if ($iterator->number_of_entries() == 0) {
      my $genome = (! $CTL_OPT{genome}) ? $CTL_OPT{genome_gff} : $CTL_OPT{genome};
      die "ERROR:  The file $genome contains no fasta entries\n";
   }

   #--decide whether to force datastore 
   if ($iterator->number_of_entries() > 1000 && ! $CTL_OPT{datastore}) {
      warn "WARNING:  There are more than 1000 fasta entries in the input file.\n".
      "Datastore will be used to avoid overloading the data structure of\n".
      "the output directory.\n\n";
      
      $CTL_OPT{datastore} = 1;
   }

   #--decide if gff database should be created
   my @gffs = grep(/\_gff$/, @{[keys %CTL_OPT]});
   foreach my $key (@gffs) {
      if ($CTL_OPT{$key}) {
	 $CTL_OPT{go_gffdb} = 1;
	 last;
      }
   }

   #--check validity of the alternate peptide
   $CTL_OPT{alt_peptide} = uc($CTL_OPT{alt_peptide});
   if ($CTL_OPT{alt_peptide} !~ /^[ABCDEFGHIKLMNPQRSTVWXYZ]$/) {
      warn "WARNING: Invalid alternate peptide \'$CTL_OPT{alt_peptide}\'\n",
      "This will be set to the default 'C'\n";
      $CTL_OPT{alt_peptide} = 'C';
   }

   #---set up blast databases for analyisis
   create_blastdb(\%CTL_OPT, $mpi_size);

   #--set values for datastructure
   my $genome = $CTL_OPT{genome};
   $genome = $CTL_OPT{genome_gff} if (not $genome);
   ($CTL_OPT{out_name}) = $genome =~ /([^\/]+)$/;
   $CTL_OPT{out_name} =~ s/\.[^\.]+$//;
   if(! $CTL_OPT{out_base}){
      $CTL_OPT{out_base} = Cwd::cwd()."/$CTL_OPT{out_name}.maker.output";
   }
   mkdir($CTL_OPT{out_base}) if(! -d $CTL_OPT{out_base});
   die "ERROR: Could not build output directory $CTL_OPT{out_base}\n"
        if(! -d $CTL_OPT{out_base});

   #--log the control files
   generate_control_files($CTL_OPT{out_base}, %CTL_OPT);

   return %CTL_OPT;
}
#-----------------------------------------------------------------------------
#this function generates generic control files
sub generate_control_files {
   my $dir = shift || Cwd::cwd();
   my %O = (@_) ? @_ : set_defaults();

   #--build maker_opts.ctl file
   open (OUT, "> $dir/maker_opts.ctl");
   print OUT "#-----Genome (Required for De-Novo Annotations)\n";
   print OUT "genome:$O{genome} #genome sequence file in fasta format\n";
   print OUT "\n";
   print OUT "#-----Re-annotation Options\n";
   print OUT "genome_gff:$O{genome_gff} #re-annotate genome based on this gff3 file\n";
   print OUT "est_pass:$O{est_pass} #use ests in genome_gff: 1 = yes, 0 = no\n";
   print OUT "altest_pass:$O{altest_pass} #use alternate organism ests in genome_gff: 1 = yes, 0 = no\n";
   print OUT "protein_pass:$O{protein_pass} #use proteins in genome_gff: 1 = yes, 0 = no\n";
   print OUT "rm_pass:$O{protein_pass} #use repeats in genome_gff: 1 = yes, 0 = no\n";
   print OUT "model_pass:$O{protein_pass} #use gene models in genome_gff: 1 = yes, 0 = no\n";
   print OUT "pred_pass:$O{pred_pass} #use ab-initio predictions in genome_gff: 1 = yes, 0 = no\n";
   print OUT "other_pass:$O{protein_pass} #passthrough everything else in genome_gff: 1 = yes, 0 = no\n";
   print OUT "\n";
   print OUT "#-----EST Evidence (you must provide a value for at least one)\n";
   print OUT "est:$O{est} #non-redundant set of assembled ESTs in fasta format (classic EST analysis)\n";
   print OUT "est_reads:$O{est_reads} #un-assembled EST reads in fasta format (for deep nextgen mRNASeq)\n";
   print OUT "altest:$O{altest} #EST/cDNA sequence file in fasta format from an alternate organism\n";
   print OUT "est_gff:$O{est_gff} #EST evidence from a seperate gff3 file\n";
   print OUT "altest_gff:$O{altest_gff} #Alternate organism EST evidence from a seperate gff3 file\n";
   print OUT "\n";
   print OUT "#-----Protein Homology Evidence (you must provide a value for at least one)\n";
   print OUT "protein:$O{protein}  #protein sequence file in fasta format\n";
   print OUT "protein_gff:$O{protein_gff}  #protein homology evidence from a gff3 file\n";
   print OUT "\n";
   print OUT "#-----Repeat Masking (leave values blank to skip)\n";
   print OUT "model_org:$O{model_org} #model organism for RepBase masking in RepeatMasker\n";
   print OUT "repeat_protein:$O{repeat_protein} #a database of transposable element proteins in fasta format\n";
   print OUT "rmlib:$O{rmlib} #an organism specific repeat library in fasta format\n";
   print OUT "rm_gff:$O{rm_gff} #repeat elements from a gff3 file\n";
   print OUT "\n";
   print OUT "#-----Gene Prediction Options\n";
   print OUT "predictor:$O{predictor} #prediction method for annotations (seperate multiple values by ',')\n";
   print OUT "snaphmm:$O{snaphmm} #SNAP HMM model\n";
   print OUT "augustus_species:$O{augustus_species} #Augustus gene prediction model\n";
   print OUT "fgenesh_par_file:$O{fgenesh_par_file} #Fgenesh parameter file\n";
   print OUT "model_gff:$O{model_gff} #gene models from a gff3 file (annotation passthrough)\n";
   print OUT "pred_gff:$O{model_gff} #ab-initio predictions from a gff3 file\n";
   print OUT "\n";
   print OUT "#-----Other Annotation Type Options (features maker doesn't recognize)\n";
   print OUT "other_gff:$O{other_gff} #features to passthrough to final output from a gff3 file\n";
   print OUT "\n";
   print OUT "#-----External Application Specific Options\n";
   print OUT "alt_peptide:$O{alt_peptide} #amino acid used to replace non standard amino acids in blast databases\n";
   print OUT "cpus:$O{cpus} #max number of cpus to use in BLAST and RepeatMasker\n";
   print OUT "\n";
   print OUT "#-----Maker Specific Options\n";
   print OUT "max_dna_len:$O{max_dna_len} #length for dividing up contigs into chunks (larger values increase memory usage)\n";
   print OUT "min_contig:$O{min_contig} #all contigs from the input genome file below this size will be skipped\n";
   print OUT "split_hit:$O{split_hit} #length for the splitting of hits (expected max intron size for EST and protein alignments)\n";
   print OUT "pred_flank:$O{pred_flank} #length of sequence surrounding EST and protein evidence used to extend gene predictions\n";
   print OUT "single_exon:$O{single_exon} #consider EST hits aligning to single exons when generating annotations, 1 = yes, 0 = no\n";
   print OUT "keep_preds:$O{keep_preds} #keep ab-initio predictions that do not overlap a maker annotation, 1 = yes, 0 = no\n";
   print OUT "retry:$O{retry} #number of times to retry a contig if annotation fails for some reason\n";
   print OUT "clean_up:$O{clean_up} #remove theVoid directory: 1 = yes, 0 = no\n";
   close (OUT);
    
   #--build maker_bopts.ctl file
   open (OUT, "> $dir/maker_bopts.ctl");
   print OUT "#-----BLAST and Exonerate Statistics Thresholds\n";
   print OUT "blast_type:$O{blast_type} #set to 'wublast' or 'ncbi'\n";
   print OUT "\n";
   print OUT "pcov_blastn:$O{pcov_blastn} #Blastn Percent Coverage Threhold EST-Genome Alignments\n";
   print OUT "pid_blastn:$O{pid_blastn} #Blastn Percent Identity Threshold EST-Genome Aligments\n";
   print OUT "eval_blastn:$O{eval_blastn} #Blastn eval cutoff\n";
   print OUT "bit_blastn:$O{bit_blastn} #Blastn bit cutoff\n";
   print OUT "\n";
   print OUT "pcov_blastx:$O{pcov_blastx} #Blastx Percent Coverage Threhold Protein-Genome Alignments\n";
   print OUT "pid_blastx:$O{pid_blastx} #Blastx Percent Identity Threshold Protein-Genome Aligments\n";
   print OUT "eval_blastx:$O{eval_blastx} #Blastx eval cutoff\n";
   print OUT "bit_blastx:$O{bit_blastx} #Blastx bit cutoff\n";
   print OUT "\n";
   print OUT "pcov_rm_blastx:$O{pcov_rm_blastx} #Blastx Percent Coverage Threhold For Transposable Element Masking\n";
   print OUT "pid_rm_blastx:$O{pid_rm_blastx} #Blastx Percent Identity Threshold For Transposbale Element Masking\n";
   print OUT "eval_rm_blastx:$O{eval_rm_blastx} #Blastx eval cutoff for transposable element masking\n";
   print OUT "bit_rm_blastx:$O{bit_rm_blastx} #Blastx bit cutoff for transposable element masking\n";
   print OUT "\n";
   print OUT "pcov_tblastx:$O{pcov_tblastx} #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments\n";
   print OUT "pid_tblastx:$O{pid_tblastx} #tBlastx Percent Identity Threshold alt-EST-Genome Aligments\n";
   print OUT "eval_tblastx:$O{eval_tblastx} #tBlastx eval cutoff\n";
   print OUT "bit_tblastx:$O{bit_tblastx} #tBlastx bit cutoff\n";
   print OUT "\n";
   print OUT "ep_score_limit:$O{ep_score_limit} #exonerate protein percent of maximal score threshold\n";
   print OUT "en_score_limit:$O{en_score_limit} #exonerate nucleotide percent of maximal score threshold\n";
   close(OUT);
    
   #--build maker_exe.ctl file
   open (OUT, "> $dir/maker_exe.ctl");
   print OUT "#-----Location of Executables Used by Maker\n";
   print OUT "formatdb:$O{formatdb} #location of NCBI formatdb executable\n";
   print OUT "blastall:$O{blastall} #location of NCBI blastall executable\n";
   print OUT "xdformat:$O{xdformat} #location of WUBLAST xdformat executable\n";
   print OUT "blastn:$O{blastn} #location of WUBLAST blastn executable\n";
   print OUT "blastx:$O{blastx} #location of WUBLAST blastx executable\n";
   print OUT "tblastx:$O{tblastx} #location of WUBLAST tblastx executable\n";
   print OUT "RepeatMasker:$O{RepeatMasker} #location of RepeatMasker executable\n";
   print OUT "exonerate:$O{exonerate} #location of exonerate executable\n";
   print OUT "\n";
   print OUT "#-----Ab-initio Gene Prediction Algorithms\n";
   print OUT "snap:$O{snap} #location of snap executable\n";
   print OUT "augustus:$O{augustus} #location of augustus executable\n";
   print OUT "fgenesh:$O{fgenesh} #location of fgenesh executable\n";
#   print OUT "twinscan:$O{twinscan} #location of twinscan executable\n";
#   print OUT "fathom:$O{fathom} #location of fathom executable\n";
   print OUT "\n";
   print OUT "#-----Other Algorithms\n";
   print OUT "jigsaw:$O{jigsaw} #location of jigsaw executable (not yet implemented)\n";
   print OUT "qrna:$O{qrna} #location of qrna executable (not yet implemented)\n";
   close(OUT);
    
   #--build evaluator.ctl file
   open (OUT, "> $dir/evaluator.ctl");
   print OUT "#-----EVALUATOR Control Options\n";
   print OUT "eva_pcov_blastn:$O{eva_pcov_blastn} #Blastn Percent Coverage Threshold EST-Genome Alignments\n";
   print OUT "eva_pid_blastn:$O{eva_pid_blastn} #Blastn Percent Identity Threshold EST-Genome Alignments\n";
   print OUT "eva_eval_blastn:$O{eva_eval_blastn} #Blastn eval cutoff\n";
   print OUT "eva_bit_blastn:$O{eva_bit_blastn} #Blastn bit cutoff\n";
   print OUT "side_thre:$O{side_thre}\n";
   print OUT "eva_window_size:$O{eva_window_size}\n";
   print OUT "eva_split_hit:$O{eva_split_hit}\n";
   print OUT "eva_hspmax:$O{eva_hspmax}\n";
   print OUT "eva_gspmax:$O{eva_gspmax}\n";
#   print OUT "enable_fathom:$O{enable_fathom}\n";
   close (OUT);
}

1;
