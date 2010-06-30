#------------------------------------------------------------------------
#----                                GI                              ----
#------------------------------------------------------------------------
package GI;

use strict;
use vars qw(@ISA @EXPORT $VERSION $TMP $LOCK);
use FindBin;
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
use Widget::genemark; 
use Widget::fgenesh;
use Widget::xdformat;
use Widget::formatdb;
use PhatHit_utils;
use Shadower;
use polisher::exonerate::protein;
use polisher::exonerate::est;
use maker::auto_annotator;
use cluster;
use repeat_mask_seq;
use maker::sens_spec;
use File::NFSLock;
use FastaDB;
use Digest::MD5;

@ISA = qw(
	);

$TMP = tempdir("maker_XXXXXX", CLEANUP => 1, TMPDIR => 1);
print STDERR "TMP_STAT: TMP is being initialized to $TMP: PID=$$\n" if($main::dtmp);
#------------------------------------------------------------------------
#--------------------------- CLASS FUNCTIONS ----------------------------
#------------------------------------------------------------------------
sub set_global_temp {
    my $dir = shift;

    return if(! $dir);

    #remove old tempdir if user supplied a new one
    if($TMP ne $dir){
	print STDERR "\nTMP_STAT: Trying to change TMP from $TMP to $dir: PID=$$\n" if($main::dtmp); ##temp
	my $base = $dir;
	$base =~ s/[^\/]+$//;

	if(! -d $base){
	    print STDERR "TMP_STAT: base directory $base does not exist, keeping TMP as $TMP: PID=$$\n" if($main::dtmp); ##temp
	    return;
	}

	if(! -d $dir){
	    print STDERR "TMP_STAT: base directory $base exists but directory $dir does not, trying to create: PID=$$\n" if($main::dtmp); ##temp
	    mkdir($dir);
	}

	if(! -d $dir){
	    print STDERR "TMP_STAT: directory $dir does not exist, keeping TMP as $TMP: PID=$$\n\n" if($main::dtmp); ##temp
	    return;
	}

	File::Path::rmtree($TMP);

	$TMP = $dir;
	print STDERR "TMP_STAT: Success TMP is now $dir: PID=$$\n\n" if($main::dtmp); ##temp
    }
}
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
   my $fasta_index = shift @_;
   my $blast_keepers = shift @_;
   my $blast_holdovers = shift @_;
   my $the_void = shift @_;
   my %CTL_OPT = %{shift @_};
   my $type = shift @_; #blastn, blastx, or tblastx
   my $LOG = shift @_;

   PhatHit_utils::merge_hits($blast_keepers,
			     $blast_holdovers, 
			     $CTL_OPT{split_hit},
			    );
   @{$blast_holdovers} = ();

   $blast_keepers = reblast_merged_hits($fasta,
					$blast_keepers,
					$fasta_index,
					$the_void,
					$type,
					\%CTL_OPT,
					$LOG
					);

   return $blast_keepers;
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
      my $fastaObj = $db_index->get_Seq_for_hit($hit);
      
      #still no sequence? try rebuilding the index and try again
      if (not $fastaObj) {
	  print STDERR "WARNING: Cannot find> ".$hit->name.", trying to re-index the fasta.\n";
	  $db_index->reindex();
	  $fastaObj = $db_index->get_Seq_for_hit($hit);
	  if (not $fastaObj) {
	      print STDERR "stop here:".$hit->name."\n";
	      die "ERROR: Fasta index error\n";
	  }
      }
      
      #get fasta def and seq
      my $t_seq      = $fastaObj->seq();
      my $t_def      = $db_index->header_for_hit($hit);
      
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
			      \%CTL_OPT,
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
			      \%CTL_OPT,
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
			       \%CTL_OPT,
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
    my $chunk      = shift @_;
    my $split_hit  = shift @_;
    my $pred_flank = shift @_;
    my $s_flag     = shift @_; #indicates whether to treat strands independantly 
    my $groups_cfh = shift @_; #group to cluster and find holdovers
    
    my $phat_hits;
    
    foreach my $group (@{$groups_cfh}) {
	push(@{$phat_hits}, @{$group});
    }
    
    my $p_hits = [];
    my $m_hits = [];
    
    #seperate by strand or not (this makes chunk cutoffs strand independant)
    if($s_flag){
	($p_hits, $m_hits) = PhatHit_utils::seperate_by_strand('query', $phat_hits, 1); #exonerate flag set
    }
    else{
	$p_hits = $phat_hits;
    }
    
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
    
    my $p_pieces = Shadower::getPieces(\$chunk->seq, $p_coors, $pred_flank);
    $p_pieces = [sort {$b->{e} <=> $a->{e}} @{$p_pieces}];
    my $m_pieces = Shadower::getPieces(\$chunk->seq, $m_coors, $pred_flank);
    $m_pieces = [sort {$b->{e} <=> $a->{e}} @{$m_pieces}];
    
    my $cutoff = $chunk->length + $chunk->offset - $split_hit;
    my $p_cutoff = $chunk->length + $chunk->offset + 1;
    my $m_cutoff = $chunk->length + $chunk->offset + 1;

    my @keepers;
    my @holdovers;

    #no internal cutoff if this is the last contig
    $cutoff = $chunk->length + $chunk->offset + 1 if($chunk->is_last);
    
    #adjust cutoff to overlapping hits
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
    
    #too small, all are heldover for next round
    if ($p_cutoff <= 1 + $chunk->offset &&
	$m_cutoff <= 1 + $chunk->offset) {
	foreach my $g (@{$groups_cfh}){
	    push (@holdovers, $g);
	    push (@keepers, []);
	}
	return @keepers, @holdovers;
    }
    
    #seperate holdovers and keepers
    foreach my $group (@{$groups_cfh}) {
	my $group_keepers = [];
	my $group_holdovers = [];
	
	foreach my $hit (@{$group}) {
	    my $b = $hit->nB('query');
	    my $e = $hit->nE('query');
	    my $strand = $hit->strand;
	    
	    #exonerate counterpart check (blastn with flipped exonerate)
	    $strand *= -1 if ($hit->{_exonerate_flipped});
	    
	    #if stands are not being treated independantly, treat all as plus strand
	    $strand = 1 if (!$s_flag);
	    
	    ($b, $e) = ($e, $b) if $b > $e;
	    
	    if (($strand eq '1' && $e < $p_cutoff && $p_cutoff > $chunk->offset +1) ||
		($strand eq '-1' && $e < $m_cutoff && $m_cutoff > $chunk->offset +1)
		) {
		$hit->{_holdover} = 0;
		push(@{$group_keepers}, $hit);
	    }
	    else {
		$hit->{_holdover} = 1;
		push(@{$group_holdovers}, $hit);
	    }
	}
	
	push(@keepers, $group_keepers);
	push(@holdovers, $group_holdovers);
    }

    #keepers and hit holdovers are returned in same order given by user
    return @keepers, @holdovers;
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
sub maker_p_and_t_fastas {
   my $maker    = shift @_;
   my $non_over = shift @_;
   my $abinit   = shift @_;
   my $p_fastas = shift @_;
   my $t_fastas = shift @_;
   
   foreach my $an (@$maker) {
      foreach my $a (@{$an->{t_structs}}) {
	 my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	 $p_fastas->{maker} .= $p_fasta;
	 $t_fastas->{maker} .= $t_fasta;
      }
   }

   foreach my $an (@$non_over) {
       foreach my $a (@{$an->{t_structs}}) {
	   my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	   $p_fastas->{non_overlapping_ab_initio} .= $p_fasta;
	   $t_fastas->{non_overlapping_ab_initio} .= $t_fasta;
       }
   }
   
   foreach my $an (@$abinit) {
       foreach my $a (@{$an->{t_structs}}) {
	   my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	   my $source = $a->{hit}->algorithm;
	   $source =~ s/pred_gff\://;
	   $p_fastas->{$source} .= $p_fasta;
	   $t_fastas->{$source} .= $t_fasta;
       }
   }
}

#-----------------------------------------------------------------------------
sub get_p_and_t_fastas {
   my $t_struct = shift;
	
   my $t_seq  = $t_struct->{t_seq};
   my $p_seq  = $t_struct->{p_seq};
   my $t_off  = $t_struct->{t_offset};
   my $t_name = $t_struct->{t_name};
   my $AED    = $t_struct->{AED};
   my $QI     = $t_struct->{t_qi};
	
   my $p_def = '>'.$t_name.' protein AED:'.$AED.' QI:'.$QI; 
   my $t_def = '>'.$t_name.' transcript offset:'.$t_off.' AED:'.$AED.' QI:'.$QI;
	
   my $p_fasta = Fasta::toFasta($p_def, \$p_seq);
   my $t_fasta = Fasta::toFasta($t_def, \$t_seq);
	
   return($p_fasta, $t_fasta);
}
#----------------------------------------------------------------------------
sub write_p_and_t_fastas{
    my $p_fastas = shift @_;
    my $t_fastas = shift @_;
    my $safe_seq_id = shift @_;
    my $out_dir = shift @_;

    while( my $key = each %$p_fastas){
	my $name = "$out_dir/$safe_seq_id.maker";
	$name .= ".$key" unless($key eq 'maker');
	$name .= "\.proteins.fasta";

	FastaFile::writeFile(\$p_fastas->{$key},
			     $name,
	                    );
    }

    while( my $key = each %$t_fastas){
	my $name = "$out_dir/$safe_seq_id.maker";
	$name .= ".$key" unless($key eq 'maker');
	$name .= "\.transcripts.fasta";

	FastaFile::writeFile(\$t_fastas->{$key},
			     $name,
	                    );
    }
}
#----------------------------------------------------------------------------
sub create_blastdb {
   my $CTL_OPT = shift @_;
   my $mpi_size = shift@_ || 1;
       
   #rebuild all fastas when specified
   File::Path::rmtree($CTL_OPT->{out_base}."/mpi_blastdb") if ($CTL_OPT->{force} &&
							       ! $CTL_OPT->{_multi_chpc});
   
   ($CTL_OPT->{_protein}, $CTL_OPT->{p_db}) = split_db($CTL_OPT, 'protein', $mpi_size);
   ($CTL_OPT->{_est}, $CTL_OPT->{e_db}) = split_db($CTL_OPT, 'est', $mpi_size);
   ($CTL_OPT->{_est_reads},  $CTL_OPT->{d_db}) = split_db($CTL_OPT, 'est_reads', $mpi_size);
   ($CTL_OPT->{_altest},  $CTL_OPT->{a_db}) = split_db($CTL_OPT, 'altest', $mpi_size);
   ($CTL_OPT->{_repeat_protein}, $CTL_OPT->{r_db}) = split_db($CTL_OPT, 'repeat_protein', $mpi_size);
}
#----------------------------------------------------------------------------
sub concatenate_files {
    my $infiles = shift;
    my $outfile = shift;

    my ($tFH, $t_file) = tempfile(DIR => $TMP);
    foreach my $file (@$infiles){
	open(my $IN, "< $file");
	while(defined(my $line = <$IN>)){
	    print $tFH $line;
	}
	close($IN);
    }
    close($tFH);
    File::Copy::move($t_file, $outfile);
}
#----------------------------------------------------------------------------
sub split_db {
   my $CTL_OPT  = shift @_;
   my $key      = shift @_;
   my $mpi_size = shift @_ || 1;

   #always set to at least 10 for faster fasta indexing
   $mpi_size = 10 if($mpi_size < 10);

   my $file = $CTL_OPT->{$key};
   my $alt = $CTL_OPT->{alt_peptide} if($key =~ /protein/);
   
   return ('', []) if (not $file);
   
   #set up names and variables
   my $fasta_iterator = new Iterator::Fasta($file);
   my $db_size = $fasta_iterator->number_of_entries();
   my $bins = $mpi_size;
   $bins = $db_size if ($db_size < $bins);

   my ($f_name) = $file =~ /([^\/]+)$/;
   $f_name =~ s/\.fasta$//;
    
   my $d_name = "$f_name\.mpi\.$mpi_size";
   my $b_dir = $CTL_OPT->{out_base}."/mpi_blastdb";
   my $f_dir = "$b_dir/$d_name";
   my $t_dir = $TMP."/$d_name";

   #make needed output directories
   mkdir($t_dir);
   mkdir($b_dir) unless (-e $b_dir);

   if(my $lock = new File::NFSLock($f_dir, 'EX', 600, 40)){
       $lock->maintain(30);
       
       if(-e "$f_dir"){ #on multi processors check if finished
	   my @t_db = <$f_dir/*$d_name*\.fasta>;
	   
	   my @existing_db;
	   foreach my $f (@t_db) {
	       push (@existing_db, $f) if (! -d $f);
	   }
	   
	   if(@existing_db == $bins){ #use existing if right count
	       $lock->unlock;
	       return $f_dir, \@existing_db;
	   }
	   else{ #remove if there is an error
	       File::Path::rmtree($f_dir);
	     }
       }
       
       #open filehandles for  pieces on multi processors
       my @fhs;

       for (my $i = 0; $i < $bins; $i++) {
	   my $name = "$t_dir/$d_name\.$i\.fasta";
	   my $fh;
	   open ($fh, "> $name");
	   
	   push (@fhs, $fh);
       }
       
       #write fastas here
       my %alias;
       
       my $wflag = 1; #flag set so warnings gets printed only once 
       while (my $fasta = $fasta_iterator->nextEntry()) {
	   my $def = Fasta::getDef(\$fasta);
	   my $seq_id = Fasta::def2SeqID($def);
	   my $seq_ref = Fasta::getSeqRef(\$fasta);
	   
	   #fix non standard peptides
	   if (defined $alt) {
	       $$seq_ref =~ s/[\*\-]//g;
	       $$seq_ref =~ s/[^abcdefghiklmnpqrstvwyzxACDEFGHIKLMNPQRSTVWYX\-\n]/$alt/g;
	   }
	   #fix nucleotide sequences
	   elsif($key !~ /protein/){
	       #most programs use N for masking but for some reason the NCBI decided to
	       #use X to mask their sequence, which causes many many programs to fail
	       $$seq_ref =~ s/\-//g;
	       $$seq_ref =~ s/X/N/g;
	       die "ERROR: The nucleotide sequence file \'$file\'\n".
		   "appears to contain protein sequence or unrecognized characters.\n".
		   "Please check/fix the file before continuing.\n".
		   "Invalid Character: $1\n\n"
		   if($$seq_ref =~ /([^acgturykmswbdhvnxACGTURYKMSWBDHVNX\-\n])/);
	   }
	   
	   #Skip empty fasta entries
	   next if($$seq_ref eq '');
	   
	   #fix weird blast trimming error for long seq IDs by replacing them
	   if(length($seq_id) > 78){
	       warn "WARNING: The fasta file contains sequences with names longer\n".
		   "than 78 characters.  Long names get trimmed by BLAST, making\n".
		   "it harder to identify the source of an alignmnet. You might\n".
		   "want to reformat the fasta file with shorter IDs.\n".
		   "File_name:$file\n\n" if($wflag-- > 0);
	       
	       my $new_id = uri_escape(Digest::MD5::md5_base64($seq_id), "^A-Za-z0-9\-\_");

	       die "ERROR: The id $seq_id is too long for BLAST, and I can'y uniquely fix it\n"
		   if($alias{$new_id});

	       $alias{$new_id}++;
	       $def =~ s/^(>\S+)/$1 MD5_alias=$new_id/;
	   }
	   
	   #reformat fasta, just incase
	   my $fasta_ref = Fasta::toFastaRef($def, $seq_ref);
	   
	   #build part files only on multi processor
	   my $fh = shift @fhs;
	   print $fh $$fasta_ref;
	   push (@fhs, $fh);
       }
       
       #close part file handles
       foreach my $fh (@fhs) {
	   close ($fh);
       }
       
       #move finished files into place
       system("mv $t_dir $f_dir");
       
       #check if everything is ok
       if (-e $f_dir) { #multi processor
	   my @t_db = <$f_dir/*$d_name*\.fasta>;
	   
	   my @db_files;
	   foreach my $f (@t_db) {
	       push (@db_files, $f) if (! -d $f);
	   }
	   
	   die "ERROR: SplitDB not created correctly\n\n" if(@db_files != $bins); #not o

	   $lock->unlock;
	   return $f_dir, \@db_files;
       }
       else {
	   die "ERROR: Could not split db\n"; #not ok
       }
   }
   else{
       die "ERROR: Could not get lock to process fasta\n\n";
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
   push(@preds, @{genemark(@_)});
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
   
   return [] if(! grep {/snap/} @{$CTL_OPT->{_run}});
   
   my %params;
   my $out_file = "$the_void/$seq_id\.all\.snap";
   
   $LOG->add_entry("STARTED", $out_file, "");   

   my $command  = $exe;
   $command .= " $snaphmm";
   $command .= " $in_file";
   $command .= " > $out_file";
	
   my $w = new Widget::snap();
	
   if (-e $out_file) {
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
sub genemark {
   my $in_file     = shift;
   my $the_void    = shift;
   my $seq_id      = shift;
   my $CTL_OPT     = shift;
   my $LOG         = shift;

   $in_file        = shift; #temp
   return []      if(! $in_file); #temp

   #genemark sometimes fails if called directly so I built a wrapper
   my $wrap = "$FindBin::Bin/../lib/Widget/genemark/gmhmm_wrap";
   my $exe  = $CTL_OPT->{organism_type} eq 'eukaryotic' ? $CTL_OPT->{gmhmme3} : $CTL_OPT->{gmhmmp}; #genemark
   my $pro = $CTL_OPT->{probuild}; #helper exe
   my $hmm = $CTL_OPT->{gmhmm};
   
   return [] if(! grep {/genemark/} @{$CTL_OPT->{_run}});
   
   my %params;
   my $out_file = "$the_void/$seq_id\.all\.genemark";
   
   $LOG->add_entry("STARTED", $out_file, "");   


   my $command  = $wrap;
   $command .= " -m $hmm";
   $command .= " -g $exe";
   $command .= " -p $pro";
   $command .= " -o $out_file";
   #$command .= " -t $TMP";
   $command .= " $in_file";

   my $w = new Widget::genemark();
	
   if (-e $out_file) {
      print STDERR "re reading genemark report.\n" unless $main::quiet;
      print STDERR "$out_file\n" unless $main::quiet;
   }
   else {
      print STDERR "running  genemark.\n" unless $main::quiet;
      $w->run($command);
   }
	
   $params{min_exon_score}  = -100000; #-10000;
   $params{min_gene_score}  = -100000; #0;
		
   my $keepers = Widget::genemark::parse($out_file,
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

   return [] if(! grep {/augustus/} @{$CTL_OPT->{_run}});

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.augustus";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = $exe;
   $command .= " --species=$org";
   $command .= " --UTR=off"; #added 3/19/2009
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::augustus();

   if (-e $out_file) {
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

   my $wrap = "$FindBin::Bin/../lib/Widget/fgenesh/fgenesh_wrap"; #fgenesh wrapper
   my $exe = $CTL_OPT->{fgenesh};
   my $org = $CTL_OPT->{fgenesh_par_file};

   return [] if(! grep {/fgenesh/} @{$CTL_OPT->{_run}});

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.fgenesh";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = "$wrap $exe";
   #$command .= " -tmp $TMP";
   $command .= " $org";
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::fgenesh();

   if (-e $out_file) {
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

   return [] if(! grep {/twinscan/} @{$CTL_OPT->{_run}});

   my %params;
   my $out_file    = "$the_void/$seq_id\.all\.twinscan";

   $LOG->add_entry("STARTED", $out_file, ""); 

   my $command  = $exe;
   $command .= ' --species='."$org";
   $command .= " $in_file";
   $command .= " > $out_file";

   my $w = new Widget::twinscan();

   if (-e $out_file) {
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
    my $g_fasta     = shift;
    my $phat_hits   = shift;
    my $db_index    = shift;
    my $the_void    = shift;
    my $type        = shift;
    my $exonerate   = shift;
    my $pcov        = shift;
    my $pid         = shift;
    my $score_limit = shift;
    my $matrix      = shift;
    my $est_forward = shift;
    my $LOG         = shift;
    
    my $def = Fasta::getDef($g_fasta);
    my $seq = Fasta::getSeqRef($g_fasta);
    
    my $exe = $exonerate;
    
    my @exonerate_data;
    
    foreach my $hit (@{$phat_hits}) {
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
	
	my $fastaObj = $db_index->get_Seq_for_hit($hit);
	
	#still no sequence? try rebuilding the index and try again
	if (not $fastaObj) {
	    print STDERR "WARNING: Cannot find> ".$hit->name.", trying to re-index the fasta.\n";
	    $db_index->reindex();
	    $fastaObj = $db_index->get_Seq_for_hit($hit);
	    
	    if (not $fastaObj) {
		print STDERR "stop here:".$hit->name."\n";
		die "ERROR: Fasta index error\n";
	    }
	}
	
	my $seq      = $fastaObj->seq();
	my $def      = $db_index->header_for_hit($hit);
	
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
					 $LOG
					 );
	
	foreach my $exonerate_hit (@{$exonerate_hits}) {
	    next if(! defined $exonerate_hit);

	    #fix flipped hits when mapping ESTs to gene models as is
	    if($est_forward && $exonerate_hit->num_hsps == 1 && $exonerate_hit->{_was_flipped}){
		$exonerate_hit = PhatHit_utils::copy($exonerate_hit, 'both');
		$exonerate_hit->{_was_flipped} = 0;
	    }

	    if (exonerate_okay($exonerate_hit)) {
		#tag the source blastn hit to let you know the counterpart
		#exonerate hit was flipped to the other strand
		$hit->{_exonerate_flipped} = 1 if($exonerate_hit->{_was_flipped});
		$hit->type("exonerate:$type"); #set hit type (exonerate only)
		
		push(@exonerate_data, $exonerate_hit);	       
	    }
	}
    }
    
    return \@exonerate_data;
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
    
    return 1;
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
   my $LOG   = shift;

   if ($type eq 'p') {
      return polisher::exonerate::protein::polish($d_file,
						  $t_file,
						  $the_void,
						  $offset,
						  $exe,
						  $score_limit,
						  $matrix,
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
	 my $fastaObj = $index->get_Seq_for_hit($hit);
	 
	 #still no sequence? try rebuilding the index and try again
	 if (not $fastaObj) {
	     print STDERR "WARNING: Cannot find> ".$hit->name.", trying to re-index the fasta.\n";
	     $index->reindex();
	     $fastaObj = $index->get_Seq_for_hit($hit);

	     if (not $fastaObj) {
		 print STDERR "stop here:".$hit->name."\n";
		 die "ERROR: Fasta index error\n";
	     }
	 }

	 my $seq      = $fastaObj->seq(); 
	 my $def      = $index->header_for_hit($hit);
	 my $fasta    = Fasta::toFasta('>'.$def, \$seq);
	 $fastas     .= $$fasta; 
      }
   }
   return \$fastas;
}
#-----------------------------------------------------------------------------
sub build_fasta_index {
   my $db = shift;

   my $index = new FastaDB($db);

   return $index;
}
#-----------------------------------------------------------------------------
sub build_all_indexes {
   my $CTL_OPT = shift;

   my @dbs = ($CTL_OPT->{_est},
	      $CTL_OPT->{_protein},
	      $CTL_OPT->{_repeat_protein},
	      $CTL_OPT->{_est_reads},
	      $CTL_OPT->{_altest}
	     );

   foreach my $db (@dbs){
       next if(! $db);
       my $index = build_fasta_index($db);
       $index->reindex() if($CTL_OPT->{force} && !$CTL_OPT->{_multi_chpc});
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

   my $lock;
   unless($lock = new File::NFSLock("$file.dbformat", 'EX', 600, 40)){
       die "ERROR:  Could not obtain lock to format database\n\n";
   }

   if ($command =~ /xdformat/) {
      if (($type eq 'blastn' && ! -e $file.'.xnd') ||
	  ($type eq 'blastx' && ! -e $file.'.xpd') ||
	  ($type eq 'tblastx' && ! -e $file.'.xnd')
	 ) {
	 $command .= " -p" if($type eq 'blastx');
	 $command .= " -n" if($type eq 'blastn' || $type eq 'tblastx');
	 $command .= " $file";

	 $lock->maintain(30);
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

	 $lock->maintain(30);
	 my $w = new Widget::formatdb();
	 print STDERR "formating database...\n" unless $main::quiet;
	 $w->run($command);
      }
   }
   else {
      die "ERROR: databases can only be formated by xdformat or formatdb not \'$command\'\n";
   }

   $lock->unlock;
}
#-----------------------------------------------------------------------------
sub blastn_as_chunks {
   my $chunk      = shift;
   my $db         = shift;
   my $old_db     = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $rank       = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;

   my $blastn      = $CTL_OPT->{_blastn};
   my $bit_blastn  = $CTL_OPT->{bit_blastn};
   my $eval_blastn = $CTL_OPT->{eval_blastn};
   my $pcov_blastn = $CTL_OPT->{pcov_blastn};
   my $pid_blastn  = $CTL_OPT->{pid_blastn};
   my $split_hit   = $CTL_OPT->{split_hit};
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

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
   my $tmp_db = "$TMP/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #copy db to local tmp dir and run xdformat or formatdb
   if (((! @{[<$tmp_db.x?d*>]} && $formater =~ /xdformat/) ||
        (! @{[<$tmp_db.?sq*>]} && $formater =~ /formatdb/)) &&
       (! -e $blast_finished)
      ){
       if(my $lock = new File::NFSLock("$tmp_db.copy", 'EX', 600, 40)){
	   $lock->maintain(30);
	   copy($db, $tmp_db) if(! -e $tmp_db);
	   dbformat($formater, $tmp_db, 'blastn');
	   $lock->unlock;
       }
       else{
	   die "ERROR: Could not get lock.\n\n";
       }
   }
   elsif (-e $blast_finished) {
      print STDERR "re reading blast report.\n" unless ($main::quiet || !$LOG_FLAG);
      print STDERR "$blast_finished\n" unless ($main::quiet || !$LOG_FLAG);
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
	     $org_type,
	     $softmask
	    );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_blastn{
   my $chunk     = shift;
   my $blast_dir = shift;
   my $CTL_OPT   = shift;
   my $LOG       = shift;

   my $eval_blastn = $CTL_OPT->{eval_blastn};
   my $bit_blastn  = $CTL_OPT->{bit_blastn},
   my $pcov_blastn = $CTL_OPT->{pcov_blastn};
   my $pid_blastn  = $CTL_OPT->{pid_blastn};
   my $split_hit   = $CTL_OPT->{split_hit};

   my $blast_finished = $blast_dir;
   $blast_finished =~ s/\.temp_dir$//;

   #merge blast reports
   if (! -e $blast_finished) {
      system ("cat $blast_dir/*blast* > $blast_finished");
      File::Path::rmtree ("$blast_dir");
   }

   my %params;
   $params{significance}  = $eval_blastn;
   $params{hsp_bit_min}   = $bit_blastn;
   $params{percov}        = $pcov_blastn;
   $params{percid}        = $pid_blastn;
   $params{split_hit}     = $split_hit;

   $LOG->add_entry("FINISHED", $blast_finished, "");

   my $chunk_keepers = Widget::blastn::parse($blast_finished,
					     \%params,
					    );
   
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
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;

   my $blastn      = $CTL_OPT->{_blastn};
   my $bit_blastn  = $CTL_OPT->{bit_blastn};
   my $eval_blastn = $CTL_OPT->{eval_blastn};
   my $pcov_blastn = $CTL_OPT->{pcov_blastn};
   my $pid_blastn  = $CTL_OPT->{pid_blastn};
   my $split_hit   = $CTL_OPT->{split_hit};
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

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
	     $org_type,
	     $softmask
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
   my $q_file     = shift;
   my $db         = shift;
   my $out_file   = shift;
   my $blast      = shift;
   my $eval_blast = shift;
   my $split_hit  = shift;
   my $cpus       = shift;
   my $org_type   = shift;
   my $softmask   = shift;

   my $command  = $blast;
   if ($command =~ /blastn$/) {
      $command .= " $db $q_file B=100000 V=100000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";;
      $command .= " R=3";
      $command .= " W=15";
      $command .= " M=1";
      $command .= " N=-3";
      $command .= " Q=3";
      $command .= " Z=1000";
      $command .= ($org_type eq 'eukaryotic') ? " Y=500000000" : " Y=20000000";
      $command .= " cpus=$cpus";	
      $command .= ($org_type eq 'eukaryotic') ? " topcomboN=1" : "";
      $command .= ($org_type eq 'eukaryotic') ? " hspmax=100" : " hspmax=5";
      $command .= ($org_type eq 'eukaryotic') ? " gspmax=100" : " gspmax=5";
      $command .= ($org_type eq 'eukaryotic') ? " hspsepqmax=$split_hit" : "";
      $command .= " lcmask";
      $command .= " maskextra=10";
      $command .= " gi";
      $command .= " warnings"; #suppress certain warnings
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      $command .= ($org_type eq 'eukaryotic') ? "" : " kap";
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
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= ($org_type eq 'eukaryotic') ? " -K 100" : " -K 5";
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
   if (-e $out_file) {
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
   my $chunk      = shift;
   my $db         = shift;
   my $old_db     = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $rank       = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;

   my $rflag = 1 if($old_db && $CTL_OPT->{repeat_protein} eq $old_db); #am I running repeat data?

   my $blastx      = $CTL_OPT->{_blastx};
   my $bit_blastx  = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{bit_blastx};
   my $eval_blastx = ($rflag) ? $CTL_OPT->{eval_rm_blastx} : $CTL_OPT->{eval_blastx};
   my $pcov_blastx = ($rflag) ? $CTL_OPT->{pcov_rm_blastx} : $CTL_OPT->{pcov_blastx};
   my $pid_blastx  = ($rflag) ? $CTL_OPT->{pid_rm_blastx} : $CTL_OPT->{pid_blastx};
   my $split_hit   = $CTL_OPT->{split_hit}; #repeat proteins get shatttered later anyway
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = ($rflag) ? 1 : $CTL_OPT->{softmask}; #always on for repeats
   my $org_type    = $CTL_OPT->{organism_type};

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
   my $tmp_db = "$TMP/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG);

   #copy db to local tmp dir and run xdformat or format db 
   if (((! @{[<$tmp_db.x?d*>]} && $formater =~ /xdformat/) ||
	(! @{[<$tmp_db.?sq*>]} && $formater =~ /formatdb/)) &&
       (! -e $blast_finished)
      ){
       if(my $lock = new File::NFSLock("$tmp_db.copy", 'EX', 600, 40)){
	   $lock->maintain(30);
           copy($db, $tmp_db) if(! -e $tmp_db);
	   dbformat($formater, $tmp_db, 'blastx');
           $lock->unlock;
       }
       else{
	   die "ERROR: Could not get lock.\n\n";
       }
   }
   elsif (-e $blast_finished) {
      print STDERR "re reading blast report.\n" unless ($main::quiet || !$LOG_FLAG);
      print STDERR "$blast_finished\n" unless ($main::quiet || !$LOG_FLAG);
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
	     $org_type,
	     $softmask
	    );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_blastx{
   my $chunk     = shift;
   my $blast_dir = shift;
   my $CTL_OPT   = shift;
   my $LOG       = shift;

   my $eval_blastx = $CTL_OPT->{eval_blastx};
   my $bit_blastx  = $CTL_OPT->{bit_blastx},
   my $pcov_blastx = $CTL_OPT->{pcov_blastx};
   my $pid_blastx  = $CTL_OPT->{pid_blastx};
   my $split_hit   = $CTL_OPT->{split_hit};

   my $blast_finished = $blast_dir;
   $blast_finished =~ s/\.temp_dir$//;

   #merge blast reports
   if (! -e $blast_finished) {
      system ("cat $blast_dir/*blast* > $blast_finished");
      rmtree ("$blast_dir");
   }

   my %params;
   $params{significance}  = $eval_blastx;
   $params{hsp_bit_min}   = $bit_blastx;
   $params{percov}        = $pcov_blastx;
   $params{percid}        = $pid_blastx;
   $params{split_hit}     = $split_hit;

   $LOG->add_entry("FINISHED", $blast_finished, "");
   
   my $chunk_keepers = Widget::blastx::parse($blast_finished,
					     \%params,
					    );

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
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;

   my $rflag = 1 if($db && $CTL_OPT->{repeat_protein} eq $db); #am I running repeat data?

   my $blastx      = $CTL_OPT->{_blastx};
   my $bit_blastx  = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{bit_blastx};
   my $eval_blastx = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{eval_blastx};
   my $pcov_blastx = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{pcov_blastx};
   my $pid_blastx  = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{pid_blastx};
   my $split_hit   = $CTL_OPT->{split_hit}; #repeat proteins get shatttered later anyway
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = ($rflag) ? 1 : $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

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
	     $org_type,
	     $softmask
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
   my $org_type = shift;
   my $softmask = shift;


   my $command  = $blast;
   if ($command =~ /blastx$/) {
      $command .= " $db $q_file B=10000 V=10000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";
      #$command .= " T=20";
      #$command .= " W=5";
      #$command .= " wink=5";
      $command .= " Z=300";
      $command .= ($org_type eq 'eukaryotic') ? " Y=500000000" : " Y=20000000";
      $command .= ($org_type eq 'eukaryotic') ? " hspmax=100" : " hspmax=5";
      $command .= " cpus=$cpus";
      $command .= ($org_type eq 'eukaryotic') ? " gspmax=100" : " gspmax=5";
      #$command .= " hspsepqmax=10000";
      $command .= " lcmask";
      $command .= " maskextra=10";
      $command .= " kap";
      $command .= " gi";
      $command .= " warnings"; #suppress certain warnings
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p blastx";
      $command .= " -d $db -i $q_file -b 100000 -v 100000 -e $eval_blast";
      $command .= " -z 300";
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= ($org_type eq 'eukaryotic') ? " -K 100" : " -K 5";
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

   if (-e $out_file) {
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
   my $old_db     = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $rank       = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;

   my $tblastx      = $CTL_OPT->{_tblastx};
   my $bit_tblastx  = $CTL_OPT->{bit_tblastx};
   my $eval_tblastx = $CTL_OPT->{eval_tblastx};
   my $pcov_tblastx = $CTL_OPT->{pcov_tblastx};
   my $pid_tblastx  = $CTL_OPT->{pid_tblastx};
   my $split_hit    = $CTL_OPT->{split_hit};
   my $cpus         = $CTL_OPT->{cpus};
   my $formater     = $CTL_OPT->{_formater};
   my $softmask     = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

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
   my $tmp_db = "$TMP/$1";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #copy db to local tmp dir and run xdformat or formatdb
   if (((! @{[<$tmp_db.x?d*>]} && $formater =~ /xdformat/) ||
        (! @{[<$tmp_db.?sq*>]} && $formater =~ /formatdb/)) &&
       (! -e $blast_finished)
      ){
       if(my $lock = new File::NFSLock("$tmp_db.copy", 'EX', 600, 40)){
	   $lock->maintain(30);
           copy($db, $tmp_db) if(! -e $tmp_db);
	   dbformat($formater, $tmp_db, 'tblastx');
           $lock->unlock;
       }
       else{
	   die "ERROR: Could not get lock.\n\n";
       }
   }
   elsif (-e $blast_finished) {
      print STDERR "re reading blast report.\n" unless ($main::quiet || !$LOG_FLAG);
      print STDERR "$blast_finished\n" unless ($main::quiet || !$LOG_FLAG);
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
	      $org_type,
	      $softmask
	     );

   $chunk->erase_fasta_file();

   return $blast_dir;
}
#-----------------------------------------------------------------------------
sub collect_tblastx{
   my $chunk     = shift;
   my $blast_dir = shift;
   my $CTL_OPT   = shift;
   my $LOG       = shift;

   my $eval_tblastx = $CTL_OPT->{eval_tblastx};
   my $bit_tblastx  = $CTL_OPT->{bit_tblastx},
   my $pcov_tblastx = $CTL_OPT->{pcov_tblastx};
   my $pid_tblastx  = $CTL_OPT->{pid_tblastx};
   my $split_hit   = $CTL_OPT->{split_hit};

   my $blast_finished = $blast_dir;
   $blast_finished =~ s/\.temp_dir$//;

   #merge blast reports
   if (! -e $blast_finished) {
      system ("cat $blast_dir/*blast* > $blast_finished");
      File::Path::rmtree ("$blast_dir");
   }

   my %params;
   $params{significance}  = $eval_tblastx;
   $params{hsp_bit_min}   = $bit_tblastx;
   $params{percov}        = $pcov_tblastx;
   $params{percid}        = $pid_tblastx;
   $params{split_hit}     = $split_hit;

   $LOG->add_entry("FINISHED", $blast_finished, "");

   my $chunk_keepers = Widget::tblastx::parse($blast_finished,
					      \%params,
					     );
   
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
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;

   my $tblastx      = $CTL_OPT->{_tblastx};
   my $bit_tblastx  = $CTL_OPT->{bit_tblastx};
   my $eval_tblastx = $CTL_OPT->{eval_tblastx};
   my $pcov_tblastx = $CTL_OPT->{pcov_tblastx};
   my $pid_tblastx  = $CTL_OPT->{pid_tblastx};
   my $split_hit    = $CTL_OPT->{split_hit};
   my $cpus         = $CTL_OPT->{cpus};
   my $formater     = $CTL_OPT->{_formater};
   my $softmask     = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

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
	      $org_type,
	      $softmask
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

   return $chunk_keepers;
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
   my $org_type = shift;
   my $softmask = shift;


   my $command  = $blast;
   if ($command =~ /tblastx$/) {
      $command .= " $db $q_file B=100000 V=100000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";
      #$command .= " W=15";
      $command .= " Z=1000";
      $command .= ($org_type eq 'eukaryotic') ? " Y=500000000" : " Y=20000000";
      $command .= " cpus=$cpus";	
      $command .= ($org_type eq 'eukaryotic') ? " topcomboN=1" : "";
      $command .= ($org_type eq 'eukaryotic') ? " hspmax=100" : " hspmax=5";
      $command .= ($org_type eq 'eukaryotic') ? " gspmax=100" : " gspmax=5";
      $command .= ($org_type eq 'eukaryotic') ? " hspsepqmax=$split_hit" : "";
      $command .= " lcmask";
      $command .= " maskextra=10";
      $command .= " gi";
      $command .= " warnings"; #suppress certain warnings
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      $command .= ($org_type eq 'eukaryotic') ? "" : " kap";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p tblastx";
      $command .= " -d $db -i $q_file -b 100000 -v 100000 -e $eval_blast";
      #$command .= " -W 15";
      $command .= " -z 1000";
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= ($org_type eq 'eukaryotic') ? " -K 100" : " -K 5";
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
   if (-e $out_file) {
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
   if (-e $o_file) {
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
   my $user_default = shift; #hash ref
   
   if ($type !~ /^all$|^opts$|^bopts$|^exe$|^menus$|^server$/) {
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
      $CTL_OPT{'model_pass'} = 0;
      $CTL_OPT{'pred_pass'} = 0;
      $CTL_OPT{'other_pass'} = 0;
      $CTL_OPT{'est'} = '';
      $CTL_OPT{'est_reads'} = '';
      $CTL_OPT{'altest'} = '';
      $CTL_OPT{'est_gff'} = '';
      $CTL_OPT{'altest_gff'} = '';
      $CTL_OPT{'protein'} = '';
      $CTL_OPT{'protein_gff'} = '';
      $CTL_OPT{'model_org'} = 'all';
      $CTL_OPT{'model_org'} .= '=STATIC' if($main::server);
      $CTL_OPT{'repeat_protein'} = Cwd::abs_path("$FindBin::Bin/../data/te_proteins.fasta");
      $CTL_OPT{'repeat_protein'} .= "=STATIC" if($main::server);
      $CTL_OPT{'rmlib'} = '';
      $CTL_OPT{'rm_gff'} = '';
      $CTL_OPT{'organism_type'} = 'eukaryotic';
      $CTL_OPT{'predictor'} = '';
      $CTL_OPT{'predictor'} = 'model_gff' if($main::eva);
      $CTL_OPT{'snaphmm'} = '';
      $CTL_OPT{'gmhmm'} = '';
      $CTL_OPT{'gmhmm_e'} = '' if($main::server);
      $CTL_OPT{'gmhmm_p'} = '' if($main::server);
      $CTL_OPT{'augustus_species'} = '';
      $CTL_OPT{'fgenesh_par_file'} = '';
      $CTL_OPT{'model_gff'} = '';
      $CTL_OPT{'pred_gff'} = '';
      $CTL_OPT{'other_gff'} = '';
      $CTL_OPT{'domain'} = '0';
      $CTL_OPT{'function'} = '0';
      $CTL_OPT{'short_name'} = '';
      $CTL_OPT{'alt_peptide'} = 'C';
      $CTL_OPT{'cpus'} = 1;
      $CTL_OPT{'cpus'} .= '=DISABLED' if($main::server);
      $CTL_OPT{'evaluate'} = 0;
      $CTL_OPT{'evaluate'} = 1 if($main::eva);
      $CTL_OPT{'max_dna_len'} = 100000;
      $CTL_OPT{'min_contig'} = 1;
      $CTL_OPT{'split_hit'} = 10000;
      $CTL_OPT{'softmask'} = 1;
      $CTL_OPT{'pred_flank'} = 200;
      $CTL_OPT{'single_exon'} = 0;
      $CTL_OPT{'single_length'} = 250;
      $CTL_OPT{'min_protein'} = 0;
      $CTL_OPT{'AED_threshold'} = 1;
      $CTL_OPT{'keep_preds'} = 0;
      $CTL_OPT{'map_forward'} = 0;
      $CTL_OPT{'retry'} = 1;
      $CTL_OPT{'clean_try'} = 0;
      $CTL_OPT{'TMP'} = '';
      $CTL_OPT{'TMP'} .= '=DISABLED' if($main::server);
      $CTL_OPT{'run'} = ''; #hidden option
      $CTL_OPT{'unmask'} = 0;
      $CTL_OPT{'clean_up'} = 0;
      $CTL_OPT{'clean_up'} = 1 if($main::server);
      #evaluator below here
      $CTL_OPT{'side_thre'} = 5;
      $CTL_OPT{'eva_window_size'} = 70;
      $CTL_OPT{'eva_split_hit'} = 1;
      $CTL_OPT{'eva_hspmax'} = 100;
      $CTL_OPT{'eva_gspmax'} = 100;
      $CTL_OPT{'enable_fathom'} = 0;
      $CTL_OPT{'enable_fathom'} = 1 if($main::eva);
   }

   #maker_bopts
   if ($type eq 'all' || $type eq 'bopts') {
      $CTL_OPT{'blast_type'} = 'wublast';
      $CTL_OPT{'blast_type'} .= '=DISABLED' if($main::server);
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
      #evaluator below here
      $CTL_OPT{'eva_pcov_blastn'} = 0.80;
      $CTL_OPT{'eva_pid_blastn'} = 0.85;
      $CTL_OPT{'eva_eval_blastn'} = 1e-10;
      $CTL_OPT{'eva_bit_blastn'} = 40;
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
		  'gmhmme3',
		  'gmhmmp',
		  'augustus',
		  'fgenesh',
		  'twinscan',
		  'jigsaw',
		  'qrna',
		  'fathom',
		  'probuild'
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

   #server
   if ($type eq 'server') {
      $CTL_OPT{'DBI'} = 'SQLite';
      $CTL_OPT{'dbname'} = 'makerweb';
      $CTL_OPT{'host'} = '';
      $CTL_OPT{'port'} = '';
      $CTL_OPT{'username'} = 'maker';
      $CTL_OPT{'password'} = '';
      $CTL_OPT{'admin_email'} = '';
      $CTL_OPT{'smtp_server'} = '';
      $CTL_OPT{'MPI'} = 0;
      $CTL_OPT{'mpiexec'} = 'mpiexec';
      $CTL_OPT{'max_cpus'} = 1;
      $CTL_OPT{'job_cpus'} = 1;
      $CTL_OPT{'use_login'} = 1;
      $CTL_OPT{'allow_guest'} = 1;
      $CTL_OPT{'allow_register'} = 1;
      $CTL_OPT{'tutorials'} = 1;
      $CTL_OPT{'max_submit_user'} = 2000000; #length in base pairs
      $CTL_OPT{'max_submit_guest'} = 200000; #length in base pairs
      $CTL_OPT{'persist_user'} = 336; #in hours
      $CTL_OPT{'persist_guest'} = 72; #in hours
      $CTL_OPT{'inactive_user'} = 0; #in days
      $CTL_OPT{'inactive_guest'} = 14; #in days
      $CTL_OPT{'cgi_dir'} = '/var/www/cgi-bin';
      $CTL_OPT{'cgi_dir'} = '/Library/WebServer/CGI-Executables' if(! -d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_dir'} = '/usr/lib/cgi-bin' if(! -d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_dir'} = '' if(! -d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_dir'} .= '/maker' if(-d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_web'} = '/cgi-bin/maker';
      $CTL_OPT{'html_dir'} = '/var/www/html';
      $CTL_OPT{'html_dir'} = '/Library/WebServer/Documents' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} = '/var/www' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} = '' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} .= '/maker' if(-d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_web'} = '/maker';
      $CTL_OPT{'data_dir'} = '';
      $CTL_OPT{'data_dir'} = "$CTL_OPT{html_dir}/data" if($CTL_OPT{html_dir});
      $CTL_OPT{'web_address'} = 'http://'.[`hostname` =~ /^([^\n]+)/]->[0];
      $CTL_OPT{'apache_user'} = '';
      $CTL_OPT{'apache_user'} = 'apache' if(@{[getpwnam('apache')]});
      $CTL_OPT{'apache_user'} = 'www' if(@{[getpwnam('www')]});
      $CTL_OPT{'apache_user'} = 'www-data' if(@{[getpwnam('www-data')]});
      $CTL_OPT{'font_file'} = '/usr/share/fonts/bitstream-vera/VeraMono.ttf';
      $CTL_OPT{'font_file'} = '/Library/Fonts/Verdana.ttf' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'font_file'} = '/usr/share/fonts/truetype/freefont/FreeMono.ttf' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'font_file'} = '' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'soba_url'} = 'http://www.sequenceontology.org/cgi-bin/soba.cgi';
      $CTL_OPT{'APOLLO_ROOT'} = $ENV{APOLLO_ROOT} || '';
      $CTL_OPT{'JBROWSE_ROOT'} = '';
      $CTL_OPT{'GBROWSE_MASTER'} = '/etc/gbrowse/GBrowse.conf';
      $CTL_OPT{'GBROWSE_MASTER'} = '' if(! -f $CTL_OPT{'GBROWSE_MASTER'});
   }

   #server menus
   if ($type eq 'menus') {
      #this step is required since some defaults are dependent on server setting dependant
      my %server_ctl = set_defaults('server', $user_default); 

      #now add static defaults
      $CTL_OPT{'genome'}           = {'D. melanogaster : example contig' => "$server_ctl{data_dir}/maker/MWAS/../data/dpp_contig.fasta",
				      'De novo Annotation : example contig' => "$server_ctl{data_dir}/maker/MWAS/data/pyu-contig.fasta",
				      'Pass-through : example contig' => "$server_ctl{data_dir}/maker/MWAS/data/pass-contig.fasta",
				      'Legacy Annotation : example contig' => "$server_ctl{data_dir}/maker/MWAS/data/legacy-contig.fasta",
				      'E. coli : example contig' => "$server_ctl{data_dir}/maker/MWAS/data/ecoli-contig.fasta"};
      $CTL_OPT{'snaphmm'}          = {'P. ultimum' => "$server_ctl{data_dir}/maker/MWAS/data/pyu.hmm"};
      $CTL_OPT{'augustus_species'} = {};
      $CTL_OPT{'fgenesh_par_file'} = {};
      $CTL_OPT{'gmhmm_e'}          = {'P. ultimum' => "$server_ctl{data_dir}/maker/MWAS/data/pyu.mod"};
      $CTL_OPT{'gmhmm_p'}          = {'E. coli' => "$server_ctl{data_dir}/maker/MWAS/data/ecoli.mod"};
      $CTL_OPT{'est'}              = {'D. melanogaster : example cDNA' => "$server_ctl{data_dir}/maker/MWAS/../data/dpp_transcripts.fasta",
				      'De novo/Legacy/Pass-through : example ESTs' => "$server_ctl{data_dir}/maker/MWAS/data/pyu-est.fasta",
				      'E. coli : example ESTs' => "$server_ctl{data_dir}/maker/MWAS/data/ecoli-est.fasta"};
      $CTL_OPT{'altest'}           = {};
      $CTL_OPT{'protein'}          = {'D. melanogaster : example proteins' => "$server_ctl{data_dir}/maker/MWAS/../data/dpp_proteins.fasta",
				      'E. coli : example proteins' => "$server_ctl{data_dir}/maker/MWAS/data/ecoli-protein.fasta",
				      'De novo/Legacy/Pass-through : example proteins' => "$server_ctl{data_dir}/maker/MWAS/data/pyu-protein.fasta"};
      $CTL_OPT{'est_gff'}          = {'Pass-through : example mRNAseq' => "$server_ctl{data_dir}/maker/MWAS/data/pass-mRNAseq.gff"};
      $CTL_OPT{'altest_gff'}       = {};
      $CTL_OPT{'protein_gff'}      = {};
      $CTL_OPT{'pred_gff'}         = {};
      $CTL_OPT{'model_gff'}        = {'Legacy Annotation : example model set 1' => "$server_ctl{data_dir}/maker/MWAS/data/legacy-set1.gff",
				      'Legacy Annotation : example model set 2' => "$server_ctl{data_dir}/maker/MWAS/data/legacy-set2.gff"};
      $CTL_OPT{'repeat_gff'}       = {};
      $CTL_OPT{'rmlib'}            = {};
      $CTL_OPT{'repeat_protein'}   = {'RepeatRunner te_proteins' => "$server_ctl{data_dir}/maker/MWAS/../data/te_proteins.fasta"};
      $CTL_OPT{'model_org'}        = {'All species' => 'all',
				      'Fungi' => 'fungi',
				      'Deuterostomes' => 'deuterostomes',
				      'Protostomes' => 'protostomes',
				      'Drosophila' => 'drosophila',
				      'Human' => 'human',
				      'Mouse' => 'mouse',
				      'Nematode' => 'nematode',
				      'Vertibrates' => 'vertibrates',
				      'Plants' => 'plants'};


      #auto add uniprot if user downloaded it into data directory
      $CTL_OPT{'protein'}{'UniProt'} = "$server_ctl{data_dir}/maker/MWAS/data/uniprot_sprot.fasta"
	  if(-f "$server_ctl{data_dir}/maker/MWAS/data/uniprot_sprot.fasta");

      #this step is required since menu defaults are exe dependant
      my %exe_ctl = set_defaults('exe', $user_default); 
      my %hmm_ctl = %{collect_hmms(\%exe_ctl)};
      while(my $key = each %hmm_ctl){
	  $CTL_OPT{$key} = {} if(! $CTL_OPT{$key});
	  %{$CTL_OPT{$key}} = (%{$CTL_OPT{$key}}, %{$hmm_ctl{$key}}); #add exe dependant values
      }

      #restore any user supplied values
#      if($user_default->{menus}){
#	  my %user_ctl = %{$user_default->{menus}};
#	  while(my $key = each %user_ctl){
#	      $CTL_OPT{$key} = {} if(! $CTL_OPT{$key});
#	      %{$CTL_OPT{$key}} = (%{$CTL_OPT{$key}}, %{$user_ctl{$key}});
#	  }
#      }
   }
   #reset values with user supplied defaults
   if($user_default && $type ne 'menus'){
       while(my $key = each %$user_default){
	   #will ignore invalid/inappropriate entries
	   $CTL_OPT{$key} = $user_default->{$key} if(exists $CTL_OPT{$key});
       }
   }

   return %CTL_OPT;
}
#-----------------------------------------------------------------------------
#this function will collect HMM file names and there locations for different
#prediction algorithms, and it returns the data in a hash reference
sub collect_hmms {
    my %exes = %{shift @_}; #make a copy of the location of executables

    my %hmms; #hold the hmms for each algorithm

    #find augustus HMMs
    if(-f $exes{augustus}){
	open(my $EXE, "$exes{augustus} --species=help 2>&1|");
	my $flag;
	while(my $line = <$EXE>){
	    if($flag){
		chomp $line;
		my ($value, $name) = split(/\|/, $line);
		next if(!$value || !$name);

		$name =~ s/^[\s\t\n]+|[\s\t\n]+$//g;
		$value =~ s/^[\s\t\n]+|[\s\t\n]+$//g;

		next if( $value =~ /^\(/ );

		$hmms{augustus_species}{$name} = $value;
	    }

	    $flag++ if($line =~ /----------------/);
	}
	close($EXE);
    }

    #find snap HMMs
    if(defined $ENV{ZOE} && -d "$ENV{ZOE}/HMM/"){
	$exes{snap} = "$ENV{ZOE}/HMM/";
    }
    elsif($exes{snap}){
	$exes{snap} =~ s/[^\/]+$//;
	$exes{snap} = "$exes{snap}/HMM/";
    }
    if(-d $exes{snap}){
	foreach my $file (grep {!/README/} <$exes{snap}/*>){
	    my ($name) = $file =~ /([^\/]+)$/;
	    $name =~ s/\.hmm$//;
	    $name =~ s/\./\. /;
	    my $value = Cwd::abs_path("$file");

	    next if(!$name || !$value);

	    $hmms{snaphmm}{$name} = $value;
	}
    }

    #find genemark Eukaryotic HMMs
    if($exes{gmhmme3}){
	$exes{gmhmme3} =~ s/[^\/]+$//;
	$exes{gmhmme3} = "$exes{gmhmme3}/HMM/";
    }
    if(-d $exes{gmhmme3}){
	foreach my $file (<$exes{gmhmme3}/*.mod>){
	    my ($name) = $file =~ /([^\/]+)\.mod$/;
	    $name =~ s/_/\. /;
	    $name = ucfirst($name);
	    my $value = Cwd::abs_path("$file");

	    next if(!$name || !$value);

	    $hmms{gmhmm_e}{$name} = $value;
	}
    }

    #find fgenesh HMMs
    if($exes{fgenesh}){
	$exes{fgenesh} =~ s/[^\/]+$//;
    }
    if(-d $exes{fgenesh}){
	foreach my $file (<$exes{fgenesh}/*>){
	    my ($name) = $file =~ /([^\/]+)$/;
	    my $value =Cwd::abs_path("$file");

	    next if(!$name || !$value);

	    my $filesize = [stat($file)]->[7]; #size in bytes

        #rough size of all parameter files
	    $hmms{fgenesh_par_file}{$name} = $value if($filesize >= 200000);
	}
    }

    return \%hmms;
}
#-----------------------------------------------------------------------------
#this function parses the control files and does no error checking
#this is dirty method to load subsections of control files
sub parse_ctl_files {
    my @ctlfiles = @{shift @_};

    my %CTL_OPT;
    #--load values from control files
    foreach my $ctlfile (@ctlfiles) {
	open (CTL, "< $ctlfile") or die"ERROR: Could not open the control file \"$ctlfile\".\n";
	
	while (my $line = <CTL>) {
	    chomp($line);
	   
	    if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:]+)\:([^\n\#]*)/) {
		my $key = $1;
		my $value = $2;
		my $stat;

		#remove preceding and trailing whitespace
		$value =~ s/^[\s\t]+|[\s\t]+$//g;

		($value, $stat) = split("=", $value);

		#set value
		$CTL_OPT{$key} = defined($value) ? $value : '';
		$CTL_OPT{STAT}{$key} = defined($stat) ? $stat : ''
		    if($stat);
		
	    }#now load menus
	    elsif($line =~ /^\#\#Menus from Data::Dumper/){ #only non-sandard format control file
		my $data = join('', <CTL>);
		my $menus; #will be set by data
		eval $data;

		while(my $key = each %$menus){
		    #set value
		    $CTL_OPT{menus}{$key} = $menus->{$key};
		}
	    }
	}
    }

    return %CTL_OPT;
}
#-----------------------------------------------------------------------------
#this function parses the control files and sets up options for each maker run
#error checking for starup occurs here
sub load_server_files {
    my @ctlfiles = @{shift @_};

    #make list of permited values
    my %CTL_OPT  = set_defaults('server');
    my %defaults = set_defaults('all');
    my %menus    = set_defaults('menus');

    #initialize all display status flags to empty
    while (my $key = each %defaults){
	$CTL_OPT{STAT}{$key} = '';
	$CTL_OPT{$key} = $defaults{$key};
    }

    #add menus to the control options
    $CTL_OPT{menus} = \%menus;

    #--load values from control files
    foreach my $ctlfile (@ctlfiles) {
	open (CTL, "< $ctlfile") or die"ERROR: Could not open the control file \"$ctlfile\".\n";
	
	while (my $line = <CTL>) {
	    chomp($line);
	   
	    if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:]+)\:([^\n\#]*)/) {
		my $key = $1;
		my $value = $2;
		my $stat;

		#remove preceding and trailing whitespace
		$value =~ s/^[\s\t]+|[\s\t]+$//g;

		($value, $stat) = split("=", $value);
		
		if(exists $CTL_OPT{STAT}{$key}){
		    #set value
		    $CTL_OPT{$key} = defined($value) ? $value : '';
		    $CTL_OPT{STAT}{$key} = defined($stat) ? $stat : '';
		}
		elsif (exists $CTL_OPT{$key}) { #should already exist or is a bad value
		    #fix database to not be file location on sqlite
		    if($key eq 'dbname' && defined $value){
			$value =~ /([^\/]+)$/;
			$value = $1;
		    }
		    #set value
		    $CTL_OPT{$key} = defined($value) ? $value : '';
		}
		else {
		    warn "WARNING: Invalid option \'$key\' in control file $ctlfile\n\n";
		}
	    }#now load menus
	    elsif($line =~ /^\#\#Menus from Data::Dumper/){ #only non-sandard format control file
		my $data = join('', <CTL>);
		my $menus; #will be set by data
		eval $data;

		while(my $key = each %$menus){
		    if (exists $CTL_OPT{menus}{$key}) { #should already exist or is a bad value
			#set value
			$CTL_OPT{menus}{$key} = $menus->{$key};
		    }
		    else {
			warn "WARNING: Invalid option \'$key\' in control file $ctlfile\n\n";
		    }
		}
	    }
	}
    }

    #error correct status values
    while(my $key = each %{$CTL_OPT{STAT}}){
	#empty static values are the same as disabled values
	$CTL_OPT{STAT}{$key} = 'DISABLED' if($CTL_OPT{$key} eq '' && $CTL_OPT{STAT}{$key} eq 'STATIC');

	#disabling some options must automatically disable others
	if($key eq 'predictor' && $CTL_OPT{STAT}{$key} eq 'DISABLED'){
	    $CTL_OPT{STAT}{snaphmm} = 'DISABLED';
	    $CTL_OPT{STAT}{gmhmm} = 'DISABLED';
	    $CTL_OPT{STAT}{augustus_species} = 'DISABLED';
	    $CTL_OPT{STAT}{fgenesh_par_file} = 'DISABLED';
	    $CTL_OPT{STAT}{self_train} = 'DISABLED';
	}

	if($key eq 'gmhmm' && $CTL_OPT{STAT}{$key} eq 'DISABLED'){
	    $CTL_OPT{STAT}{self_train} = 'DISABLED';
	}
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
      open (CTL, "< $ctlfile") or die"ERROR: Could not open the control file \"$ctlfile\".\n";
	
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
		  chomp $value;
	       }

	       #require numerical values for certain options
	       if ($CTL_OPT{$key} =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/ &&
		   $value !~  /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/
		  ) {
		  $error .= "ERROR: Invalid setting for the option \'$key\'. The value must be numerical.\n\n"
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
   my @OK = qw(genome
	       protein
	       genome_gff
	       model_gff
	       force
	       predictor
	       retry
	       cpus
	       clean_try
	       again
	       est
	       est_forward
	       single_exon
	       single_length
	       pcov_blastn
	       pid_blastn
	       en_score_limit
	       out_name
	       datastore
	       );

   foreach my $key (@OK){
       $CTL_OPT{$key} = $OPT{$key} if (defined $OPT{$key});
   }

   #check organism type values
   if ($CTL_OPT{organism_type} !~ /^eukaryotic$|^prokaryotic$/) {
      $error .=  "ERROR: organism_type must be set to \'eukaryotic\' or \'prokaryotic\'.\n".
      "The value $CTL_OPT{organism_type} is invalid.\n\n";

      #set the default just to assist in populating remaining errors
      $CTL_OPT{organism_type} = 'eukaryotic';
   }

   #skip repeat masking command line option
   if ($OPT{R} || $CTL_OPT{organism_type} eq 'prokaryotic') {
      $CTL_OPT{model_org} = '';
      $CTL_OPT{repeat_protein} = '';
      $CTL_OPT{rmlib} = '';
      $CTL_OPT{rm_gff} = '';
      $CTL_OPT{rm_pass} = 0;

      print STDERR "INFO: " unless($main::qq);
      print STDERR "No repeats expected in prokaryotic organisms.\n"
	  if($CTL_OPT{organism_type} eq 'prokaryotic' && !$main::qq);
      print STDERR "All repeat masking options will be skipped.\n\n" unless($main::qq);
   }

   #if no repeat masking options are set don't run masking dependent methods
   #i.e. unmasked predictions
   if ($CTL_OPT{model_org} eq '' &&
       $CTL_OPT{repeat_protein} eq '' &&
       $CTL_OPT{rmlib} eq '' &&
       $CTL_OPT{rm_gff} eq '' &&
       ($CTL_OPT{rm_pass} == 0 ||
	$CTL_OPT{genome_gff} eq '')
      ) {
       warn "WARNING: There are no masking options set, yet unmask is set to 1.\n".
	   "This is not valid. The value for unmask will be set to 0.\n\n"
	   if($CTL_OPT{unmask} == 1);

       $CTL_OPT{unmask} = 0;
       $CTL_OPT{_no_mask} = 1; #no masking options found
   }

   #required for evaluator to work
   $CTL_OPT{predictor} = 'model_gff' if($main::eva);
   $CTL_OPT{model_pass} = 1 if($main::eva);
   $CTL_OPT{evaluate} = 1 if($main::eva);

   #evaluator only handles one or the other
   if($main::eva && $CTL_OPT{genome_gff} && $CTL_OPT{model_gff}){
       $error .= "ERROR: In EVALUATOR you can have either models from a MAKER\n".
	         "produced GFF3 file or models from an external GFF3 file, but\n".
		 "not both (Check 'genome_gff' and 'model_gff' options)!\n\n"
   }

   #parse predictor and error check
   $CTL_OPT{predictor} =~ s/\s+//g;
   my @predictors = split(',', $CTL_OPT{predictor});

   $CTL_OPT{_predictor} = {}; #temporary hash
   $CTL_OPT{_run} = {}; #temporary hash
   foreach my $p (@predictors) {
       if ($p !~ /^snap$|^augustus$|^est2genome$|^protein2genome$|^fgenesh$/ &&
	   $p !~ /^genemark$|^jigsaw$|^model_gff$|^pred_gff$/
	   ) {
	   $error .= "FATAL: Invalid predictor defined: $p\n".
	       "Valid entries are: est2genome, model_gff, pred_gff,\n".
	       "snap, genemark, augustus, and fgenesh\n\n";
	   next;
       }
       if($CTL_OPT{organism_type} eq 'prokaryotic' &&
	  $p =~ /^snap$|^augustus$|^fgenesh$|^jigsaw$/
	  ){
	   warn "WARNING: the predictor $p does not support prokaryotic organisms\n".
	       "and will be ignored.\n\n";
	   next;
       }
       if($p =~ /^protein2genome$/ && $CTL_OPT{organism_type} eq 'eukaryotic'){
	   warn "WARNING: protein2genome is currently not supported for eukaryotic organisms\n".
	       "and will be ignored.\n\n";
	   next;
       }

       $CTL_OPT{_predictor}{$p}++;
       $CTL_OPT{_run}{$p}++ unless($p =~ /est2genome|protein2genome|model_gff|pred_gff/);
   }
   $CTL_OPT{_predictor} = [keys %{$CTL_OPT{_predictor}}]; #convert to array
   $CTL_OPT{predictor} = join(",", @{$CTL_OPT{_predictor}}); #reset value for log

   #parse run and error check
   $CTL_OPT{run} =~ s/\s+//g;
   my @run = split(',', $CTL_OPT{run});

   foreach my $p (@run) {
      if ($p !~ /^snap$|^augustus$|^fgenesh$|^genemark$|^jigsaw$/) {
	 $error .= "ERROR: Invalid value defined for run: $p\n".
	 "Valid entries are: snap, augustus, genemark, or fgenesh\n\n";
	 next;
      }
      if($CTL_OPT{organism_type} eq 'prokaryotic' &&
          $p =~ /^snap$|^augustus$|^fgenesh$|^jigsaw$/
	 ){
           warn "WARNING: the predictor $p does not support prokaryotic organisms\n".
               "and will be ignored in option 'run'.\n\n";
           next;
      }
      $CTL_OPT{_run}{$p}++;
   }
   $CTL_OPT{_run} = [keys %{$CTL_OPT{_run}}]; #convert to array
   $CTL_OPT{run} = join(",", @{$CTL_OPT{_run}}); #reset value for log

   #check blast type validity and related values (NCBI vs. WUBLAST)
   if ($CTL_OPT{blast_type} !~ /^wublast$|^ncbi$/) {
      warn "WARNING: blast_type must be set to \'wublast\' or \'ncbi\'.\n",
      "The value $CTL_OPT{blast_type} is invalid.\n",
      "This will now be reset to the default 'wublast'.\n\n";
      
      $CTL_OPT{blast_type} = 'wublast';
   }
   
   if ($CTL_OPT{blast_type} =~ /^wublast$/ &&
       ! -f $CTL_OPT{blastn} &&
       ! -f $CTL_OPT{blastx} &&
       ! -f $CTL_OPT{tblastx} &&
       -f $CTL_OPT{blastall}
      ) {
      warn "WARNING: blast_type is set to \'wublast\' but wublast executables\n",
      "can not be located.  NCBI blast will be used instead.\n\n";

      $CTL_OPT{blast_type} = 'ncbi';
   }

   if ($CTL_OPT{blast_type} =~ /^ncbi$/ &&
       ! -f $CTL_OPT{blastall} &&
       -f $CTL_OPT{blastn} &&
       -f $CTL_OPT{blastx} &&
       -f $CTL_OPT{tblastx}
      ) {
      warn "WARNING: blast_type is set to \'ncbi\' but ncbi executables\n",
      "can not be located.  WUBLAST blast will be used instead.\n\n";

      $CTL_OPT{blast_type} = 'wublast';
   }
   
   #use standard value to refer to both NCBI and WUBLAST
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
   push (@infiles, 'est') if($CTL_OPT{est}); 
   push (@infiles, 'protein') if($CTL_OPT{protein}); 
   push (@infiles, 'altest') if($CTL_OPT{altest}); 
   push (@infiles, 'est_reads') if($CTL_OPT{est_reads}); 
   push (@infiles, 'probuild') if (grep {/genemark/} @{$CTL_OPT{_run}});
   push (@infiles, 'fathom') if ($CTL_OPT{enable_fathom} && $CTL_OPT{evaluate});
   push (@infiles, 'est_gff') if($CTL_OPT{est_gff});
   push (@infiles, 'protein_gff') if($CTL_OPT{protein_gff});
   push (@infiles, 'genome_gff') if($CTL_OPT{genome_gff});
   push (@infiles, 'pred_gff') if($CTL_OPT{pred_gff});
   push (@infiles, 'model_gff') if ($CTL_OPT{model_gff});
   push (@infiles, 'snap') if (grep {/snap/} @{$CTL_OPT{_run}});
   push (@infiles, 'augustus') if (grep {/augustus/} @{$CTL_OPT{_run}}); 
   push (@infiles, 'fgenesh') if (grep {/fgenesh/} @{$CTL_OPT{_run}});
   push (@infiles, 'twinscan') if (grep {/twinscan/} @{$CTL_OPT{_run}});
   push (@infiles, 'jigsaw') if (grep {/jigsaw/} @{$CTL_OPT{_run}});
   push (@infiles, 'repeat_protein') if ($CTL_OPT{repeat_protein});
   push (@infiles, 'RepeatMasker') if($CTL_OPT{rmlib});
   push (@infiles, 'RepeatMasker') if($CTL_OPT{model_org});
   push (@infiles, 'rm_gff') if($CTL_OPT{rm_gff});
   push (@infiles, 'rmlib') if ($CTL_OPT{rmlib});
   
   if($CTL_OPT{organism_type} eq 'eukaryotic'){
       push (@infiles, 'exonerate') if($CTL_OPT{est}); 
       push (@infiles, 'exonerate') if($CTL_OPT{protein});
       push (@infiles, 'gmhmme3') if (grep {/genemark/} @{$CTL_OPT{_run}});
   }
   elsif($CTL_OPT{organism_type} eq 'prokaryotic'){
       push (@infiles, 'gmhmmp') if (grep {/genemark/} @{$CTL_OPT{_run}});
   }

   #uniq the array
   my %uniq;
   @uniq{@infiles} = ();
   @infiles = keys %uniq;

   #verify existence of required values
   foreach my $in (@infiles) {
      if (not $CTL_OPT{$in}) {
	 $error .= "ERROR: You have failed to provide a value for \'$in\' in the control files.\n\n";
	 next;
      }
      elsif ((my @files = split(/\,/, $CTL_OPT{$in})) > 1){#handle comma seperated list
	  my %uniq; #make files uniq
	  @files = grep {! $uniq{$_}++} @files;
	  my @non = grep {! -f $_} @files;
	  $error .= "ERROR: The \'$in\' files ".join(', ', @files)." do not exist.\n".
	      "Please check settings in the control files.\n\n"if(@non);

	  @files = map {Cwd::abs_path($_)} @files unless ($in =~ /^_blastn$|^_blastx$|^_tblastx$|^_formater$/);
	  $CTL_OPT{$in} = join(',', @files); #order entries for logging (lets run log work as is)
      }
      elsif (! -f $CTL_OPT{$in}) {
	 $error .= "ERROR: The \'$in\' file $CTL_OPT{$in} does not exist.\n".
	 "Please check settings in the control files.\n\n";
	 next;
      }
      else{
	  #set the absolute path to the file to reduce ambiguity
	  #this breaks blast which requires symbolic links
	  $CTL_OPT{$in} = Cwd::abs_path($CTL_OPT{$in}) unless ($in =~ /^_blastn$|^_blastx$|^_tblastx$|^_formater$/);
      }
   }

   #--error check sometimes required values
   if ((grep {/model_gff/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{model_gff} &&
       (!$CTL_OPT{genome_gff} || !$CTL_OPT{model_pass})
      ){
       $error .= "ERROR: You must provide gene models in a GFF3 file to use model_gff as a predictor.\n\n";
   }
   if ((grep {/pred_gff/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{pred_gff} &&
       (!$CTL_OPT{genome_gff} || !$CTL_OPT{pred_pass})
      ){
       $error .= "ERROR: You must provide gene predictions in a GFF3 file to use pred_gff as a predictor.\n\n";
   }
   if ((grep {/est2genome/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{est} &&
       !$CTL_OPT{altest} &&
       !$CTL_OPT{est_gff} &&
       !$CTL_OPT{altest_gff} &&
       (!$CTL_OPT{est_pass} || !$CTL_OPT{genome_gff})
      ){
       $error .= "ERROR: You must provide some form of EST evidence to use est2genome as a predictor.\n\n";
   } 
   if ((grep {/protein2genome/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{protein} &&
       !$CTL_OPT{protein_gff} &&
       (!$CTL_OPT{protein_pass} || !$CTL_OPT{genome_gff})
       ){
       $error .= "ERROR: You must provide some form of protein evidence to use protein2genome as a predictor.\n\n";
   }
			       
   #--error check that values are meaningful
   if ((grep {/^augustus$/} @infiles) && not $CTL_OPT{augustus_species}) {
       $error .= "ERROR: There is no species specified for Augustus (augustus_species).\n\n";
   }
   if ((grep {/^augustus$/} @infiles) &&
       (! $ENV{AUGUSTUS_CONFIG_PATH} || ! -f "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg")
      ) {
      $error .= "ERROR: The environmental variable AUGUSTUS_CONFIG_PATH has not been set\n".
      "or is not set correctly Please set this in your profile per Augustus\n".
      "installation instructions then try again.\n\n";
   }
   if ((grep {/^snaps$|^fathom$/} @infiles) && not $CTL_OPT{snaphmm}) {
       $error .= "ERROR: There is no HMM specified for for SNAP/Fathom (snaphmm).\n\n";
   }
   elsif ((grep {/^snap$|^fathom$/} @infiles) && ! -f $CTL_OPT{snaphmm} &&
       (! exists $ENV{ZOE} || ! -f $ENV{ZOE}."/HMM/".$CTL_OPT{snaphmm})
      ) {
      $error .= "ERROR: The snaphmm specified for SNAP/Fathom does not exist.\n\n";
   }
   if ((grep {/^gmhmme3$/} @infiles) && not $CTL_OPT{gmhmm}) {
      $ error .=  "ERROR: There is no HMM specified for for GeneMark (gmhmm).\n\n";
   }
   elsif ((grep {/^gmhmme3$/} @infiles) && ! -f $CTL_OPT{gmhmm}) {
      $error .= "ERROR: The HMM specified for GeneMark does not exist.\n\n";
   }
   if (grep {/^fgenesh$/} @infiles) {
      if (! $CTL_OPT{fgenesh_par_file}) {
	  $error .= "ERROR: There is no parameter file secified for FgenesH (fgenesh_par_file)\n\n";
      }
      elsif (! -f $CTL_OPT{fgenesh_par_file}) {
	 $error .= "ERROR: The parameter file specified for fgenesh does not exist.\n\n";
      }
   }
   if ($CTL_OPT{max_dna_len} < 50000) {
      warn "WARNING: \'max_dna_len\' is set too low.  The minimum value permited is 50,000.\n".
      "max_dna_len will be reset to 50,000\n\n";
      $CTL_OPT{max_dna_len} = 50000;
   }
   if ($CTL_OPT{split_hit} > 0 && $CTL_OPT{organism_type} eq 'prokaryotic') {
      $CTL_OPT{split_hit} = 0;
   }
   if ($CTL_OPT{single_exon} == 0 && $CTL_OPT{organism_type} eq 'prokaryotic') {
      warn "WARNING: \'single_exon\' is required for prokaryotic genomes and will be set to 1.\n\n" unless($main::qq);
      $CTL_OPT{single_exon} = 1;
   }
   if ($CTL_OPT{min_contig} <= 0) {
      warn "WARNING: \'min_contig\' must be set to 1 or higher.\n".
      "min_contig will be reset to 1\n\n";
      $CTL_OPT{min_contig} = 1;
   }
   if ($CTL_OPT{min_contig} < 0) {
      warn "WARNING: \'min_protein\' must be set to 0 or higher.\n".
      "min_protein will be reset to 0\n\n";
      $CTL_OPT{min_protein} = 0;
   }
   if ($CTL_OPT{AED_threshold} < 0 || $CTL_OPT{AED_threshold} > 1) {
      warn "WARNING: \'AED_threshold\' must be set to a value betweeb 0 and 1.\n".
      "AED_threshold will be reset to 1\n\n";
      $CTL_OPT{AED_threshold} = 1;
   }
   if ($CTL_OPT{retry} < 0) {
      warn "WARNING: \'retry\' must be set to 0 or greater.\n".
	   "It will now be set to 0\n\n";
      $CTL_OPT{retry} = 0;
   } 
   if($CTL_OPT{TMP} && ! -d $CTL_OPT{TMP}){
       $error .= "The TMP value \'$CTL_OPT{TMP}\' is not a directory or does not exist.\n\n";
   }
   if($main::eva && $CTL_OPT{genome_gff} && $CTL_OPT{model_gff}){ #only for evaluator
       $error .= "You can only specify a GFF3 file for genome_gff or model_gff no both!!\n\n";
   }
   
   die $error if ($error);   

   #--check genome fasta file
   my $fasta_gff = ($CTL_OPT{genome_gff}) ? $CTL_OPT{genome_gff} : $CTL_OPT{model_gff};
   my $iterator = new Iterator::Any( -fasta => $CTL_OPT{genome},
				     -gff => $fasta_gff
				   );

   if ($iterator->number_of_entries() == 0) {
      my $genome = (! $CTL_OPT{genome}) ? $fasta_gff : $CTL_OPT{genome};
      die "ERROR:  The file $genome contains no fasta entries\n\n";
   }

   #--decide whether to force datastore, datastore will already be defined if selected by user 
   if(! defined $CTL_OPT{datastore}){
       if($iterator->number_of_entries() > 1000) {
	   warn "WARNING:  There are more than 1000 fasta entries in the input file.\n".
	       "A two depth datastore will be used to avoid overloading the data structure of\n".
	       "the output directory.\n\n" unless($main::qq);
	   
	   $CTL_OPT{datastore} = 1;
       }
       else{
	   $CTL_OPT{datastore} = 0;
       }
   }

   #--decide if gff database should be created
   my @gffs = grep {/\_gff$/} @{[keys %CTL_OPT]};
   foreach my $key (@gffs) {
      if ($CTL_OPT{$key}) {
	 $CTL_OPT{go_gffdb} = 1;
	 last;
      }
   }

   #--check validity of the alternate peptide
   $CTL_OPT{alt_peptide} = uc($CTL_OPT{alt_peptide});
   if ($CTL_OPT{alt_peptide} !~ /^[ACDEFGHIKLMNPQRSTVWXY]$/) {
      warn "WARNING: Invalid alternate peptide \'$CTL_OPT{alt_peptide}\'.\n",
      "This will be set to the default 'C'.\n\n";
      $CTL_OPT{alt_peptide} = 'C';
   }

   #--set values for datastructure
   my $genome = $CTL_OPT{genome};
   $genome = $CTL_OPT{genome_gff} if (not $genome);
   $CTL_OPT{CWD} = Cwd::cwd();
   if(! $CTL_OPT{out_name}){
       ($CTL_OPT{out_name}) = $genome =~ /([^\/]+)$/;
       $CTL_OPT{out_name} =~ s/\.[^\.]+$//;
   }
   if(! $CTL_OPT{out_base}){
      $CTL_OPT{out_base} = $CTL_OPT{CWD}."/$CTL_OPT{out_name}.maker.output";
   }
   mkdir($CTL_OPT{out_base}) if(! -d $CTL_OPT{out_base});
   die "ERROR: Could not build output directory $CTL_OPT{out_base}\n"
        if(! -d $CTL_OPT{out_base});

   #shared database for logging encounters of IDs or other values
   $CTL_OPT{SEEN_file} = "$CTL_OPT{out_base}/seen.dbm";   
   
   #--set up optional global TMP (If TMP is not accessible from other nodes
   #they will default back to /tmp)
   if($CTL_OPT{TMP}){
       print STDERR "\nTMP_STAT: User specified TMP is $CTL_OPT{TMP}: PID=$$-root\n" if($main::dtmp); ##temp
       $CTL_OPT{_TMP} = tempdir("maker_XXXXXX", CLEANUP => 1, DIR => $CTL_OPT{TMP});
   }
   else{
       $CTL_OPT{_TMP} = $TMP;
   }

   set_global_temp($CTL_OPT{_TMP});

   #--exit with status of 0 if just checking control files with -check flag
   exit(0) if($OPT{check});

   #--set an initialization lock so steps are locked to a single process
   my $i_lock; #init lock, it is only a temporary blocking lock
   unless($i_lock = new File::NFSLock($CTL_OPT{out_base}."/.init_lock", 'EX', 40, 45)){
       die "ERROR: Cannot get initialization lock.\n\n";
   }

   #--check if another instance of maker is running, and lock the directory
   #lock must be global or it will be destroyed outside of block
   if($LOCK = new File::NFSLock($CTL_OPT{out_base}."/.gi_lock", 'SH', 40, 40)){
       $LOCK->maintain(30);
   }
   else{
       die "ERROR: The directory is locked.  Perhaps by an instance of MAKER or EVALUATOR.\n\n";
   }

   #check who else is also sharing the lock and if running same settings
   my $app = ($main::eva) ? "eval" : "maker";

   if($LOCK->owners() == 1){ #I am only holder of the lock
       #log the control files
       generate_control_files($CTL_OPT{out_base}, 'all', \%CTL_OPT, 1);
       unlink($CTL_OPT{SEEN_file}) if (-f $CTL_OPT{SEEN_file});
   }
   else{
       #compare current control files to logged files
       my @ctl_logs = ($CTL_OPT{out_base}."/$app\_opts.log",
		       $CTL_OPT{out_base}."/$app\_bopts.log",
		       $CTL_OPT{out_base}."/$app\_exe.log"
		       );

       my @ctl_news = ($TMP."/$app\_opts.log",
		       $TMP."/$app\_bopts.log",
		       $TMP."/$app\_exe.log"
		       );

       unless (-f $ctl_logs[0] && -f $ctl_logs[1] && -f $ctl_logs[2]){
	   die "ERROR: Could not query control option logs\n\n";
       }

       #log the current control files for comparison
       generate_control_files($TMP, 'all', \%CTL_OPT, 1);

       my $log_data;
       my $new_data;
       foreach my $ctl (@ctl_logs){
	   open(my $FH, "< $ctl");
	   $log_data .= join('', <$FH>);
	   close($FH);
       }

       foreach my $ctl (@ctl_news){
	   open(my $FH, "< $ctl");
	   $new_data .= join('', <$FH>);
	   close($FH);
       }

       #should be exactly identical
       if($log_data ne $new_data){
	   die "ERROR: Cannot start process. MAKER/EVALUATOR already running\n".
	       "with different settings in this same directory.\n\n";
       }
       else{#start a second MAKER process, but give a warning
	   warn "WARNING: Multiple MAKER processes have been started in the\n".
	        "same directory.\n\n";

	   $CTL_OPT{_multi_chpc}++;
       }
   }

   $i_lock->unlock; #release init lock

   #---set up blast databases and indexes for analyisis
   create_blastdb(\%CTL_OPT, $mpi_size);
   build_all_indexes(\%CTL_OPT);

   return %CTL_OPT;
}
#-----------------------------------------------------------------------------
#this function generates generic control files
sub generate_control_files {
   my $dir = shift @_ || Cwd::cwd();
   my $type = shift @_ || 'all';
   my %O = set_defaults($type, shift @_);
   my $log = shift;
   my $ev = 1 if($main::eva);

   my $app = ($ev) ? "eval" : "maker"; #extension
   my $ext = ($log) ? "log" : "ctl"; #extension

   if ($type !~ /^all$|^opts$|^bopts$|^exe$|^menus$|^server$/) {
       warn "WARNING: Invalid type \'$type\' in GI::generate_control_files";
       $type = 'all';
   }

   #--build opts.ctl file
   if($type eq 'all' || $type eq 'opts'){
       open (OUT, "> $dir/$app\_opts.$ext");
       print OUT "#-----Genome (Required for De-Novo Annotation)\n" if(!$ev);
       print OUT "#-----Genome (Required if not internal to GFF3 file)\n" if($ev);
       print OUT "genome:$O{genome} #genome sequence file in fasta format\n";
       print OUT "\n";
       print OUT "#-----Re-annotation Options (Only MAKER derived GFF3)\n" if(!$ev);
       print OUT "#-----MAKER Derived GFF3 Annotations to Evaluate (genome fasta is internal to GFF3)\n" if($ev);
       print OUT "genome_gff:$O{genome_gff} #re-annotate genome based on this gff3 file\n" if(!$ev);
       print OUT "genome_gff:$O{genome_gff} #MAKER derived gff3 file\n" if($ev);
       print OUT "est_pass:$O{est_pass} #use ests in genome_gff: 1 = yes, 0 = no\n";
       print OUT "altest_pass:$O{altest_pass} #use alternate organism ests in genome_gff: 1 = yes, 0 = no\n";
       print OUT "protein_pass:$O{protein_pass} #use proteins in genome_gff: 1 = yes, 0 = no\n";
       print OUT "rm_pass:$O{rm_pass} #use repeats in genome_gff: 1 = yes, 0 = no\n";
       print OUT "model_pass:$O{model_pass} #use gene models in genome_gff: 1 = yes, 0 = no\n" if(!$ev);
       print OUT "pred_pass:$O{pred_pass} #use ab-initio predictions in genome_gff: 1 = yes, 0 = no\n";
       print OUT "other_pass:$O{other_pass} #passthrough everything else in genome_gff: 1 = yes, 0 = no\n" if(!$ev);
       print OUT "\n";
       print OUT "#-----External GFF3 Annotations to Evaluate\n" if($ev);
       print OUT "model_gff:$O{model_gff} #gene models from an external gff3 file\n" if($ev);
       print OUT "\n"if($ev);
       print OUT "#-----EST Evidence (you should provide a value for at least one)\n";
       print OUT "est:$O{est} #non-redundant set of assembled ESTs in fasta format (classic EST analysis)\n";
       print OUT "est_reads:$O{est_reads} #unassembled nextgen mRNASeq in fasta format (not fully implemented)\n";
       print OUT "altest:$O{altest} #EST/cDNA sequence file in fasta format from an alternate organism\n";
       print OUT "est_gff:$O{est_gff} #EST evidence from an external gff3 file\n";
       print OUT "altest_gff:$O{altest_gff} #Alternate organism EST evidence from a seperate gff3 file\n";
       print OUT "\n";
       print OUT "#-----Protein Homology Evidence (you should provide a value for at least one)\n";
       print OUT "protein:$O{protein}  #protein sequence file in fasta format\n";
       print OUT "protein_gff:$O{protein_gff}  #protein homology evidence from an external gff3 file\n";
       print OUT "\n";
       print OUT "#-----Repeat Masking (leave values blank to skip)\n";
       print OUT "model_org:$O{model_org} #model organism for RepBase masking in RepeatMasker\n";
       print OUT "repeat_protein:$O{repeat_protein} #a database of transposable element proteins in fasta format\n";
       print OUT "rmlib:$O{rmlib} #an organism specific repeat library in fasta format\n";
       print OUT "rm_gff:$O{rm_gff} #repeat elements from an external gff3 file\n";
       print OUT "\n";
       print OUT "#-----Gene Prediction Options\n" if(!$ev);
       print OUT "#-----EVALUATOR Ab-Initio Comparison Options\n" if($ev);
       print OUT "organism_type:$O{organism_type} #eukaryotic or prokaryotic. Default is eukaryotic\n";
       print OUT "run:$O{run} #ab-initio methods to run (seperate multiple values by ',')\n" if($ev);
       print OUT "predictor:$O{predictor} #prediction methods for annotations (seperate multiple values by ',')\n" if(!$ev);
       print OUT "unmask:$O{unmask} #Also run ab-initio methods on unmasked sequence, 1 = yes, 0 = no\n";
       print OUT "snaphmm:$O{snaphmm} #SNAP HMM model\n";
       print OUT "gmhmm:$O{gmhmm} #GeneMark HMM model\n";
       print OUT "augustus_species:$O{augustus_species} #Augustus gene prediction model\n";
       print OUT "fgenesh_par_file:$O{fgenesh_par_file} #Fgenesh parameter file\n";
       print OUT "model_gff:$O{model_gff} #gene models from an external gff3 file (annotation pass-through)\n" if(!$ev);
       print OUT "pred_gff:$O{pred_gff} #ab-initio predictions from an external gff3 file\n";
       print OUT "\n";
       print OUT "#-----Other Annotation Type Options (features maker doesn't recognize)\n" if(!$ev);
       print OUT "other_gff:$O{other_gff} #features to pass-through to final output from an extenal gff3 file\n" if(!$ev);
       print OUT "\n" if(!$ev);
       print OUT "#-----External Application Specific Options\n";
       print OUT "alt_peptide:$O{alt_peptide} #amino acid used to replace non standard amino acids in blast databases\n";
       print OUT "cpus:$O{cpus} #max number of cpus to use in BLAST and RepeatMasker\n";
       print OUT "\n";
       print OUT "#-----MAKER Specific Options\n";
       print OUT "evaluate:$O{evaluate} #run EVALUATOR on all annotations, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "max_dna_len:$O{max_dna_len} #length for dividing up contigs into chunks (larger values increase memory usage)\n";
       print OUT "min_contig:$O{min_contig} #all contigs from the input genome file below this size will be skipped\n" if(!$ev);
       print OUT "min_protein:$O{min_protein} #all gene annotations must produce a protein of at least this many amino acids in length\n" if(!$ev);
       print OUT "AED_threshold:$O{AED_threshold} #Maximum Annotation Edit Distance allowed for annotations (bound by 0 and 1)\n" if(!$ev);
       print OUT "softmask:$O{softmask} #use soft-masked rather than hard-masked seg filtering for wublast\n";
       print OUT "split_hit:$O{split_hit} #length for the splitting of hits (expected max intron size for evidence alignments)\n";
       print OUT "pred_flank:$O{pred_flank} #length of sequence surrounding EST and protein evidence used to extend gene predictions\n";
       print OUT "single_exon:$O{single_exon} #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no\n";
       print OUT "single_length:$O{single_length} #min length required for single exon ESTs if \'single_exon\ is enabled'\n";
       print OUT "keep_preds:$O{keep_preds} #Add non-overlapping ab-inito gene prediction to final annotation set, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "map_forward:$O{map_forward} #try to map names and attributes forward from gff3 annotations, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "retry:$O{retry} #number of times to retry a contig if there is a failure for some reason\n";
       print OUT "clean_try:$O{clean_try} #removeall data from previous run before retrying, 1 = yes, 0 = no\n";
       print OUT "clean_up:$O{clean_up} #removes theVoid directory with individual analysis files, 1 = yes, 0 = no\n";
       print OUT "TMP:$O{TMP} #specify a directory other than the system default temporary directory for temporary files\n";
       print OUT "\n";
       print OUT "#-----EVALUATOR Control Options\n";
       print OUT "side_thre:$O{side_thre}\n";
       print OUT "eva_window_size:$O{eva_window_size}\n";
       print OUT "eva_split_hit:$O{eva_split_hit}\n";
       print OUT "eva_hspmax:$O{eva_hspmax}\n";
       print OUT "eva_gspmax:$O{eva_gspmax}\n";
       print OUT "enable_fathom:$O{enable_fathom}\n";
       close (OUT);
   }
    
   #--build bopts.ctl file
   if($type eq 'all' || $type eq 'bopts'){
       open (OUT, "> $dir/$app\_bopts.$ext");
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
       print OUT "eva_pcov_blastn:$O{eva_pcov_blastn} #EVALUATOR Blastn Percent Coverage Threshold EST-Genome Alignments\n";
       print OUT "eva_pid_blastn:$O{eva_pid_blastn} #EVALUATOR Blastn Percent Identity Threshold EST-Genome Alignments\n";
       print OUT "eva_eval_blastn:$O{eva_eval_blastn} #EVALUATOR Blastn eval cutoff\n";
       print OUT "eva_bit_blastn:$O{eva_bit_blastn} #EVALUATOR Blastn bit cutoff\n";
       print OUT "\n";
       print OUT "ep_score_limit:$O{ep_score_limit} #Exonerate protein percent of maximal score threshold\n";
       print OUT "en_score_limit:$O{en_score_limit} #Exonerate nucleotide percent of maximal score threshold\n";
       close(OUT);
   }

   #--build maker_exe.ctl file
   if($type eq 'all' || $type eq 'exe'){
       open (OUT, "> $dir/$app\_exe.$ext");
       print OUT "#-----Location of Executables Used by MAKER/EVALUATOR\n";
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
       print OUT "gmhmme3:$O{gmhmme3} #location of eukaryotic genemark executable\n";
       print OUT "gmhmmp:$O{gmhmmp} #location of prokaryotic genemark executable\n";
       print OUT "augustus:$O{augustus} #location of augustus executable\n";
       print OUT "fgenesh:$O{fgenesh} #location of fgenesh executable\n";
       print OUT "twinscan:$O{twinscan} #location of twinscan executable (not yet implemented)\n";
       print OUT "\n";
       print OUT "#-----Other Algorithms\n";
       print OUT "jigsaw:$O{jigsaw} #location of jigsaw executable (not yet implemented)\n";
       print OUT "qrna:$O{qrna} #location of qrna executable (not yet implemented)\n";
       print OUT "fathom:$O{fathom} #location of fathom executable (not yet implemented)\n";
       print OUT "probuild:$O{probuild} #location of probuild executable (required for genemark)\n";
       close(OUT);    
   }

   #--build server.ctl file
   if($type eq 'server'){
       open (OUT, "> $dir/server.$ext");
       print OUT "#-----Database Setup\n";
       print OUT "DBI:$O{DBI} #interface type to database\n";
       print OUT "dbname:$O{dbname} #database name\n";
       print OUT "host:$O{host} #host on which database is found\n";
       print OUT "port:$O{port} #port on host to access database\n";
       print OUT "username:$O{username} #username to connect to database\n";
       print OUT "password:$O{password} #password to connect to database\n";
       print OUT "\n";
       print OUT "#-----Communication Options\n";
       print OUT "admin_email:$O{admin_email} #Address for sending error and status information\n";
       print OUT "smtp_server:$O{smtp_server} #Outgoing e-mail server\n";
       print OUT "\n";
       print OUT "#-----Web Setup\n";
       print OUT "apache_user:$O{apache_user} #username apache runs as\n";
       print OUT "web_address:$O{web_address} #base web address to server hosting MWAS\n";
       print OUT "cgi_dir:$O{cgi_dir} #web accesible directory to house MWAS CGI content\n";
       print OUT "cgi_web:$O{cgi_web} #url to cgi_dir above (can be relative)\n";
       print OUT "html_dir:$O{html_dir} #web accesible directory to house MWAS HTML conent\n";
       print OUT "html_web:$O{html_web} #url to html_dir (can be relative)\n";
       print OUT "data_dir:$O{data_dir} #directory for saving user uploaded files, running jobs, and storing results\n";
       print OUT "font_file:$O{font_file} #font file for webpage CAPTCHA\n";
       print OUT "\n";
       print OUT "#-----External Viewer Setup\n";
       print OUT "soba_url:$O{soba_url} #url to Sequence Ontology SOBA CGI script\n";
       print OUT "APOLLO_ROOT:$O{APOLLO_ROOT} #base directory for Apollo installation.\n";
       print OUT "JBROWSE_ROOT:$O{JBROWSE_ROOT} #base directory for JBrowse installation.\n";
       print OUT "GBROWSE_MASTER:$O{GBROWSE_MASTER} #path to GBrowse.conf file.\n";
       print OUT "\n";
       print OUT "#-----MAKER Server Specific Options\n";
       print OUT "use_login:$O{use_login} #whether to require login to access the web interface, 1 = yes, 0 = no\n";
       print OUT "allow_guest:$O{allow_guest} #enable guest accounts on the server, 1 = yes, 0 = no\n";
       print OUT "allow_register:$O{allow_register} #allow users to register themselves, 1 = yes, 0 = no\n";
       print OUT "tutorials:$O{tutorials} #show example data on \"New Job\" screen, 1 = yes, 0 = no\n";
       print OUT "max_cpus:$O{max_cpus} #maximum number of cpus that can be dedicated to all MAKER jobs\n";
       print OUT "job_cpus:$O{job_cpus} #maximum number of cpus that can be used by a single MAKER job\n";
       print OUT "max_submit_user:$O{max_submit_user} #maximum submission size for registered users (0 = no limit)\n";
       print OUT "max_submit_guest:$O{max_submit_guest} #maximum submission size for guest users (0 = no limit)\n";
       print OUT "persist_user:$O{persist_user} #time results persist for registered users, in hours (0 = no limit)\n";
       print OUT "persist_guest:$O{persist_guest} #time results persist for guest users, in hours (0 = no limit)\n";
       print OUT "inactive_user:$O{inactive_user} #time user account can be inactive, in days (0 = no limit)\n";
       print OUT "inactive_guest:$O{inactive_guest} #time guest account can be inactive, in days (0 = no limit)\n";
       print OUT "MPI:$O{MPI} #use mpi_maker instead of maker\n";
       print OUT "mpiexec:$O{mpiexec} #mpiexec command line for running MPI\n";

       close(OUT);    
   }

   #--build menus.ctl file
   if($type eq 'menus'){
       open (OUT, "> $dir/menus.$ext");
       print OUT "##Menus from Data::Dumper\n";
       print OUT Data::Dumper->Dump([\%O], [qw(menus)]);
       close(OUT);    
   }
}

1;
