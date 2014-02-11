#------------------------------------------------------------------------
#----                                GI                              ----
#------------------------------------------------------------------------
package GI;

use strict;
use vars qw(@ISA @EXPORT $VERSION $TMP $RANK $LOCK $HOST);
use FindBin;
use Exporter;
use FileHandle;
use File::Temp qw(tempfile tempdir);
use File::Path;
use File::Copy;
use Carp;
use Cwd qw(cwd);
use File::NFSLock;
use File::Which;
use File::Glob;
use Dumper::GFF::GFFV3;
use Dumper::XML::Game;
use Error qw(:try);
use Error::Simple;
use URI::Escape;
use Data::Dumper;
use Getopt::Long;
use FileHandle;
use PostData;
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
use Widget::trnascan;
use Widget::snoscan;
use Widget::formater;
use PhatHit_utils;
use Shadower;
use polisher::exonerate::protein;
use polisher::exonerate::est;
use polisher::exonerate::altest;
use maker::auto_annotator;
use cluster;
use repeat_mask_seq;
use maker::sens_spec;
use Digest::MD5;
use FastaDB;
use List::Util;
use Sys::Hostname;

@ISA = qw(
	);
$RANK = 0;
$TMP = tempdir("maker_XXXXXX", CLEANUP => 1, TMPDIR => 1);
$HOST = Sys::Hostname::hostname();
print STDERR "TMP_STAT: TMP is being initialized to $TMP: PID=$$\n" if($main::dtmp);
#------------------------------------------------------------------------
#--------------------------- CLASS FUNCTIONS ----------------------------
#------------------------------------------------------------------------
sub version {
    $VERSION = '2.31';
    return $VERSION;
}
#------------------------------------------------------------------------
sub new_instance_temp {
    my $base = shift;

    if($base && -d $base){
	my $dir = tempdir("maker_XXXXXX", CLEANUP => 1, DIR => $base);
	set_global_temp($dir);
    }

    return get_global_temp();
}
#------------------------------------------------------------------------
{
my %mounts; #persistent list of mounts to avoid calling df so often
my $host = Sys::Hostname::hostname();
sub mount_check {
    my $path = shift;

    return undef if(! -e $path);

    $path = Cwd::abs_path($path);

    #it's a file and not a directory
    my $f_name;
    if(-f $path){
	($path, $f_name) = $path =~ /(.*)\/([^\/]+)$/;
    }
    if($mounts{$path}){
	return length($f_name) ? "$mounts{$path}/$f_name" : $mounts{$path};
    }

    #get only line with mount point
    (my $p = $path) =~ s/([^A-Za-z_0-9\.\-\/\~])/\\$1/g;
    open(DF, "df -P $p |");
    my @F;
    while(my $line = <DF>){
        chomp $line;
        @F = split(/\s+/, $line);
        last if grep {/\%/} @F;
    }
    close(DF);

    #accept mount points with spaces
    my $mounted_on = '';
    my $filesystem = '';
    for(my $i = $#F; $i >= 0; $i--){
        if($F[$i] =~ /^\d+\%$/ &&
           -e $mounted_on &&
	   Cwd::abs_path($mounted_on) eq $mounted_on
          ){
            $filesystem = join(' ', @F[0..$i-4]);
            last;
        }

        $mounted_on = ($mounted_on) ? "$F[$i] $mounted_on" : $F[$i];
    }
    $mounted_on = Cwd::abs_path($mounted_on);

    #build new mount point
    if($filesystem =~ /^fhgfs/){
	$mounts{$path} = "$filesystem:$path";
    }
    elsif($filesystem =~ /\:/){
        $mounted_on = quotemeta($mounted_on);
	my $mpath = $path;
        $mpath =~ s/^$mounted_on//;
        $mpath = "$filesystem/$mpath";
        $mpath =~ s/\/\/+/\//g;
	$mounts{$path} = $mpath;
    }
    else{
	$mounts{$path} = "$host:$path";
    }

    return length($f_name) ? "$mounts{$path}/$f_name" : $mounts{$path};
}
#------------------------------------------------------------------------
sub is_NFS_mount {
    my $path = shift;

    $path = mount_check($path);
    return undef if(! $path);
    my $safe = quotemeta($host);
    my $is_NFS = ($path =~ /^$safe\:/) ? 0 : 1;
    $is_NFS = 2 if($path =~ /^fhgfs/ && $path !~ /\:/);

    return $is_NFS;
}
}
#------------------------------------------------------------------------
sub set_global_temp {
    my $dir = shift;
    
    return if(! $dir || ! -d $dir);

    #remove old tempdir if user supplied a new one
    if(Cwd::abs_path($TMP) ne Cwd::abs_path($dir)){
	File::Path::rmtree($TMP);
    }

    $TMP = $dir;
    mkdir "$TMP/$RANK" unless(-d "$TMP/$RANK"); #for shared MPI
}
#------------------------------------------------------------------------
sub get_global_temp {
    
    return $TMP;
}
#------------------------------------------------------------------------
sub RANK {
    if(@_){
	$RANK = shift @_;
	mkdir "$TMP/$RANK" unless(-d "$TMP/$RANK"); #for shared MPI
    }
    return $RANK;
}
#------------------------------------------------------------------------
sub LOCK {
    return $LOCK;
}
#------------------------------------------------------------------------
sub s_abs_path {
    my $path = shift;

    if(-l $path){
	my ($dir, $file) = $path =~ /(.*\/)?([^\/]+)$/;
	$dir = '.' if(!$dir);
	$path = Cwd::abs_path($dir)."/$file";
	$path =~ s/\/+/\//g;
	return $path;
    }
    else{
	return Cwd::abs_path($path);
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
   my $q_seq_obj = shift @_;
   my $break = shift @_;
   my $def = shift @_;
   my $s_length = shift @_;
   my $index_files = shift @_;
   my $blast_keepers = shift @_;
   my $blast_holdovers = shift @_;
   my $the_void = shift @_;
   my %CTL_OPT = %{shift @_};
   my $type = shift @_; #blastn, blastx, or tblastx
   my $LOG = shift @_;

   print STDERR "merging blast reports...\n" unless($main::quiet);

   my %sets;
   foreach my $h (@$blast_keepers, @$blast_holdovers){
       push (@{$sets{$h->name}}, $h);
   }

   my $low;
   my $high;
   my @to_reblast;
   my $no_merge = [];
   while(my $key = each %sets){
       if(@{$sets{$key}} == 1){
           push(@$no_merge, @{$sets{$key}}) unless($sets{$key}[0]{_remove});
           next;
       }
       
       my $clusters = SimpleCluster::cluster_hits($sets{$key}, $CTL_OPT{split_hit}, 1); #strand ignore flag set
       
       foreach my $c (@$clusters){
           if(@$c == 1){
               push(@$no_merge, @$c) unless($c->[0]{_remove});
	       next;
           }

	   my $start; 
	   my $end;
	   foreach my $h (@$c){
	       $start  = $h->start('query') if(!$start || $h->start('query') < $start);
	       $end = $h->end('query') if(!$end || $h->end('query') > $end);
	   }
	   
	   if($start <= $break && $break <= $end){
	       foreach my $h (@$c){
		   $h->{'_hsps'} = []; #memory optimization
		   $h->{'_reblast'} = 1;
	       }
	       push(@to_reblast, $c->[0]); #doesn't matter which one I just need the ID and DB

	       $low = $start if(!$low || $start < $low);
	       $high = $end if(!$high || $end > $high);
	   }
	   else{
	       foreach my $h (@$c){
		   push(@$no_merge, $h) unless($h->{_remove});
	       }
	   }
       }
   }

   my $merged = reblast_merged_hits($q_seq_obj,
				    $def,
				    $low,
				    $high,
				    $s_length,
				    \@to_reblast,
				    $index_files,
				    $the_void,
				    $type,
				    \%CTL_OPT,
				    $LOG);

   $blast_keepers = $no_merge;
   foreach my $h (@$merged){
       if($h->start <= $break && $break <= $h->end){
	   push(@$blast_keepers, $h);
       }
       elsif(abs($break - $h->start) <= $CTL_OPT{split_hit}){
	   push(@$blast_keepers, $h);
       }
       elsif(abs($break - $h->end) <= $CTL_OPT{split_hit}){
	   push(@$blast_keepers, $h);
       }
   }

   return $blast_keepers;
}
#-----------------------------------------------------------------------------
sub reblast_merged_hits {
   my $par_seq     = shift @_;
   my $par_def     = shift @_;
   my $nB          = shift @_;
   my $nE          = shift @_;
   my $seq_length  = shift @_;
   my $hits        = shift @_;
   my $index_files = shift @_;
   my $the_void    = shift @_;
   my $type        = shift @_;
   my %CTL_OPT = %{shift @_};
   my $LOG         = shift @_;

   return [] if(!@$hits);

   my @to_blast;
   my @blast_keepers;
   foreach my $hit (@{$hits}) {
       #if not a merged hit take as is
       if (! $hit->{'_sequences_was_merged'} && ! $hit->{'_reblast'}) {
	   die "ERROR: Logic error in reblast_merged_hits\n";
	   #push (@blast_keepers, $hit);
       }
       else{
	   push (@to_blast, $hit);
       }
   }
   return \@blast_keepers if(!@to_blast);

   #==get data from parent fasta
   #parent fasta get def and seq
   my $p_id  = Fasta::def2SeqID($par_def);
   my $p_safe_id = Fasta::seqID2SafeID($p_id);

   #== excise region of query fasta and build new chunk object for blast
   ($nB, $nE) = ($nE, $nB) if($nB > $nE);
   my $F = ($nB - $CTL_OPT{split_hit} < 1) ? 1 : $nB - $CTL_OPT{split_hit};
   my $L = ($nE + $CTL_OPT{split_hit} > $seq_length) ? $seq_length : $nE + $CTL_OPT{split_hit};

   #get input and output blast file names
   my $name = "db.$F-$L.for_$type";
   my $o_file = get_blast_finished_name(0, $name, $the_void, "$p_safe_id.$F.$L", $type);

   my $tmp = get_global_temp();
   my $rank = RANK();
   my $t_dir = "$tmp/$rank";
   my $t_file = "$t_dir/$name.fasta";

   #extract sequence
   my $piece_seq;
   if(!ref($par_seq) || ref($par_seq) eq 'SCALAR'){
       my @coors = [$nB, $nE];
       my $piece = Shadower::getPieces($par_seq, \@coors, $CTL_OPT{split_hit});
       $piece_seq = $piece->[0]->{piece};
   }
   else{ #seq object
       $piece_seq = $par_seq->subseq($F, $L);
   }
   
   #get piece fasta def, seq, and offset
   my $piece_def = $par_def." ".$F." ".$L;
   my $offset = $F - 1;
   
   #will require chunk to blast
   my $chunk = new FastaChunk();
   $chunk->number(0);
   $chunk->seq(\$piece_seq);
   $chunk->def($piece_def);
   $chunk->parent_def($par_def);
   $chunk->size($seq_length); #the max size of a chunk
   $chunk->length($L - $F + 1); #the actual size of a chunk
   $chunk->offset($offset);
   $chunk->is_last(1);
   $chunk->is_first(1);
   $chunk->parent_seq_length($seq_length);

   if(! -f $o_file){                   
       #build hit db
       if(! -f $t_file){
	   my $db_index = build_fasta_index($index_files);
	   open(my $FA, ">$t_file.tmp");
	   foreach my $hit (@to_blast) {
	       #build a safe name for file names from the sequence identifier
	       my $t_safe_id = Fasta::seqID2SafeID($hit->name());
	       
	       #==build new fasta and db for blast search from hit name and db index
	       
	       #search db index
	       my $fastaObj = $db_index->get_Seq_for_hit($hit);
	       
	       #still no sequence? try rebuilding the index and try again
	       if(!$fastaObj) {
		   for(my $i = 0; $i < 2 && !$fastaObj; $i++){
		       sleep 5;
		       print STDERR "WARNING: Cannot find >".$hit->name.", trying to re-index the fasta.\n" if($i);
		       $db_index->reindex($i);
		       $fastaObj = $db_index->get_Seq_for_hit($hit);
		   }
		   
		   if (not $fastaObj) {
		       print STDERR "stop here:".$hit->name."\n";
		       confess "ERROR: Fasta index error\n";
		   }
	       }
	       
	       #get fasta def and seq
	       my $t_seq      = $fastaObj->seq();
	       my $t_def      = $db_index->header_for_hit($hit);
	       
	       #write fasta file
	       my $fasta = Fasta::toFastaRef('>'.$t_def, \$t_seq);
	       print $FA $$fasta;
	   }
	   close($FA);
	   File::Copy::move("$t_file.tmp", $t_file);
       }
   }

   #save labels
   my %labels;
   foreach my $hit (@{$hits}) {
       $labels{$hit->name} = $hit->{_label} if($hit->{_label});
   }
   
   #==run the blast search
   my $entry = $t_file;
   if ($type eq 'blastx') {
       print STDERR "re-running blast for edge cases...\n" unless $main::quiet;
       my $keepers = blastx($chunk,
			    $t_file,
			    $the_void,
			    "$p_safe_id.$F.$L",
			    \%CTL_OPT,
			    $LOG);
       
       push(@blast_keepers, @{$keepers});
       print STDERR "...finished\n" unless $main::quiet;
   }
   elsif ($type eq 'blastn') {
       print STDERR "re-running blast for edge cases...\n" unless $main::quiet;
       my $keepers = blastn($chunk, 
			    $t_file,
			    $the_void,
			    "$p_safe_id.$F.$L",
			    \%CTL_OPT,
			    $LOG);
       
       push(@blast_keepers, @{$keepers});
       print STDERR "...finished\n" unless $main::quiet;
   }
   elsif ($type eq 'tblastx') {
       print STDERR "re-running blast for edge cases...\n" unless $main::quiet;
       my $keepers = tblastx($chunk,
			     $t_file,
			     $the_void,
			     "$p_safe_id.$F.$L",
			     \%CTL_OPT,
			     $LOG);
       
       push(@blast_keepers, @{$keepers});
       print STDERR "...finished\n" unless $main::quiet;
   }
   else {
       confess "ERROR: Invalid type \'$type\' in maker::reblast_merged_hit\n";
   }
   unlink($t_file);
   
   #restore labels
   foreach my $hit (@blast_keepers){
       if($labels{$hit->name}){
	   $hit->{_label} = $labels{$hit->name};
	   $_->{_label} = $labels{$hit->name} foreach($hit->hsps);
       }
   }
   
   #mark merged hits as having been heldover once
   foreach my $h (@blast_keepers){
      $h->{_holdover} = 1;
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
    my $a_flag     = shift @_; #indicates whether to have a data adjusted floating cutoff
    my $b_flag     = shift @_; #indicates which end to have cutoff on (-1 => start, 0 => both, 1 => end)
    my $groups_cfh = shift @_; #group to cluster and find holdovers

    $b_flag = 1 if($chunk->is_first && $b_flag == 0);
    $b_flag = -1 if($chunk->is_last && $b_flag == 0);

    my $phat_hits;
    
    foreach my $group (@{$groups_cfh}) {
	push(@{$phat_hits}, @{$group});
    }
    
    my $p_hits = [];
    my $m_hits = [];
    
    #separate by strand or not (this makes chunk cutoffs strand independant)
    if($s_flag){
	($p_hits, $m_hits) = PhatHit_utils::separate_by_strand('query', $phat_hits, 1); #exonerate flag set
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
    
    my @keepers;
    my @holdovers;

    #no internal cutoff if this is the last contig, return everything
    if($chunk->is_last && $b_flag == 1) {
        foreach my $g (@{$groups_cfh}){
            push (@holdovers, []);
            push (@keepers, $g);
        }
        return @keepers, @holdovers;
    }

    #no internal cutoff if this is the first contig, return everything
    if($chunk->is_first && $b_flag == -1) {
        foreach my $g (@{$groups_cfh}){
            push (@holdovers, []);
            push (@keepers, $g);
        }
        return @keepers, @holdovers;
    }

    #set up cutoffs
    my $cutoff = $chunk->length + $chunk->offset - $split_hit;
    my $p_cutoff = $chunk->length + $chunk->offset + 1;
    my $m_cutoff = $chunk->length + $chunk->offset + 1;
    my $scutoff = $chunk->offset + 1 + $split_hit;
    my $p_scutoff = $chunk->offset + 1 - 1;
    my $m_scutoff = $chunk->offset + 1 - 1;

    #adjust cutoff to overlapping hits
    if($a_flag){
       my $p_pieces = Shadower::getVectorPieces($p_coors, $pred_flank);
       $p_pieces = [sort {$b->{e} <=> $a->{e}} @{$p_pieces}];
       my $m_pieces = Shadower::getVectorPieces($m_coors, $pred_flank);
       $m_pieces = [sort {$b->{e} <=> $a->{e}} @{$m_pieces}];
       
       foreach my $p_piece (@{$p_pieces}) {
	  if ($b_flag >= 0 && $p_piece->{e} + $chunk->offset >= $cutoff) {
	     $p_cutoff = $p_piece->{b} + $chunk->offset;
	  }
	  if($b_flag <= 0 && $p_piece->{b} + $chunk->offset <= $scutoff){
	     $p_scutoff = $p_piece->{e} + $chunk->offset;
	  }
       }

       foreach my $m_piece (@{$m_pieces}) {
	  if ($b_flag >= 0 && $m_piece->{e} + $chunk->offset >= $cutoff) {
	     $m_cutoff = $m_piece->{b} + $chunk->offset;
	  }
	  if ($b_flag <= 0 && $m_piece->{b} + $chunk->offset <= $cutoff) {
	     $m_scutoff = $m_piece->{e} + $chunk->offset;
	  }
       }
    }
    else{
	  $p_cutoff = $cutoff;
	  $m_cutoff = $cutoff;
	  $p_scutoff = $scutoff;
	  $m_scutoff = $scutoff;
    }

    #too small, all are heldover for next round
    if ($b_flag == 1 &&
	$p_cutoff <= 1 + $chunk->offset &&
	$m_cutoff <= 1 + $chunk->offset
	) {
       foreach my $g (@{$groups_cfh}){
	  push (@holdovers, $g);
	  push (@keepers, []);
       }
    }
    #too small, all are heldover for next round
    elsif ($b_flag == -1 &&
	   $p_scutoff >= $chunk->length + $chunk->offset &&
	   $m_scutoff >= $chunk->length + $chunk->offset
	   ) {
       foreach my $g (@{$groups_cfh}){
	  push (@holdovers, $g);
	  push (@keepers, []);
       }
    }
    else{
       #separate holdovers and keepers
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
	     	     
	     $hit->{_holdover} = 0; #initialize to 0
	     if ($b_flag >= 0 &&
		 (($strand eq '1' && ($e >= $p_cutoff || $p_cutoff <= $chunk->offset +1)) ||
		  ($strand eq '-1' && ($e >= $m_cutoff || $m_cutoff <= $chunk->offset +1)))
		){
		$hit->{_holdover} += 2;
		push(@{$group_holdovers}, $hit);
	     }

	     if($b_flag <= 0 &&
		   (($strand eq '1' && ($b <= $p_scutoff || $p_scutoff >= $chunk->offset + $chunk->length)) ||
		    ($strand eq '-1' && ($b <= $m_cutoff || $m_cutoff >= $chunk->offset + $chunk->length)))
		  ){
		$hit->{_holdover} += 1;
                push(@{$group_holdovers}, $hit) unless($hit->{_holdover} != 1); #spanners will != 1
	     }	     
	     if($hit->{_holdover} == 0){
		push(@{$group_keepers}, $hit);
	     }
	  }
	  
	  push(@keepers, $group_keepers);
	  push(@holdovers, $group_holdovers);
       }
    }

    #keepers and hit holdovers are returned in same order given by user
    return @keepers, @holdovers;
}
#-----------------------------------------------------------------------------
sub maker_p_and_t_fastas {
   my $maker    = shift @_;
   my $non_over = shift @_;
   my $all      = shift @_;
   my $p_fastas = shift @_;
   my $t_fastas = shift @_;

   my $abinit = [];
   my $ncrna = [];
   my $mod_gff = [];
   my @ab_keys = grep {/^pred_gff|_abinit$/} keys %$all;
   my @mod_keys = grep {/^model_gff/} keys %$all;
   my @nc_keys = grep {/^ncrna_|_ncrna$/} keys %$all;

   foreach my $key (@ab_keys){
       push(@$abinit, @{$all->{$key}});
   }   

   foreach my $key (@mod_keys){
       push(@$mod_gff, @{$all->{$key}});
   }

   foreach my $key (@nc_keys){
       push(@$ncrna, @{$all->{$key}});
   }

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
	   my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a, 'abinit');
	   my $source = $a->{hit}->algorithm;
	   $source = uri_escape($source, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
	   $p_fastas->{$source} .= $p_fasta;
	   $t_fastas->{$source} .= $t_fasta;
       }
   }

   foreach my $an (@$mod_gff) {
       foreach my $a (@{$an->{t_structs}}) {
	   my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a, 'model_gff');
	   my $source = $a->{hit}->algorithm;
	   $source = uri_escape($source, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
	   $p_fastas->{$source} .= $p_fasta;
	   $t_fastas->{$source} .= $t_fasta;
       }
   }

   foreach my $an (@$ncrna) {
       foreach my $a (@{$an->{t_structs}}) {
	   my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
	   my $source = $a->{hit}->algorithm;
	   $source = uri_escape($source, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
	   $t_fastas->{$source} .= $t_fasta;
       }
   }
}

#-----------------------------------------------------------------------------
sub get_p_and_t_fastas {
   my $t_struct = shift;
   my $type     = shift || 'maker';
	
   my $t_seq  = $t_struct->{t_seq};
   my $p_seq  = $t_struct->{p_seq};
   my $t_name = $t_struct->{t_name} || $t_struct->{t_id};
   my $t_off  = "offset:".$t_struct->{t_offset};
   my $AED    = "AED:".sprintf('%.2f', $t_struct->{AED});
   my $eAED   = "eAED:".sprintf('%.2f', $t_struct->{eAED});
   my $QI     = "QI:".$t_struct->{t_qi};

   #for ab initio annotations use the stats from the unmodified model 
   if($type ne 'maker'){
       my $p_struct = ($type eq 'abinit') ? $t_struct->{p_struct} : $t_struct;
       $t_name = $p_struct->{t_name} || $p_struct->{t_id};
       $t_seq = $p_struct->{t_seq};
       $p_seq = $p_struct->{p_seq};
       $t_off = "offset:".$p_struct->{t_offset};
       $AED = (defined $p_struct->{AED}) ? "AED:".sprintf('%.2f', $p_struct->{AED}) : '';
       $eAED = (defined $p_struct->{eAED}) ? "eAED:".sprintf('%.2f', $p_struct->{eAED}) : '';
       $QI = ($p_struct->{t_qi}) ? "QI:".$p_struct->{t_qi} : '';
   }

   my $p_def = ">$t_name protein $AED $eAED $QI";
   my $t_def = ">$t_name transcript $t_off $AED $eAED $QI";

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
sub get_split_args {
    my $CTL_OPT = shift;
    my @set = @_ || qw(genome est altest protein repeat_protein);

    my @to_split;
    foreach my $in (@set){
	my @files = split(/\,/, $CTL_OPT->{$in});
	my %uniq = map {/^([^\:]+)\:?(.*)?/} @files;
	@files = map {($uniq{$_}) ? s_abs_path($_).":$_" : $_} keys %uniq;
	
	foreach my $file (@files){
	    my $bins = ($in eq 'genome') ? 1 : 10;
	    my $alt  = ($in =~ /protein/) ? $CTL_OPT->{mpi_blastdb} : undef;
	    my @args = ($file, $CTL_OPT->{mpi_blastdb}, $bins, $alt);
	    push(@to_split, \@args);
	}
    }

    return @to_split;
}
#----------------------------------------------------------------------------
{
my %localized; #persistent list of already localized files
sub localize_file {
    my $file = shift;
    my $name = shift;

    die "ERROR: Cannot localize non-existant file $file\n" if(! -f $file);

    $file = Cwd::abs_path($file);
    if(! length($name)){
       ($name) = $file =~ /([^\/]+)$/;
    }
    my $tmp = GI::get_global_temp();
    my $rank = GI::RANK();

    if($localized{$file} && $localized{$file} eq "$tmp/$name" && -f "$tmp/$name"){
	return $localized{$file};
    }
    elsif($localized{$file} && $localized{$file} ne "$tmp/$name" && -f $localized{$file}){
       symlink($localized{$file}, "$tmp/$name") if(! -f "$tmp/$name");
       return "$tmp/$name";
    }

    while(!-f "$tmp/$name"){
	my $lock = new File::NFSLock("$tmp/$name", 'EX', 1800, 60);
	next if(!$lock);

	#always check again after locking
	last if(-f "$tmp/$name");
	next unless($lock->maintain(30));

	File::Copy::copy($file, "$tmp/$rank/$name.tmp");
	File::Copy::move("$tmp/$rank/$name.tmp", "$tmp/$name");
	$lock->unlock if($lock);
    }    
    
    $localized{$file} = "$tmp/$name";
    return $localized{$file};
}
}
#----------------------------------------------------------------------------
sub create_blastdb {
   my $CTL_OPT = shift @_;

   my $bdir = $CTL_OPT->{mpi_blastdb};
   mkdir($bdir) if(! -d $bdir);

   my %source2dest = (genome         => '_g_db',
		      protein        => '_p_db',
		      est            => '_e_db',
		      altest         => '_a_db',
		      repeat_protein => '_r_db');

   foreach my $in (List::Util::shuffle(keys %source2dest)){
       $CTL_OPT->{$source2dest{$in}} = [];
       my @files = split(/\,/, $CTL_OPT->{$in});

       my $key  = ($in =~ /protein/) ? 'protein' : 'nucleotide';
       my $bins = ($in eq 'genome') ? 1 : 10;
       my $bdir = $CTL_OPT->{mpi_blastdb};
       my $alt  = $CTL_OPT->{alt_peptide};

       my @to_do = map {[$_, $key, $bins, $bdir, $alt]} @files;
       foreach my $args (@to_do){
	   my $db = split_db(@$args);
	   push(@{$CTL_OPT->{$source2dest{$in}}}, @$db);
       }
   }

   #handle maker_coor hints given in fasta headers
   if($CTL_OPT->{est_forward}){
       my %to_exonerate;
       foreach my $db (@{$CTL_OPT->{_e_db}}){
	   open(IN, "< $db");
	   while(my $line = <IN>){
	       next unless($line =~ /maker_coor=([^\;\n]+)/);
	       my @coors = split(/,/, $1);
	       my ($id, $def) = $line =~ /^>([^\s]+)\s*([^\n]*)/;
	       $def =~ s/maker_coor=[^\;\n]+\;?//g;
	       
	       foreach my $coor (@coors){
		   $coor =~ /^([^\:]+)/;
		   push(@{$to_exonerate{$1}}, ">$id maker_coor=$coor\; $def");
	       }
	   }
	   close(IN);
       }

       $CTL_OPT->{_to_exonerate} = \%to_exonerate;
   }
}
#----------------------------------------------------------------------------
sub concatenate_files {
    my $infiles = shift;
    my $outfile = shift;

    my $tmp = get_global_temp();
    my $rank = RANK();
    my ($tFH, $t_file) = tempfile(DIR => "$tmp/$rank");
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
    my $entry = shift;
    my $key   = shift;
    my $bins  = shift;
    my $b_dir = shift;
    my $alt   = shift;

    if(! -d $b_dir){
       mkdir($b_dir);
       confess "ERROR: Could not create $b_dir\n-->$!" if(! -d $b_dir);
    }

    my $tmp = get_global_temp();
    my $rank = RANK();

    my @db_files;
    my ($file, $label) = $entry =~ /^([^\:]+)\:?(.*)/;    
    my ($f_name) = $file =~ /([^\/]+)$/;
    $f_name = uri_escape($f_name, '^a-zA-Z0-9\-\_');
    my $d_name = "$f_name\.mpi\.$bins";
    my $f_dir = "$b_dir/$d_name";
    my $t_dir = "$tmp/$rank/$d_name";
    $bins = ($bins > 1 && $bins < 30 && -s $file > 1000000000) ? 30 : $bins;

    #check if already finished
    my @t_db = map {($label) ? "$_:$label" : $_} grep {-f $_ && /$d_name\.\d+$/} File::Glob::bsd_glob("$f_dir/$d_name\.*");
    if(@t_db == $bins){ #use existing if right count
	push(@db_files, @t_db);
	return \@db_files;
    }

    #lock and check again
    my $lock;
    while(! $lock){
	$lock = new File::NFSLock($f_dir, 'EX', 1800, 60);

	#check if already finished
	my @t_db = map {($label) ? "$_:$label" : $_} grep {-f $_ && /$d_name\.\d+$/} File::Glob::bsd_glob("$f_dir/$d_name\.*");
	if(@t_db == $bins){ #use existing if right count
	    push(@db_files, @t_db);
	    carp "Calling File::NFSLock::unlock" if($main::debug);
	    $lock->unlock if($lock);
	    return \@db_files;
	}
	elsif($lock->maintain(30)){
	    carp "Calling File::Path::rmtree" if($main::debug);
	    File::Path::rmtree($f_dir) if(-d $f_dir);
	    last;
	}
	else{
	    undef $lock;
	}
    }    

    #set up a new database
    carp "Calling Iterator::Any::new" if($main::debug);
    my $fasta_iterator = new Iterator::Any(-fasta => $file, -gff => $file); #handle both cases
    my $count = 0;
    carp "Calling Iterator::Any::nextDef" if($main::debug);
    while($fasta_iterator->nextDef){
	$count++;
	last if($count > $bins * 10 || $bins == 1);
    }
    die "ERROR: The fasta file $file appears to be empty.\n" if(! $count);
    my $max = ($count > 10) ? int($count / 10) : 1; #min seq per bin
    carp "Calling Iterator::Any::new" if($main::debug);
    $fasta_iterator = new Iterator::Any(-fasta => $file, -gff => $file); #rebuild
    
    if ($max < $bins){
	$bins = $max;
	$d_name = "$f_name\.mpi\.$bins";
	$f_dir = "$b_dir/$d_name";
	$t_dir = "$tmp/$rank/$d_name";
    }
   
    #make needed output directories
    if(! -d $t_dir){
	carp "Calling mkdir" if($main::debug);
	mkdir($t_dir) or confess "ERROR: Could not create $t_dir\n-->$!";
    }

    if(-d "$f_dir"){ #on multi processors check if finished
	my @t_db = map {($label) ? "$_:$label" : $_} grep {-f $_ && /$d_name\.\d+$/} File::Glob::bsd_glob("$f_dir/$d_name\.*");
	
	if(@t_db == $bins){ #use existing if right count
	    push(@db_files, @t_db);
	    $lock->unlock if($lock);
	    return \@db_files;
	}
	else{ #remove if there is an error
	    carp "Calling File::Path::rmtree" if($main::debug);
	    File::Path::rmtree($f_dir);
	}
    }
    
    #open filehandles for  pieces on multi processors
    my @fhs;    
    for (my $i = 0; $i < $bins; $i++) {
	my $name = "$t_dir/$d_name\.$i";
	my $fh;
	open ($fh, "> $name");
	push (@fhs, $fh);
    }
    
    #write fastas here
    my %alias;
    
    my $wflag = 1; #flag set so warnings gets printed only once 
    carp "Calling Iterator::Any::nextFastaRef" if($main::debug);
    while (my $fasta = $fasta_iterator->nextFastaRef()) {
	my $def = Fasta::getDef($fasta);
	my $seq_id = Fasta::def2SeqID($def);
	my $seq_ref = Fasta::fasta2seqRef($fasta);
	
	#fix non standard peptides
	if ($key =~ /protein/ && defined $alt) {
	    $$seq_ref =~ tr/[a-z]/[A-Z]/;
	    $$seq_ref =~ s/[\*\-]//g;
	    $$seq_ref =~ s/[^ACDEFGHIKLMNPQRSTVWYX\-\n]/$alt/g;
	}
	#fix nucleotide sequences
	elsif($key !~ /protein/){
	    #most programs use N for masking but for some reason the NCBI decided to
	    #use X to mask their sequence, which causes many many programs to fail
	    $$seq_ref =~ tr/[a-z]/[A-Z]/;
	    $$seq_ref =~ s/\-//g;
	    $$seq_ref =~ tr/XU/NT/;
	    $$seq_ref =~ s/[RYKMSWBDHV]/N/g if($main::fix_nucleotides);
	    die "ERROR: The nucleotide sequence file \'$file\'\n".
		"appears to contain protein sequence or unrecognized characters. Note\n".
		"the following nucleotides may be valid but are unsupported [RYKMSWBDHV]\n".
		"Please check/fix the file before continuing, or set -fix_nucleotides on\n".
		"the command line to fix this automatically.\n".
		"Invalid Character: '$1'\n\n"
		if($$seq_ref =~ /([^ATCGN\n])/i);
	}
	
	#Skip empty fasta entries
	next if($$seq_ref eq '');

	#fix weird super long headers
	if(length($def) > 10000){
	    warn "WARNING: Fasta header of length >10000.  Long headers kill BLAST.\n".
		"I will truncate it for you whether you like it or not\n\n";
	    $def = substr($def, 0, 10000);
	}

	#fix weird blast trimming error for long seq IDs by replacing them
	if(length($seq_id) > 78){
	    warn "WARNING: The fasta file contains sequences with names longer\n".
		"than 78 characters.  Long names get trimmed by BLAST, making\n".
		"it harder to identify the source of an alignmnet. You might\n".
		"want to reformat the fasta file with shorter IDs.\n".
		"File_name:$file\n\n" if($wflag-- > 0);
	    
	    carp "Calling Digest::MD5::md5_base6" if($main::debug);
	    carp "Calling uri_escape" if($main::debug);
	    my $new_id = uri_escape(Digest::MD5::md5_base64($seq_id), "^A-Za-z0-9\-\_");
	    
	    die "ERROR: The id $seq_id is too long for BLAST, and I can't uniquely fix it\n"
		if($alias{$new_id});
	    
	    $alias{$new_id}++;
	    $def =~ s/^(>\S+)/$1 MD5_alias=$new_id/;
	}
	
	#reformat fasta, just incase
	my $fasta_ref = Fasta::seq2fastaRef($def, $seq_ref);
	
	#build part files
	my $fh = shift @fhs;
	print $fh $$fasta_ref;
	push (@fhs, $fh);
    }
    
    #close part file handles
    foreach my $fh (@fhs) {
	close ($fh);
    }
    
    #move finished file directory into place
    #File::Copy::move cannot move directories
    confess "ERROR: logic problem. $f_dir already exists. Cannot replace.\n" if(-d $f_dir);
    carp "Calling system" if($main::debug);
    mkdir($f_dir);
    mkdir("$f_dir.$rank"); #temporary directory to hold active move
    system("mv $t_dir/* $f_dir.$rank/"); #move into temporary (slow)
    confess "ERROR: mv failed: $!\n" if($? == -1);
    system("mv $f_dir.$rank/* $f_dir/"); #move into active (fast)
    confess "ERROR: mv failed: $!\n" if($? == -1);
    File::Path::rmtree("$f_dir.$rank");
    
    #check if everything is ok
    if (-e $f_dir) { #multi processor
	my @t_db = map {($label) ? "$_:$label" : $_} grep {-f $_ && /$d_name\.\d+$/} File::Glob::bsd_glob("$f_dir/$d_name\.*");
	
	confess "ERROR: SplitDB not created correctly\n\n" if(@t_db != $bins);
	
	push(@db_files, @t_db);
    }
    else {
	confess "ERROR: Could not split db\n"; #not ok
    }
    
    carp "Calling File::NFSLock::unlock" if($main::debug);
    $lock->unlock() if $lock;
    return \@db_files;
}
#-----------------------------------------------------------------------------
sub flatten {
   my $clusters = shift;
   my $type     = shift;
   my @hits;
   foreach my $c (@{$clusters}) {
      if(ref($c) ne 'ARRAY'){
	 push(@hits, $c);
	 next;
      }
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
sub parse_abinit_file {
   my $file  = shift;
   my $entry = shift;
   my $chunk = shift;
   my $type  = shift;

   ($type) = $file =~ /\.([^\.]+)$/ if(!$type);

   my %params;
   my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my $keepers = "Widget::$type"->parse($file,
					\%params,
					$chunk);
   

   PhatHit_utils::add_offset($keepers,
                             $chunk->offset_w_flank());

   my $masked = 1 if($file =~ /\.abinit\_masked\.\d+\..*.$type$/);

   foreach my $h (@$keepers){
       $h->{_HMM}   = $hmm;
       if($label){
	   $h->{_label} = $label;
	   $_->{_label} = $label foreach($h->hsps);
       }
       if($masked){
	   my $alg = $h->algorithm();
	   $h->algorithm("$alg\_masked");
	   $_->algorithm("$alg\_masked") foreach($h->hsps);
       }
   }

   return $keepers;
}
#-----------------------------------------------------------------------------
sub snap {
   my $in_file     = shift;
   my $the_void    = shift;
   my $CTL_OPT     = shift;
   my $LOG         = shift;

   my $exe    = $CTL_OPT->{snap};
   my @entries = split(',', $CTL_OPT->{snaphmm});

   #make sure ZOE is set or snap can fail
   $ENV{ZOE} = $CTL_OPT->{ZOE} if($CTL_OPT->{ZOE} && -d $CTL_OPT->{ZOE});
   if(!$ENV{ZOE} || ! -d $ENV{ZOE}){
       #try and find it
       my ($path) = Cwd::abs_path($CTL_OPT->{snap});
       $path =~ s/snap$//;
       $ENV{ZOE} = $path;
   }

   my @out_files;
   foreach my $entry (@entries){
       my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
       my ($hmm_n) = $hmm =~ /([^\/]+)$/;
       $hmm_n = uri_escape($hmm_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
       my $out_file = "$in_file\.$hmm_n\.snap";
       (my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
       $LOG->add_entry("STARTED", $backup, "");
       
       my $command  = $exe;
       $command .= " $hmm";
       $command .= " $in_file";
       $command .= " > $out_file";
       
       my $w = new Widget::snap();
       if (-f $backup) {
	   print STDERR "using existing snap report.\n" unless $main::quiet;
	   print STDERR "$backup\n" unless $main::quiet;
	   $out_file = $backup;
       }
       else {
	   print STDERR "running  snap.\n" unless $main::quiet;
	   $w->run($command);
	   File::Copy::copy($out_file, $backup);
	   unlink($out_file);
       }
       $LOG->add_entry("FINISHED", $backup, "");

       push(@out_files, [$backup, $entry]);
   }

   return \@out_files;
}
#-----------------------------------------------------------------------------
sub genemark {
   my $in_file     = shift;
   my $the_void    = shift;
   my $CTL_OPT     = shift;
   my $LOG         = shift;

   #genemark sometimes fails if called directly so I built a wrapper
   my $wrap = "$FindBin::Bin/../lib/Widget/genemark/gmhmm_wrap";
   my $exe  = $CTL_OPT->{organism_type} eq 'eukaryotic' ? $CTL_OPT->{gmhmme3} : $CTL_OPT->{gmhmmp}; #genemark
   my $pro = $CTL_OPT->{probuild}; #helper exe

   my @entries = split(',', $CTL_OPT->{gmhmm});

   my @out_files;
   foreach my $entry (@entries){
       my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
       my ($hmm_n) = $hmm =~ /([^\/]+)$/;
       $hmm_n = uri_escape($hmm_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
       my $out_file = "$in_file\.$hmm_n\.genemark";
       (my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
       $LOG->add_entry("STARTED", $backup, "");

       my $command  = "$^X $wrap";
       $command .= " -m $hmm";
       $command .= " -g $exe";
       $command .= " -p $pro";
       $command .= " -o $out_file";
       #$command .= " -t $TMP";
       $command .= " $in_file";
       
       my $w = new Widget::genemark();
       if (-f $backup) {
           print STDERR "using existing genemark report.\n" unless $main::quiet;
           print STDERR "$backup\n" unless $main::quiet;
	   $out_file = $backup;
       }
       else {
	   print STDERR "running  genemark.\n" unless $main::quiet;
	   $w->run($command);
	   File::Copy::copy($out_file, $backup);
	   unlink($out_file);
       }
       $LOG->add_entry("FINISHED", $backup, "");

       push(@out_files, [$backup, $entry]);
   }

   return \@out_files;
}
#-----------------------------------------------------------------------------
sub augustus {
   my $in_file     = shift;
   my $the_void    = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;

   my $exe = $CTL_OPT->{augustus};
   my @entries = split(',', $CTL_OPT->{augustus_species});

   #make sure AUGUSTUS_CONFIG_PATH is set or augustus can fail
   $ENV{AUGUSTUS_CONFIG_PATH} = $CTL_OPT->{AUGUSTUS_CONFIG_PATH} if($CTL_OPT->{AUGUSTUS_CONFIG_PATH} &&
								    ! -f $CTL_OPT->{AUGUSTUS_CONFIG_PATH}."/extrinsic/extrinsic.MPE.cfg");
   if (! $ENV{AUGUSTUS_CONFIG_PATH} || ! -f "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg") {
       #try and find it
       my ($path) = Cwd::abs_path($CTL_OPT->{augustus});
       $path =~ s/bin\/augustus$/config/;
       $ENV{AUGUSTUS_CONFIG_PATH} = $path;
   }

   my @out_files;
   foreach my $entry (@entries){
       my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
       my ($hmm_n) = $hmm =~ /([^\/]+)$/;
       $hmm_n = uri_escape($hmm_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
       my $out_file = "$in_file\.$hmm_n\.augustus";
       (my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
       $LOG->add_entry("STARTED", $backup, ""); 
       
       my $command  = $exe;
       $command .= " --species=$hmm";
       $command .= " --UTR=off"; #added 3/19/2009
       $command .= " $in_file";
       $command .= " > $out_file";
       
       my $w = new Widget::augustus();
       if (-f $backup) {
           print STDERR "using existing augustus report.\n" unless $main::quiet;
           print STDERR "$backup\n" unless $main::quiet;
	   $out_file = $backup;
       }
       else {
	   print STDERR "running  augustus.\n" unless $main::quiet;
	   $w->run($command);
	   File::Copy::copy($out_file, $backup);
	   unlink($out_file);
       }
       $LOG->add_entry("FINISHED", $backup, "");

       push(@out_files, [$backup, $entry]);
   }

   return \@out_files;
}
#-----------------------------------------------------------------------------
sub fgenesh {
   my $in_file     = shift;
   my $the_void    = shift;
   my $CTL_OPT = shift;
   my $LOG         = shift;

   my $wrap = "$FindBin::Bin/../lib/Widget/fgenesh/fgenesh_wrap"; #fgenesh wrapper
   my $exe = $CTL_OPT->{fgenesh};

   my @entries = split(',', $CTL_OPT->{fgenesh_par_file});

   my @out_files;
   foreach my $entry (@entries){
       my ($hmm, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
       my ($hmm_n) = $hmm =~ /([^\/]+)$/;
       $hmm_n = uri_escape($hmm_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
       my $out_file = "$in_file\.$hmm_n\.fgenesh";       
       (my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
       $LOG->add_entry("STARTED", $backup, "");

       my $command  = "$^X $wrap $exe";
       $command .= " $hmm";
       $command .= " $in_file";
       $command .= " > $out_file";
       
       my $w = new Widget::fgenesh();
       if (-f $backup) {
           print STDERR "using existing fgenesh report.\n" unless $main::quiet;
           print STDERR "$backup\n" unless $main::quiet;
	   $out_file = $backup;
       }
       else {
	   print STDERR "running  fgenesh.\n" unless $main::quiet;
	   $w->run($command);
	   File::Copy::copy($out_file, $backup);
	   unlink($out_file);
       }
       $LOG->add_entry("FINISHED", $backup, "");

       push(@out_files, [$backup, $entry]);
   }

   return \@out_files;
}
#-----------------------------------------------------------------------------
sub snoscan {
    my $in_file     = shift;
    my $the_void    = shift;
    my $CTL_OPT     = shift;
    my $LOG         = shift;

    my $exe         = $CTL_OPT->{snoscan};
    my @entries = split(',', $CTL_OPT->{snoscan_rrna});

    my @out_files;
    foreach my $entry (@entries){
        my ($rna_file, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
        my ($rna_n) = $rna_file =~ /([^\/]+)$/;
        $rna_n = uri_escape($rna_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
        my $out_file = "$in_file\.$rna_n\.snoscan";
	(my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
	$LOG->add_entry("STARTED", $backup, "");

	my $command  = $exe;
	$command .= " $rna_file";
	$command .= " $in_file";
	$command .= " > $out_file";

	my $w = new Widget::snoscan();
	if (-f $backup) {
	    print STDERR "using existing snoscan report.\n" unless $main::quiet;
	    print STDERR "$backup\n" unless $main::quiet;
	    $out_file = $backup;
	}
	else {
	    print STDERR "running  snoscan.\n" unless $main::quiet;
	    $w->run($command);
	    File::Copy::copy($out_file, $backup) unless($CTL_OPT->{clean_up});
	}
	$LOG->add_entry("FINISHED", $backup, "");
	
	push(@out_files, [$out_file, $entry]);
    }
    
    return \@out_files;
}
#-----------------------------------------------------------------------------
sub trnascan {
    my $in_file     = shift;
    my $the_void    = shift;
    my $CTL_OPT     = shift;
    my $LOG         = shift;
    
    my $exe    = $CTL_OPT->{'tRNAscan-SE'};

    my @outfiles;
    my $entry = $CTL_OPT->{organism_type};
    my $out_file = "$in_file.$entry.trnascan";
    (my $backup = $out_file) =~ s/.*\/([^\/]+)$/$the_void\/$1/;
    $LOG->add_entry("STARTED", $backup, "");

    my $command  = $exe;
    $command .= " -b -q";
    $command .= " -B" if($CTL_OPT->{organism_type} ne 'eukaryotic');
    $command .= " $in_file";
    $command .= " > $out_file";

    my $w = new Widget::trnascan();
    if (-f $backup) {
        print STDERR "using existing trnascan report.\n" unless $main::quiet;
        print STDERR "$backup\n" unless $main::quiet;
        $out_file = $backup;
    }
    else {
        print STDERR "running  trnascan.\n" unless $main::quiet;
        $w->run($command);
	File::Copy::copy($out_file, $backup) unless($CTL_OPT->{clean_up});
    }
    $LOG->add_entry("FINISHED", $backup, "");
    push(@outfiles, [$out_file, $entry]);
    
    return \@outfiles;
}
#-----------------------------------------------------------------------------
sub polish_exonerate {
    my $chunk        = shift;
    my $q_seq_obj    = shift;
    my $length       = shift;
    my $def          = shift;
    my $phat_hits    = shift;
    my $db_set       = shift;
    my $the_void     = shift;
    my $type         = shift;
    my $exonerate    = shift;
    my $pcov         = shift;
    my $pid          = shift;
    my $score_limit  = shift;
    my $max_intron   = shift;
    my $matrix       = shift;
    my $pred_flank   = shift;
    my $est_forward  = shift;
    my $LOG          = shift;

    return [] if(! @{$phat_hits});

    my $db_index = GI::build_fasta_index($db_set);
    my $name =  Fasta::def2SeqID($def);
    my $safe_name = Fasta::seqID2SafeID($name);
    my $exe = $exonerate;
    
    my @exonerate_data;
    
    my %uniq; #make sure this same exonerate hit does not already exist
    foreach my $hit (@{$phat_hits}) {
	my $h_name;
	my $h_description;
	my $B;
	my $E;
	my $min_intron = 20; #constant

	#scalar for exonerate
	if(ref($hit) eq '' && $hit =~ />([^\s]+)\s+(.*)/){
	    $h_name = $1;
	    $h_description = $2;
	    $B = $chunk->offset + 1;
	    $E = $chunk->offset + $chunk->length;
	}
	#hit from blastn
	else{
	    $h_name = $hit->name;
	    $h_description = $hit->description;

	    next if($hit->pAh < $pcov/2);
	    next if($hit->hsp('best')->frac_identical < $pid &&
		    $hit->frac_identical < $pid);

	    ($B, $E) = PhatHit_utils::get_span_of_hit($hit,'query');
	    ($B, $E) = ($E, $B) if($B > $E);
	}

	#gene id specified for est_forward
	my $gene_id;
	if($est_forward && $h_description =~ /gene_id\=([^\s\;]+)/){
	    $gene_id = $1;
	}

	#check if fasta contains coordinates for maker
	if($h_description =~ /maker_coor\=([^\s\;]+)/){
	    my $go;
	    $min_intron = 1;
	    $max_intron = 200000;
	    foreach my $coor (split(',', $1)){
		my ($bName, $bB, $bE);
		if(($bName, $bB, $bE) = $coor =~ /^([^\s\;]+)\:(\d+)\-(\d+)$/){
		    next if($type eq 'e' && ref($hit));
		    ($bB, $bE) = ($bE, $bB) if($bB > $bE);
		    
		    next if($name ne $bName);
		    next if(compare::compare($B, $E, $bB, $bE) eq '0');
		    next if(ref($hit) eq '' && ($bB < $B || $bB > $E));
		    
		    ($B, $E) = ($bB, $bE); #switch to specified coordinates
		    $go++;
		    last;
		}
		elsif((($bName) = $coor =~ /^([^\s\;]+)$/) && (ref($hit) ne '')){
		    next if($name ne $bName);
		    $go++;
		    last;
		}
	    }
	    next if(!$go); #skip because the tag limits the seq to only these coors
	}

	my $rank = GI::RANK();
	my $id      = $h_name;
	my $safe_id = Fasta::seqID2SafeID($id);
	my $F = ($B - 2*$pred_flank > 0) ? $B - 2*$pred_flank : 1;
	my $L = ($E + 2*$pred_flank > $length) ? $length : $E + 2*$pred_flank;
	my $offset  = $F - 1;
	my $tmp = get_global_temp();
        my $backup  = "$the_void/$safe_name.$F-$L.$safe_id";
	$backup    .= '.p_exonerate' if($type eq 'p');
	$backup    .= '.est_exonerate' if($type eq 'e');
	$backup    .= '.alt_exonerate' if($type eq 'a');
	my $o_tfile = "$tmp/$rank/$safe_name.$F-$L.$safe_id.$type.exonerate";
	my $t_file  = "$tmp/$rank/$safe_id.for.$F-$L.$rank.fasta";
	my $d_file  = "$tmp/$rank/$safe_name.$F-$L.$rank.fasta";

	my $d_len = abs($L - $F) + 1;
	my $t_len = (ref($hit)) ? $hit->length : undef;
	if(! -f $backup || ! ref($hit)){
	    #get fasta for EST/protein
	    my $fastaObj = $db_index->get_Seq_for_hit($hit);

	    #still no sequence? try rebuilding the index and try again
            if(!$fastaObj) {
                for(my $i = 0; $i < 2 && !$fastaObj; $i++){
                    sleep 5;
                    print STDERR "WARNING: Cannot find >".$h_name.", trying to re-index the fasta.\n" if($i);
                    $db_index->reindex($i);
                    $fastaObj = $db_index->get_Seq_for_hit($hit);
                }

                if (not $fastaObj) {
                    print STDERR "stop here:".$h_name."\n";
                    confess "ERROR: Fasta index error\n";
                }
            }
	    
	    my $seq     = $fastaObj->seq();
	    $t_len      = length($seq) if(!$t_len);
	    my $header  = $db_index->header_for_hit($hit);
	    my $fasta   = Fasta::toFastaRef('>'.$header, \$seq);
	    FastaFile::writeFile($fasta, $t_file);	
	}
	if(! -f $backup){
	    #get substring fasta of contig
	    my $p_seq = $q_seq_obj->subseq($F, $L);
	    my $p_def = $def." ".$F." ".$L;
	    my $p_fasta = Fasta::toFasta($p_def, \$p_seq);
	    FastaFile::writeFile($p_fasta, $d_file);

	    unlink($o_tfile) if(-f $o_tfile);
	}

	#run exonerate
	unlink($o_tfile) if(-f $o_tfile);
	$o_tfile = $backup if(-f $backup);
	$LOG->add_entry("STARTED", $backup, "") if(defined $LOG);
	my $exonerate_hits = to_polisher($d_file,
					 $t_file,
					 $o_tfile,
					 $d_len,
					 $t_len,
					 $the_void,
					 $offset,
					 $type,
					 $exe,
					 $score_limit,
					 $min_intron,
					 $max_intron,
					 $matrix
					 );

	#temp
	#make backup except on TACC cluster
	if($o_tfile ne $backup && $HOST !~ /tacc\.utexas\.edu/){
	    #File::Copy::move($o_tfile, $backup);
	}
	$LOG->add_entry("FINISHED", $backup, "") if(defined $LOG);

	#delete fastas
	unlink($d_file, $t_file, $o_tfile);

	#evaluate exonerate results
	my @keepers;
	foreach my $e (@{$exonerate_hits}) {
	    next if(! defined $e);
	    next if $e->pAh < $pcov;

	    #fix flipped hits when mapping ESTs to gene models as is
	    if($type eq 'e' && $est_forward && $e->strand('hit') == -1){
		$e = PhatHit_utils::copy($e, 'both');
	    }

	    #double check was_flipped
	    $e->{_was_flipped} = (ref($hit) && $e->strand ne $hit->strand) ? 1 : 0;

	    #add gene_id if specified
	    $e->{gene_id} = $gene_id if($gene_id);

	    #trim poly-A tails
	    if($type eq 'e' && ! $est_forward){
		my @hsps = reverse @{$e->sortedHSPs()};
		while(my $hsp = shift @hsps){
		    my $eseq = $hsp->seq('hit')->seq;
		    my $len = length($eseq);
		    my $a_count = $eseq =~ tr/Aa/Aa/;
		    if($a_count/$len >= 0.8){
			$e->hsps(\@hsps); #make new referece rather than direct reference
			next;
		    }
		    else{
			last;
		    }
		}
	    }

	    #make sure hit overlaps blast input (for large pred_flanks)
	    my ($eB, $eE) = PhatHit_utils::get_span_of_hit($e,'query');
	    ($eB, $eE) = ($eE, $eB) if($eB > $eE);

	    if (exonerate_okay($e) && compare::compare($B, $E, $eB, $eE)) {
		if(ref($hit)){
		    #tag the source blastn hit to let you know the counterpart
		    #exonerate hit was flipped to the other strand
		    $hit->{_exonerate_flipped} = $e->{_was_flipped};
		    $hit->{_keep} = 1;
		    $hit->type("exonerate:$type"); #set hit type (exonerate only)
		    $e->{_label} = $hit->{_label} if($hit->{_label});
		    map{$_->{_label} = $hit->{_label}} $e->hsps if($hit->{_label});
		}
		else{
		    $e->{_from_ref} = 1;
		}

		#uniq structure string to keep from adding same EST multiple times
		my $u_string = join('', (map {$_->nB('query').'..'.$_->nE('query').'..'} $e->hsps), "ID=".$e->name);
		next if($uniq{$u_string});
		$uniq{$u_string}++;

		if($est_forward){
		   my $score = $e->frac_identical * $e->pAh * 100;
		   $e->score($score);
		}

		push(@keepers, $e);
	    }
	}
	
	#remove ambiguous alternate alignments when hints are given (only keep 1)
	if($h_description =~ /maker_coor\=([^\s\;]+)/){
	    my $coor = $1;
	    if($coor =~ /$name\:(\d+)-(\d+)$/){
		($B, $E) = ($1, $2);
		my @perfect = grep {$_->start == $B && $_->end == $E} @keepers;
		@keepers = @perfect if(@perfect);
	    }
	    @keepers = (sort {($b->frac_identical * $b->pAh) <=> ($a->frac_identical * $a->pAh)} @keepers)[0];
	    $hit->{_exonerate_flipped} = $keepers[0]->{_was_flipped} if(@keepers && ref($hit));
	}

	push(@exonerate_data, @keepers);
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
	    #sleep 4;
	    return 0;
	}
	elsif ($q_str =~ /Target Intron/) {
	    print STDERR "BADDD EXONERATE!\n";
	    #sleep 4;
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
   my $o_file   = shift;
   my $d_len    = shift;
   my $t_len    = shift;
   my $the_void = shift;
   my $offset   = shift;
   my $type     = shift;
   my $exe      = shift;
   my $score_limit = shift;
   my $min_intron = shift;
   my $max_intron = shift;
   my $matrix = shift;

   if ($type eq 'p') {
      return polisher::exonerate::protein::polish($d_file,
						  $t_file,
						  $o_file,
						  $d_len,
						  $t_len,
						  $the_void,
						  $offset,
						  $exe,
						  $score_limit,
						  $min_intron,
						  $max_intron,
						  $matrix);
   }
   elsif ($type eq 'e') {
      return polisher::exonerate::est::polish($d_file,
					      $t_file,
					      $o_file,
					      $d_len,
					      $t_len,
					      $the_void,
					      $offset,
					      $exe,
					      $score_limit,
					      $min_intron,
					      $max_intron,
					      $matrix);
   }
   elsif ($type eq 'a') {
      return polisher::exonerate::altest::polish($d_file,
						 $t_file,
						 $o_file,
						 $d_len,
						 $t_len,
						 $the_void,
						 $offset,
						 $exe,
						 $score_limit,
						 $min_intron,
						 $max_intron,
						 $matrix);
   }
   else {
      confess "unknown type:$type in sub to_polisher.\n";
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
	 if(!$fastaObj) {
	     for(my $i = 0; $i < 2 && !$fastaObj; $i++){
		 sleep 5;
		 print STDERR "WARNING: Cannot find >".$hit->name.", trying to re-index the fasta.\n" if($i);
		 $index->reindex($i);
		 $fastaObj = $index->get_Seq_for_hit($hit);
	     }

	     if (not $fastaObj) {
		 print STDERR "stop here:".$hit->name."\n";
		 confess "ERROR: Fasta index error\n";
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
   my $dbs = shift;
   $dbs = [$dbs] if(! ref($dbs));

   my @files;
   foreach my $db (@$dbs){
       my ($file) = split(':', $db);
       push(@files, $file);
   }

   carp "Calling FastaDB::new" if($main::debug);
   my $index = new FastaDB(\@files);

   return $index;
}
#-----------------------------------------------------------------------------
sub build_all_indexes {
   my $CTL_OPT = shift;

   my @dbs = (@{$CTL_OPT->{_e_db}},
	      @{$CTL_OPT->{_p_db}},
	      @{$CTL_OPT->{_a_db}},
	      @{$CTL_OPT->{_g_db}}
	     );

   my $index = build_fasta_index(\@dbs);
}
#-----------------------------------------------------------------------------
sub dbformat {
   my $exe = shift;
   my ($file) = shift =~ /^([^\:]+)\:?(.*)/; #peal off label
   my $type = shift;

   confess "ERROR: Can not find xdformat, formatdb, or makeblastdb executable\n" if(! -e $exe);
   confess "ERROR: Can not find the db file $file\n" if(! -e $file);
   confess "ERROR: You must define a type (blastn|blastx|tblastx)\n" if(! $type);
   
   my $tmp = get_global_temp();
   my $rank = RANK();
   my ($name) = $file =~ /([^\/]+)$/;
   my $t_dir = "$tmp/$rank/blastprep";
   my $t_file = "$t_dir/$name";
   File::Path::rmtree($t_dir) if(-d $t_dir);
   
   my $lock;
   while(1){
       my $command = $exe;
       my $run;
       if ($exe =~ /xdformat/) {
	   if (($type eq 'blastn' && ! -e $file.'.xnd') ||
	       ($type eq 'blastx' && ! -e $file.'.xpd') ||
	       ($type eq 'tblastx' && ! -e $file.'.xnd')
	       ) {
	       $command .= " -p" if($type eq 'blastx');
	       $command .= " -n" if($type eq 'blastn' || $type eq 'tblastx');
	       $command .= " $t_file";
	       $run++;
	   }
       }
       elsif ($exe =~ /formatdb/) {
	   if (($type eq 'blastn' && ! -e $file.'.nsq') ||
	       ($type eq 'blastx' && ! -e $file.'.psq') ||
	       ($type eq 'tblastx' && ! -e $file.'.nsq')
	       ) {
	       $command .= " -p T" if($type eq 'blastx');
	       $command .= " -p F" if($type eq 'blastn' || $type eq 'tblastx');
	       $command .= " -i $t_file";
	       $run++;
	   }
       }
       elsif ($exe =~ /makeblastdb/) {
	   if (($type eq 'blastn' && ! -e $file.'.nsq') ||
	       ($type eq 'blastx' && ! -e $file.'.psq') ||
	       ($type eq 'tblastx' && ! -e $file.'.nsq')
	       ) {
	       $command .= " -dbtype prot" if($type eq 'blastx');
	       $command .= " -dbtype nucl" if($type eq 'blastn' || $type eq 'tblastx');
	       $command .= " -in $t_file";
	       $run++;
	   }
       }
       else {
	   confess "ERROR: databases can only be formated by xdformat, formatdb, or makeblastdb, not \'$exe\'\n";
       }

       if($run){
	   #use symlink in rank specific dir incase of broken NFSLocking
	   mkdir($t_dir) if(!-d $t_dir);
	   symlink($file, $t_file) if(!-f $t_file);

	   if(! $lock){
	       $lock = new File::NFSLock("$file.dbformat", 'EX', 1800, 60);
	       confess "ERROR:  Could not obtain lock to format database\n\n" if(!$lock);
	       next;
	   }
	   elsif($lock->maintain(30)){
	       my $w = new Widget::formater();
	       print STDERR "formating database...\n" unless $main::quiet;
	       $w->run($command);

	       #move the BLAST indexes into place
	       unlink($t_file);
	       my @files = File::Glob::bsd_glob("$t_dir/*");
	       File::Copy::move($_, $tmp) foreach(@files);
	       last;
	   }
	   else{
	       undef $lock;
	       next;
	   }
       }
       else{
	   last;
       }
   }

   $lock->unlock if($lock);

   return;
}
#-----------------------------------------------------------------------------
sub get_blast_finished_name {
    my $chunk_num  = shift;
    my $entry      = shift || '';
    my $the_void   = shift || '';
    my $seq_id     = shift || '';
    my $type       = shift || '';

    #peal off label and build names for files to use and copy
    my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
    my ($db_n) = $db =~ /([^\/]+)$/;
    $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\+');
    my $db_old_n = $db_n;
    $db_old_n =~ s/([^\/]+)\.mpi\.\d+\.\d+$/$1/;
    
    return "$the_void/$seq_id\.$chunk_num\.$db_old_n\.$type";
}
#-----------------------------------------------------------------------------
sub blastn_as_chunks {
   my $chunk      = shift;
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;
   my $retry      = shift;

   $retry = 1 if(! defined $retry);

   my $blast       = $CTL_OPT->{_blastn};
   my $depth_blast = $CTL_OPT->{depth_blastn};
   my $bit_blast   = $CTL_OPT->{bit_blastn};
   my $eval_blast  = $CTL_OPT->{eval_blastn};
   my $pcov_blast  = $CTL_OPT->{pcov_blastn};
   my $pid_blast   = $CTL_OPT->{pid_blastn};
   my $split_hit   = $CTL_OPT->{split_hit};
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

   #finished blast report name
   my $chunk_number = $chunk->number();
   my $blast_finished = get_blast_finished_name($chunk_number, $entry, $the_void, $seq_id, 'blastn');

   #peal off label and build names for files to use and copy
   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\+');

   my $tmp = get_global_temp();
   my $rank = RANK();
   my $t_dir = "$tmp/$rank";

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.blastn";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #parse blast
   if (-e $blast_finished) {
      return ([], $blast_dir) if(! $LOG_FLAG);

      print STDERR "re reading blast report.\n" unless ($main::quiet);
      print STDERR "$blast_finished\n" unless ($main::quiet);

      $o_file = $blast_finished;
   }
   else{
       #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'blastn');
   
       #call blast executable
       $chunk->write_file_w_flank($t_file_name);
       
       runBlastn($t_file_name,
		 $tmp_db,
		 $o_file,
		 $blast,
		 $eval_blast,
		 $split_hit,
		 $cpus,
		 $org_type,
		 $softmask);
       
       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = $split_hit;
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try{
      $chunk_keepers = Widget::blastn::parse($o_file,
					     \%params,
					     );
   }
   catch Error::Simple with {
      my $E = shift;
      
      if($retry){
	 unlink($o_file);
	 return blastn_as_chunks( $chunk,
				  $entry,
				  $the_void,
				  $seq_id,
				  $CTL_OPT,
				  $LOG,
				  $LOG_FLAG,
				  0 );
      }
      else{
	 throw $E;
      }
   };
   
   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset_w_flank(),
			    );

   #filter out hits that are not really on this chunk, just on flank
   if($chunk->length != $chunk->length_w_flank){
       my @keepers;
       foreach my $hit (@$chunk_keepers){
	   next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	   push(@keepers, $hit);
	   $chunk_keepers = \@keepers;
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

   return ($chunk_keepers, $blast_dir);
}
#-----------------------------------------------------------------------------
sub blastn {
   my $chunk      = shift;
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $retry      = shift;

   $retry = 1 if(! defined $retry);

   my $blast      = $CTL_OPT->{_blastn};
   my $depth_blast = $CTL_OPT->{depth_blastn};
   my $bit_blast  = $CTL_OPT->{bit_blastn};
   my $eval_blast = $CTL_OPT->{eval_blastn};
   my $pcov_blast = $CTL_OPT->{pcov_blastn};
   my $pid_blast  = $CTL_OPT->{pid_blastn};
   my $split_hit   = $CTL_OPT->{split_hit};
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
	
   my $chunk_number = $chunk->number();
   my $q_length = $chunk->parent_seq_length();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.blastn";

   $LOG->add_entry("STARTED", $o_file, ""); 

   #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
   if(! -e $o_file){
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'blastn');
       
       $chunk->write_file_w_flank($file_name);
       runBlastn($file_name,
		 $tmp_db,
		 $o_file,
		 $blast,
		 $eval_blast,
		 $split_hit,
		 $cpus,
		 $org_type,
		 $softmask);

       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = 0; #don't use holdover filter
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try{
      $chunk_keepers = Widget::blastn::parse($o_file,
					     \%params,
					     );
   }
   catch Error::Simple with {
      my $E = shift;

      if($retry){
	 unlink($o_file);
         return blastn( $chunk,
			$entry,
			$the_void,
			$seq_id,
			$CTL_OPT,
			$LOG,
			0 );
      }
      else{
         throw $E;
      }
   };

   $LOG->add_entry("FINISHED", $o_file, "");

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset_w_flank(),
			    );

   #filter out hits that are not really on this chunk, just on flank
   if($chunk->length != $chunk->length_w_flank){
       my @keepers;
       foreach my $hit (@$chunk_keepers){
	   next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	   push(@keepers, $hit);
	   $chunk_keepers = \@keepers;
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

   return $chunk_keepers;
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

   my $tmp = get_global_temp();
   my $command  = $blast;
   if ($command =~ /blasta$/) {
      symlink($blast, "$tmp/blastn") if(! -e "$tmp/blastn"); #handle blasta linking
      $command = "$tmp/blastn";
      $command .= " $db $q_file B=10000 V=10000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";
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
      $command .= " notes"; #suppress certain notes
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " nonnegok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      $command .= ($org_type eq 'eukaryotic') ? "" : " kap";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p blastn";
      $command .= " -d $db -i $q_file -b 10000 -v 10000 -e $eval_blast";
      $command .= " -E 3";
      $command .= " -W 15";
      $command .= " -r 1";
      $command .= " -q -3";
      $command .= " -G 3";
      $command .= " -z 1000";
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I T";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastn$/) {
      #NCBI BLAST+ needs to be more restrictive or it overaligns
      #$command .= " -task blastn";
      $command .= " -db $db -query $q_file";
      $command .= " -num_alignments 10000 -num_descriptions 10000 -evalue $eval_blast";
      $command .= " -word_size 28";
      $command .= " -reward 1";
      $command .= " -penalty -5";
      $command .= " -gapopen 5";
      $command .= " -gapextend 5";
      $command .= " -dbsize 1000";
      $command .= ($org_type eq 'eukaryotic') ? " -searchsp 500000000" : " -searchsp 20000000";
      $command .= " -num_threads $cpus";
      $command .= " -lcase_masking";
      $command .= " -dust yes";
      $command .= ($softmask) ? " -soft_masking true" : " -soft_masking false";
      $command .= " -show_gis";
      #$command .= " -outfmt 6"; # remove for full report
      $command .= " -out $out_file";
   }
   else{
      confess "ERROR: Must be a blastn executable";  
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
sub repeatrunner_as_chunks {
    blastx_as_chunks(@_, 1); #run blastx_as_chunks but with repeat flag
}
#-----------------------------------------------------------------------------
sub blastx_as_chunks {
   my $chunk      = shift;
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;
   my $rflag      = shift; #am I running repeatrunner?
   my $retry      = shift;

   $retry = 1 if(! defined $retry);

   my $blast       = $CTL_OPT->{_blastx};
   my $depth_blast = ($rflag) ? undef : $CTL_OPT->{depth_blastx};
   my $bit_blast   = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{bit_blastx};
   my $eval_blast  = ($rflag) ? $CTL_OPT->{eval_rm_blastx} : $CTL_OPT->{eval_blastx};
   my $pcov_blast  = ($rflag) ? $CTL_OPT->{pcov_rm_blastx} : $CTL_OPT->{pcov_blastx};
   my $pid_blast   = ($rflag) ? $CTL_OPT->{pid_rm_blastx} : $CTL_OPT->{pid_blastx};
   my $split_hit   = ($rflag) ? 0 : $CTL_OPT->{split_hit}; #repeat proteins get shatttered later anyway
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = ($rflag) ? 1 : $CTL_OPT->{softmask}; #always on for repeats
   my $org_type    = $CTL_OPT->{organism_type};

   #finished blast report name
   my $type = ($rflag) ? 'repeatrunner': 'blastx';
   my $chunk_number = $chunk->number();
   my $blast_finished = get_blast_finished_name($chunk_number, $entry, $the_void, $seq_id, $type);

   #peal off label and build names for files to use and copy
   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\+');

   my $tmp = get_global_temp();
   my $rank = GI::RANK();
   my $t_dir = "$tmp/$rank";

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.$type";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG);

   #parse blast
   if (-e $blast_finished) {
      return ([], $blast_dir) if(! $LOG_FLAG);

      print STDERR "re reading blast report.\n" unless ($main::quiet);
      print STDERR "$blast_finished\n" unless ($main::quiet);

      $o_file = $blast_finished;
   }
   else{
       #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'blastx');

       #call blast executable
       if($rflag){
	   $chunk->write_file($t_file_name);  
       }
       else{
	   $chunk->write_file_w_flank($t_file_name);  
       }

       runBlastx($t_file_name,
		 $tmp_db,
		 $o_file,
		 $blast,
		 $eval_blast,
		 $split_hit,
		 $cpus,
		 $org_type,
		 $softmask);
       
       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = $split_hit;
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try{
      $chunk_keepers = Widget::blastx::parse($o_file,
					     \%params,
					     );
   }
   catch Error::Simple with {
      my $E = shift;

      if($retry){
	 unlink($o_file);
         return blastx_as_chunks( $chunk,
				  $entry,
				  $the_void,
				  $seq_id,
				  $CTL_OPT,
				  $LOG,
				  $LOG_FLAG,
				  $rflag,
				  0 );
      }
      else{
         throw $E;
      }
   };   

   if($rflag){
       PhatHit_utils::add_offset($chunk_keepers,
				 $chunk->offset());
   }
   else{
       PhatHit_utils::add_offset($chunk_keepers,
				 $chunk->offset_w_flank());
       
       #filter out hits that are not really on this chunk, just on flank
       if($chunk->length != $chunk->length_w_flank){
	   my @keepers;
	   foreach my $hit (@$chunk_keepers){
	       next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	       push(@keepers, $hit);
	       $chunk_keepers = \@keepers;
	   }
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

   return ($chunk_keepers, $blast_dir);
}
#-----------------------------------------------------------------------------
sub combine_blast_report {
   my $chunk     = shift;
   my $hits      = shift;
   my $blast_dir = shift;
   my $LOG       = shift;

   my %uniq;
   @uniq{@$blast_dir} = map {1} (1..@$blast_dir);
   foreach my $dir (keys %uniq){
       my $blast_finished = $dir;
       $blast_finished =~ s/\.temp_dir$//;
       my ($ext) = $blast_finished =~ /([^\.]+)$/;

       #merge blast reports
       if ( -d $dir && ! -e $blast_finished) {
	   system ("cat $dir/\*\.$ext > $blast_finished")
	       && confess "ERROR: Could not colapse BLAST reports\n";
	   File::Path::rmtree ("$dir");
       }
       $LOG->add_entry("FINISHED", $blast_finished, ""); 
   }

   my @flat;
   foreach my $set (@$hits){
      if(ref($set) eq 'ARRAY'){
	 push(@flat, @$set);
      }
      else{
	 push(@flat, $set);
      }
   }

   return \@flat;
}
#-----------------------------------------------------------------------------
sub repeatrunner {
    blastx(@_, 1); #run with repeat flag
}
#-----------------------------------------------------------------------------
sub blastx {
   my $chunk      = shift;
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $rflag      = shift;
   my $retry      = shift;

   $retry = 1 if(! defined $retry); 

   my $blast      = $CTL_OPT->{_blastx};
   my $depth_blast = ($rflag) ? undef : $CTL_OPT->{depth_blastx};
   my $bit_blast  = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{bit_blastx};
   my $eval_blast = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{eval_blastx};
   my $pcov_blast = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{pcov_blastx};
   my $pid_blast  = ($rflag) ? $CTL_OPT->{bit_rm_blastx} : $CTL_OPT->{pid_blastx};
   my $split_hit   = ($rflag) ? 0 : $CTL_OPT->{split_hit}; #repeat proteins get shatttered later anyway
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = ($rflag) ? 1 : $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');

   my $q_length = $chunk->parent_seq_length();
   my $chunk_number = $chunk->number();
		
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.";
   $o_file .= ($rflag) ? 'repeatrunner': 'blastx';

   $LOG->add_entry("STARTED", $o_file, ""); 

   #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
   if(! -e $o_file){
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'blastx');
       
       if($rflag){
	   $chunk->write_file($file_name);
       }
       else{
	   $chunk->write_file_w_flank($file_name);
       }

       runBlastx($file_name,
		 $tmp_db,
		 $o_file,
		 $blast,
		 $eval_blast,
		 $split_hit,
		 $cpus,
		 $org_type,
		 $softmask);

       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = 0; #don't use holdover filter
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try{
      $chunk_keepers = Widget::blastx::parse($o_file,
					     \%params,
					     );
   }
   catch Error::Simple with {
      my $E = shift;
      
      if($retry){
	 unlink($o_file);
         return blastx( $chunk,
                        $entry,
                        $the_void,
                        $seq_id,
                        $CTL_OPT,
                        $LOG,
			$rflag,
                        0 );
      }
      else{
         throw $E;
      }
   };

   $LOG->add_entry("FINISHED", $o_file, "");

   if($rflag){
       PhatHit_utils::add_offset($chunk_keepers,
				 $chunk->offset());
   }
   else{
       PhatHit_utils::add_offset($chunk_keepers,
				 $chunk->offset_w_flank());
       
       #filter out hits that are not really on this chunk, just on flank
       if($chunk->length != $chunk->length_w_flank){
	   my @keepers;
	   foreach my $hit (@$chunk_keepers){
	       next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	       push(@keepers, $hit);
	       $chunk_keepers = \@keepers;
	   }
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

   return $chunk_keepers
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

   my $tmp = get_global_temp();
   my $command  = $blast;
   if ($command =~ /blasta$/) {
      symlink($blast, "$tmp/blastx") if(! -e "$tmp/blastx"); #handle blasta linking
      $command = "$tmp/blastx";
      $command .= " $db $q_file B=10000 V=10000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";
      $command .= " Z=300";
      $command .= ($org_type eq 'eukaryotic') ? " Y=500000000" : " Y=20000000";
      $command .= ($org_type eq 'eukaryotic') ? " hspmax=100" : " hspmax=5";
      $command .= " cpus=$cpus";
      $command .= ($org_type eq 'eukaryotic') ? " gspmax=100" : " gspmax=5";
      $command .= " lcmask";
      $command .= " maskextra=10";
      $command .= " kap";
      $command .= " gi";
      $command .= " warnings"; #suppress certain warnings
      $command .= " notes"; #suppress certain notes
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " nonnegok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p blastx";
      $command .= " -d $db -i $q_file -b 10000 -v 10000 -e $eval_blast";
      $command .= " -z 300";
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I T";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastx$/) {
      $command .= " -db $db -query $q_file";
      $command .= " -num_alignments 10000 -num_descriptions 10000 -evalue $eval_blast";
      $command .= " -dbsize 300";
      $command .= ($org_type eq 'eukaryotic') ? " -searchsp 500000000" : " -searchsp 20000000";
      $command .= " -num_threads $cpus";
      $command .= " -seg yes";
      $command .= ($softmask) ? " -soft_masking true" : " -soft_masking false";
      $command .= " -lcase_masking";
      $command .= " -show_gis";
      #$command .= " -outfmt 6"; # remove for full report
      $command .= " -out $out_file";
   }
   else{
      confess "ERROR: Must be a blastx executable";  
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
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $LOG_FLAG   = shift;
   my $retry      = shift;

   $retry = 1 if(! defined $retry);

   my $blast       = $CTL_OPT->{_tblastx};
   my $depth_blast = $CTL_OPT->{depth_tblastx};
   my $bit_blast   = $CTL_OPT->{bit_tblastx};
   my $eval_blast  = $CTL_OPT->{eval_tblastx};
   my $pcov_blast  = $CTL_OPT->{pcov_tblastx};
   my $pid_blast   = $CTL_OPT->{pid_tblastx};
   my $split_hit   = $CTL_OPT->{split_hit};
   my $cpus        = $CTL_OPT->{cpus};
   my $formater    = $CTL_OPT->{_formater};
   my $softmask    = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

   #finished blast report name
   my $chunk_number = $chunk->number();
   my $blast_finished = get_blast_finished_name($chunk_number, $entry, $the_void, $seq_id, 'tblastx');

   #peal off label and build names for files to use and copy
   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\+');

   my $tmp = get_global_temp();
   my $rank = GI::RANK();
   my $t_dir = "$tmp/$rank";

   my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
   my $blast_dir = "$blast_finished\.temp_dir";
   my $o_file    = "$blast_dir/$db_n\.tblastx";

   $LOG->add_entry("STARTED", $blast_finished, "") if($LOG_FLAG); 

   #parse blast
   if (-e $blast_finished) {
      return ([], $blast_dir) if(! $LOG_FLAG);

      print STDERR "re reading blast report.\n" unless ($main::quiet);
      print STDERR "$blast_finished\n" unless ($main::quiet);

      $o_file = $blast_finished;
   }
   else{
       #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'tblastx');

       #call blast executable
       $chunk->write_file_w_flank($t_file_name);  
       
       runtBlastx($t_file_name,
		  $tmp_db,  
		  $o_file, 
		  $blast,
		  $eval_blast,
		  $split_hit,
		  $cpus,
		  $org_type,
		  $softmask);
       
       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = $split_hit;
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try {
      $chunk_keepers = Widget::tblastx::parse($o_file,
					      \%params,
					      );
   }
   catch Error::Simple with {
      my $E = shift;
      
      if($retry){
	 unlink($o_file);
         return tblastx_as_chunks( $chunk,
				   $entry,
				   $the_void,
				   $seq_id,
				   $CTL_OPT,
				   $LOG,
				   $LOG_FLAG,
				   0 );
      }
      else{
	 throw $E;
      }
   };
   
   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset_w_flank(),
			    );

   #filter out hits that are not really on this chunk, just on flank
   if($chunk->length != $chunk->length_w_flank){
       my @keepers;
       foreach my $hit (@$chunk_keepers){
	   next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	   push(@keepers, $hit);
	   $chunk_keepers = \@keepers;
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

   return ($chunk_keepers, $blast_dir);
}
#-----------------------------------------------------------------------------
sub tblastx {
   my $chunk      = shift;
   my $entry      = shift;
   my $the_void   = shift;
   my $seq_id     = shift;
   my $CTL_OPT    = shift;
   my $LOG        = shift;
   my $retry      = shift;

   $retry = 1 if(! defined $retry);

   my $blast      = $CTL_OPT->{_tblastx};
   my $depth_blast = $CTL_OPT->{depth_tblastx};
   my $bit_blast  = $CTL_OPT->{bit_tblastx};
   my $eval_blast = $CTL_OPT->{eval_tblastx};
   my $pcov_blast = $CTL_OPT->{pcov_tblastx};
   my $pid_blast  = $CTL_OPT->{pid_tblastx};
   my $split_hit    = $CTL_OPT->{split_hit};
   my $cpus         = $CTL_OPT->{cpus};
   my $formater     = $CTL_OPT->{_formater};
   my $softmask     = $CTL_OPT->{softmask};
   my $org_type    = $CTL_OPT->{organism_type};

   my ($db, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $db =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');
	
   my $chunk_number = $chunk->number();
   my $q_length = $chunk->parent_seq_length();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   my $o_file    = "$the_void/$seq_id\.$chunk_number\.$db_n\.tblastx";

   $LOG->add_entry("STARTED", $o_file, ""); 

   if(! -e $o_file){
       #copy db to local tmp dir and run xdformat, formatdb, or makeblastdb
       my $tmp_db = localize_file($db);
       dbformat($formater, $tmp_db, 'tblastx');
       
       $chunk->write_file_w_flank($file_name);

       runtBlastx($file_name,
		  $tmp_db,
		  $o_file,
		  $blast,
		  $eval_blast,
		  $split_hit,
		  $cpus,
		  $org_type,
		  $softmask);

       $chunk->erase_fasta_file();
   }

   my %params;
   $params{significance}  = $eval_blast;
   $params{hsp_bit_min}   = $bit_blast;
   $params{percov}        = $pcov_blast;
   $params{percid}        = $pid_blast;
   $params{split_hit}     = 0; #don't use holdover filter
   $params{depth}         = $depth_blast;
   $params{is_first}      = $chunk->is_first;
   $params{is_last}       = $chunk->is_last;

   my $chunk_keepers;
   try {
      $chunk_keepers = Widget::tblastx::parse($o_file,
					      \%params,
					      );
   }
   catch Error::Simple with {
      my $E = shift;
      
      if($retry){
	 unlink($o_file);
	 return tblastx( $chunk,
			 $entry,
			 $the_void,
			 $seq_id,
			 $CTL_OPT,
			 $LOG,
			 0 );
      }
      else{
	 throw $E;
      }
   };

   $LOG->add_entry("FINISHED", $o_file, "");

   PhatHit_utils::add_offset($chunk_keepers,
			     $chunk->offset_w_flank(),
			    );

   #filter out hits that are not really on this chunk, just on flank
   if($chunk->length != $chunk->length_w_flank){
       my @keepers;
       foreach my $hit (@$chunk_keepers){
	   next if($hit->end < $chunk->start || $hit->start > $chunk->end);
	   push(@keepers, $hit);
	   $chunk_keepers = \@keepers;
       }
   }

   #add user defined labels
   foreach my $hit (@$chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }

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

   my $tmp = get_global_temp();
   my $command  = $blast;
   if ($command =~ /blasta$/) {
      symlink($blast, "$tmp/tblastx") if(! -e "$tmp/tblastx"); #handle blasta linking
      $command = "$tmp/tblastx";
      $command .= " $db $q_file B=10000 V=10000 E=$eval_blast";
      $command .= ($softmask) ? " wordmask=seg" : " filter=seg";
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
      $command .= " notes"; #suppress certain notes
      $command .= " novalidctxok"; #fixes failure related to short and masked sequence
      $command .= " nonnegok"; #fixes failure related to short and masked sequence
      $command .= " shortqueryok"; #fixes failure related to very short sequence
      $command .= ($org_type eq 'eukaryotic') ? "" : " kap";
      #$command .= " mformat=2"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /blastall$/) {
      $command .= " -p tblastx";
      $command .= " -d $db -i $q_file -b 10000 -v 10000 -e $eval_blast";
      $command .= " -z 1000";
      $command .= ($org_type eq 'eukaryotic') ? " -Y 500000000" : " -Y 20000000";
      $command .= " -a $cpus";	
      $command .= " -U";
      $command .= " -F T";
      $command .= " -I T";
      #$command .= " -m 8"; # remove for full report
      $command .= " -o $out_file";
   }
   elsif ($command =~ /tblastx$/) {
      $command .= " -db $db -query $q_file";
      $command .= " -num_alignments 10000 -num_descriptions 10000 -evalue $eval_blast";
      $command .= " -dbsize 1000";
      $command .= ($org_type eq 'eukaryotic') ? " -searchsp 500000000" : " -searchsp 20000000";
      $command .= " -num_threads $cpus";
      $command .= " -lcase_masking";
      $command .= " -seg yes";
      $command .= ($softmask) ? " -soft_masking true" : " -soft_masking false";
      $command .= " -show_gis";
      #$command .= " -outfmt 6"; # remove for full report
      $command .= " -out $out_file";
   }
   else{
      confess "ERROR: Must be a tblastx executable";  
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
sub clean_blast_hits{
    my $hits = shift;
    my $pcov = shift || 0;
    my $pid  = shift || 0;
    my $sig  = shift || 1000;
    my $con  = shift || 0; #contiguity flag

    my @keepers;

    foreach my $hit (@$hits){
	my $significance = $hit->significance();
	$significance = "1".$significance if($significance =~ /^e/);
	$significance =~ s/\.$//; # 0.

	#exonerate runs before clean because there may be a good exonerate hit
	#after extending.  So we want to keep the hit if there is an exonerate pair
	if($hit->{_keep}){ #set in polisher
	    push(@keepers, $hit);
	    next;
	}

	next unless ($significance < $sig);
	next unless ($hit->pAh > $pcov);
	next unless ($hit->hsp('best')->frac_identical() > $pid ||
		     $hit->frac_identical() > $pid);
	next unless (!$con || PhatHit_utils::is_contigous($hit));
	push(@keepers, $hit);
    }
    
    return \@keepers;
}
#-----------------------------------------------------------------------------
sub repeatmask {
   my $chunk        = shift;
   my $the_void     = shift;
   my $seq_id       = shift;
   my $model_org    = shift;
   my $RepeatMasker = shift;
   my $entry        = shift; #rmlib
   my $cpus         = shift;
   my $LOG = shift;

   #peal off label and build names for files to use and copy
   my ($rmlib, $label) = $entry =~ /^([^\:]+)\:?(.*)/;
   my ($db_n) = $rmlib =~ /([^\/]+)$/;
   $db_n = uri_escape($db_n, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.\+');

   my $chunk_number = $chunk->number();
   my $file_name = "$the_void/$seq_id\.$chunk_number";
   $file_name   .= ($rmlib) ? "\.$db_n\.specific" : "\.$model_org\.rb";
   my $o_file    = "$file_name\.out";

   my $q_length  = $chunk->parent_seq_length();
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
		   $cpus); # -no_low

   $chunk->erase_fasta_file();
   
   my $rm_chunk_keepers = Widget::RepeatMasker::parse($o_file, 
						      $seq_id, 
						      $q_length
						     );

   #delete other unneeded output files
   my @files = map {File::Glob::bsd_glob("$file_name.$_")} qw(ref tbl cat masked);
   unlink(@files);

   $LOG->add_entry("FINISHED", $o_file, ""); 
   
   PhatHit_utils::add_offset($rm_chunk_keepers, 
			     $chunk->offset(),
			    );

   #add user defined labels
   foreach my $hit (@$rm_chunk_keepers){
       $hit->{_label} = $label if($label);
       map{$_->{_label} = $label} $hit->hsps if($label);
   }
	
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

   if (-e $o_file) {
      print STDERR "re reading repeat masker report.\n" unless $main::quiet;
      print STDERR "$o_file\n" unless $main::quiet;
      return;
   }

   my $t_file;
   my $tmp = get_global_temp();
   my $command  = "cd $tmp; $RepeatMasker";

   if ($rmlib) {
      $command .= " $q_file -dir $dir -pa $cpus -lib $rmlib";
   }
   elsif($species eq 'simple'){
       my $lib = "$tmp/simple.lib";
       if(!-f $lib){
	   (my $tFH, $t_file) = tempfile(DIR => $tmp);
	   print $tFH ">(N)n#Dummy_repeat \@root  [S:25]\nnnnnnnnnnnnnnnnnnn\n";
	   close($tFH);
	   File::Copy::move($t_file, $lib);
       }
       $command .= " $q_file -dir $dir -pa $cpus -lib $lib";
   }
   else {
      $command .= " $q_file -species $species -dir $dir -pa $cpus";
   }
   $command .= " -nolow" if defined($no_low);
	
   my $w = new Widget::RepeatMasker();

   print STDERR "running  repeat masker.\n" unless $main::quiet;
   $w->run($command);
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
      $CTL_OPT{'maker_gff'} = '';
      $CTL_OPT{'est_pass'} = 0;
      $CTL_OPT{'altest_pass'} = 0;
      $CTL_OPT{'protein_pass'} = 0;
      $CTL_OPT{'rm_pass'} = 0;
      $CTL_OPT{'model_pass'} = 0;
      $CTL_OPT{'pred_pass'} = 0;
      $CTL_OPT{'other_pass'} = 0;
      $CTL_OPT{'est'} = '';
      $CTL_OPT{'est_reads'} = ''; #depricated
      $CTL_OPT{'altest'} = '';
      $CTL_OPT{'est_gff'} = '';
      $CTL_OPT{'altest_gff'} = '';
      $CTL_OPT{'protein'} = '';
      $CTL_OPT{'protein_gff'} = '';
      $CTL_OPT{'model_org'} = 'all';
      $CTL_OPT{'repeat_protein'} = Cwd::abs_path("$FindBin::Bin/../data/te_proteins.fasta") || '';
      $CTL_OPT{'rmlib'} = '';
      $CTL_OPT{'rm_gff'} = '';
      $CTL_OPT{'organism_type'} = 'eukaryotic';
      $CTL_OPT{'prok_rm'} = 0;
      $CTL_OPT{'predictor'} = '';
      $CTL_OPT{'predictor'} = 'model_gff' if($main::eva);
      $CTL_OPT{'snaphmm'} = '';
      $CTL_OPT{'gmhmm'} = '';
      $CTL_OPT{'gmhmm_e'} = '' if($main::server);
      $CTL_OPT{'gmhmm_p'} = '' if($main::server);
      $CTL_OPT{'augustus_species'} = '';
      $CTL_OPT{'fgenesh_par_file'} = '';
      $CTL_OPT{'snoscan_rrna'} = '';
      $CTL_OPT{'trna'} = 0;
      $CTL_OPT{'model_gff'} = '';
      $CTL_OPT{'pred_gff'} = '';
      $CTL_OPT{'est2genome'} = 0;
      $CTL_OPT{'altest2genome'} = 0;
      $CTL_OPT{'protein2genome'} = 0;
      $CTL_OPT{'other_gff'} = '';
      $CTL_OPT{'domain'} = '0';
      $CTL_OPT{'function'} = '0';
      $CTL_OPT{'short_name'} = '';
      $CTL_OPT{'snap_train'} = '0';
      $CTL_OPT{'alt_peptide'} = 'C';
      $CTL_OPT{'cpus'} = 1;
      $CTL_OPT{'cpus'} .= '=DISABLED' if($main::server);
      $CTL_OPT{'evaluate'} = 0;
      $CTL_OPT{'evaluate'} = 1 if($main::eva);
      $CTL_OPT{'max_dna_len'} = 100000;
      $CTL_OPT{'min_contig'} = 1;
      $CTL_OPT{'min_protein'} = 0;
      $CTL_OPT{'AED_threshold'} = 1;
      $CTL_OPT{'map_forward'} = 0;
      $CTL_OPT{'est_forward'} = 0; #only used to map old annotations to new assembly
      $CTL_OPT{'correct_est_fusion'} = 0;      
      $CTL_OPT{'always_complete'} = 0;
      $CTL_OPT{'pred_flank'} = 200;
      $CTL_OPT{'pred_stats'} = 0;
      $CTL_OPT{'keep_preds'} = 0;
      $CTL_OPT{'split_hit'} = 10000;
      $CTL_OPT{'softmask'} = 1;
      $CTL_OPT{'single_exon'} = 0;
      $CTL_OPT{'single_length'} = 250;
      $CTL_OPT{'tries'} = 2;
      $CTL_OPT{'retry'} = ''; #depricated
      $CTL_OPT{'clean_try'} = 0;
      $CTL_OPT{'TMP'} = '';
      $CTL_OPT{'TMP'} .= '=DISABLED' if($main::server);
      $CTL_OPT{'mpi_blastdb'} = ''; #hidden option
      $CTL_OPT{'run'} = ''; #hidden option
      $CTL_OPT{'unmask'} = 0;
      $CTL_OPT{'alt_splice'} = 0;
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
      $CTL_OPT{'datastore'} = 1;
      $CTL_OPT{'maker_v'} = version();
   }

   #maker_bopts
   if ($type eq 'all' || $type eq 'bopts') {
      $CTL_OPT{'blast_type'} = 'ncbi+';
      $CTL_OPT{'blast_type'} .= '=DISABLED' if($main::server);
      $CTL_OPT{'pcov_blastn'} = 0.80;
      $CTL_OPT{'pid_blastn'} = 0.85;
      $CTL_OPT{'eval_blastn'} = 1e-10;
      $CTL_OPT{'bit_blastn'} = 40;
      $CTL_OPT{'depth_blastn'} = 0;
      $CTL_OPT{'pcov_blastx'} = 0.50;
      $CTL_OPT{'pid_blastx'} = 0.40;
      $CTL_OPT{'eval_blastx'} = 1e-6;
      $CTL_OPT{'bit_blastx'} = 30;
      $CTL_OPT{'depth_blastx'} = 0;
      $CTL_OPT{'pcov_rm_blastx'} = 0.50;
      $CTL_OPT{'pid_rm_blastx'} = 0.40;
      $CTL_OPT{'eval_rm_blastx'} = 1e-6;
      $CTL_OPT{'bit_rm_blastx'} = 30;
      $CTL_OPT{'pcov_tblastx'} = 0.80;
      $CTL_OPT{'pid_tblastx'} = 0.85;
      $CTL_OPT{'eval_tblastx'} = 1e-10;
      $CTL_OPT{'bit_tblastx'} = 40;
      $CTL_OPT{'depth_tblastx'} = 0;
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
		  'makeblastdb',
		  'blastall',
		  'blasta',
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
		  'fathom',
		  'probuild',
		  'twinscan',
		  'qrna',
		  'jigsaw',
		  'tRNAscan-SE',
		  'snoscan'
		 );

      #get MAKER overriden exe locations
      my @all_alts = grep {-f $_ && -x $_} (File::Glob::bsd_glob("{$FindBin::Bin/../exe/*/*,$FindBin::Bin/../exe/*/bin/*}"));
      foreach my $exe (@exes) {
	  my @alts = grep {/\/$exe$/} @all_alts;
	  my $loc = shift @alts || File::Which::which($exe) || '';
	  if(! $loc && $exe eq 'blasta'){ #find blasta using blastx
	      ($loc) = (map {Cwd::abs_path($_)} grep {Cwd::abs_path($_) =~ /blasta$/} (File::Which::where('blastx')));
	  }
	  elsif($loc && $exe =~ /^.?blast[nx]$/ && Cwd::abs_path($loc) =~ /blasta$/){ #verify not blasta
	      ($loc) = (grep {Cwd::abs_path($_) !~ /blasta$/} (@alts, File::Which::where($exe)));
	  }
	  $CTL_OPT{$exe} = $loc || '';
      }
   }

   #server
   if ($type eq 'server') {
      $CTL_OPT{'DBI'} = 'SQLite';
      $CTL_OPT{'dbname'} = 'mwas_db';
      $CTL_OPT{'host'} = '';
      $CTL_OPT{'port'} = '';
      $CTL_OPT{'username'} = 'mwas';
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
      $CTL_OPT{'cgi_dir'} = '/data/var/www/cgi-bin' if(! -d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_dir'} = '' if(! -d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_dir'} .= '/mwas' if(-d $CTL_OPT{'cgi_dir'});
      $CTL_OPT{'cgi_web'} = '/cgi-bin/mwas';
      $CTL_OPT{'html_dir'} = '/var/www/html';
      $CTL_OPT{'html_dir'} = '/Library/WebServer/Documents' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} = '/var/www' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} = '/data/var/www/html' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} = '' if(! -d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_dir'} .= '/mwas' if(-d $CTL_OPT{'html_dir'});
      $CTL_OPT{'html_web'} = '/mwas';
      $CTL_OPT{'data_dir'} = '';
      $CTL_OPT{'data_dir'} = "/var/lib/mwas";
      $CTL_OPT{'web_address'} = 'http://'.Sys::Hostname::hostname();
      $CTL_OPT{'apache_user'} = '';
      $CTL_OPT{'apache_user'} = 'apache' if(@{[getpwnam('apache')]});
      $CTL_OPT{'apache_user'} = 'www' if(@{[getpwnam('www')]});
      $CTL_OPT{'apache_user'} = 'www-data' if(@{[getpwnam('www-data')]});
      $CTL_OPT{'font_file'} = '/usr/share/fonts/bitstream-vera/VeraMono.ttf';
      $CTL_OPT{'font_file'} = '/Library/Fonts/Verdana.ttf' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'font_file'} = '/Library/Fonts/Georgia.ttf' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'font_file'} = '/usr/share/fonts/truetype/freefont/FreeMono.ttf' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'font_file'} = '' if(! -f $CTL_OPT{'font_file'});
      $CTL_OPT{'soba_url'} = 'http://www.sequenceontology.org/cgi-bin/soba.cgi';
      $CTL_OPT{'JBROWSE_ROOT'} = "$FindBin::Bin/../exe/jbrowse";
      $CTL_OPT{'JBROWSE_ROOT'} = '/var/www/html/jbrowse' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'JBROWSE_ROOT'} = '/Library/WebServer/Documents/jbrowse' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'JBROWSE_ROOT'} = '/var/www/jbrowse' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'JBROWSE_ROOT'} = '/usr/local/gmod/jbrowse' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'JBROWSE_ROOT'} = '/data/var/www/jbrowse' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'JBROWSE_ROOT'} = '' if(! -d $CTL_OPT{'JBROWSE_ROOT'});
      $CTL_OPT{'GBROWSE_MASTER'} = '/etc/gbrowse/GBrowse.conf';
      $CTL_OPT{'GBROWSE_MASTER'} = '/etc/gbrowse2/GBrowse.conf' if(! -f $CTL_OPT{'GBROWSE_MASTER'});
      $CTL_OPT{'GBROWSE_MASTER'} = '' if(! -f $CTL_OPT{'GBROWSE_MASTER'});
      $CTL_OPT{'APOLLO_ROOT'} = "$FindBin::Bin/../exe/apollo";
      $CTL_OPT{'APOLLO_ROOT'} = $ENV{APOLLO_ROOT} if(! -d $CTL_OPT{'APOLLO_ROOT'} &&
						     $ENV{APOLLO_ROOT} &&
						     -d $ENV{APOLLO_ROOT});
      $CTL_OPT{'APOLLO_ROOT'} = (File::Which::which('apollo') || '') =~ /^([^\n]+)\/bin\/apollo\n?$/
	  if(! -d $CTL_OPT{'APOLLO_ROOT'});
      $CTL_OPT{'APOLLO_ROOT'} = '/usr/local/gmod/apollo' if(! -d $CTL_OPT{'APOLLO_ROOT'});
      $CTL_OPT{'APOLLO_ROOT'} = '' if(! -d $CTL_OPT{'APOLLO_ROOT'});
      $CTL_OPT{'ZOE'} = $ENV{'ZOE'} || '';
      $CTL_OPT{'AUGUSTUS_CONFIG_PATH'} = $ENV{'AUGUSTUS_CONFIG_PATH'} || '';
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
      $CTL_OPT{'est'}              = {'D. melanogaster : example cDNA' => "$server_ctl{data_dir}/maker/MWAS/../data/dpp_est.fasta",
				      'De novo/Legacy/Pass-through : example ESTs' => "$server_ctl{data_dir}/maker/MWAS/data/pyu-est.fasta",
				      'E. coli : example ESTs' => "$server_ctl{data_dir}/maker/MWAS/data/ecoli-est.fasta"};
      $CTL_OPT{'altest'}           = {};
      $CTL_OPT{'protein'}          = {'D. melanogaster : example proteins' => "$server_ctl{data_dir}/maker/MWAS/../data/dpp_protein.fasta",
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
      if($user_default->{menus}){
	  my %user_ctl = %{$user_default->{menus}};
	  while(my $key = each %user_ctl){
	      $CTL_OPT{$key} = {} if(! $CTL_OPT{$key});
	      %{$CTL_OPT{$key}} = %{$user_ctl{$key}} if($user_ctl{$key});
	  }
      }
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
	if(!$ENV{AUGUSTUS_CONFIG_PATH}){
	    $exes{augustus} = Cwd::abs_path($exes{augustus});
	    ($ENV{AUGUSTUS_CONFIG_PATH} = $exes{augustus}) =~ s/\/(bin|src)\/augustus/\/config/;
	}

	open(my $EXE, quotemeta($exes{augustus})." --species=help 2>&1|");
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
	$exes{snap} = Cwd::abs_path($exes{snap});
	$exes{snap} =~ s/[^\/]+$//;
	$exes{snap} = "$exes{snap}/HMM/";
    }

    if($exes{snap} && -d $exes{snap}){
	foreach my $file (grep {!/README/} File::Glob::bsd_glob("$exes{snap}/*")){
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
	$exes{gmhmme3} = Cwd::abs_path($exes{gmhmme3});
	$exes{gmhmme3} =~ s/[^\/]+$//;
	$exes{gmhmme3} = "$exes{gmhmme3}/HMM/";
 
	if(-d $exes{gmhmme3}){
	    foreach my $file (File::Glob::bsd_glob("$exes{gmhmme3}/*.mod")){
		my ($name) = $file =~ /([^\/]+)\.mod$/;
		$name =~ s/_/\. /;
		$name = ucfirst($name);
		my $value = Cwd::abs_path("$file");
		
		next if(!$name || !$value);
		
		$hmms{gmhmm_e}{$name} = $value;
	    }
	}
    }

    #find fgenesh HMMs
    if($exes{fgenesh}){
	$exes{fgenesh} = Cwd::abs_path($exes{fgenesh});
	$exes{fgenesh} =~ s/[^\/]+$//;
    
	if(-d $exes{fgenesh}){
	    foreach my $file (File::Glob::bsd_glob("$exes{fgenesh}/*")){
		my ($name) = $file =~ /([^\/]+)$/;
		my $value = Cwd::abs_path("$file");
		
		next if(!$name || !$value);
		next if($name =~ /\.zip$/);

		my $filesize = [stat($file)]->[7]; #size in bytes
		next unless(200000 <= $filesize && $filesize <= 400000); #rough size of all parameter files

		my $ok;
		open(IN, "< $file");
		while(my $line = <IN>){
		    if($line =! /Fgenesh_data/){
			$ok++;
			last;
		    }
		}
		close(IN);
		next unless($ok);

		#rough size of all parameter files
		$hmms{fgenesh_par_file}{$name} = $value;
	    }
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
	open (CTL, "< $ctlfile") or die "ERROR: Could not open the control file \"$ctlfile\".\n";
	
	while (my $line = <CTL>) {
	    chomp($line);
	   
	    if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:\=]+)[\:\=]([^\n\#]*)/) {
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
	open (CTL, "< $ctlfile") or die "ERROR: Could not open the control file \"$ctlfile\".\n";
	
	while (my $line = <CTL>) {
	    chomp($line);
	   
	    if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:\=]+)[\:\=]([^\n\#]*)/) {
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
	    elsif($line =~ /^\#\#Menus from Data::Dumper/){ #only non-standard format control file
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

	if($CTL_OPT{$key} eq '' && $key =~ /^(snap|fgenesh|augustus|gmhmme3|gmhmmp|tRNAscan|snoscan)$/){
	   $CTL_OPT{STAT}{$key} = 'DISABLED';
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
      open (CTL, "< $ctlfile") or die "ERROR: Could not open the control file \"$ctlfile\".\n";
	
      while (my $line = <CTL>) {
	 chomp($line);
	    
	 if ($line !~ /^[\#\s\t\n]/ && $line =~ /^([^\:\=]+)[\:\=]([^\n\#]*)/) {
	    my $key = $1;
	    my $value = $2;

	    #backwards compatability fix
	    $key = 'maker_gff' if($key eq 'genome_gff');

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
	       $CTL_OPT{retry} = '' if($key eq 'tries');
	    }
	    else {
	       warn "WARNING: Invalid option \'$key\' in control file $ctlfile\n\n";
	    }
	 }
      }
   }

   #--load command line options
   while (my $key = each %OPT){
       $CTL_OPT{$key} = $OPT{$key} if(defined $OPT{$key});
       $CTL_OPT{retry} = '' if($key eq 'tries');
   }

   #check organism type values
   if ($CTL_OPT{organism_type} !~ /^eukaryotic$|^prokaryotic$/) {
      $error .=  "ERROR: organism_type must be set to \'eukaryotic\' or \'prokaryotic\'.\n".
      "The value $CTL_OPT{organism_type} is invalid.\n\n";

      #set the default just to assist in populating remaining errors
      $CTL_OPT{organism_type} = 'eukaryotic';
   }

   #skip repeat masking command line option
   if ($OPT{R}) {
      $CTL_OPT{model_org} = '';
      $CTL_OPT{repeat_protein} = '';
      $CTL_OPT{rmlib} = '';
      $CTL_OPT{rm_gff} = '';
      $CTL_OPT{rm_pass} = 0;

      print STDERR "INFO: All repeat masking options will be skipped.\n\n" unless($main::qq);
   }

   if ($CTL_OPT{organism_type} eq 'prokaryotic' && ! $CTL_OPT{prok_rm}) {
      $CTL_OPT{model_org} = '';
      $CTL_OPT{repeat_protein} = '';
      $CTL_OPT{rmlib} = '';
      $CTL_OPT{rm_gff} = '';
      $CTL_OPT{rm_pass} = 0;

      print STDERR "WARNING: No repeats are expected in prokaryotic organisms. All repeat\n".
	           "masking options will be skipped. You can override this behavior by\n".
		   "setting prok_rm to 1 in the control files.\n\n" unless($main::qq);
   }

   #if no repeat masking options are set don't run masking dependent methods
   #i.e. unmasked predictions
   if ($CTL_OPT{model_org} eq '' &&
       $CTL_OPT{repeat_protein} eq '' &&
       $CTL_OPT{rmlib} eq '' &&
       $CTL_OPT{rm_gff} eq '' &&
       ($CTL_OPT{rm_pass} == 0 ||
	$CTL_OPT{maker_gff} eq '')
      ) {
       $CTL_OPT{unmask} = 0;
       $CTL_OPT{_no_mask} = 1; #no masking options found
   }

   #required for evaluator to work
   if($main::eva){
       $CTL_OPT{predictor} = 'model_gff';
       $CTL_OPT{model_pass} = 1;
       $CTL_OPT{evaluate} = 1;

       #evaluator only handles one or the other
       if($CTL_OPT{maker_gff} && $CTL_OPT{model_gff}){
	   $error .= "ERROR: In EVALUATOR you can have either models from a MAKER\n".
	       "produced GFF3 file or models from an external GFF3 file, but\n".
	       "not both (Check 'maker_gff' and 'model_gff' options)!\n\n";
       }
   }

   #parse predictor and error check
   $CTL_OPT{predictor} =~ s/\s+//g;
   my @predictors = split(',', $CTL_OPT{predictor});
   my @run;

   push(@run, 'snap') if($CTL_OPT{snaphmm});
   push(@run, 'genemark') if($CTL_OPT{gmhmm});
   push(@run, 'augustus') if($CTL_OPT{augustus_species});
   push(@run, 'fgenesh') if($CTL_OPT{fgenesh_par_file});
   push(@run, 'trnascan') if($CTL_OPT{trna});
   push(@run, 'snoscan') if($CTL_OPT{snoscan_rrna});

   if(! @predictors){ #build predictors if not provided
       push(@predictors, @run);
       push(@predictors, 'pred_gff') if($CTL_OPT{pred_gff} || ($CTL_OPT{maker_gff} && $CTL_OPT{pred_pass}));
       push(@predictors, 'model_gff') if($CTL_OPT{model_gff} || ($CTL_OPT{maker_gff} && $CTL_OPT{model_pass}));
   }
   push(@predictors, 'est2genome') if($CTL_OPT{est2genome} && ! $main::eva);
   push(@predictors, 'protein2genome') if($CTL_OPT{protein2genome} && ! $main::eva);

   $CTL_OPT{_predictor} = {}; #temporary hash
   $CTL_OPT{_run} = {}; #temporary hash
   foreach my $p (@predictors) {
       if ($p !~ /^snap$|^augustus$|^est2genome$|^protein2genome$|^fgenesh$/ &&
	   $p !~ /^genemark$|^model_gff$|^pred_gff$|^trnascan$|^snoscan$/
	   ) {
	   $error .= "FATAL: Invalid predictor defined: $p\n".
	       "Valid entries are: est2genome, model_gff, pred_gff,\n".
	       "snap, genemark, augustus, and fgenesh\n\n";
	   next;
       }
       if($CTL_OPT{organism_type} eq 'prokaryotic' &&
	  $p =~ /^snap$|^augustus$|^fgenesh$/
	  ){
	   warn "WARNING: the predictor $p does not support prokaryotic organisms\n".
	       "and will be ignored.\n\n";
	   next;
       }

       $CTL_OPT{_predictor}{$p}++;
       $CTL_OPT{_run}{$p}++ unless($p =~ /est2genome|protein2genome|model_gff|pred_gff/);
   }
   $CTL_OPT{_predictor} = [keys %{$CTL_OPT{_predictor}}]; #convert to array
   $CTL_OPT{predictor} = join(",", @{$CTL_OPT{_predictor}}); #reset value for log
   $CTL_OPT{_run} = [keys %{$CTL_OPT{_run}}]; #convert to array
   $CTL_OPT{run} = join(",", @{$CTL_OPT{_run}}); #reset value for log

   #check blast type validity and related values (NCBI vs. WUBLAST)
   $CTL_OPT{blast_type} = lc($CTL_OPT{blast_type});
   if ($CTL_OPT{blast_type} !~ /^wublast$|^ncbi$|^ncbi\+$/) {
      warn "WARNING: blast_type must be set to 'wublast', 'ncbi', or 'ncbi+'.\n",
      "The value $CTL_OPT{blast_type} is invalid.\n",
      "This will now be reset to the default 'ncbi+'.\n\n" unless($main::qq);
      
      $CTL_OPT{blast_type} = 'ncbi+';
   }
   
   if (($CTL_OPT{blast_type} =~ /^wublast$/ && ! -f $CTL_OPT{xdformat}) ||
       ($CTL_OPT{blast_type} =~ /^ncbi$/ && ! -f $CTL_OPT{formatdb}) ||
       ($CTL_OPT{blast_type} =~ /^ncbi\+$/ && ! -f $CTL_OPT{makeblastdb})
      ) {
       my $new;
       if(-f $CTL_OPT{makeblastdb}){
	   $new = 'ncbi+';
       }
       elsif(-f $CTL_OPT{formatdb}){
	   $new = 'ncbi';
       }
       elsif(-f $CTL_OPT{xdformat}){
	   $new = 'wublast';
       }

       warn "WARNING: blast_type is set to '$CTL_OPT{blast_type}' but executables cannot be located\n" unless($main::qq);
       warn "The blast_type '$new' will be used instead.\n\n" if($new && !$main::qq);
       $error .= "ERROR: Please provide a valid locaction for a BLAST algorithm in the control files.\n\n" if(!$new);

       $CTL_OPT{blast_type} = $new;
   }
   
   #use standard value to refer to both NCBI and WUBLAST
   if ($CTL_OPT{blast_type} =~ /^wublast$/i) {
      $CTL_OPT{_formater} = $CTL_OPT{xdformat};
      $CTL_OPT{_blastn} = $CTL_OPT{blasta};
      $CTL_OPT{_blastx} = $CTL_OPT{blasta};
      $CTL_OPT{_tblastx} = $CTL_OPT{blasta};
   }
   elsif ($CTL_OPT{blast_type} =~ /^ncbi$/i) {
      $CTL_OPT{_formater} = $CTL_OPT{formatdb};
      $CTL_OPT{_blastn} = $CTL_OPT{blastall};
      $CTL_OPT{_blastx} = $CTL_OPT{blastall};
      $CTL_OPT{_tblastx} = $CTL_OPT{blastall};
   }
   elsif ($CTL_OPT{blast_type} =~ /^ncbi\+$/i) {
      $CTL_OPT{_formater} = $CTL_OPT{makeblastdb};
      $CTL_OPT{_blastn} = $CTL_OPT{blastn};
      $CTL_OPT{_blastx} = $CTL_OPT{blastx};
      $CTL_OPT{_tblastx} = $CTL_OPT{tblastx};
   }
   
   #--validate existence of required values from control files
   my @infiles;
   if($CTL_OPT{blast_type} =~ /^wublast$/i){
       push (@infiles, 'blasta', 'xdformat') if($CTL_OPT{est});
       push (@infiles, 'blasta', 'xdformat') if($CTL_OPT{protein}); 
       push (@infiles, 'blasta', 'xdformat') if($CTL_OPT{repeat_protein}); 
       push (@infiles, 'blasta', 'xdformat') if($CTL_OPT{altest});
   }
   elsif($CTL_OPT{blast_type} =~ /^ncbi$/i){
       push (@infiles, 'blastall', 'formatdb') if($CTL_OPT{est});
       push (@infiles, 'blastall', 'formatdb') if($CTL_OPT{protein}); 
       push (@infiles, 'blastall', 'formatdb') if($CTL_OPT{repeat_protein}); 
       push (@infiles, 'blastall', 'formatdb') if($CTL_OPT{altest});
   }
   elsif($CTL_OPT{blast_type} =~ /^ncbi\+$/i){
       push (@infiles, 'blastn', 'makeblastdb') if($CTL_OPT{est});
       push (@infiles, 'blastx', 'makeblastdb') if($CTL_OPT{protein}); 
       push (@infiles, 'blastx', 'makeblastdb') if($CTL_OPT{repeat_protein}); 
       push (@infiles, 'tblastx', 'makeblastdb') if($CTL_OPT{altest});
   }

   push (@infiles, 'genome');
   push (@infiles, 'est') if($CTL_OPT{est}); 
   push (@infiles, 'protein') if($CTL_OPT{protein}); 
   push (@infiles, 'altest') if($CTL_OPT{altest}); 
   push (@infiles, 'probuild') if (grep {/genemark/} @{$CTL_OPT{_run}});
   push (@infiles, 'fathom') if ($CTL_OPT{enable_fathom} && $CTL_OPT{evaluate});
   push (@infiles, 'est_gff') if($CTL_OPT{est_gff});
   push (@infiles, 'protein_gff') if($CTL_OPT{protein_gff});
   push (@infiles, 'maker_gff') if($CTL_OPT{maker_gff});
   push (@infiles, 'pred_gff') if($CTL_OPT{pred_gff});
   push (@infiles, 'model_gff') if ($CTL_OPT{model_gff});
   push (@infiles, 'snap') if (grep {/snap/} @{$CTL_OPT{_run}});
   push (@infiles, 'augustus') if (grep {/augustus/} @{$CTL_OPT{_run}}); 
   push (@infiles, 'fgenesh') if (grep {/fgenesh/} @{$CTL_OPT{_run}});
   push (@infiles, 'repeat_protein') if ($CTL_OPT{repeat_protein});
   push (@infiles, 'RepeatMasker') if($CTL_OPT{rmlib});
   push (@infiles, 'RepeatMasker') if($CTL_OPT{model_org});
   push (@infiles, 'rm_gff') if($CTL_OPT{rm_gff});
   push (@infiles, 'rmlib') if ($CTL_OPT{rmlib});
   push (@infiles, 'tRNAscan-SE') if (grep {/trnascan/} @{$CTL_OPT{_run}});
   push (@infiles, 'snoscan') if (grep {/snoscan/} @{$CTL_OPT{_run}});
   
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
      if (! $CTL_OPT{$in}) {
	 $error .= "ERROR: You have failed to provide a value for \'$in\' in the control files.\n\n";
	 next;
      }
      else{
	  my @files = split(/\,/, $CTL_OPT{$in}); #handle comma separated list
	  foreach my $file (@files){
	      $file =~ s/^\s+|\s+$//g;
	      $file =~ s/\s+\:|\:\s+$//g;
	  }
	  my %uniq = map {/^([^\:]+)\:?(.*)?/} @files;
	  my @non = grep {! -f $_} keys %uniq;
	  $error .= "ERROR: The \'$in\' file[s] ".join(', ', @non)." do not exist.\n".
	      "Please check settings in the control files.\n\n"if(@non);

	  @files = map {($uniq{$_}) ? s_abs_path($_).":$uniq{$_}" : s_abs_path($_)} keys %uniq;
	  $CTL_OPT{$in} = join(',', @files); #order entries for logging
      }
   }

   #--error check sometimes required values
   if ((grep {/model_gff/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{model_gff} &&
       (!$CTL_OPT{maker_gff} || !$CTL_OPT{model_pass})
      ){
       $error .= "ERROR: You must provide gene models in a GFF3 file to use model_gff as a predictor.\n\n";
   }
   if ((grep {/pred_gff/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{pred_gff} &&
       (!$CTL_OPT{maker_gff} || !$CTL_OPT{pred_pass})
      ){
       $error .= "ERROR: You must provide gene predictions in a GFF3 file to use pred_gff as a predictor.\n\n";
   }
   if ((grep {/est2genome/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{est} &&
       !$CTL_OPT{altest} &&
       !$CTL_OPT{est_gff} &&
       !$CTL_OPT{altest_gff} &&
       (!$CTL_OPT{est_pass} || !$CTL_OPT{maker_gff})
      ){
       $error .= "ERROR: You must provide some form of EST evidence to use est2genome as a predictor.\n\n";
   } 
   if ((grep {/protein2genome/} @{$CTL_OPT{_predictor}}) &&
       !$CTL_OPT{protein} &&
       !$CTL_OPT{protein_gff} &&
       (!$CTL_OPT{protein_pass} || !$CTL_OPT{maker_gff})
       ){
       $error .= "ERROR: You must provide some form of protein evidence to use protein2genome as a predictor.\n\n";
   }
			       
   #--error check that values are meaningful
   if (grep {/augustus/} @{$CTL_OPT{_run}}){
       if (! $ENV{AUGUSTUS_CONFIG_PATH} || ! -f "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg") {
	   #try and find it
	   my ($path) = Cwd::abs_path($CTL_OPT{augustus});
	   $path =~ s/bin\/augustus$/config/;
	   $ENV{AUGUSTUS_CONFIG_PATH} = $path;
	   
	   if(-f $CTL_OPT{augustus} && ! -f "$ENV{AUGUSTUS_CONFIG_PATH}/extrinsic/extrinsic.MPE.cfg"){
	       $error .= "ERROR: The environmental variable AUGUSTUS_CONFIG_PATH has not been set\n".
		   "or is not set correctly. Please set this in your profile per Augustus\n".
		   "installation instructions then try running MAKER again.\n\n";
	   }
       }

       if(!$CTL_OPT{augustus_species}) {
	   $error .= "ERROR: There is no species specified for Augustus (augustus_species).\n\n";
       }
   }
   if ($CTL_OPT{blast_type} =~ /^wublast$/i && -f $CTL_OPT{blasta}) {
       (my $base = Cwd::abs_path($CTL_OPT{blasta})) =~ s/\/blasta$//;

       #set both variable types
       my $matrix = $ENV{WUBLASTMAT} || $ENV{BLASTMAT};
       my $filter = $ENV{WUBLASTFILTER} || $ENV{BLASTFILTER};

       if(!$matrix || ! -d $matrix){
	   $matrix = "$base/matrix";
       }
       if(!$filter || ! -d $filter){
	   $filter = "$base/filter";
       }

       $error.= "ERROR: You must properly set the BLASTFILTER/WUBLASTFILTER and\n".
	        "BLASTMAT/WUBLASTMAT environmental variables or WUBLAST will not\n".
	        "function\n" if(! -d $matrix || ! -d $filter);

       $ENV{BLASTMAT} = $ENV{WUBLASTMAT} = $matrix;
       $ENV{BLASTFILTER} = $ENV{WUBLASTFILTER} = $filter;
   }
   if (grep {/^snap$|^fathom$/} @{$CTL_OPT{_run}}){
       if(! -d $ENV{ZOE}){
	   #try and find it
	   my ($path) = Cwd::abs_path($CTL_OPT{snap});
	   $path =~ s/snap$//;
	   $ENV{ZOE} = $path;
       }

       if(! $CTL_OPT{snaphmm}) {
	   $error .= "ERROR: There is no HMM specified for for SNAP/Fathom (snaphmm).\n\n";
       }
       else{
	   my @files = split(/\,/, $CTL_OPT{snaphmm});
	   my %uniq;
	   @files = grep {! $uniq{$_}++} @files;
	   my @non = grep {/^([^\:]+)/ && ! -f $1 && ! -f $ENV{ZOE}."/HMM/".$1} @files;
	   $error .= "ERROR: The HMM file[s] provided for snaphmm do not exist:\n".
	       "\t".join("\n\t", @non)."\n\n" if(@non);
       }
   }
   if (grep {/^genemark$/} @{$CTL_OPT{_run}}) {
       if(! $CTL_OPT{gmhmm}) {
	   $error .= "ERROR: There is no HMM specified for for GeneMark (gmhmm).\n\n";
       }
       else{
	   my @files = split(/\,/, $CTL_OPT{gmhmm});
	   my %uniq;
	   @files = grep {! $uniq{$_}++} @files;
	   my @non = grep {/^([^\:]+)/ && ! -f $1} @files;
	   $error .= "ERROR: The HMM file[s] provided for gmhmm do not exist:\n".
	       "\t".join("\n\t", @non)."\n\n" if(@non);
       }
   }
   if (grep {/^fgenesh$/} @{$CTL_OPT{_run}}) {
      if (! $CTL_OPT{fgenesh_par_file}) {
	  $error .= "ERROR: There is no parameter file secified for FgenesH (fgenesh_par_file)\n\n";
      }
      else {
	  my @files = split(/\,/, $CTL_OPT{gmhmm});
	  my %uniq;
	  @files = grep {! $uniq{$_}++} @files;
	  my @non = grep {/^([^\:]+)/ && ! -f $1} @files;
           $error .= "ERROR: The parameter file[s] provided for fgenesh_par_file do not exist:\n".
               "\t".join("\n\t", @non)."\n\n" if(@non);
      }
   }
   if ($CTL_OPT{max_dna_len} < 50000) {
      warn "WARNING: \'max_dna_len\' is set too low.  The minimum value permited is 50,000.\n".
      "max_dna_len will be reset to 50,000\n\n";
      $CTL_OPT{max_dna_len} = 50000;
   }
   if ($CTL_OPT{split_hit} < 50 || $CTL_OPT{organism_type} eq 'prokaryotic') {
      $CTL_OPT{split_hit} = 50; #important or hits will not be merged across chunk junctions
   }
   if ($CTL_OPT{split_hit} > $CTL_OPT{max_dna_len}/3){
       $error .= "ERROR: split_hit cannot me more than 1/3 the value of max_dna_len\n".
   	         "Try raising max_dna_len or lowering split_hit\n";
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
   #catch depricated parameter
   if ($CTL_OPT{retry} ne ''){
       $CTL_OPT{tries} = 1 + $CTL_OPT{retry};
       delete($CTL_OPT{retry});
   }
   else{
       delete($CTL_OPT{retry});
   }
   if ($CTL_OPT{tries} < 0) {
      warn "WARNING: \'tries\' must be set to 0 or greater.\n".
	   "It will now be set to 0\n\n";
      $CTL_OPT{tries} = 1;
   }
   if($CTL_OPT{TMP} && ! -d $CTL_OPT{TMP}){
       $error .= "ERROR: The TMP value \'$CTL_OPT{TMP}\' is not a directory or does not exist.\n\n";
   }
   my $TMPDIR = $CTL_OPT{TMP} || File::Spec->tmpdir();
   if(is_NFS_mount($TMPDIR)){
       if($CTL_OPT{ignore_nfs_tmp}){
	   warn "WARNING: Temporary directory set to an NFS location.\n".
	        "TMP=$TMPDIR\n".
		"The temporary directory in MAKER is specifically for\n".
		"operations that are not NFS-safe, but you have chosen\n".
		"to ignore this error. If you experience seemly random\n".
		"freezing and failures, the TMP directory is the cause.\n\n";
       }
       else{
	   $error .= "ERROR: Temporary directory set to an NFS location.\n".
	             "TMP=$TMPDIR\n".
		     "The temporary directory in MAKER is specifically for\n".
		     "operations that are not NFS-safe. You must set TMP\n".
		     "to a locally mounted directory such as /tmp or add\n".
		     "--ignore_nfs_tmp to the maker command line to\n".
		     "override this error message.\n\n";
       }
   }   
   if($main::eva && $CTL_OPT{maker_gff} && $CTL_OPT{model_gff}){ #only for evaluator
       $error .= "You can only specify a GFF3 file for maker_gff or model_gff no both!!\n\n";
   }
   if($CTL_OPT{out_name} && $CTL_OPT{out_name} =~ /\//){
       $error .= "The 'base' option is an ID used in generating names for\n".
	         "multiple files and directories. It is not a directory path.\n".
		 "You must remove all '/' characters from '$CTL_OPT{out_name}'\n\n";
   }

   #--if just parsing without error check stop here
   return %CTL_OPT if($OPT{parse});
   
   #--report errors
   die $error if ($error);   

   #--check genome file for fasta entries
   my $iterator = new Iterator::Any( -fasta => $CTL_OPT{genome},
				     -gff => $CTL_OPT{genome}
				   );
   
   unless($iterator->nextDef) {
      my $genome = $CTL_OPT{genome};
      die "ERROR:  The file $genome contains no fasta entries\n\n";
   }

   #--decide whether to force datastore, datastore will already be defined if selected by user 
   $CTL_OPT{datastore} = 1 if(! defined $CTL_OPT{datastore}); #on by default

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
   $CTL_OPT{CWD} = Cwd::cwd();
   if(! $CTL_OPT{out_name}){
       ($CTL_OPT{out_name}) = $genome =~ /([^\/]+)$/;
       $CTL_OPT{out_name} =~ s/\.[^\.]+$//;
   }
   if(! $CTL_OPT{out_base}){
      $CTL_OPT{out_base} = $CTL_OPT{CWD}."/$CTL_OPT{out_name}.maker.output";
   }
   if(! $CTL_OPT{mpi_blastdb}){
      $CTL_OPT{mpi_blastdb} = $CTL_OPT{out_base}."/mpi_blastdb";
   }
   else{
       $CTL_OPT{mpi_blastdb} = s_abs_path($CTL_OPT{mpi_blastdb});
   }
   mkdir($CTL_OPT{out_base}) if(! -d $CTL_OPT{out_base});
   confess "ERROR: Could not build output directory $CTL_OPT{out_base}\n"
        if(! -d $CTL_OPT{out_base});

   #--make sure repbase is installed
   if($CTL_OPT{model_org}){
       my $exe = Cwd::abs_path($CTL_OPT{RepeatMasker});
       my ($lib) = $exe =~ /(.*\/)RepeatMasker$/;
       die "ERROR: Could not determine if RepBase is installed\n" if(! $lib);

       $lib .= "Libraries/RepeatMaskerLib.embl";
       die "ERROR: Could not determine if RepBase is installed\n" if(! -f $lib);

       open(my $IN, "< $lib");
       my $rb_flag;
       for(my $i = 0; $i < 20; $i++){
           my $line = <$IN>;
           if($line =~ /RELEASE \d+(\-min)?\;/){
	       $rb_flag = ($1 && $1 eq '-min') ? 0 : 1;
	       last;
           }
       }
       close($IN);

       if(! $rb_flag){
	   warn "WARNING: RepBase is not installed for RepeatMasker. This limits\n".
	       "RepeatMasker's functionality and makes the model_org option in the\n".
	       "control files virtually meaningless. MAKER will now reconfigure\n".
	       "for simple repeat masking only.\n";
	   $CTL_OPT{model_org} = 'simple';
       }
   }

   #--set an initialization lock so steps are locked to a single process
   #--take extra steps since lock is vulnerable to race conditions
   my $i_lock; #init lock, it is only a temporary blocking lock
   my $is_NFS = is_NFS_mount($CTL_OPT{out_base});
   $CTL_OPT{_is_NFS} = $is_NFS;
   if($CTL_OPT{_no_lock}){
       warn "WARNING: You have chosen to turn locking off which may create\n".
	    "race conditions if running in parallel.  You have been warned.\n\n";
       $File::NFSLock::TYPE = 0;
   }
   elsif($is_NFS == 2){ #FhGFS type global storage
       warn "WARNING: The output directory if FhGFS which does not support\n".
	    "hardlinks and may have broken posix locks, but I will try\n".
	    "posix locking anyway (hopefully no race conditions emerge).\n\n";
       $File::NFSLock::TYPE = 'POSIX';
   }
   my $tries = 0;
   while(!$i_lock){
       $tries++;
       $i_lock = new File::NFSLock($CTL_OPT{out_base}."/init_lock", 'EX', 150, 120);
       my $is_mine = $i_lock->still_mine if($i_lock);
       if(!$is_mine && $tries >= 5){ #try up to 5 times
	   my $err = "ERROR: Cannot get initialization lock.\n".
	             "If you are running maker in parallel or via MPI\n".
	             "You may be facing a race condition.\n\n";
	   $err .= "The output directory if NFS and may be a factor\n\n" if($is_NFS == 1);
	   die $err;
       }
       elsif(!$is_mine){ #sleep with random wait on failure
	   undef $i_lock;
	   warn "WARNING: Could not get initialization lock. Trying Again...\n";
	   my $time = int(rand(20))+1;
	   sleep $time;
       }
       elsif($is_mine){ #success
	   last;
       }
   }

   #--check if MAKER is already running and lock the directory
   #lock must be global or it will be destroyed outside of block
   unless(($LOCK = new File::NFSLock($CTL_OPT{out_base}."/gi_lock", 'SH', 40, 40)) && $LOCK->maintain(30)){
       die "ERROR: The directory is locked.  Perhaps by an instance of MAKER.\n\n";
   }

   #compare current control files to logged files
   my $app = ($main::eva) ? "eval" : "maker";
   my @ctl_logs = ($CTL_OPT{out_base}."/$app\_opts.log",
		   $CTL_OPT{out_base}."/$app\_bopts.log",
		   $CTL_OPT{out_base}."/$app\_exe.log"
		   );

   #check who else is also sharing the lock and if running same settings
   $CTL_OPT{_step} = $LOCK->owners() || 1;
   $CTL_OPT{_shared_id} = $LOCK->shared_id();
   if($CTL_OPT{_step} == 1){ #I am only/first holder of the lock
       #control files already exist see if I should remove the datastore log
       if(-f $ctl_logs[0] && -f $ctl_logs[1] && -f $ctl_logs[2]){
	   $CTL_OPT{_resume}++ if(runlog::are_same_opts(\%CTL_OPT, \@ctl_logs));
	   $CTL_OPT{_resume} = 0 if($CTL_OPT{force} || $CTL_OPT{again});
       }

       #log the control files
       generate_control_files($CTL_OPT{out_base}, 'all', \%CTL_OPT, 1);

       $i_lock->unlock if($i_lock); #release init lock
   }
   else{
       $i_lock->unlock if($i_lock); #release init lock

       unless (-f $ctl_logs[0] && -f $ctl_logs[1] && -f $ctl_logs[2]){
	   $LOCK->unlock if($LOCK);
	   confess "ERROR: Could not query control option logs\n\n";
       }

       #should be same
       if(! runlog::are_same_opts(\%CTL_OPT, \@ctl_logs)){
	   $LOCK->unlock if($LOCK);
	   die "ERROR: Cannot start process. MAKER/EVALUATOR already running\n".
	       "with different settings in this same directory.\n\n";
       }
       else{#start a second MAKER process, but give a warning
	   warn "WARNING: Multiple MAKER processes have been started in the\n".
	        "same directory.\n\n" unless($main::qq);

	   $CTL_OPT{_multi_chpc}++; #multi process flag
       }
   }

   #--exit with status of 0 if just checking control files with -check fla
   exit(0) if($OPT{check});

   #---set up blast databases and indexes for analyisis
   $CTL_OPT{_mpi_size} = $mpi_size;

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
       open (OUT, "> $dir/$app\_opts.$ext") or
	   die "ERROR: Could not create $dir/$app\_opts.$ext\n";
       print OUT "#-----Genome (these are always required)\n";
       print OUT "genome=$O{genome} #genome sequence (fasta file or fasta embeded in GFF3 file)\n";
       print OUT "organism_type=$O{organism_type} #eukaryotic or prokaryotic. Default is eukaryotic\n";
       print OUT "\n";
       print OUT "#-----Re-annotation Using MAKER Derived GFF3\n" if(!$ev);
       print OUT "#-----MAKER Derived GFF3 Annotations to Evaluate\n" if($ev);
       print OUT "maker_gff=$O{maker_gff} #MAKER derived GFF3 file\n";
       print OUT "est_pass=$O{est_pass} #use ESTs in maker_gff: 1 = yes, 0 = no\n";
       print OUT "altest_pass=$O{altest_pass} #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no\n";
       print OUT "protein_pass=$O{protein_pass} #use protein alignments in maker_gff: 1 = yes, 0 = no\n";
       print OUT "rm_pass=$O{rm_pass} #use repeats in maker_gff: 1 = yes, 0 = no\n";
       print OUT "model_pass=$O{model_pass} #use gene models in maker_gff: 1 = yes, 0 = no\n" if(!$ev);
       print OUT "pred_pass=$O{pred_pass} #use ab-initio predictions in maker_gff: 1 = yes, 0 = no\n";
       print OUT "other_pass=$O{other_pass} #passthrough anyything else in maker_gff: 1 = yes, 0 = no\n" if(!$ev);
       print OUT "\n";
       print OUT "#-----External GFF3 Annotations to Evaluate\n" if($ev);
       print OUT "model_gff=$O{model_gff} #gene models from an external GFF3 file\n" if($ev);
       print OUT "\n"if($ev);
       print OUT "#-----EST Evidence (for best results provide a file for at least one)\n";
       print OUT "est=$O{est} #set of ESTs or assembled mRNA-seq in fasta format\n";
       print OUT "altest=$O{altest} #EST/cDNA sequence file in fasta format from an alternate organism\n";
       print OUT "est_gff=$O{est_gff} #aligned ESTs or mRNA-seq from an external GFF3 file\n";
       print OUT "altest_gff=$O{altest_gff} #aligned ESTs from a closly relate species in GFF3 format\n";
       print OUT "\n";
       print OUT "#-----Protein Homology Evidence (for best results provide a file for at least one)\n";
       print OUT "protein=$O{protein}  #protein sequence file in fasta format (i.e. from mutiple oransisms)\n";
       print OUT "protein_gff=$O{protein_gff}  #aligned protein homology evidence from an external GFF3 file\n";
       print OUT "\n";
       print OUT "#-----Repeat Masking (leave values blank to skip repeat masking)\n";
       print OUT "model_org=$O{model_org} #select a model organism for RepBase masking in RepeatMasker\n";
       print OUT "rmlib=$O{rmlib} #provide an organism specific repeat library in fasta format for RepeatMasker\n";
       print OUT "repeat_protein=$O{repeat_protein} #provide a fasta file of transposable element proteins for RepeatRunner\n";
       print OUT "rm_gff=$O{rm_gff} #pre-identified repeat elements from an external GFF3 file\n";
       print OUT "prok_rm=$O{prok_rm} #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no\n";
       print OUT "softmask=$O{softmask} #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)\n";
       print OUT "\n";
       print OUT "#-----Gene Prediction\n" if(!$ev);
       print OUT "#-----EVALUATOR Ab-Initio Comparison Options\n" if($ev);
       print OUT "snaphmm=$O{snaphmm} #SNAP HMM file\n";
       print OUT "gmhmm=$O{gmhmm} #GeneMark HMM file\n";
       print OUT "augustus_species=$O{augustus_species} #Augustus gene prediction species model\n";
       print OUT "fgenesh_par_file=$O{fgenesh_par_file} #FGENESH parameter file\n";
       print OUT "pred_gff=$O{pred_gff} #ab-initio predictions from an external GFF3 file\n";
       print OUT "model_gff=$O{model_gff} #annotated gene models from an external GFF3 file (annotation pass-through)\n" if(!$ev);
       print OUT "est2genome=$O{est2genome} #infer gene predictions directly from ESTs, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "protein2genome=$O{protein2genome} #infer predictions from protein homology, 1 = yes, 0 = no\n"  if(!$ev);
       print OUT "trna=$O{trna} #find tRNAs with tRNAscan, 1 = yes, 0 = no\n";
       print OUT "snoscan_rrna=$O{snoscan_rrna} #rRNA file to have Snoscan find snoRNAs\n";
       print OUT "unmask=$O{unmask} #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no\n";
       print OUT "\n";
       print OUT "#-----Other Annotation Feature Types (features MAKER doesn't recognize)\n" if(!$ev);
       print OUT "other_gff=$O{other_gff} #extra features to pass-through to final MAKER generated GFF3 file\n" if(!$ev);
       print OUT "\n" if(!$ev);
       print OUT "#-----External Application Behavior Options\n";
       print OUT "alt_peptide=$O{alt_peptide} #amino acid used to replace non-standard amino acids in BLAST databases\n";
       print OUT "cpus=$O{cpus} #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)\n";
       print OUT "\n";
       print OUT "#-----MAKER Behavior Options\n";
       print OUT "max_dna_len=$O{max_dna_len} #length for dividing up contigs into chunks (increases/decreases memory usage)\n";
       print OUT "min_contig=$O{min_contig} #skip genome contigs below this length (under 10kb are often useless)\n" if(!$ev);
       print OUT "\n";
       print OUT "pred_flank=$O{pred_flank} #flank for extending evidence clusters sent to gene predictors\n";
       print OUT "pred_stats=$O{pred_stats} #report AED and QI statistics for all predictions as well as models\n";
       print OUT "AED_threshold=$O{AED_threshold} #Maximum Annotation Edit Distance allowed (bound by 0 and 1)\n" if(!$ev);
       print OUT "min_protein=$O{min_protein} #require at least this many amino acids in predicted proteins\n" if(!$ev);
       print OUT "alt_splice=$O{alt_splice} #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "always_complete=$O{always_complete} #extra steps to force start and stop codons, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "map_forward=$O{map_forward} #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no\n" if(!$ev);
       print OUT "est_forward=$O{est_forward} #reserve flag for map2assembly\n" if($O{est_forward});
       print OUT "keep_preds=$O{keep_preds} #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)\n" if(!$ev);
       print OUT "\n";
       print OUT "split_hit=$O{split_hit} #length for the splitting of hits (expected max intron size for evidence alignments)\n";
       print OUT "single_exon=$O{single_exon} #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no\n";
       print OUT "single_length=$O{single_length} #min length required for single exon ESTs if \'single_exon\ is enabled'\n";
       print OUT "correct_est_fusion=$O{correct_est_fusion} #limits use of ESTs in annotation to avoid fusion genes\n";
       print OUT "\n";
       print OUT "tries=$O{tries} #number of times to try a contig if there is a failure for some reason\n";
       print OUT "clean_try=$O{clean_try} #remove all data from previous run before retrying, 1 = yes, 0 = no\n";
       print OUT "clean_up=$O{clean_up} #removes theVoid directory with individual analysis files, 1 = yes, 0 = no\n";
       print OUT "TMP=$O{TMP} #specify a directory other than the system default temporary directory for temporary files\n";
       print OUT "\n" if($ev);
       print OUT "#-----EVALUATOR Control Options (ignore these)\n" if($ev);
       #print OUT "evaluate=$O{evaluate} #run EVALUATOR on all annotations (very experimental), 1 = yes, 0 = no\n";
       print OUT "side_thre=$O{side_thre}\n" if($ev);
       print OUT "eva_window_size=$O{eva_window_size}\n" if($ev);
       print OUT "eva_split_hit=$O{eva_split_hit}\n" if($ev);
       print OUT "eva_hspmax=$O{eva_hspmax}\n" if($ev);
       print OUT "eva_gspmax=$O{eva_gspmax}\n" if($ev);
       print OUT "enable_fathom=$O{enable_fathom}\n" if($ev);
       close (OUT);
   }
    
   #--build bopts.ctl file
   if($type eq 'all' || $type eq 'bopts'){
       open (OUT, "> $dir/$app\_bopts.$ext") or
	   die "ERROR: Could not create $dir/$app\_bopts.$ext\n";
       print OUT "#-----BLAST and Exonerate Statistics Thresholds\n";
       print OUT "blast_type=$O{blast_type} #set to 'ncbi+', 'ncbi' or 'wublast'\n";
       print OUT "\n";
       print OUT "pcov_blastn=$O{pcov_blastn} #Blastn Percent Coverage Threhold EST-Genome Alignments\n";
       print OUT "pid_blastn=$O{pid_blastn} #Blastn Percent Identity Threshold EST-Genome Aligments\n";
       print OUT "eval_blastn=$O{eval_blastn} #Blastn eval cutoff\n";
       print OUT "bit_blastn=$O{bit_blastn} #Blastn bit cutoff\n";
       print OUT "depth_blastn=$O{depth_blastn} #Blastn depth cutoff (0 to disable cutoff)\n";
       print OUT "\n";
       print OUT "pcov_blastx=$O{pcov_blastx} #Blastx Percent Coverage Threhold Protein-Genome Alignments\n";
       print OUT "pid_blastx=$O{pid_blastx} #Blastx Percent Identity Threshold Protein-Genome Aligments\n";
       print OUT "eval_blastx=$O{eval_blastx} #Blastx eval cutoff\n";
       print OUT "bit_blastx=$O{bit_blastx} #Blastx bit cutoff\n";
       print OUT "depth_blastx=$O{depth_blastx} #Blastx depth cutoff (0 to disable cutoff)\n";
       print OUT "\n";
       print OUT "pcov_tblastx=$O{pcov_tblastx} #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments\n";
       print OUT "pid_tblastx=$O{pid_tblastx} #tBlastx Percent Identity Threshold alt-EST-Genome Aligments\n";
       print OUT "eval_tblastx=$O{eval_tblastx} #tBlastx eval cutoff\n";
       print OUT "bit_tblastx=$O{bit_tblastx} #tBlastx bit cutoff\n";
       print OUT "depth_tblastx=$O{depth_tblastx} #tBlastx depth cutoff (0 to disable cutoff)\n";
       print OUT "\n";
       print OUT "pcov_rm_blastx=$O{pcov_rm_blastx} #Blastx Percent Coverage Threhold For Transposable Element Masking\n";
       print OUT "pid_rm_blastx=$O{pid_rm_blastx} #Blastx Percent Identity Threshold For Transposbale Element Masking\n";
       print OUT "eval_rm_blastx=$O{eval_rm_blastx} #Blastx eval cutoff for transposable element masking\n";
       print OUT "bit_rm_blastx=$O{bit_rm_blastx} #Blastx bit cutoff for transposable element masking\n";
       print OUT "\n";
       print OUT "eva_pcov_blastn=$O{eva_pcov_blastn} #EVALUATOR Blastn Percent Coverage Threshold EST-Genome Alignments\n" if($ev);
       print OUT "eva_pid_blastn=$O{eva_pid_blastn} #EVALUATOR Blastn Percent Identity Threshold EST-Genome Alignments\n" if($ev);
       print OUT "eva_eval_blastn=$O{eva_eval_blastn} #EVALUATOR Blastn eval cutoff\n" if($ev);
       print OUT "eva_bit_blastn=$O{eva_bit_blastn} #EVALUATOR Blastn bit cutoff\n" if($ev);
       print OUT "\n" if($ev);
       print OUT "ep_score_limit=$O{ep_score_limit} #Exonerate protein percent of maximal score threshold\n";
       print OUT "en_score_limit=$O{en_score_limit} #Exonerate nucleotide percent of maximal score threshold\n";
       close(OUT);
   }

   #--build maker_exe.ctl file
   if($type eq 'all' || $type eq 'exe'){
       open (OUT, "> $dir/$app\_exe.$ext") or
	   die "ERROR: Could not create $dir/$app\_exe.$ext\n";
       print OUT "#-----Location of Executables Used by MAKER/EVALUATOR\n";
       print OUT "makeblastdb=$O{makeblastdb} #location of NCBI+ makeblastdb executable\n";
       print OUT "blastn=$O{blastn} #location of NCBI+ blastn executable\n";
       print OUT "blastx=$O{blastx} #location of NCBI+ blastx executable\n";
       print OUT "tblastx=$O{tblastx} #location of NCBI+ tblastx executable\n";
       print OUT "formatdb=$O{formatdb} #location of NCBI formatdb executable\n";
       print OUT "blastall=$O{blastall} #location of NCBI blastall executable\n";
       print OUT "xdformat=$O{xdformat} #location of WUBLAST xdformat executable\n";
       print OUT "blasta=$O{blasta} #location of WUBLAST blasta executable\n";
       print OUT "RepeatMasker=$O{RepeatMasker} #location of RepeatMasker executable\n";
       print OUT "exonerate=$O{exonerate} #location of exonerate executable\n";
       print OUT "\n";
       print OUT "#-----Ab-initio Gene Prediction Algorithms\n";
       print OUT "snap=$O{snap} #location of snap executable\n";
       print OUT "gmhmme3=$O{gmhmme3} #location of eukaryotic genemark executable\n";
       print OUT "gmhmmp=$O{gmhmmp} #location of prokaryotic genemark executable\n";
       print OUT "augustus=$O{augustus} #location of augustus executable\n";
       print OUT "fgenesh=$O{fgenesh} #location of fgenesh executable\n";
       print OUT "tRNAscan-SE=$O{'tRNAscan-SE'} #location of trnascan executable\n";
       print OUT "snoscan=$O{snoscan} #location of snoscan executable\n";
       print OUT "\n";
       print OUT "#-----Other Algorithms\n";
       print OUT "fathom=$O{fathom} #location of snap's fathom executable (experimental)\n" if($ev);
       print OUT "probuild=$O{probuild} #location of probuild executable (required for genemark)\n";
       close(OUT);
   }

   #--build server.ctl file
   if($type eq 'server'){
       open (OUT, "> $dir/server.$ext") or
	   die "ERROR: Could not create $dir/server.$ext\n";
       print OUT "#-----Database Setup\n";
       print OUT "DBI=$O{DBI} #interface type to database\n";
       print OUT "dbname=$O{dbname} #database name\n";
       print OUT "host=$O{host} #host on which database is found\n";
       print OUT "port=$O{port} #port on host to access database\n";
       print OUT "username=$O{username} #username to connect to database\n";
       print OUT "password=$O{password} #password to connect to database\n";
       print OUT "\n";
       print OUT "#-----Communication Options\n";
       print OUT "admin_email=$O{admin_email} #Address for sending error and status information\n";
       print OUT "smtp_server=$O{smtp_server} #Outgoing e-mail server\n";
       print OUT "\n";
       print OUT "#-----Web Setup\n";
       print OUT "apache_user=$O{apache_user} #username apache runs as\n";
       print OUT "web_address=$O{web_address} #base web address to server hosting MWAS\n";
       print OUT "cgi_dir=$O{cgi_dir} #web accesible directory to house MWAS CGI content\n";
       print OUT "cgi_web=$O{cgi_web} #url to cgi_dir above (can be relative)\n";
       print OUT "html_dir=$O{html_dir} #web accesible directory to house MWAS HTML conent\n";
       print OUT "html_web=$O{html_web} #url to html_dir (can be relative)\n";
       print OUT "data_dir=$O{data_dir} #directory for saving user uploaded files, running jobs, and storing results\n";
       print OUT "font_file=$O{font_file} #font file for webpage CAPTCHA\n";
       print OUT "\n";
       print OUT "#-----External Viewer Setup\n";
       print OUT "soba_url=$O{soba_url} #url to Sequence Ontology SOBA CGI script\n";
       print OUT "APOLLO_ROOT=$O{APOLLO_ROOT} #base directory for Apollo installation.\n";
       print OUT "JBROWSE_ROOT=$O{JBROWSE_ROOT} #base directory for JBrowse installation.\n";
       print OUT "GBROWSE_MASTER=$O{GBROWSE_MASTER} #path to GBrowse.conf file.\n";
       print OUT "\n";
       print OUT "#-----Environmental Variables for Dependencies\n";
       print OUT "ZOE=$O{ZOE} #required by SNAP, see SNAP documentation.\n";
       print OUT "AUGUSTUS_CONFIG_PATH=$O{AUGUSTUS_CONFIG_PATH} #required by AUGUSTUS, see AUGUSTUS documentation.\n";
       print OUT "\n";
       print OUT "#-----MAKER Server Specific Options\n";
       print OUT "use_login=$O{use_login} #whether to require login to access the web interface, 1 = yes, 0 = no\n";
       print OUT "allow_guest=$O{allow_guest} #enable guest accounts on the server, 1 = yes, 0 = no\n";
       print OUT "allow_register=$O{allow_register} #allow users to register themselves, 1 = yes, 0 = no\n";
       print OUT "tutorials=$O{tutorials} #show example data on \"New Job\" screen, 1 = yes, 0 = no\n";
       print OUT "max_cpus=$O{max_cpus} #maximum number of cpus that can be dedicated to all MAKER jobs\n";
       print OUT "job_cpus=$O{job_cpus} #maximum number of cpus that can be used by a single MAKER job\n";
       print OUT "max_submit_user=$O{max_submit_user} #maximum submission size for registered users (0 = no limit)\n";
       print OUT "max_submit_guest=$O{max_submit_guest} #maximum submission size for guest users (0 = no limit)\n";
       print OUT "persist_user=$O{persist_user} #time results persist for registered users, in hours (0 = no limit)\n";
       print OUT "persist_guest=$O{persist_guest} #time results persist for guest users, in hours (0 = no limit)\n";
       print OUT "inactive_user=$O{inactive_user} #time user account can be inactive, in days (0 = no limit)\n";
       print OUT "inactive_guest=$O{inactive_guest} #time guest account can be inactive, in days (0 = no limit)\n";
       print OUT "MPI=$O{MPI} #use mpi_maker instead of maker\n";
       print OUT "mpiexec=$O{mpiexec} #mpiexec command line for running MPI\n";

       close(OUT);    
   }

   #--build menus.ctl file
   if($type eq 'menus'){
       open (OUT, "> $dir/menus.$ext") or
	   die "ERROR: Could not create $dir/menus.$ext\n";
       print OUT "##Menus from Data::Dumper\n";
       print OUT Data::Dumper->Dump([\%O], [qw(menus)]);
       close(OUT);    
   }
}

1;
