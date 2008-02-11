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
    my $vars       = shift; #this should be an array reference
    $self->{ID}    = shift || 0;
    $self->{RANK} = 0;
    
    local $Storable::forgive_me = 1; #I hate code references
    $self->{VARS} = freeze($vars);
}

#--------------------------------------------------------------
sub run {
    my $self = shift;
    $self->{RANK} = shift;

    if (exists $self->{RESULT}){
        return;
    }

    my $level = $self->{LEVEL};
    my @vars = @{thaw($self->{VARS})};
    undef $self->{VARS};
    my @results;

    #--redirect STDERR to a log file
    open (OLDERR, ">&STDERR");
    my (undef, $t_name) = tempfile();
    close(STDERR);
    open(STDERR, "> $t_name");

    if ($level == 0) {
	#------------------------ARGS_IN
	my $fasta       = shift @vars;
	my $the_void    = shift @vars;
	my $seq_id      = shift @vars;
	my $query_seq   = shift @vars;
	my $query_def   = shift @vars;
	my %CTL_OPTIONS = %{shift @vars};
	my $opt_f       = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- repeatmask the input file
	my ($temp_masked_fasta, $rma_keepers) = repeatmask($fasta, 
							   $the_void,
							   length($$query_seq),
							   $CTL_OPTIONS{'max_dna_len'},
							   $seq_id,
							   $query_seq,
							   $query_def,
							   $CTL_OPTIONS{'model_org'},
							   $CTL_OPTIONS{'RepeatMasker'},
							   $CTL_OPTIONS{'rmlib'},
							   $opt_f
							  );
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($temp_masked_fasta, $rma_keepers);
	#------------------------RESULTS
    }
    elsif ($level == 1) {
	#------------------------ARGS_IN
	my $temp_masked_fasta = shift @vars;
	my $repeat_protein    = shift @vars;
	my $the_void          = shift @vars;
	my $query_seq         = shift @vars;
	my $seq_id            = shift @vars;
	my %CTL_OPTIONS       = %{shift @vars};
	my $opt_f             = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- blastx agains a repeat library (for better masking)
	my $repeat_blastx_keepers = blastx($temp_masked_fasta,
					   $repeat_protein,
					   $the_void,
					   length($$query_seq),
					   $CTL_OPTIONS{'max_dna_len'},
					   $seq_id,
					   $CTL_OPTIONS{blastx},
					   $CTL_OPTIONS{eval_blastx},
					   $CTL_OPTIONS{bit_blastx},
					   $CTL_OPTIONS{percov_blastx},
					   $CTL_OPTIONS{percid_blastx},
					   $CTL_OPTIONS{cpus},
					   $CTL_OPTIONS{old_repeat_protein},
					   $CTL_OPTIONS{split_hit},
					   $CTL_OPTIONS{setdb},
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
	my $query_seq             = shift @vars;
	my $rma_keepers           = shift @vars;
	my $repeat_blastx_keepers = shift @vars;
	my $query_def             = shift @vars;
	my $the_void              = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	my($f_util) = File::Util->new();
	my(@dirs) = $f_util->list_dir($the_void, '--dirs-only');
	@dirs = grep (/\.blastx\.temp_dir$/, @dirs);

	foreach my $dir (@dirs){
	    my $blast_finished = $dir;
	    $blast_finished =~ s/\.temp_dir$//;
	    system ("cat $the_void/$dir/*.blastx > $the_void/$blast_finished");
	    rmtree ("$the_void/$dir");
	}

	#-- now use the repeatmasker data + the blastx repeat data
	#-- to generate a better masked fasta
	my ($masked_seq, $rm_keepers) = repeat_mask_seq::process($query_seq, 
								 $rma_keepers, 
								 $repeat_blastx_keepers
								);
	
	my $masked_fasta = Fasta::toFasta($query_def.' masked', \$masked_seq);
	#-------------------------CHUNK
	    
	#------------------------RESULTS
	@results = ($rm_keepers, $masked_fasta);
	#------------------------RESULTS
    }
    elsif ($level == 3) {
	#------------------------ARGS_IN
	my $masked_fasta = shift @vars;
	my $the_void     = shift @vars;
	my $out_dir      = shift @vars;
	my $seq_out_name = shift @vars;
	my $rm_keepers   = shift @vars;
	my %CTL_OPTIONS  = %{shift @vars};
	#------------------------ARGS_IN

	#-------------------------CHUNK
	FastaFile::writeFile($masked_fasta ,$the_void."/query.masked.fasta");
    
	#$IOG= new Dumper::GFF::GFF3->init(FastaFile=>$CTL_OPTIONS{'genome'},
	#				  WriteFile=>$out_dir."/".$seq_out_name.".gff");
    
	my $IOX= new Dumper::XML::Game->init(FastaFile=>$CTL_OPTIONS{'genome'},
					     WriteFile=>$out_dir."/".$seq_out_name.".xml");
	my $i=0;
	my $j=0;
    
	#--add repeat masked hits to be included
	$IOX->add_hits($rm_keepers);

    	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($IOX);
	#------------------------RESULTS
    }
    elsif ($level == 4) {
	#------------------------ARGS_IN
	my $masked_fasta = shift @vars;
	my $proteins     = shift @vars;
	my $the_void     = shift @vars;
	my $query_seq    = shift @vars;
	my $seq_id       = shift @vars;
	my %CTL_OPTIONS  = %{shift @vars};
	my $opt_f        = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- blastx search  the masked input file
	my $blastx_keepers = blastx($$masked_fasta, 
				    $proteins,
				    $the_void,
				    length($$query_seq),
				    $CTL_OPTIONS{'max_dna_len'},
				    $seq_id,
				    $CTL_OPTIONS{blastx},
				    $CTL_OPTIONS{eval_blastx},
				    $CTL_OPTIONS{bit_blastx},
				    $CTL_OPTIONS{percov_blastx},
				    $CTL_OPTIONS{percid_blastx},
				    $CTL_OPTIONS{cpus},
				    $CTL_OPTIONS{old_protein},
				    $CTL_OPTIONS{split_hit},
				    $CTL_OPTIONS{setdb},
				    $self->id(),
				    $self->{RANK},
				    $opt_f
				   );
	#-------------------------CHUNK
	
	#------------------------RESULTS
	@results = ($blastx_keepers);
	#------------------------RESULTS
    }
    elsif ($level == 5) {
	#------------------------ARGS_IN
	my $masked_fasta = shift @vars;
	my $cests        = shift @vars;
	my $the_void     = shift @vars;
	my $query_seq    = shift @vars;
	my $seq_id       = shift @vars;
	my %CTL_OPTIONS  = %{shift @vars};
	my $opt_f        = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	my $tblastx_keepers = [];

	if ($cests) {
	    $tblastx_keepers = tblastx($$masked_fasta, 
				       $cests,
				       $the_void,
				       length($$query_seq),
				       $CTL_OPTIONS{'max_dna_len'},
				       $seq_id,
				       $CTL_OPTIONS{tblastx},
				       $CTL_OPTIONS{eval_tblastx},
				       $CTL_OPTIONS{bit_tblastx},
				       $CTL_OPTIONS{percov_tblastx},
				       $CTL_OPTIONS{percid_tblastx},
				       $CTL_OPTIONS{cpus},
				       $CTL_OPTIONS{old_alt_est},
				       $CTL_OPTIONS{pressdb},
				       $self->id(),
				       $self->{RANK},
				       $opt_f
				      );
	}
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($tblastx_keepers);
	#------------------------RESULTS
    }
    elsif ($level == 6) {
	#------------------------ARGS_IN
	my $masked_fasta = shift @vars;
	my $transcripts  = shift @vars;
	my $the_void     = shift @vars;
	my $query_seq    = shift @vars;
	my $seq_id       = shift @vars;
	my %CTL_OPTIONS  = %{shift @vars};
	my $opt_f        = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- blastn search the file against ESTs
	my $blastn_keepers = blastn($$masked_fasta,
				    $transcripts,
				    $the_void,
				    length($$query_seq),
				    $CTL_OPTIONS{'max_dna_len'},
				    $seq_id,
				    $CTL_OPTIONS{blastn},
				    $CTL_OPTIONS{eval_blastn},
				    $CTL_OPTIONS{bit_blastn},
				    $CTL_OPTIONS{percov_blastn},
				    $CTL_OPTIONS{percid_blastn},
				    $CTL_OPTIONS{cpus},
				    $CTL_OPTIONS{old_est},
				    $CTL_OPTIONS{pressdb},
				    $self->id(),
				    $self->{RANK},
				    $opt_f
				   );
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($blastn_keepers);
	#------------------------RESULTS
    }
    elsif ($level == 7) {
	#------------------------ARGS_IN
	my $blastx_keepers  = shift @vars;
	my $tblastx_keepers = shift @vars;
	my $blastn_keepers  = shift @vars;
	my $IOX             = shift @vars;
	my $query_seq       = shift @vars;
	my $the_void        = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	my($f_util) = File::Util->new();
	my(@dirs) = $f_util->list_dir($the_void, '--dirs-only');
	@dirs = grep (/blast.\.temp_dir$/, @dirs);

	foreach my $dir (@dirs){
	    my $blast_finished = $dir;
	    $blast_finished =~ s/\.temp_dir$//;
	    system ("cat $the_void/$dir/*blast* > $the_void/$blast_finished");
	    rmtree ("$the_void/$dir");
	}

	#--add blastx and blastn hits 
	$IOX->add_hits($blastx_keepers);
	$IOX->add_hits($tblastx_keepers);
	$IOX->add_hits($blastn_keepers);
    
	#-- cluster the blastx hits
	print STDERR "cleaning blastx...\n";
    
	my $blastx_clusters = cluster::clean_and_cluster($blastx_keepers,
							 $query_seq,
							 10);
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($IOX, $blastx_clusters);
	#------------------------RESULTS
    }
    elsif ($level == 8) {
	#------------------------ARGS_IN
	my $fasta           = shift @vars;
	my $blastx_clusters = shift @vars;
	my $proteins        = shift @vars;
	my $the_void        = shift @vars;
	my %CTL_OPTIONS     = %{shift @vars};
	my $opt_f           = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- build an index of the databases
	my $fasta_p_index     = build_fasta_index($proteins);

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
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($exonerate_p_clusters);
	#------------------------RESULTS
    }
    elsif ($level == 9) {
	#------------------------ARGS_IN
	my $exonerate_p_clusters = shift @vars;
	my $IOX                  = shift @vars;
	my $tblastx_keepers      = shift @vars;
	my $query_seq            = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#--add hits
	foreach my $pp (@$exonerate_p_clusters) {
	    $IOX->add_hits($pp);
	}
	#-- cluster the tblastx hits
	print STDERR "cleaning blastx...\n";
    
	my $tblastx_clusters = cluster::clean_and_cluster($tblastx_keepers,
							  $query_seq,
							  10);
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($IOX, $tblastx_clusters);
	#------------------------RESULTS
    }
    elsif ($level == 10) {
	#------------------------ARGS_IN
	my $blastn_keepers = shift @vars;
	my $query_seq      = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- make a multi-fasta of the seqs in the tblastx_clusters 
	#-- polish the tblastx hits with exonerate
    
	#####There is no exonerate with Tblastx#####
    
	#-- Cluster the blastn hits
	print STDERR "cleaning blastn...\n";
	my $blastn_clusters = cluster::clean_and_cluster($blastn_keepers,
							 $query_seq,
							 10);
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($blastn_clusters);
	#------------------------RESULTS
    }
    elsif ($level == 11) {
	#------------------------ARGS_IN
	my $fasta           = shift @vars;
	my $blastn_clusters = shift @vars;
	my $transcripts     = shift @vars;
	my $the_void        = shift @vars;
	my %CTL_OPTIONS     = %{shift @vars};
	my $opt_f           = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- build an index of the databases
	my $fasta_t_index     = build_fasta_index($transcripts);

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
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($exonerate_e_clusters);
	#------------------------RESULTS
    }
    elsif ($level == 12) {
	#------------------------ARGS_IN
	my $exonerate_e_clusters = shift @vars;
	my $IOX                  = shift @vars;
	my $masked_fasta         = shift @vars;
	my $the_void             = shift @vars;
	my $seq_id               = shift @vars;
	my %CTL_OPTIONS          = %{shift @vars};
	my $opt_f                = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#--add hits
	foreach my $ee (@$exonerate_e_clusters) {
	    $IOX->add_hits($ee);
	}
	#####working here###########
	#-- gene_find the input file
	my $snaps = snap($$masked_fasta, $the_void, length($$masked_fasta), $CTL_OPTIONS{'max_dna_len'},
			 $seq_id,$CTL_OPTIONS{snap}, $CTL_OPTIONS{'snaphmm'}, $opt_f
			);

	$IOX->add_hits($snaps);
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($IOX, $snaps);
	#------------------------RESULTS
    }
    elsif ($level == 13) {
	#------------------------ARGS_IN
	my $fasta                = shift @vars;
	my $masked_fasta         = shift @vars;
	my $exonerate_p_clusters = shift @vars;
	my $exonerate_e_clusters = shift @vars;
	my $blastx_clusters      = shift @vars;
	my $snaps                = shift @vars;
	my $the_void             = shift @vars;
	my %CTL_OPTIONS          = %{shift @vars};
	my $opt_f                = shift @vars;
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#-- auto-annotate the input file
	my $snap_command = $CTL_OPTIONS{snap}.' '.$CTL_OPTIONS{snaphmm};
	my $snap_flank   = $CTL_OPTIONS{snap_flank};

	my $annotations = auto_annotate($fasta,
					$$masked_fasta,
					$exonerate_p_clusters,
					$exonerate_e_clusters,
					$blastx_clusters,
					$snaps,
					$the_void,
					$snap_command,
					$snap_flank,
					$CTL_OPTIONS{'single_exon'},
					$opt_f
				       );
	#-------------------------CHUNK

	#------------------------RESULTS
	@results = ($annotations);
	#------------------------RESULTS
    }
    elsif ($level == 14) {
	#------------------------ARGS_IN
	my $blastx_clusters      = shift @vars;
	my $blastn_clusters      = shift @vars;
	my $tblastx_clusters     = shift @vars;
	my $exonerate_p_clusters = shift @vars;
	my $exonerate_e_clusters = shift @vars;
	my $query_seq            = shift @vars;
	my $seq_id               = shift @vars;
	my $annotations          = shift @vars;
	my $rm_keepers           = shift @vars;
	my $snaps                = shift @vars;
	my $out_dir              = shift @vars;
	my $seq_out_name         = shift @vars;
	my $IOX                  = shift @vars;
	my $the_void             = shift @vars;
	my %CTL_OPTIONS          = %{shift @vars};
	#------------------------ARGS_IN

	#-------------------------CHUNK
	#--- new GFF3
	my $blastx_data      = flatten($blastx_clusters);
	my $blastn_data      = flatten($blastn_clusters);
	my $tblastx_data     = flatten($tblastx_clusters);
	my $exonerate_p_data = flatten($exonerate_p_clusters, 'exonerate:p');
	my $exonerate_e_data = flatten($exonerate_e_clusters, 'exonerate:e');

	my $GFF3 = new Dumper::GFF::GFFV3();

	$GFF3->seq($query_seq);
	$GFF3->seq_id($seq_id);


	$GFF3->genes($annotations);
	$GFF3->phat_hits($rm_keepers);
	$GFF3->phat_hits($blastx_data);
	$GFF3->phat_hits($blastn_data);
	$GFF3->phat_hits($tblastx_data);
	$GFF3->phat_hits($exonerate_p_data);
	$GFF3->phat_hits($exonerate_e_data);
	$GFF3->predictions($snaps);
	$GFF3->print($out_dir."/".$seq_out_name.".gff");
	#$GFF3->print('foo.gff');
	
	#--adding auto annotations
	my $p_fastas = '';
	my $t_fastas = '';
	my @quality_indices;
	foreach my $an (@$annotations) {
	    my $g_name     = $an->{g_name};
	    my $g_s        = $an->{g_start};
	    my $g_e        = $an->{g_end};
	    my $g_strand   = $an->{g_strand};
	
	    #$IOG->add_auto_gene($an);
	    my @temp_ant;
	    foreach my $a (@{$an->{t_structs}}) {
		push(@quality_indices, [$a->{t_name}, $a->{t_qi}]);
		
		my ($p_fasta, $t_fasta) = get_p_and_t_fastas($a);
		
		$t_fastas .= $$t_fasta;
		$p_fastas .= $$p_fasta;
		
		#$IOG->add_autoant($a) if defined($a->{hit});
		push (@temp_ant,$a) if defined($a->{hit});
		$IOX->add_autoant($a) if defined($a->{hit});
	    }
	    #$IOG->add_autoant(\@temp_ant); 
	}
	
	#Write the quality index of the mRNAs to a separate file.
	#write_quality_data(\@quality_indices, $seq_id);

	FastaFile::writeFile(\$p_fastas ,"$out_dir\/$seq_out_name\.maker.proteins.fasta");
	FastaFile::writeFile(\$t_fastas ,"$out_dir\/$seq_out_name\.maker.transcripts.fasta");

	#$IOG->GFF3();
	$IOX->Game();

	my ($p_snap_fastas, $t_snap_fastas) = get_snap_p_and_t_fastas($query_seq, $snaps);

	FastaFile::writeFile(\$p_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.proteins.fasta");
	FastaFile::writeFile(\$t_snap_fastas ,"$out_dir\/$seq_out_name\.maker.snap.transcript.fasta");

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

    local $Storable::forgive_me = 1; #I hate code references
    $self->{RESULT} = freeze(\@results);
}

#--------------------------------------------------------------
sub result {
    my $self = shift;
    
    return @{thaw($self->{RESULT})};
}

#--------------------------------------------------------------
sub error {
    my $self = shift;
    
    return $self->{ERROR};
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
sub auto_annotate {
    my $virgin_fasta         = shift;
    my $masked_fasta         = shift;
    my $exonerate_p_clusters = shift;
    my $exonerate_e_clusters = shift;
    my $blastx_clusters      = shift;
    my $snaps                = shift;
    my $the_void             = shift;
    my $snap_command         = shift;
    my $snap_flank           = shift;
    my $single_exon          = shift;
    my $opt_f                = shift;

    my $blastx_hits      = flatten($blastx_clusters);
    my $exonerate_p_hits = flatten($exonerate_p_clusters, 'exonerate:p');
    my $exonerate_e_hits = flatten($exonerate_e_clusters, 'exonerate:e');

	
    my $annotations = maker::auto_annotator::annotate($virgin_fasta,
						      $masked_fasta,
						      $exonerate_p_hits,
						      $exonerate_e_hits,
						      $blastx_hits,
						      $snaps,
						      $the_void,
						      $snap_command,
						      $snap_flank,
						      $single_exon,
						      $opt_f
						     );
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
    my $q_length   = shift;
    my $chunk_size = shift;
    my $seq_id     = shift;
    my $snap = shift;
    my $snaphmm = shift;
    my $opt_f = shift;
	
    my $fasta_chunker = new FastaChunker();
    $fasta_chunker->parent_fasta($fasta);
    $fasta_chunker->chunk_size($chunk_size);
    $fasta_chunker->load_chunks();
	
    my %params;
    my $snap_keepers = [];
    my $i = 0;
    while (my $chunk = $fasta_chunker->get_chunk($i)) {
	my $chunk_number = $chunk->number();
	my $file_name = "$the_void/$seq_id\.$chunk_number";
	my $o_file    = "$the_void/$seq_id\.$chunk_number\.snap";

	$chunk->write_file($file_name);
		
	runSnap($file_name, $o_file, $snap, $snaphmm, $opt_f);
		
	$params{min_exon_score}  = -100000;	    #-10000;
	$params{min_gene_score}  = -100000;	    #0;
		
	my $chunk_keepers =
	Widget::snap::parse($o_file,
			    \%params,
			    $file_name,
			   );
	PhatHit_utils::add_offset($chunk_keepers,
				  $chunk->offset(),
				 );
	PhatHit_utils::merge_hits($snap_keepers,
				  $chunk_keepers,
				  10000,
				 );
	$chunk->erase_fasta_file();
	$i++;
    }
    return $snap_keepers;
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

	    if ($type eq 'e') {
		next if $hit->pAh < $percov;
		next if $hit->hsp('best')->frac_identical < $percid;
	    }
	    elsif ($type eq 'p') {
		next if $hit->pAh < $percov;
		next if $hit->hsp('best')->frac_identical < $percid;
	    }
	    my ($nB, $nE) =
	    PhatHit_utils::get_span_of_hit($hit,'query');
	    my @coors = [$nB, $nE];
	    my $p = Shadower::getPieces($seq, \@coors, 50);
	    my $p_def = $def." ".$p->[0]->{b}." ".$p->[0]->{e};
	    my $p_fasta = Fasta::toFasta($p_def, \$p->[0]->{piece});
	    my ($name) = $p_def =~ />([^\s\t\n]+)/;
	    my $safe_name = uri_escape($name,  #build a safe name for file names from the sequence identifier
				       '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				       );
	    $safe_name .= '.fasta';
	    my $d_file = $the_void."/".$safe_name.'.'.$i.'.'.$n;
	    FastaFile::writeFile($p_fasta, $d_file);
	    my $offset = $p->[0]->{b};
	    my $id  = $hit->name();
	    $id =~ s/\s+/_/g;
	    $id =~ s/\|/_/g;
	    my $fastaObj = $db_index->get_Seq_by_id($hit->name);
	    if (not $fastaObj) {
		print "stop here:".$hit->name."\n";
		die;
	    }
	    my $seq2      = $fastaObj->seq();
	    my $def2      = $db_index->header($hit->name);
	    $def2 =~ s/\|/_/g;
	    my $fasta    = Fasta::toFasta('>'.$def2, \$seq2);
	    my $safe_id = uri_escape($id,  #build a safe name for file names from the sequence identifier
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
    my $fasta      = shift;
    my $the_void    = shift;
    my $q_length   = shift;
    my $chunk_size = shift;
    my $seq_id     = shift;
    my $query_seq = shift;
    my $query_def = shift;
    my $model_org = shift;
    my $RepeatMasker = shift;
    my $rmlib = shift;
    my $opt_f = shift;

    my $fasta_chunker = new FastaChunker();
    $fasta_chunker->parent_fasta($fasta);
    $fasta_chunker->chunk_size($chunk_size);
    $fasta_chunker->load_chunks();
	
    my $rm_keepers = [];
    my $i = 0;
    while (my $chunk = $fasta_chunker->get_chunk($i)) {

	my $chunk_number = $chunk->number();
	my $file_name = "$the_void/$seq_id\.$chunk_number";
	my $o_file    = "$the_void/$seq_id\.$chunk_number\.out";
	$chunk->write_file($file_name);
		
	runRepeatMasker($file_name, 
			$model_org, 
			$the_void, 
			$o_file,
			$RepeatMasker,
			$rmlib,
			$opt_f,
		       );	# -no_low
		
	my $rm_chunk_keepers = 
	Widget::RepeatMasker::parse($o_file, 
				    $seq_id, 
				    $q_length,
				   );
	PhatHit_utils::add_offset($rm_chunk_keepers, 
				  $chunk->offset(),
				 );
	PhatHit_utils::merge_hits($rm_keepers,  
				  $rm_chunk_keepers, 
				  20,
				 );
	$chunk->erase_fasta_file();
	$i++;
    }
	
    my ($tes, $lcs) = repeat_mask_seq::seperate_types($rm_keepers);
	
    my $masked_seq = repeat_mask_seq::mask_seq($query_seq, $tes, $lcs);
	
    my $masked_fasta = Fasta::toFasta($query_def.' masked', \$masked_seq);
	
    return ($$masked_fasta, $rm_keepers);
}

#-----------------------------------------------------------------------------
sub blastn {
    my $fasta      = shift;
    my $db         = shift;
    my $the_void    = shift;
    my $q_length   = shift;
    my $chunk_size = shift;
    my $seq_id     = shift;
    my $blastn = shift;
    my $eval_blastn = shift;
    my $bit_blastn = shift,
    my $percov_blastn = shift;
    my $percid_blastn = shift;
    my $cpus = shift;
    my $old_db = shift;
    my $pressdb = shift;
    my $id = shift;
    my $rank = shift;
    my $opt_f = shift;

    my ($db_n) = $db =~ /([^\/]+)$/;
    $db_n  =~ s/\.fasta$//;
	
    my $fasta_chunker = new FastaChunker();
    $fasta_chunker->parent_fasta($fasta);
    $fasta_chunker->chunk_size($chunk_size);
    $fasta_chunker->load_chunks();
	
    my $blastn_keepers = [];
    my $i = 0;
    while (my $chunk = $fasta_chunker->get_chunk($i)) {
	my $chunk_number = $chunk->number();

	my ($db_old_n) = $old_db =~ /([^\/]+)$/;
	$db_old_n  =~ s/\.fasta$//;
	my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.blastn";

	my $t_dir = "/tmp/rank".$rank;
	mkpath($t_dir);

	my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
	my $o_file    = "$blast_finished\.temp_dir/$db_n\.blastn";

	$db =~ /([^\/]+$)/;
	my $tmp_db = "$t_dir/$1";

	if (-e $blast_finished && ! $opt_f){
	    $o_file = $blast_finished;
	    
	    return [] if ($id !~ /\:0$/);
	}
	elsif (! -e $blast_finished && ! `ls $tmp_db\.n*`) {
	    system("cp $db $tmp_db");
	    system ("$pressdb $tmp_db");
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
	my $chunk_keepers =
	Widget::blastn::parse($o_file,
			      \%params,
			     );
	PhatHit_utils::add_offset($chunk_keepers,
				  $chunk->offset(),
				 );
	PhatHit_utils::merge_hits($blastn_keepers,
				  $chunk_keepers,
				  10000,
				 );
	$chunk->erase_fasta_file();
	$i++;
    }
    
    my @purge;
    
    foreach my $hit (@{$blastn_keepers}) {
	next unless $hit->pAh > $percov_blastn;
	next unless $hit->hsp('best')->frac_identical() > $percid_blastn;
	next unless PhatHit_utils::is_contigous($hit);
	push(@purge, $hit);
    }
    
    my $a = @{$blastn_keepers};
    my $b = @purge;
    my $diff = $a - $b;
    print STDERR "purging blastns deleted $diff hits!\n";
    sleep 1;
    return \@purge;
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
    my $fasta      = shift;
    my $db         = shift;
    my $the_void    = shift;
    my $q_length   = shift;
    my $chunk_size = shift;
    my $seq_id     = shift;
    my $blastx = shift;
    my $eval_blastx = shift;
    my $bit_blastx = shift;
    my $percov_blastx = shift;
    my $percid_blastx = shift;
    my $cpus = shift;
    my $old_db = shift;
    my $split_hit = shift;
    my $setdb = shift;
    my $id = shift;
    my $rank = shift;
    my $opt_f = shift;
	
    my ($db_n) = $db =~ /([^\/]+)$/;
    $db_n  =~ s/\.fasta$//;
	
    my $fasta_chunker = new FastaChunker();
    $fasta_chunker->parent_fasta($fasta);
    $fasta_chunker->chunk_size($chunk_size);
    $fasta_chunker->load_chunks();
	
    my $blastx_keepers = [];
    my $i = 0;
    while (my $chunk = $fasta_chunker->get_chunk($i)) {
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

	if (-e $blast_finished && ! $opt_f){
	    $o_file = $blast_finished ;
	    
	    return [] if ($id !~ /\:0$/);
	}
	elsif (! `ls $tmp_db\.n*`) {
	    system("cp $db $tmp_db");
	    system ("$setdb $tmp_db");
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
	$params{significance} = $eval_blastx;
	$params{hsp_bit_min}  = $bit_blastx;
	
	my $chunk_keepers =
	Widget::blastx::parse($o_file,
			      \%params,
			     );
	PhatHit_utils::add_offset($chunk_keepers,
				  $chunk->offset(),
				 );
	PhatHit_utils::merge_hits($blastx_keepers,
				  $chunk_keepers,
				  10000,
				 );
	$chunk->erase_fasta_file();
	$i++;
    }
    my @purge;
    foreach my $hit (@{$blastx_keepers}) {
	#my $new_hit = PhatHit_utils::normalize($hit, 'hit');
	my $split_hits = PhatHit_utils::split_hit($hit, $split_hit);
							  
	foreach my $s_hit (@{$split_hits}) {
	    push(@purge, $s_hit);
	}
							  
							  
	#if (!PhatHit_utils::is_contigous($hit)){
	#	my $shatter_hits = PhatHit_utils::shatter_hit($hit);
	#	push(@purge, @{$shatter_hits});
	#}
	#else {
	#	push(@purge, $hit);
	#}
    }
    my $a = @{$blastx_keepers};
    my $b = @purge;
    my $diff = $a - $b;
    print STDERR "purging blastxs deleted $diff hits!\n";
    sleep 1;
    return \@purge;
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
sub tblastx {
    my $fasta      = shift;
    my $db         = shift;
    my $the_void   = shift;
    my $q_length   = shift;
    my $chunk_size = shift;
    my $seq_id     = shift;
    my $tblastx = shift;
    my $eval_tblastx = shift;
    my $bit_tblastx = shift;
    my $percov_tblastx = shift;
    my $percid_tblastx = shift;
    my $cpus = shift;
    my $old_db = shift;
    my $setdb = shift;
    my $id = shift;
    my $rank = shift;
    my $opt_f = shift;
	
    my ($db_n) = $db =~ /([^\/]+)$/;
    $db_n  =~ s/\.fasta$//;
	
    my $fasta_chunker = new FastaChunker();
    $fasta_chunker->parent_fasta($fasta);
    $fasta_chunker->chunk_size($chunk_size);
    $fasta_chunker->load_chunks();
	
    my $tblastx_keepers = [];
    my $i = 0;
    while (my $chunk = $fasta_chunker->get_chunk($i)) {
	my $chunk_number = $chunk->number();

	my $t_dir = "/tmp/rank".$rank;
        mkpath($t_dir);
	
	my ($db_old_n) = $old_db =~ /([^\/]+)$/;
	$db_old_n  =~ s/\.fasta$//;
	my $blast_finished = "$the_void/$seq_id\.$chunk_number\.$db_old_n\.tblastx";
		
	my $t_file_name = "$t_dir/$seq_id\.$chunk_number";
	my $o_file    = "$blast_finished\.temp_dir/$db_n\.tblastx";

	$db =~ /([^\/]+$)/;
	my $tmp_db = "$t_dir/$1";

	if (-e $blast_finished && ! $opt_f){
	    $o_file = $blast_finished ;
	    
	    return [] if ($id !~ /\:0$/);
	}	
	elsif (! `ls $tmp_db\.n*`) {
	    system("cp $db $tmp_db");
	    system ("$setdb $tmp_db");
	}

	$chunk->write_file($t_file_name);	

	runTblastx($t_file_name,
		   $tmp_db,
		   $o_file,
		   $tblastx,
		   $eval_tblastx,
		   $cpus,
		   $opt_f
		  );
	my %params;
	$params{significance} = $eval_tblastx;
	$params{hsp_bit_min}  = $bit_tblastx;
	
	my $chunk_keepers =
	Widget::tblastx::parse($o_file,
			       \%params,
			      );
	PhatHit_utils::add_offset($chunk_keepers,
				  $chunk->offset(),
				 );
	PhatHit_utils::merge_hits($tblastx_keepers,
				  $chunk_keepers,
				  10000,
				 );
	$chunk->erase_fasta_file();
	$i++;
    }
    my @purge;
    foreach my $hit (@{$tblastx_keepers}) {
	if (!PhatHit_utils::is_contigous($hit)) {
	    my $shatter_hits = PhatHit_utils::shatter_hit($hit);
	    push(@purge, @{$shatter_hits});
	}
	else {
	    push(@purge, $hit);
	}
    }
    my $a = @{$tblastx_keepers};
    my $b = @purge;
    my $diff = $a - $b;
    print STDERR "purging tblastxs deleted $diff hits!\n";
    sleep 1;
    return \@purge;
}

#-----------------------------------------------------------------------------
sub runTblastx {
    my $q_file   = shift;
    my $db       = shift;
    my $out_file = shift;
    my $tblastx = shift;
    my $eval_tblastx = shift;
    my $cpus = shift;
    my $opt_f = shift;

    my $command  = $tblastx;
    $command .= " $db $q_file B=10000 V=10000 E=$eval_tblastx";
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
    $command .= " > $out_file";
    my $w = new Widget::tblastx();
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
sub runRepeatMasker {
    my $q_file   = shift;
    my $species  = shift;
    my $dir      = shift;
    my $o_file   = shift;
    my $RepeatMasker = shift;
    my $rmlib = shift;
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
    my $run_id  = shift;
    my $out_dir = shift;

    my $vid = "theVoid\.$seq_id\.$run_id";   
    my $the_void = "$out_dir"."$vid";
    mkpath ($the_void);

    return $the_void;
}

#-----------------------------------------------------------------------------
1;
