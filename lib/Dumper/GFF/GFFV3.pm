#------------------------------------------------------------------------
#----                         Dumper::GFF::GFFV3                     ---- 
#------------------------------------------------------------------------
package Dumper::GFF::GFFV3;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use PhatHit_utils;
use File::Copy;
use URI::Escape;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- METHODS   ----------------------------------
#------------------------------------------------------------------------
sub new {
   my $class = shift;
   my $filename = shift;
   my $build = shift || "Build1";;
   my $t_dir = shift || "/tmp";

   my $self = {};

   bless ($self, $class);

   $self->_initialize($filename, $build, $t_dir);

   return $self;
}
#------------------------------------------------------------------------
sub _initialize {
   my $self = shift;
   my $filename = shift;
   my $build = shift;
   my $t_dir = shift;

   $self->{build} = $build;

   my $gff_file = "$filename";
   my ($name) = $filename =~ /([^\/]+)$/;
   my $ann_file = "$t_dir/$name.ann";
   my $seq_file = "$t_dir/$name.seq";

   open(my $ANN, "> $ann_file") || die "ERROR: Could not open file: $ann_file\n";
   print_txt($ANN, $self->header."\n");
   close($ANN);

   open(my $SEQ, "> $seq_file") || die "ERROR: Could not open file: $seq_file\n";
   print_txt($SEQ, "##FASTA\n");
   close($SEQ);

   $self->{gff_file} = $gff_file;
   $self->{ann_file} = $ann_file;
   $self->{seq_file} = $seq_file;
}
#------------------------------------------------------------------------
sub fasta {
    my $self  = shift;

    my $fasta_ref = Fasta::toFastaRef('>'.$self->{seq_id}, \uc(${$self->seq}));

    return $fasta_ref;
}
#------------------------------------------------------------------------
sub resolved_flag {
    my $self   = shift;

    open(my $ANN, ">>", $self->{ann_file}) || die "ERROR: Can't open annotation file\n\n";
    print_txt($ANN, "###\n");
    close($ANN);
}
#------------------------------------------------------------------------
sub set_current_contig {
    my $self   = shift;

    my $flag = 0;
    $flag = 1 if (defined $self->{seq_id});
    $self->{seq_id} = shift;
    $self->{seq} = shift;

    $self->{seq_id} = uri_escape($self->{seq_id}, "^a-zA-Z0-9\.\:\^\*\\\$\@\!\+\_\?\\-\|"); #per gff standards

    open(my $ANN, ">>", $self->{ann_file}) || die "ERROR: Can't open annotation file\n\n";
    print_txt($ANN, "###\n") if($flag);

    #skip adding the optional sequence-region line because many programs
    #do not handle it correctly or consistently - 06/05/2010
    #print_txt($ANN, $self->contig_comment."\n");

    print_txt($ANN, $self->contig_line."\n"); 
    close($ANN);

    open(my $SEQ, ">>", $self->{seq_file}) || die "ERROR: opening fasta for GFF3\n\n";
    print_txt($SEQ, ${$self->fasta});
    close($SEQ);
}
#------------------------------------------------------------------------
sub seq {
    my $self   = shift;

    return $self->{seq} || undef;
}

#------------------------------------------------------------------------
sub seq_id {
    my $self   = shift;

    return $self->{seq_id} || undef;
}
#------------------------------------------------------------------------
sub finalize {
    my $self = shift;

    my $gff_f = $self->{gff_file};
    my $ann_f = $self->{ann_file};
    my $seq_f = $self->{seq_file};

    open(my $SEQ, "< $seq_f")|| die "ERROR: Can't open seq file\n\n";
    open(my $ANN, ">> $ann_f")|| die "ERROR: Can't open annotation file\n\n";
    my $line = <$SEQ>;
    if($line !~ /^\#\#FASTA/){
	die "ERROR: There was a problem in the writing the fasta entry\n".
	    "Either no sequence was given, or there was an error in writing\n\n";
    }
    print $ANN $line;
    $line = <$SEQ>;
    if($line !~ /^>/){
	die "ERROR: There was a problem in the writing the fasta entry\n".
	    "Either no sequence was given, or there was an error in writing\n\n";
    }
    print $ANN $line;
    while(defined($line = <$SEQ>)){
	print $ANN $line;
    }
    close($ANN);
    close($SEQ);

    unlink($seq_f);
    move($ann_f, $gff_f);

    die "ERROR: GFF3 file not created\n" if(! -e $gff_f);
}
#------------------------------------------------------------------------
sub contig_line {
    my $self = shift;
    
    die "no contig seq in Dumper::GFFV3::contig_line\n"
	unless defined($self->seq);
    
    my $seq = $self->seq();
    
    my $length = length($$seq);
    my $id     = $self->seq_id();
    my $name   = $id;
    if ($id =~ m/gnl\%7Cv3\%7C(Contig\d+)/) {
	$name = $1;
    }
    my @data;
    push(@data, $id, '.', 'contig', 1, $length, '.','.','.');
    push(@data, 'ID='.$id.';Name='.$name.';');
    
    my $l = join("\t", @data);
    return $l;
}
#------------------------------------------------------------------------
sub contig_comment {
    my $self = shift;

    die "no contig seq in Dumper::GFFV3::contig_comment\n"
        unless defined($self->seq);

    my $seq = $self->seq;

    my $length = length($$seq);
    my $id     = $self->seq_id;
    my $name   = $id;
    if ($id =~ m/gnl\%7Cv3\%7C(Contig\d+)/) {
        $name = $1;
    }
    my @data;
    push(@data, "##sequence-region", $id, 1, $length);

    my $l = join(" ", @data);
    return $l;
}
#------------------------------------------------------------------------
sub print_txt {
    my $fh  = shift;
    my $str = shift;;
    
    if (defined($fh)){
	print $fh $str;
    }
    else {
	print $str;
    }
}
#------------------------------------------------------------------------
sub header {
	my $self = shift;

	my $build = $self->{build};
	my $h = "##gff-version 3";
	#removed for gmod
	#$h .= "\n##genome-build maker $build" if(defined $build);

	return $h;
    }
#------------------------------------------------------------------------
sub add_predictions {
    my $self  = shift;
    my $hits  = shift;
    
    open(my $ANN, '>>', $self->{ann_file})|| die "ERROR: Can't open annotation file\n\n";
    foreach my $p (@{$hits}){
       print_txt($ANN, hit_data($p, $self->seq_id));
    }
    close($ANN);
}
#------------------------------------------------------------------------
sub add_phathits {
   my $self  = shift;
   my $hits  = shift;
   
   open(my $ANN, '>>', $self->{ann_file})|| die "ERROR: Can't open annotation file\n\n";
    foreach my $h (@{$hits}){
       print_txt($ANN, hit_data($h, $self->seq_id));
    }
    close($ANN);
}
#------------------------------------------------------------------------
sub add_repeat_hits {
   my $self  = shift;
   my $hits  = shift;
   
   open(my $ANN, '>>', $self->{ann_file}) || die "ERROR: Can't open annotation file\n\n";
    foreach my $r (@{$hits}){
       print_txt($ANN, repeat_data($r, $self->seq_id));
    }
    close($ANN);
}
#------------------------------------------------------------------------
sub add_genes {
    my $self  = shift;
    my $genes = shift;
    
    open(my $ANN, '>>', $self->{ann_file}) || die "ERROR: Can't open annotation file\n\n";
    foreach my $g (@{$genes}){
       print_txt($ANN, gene_data($g, $self->seq_id));
    }
    close($ANN);
}
#------------------------------------------------------------------------
#------------------------------------------------------------------------
#----------------------------- FUNCTIONS --------------------------------
#------------------------------------------------------------------------
{
        my $id = -1;
sub get_id_cds {
        $id++;
        return 'cds:'.$id;
}
}
#------------------------------------------------------------------------

{
        my $id = -1;
sub get_id_hit {
        $id++;
        return 'hit:'.$id;
}
}
#------------------------------------------------------------------------
{
        my $id = -1;
sub get_id_hsp {
        $id++;
        return 'hsp:'.$id;
}
}
#------------------------------------------------------------------------

{
	my $id = -1;
sub get_id_gene {
	$id++;
	return 'gene:'.$id;
}
}
#------------------------------------------------------------------------
{
        my $id = -1;
sub get_id_mRNA {
        $id++;
        return 'mRNA:'.$id;
}
}
#------------------------------------------------------------------------
{
        my $id = -1;
sub get_id_exon {
        $id++;
        return 'exon:'.$id;
}
}
#------------------------------------------------------------------------
sub gene_data {
    my $g      = shift;
    my $seq_id = shift;
    
    my $g_name     = $g->{g_name};
    my $g_id       = $g->{g_id} || $g_name;
    my $g_s        = $g->{g_start};
    my $g_e        = $g->{g_end};
    my $g_strand   = $g->{g_strand} ==1 ? '+' : '-';
    my $t_structs  = $g->{t_structs};
    
    my @g_data;
    push(@g_data, $seq_id, 'maker', 'gene', $g_s, $g_e, '.', $g_strand, '.');
    my $attributes = 'ID='.$g_id.';Name='.$g_name.';';
    $attributes .= $g->{g_attrib} if($g->{g_attrib});
    $attributes .= ';' if($attributes !~ /\;$/);
    push(@g_data, $attributes); 
    
    my $g_l = join("\t", @g_data)."\n";
    
    my @transcripts = @{$g->{t_structs}};
    
    my %epl;
    my %cdss;
    foreach my $t (@transcripts){
	my $t_id = (split(/\s+/, $t->{t_id}))[0] || (split(/\s+/, $t->{t_name}))[0];
	my $t_l = get_transcript_data($t, $seq_id, $g_id);
	$g_l .= $t_l."\n";
	
	grow_exon_data_lookup($t->{hit}, $t_id, \%epl);
	grow_cds_data_lookup($t->{hit}, $t_id, \%cdss, $t->{t_offset}, $t->{t_end});
    }
    my $e_l = get_exon_data($seq_id, \%epl, 'maker');
    
    $g_l .= $e_l;
    
    my $c_l = get_cds_data($seq_id, \%cdss, 'maker');
    
    $g_l .= $c_l;
    return $g_l;
}
#------------------------------------------------------------------------
sub hit_data {
   my $h      = shift;
   my $seq_id = shift;

   my $h_str = $h->strand('query') == 1 ? '+' : '-';
   
   my ($h_s, $h_e) = PhatHit_utils::get_span_of_hit($h, 'query');
   
   ($h_s, $h_e) = ($h_e, $h_s) if $h_s > $h_e; 
   
   
   my ($t_s, $t_e) = PhatHit_utils::get_span_of_hit($h, 'hit');
   
   my $t_strand = $h->strand('hit') == -1  ? '-' : '+';
   
   ($t_s, $t_e) = ($t_e, $t_s) if $t_s > $t_e;
   
   my ($class, $type) = get_class_and_type($h, 'hit');
   
   my $h_n = $h->name();
   
   my $name = $h_n;
   $name =~ s/\s/_/g;

   my $AED .= sprintf '%.2f', $h->{_AED} if($h->{_AED});
   my $QI .= $h->{_QI} if($h->{_QI});

   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id);
   my $score = $h->score() || '.';
   $score .= 0 if $score eq '0.';
   
   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str, '.');
   my $attributes = 'ID='.$h_id.';Name='.$name.';';#.'Target='.$h_n.' '.$t_s.' '.$t_e.' '.$t_strand.';';
   $attributes .= '_AED='.$AED.';' if(defined($AED));
   $attributes .= '_QI='.$QI.';' if(defined($QI));
   $attributes .= $h->{-attrib} if($h->{-attrib});
   $attributes .= ';' if($attributes !~ /\;$/);
   my $h_l = join("\t", @h_data, $attributes)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id);
      $hsp_id =~ s/\s/_/g;
      my $hsp_l =
      get_hsp_data($hsp, $hsp_id, $seq_id, $h_id, $h_n);
      
      $h_l .= $hsp_l."\n";
      
   }
   
   return $h_l;
}
#------------------------------------------------------------------------
sub repeat_data {
   my $h      = shift;
   my $seq_id = shift;

   my $h_str = $h->strand('query') == 1 ? '+' : '-';
   #my $h_str = '+';
   
   my ($h_s, $h_e) = PhatHit_utils::get_span_of_hit($h, 'query');
   
   ($h_s, $h_e) = ($h_e, $h_s) if $h_s > $h_e; 
   
   
   my ($t_s, $t_e) = PhatHit_utils::get_span_of_hit($h, 'hit');
   
   my $t_strand = $h->strand('hit') == -1  ? '-' : '+';
   
   ($t_s, $t_e) = ($t_e, $t_s) if $t_s > $t_e;
   
   my ($class, $type) = get_class_and_type($h, 'hit');

   $class = "repeatrunner" if ($class eq 'blastx');
   
   my $h_n = $h->name();
   $h_n   =~ s/\s/_/g;
   
   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id);
   my $score = $h->score() || '.';
   $score .= 0 if $score eq '0.';

   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str, '.');
   my $attributes = 'ID='.$h_id.';Name='.$h_n.';Target='.$h_n.' '.$t_s.' '.$t_e.' '.$t_strand.';';
   $attributes .= $h->{-attrib} if($h->{-attrib});
   $attributes .= ';' if($attributes !~ /\;$/);
   my $h_l = join("\t", @h_data, $attributes)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id);
      $hsp_id =~ s/\s/_/g;
      my $hsp_l =
      get_repeat_hsp_data($hsp, $hsp_id, $seq_id, $h_id, $h_n);
      
      $h_l .= $hsp_l."\n";
      
   }
   
   return $h_l;
}
#------------------------------------------------------------------------
sub get_class_and_type {
    my $h = shift;
    my $k = shift;
    
    my ($class) = lc($h->algorithm);
    $class =~ s/^exonerate\:*\_*//;

    my $type;
    if    ($class =~ /^blastx$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif ($class =~ /^protein2genome$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part'; 
    }
    elsif($class =~ /^protein_gff\:/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif    ($class =~ /^tblastx$/i){
	$type = $k eq 'hit' ? 'translated_nucleotide_match' : 'match_part';
    }
    elsif    ($class =~ /^altest_gff/i){
	$type = $k eq 'hit' ? 'translated_nucleotide_match' : 'match_part';
    }
    elsif ($class =~ /^est2genome$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
    }
    elsif ($class =~ /^blastn$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part' ;
    }
    elsif ($class =~ /^est_gff\:/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
    }
    elsif ($class =~ /^snap_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^genemark_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^augustus_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^fgenesh_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^twinscan_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^jigsaw_*/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^pred_gff\:/i){
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif ($class =~ /^repeat_gff\:/i){
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    elsif ($class =~ /^blastx\:repeat/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif ($class =~ /^repeatmasker$/i){
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    elsif ($class =~ /^maker$/i){ #pasthrough maker annotation as evidence
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    else {
	die "unknown class in GFFV3::get_class_and_type $class ".ref($h)."\n";
    }
    
    return ($class, $type);
}
#------------------------------------------------------------------------
sub get_exon_data {
	my $seq_id = shift;
	my $epl    = shift;
	my $source = shift;

	my @uniques;
	my %seen;
	foreach my $e (@{$epl->{exons}}){
		my $nB = $e->nB('query');
		my $nE = $e->nE('query');
		push(@uniques, $e) unless defined($seen{$nB}{$nE});
		$seen{$nB}{$nE}++;
	}

	my $e_l = '';
	foreach my $e (@uniques){
		my $nB = $e->nB('query');
                my $nE = $e->nE('query');
		
		my $e_id = get_id_exon();
		my $e_n  = $e->name();
		
		my @t_ids = @{$epl->{t_ids}->{$nB}->{$nE}};
		
		my $p = join(',', @t_ids);
		$e_id = join(":", $t_ids[0], $e_id);
		
		($nB, $nE) = ($nE, $nB) if $nB > $nE;
		
		my $strand = $e->strand('query') == 1 ? '+': '-';

		my $score = $e->score() || '.';
		$score .= 0 if $score eq '0.';

		my @data;
		push(@data, $seq_id, $source, 'exon', $nB, $nE, $score, $strand, '.');
		my $nine = 'ID='.$e_id.';Parent='.$p .';';
		   $nine .= $e->{-attrib} if($e->{-attrib});
		   $nine .= ';' if($nine !~ /\;$/);
		push(@data, $nine); 

		$e_l .= join("\t", @data)."\n";
	}

	return $e_l;
}
#------------------------------------------------------------------------
sub get_cds_data {
        my $seq_id = shift;
        my $cdss    = shift;
        my $source = shift;

	#PostData($cdss);
	#die "UUUUUUUUUUUUUUUUUUU\n";
        my @uniques;
        my %seen;
        foreach my $e (@{$cdss->{cds}}){
                my $nB = $e->[0];
                my $nE = $e->[1];
                push(@uniques, $e) unless defined($seen{$nB}{$nE});
                $seen{$nB}{$nE}++;
        }
	
	#reverse for minus strand, important for phase calculation
	@uniques = reverse(@uniques) if(@uniques && $uniques[0][2]==-1);

        my $c_l = '';
	my $phase = 0;# + $fix;
        foreach my $e (@uniques){
                my $nB = $e->[0];
                my $nE = $e->[1];

                my $e_id = get_id_cds();
		my $e_n  = $e->[3];
		
                my @t_ids = @{$cdss->{t_ids}->{$nB}->{$nE}};

                my $p = join(',', @t_ids);
		$e_id = join(":", $t_ids[0], $e_id);
                ($nB, $nE) = ($nE, $nB) if $nB > $nE;

                my $strand = $e->[2] == 1 ? '+' : '-';

                my $score = '.';

                my @data;
                push(@data, $seq_id, $source, 'CDS', $nB, $nE, $score, $strand, $phase);
		my $nine = 'ID='.$e_id.';Parent='.$p.';';
		  #$nine .= $e->{-attrib} if($e->{-attrib});
		   $nine .= ';' if($nine !~ /\;$/);
                push(@data, $nine);

		# $phase = (3 - (($nE - $nB + 1) % 3)) % 3;
		$phase = ($phase - ($nE - $nB + 1)) % 3;
		
		#make sure CDS lines are ordered along plus strand
		#needed because of previous reveral for minus strand
		if($strand eq '+'){
		    $c_l .= join("\t", @data)."\n";
		}
		else{
		    $c_l = join("\t", @data)."\n".$c_l;
		}
        }

        return $c_l;
}
#------------------------------------------------------------------------
sub grow_exon_data_lookup {
	my $phat_hit    = shift;
	my $id          = shift;
	my $epl         = shift;
	
	foreach my $hsp ($phat_hit->hsps){
		my $nB = $hsp->nB('query');
		my $nE = $hsp->nE('query');

		push (@{$epl->{t_ids}->{$nB}->{$nE}}, $id);

		push(@{$epl->{exons}},  $hsp);
	}
}
#------------------------------------------------------------------------
sub grow_cds_data_lookup {
        my $phat_hit    = shift;
        my $id          = shift;
        my $cdss        = shift;
	my $offset      = shift;
	my $transl_end  = shift;

	my $hsp_start = 1;
	my $hsps = PhatHit_utils::sort_hits($phat_hit);

        foreach my $hsp (@$hsps){

                my $nB = $hsp->nB('query');
                my $nE = $hsp->nE('query');

		my $q_strand = $phat_hit->strand('query');

		my $hsp_end = abs($nE -$nB) + $hsp_start;

		#print "hsp_start:$hsp_start hsp_end:$hsp_end\n";
		#print "offset:$offset transl_end:$transl_end\n";

		my ($b, $e);
		if ($hsp_start <= $offset && $hsp_end > $offset){
			my $d = $offset - $hsp_start;
			$b = $q_strand == 1 ? $nB + ($d + 1) : $nB - ($d + 1);	
		}
		else {
			$b = $nB;
		} 

		if ($hsp_start < $transl_end && $hsp_end >= $transl_end){
			my $d = $hsp_end - $transl_end;
			$e = $q_strand == 1 ? $nE - ($d + 1) : $nE + ($d + 1); 
		}
		else {
			$e = $nE;
		}

		my $warn;

		if     ($b > $e && $hsp->strand('query') == 1){
			$warn = 1;
		}
		elsif ($b < $e && $hsp->strand('query') == -1){
			$warn = 1;
		}
		else {
			$warn = 0;
		}

		warn  "WARNING:B > E in GFFV3.pm strand:$q_strand\n" if $warn; 
		sleep 5 if $warn;

		unless ($hsp_end <= $offset || $hsp_start >= $transl_end || $warn){
		        #print "B:$b E:$e nB:$nB nE:$nE\n";
			#$hsp->show();
		        if($q_strand == 1){ #helps keep output in sorted order
			    push(@{$cdss->{cds}},  [$b, 
						    $e, 
						    $q_strand, 
						    $phat_hit->name()
						    ]);
			} 
			else{
			    unshift(@{$cdss->{cds}},  [$b,
						       $e,
						       $q_strand,
						       $phat_hit->name()
						       ]);
			}

			push(@{$cdss->{t_ids}->{$b}->{$e}},  $id);
		}
		$hsp_start = $hsp_end + 1; 
        }
	#PostData($cdss);
}
#------------------------------------------------------------------------
sub get_transcript_data {
	my $t      = shift;
	my $seq_id = shift;
	my $g_id   = shift;

	my $t_hit   = $t->{hit};
	my $t_seq   = $t->{t_seq};
	my $p_seq   = $t->{p_seq};
	my $t_off   = $t->{t_offset};
	my $t_end   = $t->{end};
	my $t_name  = $t->{t_name};
	my $t_id    = $t->{t_id};
	my $t_qi    = $t->{t_qi};
	my $AED     = $t->{AED};
	my $score   = $t->{score};

	#format informative name for GFF3
	$AED = sprintf '%.2f', $AED; # two decimal places

	my $t_s = $t_hit->strand('query') == 1 ? '+' : '-';

	my ($t_b, $t_e) = PhatHit_utils::get_span_of_hit($t_hit, 'query');

	($t_b, $t_e) = ($t_e, $t_b) if $t_b > $t_e;

	my @data;
	push(@data, $seq_id, 'maker', 'mRNA', $t_b, $t_e, '.', $t_s, '.');
	my $nine = 'ID='.$t_id.';Parent='.$g_id.';Name='.$t_name.';';
	   $nine .= '_AED='.$AED.';' if(defined($AED));
	   $nine .= '_QI='.$t_qi.';' if(defined($t_qi));
	   $nine .= $t_hit->{-attrib} if($t_hit->{-attrib});
	   $nine .= ';' if($nine !~ /\;$/);
	if($t->{hit}->{_Alias}){
	    if($nine =~ /Alias\=([^\;\n]+)/){
		my @keepers = (@{[split(',', $1)]}, @{$t->{hit}->{_Alias}});
		my %uniq;
		@uniq{@keepers} =();
		my $alias = join(',', keys %uniq);
		$nine =~ s/Alias\=[^\;\n]+\;*/Alias=$alias\;/;
	    }
	    else{
                my %uniq;
		@uniq{@{$t->{hit}->{_Alias}}} =();
                my $alias = join(',', keys %uniq);
                $nine .= "Alias=$alias\;";
	    }
	}

	push(@data, $nine);

	my $l = join("\t", @data);
	return $l;
}
#------------------------------------------------------------------------
sub get_hsp_data {
        my $hsp      = shift;
        my $hsp_id   = shift;
        my $seq_id   = shift;
        my $hit_id   = shift;
        my $hit_n    = shift;

        my $hsp_str  = $hsp->strand('query') ==  1 ? '+' : '-';
	my $t_strand = $hsp->strand('hit')   == -1 ? '-' : '+';

	#make Gap attribute
	my $q_string = $hsp->query_string();
	my $h_string = $hsp->hit_string();
	if($t_strand eq '-'){
	    $q_string = reverse($q_string);
	    $h_string = reverse($h_string);
	}
	my @gap = $hsp->cigar_string() =~ /([A-Z]\d+)/g;

	my $score = $hsp->score() || '.';
	$score .= 0 if $score eq '0.';

	my $nB = $hsp->nB('query');
	my $nE = $hsp->nE('query');

        ($nB, $nE) = ($nE, $nB) if $nB > $nE;

	my $tB = $hsp->nB('hit');
	my $tE = $hsp->nE('hit');

	($tB, $tE) = ($tE, $tB) if $tB > $tE;

	my ($class, $type) = get_class_and_type($hsp, 'hsp');

	  my $hsp_name = $hit_n;
             $hsp_name =~ s/\s/_/g;


	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id.';Name='.$hsp_name;
	   $nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE;
	   $nine .= ' '.$t_strand if($hsp->strand('hit'));
	   $nine .= ';';
	   $nine .= 'Gap='.join(' ', @gap).';' if(@gap);
	   $nine .= $hsp->{-attrib} if($hsp->{-attrib});
	   $nine .= ';' if($nine !~ /\;$/);
        my @data;
        push(@data, $seq_id, $class, $type, $nB, $nE, $score, $hsp_str, '.');
        push(@data, $nine);

        my $l = join("\t", @data);
        return $l;
}
#------------------------------------------------------------------------
sub get_repeat_hsp_data {
        my $hsp      = shift;
        my $hsp_id   = shift;
        my $seq_id   = shift;
        my $hit_id   = shift;
        my $hit_n    = shift;

        my $hsp_str  = $hsp->strand('query') ==  1 ? '+' : '-';
	#my $hsp_str  = '+';
	my $t_strand = $hsp->strand('hit')   == -1 ? '-' : '+';

	my $score = $hsp->score() || '.';
	$score .= 0 if $score eq '0.';

	my $nB = $hsp->nB('query');
	my $nE = $hsp->nE('query');

        ($nB, $nE) = ($nE, $nB) if $nB > $nE;

	my $tB = $hsp->nB('hit');
	my $tE = $hsp->nE('hit');

	($tB, $tE) = ($tE, $tB) if $tB > $tE;

	my ($class, $type) = get_class_and_type($hsp, 'hsp');
	$class = "repeatrunner" if ($class eq 'blastx');
	
	my $hsp_name = $hit_n;
	   $hsp_name =~ s/\s/_/g;

 
        my @data;
        push(@data, $seq_id, $class, $type, $nB, $nE, $score, $hsp_str, '.');
	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id.';Name='.$hsp_name;
	   $nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE.' '.$t_strand.';';
	   $nine .= $hsp->{-attrib} if($hsp->{-attrib});
	   $nine .= ';' if($nine !~ /\;$/);
        push(@data, $nine);

        my $l = join("\t", @data);
        return $l;
}
#------------------------------------------------------------------------
1;


