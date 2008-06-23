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

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- METHODS   ----------------------------------
#------------------------------------------------------------------------
sub new {
    my $self = {};
    bless $self;
    return $self;
}
#------------------------------------------------------------------------
sub fasta {
    my $self  = shift;
    my $fasta = Fasta::toFasta('>'.$self->seq_id, \(uc(${$self->seq})));
    return $$fasta;
}
#------------------------------------------------------------------------
sub seq_id {
    my $self   = shift;
    my $seq_id = shift;
    if (defined($seq_id)) {
	$self->{seq_id} = $seq_id;
    }
    else {
	return $self->{seq_id};
    }
}
#------------------------------------------------------------------------
sub seq {
    my $self  = shift;
    my $seq = shift;
    
    if (defined($seq)) {
	$self->{seq} = $seq;
    }
    else {
	return $self->{seq};
    }
}
#------------------------------------------------------------------------
sub print {
    my $self = shift;
    my $file = shift;
    my $fh;
    if (defined($file)){
	$fh = new FileHandle();
	$fh->open(">$file");		
    }
    print_txt($fh, $self->header."\n");
    print_txt($fh, $self->contig_comment."\n");
    print_txt($fh, $self->contig_line."\n"); 

    print_txt($fh, $self->genes); 
    print_txt($fh, $self->predictions);
    print_txt($fh, $self->repeat_hits);
    print_txt($fh, $self->phat_hits);   
    
    print_txt($fh, "##FASTA\n");
    print_txt($fh, $self->fasta);
    $fh->close() if defined($file);
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
    push(@data, 'ID='.$id.';Name='.$name);
    
    my $l = join("\t", @data);
    return $l;
}
#------------------------------------------------------------------------
sub contig_comment {
    my $self = shift;

    die "no contig seq in Dumper::GFFV3::contig_comment\n"
        unless defined($self->seq);

    my $seq = $self->seq();

    my $length = length($$seq);
    my $id     = $self->seq_id();
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
sub print_fld {
    my $data = shift;
    my $fh   = shift;
    
    die "LLLLL\n";
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
	
	my $h = "##gff-version 3";
	
	return $h;
    }
#------------------------------------------------------------------------
sub predictions {
    my $self  = shift;
    my $hits  = shift;
    
    if (defined($hits) && defined($hits->[0])){
       foreach my $p (@{$hits}){
	  #$self->{predictions} .= pred_data($p, $self->seq_id);
	  $self->{predictions} .= hit_data($p, $self->seq_id);
       }
    }
    else {
	return $self->{predictions} || '';
    }
}
#------------------------------------------------------------------------
sub phat_hits {
   my $self  = shift;
   my $hits  = shift;
   
   if (defined($hits) && defined($hits->[0])){
      foreach my $hit (@{$hits}){
	 $self->{hits} .= hit_data($hit, $self->seq_id);
      }
   }
   else {
      return $self->{hits} || '';
   }
}
#------------------------------------------------------------------------
sub repeat_hits {
   my $self  = shift;
   my $hits  = shift;
   
   if (defined($hits) && defined($hits->[0])){
      foreach my $hit (@{$hits}){
	 $self->{repeats} .= repeat_data($hit, $self->seq_id);
      }
   }
   else {
      return $self->{repeats} || '';
   }
}
#------------------------------------------------------------------------
sub genes {
    my $self  = shift;
    my $genes = shift;
    
    if (defined($genes) && defined($genes->[0])){
       foreach my $g (@{$genes}){
	  $self->{genes} .= gene_data($g, $self->seq_id);
       }
    }
    else {
	return $self->{genes} || '';
    }
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
sub pred_data {
    my $g      = shift;
    my $seq_id = shift;
    
    my $g_name     = $g->name();
    my $g_s        = $g->nB('query');
    my $g_e        = $g->nE('query');
    my $g_strand   = $g->strand('query') == 1 ? '+' : '-';
    
    $g_name =~ s/\s/-/g;
    my $g_id = get_id_gene();
    $g_id = join("-", "snap", $seq_id, $g_id);
    my $score = $g->score();
    
    ($g_s, $g_e) = ($g_e, $g_s) if $g_s > $g_e;


	my ($class, $type) = get_class_and_type($g, 'hit');
    
    my @g_data;
    push(@g_data, $seq_id, $class, 'gene', $g_s, $g_e, $score, $g_strand);
    push(@g_data, '.', 'ID='.$g_id.';Name='.$g_name);
    
    my $g_l = join("\t", @g_data)."\n";
    
    my @transcripts = $g;
    
    my %epl;
    foreach my $t (@transcripts){
	#my $t_id = get_id_mRNA();
	my $t_id = join(":", $t->name, get_id_mRNA());
	#-------
	my $t_s = $t->strand('query') == 1 ? '+' : '-';
	
	my ($t_b, $t_e) = PhatHit_utils::get_span_of_hit($t, 'query');
	
	($t_b, $t_e) = ($t_e, $t_b) if $t_b > $t_e;
	
	my ($p_b, $p_e) = PhatHit_utils::get_span_of_hit($t, 'hit');
	
	my $nine  = 'ID='.$t_id.';Parent='.$g_id.';Name='.$t->name;
	    $nine .= ';Target='.$t->name.' '.$p_b.' '.$p_e.' '.'+';
	
	my ($class, $type) = get_class_and_type($t, 'hsp');
	
	my @t_data;
	push(@t_data, $seq_id, $class, $type, $t_b, $t_e, '.', $t_s, '.');
        	push(@t_data, $nine);

		my $t_l = join("\t", @t_data)."\n";
	#--------
	$g_l .= $t_l;
	
	grow_exon_data_lookup($t, $t_id, \%epl);
    }
    my $e_l = get_exon_data($seq_id, \%epl, 'snap');
    
    $g_l .= $e_l;
    
    return $g_l;
}
#------------------------------------------------------------------------
sub gene_data {
    my $g      = shift;
    my $seq_id = shift;
    
    my $g_name     = join("-", "maker", $seq_id, (split("-", $g->{g_name}))[1..2]);
    my $g_s        = $g->{g_start};
        my $g_e        = $g->{g_end};
        my $g_strand   = $g->{g_strand} ==1 ? '+' : '-';
    my $t_structs  = $g->{t_structs};
    
    #my $g_id = get_id_gene();
    my $g_id = $g_name;
    my @g_data;
    push(@g_data, $seq_id, 'maker', 'gene', $g_s, $g_e, '.', $g_strand);
    push(@g_data, '.', 'ID='.$g_id.';Name='.$g_name); 
    
    my $g_l = join("\t", @g_data)."\n";
    
    my @transcripts = @{$g->{t_structs}};
    
    my %epl;
    my %cdss;
    foreach my $t (@transcripts){
	#my $t_id = get_id_mRNA();
	my $t_id = (split(/\s+/, $t->{t_name}))[0];
	my $t_l = 
	    get_transcript_data($t, $t_id, $seq_id, $g_id, $g_name);
	
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
   
   my $h_n = $class eq 'repeatmasker' && $type eq 'match' 
       ? $h->hsp('best')->name() : $h->name();
   
   $h_n   =~ s/\s/-/g;
   
   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id);
   my $score = $h->significance() || 'NA';
   $score .= 0 if $score eq '0.';
   $score = '.' if $score eq 'NA';

   my $alt_score = $h->score() || '.';
   $score = $alt_score
       if ($score eq '.' && $alt_score =~ /\d/);

   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str);
   push(@h_data, '.', 'ID='.$h_id.';Name='.$h_n.';Target='.$h_n.' '.$t_s.' '.$t_e.' '.$t_strand);
   my $h_l = join("\t", @h_data)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id);
      $hsp_id =~ s/\s/-/g;
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

   $class = "blastx:repeatmask" if ($class eq 'blastx');
   
   my $h_n = $class eq 'repeatmasker' && $type eq 'match' 
       ? $h->hsp('best')->name() : $h->name();
   
   $h_n   =~ s/\s/-/g;
   
   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id);
   my $score = $h->significance() || 'NA';
   $score .= 0 if $score eq '0.';
   $score = '.' if $score eq 'NA';

   my $alt_score = $h->score() || '.';
   $score = $alt_score
       if ($score eq '.' && $alt_score =~ /\d/);

   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str);
   push(@h_data, '.', 'ID='.$h_id.';Name='.$h_n.';Target='.$h_n.' '.$t_s.' '.$t_e.' '.$t_strand);
   my $h_l = join("\t", @h_data)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id);
      $hsp_id =~ s/\s/-/g;
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

	my ($class) = ref($h) =~ /.*::(\S+)$/;

	my $type;
	if    ($class eq 'blastx'){
		$type = $k eq 'hit' ? 'protein_match' : 'match_part';
	}elsif    ($class eq 'tblastx'){
		$type = $k eq 'hit' ? 'translated_nucleotide_match' : 'match_part';
	}
	elsif ($class eq 'protein2genome'){
		$type = $k eq 'hit' ? 'protein_match' : 'match_part'; 
	}
        elsif ($class eq 'est2genome'){
                $type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
        }
        elsif ($class eq 'blastn'){
                $type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part' ;
        }
        elsif (ref($h)  =~ /snap/){
                $type = $k eq 'hit' ? 'match' : 'match_part' ;
		$class= 'snap';
        }
        elsif (ref($h)  =~ /repeatmasker/){
                $type = $k eq 'hit' ? 'match' : 'match_part';
		$class= 'repeatmasker';
        }
	else {
		die "unknown class in GFFV3::get_class_and_type:".ref($h)."\n";
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

		my $score = $source eq 'maker' ? '.' : $e->score();

		my @data;
		push(@data, $seq_id, $source, 'exon', $nB, $nE, $score);
		push(@data, $strand, '.','ID='.$e_id.';Parent='.$p); 

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

        my $c_l = '';
	my $phase = 0;
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
                push(@data, $seq_id, $source, 'CDS', $nB, $nE, $score);
                push(@data, $strand, $phase,'ID='.$e_id.';Parent='.$p);

		$phase = (3 - (($nE - $nB + 1) % 3)) % 3;

                $c_l .= join("\t", @data)."\n";
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
        foreach my $hsp ($phat_hit->hsps){

                my $nB = $hsp->nB('query');
                my $nE = $hsp->nE('query');

		my $q_strand = $hsp->strand('query');

		my $hsp_end = abs($nE -$nB) + $hsp_start;

		#print "hsp_start:$hsp_start hsp_end:$hsp_end\n";
		#print "offset:$offset transl_end:$transl_end\n";

		my ($b, $e);
		if ($hsp_start <= $offset && $hsp_end >= $offset){
			my $d = $offset  - $hsp_start;
			$b = $q_strand == 1 ? $nB + ($d + 1) : $nB - ($d + 1);	
		}
		else {
			$b = $nB;	
		} 

		if ($hsp_start <= $transl_end && $hsp_end >= $transl_end){
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

		warn  "WARNING:B > E in GFFV3.pm strand:".$hsp->strand('query')."\n" if $warn; 
		sleep 5 if $warn;

		unless ($hsp_end < $offset || $hsp_start > $transl_end || $warn){
			#print "B:$b E:$e nB:$nB nE:$nE\n";
			#$hsp->show();
			push(@{$cdss->{cds}},  [$b, 
			                        $e, 
			                        $hsp->strand('query'), 
			                        $hsp->name(),
			                        ]); 

			push(@{$cdss->{t_ids}->{$b}->{$e}},  $id);
		}
		$hsp_start = $hsp_end + 1; 
        }
	#PostData($cdss);
}
#------------------------------------------------------------------------
sub get_transcript_data {
	my $t      = shift;
	my $t_id   = shift;
	my $seq_id = shift;
	my $g_id   = shift;
	my $g_n    = shift;

	my $t_hit   = $t->{hit};
	my $t_seq   = $t->{t_seq};
	my $p_seq   = $t->{p_seq};
	my $t_off   = $t->{t_offset};
	my $t_end   = $t->{end};
	my $t_name  = $t->{t_name};
	my $t_qi    = $t->{qi};

	my $t_s = $t_hit->strand('query') == 1 ? '+' : '-';

	my ($t_b, $t_e) = PhatHit_utils::get_span_of_hit($t_hit, 'query');

	($t_b, $t_e) = ($t_e, $t_b) if $t_b > $t_e;

	my @data;
	push(@data, $seq_id, 'maker', 'mRNA', $t_b, $t_e, '.', $t_s, '.');
	push(@data, 'ID='.$t_id.';Parent='.$g_id.';Name='.$t_name);


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

	my $score = $hsp->significance() || 'NA';
	$score .= 0 if $score eq '0.';
	$score = '.' if $score eq 'NA';
	
	my $alt_score = $hsp->score() || '.';
	$score = $alt_score
	    if ( $score eq '.' && $alt_score =~ /\d/);


	my $nB = $hsp->nB('query');
	my $nE = $hsp->nE('query');

        ($nB, $nE) = ($nE, $nB) if $nB > $nE;

	my $tB = $hsp->nB('hit');
	my $tE = $hsp->nE('hit');

	($tB, $tE) = ($tE, $tB) if $tB > $tE;

	my ($class, $type) = get_class_and_type($hsp, 'hsp');

	  my $hsp_name = $hsp->name();
             $hsp_name =~ s/\s/-/g;


	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id.';Name='.$hsp_name;
	   $nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE.' '.$t_strand;
 
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

	my $score = $hsp->significance() || 'NA';
	$score .= 0 if $score eq '0.';
	$score = '.' if $score eq 'NA';
	
	my $alt_score = $hsp->score() || '.';
	$score = $alt_score
	    if ( $score eq '.' && $alt_score =~ /\d/);


	my $nB = $hsp->nB('query');
	my $nE = $hsp->nE('query');

        ($nB, $nE) = ($nE, $nB) if $nB > $nE;

	my $tB = $hsp->nB('hit');
	my $tE = $hsp->nE('hit');

	($tB, $tE) = ($tE, $tB) if $tB > $tE;

	my ($class, $type) = get_class_and_type($hsp, 'hsp');
	$class = "blastx:repeatmask" if ($class eq 'blastx');
	
	my $hsp_name = $hsp->name();
	   $hsp_name =~ s/\s/-/g;

 
        my @data;
        push(@data, $seq_id, $class, $type, $nB, $nE, $score, $hsp_str, '.');
	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id.';Name='.$hsp_name;
	$nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE.' '.$t_strand;
        push(@data, $nine);

        my $l = join("\t", @data);
        return $l;
}
#------------------------------------------------------------------------
1;


