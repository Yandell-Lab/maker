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
use File::NFSLock;
use GI;
use Carp;

@ISA = qw(
       );
#------------------------------------------------------------------------
#----------------------------- METHODS ----------------------------------
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
   my $def_file = "$t_dir/$name.def";
   my $ann_file = "$t_dir/$name.ann";
   my $seq_file = "$t_dir/$name.seq";

   #use hard linking to force clear NFS cache
   my $link = (GI::is_NFS_mount($t_dir) == 1);
   my $tdef_file = ($link) ? "$def_file.tmp" : $def_file;
   my $tann_file = ($link) ? "$ann_file.tmp" : $ann_file;
   my $tseq_file = ($link) ? "$seq_file.tmp" : $seq_file;

   #continuous try because of weird NFS behavior with false success
   my $success = 0;
   while(! $success){
       open(my $DEF, "> $tdef_file") || confess "ERROR: Could not open file: $tdef_file\n$!\n";
       print_txt($DEF, $self->header."\n");
       close($DEF);
       #link is atomic on most NFS implementations
       unlink($def_file) if($link);
       $success = link("$tdef_file", $def_file) if($link);
       $success = -f $def_file if($success || !$link);
       $success = (stat _)[3] == 2 if($success && $link);
       unlink($tdef_file) if($link);
   }
   $success = 0;
   while(! $success){
       open(my $ANN, "> $tann_file") || confess "ERROR: Could not open file: $tann_file\n$!\n";
       close($ANN);
       #link is atomic on most NFS implementations
       unlink($ann_file) if($link);
       $success = link("$tann_file", $ann_file) if($link);
       $success = -f $ann_file  if($success || !$link);
       $success = (stat _)[3] == 2 if($success && $link);
       unlink($tann_file) if($link);
   }
   $success = 0;
   while(! $success){
       open(my $SEQ, "> $tseq_file") || confess "ERROR: Could not open file: $tseq_file\n$!\n";
       print_txt($SEQ, "##FASTA\n");
       close($SEQ);
       #link is atomic on most NFS implementations
       unlink($seq_file) if($link);
       $success = link("$tseq_file", $seq_file) if($link);
       $success = -f $seq_file if($success || !$link);
       $success = (stat _)[3] == 2 if($success && $link);
       unlink($tseq_file) if($link);
   }

   $self->{gff_file} = $gff_file;
   $self->{def_file} = $def_file;
   $self->{ann_file} = $ann_file;
   $self->{seq_file} = $seq_file;
}
#------------------------------------------------------------------------
sub resolved_flag {
    my $self   = shift;

    my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
    while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
    open(my $ANN, ">>", $self->{ann_file}) || confess "ERROR: Can't open annotation file\n\n";
    print_txt($ANN, "###\n");
    close($ANN);
    $lock->unlock;
}
#------------------------------------------------------------------------
sub set_current_contig {
    my $self   = shift;

    my $flag = 0;
    $flag = 1 if (defined $self->{seq_id});
    my $id  = shift;
    my $seq = shift; #can be left blank and no fasta or contig line will be added

    #escape seqid per gff standards
    $self->{seq_id} = uri_escape($id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|');

    if(!$self->{SEEN}{$self->{seq_id}} && $seq){
	my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
	while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}

	my $seq_ref;
	if($seq && (! ref($seq) || ref($seq) eq 'SCALAR')){
	    $seq = $$seq while(ref($seq) eq 'REF');
	    $seq_ref = (ref($seq) eq '') ? \$seq : $seq;
	    $self->{seq_length} = length($$seq_ref);
	}
	elsif($seq){
	    $self->{seq_length} = $seq->length;
	    $seq_ref = \ ($seq->seq);
	    $seq_ref = $$seq_ref while(ref($seq_ref) eq 'REF');
	}

	open(my $SEQ, ">>", $self->{seq_file}) || confess "ERROR: opening fasta for GFF3\n\n";
	print $SEQ ">".uri_unescape($self->{seq_id})."\n"; #fasta doesn't require escaped seq_ids
	my $offset = 0;
	my $width = 60;
	while($offset < $self->{seq_length}){
	    my $l = substr($$seq_ref, $offset, $width);
	    print $SEQ $l."\n";
	    $offset += $width;
	}
	close($SEQ);

	$self->{SEEN}{$self->{seq_id}}++;

	#-skip adding the optional sequence-region line because many programs
	#-do not handle it correctly or consistently - 06/05/2010
	#open(my $DEF, ">>", $self->{def_file}) || confess "ERROR: Can't open definition file\n\n";
	#print_txt($DEF, $self->contig_comment."\n");
	#close($DEF);
	
	open(my $ANN, ">>", $self->{ann_file}) || confess "ERROR: Can't open annotation file\n\n";
	print_txt($ANN, "###\n") if($flag);
	print_txt($ANN, $self->contig_line."\n"); 
	close($ANN);

	$lock->unlock;
    }
}
#------------------------------------------------------------------------
sub seq_length {
    my $self   = shift;

    return $self->{seq_length} || undef;
}

#------------------------------------------------------------------------
sub seq_id {
    my $self   = shift;

    return $self->{seq_id} || undef;
}
#------------------------------------------------------------------------
sub gff_file {
    my $self = shift;

    return $self->{gff_file};
}
#------------------------------------------------------------------------
sub finalize {
    my $self = shift;

    my $gff_f = $self->{gff_file};
    my $def_f = $self->{def_file};
    my $ann_f = $self->{ann_file};
    my $seq_f = $self->{seq_file};

    #already finalized
    return if(-f $gff_f && ! -f $def_f && ! -f $ann_f && ! -f $seq_f);
    sleep 10 if(! -f $def_f && ! -f $ann_f && ! -f $seq_f); # wait because of slow NFS

    open(my $DEF, ">> $def_f") || confess "ERROR: Can't open def file: $def_f\n$!\n\n";
    open(my $SEQ, "< $seq_f") || confess "ERROR: Can't open seq file: $seq_f\n$!\n\n";
    open(my $ANN, "< $ann_f") || confess "ERROR: Can't open annotation file: $ann_f\n$!\n\n";

    while(defined(my $line = <$ANN>)){
        print $DEF $line;
    }
    close($ANN);

    my $buf = <$SEQ>;
    if($buf !~ /^\#\#FASTA/){
	confess "ERROR: There was a problem in the writing the fasta entry\n".
	    "Either no sequence was given, or there was an error in writing\n\n";
    }
    if((my $line = <$SEQ>) =~ /^>/){
	print $DEF $buf . $line;
	while(defined($line = <$SEQ>)){
	    print $DEF $line;
	}
    }
    close($SEQ);
    close($DEF);

    unlink($seq_f);
    unlink($ann_f);
    move($def_f, $gff_f);

    confess "ERROR: GFF3 file not created\n" if(! -e $gff_f);
}
#------------------------------------------------------------------------
sub merge {
    my $self = shift;
    my $files = shift;

    return if(! @$files);

    my $def_f = $self->{def_file};
    my $ann_f = $self->{ann_file};
    my $seq_f = $self->{seq_file};

    open(my $DEF, ">> $def_f")|| confess "ERROR: Can't open def file: $def_f\n$!\n\n";
    open(my $SEQ, ">> $seq_f")|| confess "ERROR: Can't open seq file: $seq_f\n$!\n\n";
    open(my $ANN, ">> $ann_f")|| confess "ERROR: Can't open annotation file: $ann_f\n$!\n\n";

    foreach  my $file (@$files){
	my $FH = $ANN;
	open(my $IN, "< $file") || confess "ERROR: Could not open file \'$file\'\n";
	print $ANN "\###\n";
	while (defined(my $line = <$IN>)){
	    if ($line =~ /^\#\#gff-version 3/){
		next;
	    }
	    elsif($line =~ /^\#\#genome-build/){
		next;
	    }
	    elsif($line =~ /^\#\#sequence-region\s+(^[\s]+)/){
		next if($self->{SEEN}{$1});
		$self->{SEEN}{$1}++;
		print $DEF $line;
		next;
            }
	    elsif($line =~ /^\#\#FASTA/){
		$FH = $SEQ;
		next;
	    }
	    elsif($line =~ /^>/){
		$FH = $SEQ;
	    }
	    elsif($line =~ /^[^\s\t\#\>\n]+\t([^\t]+)\t/){
		$FH = $ANN;

		my ($source) = $line =~ /^([^\t]+)/;
		my ($id) = $line =~ /ID\=([^\;\n]+)/;
		#my ($parent) = $line =~ /Parent\=([^\;\n]+)/;
		if($id && $source eq $id){
		    next if($self->{SEEN}{$id});
		    $self->{SEEN}{$id}++;
		}
	    }

	    chomp $line; #chomp to remove empty lines
	    print $FH $line."\n" if($line);
	}
	close($IN);
    }
    close($DEF);
    close($ANN);
    close($SEQ);
}

#------------------------------------------------------------------------
sub contig_line {
    my $self = shift;
    
    confess "no contig seq in Dumper::GFFV3::contig_line\n"
	unless defined($self->seq_id);
    
    my $length = $self->seq_length();
    my $id     = $self->seq_id();
    my $name   = $id;
    my @data;
    push(@data, $id, '.', 'contig', 1, $length, '.','.','.');
    push(@data, 'ID='.$id.';Name='.$name);
    
    my $l = join("\t", @data);
    return $l;
}
#------------------------------------------------------------------------
sub contig_comment {
    my $self = shift;

    confess "no contig seq in Dumper::GFFV3::contig_comment\n"
        unless defined($self->seq_id);

    my $length = $self->seq_length();
    my $id     = $self->seq_id();
    my $name   = $id;
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
sub add_lines {
    my $self  = shift;
    my $lines = shift;

    return if(! $lines || ! @$lines);

    my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
    while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
    open(my $ANN, '>>', $self->{ann_file})|| confess "ERROR: Can't open annotation file\n\n";
    print_txt($ANN, join("\n", @$lines));    
    close($ANN);
    $lock->unlock;
}
#------------------------------------------------------------------------
sub add_predictions {
    my $self  = shift;
    my $hits  = shift;
    my $uid   = shift || $$;

    return if(! $hits || ! @$hits);

    my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
    while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
    open(my $ANN, '>>', $self->{ann_file})|| confess "ERROR: Can't open annotation file\n\n";
    foreach my $p (@{$hits}){
       print_txt($ANN, hit_data($p, $self->seq_id, $uid));
    }
    close($ANN);
    $lock->unlock;
}
#------------------------------------------------------------------------
sub add_phathits {
   my $self  = shift;
   my $hits  = shift;
   my $uid   = shift || $$;

   return if(! $hits || ! @$hits);

   my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
   while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
   open(my $ANN, '>>', $self->{ann_file})|| confess "ERROR: Can't open annotation file\n\n";
   foreach my $h (@{$hits}){
      print_txt($ANN, hit_data($h, $self->seq_id, $uid));
   }
   close($ANN);
   $lock->unlock;
}
#------------------------------------------------------------------------
sub add_repeat_hits {
   my $self  = shift;
   my $hits  = shift;
   my $uid   = shift || $$;

   return if(! $hits || ! @$hits);  

   my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
   while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
   open(my $ANN, '>>', $self->{ann_file}) || confess "ERROR: Can't open annotation file\n\n";
   foreach my $r (@{$hits}){
      print_txt($ANN, repeat_data($r, $self->seq_id, $uid));
   }
   close($ANN);
   $lock->unlock;
}
#------------------------------------------------------------------------
sub add_genes {
    my $self  = shift;
    my $genes = shift;

    return if(! $genes || ! @$genes);    

    my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
    while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
    open(my $ANN, '>>', $self->{ann_file}) || confess "ERROR: Can't open annotation file\n\n";
    foreach my $g (@{$genes}){
       print_txt($ANN, gene_data($g, $self->seq_id));
    }
    close($ANN);
    $lock->unlock;
}
#------------------------------------------------------------------------
sub add_ncgenes {
    my $self  = shift;
    my $genes = shift;

    return if(! $genes || ! @$genes);    

    my $lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30);
    while(!$lock || !$lock->still_mine){$lock = new File::NFSLock($self->{ann_file}, 'EX', 1800, 30)}
    open(my $ANN, '>>', $self->{ann_file}) || confess "ERROR: Can't open annotation file\n\n";
    foreach my $g (@{$genes}){
       print_txt($ANN, ncgene_data($g, $self->seq_id));
    }
    close($ANN);
    $lock->unlock;
}
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

    #maximum of supported characters for GFF3 attributes using ASCII encoding
    #\x21-\x24\x27-\x2B\x2D-\x3A\x3C\x3E-\x7E
    #I really only allow a subset of these though for portability

    $g_name = uri_escape($g_name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
    $g_id = uri_escape($g_id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards

    my @g_data;
    push(@g_data, $seq_id, 'maker', 'gene', $g_s, $g_e, '.', $g_strand, '.');
    my $attributes = 'ID='.$g_id.';Name='.$g_name;
    $attributes .= ';'.$g->{g_attrib} if($g->{g_attrib});
    $attributes =~  s/\;$//;
    push(@g_data, $attributes); 
    
    my $g_l = join("\t", @g_data)."\n";
    
    my @transcripts = @{$g->{t_structs}};
    
    my %epl;
    my $c_l;
    foreach my $t (@transcripts){
	my $t_id = (split(/\s+/, $t->{t_id}))[0] || (split(/\s+/, $t->{t_name}))[0];
	$t_id = uri_escape($t_id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards 
	my $t_l = get_transcript_data($t, $seq_id, $g_id);
	$g_l .= $t_l."\n";
	
	my %cdss;
	grow_exon_data_lookup($t->{hit}, $t_id, \%epl);
	grow_cds_data_lookup($t->{hit}, $t_id, \%cdss, $t->{t_offset}, $t->{t_end});

	$c_l .= get_cds_data($seq_id, \%cdss, 'maker');
    }

    $g_l .= get_exon_data($seq_id, \%epl, 'maker');
    $g_l .= $c_l;

    return $g_l;
}
#------------------------------------------------------------------------
sub ncgene_data {
    my $g      = shift;
    my $seq_id = shift;
    
    my $g_name     = $g->{g_name};
    my $g_id       = $g->{g_id} || $g_name;
    my $g_s        = $g->{g_start};
    my $g_e        = $g->{g_end};
    my $g_strand   = $g->{g_strand} ==1 ? '+' : '-';
    my $t_structs  = $g->{t_structs};

    $g_name = uri_escape($g_name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
    $g_id = uri_escape($g_id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards

    my @g_data;
    push(@g_data, $seq_id, 'maker', 'gene', $g_s, $g_e, '.', $g_strand, '.');
    my $attributes = 'ID='.$g_id.';Name='.$g_name;
    $attributes .= ';'.$g->{g_attrib} if($g->{g_attrib});
    $attributes =~  s/\;$//;
    push(@g_data, $attributes); 
    
    my $g_l = join("\t", @g_data)."\n";
    
    my @transcripts = @{$g->{t_structs}};
    
    my %epl;
    my $c_l;
    foreach my $t (@transcripts){
	my $t_id = (split(/\s+/, $t->{t_id}))[0] || (split(/\s+/, $t->{t_name}))[0];
	$t_id = uri_escape($t_id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
	my $t_l = get_transcript_data($t, $seq_id, $g_id);
	$g_l .= $t_l."\n";
	
	my %cdss;
	grow_exon_data_lookup($t->{hit}, $t_id, \%epl);
	#grow_cds_data_lookup($t->{hit}, $t_id, \%cdss, $t->{t_offset}, $t->{t_end});
	#$c_l .= get_cds_data($seq_id, \%cdss, 'maker');
    }

    $g_l .= get_exon_data($seq_id, \%epl, 'maker');
    $g_l .= $c_l;

    return $g_l;
}
#------------------------------------------------------------------------
sub hit_data {
   my $h      = shift;
   my $seq_id = shift;
   my $uid    = shift || $$;

   my $h_str = $h->strand('query') == 1 ? '+' : '-';
   
   my ($h_s, $h_e) = PhatHit_utils::get_span_of_hit($h, 'query');
   
   ($h_s, $h_e) = ($h_e, $h_s) if $h_s > $h_e; 
   
   
   my ($t_s, $t_e) = PhatHit_utils::get_span_of_hit($h, 'hit');
   
   my $t_strand = $h->strand('hit') == -1  ? '-' : '+';
   
   ($t_s, $t_e) = ($t_e, $t_s) if $t_s > $t_e;
   
   my ($class, $type) = get_class_and_type($h, 'hit');
   
   my $name = $h->name();
   $name = uri_escape($name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards

   my $AED = sprintf '%.2f', $h->{_AED} if(defined ($h->{_AED}));
   my $eAED = sprintf '%.2f', $h->{_eAED} if(defined ($h->{_eAED}));
   my $QI = $h->{_QI} if(defined($h->{_QI}));

   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id, $uid);
   my $score = $h->score() || '.';
   $score .= 0 if $score eq '0.';
   
   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str, '.');
   my $attributes = 'ID='.$h_id.';Name='.$name;
   $attributes .= ';_AED='.$AED if(defined($AED));
   $attributes .= ';_eAED='.$eAED if(defined($eAED));
   $attributes .= ';_QI='.$QI if(defined($QI));
   $attributes .= ';'.$h->{-attrib} if($h->{-attrib});
   $attributes =~  s/\;$//;
   my $h_l = join("\t", @h_data, $attributes)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id, $uid);
      $hsp_id =~ s/\s/_/g;
      my $hsp_l =
      get_hsp_data($hsp, $hsp_id, $seq_id, $h_id, $name, $h);
      
      $h_l .= $hsp_l."\n";
      
   }
   
   return $h_l;
}
#------------------------------------------------------------------------
sub repeat_data {
   my $h      = shift;
   my $seq_id = shift;
   my $uid    = shift || $$;

   my $h_str = $h->strand('query') == 1 ? '+' : '-';
   #my $h_str = '+';
   
   my ($h_s, $h_e) = PhatHit_utils::get_span_of_hit($h, 'query');
   
   ($h_s, $h_e) = ($h_e, $h_s) if $h_s > $h_e; 
   
   
   my ($t_s, $t_e) = PhatHit_utils::get_span_of_hit($h, 'hit');
   
   my $t_strand = $h->strand('hit') == -1  ? '-' : '+';
   
   ($t_s, $t_e) = ($t_e, $t_s) if $t_s > $t_e;
   
   my ($class, $type) = get_class_and_type($h, 'hit');

   $class =~ s/^(blastx|rapsearch)\:repeat.*?$/repeatrunner/;
   $class =~ s/^(blastx|rapsearch)/repeatrunner/;

   my $h_n = $h->name();
   $h_n = uri_escape($h_n, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
   
   my $h_id = get_id_hit();
   $h_id = join(":", $seq_id, $h_id, $uid);
   my $score = $h->score() || '.';
   $score .= 0 if $score eq '0.';

   my @h_data;
   push(@h_data, $seq_id, $class, $type, $h_s, $h_e, $score, $h_str, '.');
   my $attributes = 'ID='.$h_id.';Name='.$h_n.';Target='.$h_n.' '.$t_s.' '.$t_e.' '.$t_strand;
   $attributes .= ';'.$h->{-attrib} if($h->{-attrib});
   $attributes =~  s/\;$//;
   my $h_l = join("\t", @h_data, $attributes)."\n";
   
   my $sorted = PhatHit_utils::sort_hits($h);
   
   foreach my $hsp (@{$sorted}){
      my $hsp_id = get_id_hsp();
      $hsp_id = join(":", $seq_id, $hsp_id, $uid);
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
    $class =~ s/^exonerate\:*\_*est2genome$/est2genome/;
    $class =~ s/^exonerate\:*\_*cdna2genome$/cdna2genome/;
    $class =~ s/^exonerate\:*\_*protein2genome$/protein2genome/;

    my $type;
    if($class =~ /^blastx(\:.*)?$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif($class =~ /^rapsearch(\:.*)?$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif($class =~ /^protein2genome(\:.*)?$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part'; 
    }
    elsif($class =~ /^protein_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif($class =~ /^tblastx(\:.*)?$/i){
	$type = $k eq 'hit' ? 'translated_nucleotide_match' : 'match_part';
    }
    elsif($class =~ /^altest_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'translated_nucleotide_match' : 'match_part';
    }
    elsif($class =~ /^est2genome(\:.*)?$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
    }
    elsif($class =~ /^cdna2genome(\:.*)?$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
    }
    elsif($class =~ /^blastn(\:.*)?$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part' ;
    }
    elsif($class =~ /^est_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'expressed_sequence_match' : 'match_part';
    }
    elsif($class =~ /^snap(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^trnascan(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^snoscan(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^genemark(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^evm(_masked)?(\:.*)?$/i){ 
        $class = lc($h->algorithm);
        $type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^augustus(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^fgenesh(_masked)?(\:.*)?$/i){
	$class = lc($h->algorithm);
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^pred_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^model_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'match' : 'match_part' ;
    }
    elsif($class =~ /^repeat_gff(\:.*)?$/i){
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    elsif($class =~ /^repeatrunner(\:.*)?$/i){
	$type = $k eq 'hit' ? 'protein_match' : 'match_part';
    }
    elsif ($class =~ /^repeatmasker(\:.*)?$/i){
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    elsif ($class =~ /^maker(\:.*)?$/i){ #pasthrough maker annotation as evidence
	$type = $k eq 'hit' ? 'match' : 'match_part';
    }
    else {
	confess "unknown class in GFFV3::get_class_and_type $class ".ref($h)."\n";
    }
    
    $class .= ':'.$h->{_label} if($h->{_label});

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
	for(my $i = 0; $i < @uniques; $i++){
	    my $e = $uniques[$i];
		my $nB = $e->nB('query');
                my $nE = $e->nE('query');
		
	        #my $e_id = get_id_exon();
	        my $e_id = $i+1; #exon 1 is labeled with ID of exon:1
		my $e_n  = $e->name();
		
		my @t_ids = @{$epl->{t_ids}->{$nB}->{$nE}};
		
		my $p = join(',', @t_ids);
		$e_id = join(":", $t_ids[0], $e_id);
		
		($nB, $nE) = ($nE, $nB) if $nB > $nE;
		
		my $strand = $e->strand('query') == 1 ? '+': '-';

		my $score = '.'; #exon scores cause chado warnings/errors

		my @data;
		push(@data, $seq_id, $source, 'exon', $nB, $nE, $score, $strand, '.');
		my $nine = 'ID='.$e_id.';Parent='.$p;
		   $nine .= ';'.$e->{-attrib} if($e->{-attrib});
		   $nine =~  s/\;$//;
		push(@data, $nine); 

		#make sure exon lines are ordered relative top plus strand
		#needed because of previous reversal for minus strand
		if($strand eq '+'){
		   $e_l .= join("\t", @data)."\n";
		}
		else{
		   $e_l = join("\t", @data)."\n".$e_l;
		}
	}

	return $e_l;
}
#------------------------------------------------------------------------
sub get_cds_data {
   my $seq_id = shift;
   my $cdss    = shift;
   my $source = shift;
   
   my @labels = qw(five_prime_UTR start_codon CDS stop_codon three_prime_UTR);
   my $c_l = '';
   my $cds_phase = 0;
   my $start_phase = 0;
   my $stop_phase = 0;
   foreach my $label (@labels){
      foreach my $e (@{$cdss->{$label}}){
	 my $nB = $e->[0];
	 my $nE = $e->[1];
	 my $strand = $e->[2] == 1 ? '+' : '-';
	 my $e_n  = $e->[3];
	 
	 my @t_ids = @{$cdss->{t_ids}{$label}{$nB}{$nE}};
	 
	 #my $p = join(',', @t_ids);
	 my $p = $t_ids[0]; #multiple parents not correct for CDS
	 my $e_id = join(":", $t_ids[0], lc($label));
	 ($nB, $nE) = ($nE, $nB) if $nB > $nE;
	 
	 my $score = '.';
	      
	 my $phase = '.';
	 my $len = abs($nE - $nB) + 1;
	 if($label eq 'CDS'){
	    $phase = $cds_phase;
	    $cds_phase = ($cds_phase - $len) % 3;
	 }
	 elsif($label eq 'start_codon'){
	    $phase = $start_phase;
	    $start_phase = ($start_phase - $len) % 3;
	 }
	 elsif($label eq 'stop_codon'){
	    $phase = $stop_phase;
	    $stop_phase = ($stop_phase - $len) % 3;
	 }
	 
	 #fix phase for weird single base exons (phase == 2 && length == 1 ?)
	 $phase = 1 if($phase ne '.' && $phase > $len);
	 
	 my @data;
	 push(@data, $seq_id, $source, $label, $nB, $nE, $score, $strand, $phase);
	 my $nine = 'ID='.$e_id.';Parent='.$p;
	 $nine =~  s/\;$//;
	 push(@data, $nine);
	 
	 $c_l .= join("\t", @data)."\n";
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

   #initialize
   $cdss->{CDS} = [];
   $cdss->{five_prime_UTR} = [];
   $cdss->{three_prime_UTR} = [];
   $cdss->{start_codon} = [];
   $cdss->{stop_codon} = [];

   foreach my $hsp (@$hsps){
      my $nB = $hsp->nB('query');
      my $nE = $hsp->nE('query');
      
      my $q_strand = $phat_hit->strand('query');      
      my $hsp_end = abs($nE -$nB) + $hsp_start;
      my ($b, $e);

      #print "hsp_start:$hsp_start hsp_end:$hsp_end\n";
      #print "offset:$offset transl_end:$transl_end\n";

      #handle split on translation start
      if ($hsp_start <= $offset && $hsp_end > $offset){
	 my $d = $offset - $hsp_start;
	 $b = $q_strand == 1 ? $nB + ($d + 1) : $nB - ($d + 1);
	 my $ue = $q_strand == 1 ? $b-1 : $b+1;
	 my $ub = $nB;

	 my $utr = [$ub,
		    $ue,
		    $q_strand,
		    $phat_hit->name()];

	 push(@{$cdss->{five_prime_UTR}}, $utr);
	 push(@{$cdss->{t_ids}{five_prime_UTR}{$ub}{$ue}}, $id);
      }
      else {
	 $b = $nB;
      } 

      #handle split on translation end
      if ($hsp_start < $transl_end && $hsp_end >= $transl_end){
	 my $d = $hsp_end - $transl_end;
	 $e = $q_strand == 1 ? $nE - ($d + 1) : $nE + ($d + 1); 
	 my $ub = $q_strand == 1 ? $e+1 : $e-1; 
	 my $ue = $nE;

	 my $utr = [$ub,
		    $ue,
		    $q_strand,
		    $phat_hit->name()];

	 push(@{$cdss->{three_prime_UTR}}, $utr);
	 push(@{$cdss->{t_ids}{three_prime_UTR}{$ub}{$ue}}, $id);
      }
      else {
	 $e = $nE;
      }

      #error check CDS/UTR
      if(($b > $e && $hsp->strand('query') == 1) ||
	 ($b < $e && $hsp->strand('query') == -1)
	 ){
	 confess  "ERROR:B > E in GFFV3.pm strand:$q_strand\n";
      }

      #decide if CDS or UTR 
      my $label;
      if($hsp_end > $offset && $hsp_start < $transl_end){
	 $label = 'CDS';
      }
      elsif($hsp_start <= $offset && $hsp_end <= $offset){
	 $label = 'five_prime_UTR';
      }
      elsif($hsp_start >= $transl_end && $hsp_end >= $transl_end){
	 $label = 'three_prime_UTR';
      }
      else{
	 confess "ERROR: logic error in calculating CDS/UTR boundaries.\n";
      }

      #add CDS/UTR
      my $f = [$b,
	       $e,
	       $q_strand,
	       $phat_hit->name()];
      
      push(@{$cdss->{$label}}, $f);
      push(@{$cdss->{t_ids}{$label}{$b}{$e}},  $id);

      $hsp_start += abs($nE - $nB)+1;
   }

   #skip adding start and stop codons for now because Apollo croaks
   return; #temp? (not really, I'll probably keep it forever)

   #get start codon
   if($phat_hit->{_HAS_START}){
      my $start_c = 0;
      for (my $i = 0; $i < @{$cdss->{CDS}}; $i++){
	 last unless($start_c < 3);
	 
	 my $cds = $cdss->{CDS}->[$i];	 
	 my ($b, $e, $q_strand) = @$cds;
	 
	 my $diff = 2 - $start_c;
	 my $cb = $b;
	 my $ce;
	 if($q_strand == 1){
	    $ce = ($b+$diff <= $e) ? $b+$diff : $e;
	 }
	 else{
	    $ce = ($b-$diff >= $e) ? $b-$diff : $e;
	 }
	 
	 my $c = [$cb,
		  $ce,
		  $q_strand,
		  $phat_hit->name()];
	 $start_c += abs($ce - $cb) + 1;
	 
	 push(@{$cdss->{start_codon}}, $c);
	 push(@{$cdss->{t_ids}{start_codon}{$cb}{$ce}}, $id);
      }
   }

   #get stop codon
   if($phat_hit->{_HAS_STOP}){
      my $stop_c = 0;
      for (my $i = @{$cdss->{CDS}} - 1; $i >= 0; $i--){
	 last unless($phat_hit->{_HAS_STOP});
	 last unless($stop_c < 3);
	 
	 my $cds = $cdss->{CDS}->[$i];	 
	 my ($b, $e, $q_strand) = @$cds;

	 my $diff = 2 - $stop_c;
	 my $cb;
	 my $ce = $e;
	 if($q_strand == 1){
	    $cb = ($e-$diff >= $b) ? $e-$diff : $b;
	 }
	 else{
	    $cb = ($e+$diff <= $b) ? $e+$diff : $b;
	 }
	 
	 my $c = [$cb,
		  $ce,
		  $q_strand,
		  $phat_hit->name()];
	 $stop_c += abs($ce - $cb) + 1;
	 
	 unshift(@{$cdss->{stop_codon}}, $c);
	 push(@{$cdss->{t_ids}{stop_codon}{$cb}{$ce}},  $id);
      }
   }
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
	my $t_end   = $t->{t_end};
	my $t_name  = $t->{t_name};
	my $t_id    = $t->{t_id};
	my $t_qi    = $t->{t_qi};
	my $AED     = $t->{AED};
	my $eAED    = $t->{eAED};
	my $score   = '.'; #transcript scores causes chado errors/warning

	if($t->{hit}->algorithm =~ /est2genome|cdna2genome|protein2genome|BLASTX|rapsearch|est_gff|prot_gff/){
	   $score = $t->{hit}->score();
	}

	$t_name = uri_escape($t_name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standard
	$t_id = uri_escape($t_id, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
	
	#format informative name for GFF3
	$AED = sprintf '%.2f', $AED; # two decimal places
	$eAED = sprintf '%.2f', $eAED; # two decimal places

	my $t_s = $t_hit->strand('query') == 1 ? '+' : '-';

	my ($t_b, $t_e) = PhatHit_utils::get_span_of_hit($t_hit, 'query');

	($t_b, $t_e) = ($t_e, $t_b) if $t_b > $t_e;

	my @data;
	my $type = 'mRNA';
	$type = 'tRNA' if($t->{hit}->algorithm =~ /trnascan/i);
	$type = 'snoRNA' if($t->{hit}->algorithm =~ /snoscan/i);
	push(@data, $seq_id, 'maker', $type, $t_b, $t_e, $score, $t_s, '.');
	my $nine = 'ID='.$t_id.';Parent='.$g_id.';Name='.$t_name;
	   $nine .= ';_AED='.$AED if(defined($AED));
	   $nine .= ';_eAED='.$eAED if(defined($eAED));
	   $nine .= ';_QI='.$t_qi if(defined($t_qi));
	   $nine .= ';'.$t_hit->{-attrib} if($t_hit->{-attrib});
	if($t->{hit}->{_Alias}){
	    if($nine =~ /Alias\=([^\;\n]+)/){
		my @keepers = (@{[split(',', $1)]}, @{$t->{hit}->{_Alias}});
		my %uniq;
		@uniq{@keepers} =();
		my $alias = join(',', keys %uniq);
		$nine =~ s/\;*Alias\=[^\;\n]+/;Alias=$alias/;
	    }
	    else{
                my %uniq;
		@uniq{@{$t->{hit}->{_Alias}}} =();
                my $alias = join(',', keys %uniq);
                $nine .= ";Alias=$alias";
	    }
	}
	$nine =~  s/\;$//;

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
	my $hit      = shift; #incase I need anything else from the hit itself

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
	$hsp_name = uri_escape($hsp_name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards

	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id;
	   $nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE;
	   $nine .= ' '.$t_strand if($hsp->strand('hit'));
	   $nine .= ';Length='.$hit->length;
	   $nine .= ';Gap='.join(' ', @gap) if(@gap);
	   $nine .= ';'.$hsp->{-attrib} if($hsp->{-attrib});
	   $nine =~  s/\;$//;
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
	$class = "repeatrunner" if ($class eq 'blastx' || $class eq 'rapsearch');
	
	my $hsp_name = $hit_n;
	$hsp_name = uri_escape($hsp_name, '^a-zA-Z0-9\.\:\^\*\$\@\!\+\_\?\-\|'); #per gff standards
 
        my @data;
        push(@data, $seq_id, $class, $type, $nB, $nE, $score, $hsp_str, '.');
	my $nine  = 'ID='.$hsp_id.';Parent='.$hit_id;
	   $nine .= ';Target='.$hsp_name.' '.$tB.' '.$tE.' '.$t_strand;
	   $nine .= ';'.$hsp->{-attrib} if($hsp->{-attrib});
	   $nine =~  s/\;$//;
        push(@data, $nine);

        my $l = join("\t", @data);
        return $l;
}
#------------------------------------------------------------------------
1;


