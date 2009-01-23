#-------------------------------------------------------------------------------
#------                          CGL::Parser::gff3_maker                 -------
#-------------------------------------------------------------------------------
package evaluator::gff3_to_phatHit::Maker;
use strict;
use lib '/home/hao/Eval_Temp/maker/lib';
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
use Bio::DB::GFF;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use Bio::FeatureIO::gff;
use Bio::Tools::CodonTable;
use CGL::Annotation;
use Iterator::Fasta;
use Datastore::MD5;
use gff3::PhatHit;
use gff3::PhatHsp;
use File::Temp qw(tempfile);

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub new {
        my $class      = shift; #evaluator::gff3_to_phatHit::Maker
	my $gff_file   = shift; #gff file, can be for whole genome
	my $fasta_file = shift; #multi-fasta file for gff3
	my $ids        = shift; #id of contig to get genes for

	#get embedded gff3 fasta if no fasta provided
	if(! defined $fasta_file){
	   open(my $IN, "<", $gff_file);
	   
	   while (defined (my $line = <$IN>)){
	      last if ($line =~ /^\#\#FASTA/);
	   }
	   
	   if ($line =~ /^\#\#FASTA/){
	      my ($fh, $name) = tempfile;
	      $fasta_file = $name;
	      while (defined (my $line = <$IN>)){
		 chomp $line;
		 print $fh $line."\n" unless($line eq ''); #print non-empty lines
	      }
	      close($fh);
	      close($IN);
	   }
	   else{
	      die "ERROR: No fasta file provided for annotations\n";
	   }
	}

	#load a gene structure
	my ($genes, $seq, $seq_id) = 
	get_genes($gff_file, $fasta_file, 'exon', 0, $ids);

	#build phathits for gene structure
	my @phat_hits;
	foreach my $g (@{$genes}){
        	push(@phat_hits, load_hits($g, $seq));
	}

	return \@phat_hits;
}
#-------------------------------------------------------------------------------
sub glean {

        my $gff3_file  = shift;

        my $fh = new FileHandle();
           $fh->open($gff3_file);

        my $i = 0;
        while (defined(my $line = <$fh>)){

		next unless $line =~ /\S+/;

                if ($i < 2){
                        print $line;
                        $i++;
                        next;
                }

                chomp($line);

                my @stuff = split("\t", $line);

		if ($line eq '##FASTA'){
			print $line."\n";

			while (my $fasta = <$fh>){
				print $fasta;
			}
		}
		elsif (wanted($stuff[1], $stuff[2])){
			print $line."\n"; 
		}

                $i++;
        }

        $fh->close();

}
#-----------------------------------------------------------------------------
sub wanted {
	my $s = shift;
        my $x = shift;

	return 0 unless $s eq 'maker';

        return 1 if $x eq 'gene';
        return 1 if $x eq 'mRNA';
        return 1 if $x eq 'exon';
        return 1 if $x eq 'CDS';
        return 1 if $x eq 'scaffold';
        return 0;
}
#-------------------------------------------------------------------------------
sub to_gff3_contig {
	my $seg_id = shift;
	my $length = shift;

	my @fields;
	push(@fields, $seg_id, '.', 'contig', 1, $length);
	push(@fields, '.', '.', '.', "ID=$seg_id;Name=$seg_id");

 
	return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_gene {
	my $seg_id = shift;
	my $g      = shift;


	my $g_id   = get_annotation($g->{f}, 'ID');
	my $g_n    = get_annotation($g->{f}, 'Name');
	my $strand = $g->{f}->strand() == 1 ? '+' : '-';

	my @fields;

        push(@fields, $seg_id, 'maker', 'gene', $g->{i_start});
        push(@fields, $g->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$g_id;Name=$g_n");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_transcript {
	my $seg_id = shift;
        my $g      = shift;
        my $t      = shift;

        my $g_id   = get_annotation($g->{f}, 'ID');
	my $t_id   = get_annotation($t->{f}, 'ID');
	
	my $t_n    = get_annotation($t->{f}, 'Name');

        my $strand = $t->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'maker', 'mRNA', $t->{i_start});
        push(@fields, $t->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$t_id;Parent=$g_id;Name=$t_n");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_exon {
        my $seg_id = shift;
        my $t      = shift;
        my $e      = shift;

        my $e_id   = get_annotation($e->{f}, 'ID');
        my $t_id   = get_annotation($t->{f}, 'ID');


        my $strand = $e->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'maker', 'exon', $e->{i_start});
        push(@fields, $e->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$e_id;Parent=$t_id");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_cds {
        my $seg_id = shift;
        my $t      = shift;
        my $c      = shift;

        my $c_id   = get_annotation($c->{f}, 'ID');
        my $t_id   = get_annotation($t->{f}, 'ID');

        my $strand = $c->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'maker', 'CDS', $c->{i_start});
        push(@fields, $c->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$c_id;Parent=$t_id");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub get_annotation {
	my ($f,$type)=@_;;
	my $annotation_collection = $f->annotation;
	my ($annotation) = $annotation_collection->get_Annotations($type);
	
	return $annotation->value();
}

#-------------------------------------------------------------------------------
sub to_fasta_seq {                      
        my $id  = shift;
        my $seq = shift;
                
        my $fasta = Fasta::toFasta('>'.$id, \$seq);

        return $$fasta;
}  

#-------------------------------------------------------------------------------
sub to_gff3_seq {
	my $id  = shift;
	my $seq = shift;

	my $fasta = Fasta::toFasta('>'.$id, \$seq);

	return $$fasta;
}
#-------------------------------------------------------------------------------
sub split_file {
	my $gff3_file = shift;
	my $fasta_file= shift;
	my $ds_root   = shift;
		

	my ($genes, $seq, $seq_id) = 
	get_genes($gff3_file,$fasta_file, 'exon', 500);
	
	my $ds;
	if ($ds_root) {
		$ds = Datastore::MD5->new("root"=>$ds_root, "depth"=>2);
	}
	
	foreach my $g (@{$genes}){
		my $g_id = get_annotation($g->{f}, 'ID');

                my $file_base = "${seq_id}_${g_id}";
                my $gff3_file;
                if ($ds_root) {
                        $ds->mkdir($file_base) || die "Unable to ds->mkdir() for $file_base";
                        $gff3_file = sprintf("%s/%s.%s",
                                        $ds->id_to_dir($file_base),
                                        $file_base,
                                        'gff3');
                }
                else {
                        $gff3_file = sprintf("%s/%s.%s",
                                        q{},
                                        $file_base,
                                        'gff3');
                }
		my $fh = new FileHandle();
		   $fh->open(">$gff3_file");

		print $fh "##gff-version   3\n";

		my $seg_id;
		if ($g->{f}->strand() == 1){
			$seg_id = $seq_id.':'.$g->{src_s}.':'.$g->{src_e};
		}
		else {
			$seg_id = $seq_id.':'.$g->{src_e}.':'.$g->{src_s};
		}

		my $to_gff3_contig_var= to_gff3_contig($seg_id, length($g->{seq}))."\n";
		print $fh $to_gff3_contig_var;
		my $to_gff3_gene_var= to_gff3_gene($seg_id, $g)."\n";
		print $fh $to_gff3_gene_var;

		foreach my $t (@{$g->{mRNAs}}){
			print $fh to_gff3_transcript($seg_id, $g, $t)."\n";
			foreach my $e (@{$t->{exons}}){
				print $fh to_gff3_exon($seg_id, $t, $e)."\n";
			}
                        foreach my $c (@{$t->{cdss}}){
                                print $fh to_gff3_cds($seg_id, $t, $c)."\n";
                        }

		}
		$fh->close();
                my $fasta_file;
                if ($ds_root) {
                        $fasta_file = sprintf("%s/%s.%s",
                                              $ds->id_to_dir($file_base),
                                              $file_base,
                                              'fasta');
                }
                else {
                        $fasta_file = sprintf("%s/%s.%s",
                                              q{},
                                              $file_base,
                                              'fasta');
                }


		$fh->open(">$fasta_file");

		print $fh to_fasta_seq($seg_id, $g->{seq});
		$fh->close();

	}

}
#-------------------------------------------------------------------------------
sub base_to_space {

	my $nbeg = shift;
	my $nend = shift;
    
	if($nbeg < $nend || $nbeg == $nend){
        	$nbeg--;
    	}
	elsif ($nbeg >  $nend){
		$nend--;
	}
	return ($nbeg, $nend);
}
#----------------------------------------------------------------------

sub get_features {
	my $file = shift;
	my $ids  = shift;

	my $parser = new Bio::Tools::GFF(-gff_version => 3,
		     			 -file        => $file,
		                         );

	my @features;
	while (my $f = $parser->next_feature()) {
		push(@features, $f) if($f->seq_id eq $ids);
	}
	return \@features;
}
#-------------------------------------------------------------------------------
sub load_annotations {
	my $contig    = shift;
	my $cgl_genes = shift;

	my $cgl = {};

	push(@{$cgl->{contigs}} , $contig);
	push(@{$cgl->{genes}}, @{$cgl_genes});

	bless $cgl, 'CGL::Annotation';

	return $cgl;
}
#-------------------------------------------------------------------------------
sub load_contig {
	my $id  = shift;
	my $seq = shift;

	my $contig = {};

	$contig->{feature_id} = $id;
	$contig->{id}         = $id;
	$contig->{name}       = $id;
	$contig->{residues}   = $$seq;
	$contig->{type}       = 'contig';
	$contig->{uniquename} = $id;

	push(@{$contig->{locations}}, 
	load_feature_location(1, length($$seq), $id));

	bless $contig, 'CGL::Annotation::Feature::Contig';

	return $contig;
}
#-------------------------------------------------------------------------------
sub get_t_offset_and_end {
        my $t  = shift;
                
        my $t_seq = $t->{seq};
        my $c_seq = $t->{cds_seq};
                
        my $t_offset = index($t_seq, $c_seq);
                
        my $t_end = $t_offset + length($c_seq) -1;
                
        die "problem in gff3_to_phatHit::get_t_offset\n"
        if $t_offset == -1;

        return ($t_offset, $t_end);
}     
#-------------------------------------------------------------------------------

sub load_hits {
        my $g   = shift;
        my $seq = shift;


	if $g->{f}->primary_tag();
        
        my $gene_id   = $g->{f}->get_Annotations('ID')->value();
        my $gene_name = $g->{f}->get_Annotations('Name')->value();
        
        
        my @phat_hits;
        foreach my $t (@{$g->{mRNAs}}){

                my $tran_id   = $t->{f}->get_Annotations('ID')->value();
                my $tran_name = $t->{f}->get_Annotations('Name')->value();

                my $description =
                "g_name=$gene_name;g_id=$gene_id;t_name=$tran_name;t_id=$tran_id";

                my $f = new gff3::PhatHit('-name'  => $tran_name,
                                   '-description'  => $description,
                                   '-algorithm'    => 'from gff3',
                                   '-length'       => length($t->{seq}),
                                   '-score'        => $t->{score},
                                 );


                $f->{gene_id} = $gene_id;
                $f->{gene_name} = $gene_name;
		$f->{_tran_id}   = $tran_id;
		$f->{maker_qi}  = $t->{maker_qi};

                my $type = $t->{f}->primary_tag;
		$f->{transcript_type}=$type;

                $f->queryLength(length($t->{seq}));

                my ($t_offset, $t_end) = get_t_offset_and_end($t);

                $f->{translation_offset} = $t_offset;
                $f->{translation_end}    = $t_end;
                $f->{seq} = $t->{seq};

                my $hsps = load_hsps($t, $seq);
                my $cdss = load_cdss($t, $seq);

                foreach my $hsp (@{$hsps}){
                        $f->add_hsp($hsp);
                }

                $f->{cdss} = $cdss;

                push(@phat_hits, $f);
        }

        return \@phat_hits;
}
#-------------------------------------------------------------------------------

sub load_cdss {
        my $t   = shift;
        my $seq = shift;

        my @hsps;
        my $hit_start = 1;

        foreach my $e (@{$t->{cdss}}){
                my @args;
                my $hit_end = $e->{f}->end - $e->{f}->start + $hit_start;

                push(@args, '-query_start');
                push(@args, $e->{f}->start);

                push(@args, '-query_seq');
                push(@args, $e->{seq});

                push(@args, '-score');
                push(@args, 100000000000000000);

                push(@args, '-homology_seq');
                push(@args, $e->{seq});

                push(@args, '-hit_start');
                push(@args, $hit_start);

                push(@args, '-hit_seq');
                push(@args, $e->{seq});

                push(@args, '-hsp_length');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-identical');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-hit_length');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-query_name');
                push(@args, $e->{f}->get_Annotations('ID')->value());

                push(@args, '-algorithm');
                push(@args, 'from gff3');

                push(@args, '-bits');
                push(@args, 2*($e->{f}->end - $e->{f}->start + 1));

                push(@args, '-evalue');
                push(@args, 0.0);

                push(@args, '-pvalue');
                push(@args, 0.0);
                push(@args, '-query_length');
                push(@args, length($e->{seq}));

                push(@args, '-query_end');
                push(@args, $e->{f}->end);

                push(@args, '-conserved');
                push(@args, length($e->{seq}));

                push(@args, '-hit_name');
                push(@args, $e->{f}->get_Annotations('ID')->value());

                push(@args, '-hit_end');
                push(@args, $hit_end);

                push(@args, '-query_gaps');
                push(@args, 0);

                push(@args, '-hit_gaps');
                push(@args, 0);

                my $hsp = new gff3::PhatHsp(@args);
                   $hsp->queryName($e->{f}->get_Annotations('ID')->value());

                #-------------------------------------------------
                # setting strand because bioperl is all messed up!
                #------------------------------------------------
                if ($e->{f}->strand == 1 ){
                        $hsp->{_strand_hack}->{query} = 1;
                        $hsp->{_strand_hack}->{hit}   = 1;
                }
                else {
                        $hsp->{_strand_hack}->{query} = -1;
                        $hsp->{_strand_hack}->{hit}   =  1;
               }

               $hit_start += $e->{f}->end - $e->{f}->start + 1;


                push(@hsps, $hsp);
        }

        return \@hsps;
}
#-------------------------------------------------------------------------------
sub load_hsps {
        my $t   = shift;
        my $seq = shift;

        my @hsps;
        my $hit_start = 1;

        foreach my $e (@{$t->{exons}}){
                my @args;
                my $hit_end = $e->{f}->end - $e->{f}->start + $hit_start;

                push(@args, '-query_start');
                push(@args, $e->{f}->start);

                push(@args, '-query_seq');
                push(@args, $e->{seq});

                push(@args, '-score');
                push(@args, 100000000000000000);

                push(@args, '-homology_seq');
                push(@args, $e->{seq});

                push(@args, '-hit_start');
                push(@args, $hit_start);

                push(@args, '-hit_seq');
                push(@args, $e->{seq});

                push(@args, '-hsp_length');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-identical');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-hit_length');
                push(@args, $e->{f}->end - $e->{f}->start + 1);

                push(@args, '-query_name');
                push(@args, $e->{f}->get_Annotations('ID')->value());

                push(@args, '-algorithm');
                push(@args, 'from gff3');

                push(@args, '-bits');
                push(@args, 2*($e->{f}->end - $e->{f}->start + 1));

                push(@args, '-evalue');
                push(@args, 0.0);

                push(@args, '-pvalue');
                push(@args, 0.0);
                push(@args, '-query_length');
                push(@args, length($e->{seq}));

                push(@args, '-query_end');
                push(@args, $e->{f}->end);

                push(@args, '-conserved');
                push(@args, length($e->{seq}));

                push(@args, '-hit_name');
                push(@args, $e->{f}->get_Annotations('ID')->value());

                push(@args, '-hit_end');
                push(@args, $hit_end);

                push(@args, '-query_gaps');
                push(@args, 0);

                push(@args, '-hit_gaps');
                push(@args, 0);

                my $hsp = new gff3::PhatHsp(@args);
                   $hsp->queryName($e->{f}->get_Annotations('ID')->value());

                #-------------------------------------------------
                # setting strand because bioperl is all messed up!
                #------------------------------------------------
                if ($e->{f}->strand == 1 ){
                        $hsp->{_strand_hack}->{query} = 1;
                        $hsp->{_strand_hack}->{hit}   = 1;
                }
                else {
                        $hsp->{_strand_hack}->{query} = -1;
                        $hsp->{_strand_hack}->{hit}   =  1;
               }

               $hit_start += $e->{f}->end - $e->{f}->start + 1;


                push(@hsps, $hsp);
        }

        return \@hsps;
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
sub load_feature_location {
	my $nbeg = shift;
	my $nend = shift;
	my $id   = shift;

        my $l = {};

	($nbeg, $nend) = base_to_space($nbeg, $nend);

        $l->{nbeg} = $nbeg;
        $l->{nend} = $nend;

        $l->{srcfeature_id} = $id;

        $l->{text} = 'converted from maker gff3';

        bless $l, 'CGL::Annotation::FeatureLocation';

        return $l;
}
#-------------------------------------------------------------------------------
sub get_l_beg_end {

        my $feature = shift;

        my $start  = $feature->{f} ? $feature->{f}->start()  : undef;
        my $end    = $feature->{f} ? $feature->{f}->end()    : undef;
        my $strand = $feature->{f} ? $feature->{f}->strand() : undef;

        if ($strand == 1){
                return ($start, $end);
        }
        elsif ($strand == -1){
                return ($end, $start);
        }
        else {
                die "dead in WormBase::get_l_beg_end\n";
        }
}
#-------------------------------------------------------------------------------
sub load_translations {
	my $t = shift;

	my $f_id = 'protein:'.get_annotation($t->{f},'ID');

	my $transl = {};

	$transl->{feature_id} = $f_id; 

	$transl->{id} = $f_id;
 
        my $sorted = sort_cdss($t);

        my $alpha =  $sorted->[0];
        my $omega =  $sorted->[-1];

        my $transl_offset = get_translation_offset($t);
        my $trans_end     = $transl_offset + length($t->{cds_seq}) - 3;
        
        my $hell = length($t->{seq});
        
        
        my ($a_l_beg, $a_l_end) = get_l_beg_end($alpha);
        my ($o_l_beg, $o_l_end) = get_l_beg_end($omega);

=cut;   
        my $nbeg = $sorted->[0]->{f}->start();
        my $nend = $sorted->[0]->{f}->end();

        ($nbeg, $nend) = ($nend, $nbeg) if $sorted->[0]->{f}->strand() == -1;

=cut;
        my $src_f_id = $t->{src_f_id};

        push(@{$transl->{locations}},
        load_feature_location($a_l_beg, $o_l_end, $src_f_id));

	my $trn = new Bio::Tools::CodonTable();

	my $p_seq = $trn->translate($t->{cds_seq});

	$p_seq =~ s/\*$//;

	$transl->{name} = $f_id;

	my $oF_id = get_annotation($t->{f},'ID');
	my $sF_id = $transl->{id};

	push(@{$transl->{relationships}},
	load_feature_relationship( $oF_id, $sF_id, 'produced_by'));

	my %props;
	$props{codon_start} = get_translation_offset($t);
	$props{gene}        = $t->{part_of}->[0];
	$props{product}     = $transl->{id};
	$props{protein_id}  = $transl->{id};

	$transl->{properties} = \%props;

	$transl->{residues} = $p_seq;

	$transl->{type} = 'protein';
	$transl->{uniquename} = $transl->{id};
	
	bless $transl, 'CGL::Annotation::Feature::Protein';

	return [$transl];

}
#-------------------------------------------------------------------------------
sub get_translation_offset {
        my $t = shift;

        my $sorted_exons = sort_exons($t);
        my $sorted_cdss  = sort_cdss($t);

        my $e_0_b =
            $sorted_exons->[0]->{f} ?
            $sorted_exons->[0]->{f}->start() :
            undef;
        my $c_0_b =
            $sorted_cdss->[0]->{f}  ?
            $sorted_cdss->[0]->{f}->start() :
            undef;

        return abs($e_0_b - $c_0_b) + 1;

}
#-------------------------------------------------------------------------------
sub load_exons {
	my $t       = shift;
	my $transcr = shift;

        my $sorted_exons = sort_exons($t);

        foreach my $e (@{$sorted_exons}){
                my $exon = load_exon($e);
                push(@{$transcr->{_sorted_exons}}, $exon);
                push(@{$transcr->{exons}}, $exon);
        }

} 
#-------------------------------------------------------------------------------
sub load_introns {
	my $t       = shift;
	my $transcr = shift;
	my $g_seq   = shift;

	my $sorted_exons = sort_exons($t);

	my $num_exons = @{$sorted_exons};

	return [] if $num_exons == 1;

	my $j = 0;
	for (my $i = 1; $i < @{$sorted_exons}; $i++){
		my $e_u = $sorted_exons->[$i - 1];
		my $e_d = $sorted_exons->[$i];

		my $intron = load_intron($e_u, $e_d, $g_seq, $t, $j);

		push (@{$transcr->{introns}}, $intron);
		$j++;
	}
}
#-------------------------------------------------------------------------------
sub load_intron {
	my $e_u = shift;
	my $e_d = shift;
	my $seq = shift;
	my $t   = shift;
	my $i   = shift;

	my $intron = {};

        my $int_b = $e_u->{f}->strand() == 1
                    ? $e_u->{f}->end() + 1
                    : $e_u->{f}->start();

        my $int_e = $e_d->{f}->strand() == 1
                    ? $e_d->{f}->start() - 1
                    : $e_d->{f}->end();


        my $length = $e_u->{f}->strand() == 1 ? abs($int_e - $int_b) + 1
	                                      : abs($int_e - $int_b) - 1;

	my $offset = $e_u->{f}->strand() == 1 ? $int_b - 1 : $int_e;
        my $int_seq = substr($$seq, $offset, $length);

        $int_seq = Fasta::revComp($int_seq)
        if $e_u->{f}->strand() == -1;

	my $s = $e_u->{f}->strand();

	$intron->{id} = "intron:$int_b:$int_e";
	$intron->{inScope} = 1;
	
        foreach my $p_id (@{$e_u->{part_of}}){
		push(@{$intron->{relationships}},
		load_feature_relationship($p_id, $intron->{id}, 'part_of'));
        }
		
	my $src_f_id = $t->{src_f_id};

	push(@{$intron->{locations}}, 
	load_feature_location($int_b, $int_e, $src_f_id));

	$intron->{residues} = $int_seq;

	$intron->{seqlen} =  length($int_seq);

	$intron->{src_id} = $t->{part_of}->[0].":$int_b:$int_e";

	$intron->{strand} = $e_u->{f}->strand();	

	$intron->{type} = 'intron';

	$intron->{name} = "intron:$i";
	
	bless $intron, 'CGL::Annotation::Feature::Intron';

	return $intron;
}
#-------------------------------------------------------------------------------
sub load_exon {
	my $e = shift;

	my $exon ={};

	$exon->{feature_id} = get_annotation($e->{f},'ID');
	$exon->{id}         = get_annotation($e->{f},'ID');
	$exon->{inScope}    = 1;


        my $nbeg = $e->{f}->start();
        my $nend = $e->{f}->end();

	my $scr_f_id = $e->{src_f_id};

        ($nbeg, $nend) = ($nend, $nbeg)
        if $e->{f}->strand() == -1;

	push(@{$exon->{locations}}, 
	load_feature_location($nbeg, $nend, $scr_f_id));


	$exon->{name} = get_annotation($e->{f},'ID');


	my $sF_id = get_annotation($e->{f},'ID');

	foreach my $oF_id (@{$e->{part_of}}){
		push(@{$exon->{relationships}},
		load_feature_relationship($oF_id, $sF_id, 'part_of')); 
	}	
	$exon->{residues} = $e->{seq};

	
	$exon->{src_id} = $e->{src_f_id}.':'.$nbeg.':'.$nend;

	$exon->{type} = 'exon';

	$exon->{uniquename} = 
	$e->{src_f_id}.':'.get_annotation($e->{f},'ID');

	bless $exon, 'CGL::Annotation::Feature::Exon';

	return $exon;
}
#-------------------------------------------------------------------------------
sub load_feature_relationship {
	my $oF_id  = shift;
	my $sF_id  = shift;
	my $logus  = shift;

        my $r = {};
           $r->{logus} = $logus;
           $r->{oF}    = $oF_id;
           $r->{sF}    = $sF_id;

        bless $r, 'CGL::Annotation::FeatureRelationship';
	
	return $r;
}
#-------------------------------------------------------------------------------
#takes multi-fasta file name and id of single contig to load
#returns def and reference to seq for that contig
sub load_seq {
        my $file = shift;
	my $ids = shift;

	my $index = new Bio::DB::Fasta($file);
	my $fastaObj = $index->get_Seq_by_id($ids);

	die "ERROR: ID not found in fasta file\n" if (! defined($fastaObj));
	
	my $seq      = \($fastaObj->seq());
	my $def      = $index->header($ids);

        return ($def, $seq);
}
#-------------------------------------------------------------------------------
sub grab {
        #Barry
        my $types    = shift;
        my $source   = shift;
        my $features = shift;
        my $c_id     = shift;

        my %booty;
        my $i = 0;

        for my $type (@{$types}) {
                foreach my $f (@{$features}){
                        my $tag_t = $f->primary_tag();
                        my $tag_s = $f->source_tag();

                        if ($tag_t eq $type && $tag_s eq  $source) {
                                my $id =get_annotation($f,'ID');
                                my $p_ids = get_p_ids($f);
                                my $score = $f->score();

				if(@{$p_ids}){
				   foreach my $p_id (@{$p_ids}){
				      
				      push(@{$booty{$p_id}}, {f        => $f,
							      id       => $id,
							      src_f_id => $c_id,
							      part_of  => $p_ids,
							      score    => $score,
							     });
				   }
				   $i++;
				}
				else{
				   push(@{$booty{'contig'}}, {f        => $f,
							      id       => $id,
							      src_f_id => $c_id,
							      part_of  => $p_ids,
							      score    => $score,
							     });
				}
                        }
                }
        }
        return \%booty;
}
#-------------------------------------------------------------------------------
sub get_p_ids {
	my $f = shift;

	my @parents = $f->annotation->get_Annotations('Parent');
        my @p_ids;

       foreach my $p (@parents){
       		push(@p_ids, $p->{value});
       }

	return \@p_ids;
}
#-------------------------------------------------------------------------------
sub get_evidence {
	my $gff3_file  = shift;
	my $fasta_file = shift;
	my $source     = shift;
	my $flank      = shift;
	my $ids	       = shift;

	my ($def, $seq)   = load_seq($fasta_file, $ids);

	my ($c_id) = $def =~ />(\S+)/;

	my $features = get_features($gff3_file, $ids);

	my $hsp_bases = grab(['match_part'], [$source], $features, $c_id);
	my @types = ('match',
		     'nucleotide_match',
		     'expressed_sequence_match',
		     'translated_nucleotide_match',
		     'EST_match',
		     'cDNA_match',
		     'cross_genome_match',
		     'protein_match'
		     
		    );
	my $hit_bases  = grab(\@types, [$source], $features, $c_id);

        my @genes;
        foreach my $f (@{$features}){
                my $tag = $f->primary_tag();
                my $source = $f->source_tag();

                if ($tag eq 'gene' && $source eq 'maker') {
			
		my $id=get_annotation($f,'ID');
	
		    push(@genes, {'f'        => $f,
				  'mRNAs'    => $mRNAs->{$id},
				  'src_f_id' => $c_id,
				  'part_of'  => [],
				  'id'       => $id,
			      });
                }
        }
        print STDERR "validating genes\n";
        my @valid_genes;
        for my $gene (@genes) {
                push @valid_genes, $gene if validate_gene($gene);
        }
        print STDERR "...finished\n";

        print STDERR " loading seqs\n";
        load_seqs(\@valid_genes, $seq, $flank);
        print STDERR "...finished\n";

        return (\@valid_genes, $seq, $c_id);

}
#-------------------------------------------------------------------------------
sub get_genes {
	my $gff3_file      = shift;
	my $fasta_file     = shift;
	my $base_type = shift;
	my $flank     = shift;
	my $ids	      = shift;

	my ($def, $seq)   = load_seq($fasta_file, $ids);

	my ($c_id) = $def =~ />(\S+)/;

	my $features = get_features($gff3_file, $ids);

	my $exons = grab([$base_type], 'maker', $features, $c_id);
	my $cdss  = grab(['CDS'], 'maker', $features, $c_id);
	my $mRNAs = grab(['mRNA'], 'maker', $features, $c_id);

	foreach my $p_id (keys %{$mRNAs}){
		for (my $i = 0; $ i < @{$mRNAs->{$p_id}}; $i++) {
			my $f  = $mRNAs->{$p_id}->[$i]->{f};
			my $id=get_annotation($f,'ID'); 

			$mRNAs->{$p_id}->[$i]->{exons} = $exons->{$id};
			$mRNAs->{$p_id}->[$i]->{cdss}  = $cdss->{$id};
			$mRNAs->{$p_id}->[$i]->{maker_qi} =
				get_maker_qi($mRNAs->{$p_id}->[$i]);
		}
	}
        my @genes;
        foreach my $f (@{$features}){
                my $tag = $f->primary_tag();
                my $source = $f->source_tag();

                if ($tag eq 'gene' && $source eq 'maker') {
			
		my $id=get_annotation($f,'ID');
	
		    push(@genes, {'f'        => $f,
				  'mRNAs'    => $mRNAs->{$id},
				  'src_f_id' => $c_id,
				  'part_of'  => [],
				  'id'       => $id,
			      });
                }
        }
        print STDERR "validating genes\n";
        my @valid_genes;
        for my $gene (@genes) {
                push @valid_genes, $gene if validate_gene($gene);
        }
        print STDERR "...finished\n";

        print STDERR " loading seqs\n";
        load_seqs(\@valid_genes, $seq, $flank);
        print STDERR "...finished\n";

        return (\@valid_genes, $seq, $c_id);

}
#-------------------------------------------------------------------------------
sub validate_gene {
        my $gene = shift;
                
        my @strands;
        
        #Check strand and fail if we don't get a valid strand value.
        my $g_strand = $gene->{f}{_location}{_strand};
        if ($g_strand != 1 && $g_strand != -1) {
                warn "Invalid strand in gene caught at " .
                    "FlyBase::validate_gene\n";
                return undef;
        }
        
        #Push the strand onto an array so that we can check all strands for
        #internal consistancy at the end.
        push @strands, $g_strand;
        
        for my $transcript (@{$gene->{transcripts}}) {
                #Check strand and fail if we don't get a valid strand value.
                my $t_strand = $transcript->{f}{_location}{_strand};
                if ($t_strand != 1 && $t_strand != -1) {
                        warn "Invalid strand in transcript caught " .
                            "at FlyBase::validate_gene\n";
                        return undef;
                }
                #Push the strand onto an array so that we can check 
                #all strands for internal consistancy at the end.
                push @strands, $t_strand;
                
                for my $cds (@{$transcript->{cdss}}) {
                        #Check strand and fail if we don't get a valid strand value.
                        my $c_strand = $transcript->{f}{_location}{_strand};
                        if ($c_strand != 1 && $c_strand != -1) {
                                warn "Invalid strand in cds caught " .
                                    "at FlyBase::validate_gene\n";
                                return undef;
                        }
                        #Push the strand onto an array so that we can check
                        #all strands for internal consistancy at the end.
                        push @strands, $c_strand;
                }

                for my $exon (@{$transcript->{exons}}) {
                        #Check strand and fail if we don't get a valid strand value.
                        my $e_strand = $transcript->{f}{_location}{_strand};
                        if ($e_strand != 1 && $e_strand != -1) {
                                warn "Invalid strand in exon caught " .
                                    "at FlyBase::validate_gene\n";
                                return undef;
                        }
                        #Push the strand onto an array so that we can check
                        #all strands for internal consistancy at the end.
                        push @strands, $e_strand;
                }
        }
        #Check that all strands are the same, and fail if they are not.
        my $first = shift @strands;
        return undef if grep {$first != $_} @strands;

        #No failures, so return success.
        return 1;
}
#-------------------------------------------------------------------------------
sub get_exon_seq {
	my $e   = shift;
	my $seq = shift;

        my $e_b = $e->{f}->start();
        my $e_e = $e->{f}->end();

        my $length = abs($e_e - $e_b) + 1;

        my $exon_seq = substr($$seq, $e_b - 1, $length);

        $exon_seq = Fasta::revComp($exon_seq)
        if $e->{f}->strand() == -1;


	return $exon_seq;
}
#-------------------------------------------------------------------------------
sub get_gene_seq {
        my $g   = shift;
        my $seq = shift;
	my $flank = shift || 0;

        my $g_b = $g->{f}->start();
        my $g_e = $g->{f}->end();

	my $start;
	my $length;
	my ($src_s, $src_e);
	if (($g_b - 1) - $flank < 0){
		$start = $g_b - 1;
		$length = abs($g_e - $g_b) + 1 + $flank;
		$src_s = 1;
		$src_e = $g_e + $flank > length($$seq) 
		? length($$seq) : $g_e + $flank; 
	}
	else {
		$start = ($g_b - 1) - $flank;
		$length = abs($g_e - $g_b) + 1 + $flank + $flank;
		$src_s = $g_b - $flank;
		$src_e = $g_e + $flank > length($$seq)
		? length($$seq) : $g_e + $flank; 
	} 

        my $g_seq = substr($$seq, $start, $length);

        return ($g_seq, $src_s, $src_e);
}
#-------------------------------------------------------------------------------
sub load_seqs {
	my $genes = shift;
	my $seq   = shift;
	my $flank = shift || 0;

	foreach my $g (@{$genes}){
		my ($g_seg_seq, $src_start, $src_end) = 
		    get_gene_seq($g, $seq, $flank);

		$g->{seq}     = $g_seg_seq;
		$g->{src_s}   = $src_start;
		$g->{src_e}   = $src_end;
		$g->{i_start} = 1;
		$g->{i_end}   = length($g_seg_seq);

        	foreach my $t (@{$g->{mRNAs}}){
                	foreach my $e (@{$t->{exons}}){
				my $e_seq = get_exon_seq($e, $seq); 
				
				$e->{seq} = $e_seq;
                		$e->{i_start} = 
				$e->{f}->start() - $g->{src_s} + 1;

				$e->{i_end}   =
                                $e->{f}->end() - $g->{src_s} + 1;  
                	}
                        foreach my $c (@{$t->{cdss}}){
                                my $c_seq = get_exon_seq($c, $seq);

				$c->{seq} = $c_seq;
                               	$c->{i_start} =
                               	$c->{f}->start() - $g->{src_s} + 1;

                               	$c->{i_end}   =
                               	$c->{f}->end() - $g->{src_s} + 1;
                        }

			my $t_seq   = get_mRNA_seq($t, $seq); 
			my $cds_seq = get_cds_seq($t, $seq);

                       $t->{seq} = $t_seq;
		       $t->{cds_seq} = $cds_seq;


              		$t->{i_start} = 
			$t->{f}->start() - $g->{src_s} + 1;
              		$t->{i_end}   = 
			$t->{f}->end()   - $g->{src_s} + 1; 

        	}
	}
}
#-------------------------------------------------------------------------------
sub get_cds_seq {
        my $t   = shift;
        my $seq = shift;

        my $sorted = sort_cdss($t);

        my $cds_seq;
        foreach my $c (@{$sorted}){
                my $c_seq = $c->{seq} || get_exon_seq($c, $seq);
                $cds_seq .= $c_seq;

        }
        return $cds_seq;
}
#-------------------------------------------------------------------------------
sub get_mRNA_seq {
	my $t   = shift;
	my $seq = shift;
	
	my $sorted = sort_exons($t);

	my $transcript;
	foreach my $e (@{$sorted}){
		my $exon_seq = $e->{seq} || get_exon_seq($e, $seq);
                $transcript .= $exon_seq;
	
	}
	return $transcript;
}
#-------------------------------------------------------------------------------
sub sort_exons {
	my $t = shift;

	my @sorted;
	if    ($t->{f}->strand() ==  1){
		@sorted = 
		sort {$a->{f}->start <=> $b->{f}->start} @{$t->{exons}}; 
	}
	elsif ($t->{f}->strand() == -1){
		@sorted = 
		sort {$b->{f}->start <=> $a->{f}->start} @{$t->{exons}}; 
	}
	else {
		die "unknown strand in GFF3::Maker::sort_exons!\n";
	}
	return \@sorted;
}
#-------------------------------------------------------------------------------
sub sort_cdss {
        my $t = shift;

        my @sorted;
        if    ($t->{f}->strand() ==  1){
                @sorted =
                sort {$a->{f}->start <=> $b->{f}->start} @{$t->{cdss}};
        }
        elsif ($t->{f}->strand() == -1){
                @sorted =
                sort {$b->{f}->start <=> $a->{f}->start} @{$t->{cdss}};
        }
        else {
                die "unknown strand in GFF3::Maker::sort_cdss!\n";
        }
        return \@sorted;
}
#-------------------------------------------------------------------------------
sub get_maker_qi {
	my $mRNA = shift;

	my $name = $mRNA->{f}->{annotation}->{'_annotation'}->{Name}->[0]->{value};

	return 'NA' unless defined $name;

	my ($maker_qi) = $name =~ /QI:([\d\|\.-]+)/;
	return $maker_qi;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "GFF::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#----------------------------------------------------------------------------

1;
