#-------------------------------------------------------------------------------
#------                     CGL::Annotation::GFF3::WormBase              -------
#-------------------------------------------------------------------------------
package CGL::Annotation::GFF3::WormBase;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
use Bio::DB::GFF;
use Bio::Tools::GFF;
use Bio::FeatureIO::gff;
use Bio::Tools::CodonTable;
use Bio::SeqIO;
use CGL::Annotation;

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub new {
        my $class      = shift;
        my $gff3_file  = shift;
        my $fasta_file = shift;

        my ($genes, $seq, $seq_id) =
        get_genes($gff3_file, $fasta_file, 'exon', 0);

        my @cgl_genes;
        foreach my $g (@{$genes}){
                push(@cgl_genes, load_gene($g, $seq));
        }

        my $contig = load_contig($seq_id, $seq);

        my $cgl    = load_annotations($contig, \@cgl_genes);


        return $cgl;
}
#-------------------------------------------------------------------------------
sub glean {

        my $gff3_file  = shift;

        my $fh = new FileHandle();
           $fh->open($gff3_file);
        
        my $i = 0;
        while (my $line = <$fh>){
        
                if ($i < 2){
                        print $line;
                        $i++;
                        next;
                }
        
                chomp($line);
        
                my @stuff = split("\t", $line);
        
                print $line."\n" if wanted($stuff[1], $stuff[2]);

                $i++;
        }

        $fh->close();

}
#-----------------------------------------------------------------------------
sub wanted {
	my $s = shift;
        my $x = shift;

	return 1 if 
	    $s eq 'Coding_transcript' &&
	    grep {$x eq $_} qw(gene mRNA exon CDS five_prime_UTR three_prime_UTR);

#        return 1 if $x eq 'gene' && $s eq 'Coding_transcript';
#        return 1 if $x eq 'mRNA' && $s eq 'Coding_transcript';
#        return 1 if $x eq 'exon' && $s eq 'Coding_transcript';
#        return 1 if $x eq 'CDS'  && $s eq 'Coding_transcript';
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

	my $g_id   = $g->{id};
	my $g_n    = $g->{f}->get_Annotations('Name')->value();
	my $strand = $g->{f}->strand() == 1 ? '+' : '-';

	my @fields;

        push(@fields, $seg_id, 'Coding_transcript', 'gene', $g->{i_start});
        push(@fields, $g->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$g_id;Name=$g_n");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_mRNA {
	my $seg_id = shift;
        my $g      = shift;
        my $t      = shift;

        my $g_id   = $g->{id};
	my $t_id   = $t->{id};
	my $t_n    = $t->{f}->get_Annotations('Name')->value();

        my $strand = $t->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'Coding_transcript', 'mRNA', $t->{i_start});
        push(@fields, $t->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$t_id;Parent=$g_id;Name=$t_n");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_exon {
        my $seg_id = shift;
        my $t      = shift;
        my $e      = shift;

        my $t_id   = $t->{id};
        my $e_id   = $e->{id};
	my $g_id   = $t->{f}->get_Annotations('Parent')->value();
        my $strand = $e->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'Coding_transcript', 'exon', $e->{i_start});
        push(@fields, $e->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$e_id;Parent=$g_id");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_gff3_cds {
        my $seg_id = shift;
        my $t      = shift;
        my $c      = shift;

        my $t_id   = $t->{id};
        my $c_id   = $c->{id};

        my $strand = $c->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'Coding_transcript', 'CDS', $c->{i_start});
        push(@fields, $c->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$c_id;Parent=$t_id");

        return join("\t", @fields);

}
#-------------------------------------------------------------------------------
sub to_fasta_seq {
	my $id  = shift;
	my $seq = shift;

	my $fasta = Fasta::toFasta('>'.$id, \$seq);

	return $$fasta;
}
#-------------------------------------------------------------------------------
sub split_file {
	my $gff3_file  = shift;
	my $fasta_file = shift;
	my $ds_root    = shift;

	my ($genes, $seq, $seq_id) = 
	get_genes($gff3_file, $fasta_file, 'exon', 500);

        my $ds;
        if ($ds_root) {
                $ds = Datastore::MD5->new("root" => $ds_root, "depth" => 2);
	}

	my $i = 0;
        foreach my $g (@{$genes}){
                my $g_id = $g->{id};

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

                print $fh to_gff3_contig($seg_id, length($g->{seq}))."\n";

                print $fh to_gff3_gene($seg_id, $g)."\n";
                foreach my $t (@{$g->{mRNAs}}){
                        print $fh to_gff3_mRNA($seg_id, $g, $t)."\n";
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

                $i++;
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

	my $parser = new Bio::Tools::GFF(-gff_version => 3,
		     			 -file        => $file,
		                         );

	my @features;
	while (my $f = $parser->next_feature()) {
		push(@features, $f);
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
sub load_gene {
	my $g   = shift;
	my $seq = shift;

	my $gene = {};

	$gene->{feature_id} = $g->{id};
	$gene->{id}         = $g->{id};

	my $src_f_id = $g->{src_f_id};

	my $nbeg = $g->{f}->start();
	my $nend = $g->{f}->end(); 

	($nbeg, $nend) = ($nend, $nbeg) if $g->{f}->strand() == -1;

	push(@{$gene->{locations}}, 
	load_feature_location($nbeg, $nend, $src_f_id));

	$gene->{name} = $g->{f}->get_Annotations('Name')->value();

	my $oF_id = $g->{id};

        foreach my $t (@{$g->{mRNAs}}){
                my $transcr = load_transcript($t, $seq);

		my $sF_id = $transcr->id();

		push(@{$transcr->{relationships}},
		load_feature_relationship($oF_id, $sF_id, 'part_of'));

		push(@{$gene->{relationships}},
        	load_feature_relationship($oF_id, $sF_id, 'part_of'));

		push(@{$gene->{transcripts}}, $transcr);

        }

	$gene->{uniquename} = $g->{f}->get_Annotations('Name')->value();
	$gene->{type}       = 'gene';

	bless $gene, 'CGL::Annotation::Feature::Gene';
	return $gene;
}
#-------------------------------------------------------------------------------
sub load_transcript {
	my $t   = shift;
	my $seq = shift;

	my $transcr = {};

	load_exons($t, $transcr);

	load_introns($t, $transcr, $seq);


	my $status = $t->{f}->get_Annotations('prediction_status')  ?
	    $t->{f}->get_Annotations('prediction_status')->value()  :
	    undef;

	$transcr->{status} = $status;


	$transcr->{feature_id} = $t->{id};
	$transcr->{id}         = $t->{id};

	$transcr->{gene} = $t->{part_of}->[0];
	
	my $nbeg = $t->{f}->start();
	my $nend = $t->{f}->end();

	($nbeg, $nend) = ($nend, $nbeg) if $t->{f}->strand() == -1;

	my $src_f_id = $t->{src_f_id};

	push(@{$transcr->{locations}},
	load_feature_location($nbeg, $nend, $src_f_id));

	my %props;
	$props{gene} = $t->{part_of}->[0];
	$props{transcript_id} = $t->{id};

	$transcr->{properties} = \%props;

	my $oF_id = $t->{id};
	foreach my $e (@{$transcr->{exons}}){
		my $sF_id =$e->{id};
		push(@{$transcr->{relationships}}, 
		load_feature_relationship($oF_id, $sF_id, 'part_of'));
	}
        foreach my $i (@{$transcr->{introns}}){
                my $sF_id =$i->{id};
                push(@{$transcr->{relationships}},
                load_feature_relationship($oF_id, $sF_id, 'part_of'));
        }

	$transcr->{residues} = $t->{seq};

	my $t_offset = get_translation_offset($t);

	$transcr->{translationStartInTranscript} = 
	{$transcr->{gene}.'protein-0' => $t_offset };

	$transcr->{translations} = load_translations($t);

        foreach my $p (@{$transcr->{translations}}){
                my $sF_id =$p->{id};
                push(@{$transcr->{relationships}},
                load_feature_relationship($oF_id, $sF_id, 'produced_by'));
        }

	$transcr->{name}       = $t->{f}->get_Annotations('Name')->value();
	$transcr->{uniquename} = $t->{f}->get_Annotations('Name')->value();

	$transcr->{type} = 'mRNA';

	bless $transcr, 'CGL::Annotation::Feature::Transcript';

	return $transcr;
}
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
sub load_translations {
	my $t = shift;

	my $f_id = 'protein:'.$t->{id};

	my $transl = {};

	$transl->{feature_id} = $f_id; 

	$transl->{id} = $f_id;
 
        my $sorted = sort_cdss($t);

	my $alpha =  $sorted->[0];
	my $omega =  $sorted->[-1];

	my $transl_offset = get_translation_offset($t);
	my $trans_end     = $transl_offset + length($t->{cds_seq}) - 3;

	my ($a_l_beg, $a_l_end) = get_l_beg_end($alpha);
	my ($o_l_beg, $o_l_end) = get_l_beg_end($omega);
	my $src_f_id = $t->{src_f_id};

	push(@{$transl->{locations}},
	load_feature_location($a_l_beg, $o_l_end, $src_f_id));

	my $trn = new Bio::Tools::CodonTable();

	my $p_seq = $trn->translate($t->{cds_seq});

	$p_seq =~ s/\*$//;

	$transl->{name} = $f_id;

	my $oF_id = $t->{id};
	my $sF_id = $transl->{id};

	push(@{$transl->{relationships}},
	load_feature_relationship( $oF_id, $sF_id, 'produced_by'));

	my %props;
	$props{codon_start} = $transl_offset;
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
	
	return 0 unless defined($t->{utr_5});

	my $sorted_utrs_5 = sort_utr_5($t);


	my $length = 0;
	foreach my $u (@{$sorted_utrs_5}){
		my ($l_b, $l_e) =  get_l_beg_end($u);
		$length += abs($l_b - $l_e) + 1; 
	}
	return $length;
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

	$exon->{feature_id} = $e->{id};
	$exon->{id}         = $e->{id};
	$exon->{inScope}    = 1;


        my $nbeg = $e->{f}->start();
        my $nend = $e->{f}->end();

	my $scr_f_id = $e->{src_f_id};

        ($nbeg, $nend) = ($nend, $nbeg)
        if $e->{f}->strand() == -1;

	push(@{$exon->{locations}}, 
	load_feature_location($nbeg, $nend, $scr_f_id));


	$exon->{name} = $e->{id};


	my $sF_id = $e->{id};

	foreach my $oF_id (@{$e->{part_of}}){
		push(@{$exon->{relationships}},
		load_feature_relationship($oF_id, $sF_id, 'part_of')); 
	}	
	$exon->{residues} = $e->{seq};

	
	$exon->{src_id} = $e->{src_f_id}.':'.$nbeg.':'.$nend;

	$exon->{type} = 'exon';

	$exon->{uniquename} = 
	$e->{src_f_id}.':'.$e->{id};

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
sub load_seq {
        my $file = shift;

        my $seqio = Bio::SeqIO->new(-file   => $file,
                                    -format => 'Fasta');
        my $seq = $seqio->next_seq;
        my $id  = $seq->display_id;

        my $sequence = $seq->seq;

        return ($id, \$sequence);
}       
#-------------------------------------------------------------------------------
sub grab {
	my $type     = shift;
	my $source   = shift;
	my $features = shift;
	my $c_id     = shift;

        my %booty;
	my $i = 0;
        foreach my $f (@{$features}){
                my $tag_t = $f->primary_tag();
                my $tag_s = $f->source_tag();

                if ($tag_s eq $source && $tag_t eq $type) {
                        my $id = 
			    $f->get_Annotations('ID') ?
			    $f->get_Annotations('ID')->value() :
			    $c_id.':'.$type.':'.$i;
                        my $p_ids = get_p_ids($f);
                        foreach my $p_id (@{$p_ids}){
                                push(@{$booty{$p_id}}, {f        => $f,
				                        id       => $id,
                                                        src_f_id => $c_id,
                                                        part_of  => $p_ids,
                                                       });
                        }
			$i++;
                }
        }
	return \%booty;
}
#-------------------------------------------------------------------------------
sub get_p_ids {
	my $f = shift;

	my @parents = $f->get_Annotations('Parent');
        my @p_ids;
       foreach my $p (@parents){
       		push(@p_ids, $p->{value});
       }

	return \@p_ids;
}
#-------------------------------------------------------------------------------
sub get_genes {
	my $gff3_file  = shift;
	my $fasta_file = shift;
	my $base_type  = shift;
	my $flank      = shift;

	my ($c_id, $seq)   = load_seq($fasta_file);

	print STDERR "loading features....\n";
	my $features = get_features($gff3_file);
	print STDERR "...finished.\n";

	print STDERR "loading $base_type....\n";
	my $exons = grab($base_type, 'Coding_transcript', $features, $c_id);
	print STDERR "...finished.\n";

	print STDERR "loading CDSs....\n";
	my $cdss  = grab('CDS','Coding_transcript', $features, $c_id);
	print STDERR "...finished.\n";

	print STDERR "loading five prime UTRs....\n";
        my $utr_5  = grab('five_prime_UTR','Coding_transcript', $features, $c_id);
        print STDERR "...finished.\n";

        print STDERR "loading three prime UTRs....\n";
        my $utr_3  = grab('three_prime_UTR','Coding_transcript', $features, $c_id);
        print STDERR "...finished.\n";

	print STDERR "loading mRNAs....\n";
	my $mRNAs = grab('mRNA','Coding_transcript', $features, $c_id);
	print STDERR "...finished.\n";
	
			

	foreach my $p_id (keys %{$mRNAs}){

		for (my $i = 0; $ i < @{$mRNAs->{$p_id}}; $i++) {
			my $f  = $mRNAs->{$p_id}->[$i]->{f};
			my $id = $mRNAs->{$p_id}->[$i]->{id};

			if (! defined $cdss->{$id}) {
				print STDERR "Warning: Gene $p_id has no CDS annotated!\n";
				next;
			}

			$mRNAs->{$p_id}->[$i]->{cdss}  = $cdss->{$id};
			$mRNAs->{$p_id}->[$i]->{utr_3} = $utr_3->{$id};
			$mRNAs->{$p_id}->[$i]->{utr_5} = $utr_5->{$id};

			$mRNAs->{$p_id}->[$i]->{exons} = manufacture_mRNA($cdss->{$id}, 
			                                                  $utr_5->{$id}, 
									  $utr_3->{$id},
									  $exons->{$p_id},
			                                                  );
		}
	}


	print STDERR "loading genes....\n";
        my @genes;
	my $j = 0;
        foreach my $f (@{$features}){
                my $tag = $f->primary_tag();
                my $source = $f->source_tag();

                if ($source eq 'Coding_transcript' && $tag eq 'gene') {
                        my $id = 
			    $f->get_Annotations('ID') ?
			    $f->get_Annotations('ID')->value() :
			    $c_id.':'.'gene'.':'.$j;
                        push(@genes, {'f'        => $f,
				      'id'       => $id,
                                      'mRNAs'    => $mRNAs->{$id},
                                      'src_f_id' => $c_id,
                                      'part_of'  => [],
                                     });
                }
		$j++;
        }

	print STDERR "...finished.\n";

	print STDERR "loading seqs ...\n";
	load_seqs(\@genes, $seq, $flank);
	print STDERR "...finished.\n";
	return (\@genes, $seq, $c_id);
}
#-------------------------------------------------------------------------------
sub get_e_coors {
	my $cdss   = shift;
	my $utrs_5 = shift;
	my $utrs_3 = shift;

	my ($alpha_l_beg, $alpha_l_end, $f_exclude_beg, $f_exclude_end) = 
	get_alpha_coors($cdss, $utrs_5);

	my ($omega_l_beg, $omega_l_end, $t_exclude_beg, $t_exclude_end) = 
	get_omega_coors($cdss, $utrs_3);

	my %exclude;

	if (defined($f_exclude_beg) && defined($f_exclude_end)){
		$exclude{$f_exclude_beg}{$f_exclude_end}++;
	}

	if (defined($t_exclude_beg) && defined($t_exclude_end)){
		$exclude{$t_exclude_beg}{$t_exclude_end}++;
	}

	my $c_size = @{$cdss};

	my %exon_coors;
	if ($c_size == 1){
		$exon_coors{$alpha_l_beg}{$omega_l_end}++;
	}
	else {
		$exon_coors{$alpha_l_beg}{$alpha_l_end}++;
		$exon_coors{$omega_l_beg}{$omega_l_end}++;

		for (my $i = 1; $i < @{$cdss} -1; $i++){

			my ($l_beg, $l_end) = get_l_beg_end($cdss->[$i]);
			$exon_coors{$l_beg}{$l_end}++;
  		}
	}

	add_utr_exons($utrs_5, \%exon_coors, \%exclude) if defined($utrs_5);
	add_utr_exons($utrs_3, \%exon_coors, \%exclude) if defined($utrs_3);
	
	return \%exon_coors;
}
#-------------------------------------------------------------------------------
sub add_utr_exons {
	my $utr_exons  = shift;
	my $exon_coors = shift;
	my $exclude    = shift;

        for (my $i = 0; $i < @{$utr_exons}; $i++){
                my ($l_beg, $l_end) = get_l_beg_end($utr_exons->[$i]);

                next if $exclude->{$l_beg}->{$l_end};

                $exon_coors->{$l_beg}->{$l_end}++;
        }

}
#-------------------------------------------------------------------------------
sub get_l_beg_end {

	my $feature = shift;

        my $start  = $feature->{f}->start();
        my $end    = $feature->{f}->end();
        my $strand = $feature->{f}->strand();

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
sub get_alpha_coors {
	my $cdss   = shift;
	my $utrs_5 = shift;


	my $sorted_cdss   = sort_quick($cdss);
        my $sorted_utrs_5 = sort_quick($utrs_5) if defined($utrs_5);

        my $first_cds = $sorted_cdss->[0];

	my ($c_l_beg, $c_l_end) = get_l_beg_end($first_cds);
	
	#warn  "no five prime UTR annotated !\n" if !defined($utrs_5);
	return ($c_l_beg, $c_l_end) if !defined($utrs_5);

	foreach my $utr_5 (@{$utrs_5}){
		my ($u_l_beg, $u_l_end) = get_l_beg_end($utr_5);

		if (abs($u_l_end - $c_l_beg) == 1){
			my $alpha_l_beg = $u_l_beg;
			my $alpha_l_end = $c_l_end;
			#warn  " ** five prime UTR found !\n";
			return ($alpha_l_beg, $alpha_l_end, $u_l_beg, $u_l_end);
		}
		else {
			#warn " u_l_end:$u_l_end c_l_beg:$c_l_beg\n";
		}
	}
	#warn  "no five prime UTR found !\n";
	#sleep 1;

	return ($c_l_beg, $c_l_end);
}
#-------------------------------------------------------------------------------
sub get_omega_coors {
        my $cdss   = shift;
        my $utrs_3 = shift;


        my $sorted_cdss   = sort_quick($cdss);
        my $sorted_utrs_5 = sort_quick($utrs_3) if defined($utrs_3);

        my $last_cds = $sorted_cdss->[-1];

        my ($c_l_beg, $c_l_end) = get_l_beg_end($last_cds);

	#warn  "no three  prime UTR annotated !\n" if !defined($utrs_3);
        return ($c_l_beg, $c_l_end) if !defined($utrs_3);

        foreach my $utr_3 (@{$utrs_3}){
                my ($u_l_beg, $u_l_end) = get_l_beg_end($utr_3);

                if (abs($u_l_beg - $c_l_end) == 1){
                        my $omega_l_beg = $c_l_beg;
                        my $omega_l_end = $u_l_end;
			#warn  " ** three prime UTR found !\n";
                        return ($omega_l_beg, $omega_l_end, $u_l_beg, $u_l_end);
                }
        }
        #warn  "no thee prime UTR found !\n";
        #sleep 1;

	return ($c_l_beg, $c_l_end);

}
#-------------------------------------------------------------------------------
sub manufacture_mRNA {
	my $cdss   = shift;
	my $utrs_5 = shift;
	my $utrs_3 = shift;
	my $exons  = shift;

#	my $g_id = 
#	    $exons->[0]->{f}->get_Annotations('Parent') ?
#	    $exons->[0]->{f}->get_Annotations('Parent')->value() :
#	    undef;

#	my $t_id = 
#	    $cdss->[0]->{f}->get_Annotations('Parent') ?
#	    $cdss->[0]->{f}->get_Annotations('Parent')->value() :
#	    undef;

	#print STDERR " XXXX new mRNA XXX $g_id $t_id XXXXXXXXXXXXX\n";

	my $exon_coors = get_e_coors($cdss, 
		                     $utrs_5, 
		                     $utrs_3,
		                     );	

	my @actual_exons;
	foreach my $e (@{$exons}){
		my ($l_beg, $l_end) = get_l_beg_end($e);

		push(@actual_exons, $e) 
		if defined($exon_coors->{$l_beg}->{$l_end});
	}

	return \@actual_exons;
}
#-------------------------------------------------------------------------------
sub debug {
        my $genes = shift;

        foreach my $g (@{$genes}) {
                print $g->{id}."\n";
                foreach my $t (@{$g->{mRNAs}}){
                        print "   ".$t->{id}."\n";
                        foreach my $e (@{$t->{exons}}){
                                print "      ".$e->{id}."\n";
                        } 
                        foreach my $c (@{$t->{cdss}}){
                                print "             ".$c->{id}."\n";
                        }
        
                }
        }
        die "DEBUG!\n";
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
sub sort_quick {
        my $array = shift;

        my @sorted;
        if    ($array->[0]->{f}->strand() ==  1){
                @sorted =
                sort {$a->{f}->start <=> $b->{f}->start} @{$array};
        }
        elsif ($array->[0]->{f}->strand() == -1){
                @sorted =
                sort {$b->{f}->start <=> $a->{f}->start} @{$array};
        }
        else {
                die "unknown strand in GFF3::WormBase::sort_quick!\n";
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
sub sort_utr_5 {
        my $t = shift;

        my @sorted;
        if    ($t->{f}->strand() ==  1){
                @sorted =
                sort {$a->{f}->start <=> $b->{f}->start} @{$t->{utr_5}};
        }
        elsif ($t->{f}->strand() == -1){
                @sorted =
                sort {$b->{f}->start <=> $a->{f}->start} @{$t->{utr_5}};
        }
        else {
                die "unknown strand in GFF3::Maker::sort_utr_5!\n";
        }
        return \@sorted;
}
#-------------------------------------------------------------------------------
sub sort_utr_3 {
        my $t = shift;

        my @sorted;
        if    ($t->{f}->strand() ==  1){
                @sorted =
                sort {$a->{f}->start <=> $b->{f}->start} @{$t->{utr_3}};
        }
        elsif ($t->{f}->strand() == -1){
                @sorted =
                sort {$b->{f}->start <=> $a->{f}->start} @{$t->{utr_3}};
        }
        else {
                die "unknown strand in GFF3::Maker::sort_utr_3!\n";
        }
        return \@sorted;
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
