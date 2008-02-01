#-------------------------------------------------------------------------------
#------                     CGL::Annotation::GFF3::DictyBase             -------
#-------------------------------------------------------------------------------
package CGL::Annotation::GFF3::DictyBase;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use FileHandle;
use Bio::DB::GFF;
use Bio::Tools::GFF;
use Bio::FeatureIO::gff;
use Bio::Tools::CodonTable;
use CGL::Annotation;

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub new {
	my $class = shift;
	my $file  = shift;

	my ($genes, $seq, $seq_id) = get_genes($file,'CDS', 0);

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
        while (defined(my $line = <$fh>)){

		next unless $line =~ /\S+/;

                if ($i < 2){
                        print $line;
                        $i++;
                        next;
                }

                chomp($line);

                my @stuff = split("\t", $line);

		if ($line =~ /##FASTA/){
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

	my $g_id   = $g->{id};
	my $g_n    = $g->{f}->get_Annotations('Name')->value();
	my $strand = $g->{f}->strand() == 1 ? '+' : '-';

	my @fields;

        push(@fields, $seg_id, 'maker', 'gene', $g->{i_start});
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

        my $t_id   = $t->{id};
        my $e_id   = $e->{id};

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

        my $t_id   = $t->{id};
        my $c_id   = $c->{id};

        my $strand = $c->{f}->strand() == 1 ? '+' : '-';

        my @fields;

        push(@fields, $seg_id, 'maker', 'CDS', $c->{i_start});
        push(@fields, $c->{i_end}, '.', $strand, '.');
        push(@fields, "ID=$c_id;Parent=$t_id");

        return join("\t", @fields);

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
	my $file = shift;
	my $root = shift;

	my ($genes, $seq, $seq_id) = get_genes($file, 'exon', 500);

	my $i = 0;
	foreach my $g (@{$genes}){
		my $g_id = $g->{id};

		my $file = $root.'/'.$seq_id.'_'.$g_id.'.gff3';
		my $fh = new FileHandle();
		   $fh->open(">$file");

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

		print $fh "##FASTA\n";
		print $fh to_gff3_seq($seg_id, $g->{seq})."\n";

		$fh->close();

		die if $i > 10;

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

        my $nbeg = $sorted->[0]->{f}->start();
        my $nend = $sorted->[0]->{f}->end();

	($nbeg, $nend) = ($nend, $nbeg) if $sorted->[0]->{f}->strand() == -1;

	my $src_f_id = $t->{src_f_id};

	push(@{$transl->{locations}},
	load_feature_location($nbeg, $nend, $src_f_id));

	my $trn = new Bio::Tools::CodonTable();

	my $p_seq = $trn->translate($t->{cds_seq});

	$p_seq =~ s/\*$//;

	$transl->{name} = $f_id;

	my $oF_id = $t->{id};
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

	my $e_0_b = $sorted_exons->[0]->{f}->start();
	my $c_0_b = $sorted_exons->[0]->{f}->start();

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
	
	my $fh = new FileHandle();
	   $fh->open($file);

	while(my $line = <$fh>){
		chomp($line);
		last if $line =~ /##FASTA/;
	}
	my $def = <$fh>;
	chomp($def);

	my $seq = '';
        while(my $line = <$fh>){
                chomp($line);
		$seq .= $line;
        }

	$fh->close();

	return ($def, \$seq);
}
#-------------------------------------------------------------------------------
sub grab {
	my $type     = shift;
	my $features = shift;
	my $c_id     = shift;

        my %booty;
        foreach my $f (@{$features}){
                my $tag_t = $f->primary_tag();
                my $tag_s = $f->source_tag();

		my $i = 0;
                if ($tag_t eq $type) {
                        my $id = $f->get_Annotations('ID') != 0 
			? $f->get_Annotations('ID')->value()
			: $c_id.':'.$type.':'.$i;
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
	my $file      = shift;
	my $base_type = shift;
	my $flank     = shift;

	my ($def, $seq)   = load_seq($file);

	my ($c_id) = $def =~ />(\S+)/;

	print STDERR "loading features...\n";
	my $features = get_features($file);
	print STDERR "... finished.\n";

	print STDERR "loading exons...\n";
	my $exons = grab($base_type, $features, $c_id);
	print STDERR "... finished.\n";
	print STDERR "loading CDSs...\n";
	my $cdss  = grab('CDS', $features, $c_id);
	print STDERR "... finished.\n";
	print STDERR "loading mRNAs...\n";
	my $mRNAs = grab('mRNA', $features, $c_id);
	print STDERR "... finished.\n";
	foreach my $p_id (keys %{$mRNAs}){
		for (my $i = 0; $ i < @{$mRNAs->{$p_id}}; $i++) {
			my $f  = $mRNAs->{$p_id}->[$i]->{f};
			my $id = $mRNAs->{$p_id}->[$i]->{id};

			$mRNAs->{$p_id}->[$i]->{exons} = $exons->{$id};
			$mRNAs->{$p_id}->[$i]->{cdss}  = $cdss->{$id};
		}
	}
	print STDERR "loading genes...\n";
        my @genes;
	my $j = 0;
        foreach my $f (@{$features}){
                my $tag = $f->primary_tag();
                my $source = $f->source_tag();

                if ($tag eq 'gene') {
                        my $id = $f->get_Annotations('ID') != 0
			? $f->get_Annotations('ID')->value()
			: $c_id.':'.'gene'.':'.$j;
                        push(@genes, {'f'        => $f,
				      'id'       => $id,
                                      'mRNAs'    => $mRNAs->{$id},
                                      'src_f_id' => $c_id,
                                      'part_of'  => [],
                                     });
                }
		$j++;
        }

	print STDERR "... finished.\n";
	print STDERR "loading seqs...\n";
	load_seqs(\@genes, $seq, $flank);
	print STDERR "... finished.\n";
	return (\@genes, $seq, $c_id);
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
