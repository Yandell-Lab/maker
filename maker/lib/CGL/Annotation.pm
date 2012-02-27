###################################################### main header begin ##

=head1 NAME

CGL::Annotation - An object for working with genome annotations.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  # make sure that there's a SO file available.
  BEGIN {
    $ENV{SO_OBO_FILE} = "sample_data/so.obo" unless $ENV{SO_OBO_FILE};
  }

  use CGL::Annotation;
  my $a;			# an annotation.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");

=for example end

=for example_testing
  isa_ok($a, "CGL::Annotation", "Check if it's the right type.");

=head1 DESCRIPTION

This is the base class for the CGL Annotation object hierarchy.  It
provides general functionality used by all of it inheritors.

=head1 USAGE

=head1 BUGS

Not yet.

=head1 AUTHOR

 Mark Yandell
 myandell@fruitfly.org
 http://www.yandell-lab.org

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

CGL::Ontology::SO

The SO_OBO_FILE environment variable, which defines the default
location for the OBO format file that describes Sequence Ontology.

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::Annotation;

use strict;
use warnings;

use CGL::Annotation::Feature;
use CGL::Annotation::Feature::Contig;
use CGL::Annotation::Feature::Exon;
use CGL::Annotation::Feature::Gene;
use CGL::Annotation::Feature::Protein;
use CGL::Annotation::Feature::Transcript;
use CGL::Annotation::Feature::Sequence_variant;
use CGL::Annotation::FeatureRelationship;
use CGL::Annotation::Iterator;
use CGL::Annotation::Trace;
use CGL::Ontology::SO;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}

my $SO = new CGL::Ontology::SO();

################################################ subroutine header begin ##

=head1 new

 Usage     :

=for example begin

  my $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");

=for example end

=for example_testing
 isa_ok($a, "CGL::Annotation");

 Purpose   : Create a new CGL::Annotation object, and optionally load
             it with data from a chaos-xml file.
 Returns   : a reference to a CGL::Annotation object.
 Argument  : the name of a chaos-xml file [optional].
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub new {
	my $class      = shift;
	my $gff_file   = shift;
	my $fasta_file = shift;
	my $format     = shift;

	my $self;

	if ($format =~ /Maker/i) {
		$self = Dumper::GFF3::GFFV3->new($gff_file, $fasta_file);
	}
	
	return $self;
}

################################################ subroutine header begin ##

=head1 transcript

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a transcript feature
  my $i;			# index into transcript list.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $i = 0;
  $t = $a->transcript($i);

=for example end

=for example_testing
  can_ok($a, qw(transcript));
  isa_ok($t, CGL::Annotation::Feature::Transcript);

 Purpose   : access a transcript from an annotation.
 Returns   : a reference to the i'th transcript feature.
 Argument  : the zero-based index identifying the transcript feature.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Transcript

=cut

################################################## subroutine header end ##

sub transcript {
	my $self = shift;
	my $i    = shift;

# XXXX does this handle an empty list correctly?
	return $self->transcripts->[$i];
}

################################################ subroutine header begin ##

=head1 translation

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $i;			# index of the translation
  my $p;			# reference to a protein feature

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $i = 0;
  $p = $a->translation($i);

=for example end

=for example_testing
  can_ok($a, qw(translation));
  isa_ok($p, CGL::Annotation::Feature::Protein);


 Purpose   : return the i'th protein feature from the annotation.
 Returns   : a scalar containing the translated sequence.
 Argument  : a zero-based integer identifying the protein feature.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Protein

=cut

################################################## subroutine header end ##

sub translation {
        my $self = shift;
        my $i    = shift;

# XXXX should this be "or undef";
        return $self->translations->[$i];
}

################################################ subroutine header begin ##

=head1 get_gene_by_id

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene object
  my $i;			# gene identifier

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");

  $i = "gene-107698";			# identifier for gene
  $g = $a->get_gene_by_id($i);

=for example end

=for example_testing
  can_ok($a, qw(get_gene_by_id));
  isa_ok($g, CGL::Annotation::Feature::Gene);

 Purpose   : return a reference to the gene object with the requested
             id from the annotation.
 Returns   : a reference to a gene object, undef if there isn't a gene
             with that identifier.
 Argument  : a string containing the gene id.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Gene

=cut

################################################## subroutine header end ##

sub get_gene_by_id {
	my $self = shift;
	my $id   = shift;

	my $i=0;
	while (my $g = $self->gene($i)){
		return $g if $g->id() eq $id;
		$i++;
	}
	return undef;
}

################################################ subroutine header begin ##

=head1 get_gene_by_name

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene feature.
  my $n;			# gene name

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $n = "At3g23010";			# name of gene of interest
  $g = $a->get_gene_by_name($n);

=for example end

=for example_testing
  can_ok($a, qw(get_gene_by_name));
  isa_ok($g, CGL::Annotation::Feature::Gene);

 Purpose   : return a reference to the gene object with the requested
             name from the annotation.
 Returns   : a reference to a gene object, undef if there isn't a gene
             with that name.
 Argument  : a string containing the gene name.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Gene

=cut

################################################## subroutine header end ##

sub get_gene_by_name {
	my $self = shift;
	my $name   = shift;

	my $i=0;
	while (my $g = $self->gene($i)){
		return $g if $g->name() eq $name;
		$i++;
	}
	return undef;
}

################################################ subroutine header begin ##

=head1 get_gene_by_uniquename

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene object
  my $u;			# unique gene name

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $u = "At3g23010-gene-107698";		# uniquename for gene of interest
  $g = $a->get_gene_by_uniquename($u);

=for example end

=for example_testing
  can_ok($a, qw(get_gene_by_uniquename));
  isa_ok($g, CGL::Annotation::Feature::Gene);

 Purpose   : return a reference to the gene object with the requested
             uniquename from the annotation.
 Returns   : a reference to a gene object, undef if there isn't a gene
             with that uniquename.
 Argument  : a string containing the gene name.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Gene

=cut

################################################## subroutine header end ##

sub get_gene_by_uniquename {
	my $self = shift;
	my $uniquename   = shift;

	my $i=0;
	while (my $g = $self->gene($i)){
		return $g if $g->uniquename() eq $uniquename;
		$i++;
	}
	return undef;
}

################################################ subroutine header begin ##

=head1 gene

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to a gene object
  my $i;			# index of gene feature.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $i = 0;			# identifier for gene of interest
  $g = $a->gene($i);

=for example end

=for example_testing
  can_ok($a, qw(gene));
  isa_ok($g, CGL::Annotation::Feature::Gene);

 Purpose   : return a reference to the i'th gene feature in the
             annotation.
 Returns   : a reference to a gene object, undef if there
             isn't an i'th gene.
 Argument  : an integer, the index into the list of genes.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Gene

=cut

################################################## subroutine header end ##
sub gene {
        my $self = shift;
        my $i    = shift;

	# XXXX should this be "or undef"?
        return $self->genes->[$i];
}

################################################ subroutine header begin ##

=head1 genes

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# reference to an array of gene objects

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->genes();

=for example end

=for example_testing
  can_ok($a, qw(genes));
  isa_ok($g, ARRAY);
  isa_ok($g->[0], CGL::Annotation::Feature::Gene);

 Purpose   : return a reference to a list of gene objects from the
             annotation.
 Returns   : a reference to a (possibly empty) list of gene objects,
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Gene
             CGL::Ontology::SO

=cut

################################################## subroutine header end ##

sub genes {
        my $self = shift;

        return $self->{genes} if defined($self->{genes});

        my @genes;
        foreach my $g ($self->featuresByType('gene')){
                my $gene = new CGL::Annotation::Feature::Gene($g);
                foreach my $part ($self->trace($g, 'parts')->features(0)){
			if ($SO->a_is_hyponym_of_b($part->type(), 'transcript')){
                                my $t = new CGL::Annotation::Feature::Transcript($part);
				foreach my $part ($self->trace($t, 'parts')->features(0)){
					if ($SO->a_is_hyponym_of_b($part->type(), 'exon')){
						my $exon = new CGL::Annotation::Feature::Exon($part);
						   $exon->_add_residues($self->contig(0));
						   $exon->_add_src_id($self->contig(0));

						$t->_add_exon($exon);
					}
				}
				foreach my $part ($self->trace($t, 'produces')->features(0)){
				        if (($part->type() eq 'protein') ||
					    ($part->type() eq 'polypeptide') ||
					    ($SO->a_is_hypomeronym_of_b($part->type(), 'protein')) ||
					    ($SO->a_is_hypomeronym_of_b($part->type(), 'polypeptide'))) {
						my $translation =
						new CGL::Annotation::Feature::Protein($part);

						$t->_add_translation($translation);
					}
				}
			        #$t->_add_residues($self->contig(0));
				$t->_add_introns($self->contig(0));
                                $gene->_add_transcript($t);
                        }
                }
                push(@genes, $gene);
        }
        $self->{genes} = \@genes;
# XXXX should it explicitly return the reference....
}

################################################ subroutine header begin ##

=head1 transcripts

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a list of transcripts

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = $a->transcripts();

=for example end

=for example_testing
  can_ok($a, qw(transcripts));
  isa_ok($t, ARRAY);
  isa_ok($t->[0], CGL::Annotation::Feature::Transcript);

 Purpose   : return a reference to a list of transcript objects
             from the annotation.
 Returns   : a reference to a (possibly empty) list of
             transcript objects.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Transcript

=cut

################################################## subroutine header end ##

sub transcripts {
        my $self = shift;


        return $self->{transcripts} if defined($self->{transcripts});

        my @transcripts;

	foreach my $t ($self->features){
		next unless $SO->a_is_hyponym_of_b($t->type(), 'transcript');
                my $transcript = new CGL::Annotation::Feature::Transcript($t);
                foreach my $part ($self->trace($t, 'parts')->features(0)){
		        if ($SO->a_is_hyponym_of_b($part->type(), 'exon')){
                                my $exon = new CGL::Annotation::Feature::Exon($part);
                                $transcript->_add_exon($exon);
                        }
                }
                foreach my $part ($self->trace($t, 'produces')->features(0)){
			if (($part->type() eq 'protein') ||
			    ($part->type() eq 'polypeptide') ||
			    ($SO->a_is_hypomeronym_of_b($part->type(), 'protein')) ||
			    ($SO->a_is_hypomeronym_of_b($part->type(), 'polypeptide'))) {
                                my $translation =
                                new CGL::Annotation::Feature::Protein($part);

                                $transcript->_add_translation($translation);
                        }
                }

                $transcript->_add_introns($self->contig(0));
                push(@transcripts, $transcript);
        }
        $self->{transcripts} = \@transcripts;
	# XXXX should we explicitly return() here?
}

################################################ subroutine header begin ##

=head1 contig

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $i;			# index of the contig of interest.
  my $c;			# reference to a contig feature object

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $i = 0;
  $c = $a->contig($i);

=for example end

=for example_testing
  can_ok($a, qw(contig));
  isa_ok($c, CGL::Annotation::Feature::Contig);

 Purpose   : return a reference to the i'th contig feature in the annotation.
 Returns   : a reference to the i'th contig feature, undef if there isn't an
             i'th contig.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Contig

=cut

################################################## subroutine header end ##

sub contig {
        my $self = shift;
        my $i    = shift;

        return $self->contigs->[$i];
}

################################################ subroutine header begin ##

=head1 contigs

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $c;			# reference to a list of contigs

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $c = $a->contigs();

=for example end

=for example_testing
  can_ok($a, qw(contigs));
  isa_ok($c, ARRAY);
  isa_ok($c->[0], CGL::Annotation::Feature::Contig);

 Purpose   : return a reference to a list of contig objects from the
             annotation.
 Returns   : a reference to a (possibly empty) list of contig objects.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Contig

=cut

################################################## subroutine header end ##

sub contigs {
        my $self = shift;

        return $self->{contigs} if defined($self->{contigs});

        my @contigs;
	foreach my $s ($self->features){
		next unless $SO->a_is_hyponym_of_b($s->type(), 'contig');
                my $contig = new CGL::Annotation::Feature::Contig($s);
                foreach my $part ($self->trace($s, 'parts')->features(0)){
			# add stuff here
                }
                push(@contigs, $contig);
        }
	# XXXX explicit return?
        $self->{contigs} = \@contigs;
}
#-------------------------------------------------------------------------------
sub sequence_variants {
        my $self = shift;

        return $self->{sequence_variants} if defined($self->{sequence_variants});

        foreach my $v ($self->featuresByType('sequence_variant')){
                my $poly = new CGL::Annotation::Feature::Sequence_variant($v);
                $self->_add_variant($poly);
        }
        return $self->{sequence_variants} || [];
}
#-------------------------------------------------------------------------------
sub _add_variant {
        my $self = shift;
        my $f    = shift;

        push(@{$self->{sequence_variants}}, $f);
}
#-------------------------------------------------------------------------------
################################################ subroutine header begin ##

=head1 translations

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# reference to a list of translations

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->translations();

=for example end

=for example_testing
  can_ok($a, qw(translations));
  isa_ok($g, ARRAY);
  isa_ok($g->[0], CGL::Annotation::Feature::Protein);

 Purpose   : return a reference to a list of translations from the
             annotation.
 Returns   : a reference to a (possibly empty) list of translations.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Translation

=cut

################################################## subroutine header end ##

sub translations {
        my $self = shift;

        return $self->{translations} if defined($self->{translations});

        my @translations;
	foreach my $p ($self->features){
		next unless ($p->type() eq 'protein' ||
			     $p->type() eq 'polypeptide' ||
			     $SO->a_is_hyponym_of_b($p->type(), 'protein')||
			     $SO->a_is_hyponym_of_b($p->type(), 'polypeptide'));
                my $translation = new CGL::Annotation::Feature::Protein($p);
                foreach my $producer ($self->trace($p, 'producers')->features(0)){
			if ($SO->a_is_hyponym_of_b($producer->type(), 'transcript')){
                                my $transcript =
				new CGL::Annotation::Feature::Transcript($producer);
                                $translation->_add_transcript($transcript);
                        }
                }
                push(@translations, $translation);
        }
	# XXXX explicit return?
        $self->{translations} = \@translations;
}

################################################ subroutine header begin ##

=head1 exons

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $l;			# a list of exon features.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $l = $a->exons();

=for example end

=for example_testing
  can_ok($a, qw(exons));
  isa_ok($l, ARRAY);
  isa_ok($l->[0], CGL::Annotation::Feature::Exon);

 Purpose   : return a reference to a list of gene objects from the
             annotation.
 Returns   : A (possibly empty) list of exon features.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Exon

=cut

################################################## subroutine header end ##

sub exons {
	my $self = shift;
	my @feats;
	my @exons;

	@feats = $self->featuresByType('exon');

	@exons = map {new CGL::Annotation::Feature::Exon($_)}  @feats;

	return \@exons;
}

################################################ subroutine header begin ##

=head1 featuresByType

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $t;			# the type, as a string.
  my @l;			# a list of the features.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $t = 'exon';
  @l = $a->featuresByType($t);

=for example end

=for example_testing
  can_ok($a, qw(featuresByType));
  isa_ok(\@l, ARRAY);           # Weird, take the ref so that isa_ok works...
  isa_ok($l[0], CGL::Annotation::Feature);
  isa_ok(new CGL::Annotation::Feature::Exon($l[0]),
         CGL::Annotation::Feature::Exon);

 Purpose   : retrieve a set of Feature objects from an annotation.
 Returns   : a list of Feature objects from the annotation.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub featuresByType {
        my $self = shift;
	my $type = shift;

	my @f;
        foreach my $f ($self->features){
                push(@f, $f) if $f->type eq $type;
        }

        return @f;
}

################################################ subroutine header begin ##

=head1 feature

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $id;			# the id of the feature.
  my $f;			# reference to a Feature object

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $id = "NC_003074.1";
  $f = $a->feature($id);

=for example end

=for example_testing
  can_ok($a, qw(feature));
  isa_ok($f, CGL::Annotation::Feature);
  is($f->id(), "NC_003074.1", "Check for correct feature id.");

 Purpose   : retrieve a feature with a particular id from the
             annotation.
 Returns   : a reference to a feature object, undef if there is
             no feature with the requested id.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub feature {
	my $self   = shift;
	my $id     = shift;

	foreach my $f ($self->features){
		return $f if $f->id() eq $id;
	}
	return(undef);
}

################################################ subroutine header begin ##

=head1 features

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @f;			# a list of features

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @f = $a->features();

=for example end

=for example_testing
  can_ok($a, qw(features));
  isa_ok(\@f, ARRAY);           # Weird, take the ref so that isa_ok works...


 Purpose   : retrieve the set of features from an annotation.
 Returns   : a (possibly empty) list of feature objects.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature

=cut

################################################## subroutine header end ##

sub features {
	my $self = shift;

	return @{$self->{features}|| []};
}

################################################ subroutine header begin ##

=head1 relationships

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @r;			# a list of the NodeRelationship objects in $a.

  $a = new CGL::Annotation;
  @r = $a->relationships();

=for example end

=for example_testing
  can_ok($a, qw(relationships));
  isa_ok(\@r, ARRAY);

 Purpose   : retrieve information about the annotation's
             feature relationships
 Returns   : a (possibly emtpy) list of NodeRelationship objects.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  : CGL::Ontology::NodeRelationship

=cut

################################################## subroutine header end ##

sub relationships {
        my $self = shift;

        return @{$self->{relationships}|| []};
}


################################################ subroutine header begin ##

=head1 trace

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $a_rab1;			# reference to a CGL::Annotation
  my $f;			# the interesting feature
  my $f_rab1;			# the interesting feature
  my $t_parts;			# the trace through the parts relationship.
  my $t_whole;			# the trace through the whole relationship.
  my $t_produces_rab1;		# the trace through the parts relationship.
  my $t_producers_rab1;		# the trace through the whole relationship.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $f = $a->feature('gene-107698');
  $t_parts = $a->trace($f, 'parts');
  $f = $a->feature('mRNA-107699');
  $t_whole = $a->trace($f, 'wholes');

  # the rab1 file uses 'derives_from' instead of 'produced_by'
  $a_rab1 = new CGL::Annotation("sample_data/Rab1.chaos.xml");
  $f_rab1 = $a_rab1->feature('mRNA:EMBL/GenBank/SwissProt:AE003734:52204:55287');
  $t_produces_rab1 = $a_rab1->trace($f_rab1, 'produces');
  $f_rab1 = $a_rab1->feature('AAF55873.1');
  $t_producers_rab1 = $a_rab1->trace($f_rab1, 'producers');

=for example end

=for example_testing
  can_ok($a, qw(trace));
  isa_ok($t_parts, CGL::Annotation::Trace, "Trace for parts.");
  ok(scalar(@{$t_parts->{features}}) > 0, "Check parts");
  isa_ok($t_whole, CGL::Annotation::Trace, "Trace for wholes.");
  ok(scalar(@{$t_whole->{features}}) > 0, "Check wholes");
  can_ok($a_rab1, qw(trace));
  isa_ok($t_produces_rab1, CGL::Annotation::Trace, "Trace rab1 for produces.");
  ok(scalar(@{$t_produces_rab1->{features}}) > 0, "Check rab1 produces");
  isa_ok($t_producers_rab1, CGL::Annotation::Trace, "Trace rab1 for produces.");
  ok(scalar(@{$t_producers_rab1->{features}}) > 0, "Check rab1 producers");

 Purpose   : trace through the relations for a feature in an annotation.
 Returns   : a reference to a trace object.
 Argument  : A Feature to trace from.
             A type of trace (as a string).  Valid values include:
               'parts'
               'producers'
               'wholes'
               'produces'
 Throws    :
 Comments  : Can also take an optional index and trace, used
           : internally for recursing
           : among the relationships.
           :
           : In newer SO releases, the typedef derived_from has replaced
           : the produced_by typedef.  This sub currently works with either.
 See Also  : CGL::Annotation::Trace

=cut

################################################## subroutine header end ##

sub trace {
  my $self  = shift;
  my $f     = shift;
  my $type  = shift;
  my $i     = shift;
  my $trace = shift;

  if (defined($i)){
    $i++;
  }
  else {
    $i = 0;
  }

  $trace = new CGL::Annotation::Trace unless defined($trace);

  foreach my $r ($f->relationships){
    if    ($type eq 'parts'){
      next unless ($r->logus eq 'part_of' ||
		   $r->logus eq 'member_of');
      next unless $r->oF eq $f->id;
      $trace->add_feature($self->feature($r->sF), $i);
      $self->trace($self->feature($r->sF), $type, $i, $trace);
    }
    elsif ($type eq 'wholes'){
      next unless ($r->logus eq 'part_of' ||
		   $r->logus eq 'member_of');
      next unless $r->sF eq $f->id;	
      $trace->add_feature($self->feature($r->oF), $i);
      $self->trace($self->feature($r->oF), $type, $i, $trace);
    }
    elsif ($type eq 'produces'){
      next unless ($r->logus eq 'produced_by' ||
		   $r->logus eq 'derives_from');
      next unless $r->oF eq $f->id;
      $trace->add_feature($self->feature($r->sF), $i);
      $self->trace($self->feature($r->sF), $type, $i, $trace);
    }
    elsif ($type eq 'producers'){
      next unless ($r->logus eq 'produced_by' ||
		   $r->logus eq 'derives_from');
      next unless $r->sF eq $f->id;
      $trace->add_feature($self->feature($r->oF), $i);
      $self->trace($self->feature($r->oF), $type, $i, $trace);
    }
    else {
      die "unknown type($type) in CGL::Annotation::trace!\n";
    }
  }

  return $trace;

}

################################################ subroutine header begin ##

=head1 metadata

 Usage     :

=for example begin

  my $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  my $export_user;		# a reference to the metadata.
  $export_user = $a->meta_data("export_user");

=for example end

=for example_testing
  can_ok($a, qw(meta_data));
  is($export_user, "cjm", "Check for correct export_user in xml file");


 Purpose   : Access chaos-xml metadata.
 Returns   : a ref
 Argument  : an optional key.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub meta_data {
	my $self = shift;
	my $key  = shift;

	if (defined($key)){
		return $self->{chaos_metadata}->[0]->{chaos_metadata}->{$key};
	}
	else {
		return $self->{chaos_metadata}->[0]->{chaos_metadata};
	}

	return undef;
}
#-------------------------------------------------------------------------------
#--------------------------------- PRIVATE -------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub _add_feature {
	my $self = shift;
	my $f    = shift;
	
	push(@{$self->{features}}, $f);
}
#-------------------------------------------------------------------------------
sub _add_relationship {
        my $self = shift;
        my $r    = shift;

        push(@{$self->{relationships}}, $r);
}
#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub transcript_is_in_scope {
        my $t = shift;

        foreach my $e (@{$t->exons}){
                return 0 unless $e->inScope();
        }

        return 1;
}
#-----------------------------------------------------------------------------
sub parse {
        my $e    = shift;
        my $tags = shift;
        my $type  = shift;
	
	my $hell = $e->nodeName();

        next unless ref($e) eq 'XML::LibXML::Element';

        foreach my $c (@{$e->childNodes()}){
		next unless ref($c) eq 'XML::LibXML::Element';
                if     ($c->nodeName() eq 'featureprop'){
                        load_featureprop($c, $tags, $type);
                        next;
                }
                elsif ($c->nodeName() eq 'featureloc'){
                        load_featureloc($c, $tags, $type);
                        next;
                }
                elsif ($c->nodeName() eq 'placeholder'){
                        #load_place_holder($c, $tags, $type);
                        #next;
                }
                else {
                        load_tag($c, $tags, $type);
                        next;
                }
                parse($c, $tags);
        }
        return $tags;
}
#-------------------------------------------------------------------------------
sub load_tag {
        my $e    = shift;
        my $tags = shift;
	my $type = shift;

	my $n = $e->textContent();

	 $n =~ s/[\s\n]//g;

	$tags->{$type}->{$e->nodeName} = $n;
}
#-------------------------------------------------------------------------------
sub _load_snp_struct {
	my $e   = shift;
	my $key = shift;

        my %data;       
	my $eName = $e->nodeName();
        foreach my $c (@{$e->childNodes()}){
        	my $pKey = $c->nodeName();
		my $pVal = $c->textContent();
		   $pVal =~ s/[\s\n]//g;
               	foreach my $d (@{$c->childNodes()}){
                       	my $pKey = $d->nodeName();
                        my $pVal = $d->textContent();
                           $pVal =~ s/[\s\n]//g;

                        next unless defined($pVal);
			next if $pKey eq 'text';

			if ($pKey eq 'snpfxn'){
				foreach my $f (@{$d->childNodes()}){
					my $pKey = $f->nodeName();
                        		my $pVal = $f->textContent();
                           		   $pVal =~ s/[\s\n]//g;

                        		next unless defined($pVal);
                        		next if $pKey eq 'text';

					$data{$pKey} = $pVal;
				}
			}
			else {
                        	$data{$pKey} = $pVal;
			}
              }
	}
	return \%data;
}
#-------------------------------------------------------------------------------
sub load_featureprop {
	my $e    = shift;
	my $tags = shift;
	my $type = shift;

	
        my $pKey = $e->getElementsByTagName('type')->[0]->textContent();
        my $pVal = $e->getElementsByTagName('value')->[0]->textContent();

	#-- dealing with snpstructs from CJM
	if ($pKey eq 'snpstruct'){
		my $data = _load_snp_struct($e, $pKey);

		# flatten the snpstruct & snpfxn contents
		foreach my $k (keys %{$data}){
			$tags->{$type}->{'properties'}->{$k} = $data->{$k};
		}
                return;

	}
	#-- end;
	$pVal =~ s/[\s\n]//g;

	$tags->{$type}->{'properties'}->{$pKey} = $pVal;

}
#-------------------------------------------------------------------------------
sub load_featureloc {
        my $e    = shift;
        my $tags = shift;
        my $type = shift;

	my %locs;
	foreach my $c (@{$e->childNodes()}){
	        my $pKey = $c->nodeName();
		my $pVal;
		if (defined($c->textContent())){
			$pVal = $c->textContent();
			$pVal =~ s/[\s\n]//g;
			$locs{$pKey} = $pVal;
		}
		else {
			$locs{$pKey} = '';
		}
	}	
	push(@{$tags->{$type}->{'locations'}}, \%locs);
}

#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;
	if ($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::AutoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}
        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#----------------------------------------------------------------------------

1; #this line is important and will help the module return a true value
__END__

