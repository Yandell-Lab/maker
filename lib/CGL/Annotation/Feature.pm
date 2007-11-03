###################################################### main header begin ##

=head1 NAME

CGL::Annotation::Feature - The CGL annotation feature object.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  # make sure that there's a SO file available.
  BEGIN {
    $ENV{SO_OBO_FILE} = "sample_data/so.obo" unless $ENV{SO_OBO_FILE};
  }

  use CGL::Annotation::Feature;
  my $f = new CGL::Annotation::Feature;

=for example end

=for example_testing
  isa_ok($f, "CGL::Annotation::Feature", "Check if it's the right type.");

=head1 DESCRIPTION

Stub documentation for this module was created by
ExtUtils::ModuleMaker.  And, then it was poked, prodded, and otherwise
massaged into it's current form by George.

Hopefully the module author wasn't negligent enough to leave the stub
unedited.

Blah blah blah.

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



=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

###################################################### main header begin ##

=head1 CGL::Annotation::Feature

CGL::Annotation::Feature - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use CGL::Annotation::Feature;
  my $foo = new CGL::Annotation::Feature;

=for example end

=for example_testing
  isa_ok($foo, "CGL::Annotation::Feature", "Check if it's the right type.");

=head1 DESCRIPTION

Stub documentation for this module was created by
ExtUtils::ModuleMaker.  And, then it was poked, prodded, and otherwise
massaged into it's current form by George.

Hopefully the module author wasn't negligent enough to leave the stub
unedited.

Blah blah blah.

=head1 USAGE

Expand on the examples from the SYNOPSIS.

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

List other relevant resources.

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::Annotation::Feature;

use strict;
use warnings;

use CGL::Annotation::FeatureLocation;
use CGL::Ontology::SO;
use CGL::Revcomp qw(revcomp);

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}



################################################ subroutine header begin ##

=head1 new

 Usage     :

=for example begin

  use CGL::Annotation::Feature;
  my $feat = new CGL::Annotation::Feature;

=for example end
=for example_testing
  isa_ok($feat, "CGL::Annotation::Feature", "Check if it's the right type.");

 Purpose   : Create a new CGL::Annotation::Feature object.
 Returns   : A reference to the object.
 Argument  : XXXX
 Throws    :
 Comments  : XXXX Never used anywhere.
           :
 See Also  :

=cut

################################################## subroutine header end ##
sub new {
	my $class = shift;
	my $hash  = shift;

	my $self = _load_feature($hash->{feature});


	bless($self, $class);

	$self->_add_id('feature_id'); # specifies binding, could change!
	return $self;
}

################################################ subroutine header begin ##

=head1 residues

 Usage     :

=for example begin

  use CGL::Annotation;

  my $a;			# an annotation.
  my $some_contig;		# a feature.
  my $residues;			# the residues.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $some_contig = $a->feature("contig-3283");
  $residues = $some_contig->residues();

=for example end

=for example_testing
  like($residues, qr/^TTTAT.*/, "Check the residues for contig-3283.");

 Purpose   : set/get the residues value.
 Returns   : The value of the residues field.
 Argument  : The value for the residues field if setting it.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub residues {
	my $self     = shift;
	my $residues = shift;

	if (defined($residues)){
		$self->{residues} = $residues;
	}
	else {
		return $self->{residues};
	}
	#XXXX default to returning residues?
}

################################################ subroutine header begin ##

=head1 nbeg

 Usage     :

=for example begin

  use CGL::Annotation;

  my $a;			# an annotation.
  my $feature;			# a feature.
  my $nbeg;			# its natural begin.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $feature = $a->feature("contig-3283");

  $nbeg = $feature->nbeg();

=for example end

=for example_testing
  is($nbeg, 8174345, "Check the natural begin of contig-3283.");


 Purpose   : set/get the nbeg value.
 Returns   : The value of the nbeg field.
 Argument  : The value for the nbeg field if setting it.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub nbeg {
	my $self = shift;
	my $nbeg = shift;

	if (defined($nbeg)){
		$self->featureLocation(0)->nbeg($nbeg);
	}
	else {
		return $self->featureLocation(0)->nbeg();
	}
	# XXXX default return?
}

################################################ subroutine header begin ##

=head1 nend

 Usage     :

=for example begin

  use CGL::Annotation;

  my $a;			# an annotation.
  my $feature;			# a feature.
  my $nend;			# its natural end.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $feature = $a->feature("contig-3283");

  $nend = $feature->nend();

=for example end

=for example_testing
  is($nend, 8177133, "Check the natural end of contig-3283.");


 Purpose   : Set/get the natural end of the feature.
 Returns   : A natural end of the feature, as a scalar.
 Argument  : The natural end to set the value, none otherwise.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub nend {
	my $self = shift;
        my $nend = shift;

        if (defined($nend)){
                $self->featureLocation(0)->nend($nend);
        }
        else {
		return $self->featureLocation(0)->nend();
        }

}

################################################ subroutine header begin ##

=head1 length

 Usage     :

=for example begin

  use CGL::Annotation;

  my $a;			# an annotation.
  my $feature;			# a feature.
  my $length;			# its length.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $feature = $a->feature("contig-3283");

  $length = $feature->length();

=for example end

=for example_testing
  is($length, 2788, "Check the length of contig-3283.");

 Purpose   : Get the length of the feature.
 Returns   : The length as a scalar.
 Argument  : none.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub length {
	my $self = shift;

	my $b = $self->nbeg();
	my $e = $self->nend();

	return abs($e - $b);
}
################################################ subroutine header begin ##

=head1 uniquename

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $u;			# the feature's unique name.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $u = $l[0]->uniquename();

=for example end

=for example_testing
  is($u, "At3g23010-gene-107698", "Check the unique name getter.");


 Purpose   : Get the Feature's uniquename.
 Returns   : The value of the residues field.
 Argument  :
 Throws    :
 Comments  : Replaces any occurences of "/" in the uniquename with
           : "-".
 See Also  :

=cut

################################################## subroutine header end ##
sub uniquename {
	my $self = shift;
	
	my $u = $self->{uniquename};

	$u =~ s/\//-/g;

	return $u;
}

################################################ subroutine header begin ##

=head1 featureLocation XXXX NOT FINISHED

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $location;			# the feature's unique name.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $location = $l[0]->featureLocation();

=for example end

=for example_testing
  isa_ok($location, "CGL::Annotation::FeatureLocation", "Check return class.");


 Purpose   : Get the Feature's uniquename.
 Returns   : The value of the residues field.
 Argument  :
 Throws    :
 Comments  : Replaces any occurences of "/" in the uniquename with
           : "-".
 See Also  :

=cut

################################################## subroutine header end ##

sub featureLocation {
	my $self = shift;
	my $i    = shift;

	if ($i > 0){
	  die "dead in CGL::Annotation::Feature::featureLocation(); i > 0 \n";
	}
        if (!defined($self->{locations})){
	  die;			# XXXX fix message	
        }

	return $self->{locations}->[$i];
}

################################################ subroutine header begin ##

=head1 name

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $n;			# the feature's unique name.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $n = $l[0]->name();

=for example end

=for example_testing
  is($n, "At3g23010", "Check the name getter.");


 Purpose   : Get the Feature's name.
 Returns   : The name of the Feature, as a scalar.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub name {
        my $self = shift;

        return $self->{name};

}

################################################ subroutine header begin ##

=head1 type

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $t;			# the feature's unique name.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $t = $l[0]->type();

=for example end

=for example_testing
  is($t, "gene", "Check the type getter.");

 Purpose   : Get the Feature's type.
 Returns   : The type of the Feature, as a scalar.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub type {
        my $self = shift;

        return $self->{type};
}

################################################ subroutine header begin ##

=head1 strand

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $s;			# the feature's unique name.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $s = $l[0]->strand();

=for example end

=for example_testing
  is($s, 1, "Check the strand getter.");

 Purpose   : Get the Feature's type.
 Returns   : The type of the Feature, as a scalar. (1 or -1)
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub strand {
        my $self = shift;


	my $s = $self->nbeg < $self->nend ? 1 : -1;

	#my $s = $self->featureLocation(0)->strand();

	die "strand undefined in CGL::Annotation::Feature::strand\n"
	unless ($s == 1 || $s == -1);

        return $s;
}

################################################ subroutine header begin ##

=head1 id

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $i;			# the feature's id.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $i = $l[0]->id();

=for example end

=for example_testing
  is($i, "gene-107698", "Check the id getter.");

 Purpose   : Get the Feature's id.
 Returns   : The type of the Feature, as a scalar.
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub id {
	my $self = shift;

	return $self->{id};
}

################################################ subroutine header begin ##

=head1 inScope

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $i;			# the feature's id.
  my $is_in_scope;
  my $is_not_in_scope;

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $l[0]->inScope(0);
  $is_not_in_scope =   $l[0]->inScope();
  $l[0]->inScope(1);
  $is_in_scope =   $l[0]->inScope();

=for example end

=for example_testing
  is($is_not_in_scope, 0, "Did setting inScope to false work?");
  is($is_in_scope, 1, "Did setting inScope to true work?");


 Purpose   : Get or set whether the feature's in scope or not.
 Returns   : 1 or -1
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub inScope {
	my $self = shift;
	my $arg  = shift;

	if (defined($arg)){
		$self->{inScope} = $arg;
	}
	else {
		return $self->{inScope};
	}
}

################################################ subroutine header begin ##

=head1 properties

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my $g;			# a gene from the annotation.
  my $note;			# a note from the featureprop.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  $g = $a->gene(0);
  $note = $g->properties('note');

=for example end

=for example_testing
  isa_ok($g, CGL::Annotation::Feature::Gene, "moose");
  is($note,
     'MXC7.4;similartoHcr25bGB:AAC78595(Lycopersiconesculentum)(PlantCell10,19151926(1998));containsPfamprofile:PF00560leucinerichrepeat(17copies)',
     "Check that we properly retrieve a note property."); 

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub properties {
	my $self = shift;
	my $pKey = shift;
	my $pVal = shift;

	if (defined($pVal)){
		$self->{properties}->{$pKey} = $pVal;
	}
	else {
		return $self->{properties}->{$pKey};
	}
}

################################################ subroutine header begin ##

=head1 relationships

 Usage     :

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my @r;			# the feature's relationships.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  @r = $l[0]->relationships();

=for example end

=for example_testing
  isa_ok($r[0], CGL::Annotation::FeatureRelationship, "Check relationships");

 Purpose   : Get the feature's relationshiops.
 Returns   : An array of the features (possibly empty).
 Argument  : None.
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub relationships {
        my $self = shift;

        return @{$self->{relationships}|| []};
}

################################################ subroutine header begin ##

=head1 _add_id

 Usage     : *private*

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $l[0]->_add_id('feature_id');


=for example end

=for example_testing
  is($l[0]->{id}, 'gene-107698', "Check that _add_id did the right thing.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub _add_id {
	my $self    = shift;
	my $binding = shift;

	$self->{id} = $self->{$binding};
}

################################################ subroutine header begin ##

=head1 _so

 Usage     : *private*

=for example begin

  use CGL::Annotation;
  my $a;			# reference to a CGL::Annotation
  my @l;			# a list of gene features.
  my $so;			# a reference to out so object.

  $a = new CGL::Annotation("sample_data/atha.sample.chaos.xml");
  @l = $a->featuresByType('gene');
  $so = $l[0]->_so();

=for example end

=for example_testing
  isa_ok($so, "CGL::Ontology::SO", "Check if it's the right type.");

 Purpose   : Provide access to a common SO object.
 Returns   : A reference to a SO object.
 Argument  :
 Throws    :
 Comments  : The shared so object is part of a closure.
           :
 See Also  : CGL::Ontology::SO

=cut

################################################## subroutine header end ##

# XXXX why a closure instead of a class variable, or....
{
  my $so;
  sub _so {
    my $self = shift;

    $so = new CGL::Ontology::SO() unless (defined($so));

    return $so;
  }
}

################################################ subroutine header begin ##

=head1 _add_residues

 Usage     : *private*

=for example begin

=for example end

=for example_testing

 Purpose   : Add residues from a contig to a feature.
 Returns   : Whether the feature is "inScope".
 Argument  : A reference to a contig.
 Throws    : die()s if either the Feature or the contig has an invalid
           : strand.
 Comments  :
 See Also  :

=cut

################################################## subroutine header end ##

sub _add_residues
 {
	my $self   = shift;
	my $contig = shift;

	my $nB = $self->nbeg();
	my $nE = $self->nend();

	if ($nB < 0 || $nE < 0){
		$self->{residues} = undef;
		$self->inScope(0);
		goto done;
	}
	
        if (CORE::length($contig->residues()) < $nB ||
	    CORE::length($contig->residues()) < $self->length()){
		$self->inScope(0);
		goto done;
        }

	if ($nB > $contig->length() || $nE > $contig->length()){
		$self->{residues} = undef;
		$self->inScope(0);
		goto done;
	}

	if    ($contig->strand() == 1 && $self->strand() == 1){
		$self->residues(substr($contig->residues,
				       $nB, $self->length()));
		$self->inScope(1);
	}
	elsif ($contig->strand() == 1 && $self->strand() == -1){
		$self->residues(revcomp(substr($contig->residues,
					       $nE, $self->length())));
		$self->inScope(1);
	}
        elsif ($contig->strand() == -1 && $self->strand() == 1){
                $self->residues(substr($contig->residues,
				       $nB, $self->length()));
                $self->inScope(1);
        }
        elsif ($contig->strand() == -1 && $self->strand() == -1){
                $self->residues(revcomp(substr($contig->residues,
					       $nE, $self->length())));
                $self->inScope(1);
        }
	else {
		die "Invalid strand in __add_residues.";
	}

      done:
	return $self->inScope();
}

################################################ subroutine header begin ##

=head1 _load_feature

 Usage     : *private*

=for example begin

=for example end

=for example_testing

 Purpose   : Load up a hash of feature information.
 Returns   : A reference to a hash full of feature info.
 Argument  : A hash full of the feature's info.
 Throws    :
 Comments  : Generally just copies the hash fields,
           : but it creates FeatureLocation objects for anything in the
           : the location field.
 See Also  : CGL::Annotation::FeatureLocation

=cut

################################################## subroutine header end ##

sub _load_feature {
  my $hash = shift;

  my $feature = {};

  foreach my $k (keys %{$hash}){
    if ($k eq 'locations'){
      foreach my $l (@{$hash->{$k}}){
	push(@{$feature->{$k}}, new CGL::Annotation::FeatureLocation($l));
      }
    }
    else {
      $feature->{$k} = $hash->{$k};	
    }
  }
  return $feature;
}

################################################ subroutine header begin ##

=head1 _add_location

 Usage     : *private*

=for example begin

=for example end

=for example_testing

 Purpose   : Add a location to the list of locations.
 Returns   :
 Argument  : The location to be added.
 Throws    :
 Comments  : XXXX appears to be unused.
           :
 See Also  : CGL::Annotation::FeatureLocation

=cut

################################################## subroutine header end ##

sub _add_location {
        my $self = shift;
        my $l    = shift;

        push(@{$self->{locations}}, $l);
}

################################################ subroutine header begin ##

=head1 _add_relationship

 Usage     : *private*

=for example begin

=for example end

=for example_testing

 Purpose   : Add a relationship to the list of relationships.
 Returns   :
 Argument  : The relationship to be added.
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::FeatureRelationships

=cut

################################################## subroutine header end ##

sub _add_relationship {
	my $self = shift;
	my $r    = shift;

	push(@{$self->{relationships}}, $r);
}

################################################ subroutine header begin ##

=head1 _get_begin_end_on_src

 Usage     : *private*

 Purpose   : Get via metaPos the begin & end of a feature on it's source.
 Returns   : An list, first element is the begin, second element is the end
           : of the feature on the source.
 Argument  : The source contig. (XXXX contig?)
 Throws    :
 Comments  :
           :
 See Also  : CGL::Annotation::Feature::Intron::_get_begin_end_on_src

=cut

################################################## subroutine header end ##

sub _get_begin_end_on_src
 {
        my $self = shift;
        my $c    = shift;

        my $b_on_c = $self->metaPos($c, 0);
        my $e_on_c = $self->metaPos($c, $self->length);

        my $src_b = $c->location(0)->nbeg();
        my $src_e = $c->location(0)->nend();

        my $b_on_s;
        my $e_on_s;

        if ($src_b < $src_e){
                $b_on_s = $b_on_c + $src_b;
                $e_on_s = $e_on_c + $src_b;
        }
        else {
                $b_on_s = $src_b - $b_on_c;
                $e_on_s = $src_b - $e_on_c;
        }
        return ($b_on_s, $e_on_s);

}

################################################ subroutine header begin ##

=head1 _add_src_id

 Usage     : *private*

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub _add_src_id {
        my $self = shift;
        my $c    = shift;

        my ($src_b, $src_e) = $self->_get_begin_end_on_src($c);

        my $src_id = $c->location(0)->srcfeature_id() .
	  ":" . $src_b . ":" . $src_e;

        $self->src_id($src_id);
}

################################################ subroutine header begin ##

=head1 AUTOLOAD

 Usage     : *private*

 Purpose   : Implments a generic getter/setter routine for attributes that
           : aren't worth defining explicitly.
 Returns   : The current value of that field in the object.
 Argument  : The new value for that field in the object.
 Throws    :
 Comments  : Explicitly passing an argument of "undef" will effectively
           : undefine that field in the hash.
 See Also  :

=cut

################################################## subroutine header end ##

sub AUTOLOAD
 {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

	if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Annotation::Feature::AutoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

