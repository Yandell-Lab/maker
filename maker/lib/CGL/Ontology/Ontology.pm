###################################################### main header begin ##

=head1 CGL::Ontology::GO

CGL::Ontology::GO - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example begin

  use CGL::Ontology::SO;
  my $so = new CGL::Ontology::SO("sample_data/so.new.obo");

=for example end

=for example_testing
  isa_ok($so, "CGL::Ontology::SO", "Check if it's the right type.");

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

package CGL::Ontology::Ontology;

use strict;
use warnings;

use CGL::Ontology::Node;
use CGL::Ontology::Trace;
use CGL::Ontology::NodeRelationship;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw ( );	# XXXX superclasses go here.
}

################################################ subroutine header begin ##

=head2 new

 Usage     : *This needs to be implemented by a subclass*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub new {

  die "Should be implemented in subclass.";
}

################################################ subroutine header begin ##

=head2 nodes

 Usage     : How to use this function/method

=for example begin

  use CGL::Ontology::SO;

  # for testing purposes
  $ENV{CGL_SO_SOURCE} = "sample_data/so.new.obo";
  my $so = new CGL::Ontology::SO;

  my $nodes1 = $so->nodes();
  my $nodes2 = $so->nodes(undef);
  my $nodes3 = $so->nodes($nodes1);

=for example end

=for example_testing
  isa_ok($nodes1, "HASH", "Check the getter.");
  is($nodes2, undef, "Check clearing.");
  isa_ok($nodes3, "HASH", "Check the getter.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nodes
 {
	my $self  = shift;

	if (@_){
		$self->{nodes} = shift;
	}

	return $self->{nodes};
}

################################################ subroutine header begin ##

=head2 find_depth

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub find_depth
 {
	my $self  = shift;
	my $a = shift;
	my $a_node = $self->node_by_name($a);
	return undef unless (defined $a_node);
	my $depth = $self->trace($a_node, 'parents')->depth;
	return $depth;
}

################################################ subroutine header begin ##

=head2 find_depth_to_leaf

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub find_depth_to_leaf
{
	my $self  = shift;
	my $a = shift;
	my $a_node = $self->node_by_name($a);
	return undef unless (defined $a_node);
	my $depth = $self->trace($a_node, 'children')->depth;
	#my @terminii = $self->trace($a_node, 'children')->terminii;
	return $depth;
}

################################################ subroutine header begin ##

=head2 name_id_index

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub name_id_index
 {
        my $self  = shift;
        my $names = shift;

        if (defined($names)){
                $self->{name_id_index} = $names;
        }
        else {
                return $self->{name_id_index};
        }

}

################################################ subroutine header begin ##

=head2 file

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub file
 {
        my $self  = shift;
        my $file = shift;

        if (defined($file)){
                $self->{file} = $file;
        }
        else {
                return $self->{file};
        }

}

################################################ subroutine header begin ##

=head2 all_b_have_a

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub all_b_have_a
 {
	my $self = shift;
	my $b    = shift;
	my $a    = shift;

	my $a_node = $self->node_by_name($a);
	my $b_node = $self->node_by_name($b);

        foreach my $part ($self->trace($b_node, 'parts')->features(0)){
                return 1 if $part->name eq $a_node->name();
        }


        my $p_trace_obj = $self->trace($b_node, 'parents');

        foreach my $parents ($p_trace_obj->features){
                foreach my $parent (@{$parents}){
                        foreach my $part ($self->trace($parent, 'parts')->features(0)){
                                return 2 if $part->name eq $a_node->name();
                        }
                }
        }
	return 0;

}

################################################ subroutine header begin ##

=head2 some_b_have_a

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub some_b_have_a
 {
        my $self = shift;
        my $b    = shift;
        my $a    = shift;

        my $a_node = $self->node_by_name($a);
        my $b_node = $self->node_by_name($b);

	return 1 if $self->all_b_have_a($b, $a);

        foreach my $part ($self->trace($b_node, 'parts')->features(0)){
                return 2 if $self->a_is_hyponym_of_b($a_node->name(), $part->name);
        }

        my $c_trace_obj = $self->trace($b_node, 'children');

        foreach my $children ($c_trace_obj->features){
                foreach my $child (@{$children}){
                        foreach my $part ($self->trace($child, 'parts')->features(0)){
                                return 3 if $part->name eq $a_node->name();
				return 4 if $self->a_is_hyponym_of_b($a_node->name(), $part->name);
                        }
                }
        }

	return 0;
}

################################################ subroutine header begin ##

=head2 some_b_could_have_a

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub some_b_could_have_a
 {
        my $self = shift;
        my $b    = shift;
        my $a    = shift;

        my $a_node = $self->node_by_name($a);
        my $b_node = $self->node_by_name($b);

	return 1 if $self->some_b_have_a($b, $a);

        my $p_trace_obj = $self->trace($b_node, 'parents');

        foreach my $parent ($p_trace_obj->features){
                foreach my $parent (@{$parent}){
                        foreach my $part ($self->trace($parent, 'parts')->features(0)){
                                return 3 if $part->name eq $a_node->name();
                                return 4 if $self->a_is_hyponym_of_b($a_node->name(), $part->name);
                        }
                }
        }

        return 0;
}

################################################ subroutine header begin ##

=head2 a_is_an_optional_part_of_b

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_an_optional_part_of_b
 {
        my $self = shift;
        my $a    = shift;
        my $b    = shift;

        my $a_node = $self->node_by_name($a);
        my $b_node = $self->node_by_name($b);

        foreach my $part ($self->trace($b_node, 'parts')->features(0)){
                return 0 if $part->name eq $a_node->name();
                return 1 if $self->a_is_hyponym_of_b($a_node->name(), $part->name);
        }
        return 0;

}

################################################ subroutine header begin ##

=head2 a_is_an_essential_part_of_some_b

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_an_essential_part_of_some_b
 {
        my $self = shift;
        my $a    = shift;
        my $b    = shift;

        my $a_node = $self->node_by_name($a);
        my $b_node = $self->node_by_name($b);

	return 0 if $self->a_is_an_essential_part_of_b($a, $b);

	return 0 if $self->a_is_an_optional_part_of_b($a, $b);

        my $p_trace_obj = $self->trace($b_node, 'parents');

        foreach my $parents ($p_trace_obj->features){
                foreach my $parent (@{$parents}){
                        foreach my $part ($self->trace($parent, 'parts')->features(0)){
                                return 1 if $part->name eq $a_node->name();
                        }
                }
        }
        return 0;

}

################################################ subroutine header begin ##

=head2 a_is_meronym_of_b_old

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_meronym_of_b_old
 {
	my $self = shift;
	my $a    = shift;
	my $b    = shift;

	my $a_f = $self->node_by_name($a);

	my $b_f = $self->node_by_name($b);



	#print "b_f:".$b_f->name."\n";
        foreach my $part ($self->trace($b_f, 'parts')->features(0)){
                #print "  part:".$part->name."\n";
		return 1 if $part->name eq $a_f->name();
        }

	my $p_trace_obj = $self->trace($b_f, 'parents');

	foreach my $parents ($p_trace_obj->features){
		foreach my $parent (@{$parents}){
			#print "parent:".$parent->name."\n";
			foreach my $part ($self->trace($parent, 'parts')->features(0)){
				#print "  part:".$part->name."\n";
				return 2 if $part->name eq $a_f->name();
			}
		}	
	}

        my $c_trace_obj = $self->trace($b_f, 'children');

        foreach my $children ($c_trace_obj->features){
                foreach my $child (@{$children}){
                        #print "child:".$child->name."\n";
                        foreach my $part ($self->trace($child, 'parts')->features(0)){
                                #print "  part:".$part->name."\n";
                                return 3 if $part->name eq $a_f->name();
                        }
                }
        }

	return 0;
}

################################################ subroutine header begin ##

=head2 a_is_meronym_of_b

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_meronym_of_b
 {
        my $self = shift;
        my $a    = shift;
        my $b    = shift;

        my $a_f = $self->node_by_name($a);

        my $b_f = $self->node_by_name($b);



        #print "b_f:".$b_f->name."\n";
        foreach my $part ($self->trace($b_f, 'parts')->features(0)){
                #print "  part:".$part->name."\n";
                return 1 if $part->name eq $a_f->name();
		return 2 if $self->a_is_hyponym_of_b($a_f->name, $part->name);
        }

        my $p_trace_obj = $self->trace($b_f, 'parents');

        foreach my $parents ($p_trace_obj->features){
                foreach my $parent (@{$parents}){
                        #print "parent:".$parent->name."\n";
                        foreach my $part ($self->trace($parent, 'parts')->features(0)){
                                #print "  part:".$part->name."\n";
				return 2 if $part->name eq $a_f->name()
				          || $self->a_is_hyponym_of_b($a_f->name(), $part->name);
                        }
                }
        }

        my $c_trace_obj = $self->trace($b_f, 'children');

        foreach my $children ($c_trace_obj->features){
                foreach my $child (@{$children}){
                        #print "child:".$child->name."\n";
                        foreach my $part ($self->trace($child, 'parts')->features(0)){
                                #print "  part:".$part->name."\n";
                                return 3 if $part->name eq $a_f->name()
                                         || $self->a_is_hyponym_of_b($a_f->name, $part->name);
                        }
                }
        }

        return 0;
}

################################################ subroutine header begin ##

=head2 a_is_hypomeronym_of_b

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_hypomeronym_of_b
 {
        my $self = shift;
        my $a    = shift;
        my $b    = shift;

        my $a_f = $self->node_by_name($a);

        my $b_f = $self->node_by_name($b);

	return 0 unless ($a_f and $b_f); # handle nodes not in Ontology.

        return 1 if $b_f->name() eq $a_f->name();
        my $p_trace_obj = $self->trace($a_f, 'wholes');

        #print "a_f:".$a_f->name."\n";
        foreach my $hypermeronyms ($p_trace_obj->features){
		foreach my $hypermeronym (@{$hypermeronyms}){
        		#print "hypermeronym:".$hypermeronym->name."\n";
               		return 2 if $hypermeronym->name eq $b_f->name();
		}
        }

        return 0;
}

################################################ subroutine header begin ##

=head2 a_is_hyponym_of_b

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub a_is_hyponym_of_b
 {
        my $self = shift;
        my $a    = shift;
        my $b    = shift;

        my $a_f = $self->node_by_name($a);

        my $b_f = $self->node_by_name($b);

	return 0 unless ($a_f and $b_f); # handle nodes not in Ontology.

	return 1 if $b_f->name() eq $a_f->name();
        my $p_trace_obj = $self->trace($a_f, 'parents');

	#print "a_f:".$a_f->name."\n";
        foreach my $parents ($p_trace_obj->features){
                foreach my $parent (@{$parents}){
                        #print "parent:".$parent->name."\n";
                        return 2 if $parent->name eq $b_f->name();
                }
        }

        return 0;
}

################################################ subroutine header begin ##

=head2 hypernym

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub hypernym
 {
	my $self = shift;
	my $name = shift;

        my $n = $self->node_by_name($name);

        my @hypernyms = $self->trace($n, 'parents')->features(0);

        my @names;
        foreach my $hypernym (@hypernyms){
                push(@names, $hypernym->name());
        }

        return join(' & ', @names) || 'ROOT';

}

################################################ subroutine header begin ##

=head2 hyponym

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub hyponym
 {
        my $self = shift;
        my $name = shift;

        my $n = $self->node_by_name($name);

	my @hyponyms = $self->trace($n, 'children')->features(0);

	my @names;
        foreach my $hyponym (@hyponyms){
                push(@names, $hyponym->name());
        }

        return join(' & ', @names) || 'terminal node';
}

################################################ subroutine header begin ##

=head2 hypermeronym

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub hypermeronym
 {
        my $self = shift;
        my $name = shift;

        my $n = $self->node_by_name($name);

        my @hypermeronyms = $self->trace($n, 'wholes')->features(0);

        my @names;
        foreach my $hypermeronym (@hypermeronyms){
                push(@names, $hypermeronym->name());
        }

        return join(' & ', @names) || 'ROOT';

}

################################################ subroutine header begin ##

=head2 hypomeronym

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub hypomeronym
 {
        my $self = shift;
        my $name = shift;

        my $n = $self->node_by_name($name);

        my @hypomeronyms = $self->trace($n, 'parts')->features(0);

        my @names;
        foreach my $hypomeronym (@hypomeronyms){
                push(@names, $hypomeronym->name());
        }

        return join(' & ', @names) || 'terminal node';
}

################################################ subroutine header begin ##

=head2 node_by_name

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub node_by_name
 {
	my $self = shift;
	my $name = shift;

	my $ids = $self->name_id_index->{$name};
	
	return  $self->node_by_id($ids->[0]);
	
}

################################################ subroutine header begin ##

=head2 node_by_id

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub node_by_id
 {
	my $self = shift;
	my $id   = shift;

	return $self->nodes->{$id};	
}

################################################ subroutine header begin ##

=head2 relationships

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub relationships
 {
        my $self = shift;

        return @{$self->{relationships}|| []};
}

################################################ subroutine header begin ##

=head2 trace

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub trace
 {
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

        $trace = new CGL::Ontology::Trace unless defined($trace);


        foreach my $r ($f->relationships){
		if    ($type eq 'parts'){
			next unless ($r->logus eq 'part_of' ||
				     $r->logus eq 'member_of');
			next unless $r->oF eq $f->id;
			$trace->add_feature($self->node_by_id($r->sF), $i);
                       	$self->trace($self->node_by_id($r->sF), $type, $i, $trace);
		}
		elsif ($type eq 'producers'){
			next unless ($r->logus eq 'produced_by' ||
				     $r->logus eq 'derives_from');
			next unless $r->sF eq $f->id;
			$trace->add_feature($self->node_by_id($r->oF), $i);
                       	$self->trace($self->node_by_id($r->oF), $type, $i, $trace);
		}
		elsif ($type eq 'wholes'){
			next unless ($r->logus eq 'part_of' ||
				     $r->logus eq 'member_of');
			next unless $r->sF eq $f->id;	
			$trace->add_feature($self->node_by_id($r->oF), $i);
                        $self->trace($self->node_by_id($r->oF), $type, $i, $trace);
		}
               elsif ($type eq 'produces'){
                        next unless ($r->logus eq 'produced_by' ||
				     $r->logus eq 'derives_from');
                        next unless $r->oF eq $f->id;
                        $trace->add_feature($self->node_by_id($r->sF), $i);
                        $self->trace($self->node_by_id($r->sF), $type, $i, $trace);
                }

              elsif ($type eq 'children'){
                        next unless $r->logus eq 'isa';
                        next unless $r->oF eq $f->id;
                        $trace->add_feature($self->node_by_id($r->sF), $i);
                        $self->trace($self->node_by_id($r->sF), $type, $i, $trace);
                }
              elsif ($type eq 'parents'){
                        next unless $r->logus eq 'isa';
                        next unless $r->sF eq $f->id;
                        $trace->add_feature($self->node_by_id($r->oF), $i);
                        $self->trace($self->node_by_id($r->oF), $type, $i, $trace);
                }
		else {
			die "unknown type($type) in CGL::Ontology::trace!\n";
		}
        }

       	return $trace;

}
#-------------------------------------------------------------------------------
#--------------------------------- PRIVATE -------------------------------------
#-------------------------------------------------------------------------------

################################################ subroutine header begin ##

=head2 _index_on_node_names

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub _index_on_node_names
 {
	my $self = shift;

	my %names;
	foreach my $id (keys %{$self->nodes}){
		my $node = $self->node_by_id($id);

		my $name = $node->name();
		push(@{$names{$name}}, $id);
		foreach my $syn (@{$node->synonyms}){
			push(@{$names{$syn}}, $id);
		}
	}

	$self->name_id_index(\%names);
}

################################################ subroutine header begin ##

=head2 source

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub source
 {
        my $self  = shift;

        if (@_){
                $self->{source} = shift;
        }

	return $self->{source};
}


################################################ subroutine header begin ##

=head2 parser

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub parser
 {
        my $self  = shift;

        if (@_){
                $self->{parser} = shift;
        }

	return $self->{parser};
}


################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub AUTOLOAD
 {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "SO::AutoLoader called for: ",
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

1; #this line is important and will help the module return a true value
__END__
