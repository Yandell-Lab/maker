###################################################### main header begin ##

=head1 CGL::Ontology::Parser::OBO

CGL::Ontology::Parser::OBO - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example begin

=for example end

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

package CGL::Ontology::Parser::OBO;

use strict;
use warnings;

use CGL::Ontology::Node;
use CGL::Ontology::Trace;
use CGL::Ontology::NodeRelationship;
#use CGL::Ontology::Parser;
use FileHandle;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
#  @ISA = qw(
#	    CGL::Ontology::Parser
#	   );
  @ISA = qw(
	   );
}

###########################################################################
#
# Once upon a time, this parser worked by setting $/ to be '[Term]'
# and then reading in chunks separated by occurences of $/.  At the
# time, the obo file format put a bunch of information at the top of
# the file, including any [Typedef] entries.  The parse routine above
# skipped over the entire mess in one swell foop by checking if the
# first line of the stanza was the format-version.  The rest of the
# entries were known to be [Term]'s and the parser worked fabulously.
#
# At some point, the obo file format changed, and [Typedef]'s were
# sprinkled throughout the file.
#
# This routine reads through the file and returns stanzas that are
# functionally identical what you'd get slurping up the old format
# with $/ = '[Term]'.  It achieves that goal by noticing if it's about
# to slurp up a [Typedef] (because the last pass through left the
# [...]  of the next thing to be parsed in $$pushback_ref) and just
# skipping over it.
#
# For a more thorough parsing of an obo file, consider: the go-perl
# portion of the go-dev distribution that can be found in the
# CVS repository at:
#     http://sourceforge.net/projects/geneontology
#
###########################################################################


################################################ subroutine header begin ##

=head2 new

 Usage     :

 Purpose   : Guarantee that an ...::OBO isn't ever directly instantiated.
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub new {
  my ($class, $file) = @_;
  my $self = bless ({}, $class);

  if ($file) {
    $self->file($file);
  }

  return ($self);
}

################################################ subroutine header begin ##

=head2 parse

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  : See tests in a subclass.
           :
 See Also  : CGL::Ontology::Parser::OBO::SO

=cut

################################################## subroutine header end ##

sub parse {
	my ($self, $file) = @_;

	my $fh = new FileHandle();
	
	$file = $self->file() unless (defined($file));

	$fh->open($file);

	my %nodes;
	my $pushback = undef;

	while (my $stanza = _get_term($fh, \$pushback)) {

		my @lines = split(/\n/, $stanza);

		next if $lines[0] =~ /format-version/;

		next unless @lines;

		my $term = _load_lines(\@lines);

		my $node =_term_to_node($term);

		my $isa_relationships = _term_to_isa_feature_relationships($term);
	
		push(@{$node->{relationships}}, @{$isa_relationships}) if @{$isa_relationships};

		my $part_of_relationships = _term_to_part_of_feature_relationships($term);

		push(@{$node->{relationships}}, @{$part_of_relationships}) if @{$part_of_relationships};


		my $f = new CGL::Ontology::Node($node);
		$nodes{$f->id} = $f;
	}
	$fh->close;

	_reverse_relationships(\%nodes);

	return \%nodes;
}

################################################ subroutine header begin ##

=head2 _get_term

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _get_term {
  my($fh, $pushback_ref) = @_;
  my($term) = undef;

  # skip over [Typedef] stanzas.
  if ($$pushback_ref && $$pushback_ref =~ m|\[Typedef\]|) {
    while(<$fh>) {
      last if (/\[Term\]/);
    }
  }

  # protect against running into the bottom of the file....
  return($term) if eof($fh);

  while (<$fh>) {
    $term .= $_;
    if ((/\[Term\]/) or (/\[Typedef\]/)) {
      $$pushback_ref = $_;
      return($term);
    }
  }

  return($term);
}

################################################ subroutine header begin ##

=head2 _load_lines

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _load_lines {
	my $lines = shift;

	my %data;
	foreach my $line (@{$lines}){
		next unless $line =~ /\S+/;
		next if     $line =~ /\[Term\]/;
		next if     $line =~ /\[Typedef\]/;

		$line =~ s|[^!]!.*$||;
		chomp $line;
		my ($key, $value) = $line =~ /^(\S+)\:\s+(.*)/;

		# XXXX die???
		die unless defined($key) && defined($value);

		$value =~ s/\"//g  if $key eq 'synonym';
		$value =~ s/\[//g  if $key eq 'synonym';
		$value =~ s/\]//g  if $key eq 'synonym';
		$value =~ s/\s+$// if $key eq 'synonym';

		push(@{$data{$key}}, $value);
	}

	return \%data;
}

################################################ subroutine header begin ##

=head2 _term_to_node

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _term_to_node {
	my $term = shift;


	my %node;

	$node{id}       = $term->{id}->[0];
	$node{name}     = $term->{name}->[0];
        $node{synonyms} = $term->{synonym};
	$node{def}      = $term->{def}->[0];

	return \%node;
}

################################################ subroutine header begin ##

=head2 _term_to_isa_feature_relationships

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _term_to_isa_feature_relationships {
        my $term = shift;

	my $s_id = $term->{id}->[0];

	my %hash;
	my @f_rs;
	foreach my $o_id (@{$term->{is_a}}){

		$hash{subject_id} = $s_id;
		$hash{object_id}  = $o_id;
		$hash{type}       = 'isa';
	
		my $f = new CGL::Ontology::NodeRelationship(\%hash);

		push(@f_rs, $f);
	}

	return \@f_rs;
}

################################################ subroutine header begin ##

=head2 _term_to_part_of_feature_relationships

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

sub _term_to_part_of_feature_relationships {
        my $term = shift;

        my $s_id = $term->{id}->[0];

        my %hash;
        my @f_rs;
        foreach my $str (@{$term->{relationship}}){

		my @data = split(/\s+/, $str);

                $hash{subject_id} = $s_id;
                $hash{object_id}  = $data[1];
                $hash{type}       = $data[0];

                my $f = new CGL::Ontology::NodeRelationship(\%hash);

                push(@f_rs, $f);
        }

        return \@f_rs;
}

################################################ subroutine header begin ##

=head2 _reverse_relationships

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub _reverse_relationships {
	my $nodes = shift;

	my %r_hash;
        foreach my $id (keys %{$nodes}){
		my $node = $nodes->{$id};
		foreach my $r ($node->relationships){

			push(@{$r_hash{$r->oF}}, $r);

		}
	}


        foreach my $id (keys %{$nodes}){
		my $node = $nodes->{$id};
		foreach my $r (@{$r_hash{$id}}){
			$node->_add_relationship($r);
		}
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

        if (@_){
                $self->{file} = shift;
        }

	return $self->{file};
}

################################################ subroutine header begin ##

=head2 AUTOLOAD

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
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
	    print STDERR "CGL::Ontology::Parser::OBO::AutoLoader called for: ",
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
