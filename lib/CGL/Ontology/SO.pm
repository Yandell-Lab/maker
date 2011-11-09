###################################################### main header begin ##

=head1 CGL::Ontology::SO

CGL::Ontology::SO - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

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

package CGL::Ontology::SO;

use strict;
use warnings;

use CGL::Ontology::Ontology;
use CGL::Ontology::Parser::OBO;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA         = qw (
		     CGL::Ontology::Ontology
		    );
}

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example begin

  use CGL::Ontology::SO;
  my $so = new CGL::Ontology::SO("sample_data/so.new.obo");

=for example end

=for example_testing
  isa_ok($so, "CGL::Ontology::SO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub new
 {
	my $class = shift;
	my $source  = shift;

	my $self = {};

	bless($self, $class);

	if(defined($source)){
		$self->source($source);
	}
	elsif (defined($ENV{'CGL_SO_SOURCE'})){
		$self->source($ENV{'CGL_SO_SOURCE'});
	}
	else {
	    (my $file = $INC{'CGL/Ontology/SO.pm'}) =~ s/Ontology\/SO\.pm/so.obo/;
	    if(-f $file){
		$self->source($file);
	    }
	    else{
		print "\n#######################################################\n";
		print "Error in CGL::Ontology::SO.pm.\n";
		print "You must either set the CGL_SO_SOURCE environment varible\n";
		print "or provide a so.obo file.\n";
		print "#######################################################\n\n";
		die;
	    }
	}

	# might want to instantiate a variety of parsers, depending on
	# the source....
	$self->parser(new CGL::Ontology::Parser::OBO($self->source));

	$self->nodes($self->parser->parse());

	$self->_index_on_node_names();

	return $self;
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
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "CGL::Ontology::SO::AutoLoader called for: ",
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

