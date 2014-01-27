###################################################### main header begin ##

=head1 NAME

CGL::TranslationMachine - An enhancement of BioPerl's CodonTable.

=head1 SYNOPSIS

=for example
  use UNIVERSAL qw( isa );

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;

=for example end

=for example_testing
  isa_ok($t, "CGL::TranslationMachine", "Check if it's the right type.");

=head1 DESCRIPTION

An extension of the Bio::Tools::CodonTable class, adding some methods
that the CGL libraries find useful.

=head1 USAGE

See the method documentation below.

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

Bio::Tools::CodonTable

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package CGL::TranslationMachine;

use strict;
use warnings;

use Bio::Tools::CodonTable;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;		# XXXX svn file version?
  @ISA         = qw (Bio::Tools::CodonTable);
}

use constant FRAMES => 3;

################################################ subroutine header begin ##

=head2 new

 Usage     :

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;

=for example end

=for example_testing
  isa_ok($t, "CGL::TranslationMachine", "Check if it's the right type.");

 Purpose   : Create a new CGL::TranslationMachine object.
 Returns   : Returns the new CGL::TranslationMachine object.
 Arguments : See Bio::Tools::CodonTable discussion of codon table
             identifiers.
 Throws    : None, but Bioperl\'s CodonTable may throw errors.
 Comments  :
           :
 See Also  : Bio::Tools::CodonTable

=cut

################################################ subroutine header end   ##

sub new {
        my $class  = shift;
        my @args   = @_;

        my $self = $class->SUPER::new(@args);

        bless ($self, $class);

	#make default codon table as M only start codon
	my $id = $self->add_table(
	    'Strict',
	    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
	    '-----------------------------------M----------------------------');
	$self->id($id);

        return $self;
}

################################################ subroutine header begin ##

=head2 name

 Usage     :

=for example
  $seq = "atgaaaaaauaa";

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;
  $t->name('Strict');
  my ($name) = $t->name();

=for example end

=for example_testing
  is($name, "Strict", "Did it find the return the codon table?");

 Purpose   : Get/set the codon table by name
 Returns   : Name of curret codon table
 Arguments : Codon talbe name
 Throws    :
 Comments  :
           :
 See Also  : Bio::Tools::CodonTable

=cut

################################################## subroutine header end ##
sub name{
    my ($self, $arg) = @_;

    if($arg){
	my %index;
	@index{@Bio::Tools::CodonTable::NAMES} =  (1..@Bio::Tools::CodonTable::NAMES);

	my $id = $index{$arg};

	die "ERROR: No codon table named '$arg' in CGL::TranaslationMachine::name\n"
	    if(!defined($id));

	$self->id($id);
    }

    my ($id) = $self->{'id'};
    return $Bio::Tools::CodonTable::NAMES[$id-1];
}
################################################ subroutine header begin ##

=head2 longest_translation

 Usage     :

=for example
  $seq = "atgaaaaaauaa";

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;
  my ($longest_orf,$offset) = $t->longest_translation($seq);

=for example end

=for example_testing
  is($longest_orf, "MKK", "Did it find the right orf?");
  is($offset, 0, "Did it get the offset right?");

 Purpose   : Find the longest open reading frame in any of three
             frames in the input sequence.
 Returns   : A list of the amino acid sequence of the longest
             open reading frame and the offset at which it begins in
             the input sequence.
 Arguments : A translatable sequence.
 Throws    :
 Comments  :
           :
 See Also  : Bio::Tools::CodonTable

=cut

################################################## subroutine header end ##

sub longest_translation {
    my $self = shift;
    my $seq  = shift;
    my ($longest,$longest_len,$offset,@frames) = ('',0,0);
    for( my $i = 0; $i < FRAMES; $i++ ) {
	my $counter = $i;
	foreach my $orf ( split(/\*/,$self->translate(substr($seq,$i))) ) {
	    my $orflen = length($orf);
	    if( $orflen > $longest_len ) {
		$longest = $orf;
		$offset  = $counter;
		$longest_len = $orflen;
            }
	    # multi x 3 for aa->codon
	    # 1 for the stop codon not included
	    # in length	
	    $counter += ($orflen + 1)*3;
	}
    }

    return ($longest,$offset);
}

################################################ subroutine header begin ##

=head2 longest_translation

 Usage     :

=for example
  $seq = "atgaaaaaauaa";

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;
  my ($longest_orf,$offset) = $t->longest_translation($seq);

=for example end

=for example_testing
  is($longest_orf, "MKK", "Did it find the right orf?");
  is($offset, 0, "Did it get the offset right?");

 Purpose   : Find the longest open reading frame in any of three
             frames in the input sequence.
 Returns   : A list of the amino acid sequence of the longest
             open reading frame and the offset at which it begins in
             the input sequence.
 Arguments : A translatable sequence.
 Throws    :
 Comments  :
           :
 See Also  : Bio::Tools::CodonTable

=cut

################################################## subroutine header end ##

sub longest_translation_plus_stop {
    my $self = shift;
    my $seq  = shift;
    my $seqlength = length($seq);

    my ($longest,$longest_len,$offset,@frames) = ('',0,0);
    for( my $i = 0; $i < FRAMES; $i++ ) {
	my $counter = $i;
	foreach my $orf ( split(/\*/,$self->translate(substr($seq,$i))) ) {
	    my $orflen = length($orf);
	    if( $orflen > $longest_len ) {
		$longest = $orf;
		$offset  = $counter;
		$longest_len = $orflen;
		my $c_start = ($orflen * 3) + $offset;
		my $codon = ($c_start < $seqlength) ? substr($seq, $c_start, 3) : '';
		if(length($codon) == 3 && $self->translate($codon) eq '*'){
		    $longest .= '*';
		}
            }
	    # multi by 3 for aa->codon
	    # 1 for the stop codon (included
	    # in seq but not in length)	
	    $counter += ($orflen + 1)*3;
	}
    }

    return ($longest,$offset);
}

################################################ subroutine header begin ##

=head2 translate_from_offset

 Usage     :

=for example
  $seq = "atgatgatgaaaaaauaa";

=for example begin

  use CGL::TranslationMachine;
  my $t = new CGL::TranslationMachine;
  my $aaseq = $t->translate_from_offset($seq, 6);

=for example end

=for example_testing
  is($aaseq, "MKK*", "Did it find the right orf?");

 Purpose   : Translate a sequence, beginning at an offset from the
             beginning of that sequence.
 Returns   : A string that contains the translated amino acid sequence,
             possibly including the stop codon.
 Arguments : An input sequence and an offset (starting at 0).
 Throws    : none, although it uses routines from Bio::Tools::CodonTable.
 Comments  :
           :
 See Also  : Bio::Tools::CodonTable::translate

=cut

################################################ subroutine header end   ##

sub translate_from_offset {
	my $self = shift;	
	my $seq  = shift;
	my $i    = shift;

	$seq = substr($seq, $i);

	return $self->translate($seq);
}

################################################ subroutine header begin ##

################################################ subroutine header end   ##

# XXXX autoloader to deal with.

sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "TranslationMachine:AutoLoader called for: ",
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

1;
__END__;
