#-------------------------------------------------------------------------------
#------               Chaos::Feature::Sequence_variant                 ---------
#-------------------------------------------------------------------------------
package CGL::Annotation::Feature::Sequence_variant;

use strict;
use warnings;

use CGL::Annotation::Feature;
use CGL::Clone qw(clone);

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
            CGL::Annotation::Feature
           );
}
#-------------------------------------------------------------------------------
#------------------------------- METHODS ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
        my $class    = shift;
        my $feature  = shift;


        my $self = clone($feature);

        bless($self, $class);

        return $self;
}
#-------------------------------------------------------------------------------
sub pos_on_translation {
        my $self = shift;
        my $s    = shift;
        my $t    = shift;
        my $p    = shift;
        my $pos  = shift;

        return undef  unless $self->_is_in_transcript($s, $t);

        my $pos_on_s = $self->metaPos($s, $pos);

        return undef  unless defined($pos_on_s);

        my $pos_on_t = $s->metaPos($t, $pos_on_s);

        return undef unless defined($pos_on_t);

        my $pos_on_p = $t->metaPos($p, $pos_on_t);

        return $pos_on_p;
}
#-------------------------------------------------------------------------------
sub get_property {
        my $self     = shift;
        my $property = shift;

        return $self->{properties}->{$property};
}
#-------------------------------------------------------------------------------
sub properties {
        my $self     = shift;

        return $self->{properties};
}
#-------------------------------------------------------------------------------
sub pos_on_transcript {
        my $self = shift;
        my $s    = shift;
        my $t    = shift;
        my $pos  = shift;

        my $pos_on_s = $self->metaPos($s, $pos);

        return undef unless defined($pos_on_s);

        my $pos_on_t = $s->metaPos($t, $pos_on_s);

        return $pos_on_t;
}
#-------------------------------------------------------------------------------
sub length {
        my $self = shift;

        return abs($self->nbeg - $self->nend);

}
#-------------------------------------------------------------------------------
sub metaPos {
        my $self          = shift;
        my $what_f        = shift;
        my $where_in_self = shift;
        my $where_in_feature;
        if ($what_f->type eq 'contig'){
                if ($self->strand() == 1){
                        return $self->nbeg + $where_in_self;
                }
                else {
                        return $self->nend - $where_in_self;
                }
        }
        else {
                die "this metaPos not yet supported!\n";
        }

        return $where_in_feature;
}
#-------------------------------------------------------------------------------
#------------------------------- PRIVATE ---------------------------------------
#-------------------------------------------------------------------------------
sub _is_in_transcript {
        my $self = shift;
        my $s    = shift;
        my $t    = shift;

        my $pos_self_b_on_s = $self->metaPos($s, 0);

        return 0 unless defined($pos_self_b_on_s);

        my $pos_t_b_on_s = $t->nbeg();
        my $pos_t_e_on_s = $t->nend();

        ($pos_t_b_on_s, $pos_t_e_on_s) = ($pos_t_e_on_s, $pos_t_b_on_s)
        if $t->strand == -1;

        return 0 if ($pos_self_b_on_s < $pos_t_b_on_s
                 ||  $pos_self_b_on_s > $pos_t_e_on_s);

        return 0 if $self->_is_in_intron($s, $t);

        return 1;
}
#-------------------------------------------------------------------------------
sub _is_in_intron {
        my $self = shift;
        my $s    = shift;
        my $t    = shift;

        my $pos_v_b_on_s = $self->metaPos($s, 0);

        my $i = 0;
        while(my $intron = $t->intron($i)){
                $i++;
                my $i_start = $intron->nbeg();
                my $i_end   = $intron->nend();
                ($i_start, $i_end) = ($i_end, $i_start) if $intron->strand == -1;
                return 1 if ($pos_v_b_on_s >= $i_start && $pos_v_b_on_s <= $i_end);
        }

        return 0;
}
#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        if($ENV{CGL_CHATTER}) {
	    print STDERR "Annotation::Feature::Sequence_variant AutoLoader called for: ",
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

1;

