#----------------------------------------------------------------------------
#----                            FastaChunker                            ---- 
#----------------------------------------------------------------------------
package FastaChunker;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use FastaChunk;
use Fasta;
use FastaFile;
@ISA = qw(
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;

        my $self = {};
        bless $self;

	return $self;
}
#-------------------------------------------------------------------------------
sub load_chunks {
	my $self = shift;

	die " you must assign a parent_fasta!\n" 
	unless defined($self->parent_fasta);

        die " you must assign a chunk_size!\n"
        unless defined($self->chunk_size);

	my $parent_def = Fasta::getDef($self->parent_fasta);
	my $parent_seq = Fasta::getSeq($self->parent_fasta);

	$self->parent_seq_length(length($$parent_seq));
        my $fasta = '';

	
	my $l = $self->chunk_size();

	my $c = 0;
        for (my $i=0; $i< length($$parent_seq);$i+=$l){
		my $def = $parent_def. " CHUNK number:$c size:$l offset:$i";

		my $chunk = new FastaChunk();
		   $chunk->seq(substr($$parent_seq, $i, $l));
		   $chunk->def($def);
		   $chunk->size($l);
		   $chunk->offset($i);
		   $chunk->number($c);

		$self->add_chunk($chunk);
		$c++;
        }

}
#-------------------------------------------------------------------------------
sub get_chunk {
	my $self = shift;
	my $i    = shift;

	return $self->{chunks}->[$i];
}
#-------------------------------------------------------------------------------
sub add_chunk {
	my $self  = shift;
	my $chunk = shift;

	push(@{$self->{chunks}}, $chunk);
}
#-------------------------------------------------------------------------------
#---------------------------  CLASS FUNCTIONS ----------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "FastaChunker::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------
1;


