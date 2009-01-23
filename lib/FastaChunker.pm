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
use POSIX qw(ceil);

@ISA = qw(
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;

        my $self = {};
        bless $self;

	$self->min_size(0);

	return $self;
}
#-------------------------------------------------------------------------------
sub load_chunks {
	my $self = shift;

	die " you must assign a parent_fasta!\n" 
	   unless defined($self->parent_fasta);

        die " you must assign a chunk_size!\n"
           unless defined($self->chunk_size);
	
	die " chunk_size must be greater than 0!\n"
           unless $self->chunk_size > 0;
	
	my $parent_def = Fasta::getDef($self->parent_fasta);
	my $parent_seq = Fasta::getSeq($self->parent_fasta);

	$self->parent_seq_length(length($parent_seq));
        my $fasta = '';

	my $t_c = ceil($self->parent_seq_length/$self->chunk_size);

	$self->total_chunks($t_c);
	my $l = $self->chunk_size();

	my $c = 0;
        for (my $i=0; $i< $self->parent_seq_length; $i+=$l){
		my $def = $parent_def. " CHUNK number:$c size:$l offset:$i";

		my $is_last = ($t_c - 1 == $c) ? 1 : 0;
		
		my $chunk = new FastaChunk();
		   $chunk->seq(substr($parent_seq, $i, $l));
		   $chunk->def($def);
		   $chunk->parent_def($parent_def);
		   $chunk->seqid(Fasta::def2SeqID($parent_def));
		   $chunk->size($l); #the max size of a chunk
		   $chunk->length(length($chunk->seq())); #the actual size of a chunk
		   $chunk->offset($i);
		   $chunk->number($c);
		   $chunk->is_last($is_last);
		   $chunk->parent_seq_length($self->parent_seq_length());

		if($chunk->length > $self->min_size || $chunk->number == 0){
		   $self->add_chunk($chunk);
		}
		else{
		   $self->last_chunk->seq($self->last_chunk->seq . $chunk->seq);
		   $self->last_chunk->length($self->last_chunk->length + $chunk->length);
		   $self->last_chunk->is_last(1);
		   $self->total_chunks($t_c - 1);
		}
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
sub last_chunk {
	my $self = shift;

	return $self->{chunks}->[-1];
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

        #print STDERR "FastaChunker::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------
1;


