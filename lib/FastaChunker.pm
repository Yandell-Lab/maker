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
        bless($self, $class);

	$self->min_size(0);

	return $self;
}
#------------------------------------------------------------------------------
sub parent_seq {
    my $self = shift;
    my $seq = shift;

    if(defined $seq){
	#always work with references
	$seq = $$seq while(ref($seq) eq 'REF');
	my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;
	
	$self->{parent_seq} = $seq_ref;
    }

    return $self->{parent_seq}; #this can be a sequence objection
}
#------------------------------------------------------------------------------
sub parent_fasta {
    my $self = shift;
    my $fasta = shift;

    if(defined $fasta){
	#always work with references
	$fasta = $$fasta while(ref($fasta) eq 'REF');
	my $fasta_ref = (ref($fasta) eq '') ? \$fasta : $fasta;

	$self->{parent_def} = Fasta::getDef($fasta_ref);
	$self->{parent_seq} = Fasta::getSeqRef($fasta_ref);
    }
}
#------------------------------------------------------------------------------
sub load_chunks {
	my $self = shift;

	die " you must assign a parent_seq or parent_fasta!\n" 
	    unless(defined($self->parent_seq));

	die " you must assign a parent_def or parent_fasta!\n" 
	    unless(defined($self->parent_def));

        die " you must assign a chunk_size!\n"
           unless defined($self->chunk_size);
	
	die " chunk_size must be greater than 0!\n"
           unless $self->chunk_size > 0;
	
	my $parent_def = $self->parent_def;
	my $parent_seq = $self->parent_seq;

	if(ref($parent_seq) eq 'SCALAR'){
	    $self->parent_seq_length(length($$parent_seq));
	}
	else{
	    $self->parent_seq_length($parent_seq->length);
	}

	#decide on number of chunks
	my $t_c = int($self->parent_seq_length/$self->chunk_size);
	$t_c++ if(!$t_c || ($self->parent_seq_length-$self->chunk_size*$t_c) > $self->min_size);

	$self->total_chunks($t_c);
	my $l = $self->chunk_size();
	
	$self->{INDEX} = 0;
        for (my $i=0; $i < $t_c; $i++){
	    my $offset = $i * $l;

	    my $is_last = ($t_c - 1 == $i) ? 1 : 0;
	    my $is_first = ($i == 0) ? 1 : 0;

	    #set start and end coordinates in 1 base space
	    my $start = $offset + 1;
	    my $end = ($is_last) ? $self->parent_seq_length : $offset + $l;

	    my $chunk = new FastaChunk();
	    $chunk->parent_def($parent_def);
	    $chunk->seqid(Fasta::def2SeqID($parent_def));
	    $chunk->size($l); #the max size of a chunk
	    $chunk->offset($offset);
	    $chunk->number($i);
	    $chunk->is_last($is_last);
	    $chunk->is_first($is_first);
	    $chunk->start($start);
	    $chunk->end($end);
	    $chunk->parent_seq_length($self->parent_seq_length());

	    if($self->flank()){
		my $flank = $self->flank();
		$chunk->flank($flank);

		#upstream
		my $B1 = ($offset+1) - $flank;
		$B1 = 1 if($B1 < 1);
		my $E1 = $offset;
		my $L1 = ($E1-$B1)+1;
		if($E1 >= $B1 && $L1 > 0){
		    $chunk->{_upstream}{start} = $B1;
		    $chunk->{_upstream}{end} = $E1;
		}

		#downstream
		my $B2 = $offset + $chunk->length + 1;
		my $E2 = $offset + $chunk->length + $flank;
		$E2 = $chunk->parent_seq_length if($E2 > $chunk->parent_seq_length);
		my $L2 = ($E2-$B2)+1;
		if($E2 >= $B2 && $L2 > 0){
                    $chunk->{_downstream}{start} = $B2;
                    $chunk->{_downstream}{end} = $E2;
		}
	    }

	    if(ref($parent_seq) eq 'SCALAR'){
		$chunk->seq(substr($$parent_seq, $chunk->offset_w_flank, $chunk->length_w_flank));
	    }
	    else{
		$chunk->seq($parent_seq);
	    }	   
	       
	    push(@{$self->{chunks}}, $chunk);
        }

	$self->{parent_fasta} = undef; #depricated
	$self->{parent_seq_ref} = undef; #depricated
}
#-------------------------------------------------------------------------------
sub get_chunk {
	my $self = shift;
	my $i    = shift;

	return $self->{chunks}->[$i] || undef;
}
#-------------------------------------------------------------------------------
sub next_chunk {
	my $self = shift;
	my $i    = $self->{INDEX}++;

	return $self->{chunks}->[$i] || undef;
}
#-------------------------------------------------------------------------------
sub last_chunk {
	my $self = shift;

	return $self->{chunks}->[-1];
}
#-------------------------------------------------------------------------------
sub all_chunks {
   my $self = shift;

   return $self->{chunks};
}
#-------------------------------------------------------------------------------
sub add_chunk {
	my $self  = shift;
	my $chunk = shift;

	push(@{$self->{chunks}}, $chunk);
}
#-------------------------------------------------------------------------------
sub replace {
   my $self = shift;
   my $chunk = shift;
   
   $self->{chunks}->[$chunk->number] = $chunk;
}
#-------------------------------------------------------------------------------
sub reset {
   my $self = shift;
   $self->{INDEX} = 0;
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


