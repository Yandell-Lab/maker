#----------------------------------------------------------------------------
#----                            FastaChunk                              ---- 
#----------------------------------------------------------------------------
package FastaChunk;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Fasta;
use FastaFile;
use FastaSeq;

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
sub seq {
    my $self = shift;
    my $arg  = shift;

    if(defined($arg)){
        $arg = $$arg while(ref($arg) eq 'REF');
        my $seq_ref = (ref($arg) eq '') ? \$arg : $arg;

	$self->{seq} = $seq_ref;
    }
    else{
	if(! ref($self->{seq}) || ref($self->{seq}) eq 'SCALAR'){
	    return $self->{seq};
	}
	else{
	    my $seq = $self->{seq}->subseq($self->start, $self->end);
	    return \$seq;
	}
    }
}
#-------------------------------------------------------------------------------
sub offset_w_flank {
    my $self = shift;
    return ($self->{_upstream}) ? $self->{_upstream}{start} - 1 : $self->{offset};
}
#-------------------------------------------------------------------------------
sub start_w_flank {
    my $self = shift;
    return ($self->{_upstream}) ? $self->{_upstream}{start} : $self->{start};
}
#-------------------------------------------------------------------------------
sub end_w_flank {
    my $self = shift;
    return ($self->{_downstream}) ? $self->{_downstream}{end} - 1 : $self->{end};
}
#-------------------------------------------------------------------------------
sub length_w_flank {
    my $self = shift;

    return abs($self->end_w_flank - $self->start_w_flank) + 1;
}
#-------------------------------------------------------------------------------
sub upstream_seq {
    my $self = shift;
    my $arg  = shift;

    return if(! $self->{flank} || ! $self->{_upstream});

    if($self->{_upstream}{seq}){
	return $self->{_upstream}{seq};
    }
    elsif($self->{seq} && ref($self->{seq}) && ref($self->{seq}) ne 'SCALAR'){
	my $B = $self->{_upstream}{start};
	my $E = $self->{_upstream}{end};
	my $L = abs($E-$B) +1;
	my $seq =  substr_o($self->{seq}, $B-1, $L);
	return \$seq;
    }
}
#-------------------------------------------------------------------------------
sub downstream_seq {
    my $self = shift;
    my $arg  = shift;

    return if(! $self->{flank} || ! $self->{_downstream});

    if($self->{_downstream}{seq}){
	return $self->{_downstream}{seq};
    }
    elsif($self->{seq} && ref($self->{seq}) && ref($self->{seq}) ne 'SCALAR'){
	my $B = $self->{_downstream}{start};
	my $E = $self->{_downstream}{end};
	my $L = abs($E-$B) +1;
	my $seq =  substr_o($self->{seq}, $B-1, $L);
	return \$seq;
    }
}
#-------------------------------------------------------------------------------
sub write_file {
	my $self      = shift;
	my $file_name = shift;

	$self->fasta_file_location($file_name);

	FastaFile::writeFile($self->fasta_ref, $file_name);
}
#-------------------------------------------------------------------------------
sub erase_fasta_file {
	my $self     = shift;
	my $location = shift;

	if    (defined($location)){
		unlink($location);
	}
	elsif (defined($self->fasta_file_location)){
		unlink($self->fasta_file_location);
	}
	else {
		print STDERR "cant find a file to erase!\n";
	}
}
#-------------------------------------------------------------------------------
sub fasta {
	my $self = shift;

	my $def = $self->def();
	my $seq = $self->seq();

	return Fasta::toFasta($def, $seq);
}
#-------------------------------------------------------------------------------
sub fasta_ref {
	my $self = shift;

	my $def = $self->def();
	my $seq = $self->seq();

	return Fasta::toFastaRef($def, $seq);
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

        #print STDERR "FastaChunk::AutoLoader called for: ",
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
