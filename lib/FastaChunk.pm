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
sub def {
    my $self = shift;
    my $i = $self->number;
    my $l = $self->length;
    my $offset = $self->offset;
    my $def = $self->parent_def. " CHUNK number:$i size:$l offset:$offset";
    return $def;
}
sub def_w_flank {
    my $self = shift;
    my $i = $self->number;
    my $l = $self->length_w_flank;
    my $offset = $self->offset_w_flank;
    my $def = $self->parent_def. " CHUNK number:$i size:$l offset:$offset";
    return $def;
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
	    if($self->{_upstream} || $self->{_downstream}){ #trim flank
		my $o = $self->start - $self->start_w_flank;
		my $l = $self->length;
		my $s = substr(${$self->{seq}}, $o, $l);
		return \$s;
	    }
	    else{ #no flank just give everything
		return $self->{seq};
	    }
	}
	else{ #give region w/o flank
	    my $seq = $self->{seq}->subseq($self->start, $self->end);
	    return \$seq;
	}
    }
}

sub seq_w_flank {
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
            my $seq = $self->{seq}->subseq($self->start_w_flank, $self->end_w_flank);
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
    return ($self->{_downstream}) ? $self->{_downstream}{end} : $self->{end};
}
#-------------------------------------------------------------------------------
sub length {
    my $self = shift;
    return abs($self->end - $self->start) + 1;
}
sub length_w_flank {
    my $self = shift;
    return abs($self->end_w_flank - $self->start_w_flank) + 1;
}
#-------------------------------------------------------------------------------
sub upstream_seq {
    my $self = shift;
    my $arg  = shift;

    return if(! $self->{flank} || ! $self->{_upstream});

    if($self->{seq} && ref($self->{seq}) && ref($self->{seq}) ne 'SCALAR'){
	my $B = $self->{_upstream}{start};
	my $E = $self->{_upstream}{end};
	my $L = abs($E-$B) +1;
	my $seq =  substr_o($self->{seq}, $B-1, $L);
	return \$seq;
    }
    elsif($self->{seq}){
	my $L = ($self->{_upstream}{end} - $self->{_upstream}{start})+1;
	my $seq =  substr($self->{seq}, 0, $L);
        return \$seq;
    }
}
#-------------------------------------------------------------------------------
sub downstream_seq {
    my $self = shift;
    my $arg  = shift;

    return if(! $self->{flank} || ! $self->{_downstream});

    if($self->{seq} && ref($self->{seq}) && ref($self->{seq}) ne 'SCALAR'){
	my $B = $self->{_downstream}{start};
	my $E = $self->{_downstream}{end};
	my $L = abs($E-$B) +1;
	my $seq =  substr_o($self->{seq}, $B-1, $L);
	return \$seq;
    }
    elsif($self->{seq}){
	my $B = $self->{_downstream}{start} - $self->offset_w_flank;
	my $L = ($self->{_downstream}{end} - $self->{_downstream}{start})+1;
	my $seq =  substr($self->{seq}, $B-1, $L);
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
sub write_file_w_flank {
	my $self      = shift;
	my $file_name = shift;

	$self->fasta_file_location($file_name);

	FastaFile::writeFile($self->fasta_ref_w_flank, $file_name);
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
sub fasta_w_flank {
	my $self = shift;

	my $def = $self->def_w_flank();
	my $seq = $self->seq_w_flank();

	return Fasta::toFasta($def, $seq);
}
#-------------------------------------------------------------------------------
sub fasta_ref {
	my $self = shift;

	my $def = $self->def();
	my $seq = $self->seq();

	return Fasta::toFastaRef($def, $seq);
}
sub fasta_ref_w_flank {
	my $self = shift;

	my $def = $self->def_w_flank();
	my $seq = $self->seq_w_flank();

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
