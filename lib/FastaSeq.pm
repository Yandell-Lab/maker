#-------------------------------------------------------------------------------
#------                            FastaSeq                               ---------
#-------------------------------------------------------------------------------
package FastaSeq;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
use Bio::PrimarySeqI;
use GI;

@ISA = qw(Bio::PrimarySeq::Fasta Exporter);
@EXPORT = qw(substr_o length_o);

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
#this object is a Bio::PrimarySeq::Fasta object that provides
#serialization/deserialization support
sub new {
	my ($class, @args) = @_;

	$class = 'FastaSeq' if(ref($class));

	my $self = $class->SUPER::new(@args);
	bless($self, $class);

	return $self;
}
#-------------------------------------------------------------------------------
sub convert {
    my $class = shift;
    my $obj = shift;

    return if(!$obj);

    bless($obj, $class);
    
    return $obj;
}
#-------------------------------------------------------------------------------
sub STORABLE_freeze {
    my ($self, $cloning) = @_;

    my $locs = [keys %{$self->{db}->{cacheseq}}];
    if(!@$locs){
	$locs = [$self->{db}->{dirname}."/".$self->{db}->{offsets}->{__file_0}];
    }
    my $id = $self->id;

    return '', $locs, \$id;
}
#-------------------------------------------------------------------------------
sub STORABLE_thaw {
    my ($self, $cloning, $serialized, $locs, $id) = @_;
    
    $id = $$id while(ref($id) eq 'REF' || ref($id) eq 'SCALAR');

    my $index = GI::build_fasta_index($locs);

    my $seq = $index->get_Seq_by_id($id);

    while(my $key = each %$seq){
	$self->{$key} = $seq->{$key};
    }
}
#-------------------------------------------------------------------------------
#convinience function to allow substring type context on seq object 
sub substr_o {
    my $seq = shift;
    my $off = shift;
    my $len = shift;

    #always work with references
    $seq = $$seq while(ref($seq) eq 'REF');
    my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

    if(ref($seq_ref) eq 'SCALAR'){
	$len = length($$seq_ref) - $off if(! defined($len));
	return substr($$seq_ref, $off, $len);
    }
    else{
	$len = $seq_ref->length - $off if(! defined($len));
	return $seq_ref->subseq($off+1, $off+$len);
    }
}
#-------------------------------------------------------------------------------
#convinience function to allow substring type context on seq object 
sub length_o {
    my $seq = shift;

    #always work with references
    $seq = $$seq while(ref($seq) eq 'REF');
    my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

    if(ref($seq_ref) eq 'SCALAR'){
	return length($$seq_ref);
    }
    else{
	return $seq_ref->length;
    }
}

1;
