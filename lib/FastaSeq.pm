#-------------------------------------------------------------------------------
#------                            FastaSeq                            ---------
#-------------------------------------------------------------------------------
package FastaSeq;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
use Bio::DB::Fasta;
use GI;
use Carp;

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
	$self->{LOCS} = [];
	bless($self, $class);

	return $self;
}
#-------------------------------------------------------------------------------
#override native method because it is so slow
sub length {
   my $self = shift;

   return ($self->{stop} - $self->{start}) + 1;
}
#-------------------------------------------------------------------------------
sub convert {
    my $class = shift;
    my $obj = shift;
    my $locs = shift;

    return if(!$obj);

    bless($obj, $class);
    
    if(!$locs){
	$locs = [keys %{$obj->{db}->{cacheseq}}];
	$locs = [$obj->{db}->{dirname}."/".$obj->{db}->{offsets}->{__file_0}] if(!@$locs);
    }
    $obj->{LOCS} = $locs;

    return $obj;
}
#-------------------------------------------------------------------------------
sub STORABLE_freeze {
    my ($self, $cloning) = @_;

    my $locs = $self->{LOCS};

    my $id = $self->id;

    return '', $locs, \$id;
}
#-------------------------------------------------------------------------------
sub STORABLE_thaw {
    my ($self, $cloning, $serialized, $locs, $id) = @_;
    
    $id = $$id while(ref($id) eq 'REF' || ref($id) eq 'SCALAR');

    my $index = GI::build_fasta_index($locs);
    my $seq = $index->get_Seq_by_id($id);

    #retry because of weird NFS
    if(! $seq->{db}){
	sleep 10;
	$index = GI::build_fasta_index($locs);
	$seq = $index->get_Seq_by_id($id);
    }

    #try one more time
    if(! $seq->{db}){
	sleep 10;
	$index = GI::build_fasta_index($locs);
	$seq = $index->get_Seq_by_id($id);
    }

    confess "ERROR: Could not reestablish DB to thaw FastaSeq for Storable\n"
	if(! $seq->{db});

    while(defined(my $key = each %$seq)){
	$self->{$key} = $seq->{$key};
    }
    $self->{LOCS} = $locs if(! $self->{LOCS});
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
