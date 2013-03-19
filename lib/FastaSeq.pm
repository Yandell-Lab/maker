#-------------------------------------------------------------------------------
#------                            FastaSeq                            ---------
#-------------------------------------------------------------------------------
package FastaSeq;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
use Bio::DB::Fasta;
use GI;
use Error qw(:try);
use Error::Simple;
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
#override native method because it can break on NFS
sub seq {
   my $self = shift;

   my $seq;
   my $fail = 0;
   my $ok = 0;
   while(!$ok && $fail < 5){
       try {
	   $seq = $self->SUPER::seq();
	   $ok++;
       }
       catch Error::Simple with {
	   my $E = shift;

	   $fail++;
	   if($fail == 5){
	       print STDERR $E->stacktrace;
	       throw $E;
	   }

	   sleep 10; #NFS?

	   if($fail == 3){
	       my $id = $self->id;
	       my $locs = $self->{LOCS};
	       my $index = GI::build_fasta_index($locs);
	       $index->drop_from_global_index();
	       $index = GI::build_fasta_index($locs);
	       my $obj = $index->get_Seq_by_id($id);

	       #retry because of weird NFS
	       if(! $obj->{db}){
		   sleep 10;
		   $index->drop_from_global_index();
		   $index = GI::build_fasta_index($locs);
		   $obj = $index->get_Seq_by_id($id);
	       }

	       #try one more time
	       if(! $obj->{db}){
		   sleep 10;
		   $index->drop_from_global_index();
		   $index = GI::build_fasta_index($locs);
		   $obj = $index->get_Seq_by_id($id);
	       }
	       
	       next if(! $obj->{db});

	       while(defined(my $key = each %$obj)){
		   $self->{$key} = $obj->{$key};
 	       }
	       $self->{LOCS} = $locs if(! $self->{LOCS});
	   }
       };
   }

   return $seq;
}
#-------------------------------------------------------------------------------
sub convert {
    my $class = shift;
    my $locs = shift;
    my $obj = shift;

    return if(!$obj);

    die "ERROR: Object is not a Bio::PrimarySeq::Fasta\n".
	"It is ".ref($obj)."\n"	if(ref($obj) ne 'Bio::PrimarySeq::Fasta');

    bless($obj, $class);
    
    if(!$locs && $obj->{db}){
	$locs = [keys %{$obj->{db}->{cacheseq}}];
	$locs = [$obj->{db}->{dirname}."/".$obj->{db}->{offsets}->{__file_0}] if(!@$locs);
    }
    $obj->{LOCS} = $locs || [];

    confess "ERROR: Attempt to create empty FastaSeq\n" if(!@{$obj->{LOCS}});

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

    if(! $seq->{db}){
	die join("\n", @$locs)."\n";
    }

    die "ERROR: Could not reestablish DB to thaw FastaSeq for Storable\n"
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

    return '' if($len == 0);

    #always work with references
    $seq = $$seq while(ref($seq) eq 'REF');
    my $seq_ref = (ref($seq) eq '') ? \$seq : $seq;

    if(ref($seq_ref) eq 'SCALAR'){
	$len = CORE::length($$seq_ref) - $off if(! defined($len));
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
	return CORE::length($$seq_ref);
    }
    else{
	return $seq_ref->length;
    }
}

1;
