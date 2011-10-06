#-------------------------------------------------------------------------------
#------                           FastaDB                              ---------
#-------------------------------------------------------------------------------
package FastaDB;
use strict;
use vars qw(@ISA @EXPORT $VERSION %G_DB);
use Exporter;
use PostData;
use FileHandle;
use URI::Escape;
use Bio::DB::Fasta;
use File::NFSLock;
use FastaSeq;

@ISA = qw(
          );
#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $class = shift;
    my $locs  = shift;

    my $self = {};

    bless ($self, $class);

    $self->{index} = [];
    $self->{stream} = [];
    $self->{locs} = $locs;
    my @args = @_;
    push (@args, ('-makeid' => \&makeid));

    my @files = grep {-f $_} @$locs;
    foreach my $dir (grep {-d $_} @$locs){
	push(@files, grep {-f $_ && /\.(fa|fasta|fast|FA|FASTA|FAST|dna|\.mpi\.\d+\.\d+)$/} <$dir/*>);
    }

    #build indexes
    while (my $file = shift @files){ 
	my ($title) = $file =~ /([^\/]+)$/;
	if($G_DB{$title}){
	    push(@{$self->{index}}, $G_DB{$title});
	    push(@{$self->{stream}}, $G_DB{$title}->get_PrimarySeq_stream);

	    #build reverse index to get the correct index based on file name
	    $self->{file2index}{$title} = $self->{index}->[-1];
	    next;
	}

	my $lock; #do non blocking lock to skip over active
	if (! -e "$file.index" || -e ".NFSLock.$file.index.NFSLock"){ #index is not ready to use
	    if(($lock = new File::NFSLock("$file.index", 'NB', undef, 50))){
		if(-e "$file.index"){ #release lock because index is ready to use
		    $lock->unlock();
		    $lock = undef;
		}
		elsif(! $lock->maintain(30)){ #can't get maintainer to lock
		    push(@files, $file);
		    $lock->unlock();
		    sleep 1;
		    next;
		}
	    }
	    else{
		push(@files, $file);
		sleep 1;
		next;
	    }
	}

	push(@{$self->{index}}, new Bio::DB::Fasta($file, @args));
	push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);

	#build reverse index to get the correct index based on file name
	$self->{file2index}{$title} = $self->{index}->[-1];
	
	$lock->unlock if($lock);
    }
    
    return $self;
}
#-------------------------------------------------------------------------------
sub reindex {
    my $self = shift;

    my $locs = $self->{locs};
    my @args = ('-reindex' => 1,
		'-makeid' => \&makeid
		);

    my @files = grep {-f $_} @$locs;
    foreach my $dir (grep {-d $_} @$locs){
	push(@files, grep {-f $_ && /\.(fa|fasta|fast|FA|FASTA|FAST|dna|\.mpi\.\d+\.\d+)$/} <$dir/*>);
    }

    #clear old index array
    $self->{index} = [];
    $self->{stream} = [];
    $self->{file2index} = {};

    #nothing to index
    return $self if(! @files);

    #rebuilt build indexes
    my $lock;
    if(($lock = new File::NFSLock("$files[0].reindex", 'NB', 0, 50)) && $lock->maintain(30)){ #reindex lock
	foreach my $file (@files){
	    my $lock;
	    if(($lock = new File::NFSLock("$file.index", 'EX', undef, 50)) && $lock->maintain(30)){ #stnd index lock
		push(@{$self->{index}}, new Bio::DB::Fasta($file, @args));
		push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);		

		#build reverse index to get the correct index based on file name
		my ($title) = $file =~ /([^\/]+)$/;
		$self->{file2index}{$title} = $self->{index}->[-1];
		$lock->unlock;
	    }
	    else{
		die "ERROR: Could not get lock for re-indexing\n\n";
	    }
	}
	$lock->unlock;
    }
    else{
	#pause and wait for other process to reindex
	my $lock = $lock = new File::NFSLock("$files[0].reindex", 'EX', undef, 40);
	$lock->unlock;
    }

    return $self;
}
#-------------------------------------------------------------------------------
sub _close_index {
    my $self = shift;

    my @index = @{$self->{index}};

    foreach my $db (@index){
	$db->_close_index;
    }

}
#-------------------------------------------------------------------------------
#uses hit info to search all indexes faster
sub get_Seq_for_hit {
    my $self = shift;
    my $hit = shift;

    my $r_ind = $self->{file2index}; #reverse index
    my $id = (ref($hit) eq '') ? ($hit =~ /^>([^\s]+)/)[0] : $hit->name;
    my $source = (ref($hit) eq '') ? undef : $hit->{_file};
    my $description = (ref($hit) eq '') ? $hit : $hit->description;

    if($description =~ /MD5_alias=(\S+)/){
	$id = $1;
    }

    my $dbf = (ref($hit) eq '') ? undef : $hit->database_name;

    return $self->get_Seq_by_id($id, $source) if(! defined $dbf);

    ($dbf) = $dbf =~ /([^\/]+)$/;

    my $fastaObj;
    if(exists $r_ind->{$dbf}){
	my $db = $r_ind->{$dbf};
	$fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
    }
    
    if(! $fastaObj){
	my @files = grep {!/^$dbf$/} keys %$r_ind; #check remaining files
	foreach my $dbf (@files){
	    my $db = $r_ind->{$dbf};
	    $fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
	    last if($fastaObj);
	}
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
#uses hit info to search all indexes faster
sub header_for_hit {
    my $self = shift;
    my $hit = shift;

    my $r_ind = $self->{file2index}; #reverse index
    my $id = (ref($hit) eq '') ? ($hit =~ /^>([^\s]+)/)[0] : $hit->name;
    my $source = (ref($hit) eq '') ? undef : $hit->{_file};
    my $description = (ref($hit) eq '') ? $hit : $hit->description;

    if($description =~ /MD5_alias=(\S+)/){
	$id = $1;
    }

    my $dbf = (ref($hit) eq '') ? undef : $hit->database_name;

    return $self->header($id, $source) if(! defined $dbf);

    ($dbf) = $dbf =~ /([^\/]+)$/;

    my $fastaObj;
    my $header;
    if(exists $r_ind->{$dbf}){
	my $db = $r_ind->{$dbf};
	$fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
	$header = $db->header($id) if($fastaObj);
    }
    
    if(! $fastaObj){
	my @files = grep {!/^$dbf$/} keys %$r_ind; #check remaining files
	foreach my $dbf (@files){
	    my $db = $r_ind->{$dbf};
	    $fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));

	    if($fastaObj){
		$header = $db->header($id);
		last
	    }
	}
    }

    return $header;
}
#-------------------------------------------------------------------------------
sub get_Seq_by_id {
    my $self = shift;
    my $id = shift;
    my $source = shift;

    my @index = @{$self->{index}};
    if($source){
	my @keys = keys %{$self->{file2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{file2index}{@keys};
    }

    my $fastaObj;
    foreach my $db (@index){
	$fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
	last if($fastaObj);
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
sub get_Seq_by_alias {
    my $self = shift;
    my $alias = shift;
    my $source = shift;

    $alias =~ /MD5_alias=(\S+)/;
    $alias = $1;

    return undef if(! defined $alias);

    my @index = @{$self->{index}};
    if($source){
	my @keys = keys %{$self->{file2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{file2index}{@keys};
    }

    my $fastaObj;
    foreach my $db (@index){
	$fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($alias));
	last if($fastaObj);
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
sub next_seq {
    my $self = shift;

    my $fastaObj;
    while(my $stream = $self->{stream}[0]){
	if($fastaObj = $stream->next_seq){
	    $fastaObj = FastaSeq->convert($self->{locs}, $fastaObj);
	    last if($fastaObj);
	}
	else{
	    shift @{$self->{stream}}
	}
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
sub get_all_ids {
    my $self = shift;

    my @index = @{$self->{index}};
    my @all_ids;
    foreach my $db (@index){
        my @ids = $db->get_all_ids;
	push(@all_ids, @ids);
    }

    return @all_ids;
}
#-------------------------------------------------------------------------------
sub header {
    my $self = shift;
    my $id = shift;
    my $source = shift;

    my @index = @{$self->{index}};
    if($source){
	my @keys = keys %{$self->{file2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{file2index}{@keys};
    }

    my $h;
    foreach my $db (@index){
	#do it this way first to avoid warnings
	my $fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
	next if (!$fastaObj);
	$h = $db->header($id);
	last;
    }

    return $h;
}
#-------------------------------------------------------------------------------
sub header_by_alias {
    my $self = shift;
    my $alias = shift;
    my $source = shift;

    $alias =~ /MD5_alias=(\S+)/;
    $alias = $1;

    return undef if(! defined $alias);

    my @index = @{$self->{index}};
    if($source){
	my @keys = keys %{$self->{file2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{file2index}{@keys};
    }

    my $h;
    foreach my $db (@index){
	#do it this way first to avoid warnings
	my $fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($alias));
	next if (!$fastaObj);
	$h = $db->header($alias);
	last;
    }

    return $h;
}
#-------------------------------------------------------------------------------
#for safe freezing of indexes using storable
sub STORABLE_freeze {
    my ($self, $cloning) = @_;

    my $locs = $self->{locs};

    return '', $locs;
}
#-------------------------------------------------------------------------------
#for safe thawing of indexes using storable
sub STORABLE_thaw {
    my ($self, $cloning, $serialized, $locs) = @_;

    $self->{locs} = $locs;

    my @args;
    push (@args, ('-makeid' => \&makeid));

    my @files = @$locs;
    foreach my $dir (grep {-d $_} @$locs){
	push(@files, grep {-f $_ && /\.(fa|fasta|fast|FA|FASTA|FAST|dna|\.mpi\.\d+\.\d+)$/} <$dir/*>);
    }
    
    #build indexes
    while (my $file = shift @files){ 
	push(@{$self->{index}}, new Bio::DB::Fasta($file, @args));
	push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);	

	#build reverse index to get the correct index based on file name
	my ($title) = $file =~ /([^\/]+)$/;
	$self->{file2index}{$title} = $self->{index}->[-1];
    }
}
#-------------------------------------------------------------------------------
sub add_to_global_index{
    my $self = shift;

    while(my $key = each %{$self->{file2index}}){
	my ($name) = $key =~ /([^\/]+)$/;
        $G_DB{$name} = $self->{file2index}->{$key};
    }
}
#-------------------------------------------------------------------------------
sub drop_from_global_index{
    my $self = shift;

    while(my $key = each %{$self->{file2index}}){
	my ($name) = $key =~ /([^\/]+)$/;
        delete($G_DB{$name});
    }
}
#-------------------------------------------------------------------------------
#------------------------------- SUBS ------------------------------------------
#-------------------------------------------------------------------------------
sub makeid {
    my $def = shift;

    my @ids;
    
    #get the standard id
    if($def =~ />(\S+)/){
	push(@ids, $1);
    }

    #get the MD5 ID if made by GI::split_db
    #otherwise just trim the standard name to get an alias
    #this is because the BLAST parser trims names
    if($def =~ /MD5_alias=(\S+)/){
	push(@ids, $1);
    }
    elsif(defined $ids[0] && length($ids[0]) > 78){
	push(@ids, substr($ids[0], 0, 78));
    }
    
    return @ids;
}
#-------------------------------------------------------------------------------
1;
