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
use Error qw(:try);
use Carp;

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

	my $lock;
	while(! $lock || ! $lock->still_mine){
	    carp "Calling File::NFSLock::new" if($main::debug);
	    $lock = new File::NFSLock("$file.index", 'EX', 1800, 60);
	    carp "Calling File::NFSLock::maintain" if($main::debug);
	    $lock->maintain(30);
	}

	push(@{$self->{index}}, _safe_new($file, @args));
	carp "Calling out to BioPerl get_PrimarySeq_stream" if($main::debug);
	push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);

	#build reverse index to get the correct index based on file name
	$self->{file2index}{$title} = $self->{index}->[-1];
	
	$lock->unlock if($lock);
    }
    
    return $self;
}
#-------------------------------------------------------------------------------
#handles version issues that arise with AnyDBM
sub _safe_new {
    my ($file, @args) = @_;

    my $db;
    try {
	local $SIG{'__WARN__'} = sub {
	    die $_[0];
	};

	carp "Calling out to BioPerl Bio::DB::Fasta::new" if($main::debug);
	$db = new Bio::DB::Fasta($file, @args);
    }
    catch Error::Simple with {
        my $E = shift;

	unlink("$file.index") if(-f "$file.index"); #DB_File
	unlink("$file.index.dir") if(-f "$file.index.dir"); #GDBM_File
	unlink("$file.index.pag") if(-f "$file.index.pag"); #GDBM_File
	carp "Calling out to BioPerl Bio::DB::Fasta::new" if($main::debug);
	$db = new Bio::DB::Fasta($file, @args);
    };

    return $db;
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
		push(@{$self->{index}}, _safe_new($file, @args));
		carp "Calling out to BioPerl get_PrimarySeq_stream" if($main::debug);
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
	carp "Calling out to BioPerl Bio::DB::Fasta::_close_index" if($main::debug);
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
	carp "Calling out to FastaSeq::convert" if($main::debug);
	$fastaObj = FastaSeq->convert($self->{locs}, $db->get_Seq_by_id($id));
    }
    
    if(! $fastaObj){
	my @files = grep {!/^$dbf$/} keys %$r_ind; #check remaining files
	foreach my $dbf (@files){
	    my $db = $r_ind->{$dbf};
	    carp "Calling out to FastaSeq::convert" if($main::debug);
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

    carp "Calling out to FastaDB::STORABLE_freeze" if($main::debug);
    my $locs = $self->{locs};

    return '', $locs;
}
#-------------------------------------------------------------------------------
#for safe thawing of indexes using storable
sub STORABLE_thaw {
    my ($self, $cloning, $serialized, $locs) = @_;

    carp "Calling out to FastaDB::STORABLE_thaw" if($main::debug);
    $self->{locs} = $locs;

    my @args;
    push (@args, ('-makeid' => \&makeid));

    my @files = @$locs;
    foreach my $dir (grep {-d $_} @$locs){
	push(@files, grep {-f $_ && /\.(fa|fasta|fast|FA|FASTA|FAST|dna|\.mpi\.\d+\.\d+)$/} <$dir/*>);
    }
    
    #build indexes
    while (my $file = shift @files){ 
	push(@{$self->{index}}, _safe_new($file, @args));
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

    carp "Calling FastaDB::makeid from BioPerl Bio::DB::Fasta hook" if($main::debug);

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
