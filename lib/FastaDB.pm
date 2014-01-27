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
use GI;

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
	push(@{$self->{index}}, _safe_new($file, @args));
	carp "Calling out to BioPerl get_PrimarySeq_stream" if($main::debug);
	push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);

	#build reverse index to get the correct index based on file name
	$self->{name2index}{$title} = $self->{index}->[-1];       
	$self->{path2index}{$file} = $self->{index}->[-1];       
    }
    
    return $self;
}
#-------------------------------------------------------------------------------
#handles version issues that arise with AnyDBM
{
my %localized;
sub _safe_new {
    my ($file, @args) = @_;

    #which DBM am I using?
    my %args = @args;
    if($G_DB{$file} && $args{'-reindex'}){
	delete($G_DB{$file});
    }
    elsif($G_DB{$file}){
	$G_DB{$file}[1]++; #add count of active references
	return $G_DB{$file}[0];
    }

    my @ext = ($AnyDBM_File::ISA[0] eq 'DB_File') ? 
	('index') : ('index.dir', 'index.pag');
    my $db;
    try {
	my $exists = grep {-f "$file.$_"} @ext;
	$exists = ($exists == @ext) ? 1 : 0;
	my $tmp = GI::get_global_temp();
	my $localize;
	if(defined($args{'-localize'})){
	    $localize = $args{'-localize'};
	    delete($args{'-localize'});
	}
	elsif(!$main::nolocal && GI::is_NFS_mount($file) && !GI::is_NFS_mount($tmp)){
	    $localize = 1;
	}

	#only lock if this is being newly created
	my $lock;
	if(!$exists || $args{'-reindex'}){
	    while(!$lock || !$lock->maintain(30)){
		carp "Calling File::NFSLock::new" if($main::debug);
		$lock = new File::NFSLock("$file.index", 'EX', 1800, 60);
		$exists = grep {-f "$file.$_"} @ext; #check again
		$exists = ($exists == @ext) ? 1 : 0;
		if($exists && !$args{'-reindex'}){
		    undef $lock;
		    last;
		}
	    }
	}

	(my $dir = $file) =~ s/([^\/]+)$//;
	my ($name) = $1;

	#set where active index should end up
	my $fdir = $dir;
	if($localize){ #copy things locally if NFS mount and TMP is available
	    if($localized{$file}){
		$fdir = $localized{$file};
	    }
	    else{
		$fdir = GI::get_global_temp();
		$fdir .= '/'.Digest::MD5::md5_hex($file);
		mkdir($fdir) unless(-f $fdir);
		$localized{$file} = $fdir;
	    }
	}
	my $ffile = "$fdir/$name";
	my $rank  = GI::RANK();

	#copy other files locally if they already exist globally
	symlink($file, $ffile) if(!-f $ffile);
	if($exists && !$args{'-reindex'}){
	    foreach my $ext (@ext){
		last if($dir eq $fdir); #must be different location
		my $gi = "$dir/$name.$ext";
		my $fi = "$fdir/$name.$ext";
		my $ti = "$fdir/$name.$ext.$rank";
		
		#always copy if newer than current
		my $gitime = (stat($gi))[9] || 0;
		my $fitime = (stat($fi))[9] || 0;
		if(! -f $fi || $fitime < $gitime){
		    File::Copy::copy($gi, $ti) or confess "ERROR: Copy failed: $!";
		    sleep 10 if(! -f $ti); #NFS
		    File::Copy::move($ti, $fi) or confess "ERROR: Move failed: $!";
		}
	    }
	}
	else{ #does not exist or must be reindexed
	    my $tdir  = "$fdir/.dbtmp$rank";
	    my $tfile = "$tdir/$name";
	    File::Path::rmtree($tdir) if(-d $tdir);
	    mkdir($tdir);
	    symlink($file, $tfile);
	    foreach my $ext (@ext){
		my $ti = "$tdir/$name.$ext";
		unlink($ti) if(-f $ti);
	    }

	    #build the index (block to localize db index failures)
	    carp "Calling out to BioPerl Bio::DB::Fasta::new" if($main::debug);
	    {
		local $SIG{'__WARN__'} = sub { die $_[0]; };
		$db = new Bio::DB::Fasta($tfile, @args);
	    }
	    untie(%{$db->{offsets}}); #untie first to flush (for new index)

	    foreach my $ext (@ext){
                my $ti = "$tdir/$name.$ext";
		my $fi = "$fdir/$name.$ext";
		sleep 10 if(! -f $ti); #NFS
		File::Copy::move($ti, $fi) or confess "ERROR: Move failed: $!";
            }
	    File::Path::rmtree($tdir);

	    #copy new files back globally
	    if($localize){
		foreach my $ext (@ext){
		    my $fi = "$fdir/$name.$ext";
		    my $ti = "$dir/$name.$ext.$rank";
		    my $gi = "$dir/$name.$ext";
		    File::Copy::copy($fi, $ti) or confess "ERROR: Copy failed: $!";
		    sleep 10 if(! -f $ti); #NFS
		    File::Copy::move($ti, $gi) or confess "ERROR: Move failed: $!";

		    #make timestamps equal (avoids iterative reindexing)
		    my $fitime = (stat($fi))[9] || 0;
		    utime($fitime, $fitime, $gi);
		}
	    }
	}

	#now mount the finished index
	delete($args{'-reindex'});
	carp "Calling out to BioPerl Bio::DB::Fasta::new" if($main::debug);
	{
	    local $SIG{'__WARN__'} = sub { die $_[0]; };
	    $db = new Bio::DB::Fasta($ffile, %args);
	}
    
	$lock->unlock if($lock);
    }
    catch Error::Simple with {
        my $E = shift;

	#try direct indexing with no lock? Dangerous, but why not it's already failing
	foreach my $ext (@ext){
	    unlink("$file.$ext") if(-f "$file.$ext"); #delete failed indexes
	}
	carp "Calling out to BioPerl Bio::DB::Fasta::new" if($main::debug);
	$db = new Bio::DB::Fasta($file, @args);
    };


    $G_DB{$file}[0] = $db; #creates a global index
    $G_DB{$file}[1]++;
    return $db;
}
}
#-------------------------------------------------------------------------------
#first reindex attempt will just try and use non localized index
#second attempt tries to rebuild localized index
#other attempts alternate between localized and non-localized index
sub reindex {
    my $self = shift;
    my $flag = shift;

    my $locs = $self->{locs};
    $flag = 1 if(!defined($flag));
    my @args = ('-makeid' => \&makeid);
    push(@args, '-reindex' => 1) if($flag);
    push(@args, '-localize' => 0) if($flag % 2 == 0);
    $self->drop_from_global_index(); #make sure I get new ones

    my @files = grep {-f $_} @$locs;
    foreach my $dir (grep {-d $_} @$locs){
	push(@files, grep {-f $_ && /\.(fa|fasta|fast|FA|FASTA|FAST|dna|\.mpi\.\d+\.\d+)$/} <$dir/*>);
    }

    #clear old index array
    $self->{index} = [];
    $self->{stream} = [];
    $self->{name2index} = {};
    $self->{path2index} = {};

    #nothing to index
    return $self if(! @files);

    #rebuilt build indexes
    my $lock;
    if(($lock = new File::NFSLock("$files[0].reindex", 'NB', 0, 50)) && $lock->maintain(30)){ #reindex lock
	foreach my $file (@files){
	    push(@{$self->{index}}, _safe_new($file, @args));
	    carp "Calling out to BioPerl get_PrimarySeq_stream" if($main::debug);
	    push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);
	    
	    #build reverse index to get the correct index based on file name
	    my ($title) = $file =~ /([^\/]+)$/;
	    $self->{name2index}{$title} = $self->{index}->[-1];
	    $self->{path2index}{$file} = $self->{index}->[-1];
	}
	$lock->unlock;
    }
    else{
	#pause and wait for other process to reindex
	my $lock = $lock = new File::NFSLock("$files[0].reindex", 'EX', undef, 50);

	#just rebuild (not destroy old index))
	foreach my $file (@files){
	    @args = ('-makeid' => \&makeid);
	    push(@{$self->{index}}, _safe_new($file, @args));
	    carp "Calling out to BioPerl get_PrimarySeq_stream" if($main::debug);
	    push(@{$self->{stream}}, $self->{index}[-1]->get_PrimarySeq_stream);
	    
	    #build reverse index to get the correct index based on file name
	    my ($title) = $file =~ /([^\/]+)$/;
	    $self->{name2index}{$title} = $self->{index}->[-1];
	    $self->{path2index}{$file} = $self->{index}->[-1];
	}	
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

    my $r_ind = $self->{name2index}; #reverse index
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

    my $r_ind = $self->{name2index}; #reverse index
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
	my @keys = keys %{$self->{name2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{name2index}{@keys};
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
	my @keys = keys %{$self->{name2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{name2index}{@keys};
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
	my @keys = keys %{$self->{name2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{name2index}{@keys};
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
	my @keys = keys %{$self->{name2index}}; #all file names
	$source =~ s/.*\/([^\/]+)$/$1/;
	@keys = grep {/$source(\.mpi\.\d+\.\d+)?$/} @keys;
	@index =  $self->{name2index}{@keys};
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
	$self->{name2index}{$title} = $self->{index}->[-1];
	$self->{path2index}{$file} = $self->{index}->[-1];
    }
}

sub drop_from_global_index{
    my $self = shift;

    my $stat = 0;
    while(my $key = each %{$self->{path2index}}){
	$stat = 1 if($G_DB{$key});
        delete($G_DB{$key});
    }

    return $stat;
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
sub DESTROY {
    my $self = shift;

    #remove global index when all pointers go out of scope
    foreach my $file (keys %{$self->{path2index}}){
	$G_DB{$file}[1]--;
	delete($G_DB{$file}) if($G_DB{$file}[1] <= 0);
    }
}
#-------------------------------------------------------------------------------
1;
