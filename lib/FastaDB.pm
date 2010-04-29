#-------------------------------------------------------------------------------
#------                           FastaDB                              ---------
#-------------------------------------------------------------------------------
package FastaDB;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use URI::Escape;
use Bio::DB::Fasta;
use File::NFSLock;

@ISA = qw(
          );
#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $class = shift;
    my $dir   = shift;

    my $self = {};

    bless ($self, $class);

    die "ERROR: Directory $dir does not exist or is not a directory.\n" if(! -d $dir);

    $self->{dirname} = $dir;
    my @args = @_;
    push (@args, ('-makeid' => \&makeid));

    my @files = <$dir/*>;

    #identify fastas
    my @keep;
    foreach my $file (@files){
	next if (! -f $file);
	next unless ($file =~ /\.(fa|fasta|fast|FA|FASTA|FAST|dna)$/);
	push(@keep, $file);
    }

    #build indexes
    foreach my $file (@keep){
	if(my $lock = new File::NFSLock("$file.index", 'EX', undef, 40)){
	    if(! -e "$file.index"){ #maintain lock because I must build index
		$lock->maintain(30);
	    }
	    else{ #release lock because index exists and is ready to use
		$lock->unlock();
		$lock = undef;
	    }

	    push(@{$self->{index}}, new Bio::DB::Fasta($file, @args));

	    #build reverse index to get the correct index based on file name
	    my ($title) = $file =~ /([^\/]+)$/;
	    $self->{file2index}{$title} = $self->{index}->[-1];

	    $lock->unlock if($lock);
	}
	else{
	    die "ERROR: Could not get lock for indexing\n\n";
	}
    }

    return $self;
}
#-------------------------------------------------------------------------------
sub reindex {
    my $self = shift;
    my $dir = $self->{dirname};

    die "ERROR: Directory $dir does not exist or is not a directory.\n" if(! -d $dir);

    my @args = ('-reindex' => 1,
		'-makeid' => \&makeid
		);

    my @files = <$dir/*>;
    
    #identify fastas
    my @keep;
    foreach my $file (@files){
	next if (! -f $file);
	next unless ($file =~ /\.(fa|fasta|fast|FA|FASTA|FAST|dna)$/);
	push(@keep, $file);
    }

    #clear old index array
    $self->{index} = [];

    #rebuilt build indexes
    foreach my $file (@keep){
	if(my $lock = new File::NFSLock("$file.index", 'EX', undef, 40)){
	    $lock->maintain(30);
	    
	    push(@{$self->{index}}, new Bio::DB::Fasta($file, @args));
	    
	    #build reverse index to get the correct index based on file name
	    my ($title) = $file =~ /([^\/]+)$/;
	    $self->{file2index}{$title} = $self->{index}->[-1];
	}
	else{
	    die "ERROR: Could not get lock for re-indexing\n\n";
	}
    }

    return $self;
}
#-------------------------------------------------------------------------------
#uses hit info to search all indexes faster
sub get_Seq_for_hit {
    my $self = shift;
    my $hit = shift;

    my $r_ind = $self->{file2index}; #reverse index
    my $id = $hit->name;

    if($hit->description =~ /MD5_alias=(\S+)/){
	$id = $1;
    }

    my $dbf = $hit->database_name;

    return $self->get_Seq_by_id($id) if(! defined $dbf);    

    ($dbf) = $dbf =~ /([^\/]+)$/;

    my $fastaObj;
    if(exists $r_ind->{$dbf}){
	my $db = $r_ind->{$dbf};
	$fastaObj = $db->get_Seq_by_id($id);
    }
    
    if(! $fastaObj){
	my @files = grep {!/^$dbf$/} keys %$r_ind; #check remaining files
	foreach my $dbf (@files){
	    my $db = $r_ind->{$dbf};
	    $fastaObj = $db->get_Seq_by_id($id);
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
    my $id = $hit->name;

    if($hit->description =~ /MD5_alias=(\S+)/){
	$id = $1;
    }

    my $dbf = $hit->database_name;

    return $self->header($id) if(! defined $dbf);

    ($dbf) = $dbf =~ /([^\/]+)$/;

    my $fastaObj;
    my $header;
    if(exists $r_ind->{$dbf}){
	my $db = $r_ind->{$dbf};
	$fastaObj = $db->get_Seq_by_id($id);
	$header = $db->header($id) if($fastaObj);
    }
    
    if(! $fastaObj){
	my @files = grep {!/^$dbf$/} keys %$r_ind; #check remaining files
	foreach my $dbf (@files){
	    my $db = $r_ind->{$dbf};
	    $fastaObj = $db->get_Seq_by_id($id);

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

    my $fastaObj;
    foreach my $db (@{$self->{index}}){
	$fastaObj = $db->get_Seq_by_id($id);
	last if($fastaObj);
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
sub get_Seq_by_alias {
    my $self = shift;
    my $alias = shift;

    $alias =~ /MD5_alias=(\S+)/;
    $alias = $1;

    return undef if(! defined $alias);

    my $fastaObj;
    foreach my $db (@{$self->{index}}){
	$fastaObj = $db->get_Seq_by_id($alias);
	last if($fastaObj);
    }

    return $fastaObj;
}
#-------------------------------------------------------------------------------
sub header {
    my $self = shift;
    my $id = shift;

    my $h;
    foreach my $db (@{$self->{index}}){
	#do it this way first to avoid warnings
	my $fastaObj = $db->get_Seq_by_id($id);
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

    $alias =~ /MD5_alias=(\S+)/;
    $alias = $1;

    return undef if(! defined $alias);

    my $h;
    foreach my $db (@{$self->{index}}){
	#do it this way first to avoid warnings
	my $fastaObj = $db->get_Seq_by_id($alias);
	next if (!$fastaObj);
	$h = $db->header($alias);
	last;
    }

    return $h;
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
