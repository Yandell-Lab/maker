#----------------------------------------------------------------------------
#----                             Iterator::Fasta                        ---- 
#----------------------------------------------------------------------------
package Iterator::Fasta;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Iterator;
use Fasta;
use Scalar::Util qw(openhandle);

@ISA = qw(Iterator);

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class = shift;
	my $arg   = shift;

        my $self = {};
        bless $self;

	$self->fileName($arg);
	$self->fileHandle($arg);
	$self->{SKIP} = {};

	return $self;
}
#-------------------------------------------------------------------------------
sub _set_number_of_entries {
	my $self = shift;
	my $arg  = $self->fileName();

	my $fh = new FileHandle();
	   $fh->open($arg);

	my $i = 0;
	my $line;
	while ($line = <$fh>){
		$i++ if $line =~ /^>/;
	}
	$fh->close();

	$self->{number_of_entries} = $i;
}
#-------------------------------------------------------------------------------
sub number_of_entries{
    my $self = shift;
    
    if(! defined ($self->{number_of_entries})){
	$self->_set_number_of_entries();
    }
    
    return $self->{number_of_entries};
}
#-------------------------------------------------------------------------------
sub skip_file {
    my $self = shift;
    my $file = shift;

    die "ERROR: Log file does not exist in Iterator::Fasta::skip_file\n" if(! -f $file);

    my $out_base = $file;
    $out_base =~ s/[^\/]+$//;
    open(IN, "< $file") || die "ERROR: Could not open the log file in Iterator::Fasta::skip_file\n";
    my %skip;
    while(my $line = <IN>){
	chomp $line;
	my @F = split("\t", $line);
	next unless(@F == 3);
	next unless($F[2] eq 'FINISHED');
	next unless(-d "$out_base/$F[1]");
	$skip{$F[0]}++;
    }
    close(IN);

    $self->{skip_file} = $file;
    $self->{SKIP} = \%skip;
}
#-------------------------------------------------------------------------------
sub find {
    die "ERROR: Iterator::Fasta::find is disabled\n";
}
#-------------------------------------------------------------------------------
sub nextEntry {
    my $self = shift;
    my $ref = $self->nextEntryRef;
    return ($ref) ? $$ref : undef;
}
#-------------------------------------------------------------------------------

{
my @SEEN;
my $COUNT = -1;
my @BUF; #buffer for pushing back meta character contamination
sub nextEntryRef {
    my $self = shift;

    my $fh = $self->fileHandle();
    
    if (! @BUF && ! openhandle($fh)){ #checks to see if file handle is open
	return undef; 
    }
    
    local $/ = "\n>";

    my $line;
    while($line = shift @BUF || <$fh>){
	$COUNT++;

	$line =~ s/>//;
	$line =~ s/>$//;
	$line = ">".$line;

	if($line =~ /^M\n?|\cM\n?/){
	    $line =~ s/^M\n?|\cM\n?/\n/g;
	    my @set = grep {$_ ne "\n" } split(/\n>/, $line);
	    foreach my $s (@set){
		$s = ">".$s if($s !~ /^>/);
	    }
	    $line = shift @set;
	    push(@BUF, @set);
	}

	#already seen so skip
	next if($SEEN[$COUNT]);

	#step forward in jumps if indicated
	if($self->{step} && $self->{step} != 1){
	    next if($COUNT < $self->{step} || $COUNT % $self->{step} != 0);
	}

	$SEEN[$COUNT]++;

	$line =~ /^>([^\s\n]+)/;
	next if(defined ($self->{SKIP}{$1}));

	local $/ = "\n";
	return \$line;
    }
    
    local $/ = "\n";

    #end of file but I was jumping using a step, so go back to start
    if($self->{step} && $self->{step} != 1){
	$self->{step} = undef; #remove step
	$COUNT = -1; #reset count
	$fh->setpos($self->startPos); #go to start of file
	$self->skip_file($self->{skip_file}) if($self->{skip_file} && -f $self->{skip_file}); #reload skip file
	return $self->nextEntryRef();
    }
    
    $fh->close();
    return undef;
}
}
#-------------------------------------------------------------------------------
sub nextFasta {#alias to nextEntry
   my $self = shift;
   my $ref = $self->nextEntryRef;

   return (ref($ref) eq '') ? $ref : $$ref;
}
#-------------------------------------------------------------------------------
sub nextFastaRef {#alias to nextEntry
   my $self = shift;
   return $self->nextEntryRef;
}
#-------------------------------------------------------------------------------
sub finished {
    my $self = shift;

    my $fh = $self->fileHandle();

    if (openhandle($fh)){ #checks to see if file handle is open                                                              
	return 0;
    }
    else{
	return 1;
    }
}

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Iterator::Fasta::AutoLoader called for: ",
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


