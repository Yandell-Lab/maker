#----------------------------------------------------------------------------
#----                             Iterator::Fasta                        ---- 
#----------------------------------------------------------------------------
package Iterator::Fasta;
use strict;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
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
	$self->{SEEN} = [];
	$self->{COUNT} = -1;
	$self->{BUF} = []; #buffer for pushing back meta character contamination
	
	return $self;
}
#-------------------------------------------------------------------------------
sub _set_number_of_entries {
	my $self = shift;
	my $arg  = $self->fileName();

	my $fh = new FileHandle();
	$fh->open($arg);
	$fh->setpos($self->startPos) if($self->startPos);
	$self->{number_of_entries} = scalar grep {/^>/} <$fh>;
	$fh->close();
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

    #die "ERROR: Log file does not exist in Iterator::Fasta::skip_file\n" if(! -f $file);

    my ($out_base, $name) = $file =~ /(.*\/)([^\/]+)$/;
    $name =~ s/master_datastore_index\.log$/datastore/;
    my %skip;
    if(-f $file){
	open(IN, "< $file") || die "ERROR: Could not open the log file in Iterator::Fasta::skip_file\n";
	while(my $line = <IN>){
	    chomp $line;
	    my @F = split("\t", $line);
	    next unless(@F == 3);
	    next unless($F[2] eq 'FINISHED' || $F[2] eq 'SKIPPED_SMALL');
	    next unless($F[1] =~ /^$name\/..\/..\/$F[0]\/?$/ || -d "$out_base/$F[1]");
	    $skip{$F[0]}++;
	}
	close(IN);
    }

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
sub nextEntryRef {
    my $self = shift;

    die "ERROR: Iterator type is: '".$self->{type}."' not 'full'\n"
	if($self->{type} && $self->{type} ne 'full');

    $self->{type} = 'full';

    my $fh = $self->fileHandle();
    
    $self->{BUF} = [] if(! $self->{BUF});
    if (! @{$self->{BUF}} && ! openhandle($fh)){ #checks to see if file handle is open
	return undef; 
    }

    my $line;
    {    
	local $/ = "\n>";
	while($line = shift @{$self->{BUF}} || <$fh>){
	    $self->{COUNT}++;
	    
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
		push(@{$self->{BUF}}, @set);
	    }
	    
	    #already seen so skip
	    next if($self->{SEEN}->[$self->{COUNT}]);
	    
	    #step forward in jumps if indicated
	    if($self->{step} && $self->{step} != 1){
		next if($self->{COUNT} < $self->{step} || $self->{COUNT} % $self->{step} != 0);
	    }
	    
	    $self->{SEEN}->[$self->{COUNT}]++;
	    
	    $line =~ /^>([^\s\n]+)/;
	    next if(defined ($self->{SKIP}{$1}));
	    
	    local $/ = "\n"; #just in case
	    return \$line;
	}
    }

    #end of file but I was jumping using a step, so go back to start
    if($self->{step} && $self->{step} != 1){
	$self->{step} = undef; #remove step
	$self->{COUNT} = -1; #reset count
	$fh->setpos($self->startPos); #go to start of file
	$self->skip_file($self->{skip_file}) if($self->{skip_file} && -f $self->{skip_file}); #reload skip file
	return $self->nextEntryRef();
    }
    
    $fh->close();
    return undef;
}
#-------------------------------------------------------------------------------
sub nextID {
    my $self = shift;

    die "ERROR: Iterator type is: '".$self->{type}."' not 'Def'\n"
	if($self->{type} && $self->{type} ne 'Def');

    $self->{type} = 'Def';

    my $fh = $self->fileHandle();
    
    if (! openhandle($fh)){ #checks to see if file handle is open
	return undef; 
    }
    
    my $line;
    while($line = <$fh>){
	next unless($line =~ /^>/);
	chomp $line;

	$self->{COUNT}++;

	#already seen so skip
	next if($self->{SEEN}->[$self->{COUNT}]);

	#step forward in jumps if indicated
	if($self->{step} && $self->{step} != 1){
	    next if($self->{COUNT} < $self->{step} || $self->{COUNT} % $self->{step} != 0);
	}

	$self->{SEEN}->[$self->{COUNT}]++;

	$line =~ /^>([^\s\n]+)/;
	next if(defined ($self->{SKIP}{$1}));

	return $1;
    }

    #end of file but I was jumping using a step, so go back to start
    if($self->{step} && $self->{step} != 1){
	$self->{step} = undef; #remove step
	$self->{COUNT} = -1; #reset count
	$fh->setpos($self->startPos); #go to start of file
	$self->skip_file($self->{skip_file}) if($self->{skip_file} && -f $self->{skip_file}); #reload skip file
	return $self->nextID();
    }
    
    $fh->close();
    return undef;
}
#-------------------------------------------------------------------------------
sub nextDef {
    my $self = shift;

    die "ERROR: Iterator type is: '".$self->{type}."' not 'Def'\n"
	if($self->{type} && $self->{type} ne 'Def');

    $self->{type} = 'Def';

    my $fh = $self->fileHandle();
    
    if (! openhandle($fh)){ #checks to see if file handle is open
	return undef; 
    }
    
    my $line;
    while($line = <$fh>){
	next unless($line =~ /^>/);
	chomp $line;
	$line =~ s/[\n\t\s\cM]+$//;

	$self->{COUNT}++;

	#already seen so skip
	next if($self->{SEEN}->[$self->{COUNT}]);

	#step forward in jumps if indicated
	if($self->{step} && $self->{step} != 1){
	    next if($self->{COUNT} < $self->{step} || $self->{COUNT} % $self->{step} != 0);
	}

	$self->{SEEN}->[$self->{COUNT}]++;

	$line =~ /^>([^\s\n]+)/;
	next if(defined ($self->{SKIP}{$1}));

	return $line;
    }

    #end of file but I was jumping using a step, so go back to start
    if($self->{step} && $self->{step} != 1){
	$self->{step} = undef; #remove step
	$self->{COUNT} = -1; #reset count
	$fh->setpos($self->startPos); #go to start of file
	$self->skip_file($self->{skip_file}) if($self->{skip_file} && -f $self->{skip_file}); #reload skip file
	return $self->nextDef();
    }
    
    $fh->close();
    return undef;
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


