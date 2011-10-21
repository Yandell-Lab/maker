#------------------------------------------------------------------------
#----           Bio::Search:HSP::PhatHSP::protein2genome                   ---- 
#------------------------------------------------------------------------
package Bio::Search::HSP::PhatHSP::protein2genome;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Bio::Search::HSP::PhatHSP::Base;
@ISA = qw(
	Bio::Search::HSP::PhatHSP::Base
       );


#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class  = shift;
	my @args   = @_;

	my $self = $class->SUPER::new(@args);

	bless ($self, $class);

	return $self;
}
#-------------------------------------------------------------------------------
sub strand {
	my $self = shift;
	my $w    = shift;

	$w = 'hit' if $w eq 'sbjct';

	return $self->{_strand_hack}->{$w};
}
#-------------------------------------------------------------------------------
sub show {
	my $self = shift;

        print "--------------------------------------\n";
        print "Bio::Search:HSP::PhatHSP::protein2genome\n";
	print "   ".$self->hit->seq_id()."\n";
        print "--------------------------------------\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDq:".$self->frac_identical('query')."\n";
	print "fIDt:".$self->frac_identical('total')."\n";
}
#-------------------------------------------------------------------------------
sub nB {
        my $self = shift;
        my $w    = shift;

        die unless defined($w);

        if ($self->strand($w) == 1){
                return $self->start($w);
        }
        else {
                return $self->end($w);
        }
}
#-------------------------------------------------------------------------------
sub nE {
        my $self = shift;
        my $w    = shift;

        if ($self->strand($w) == 1){
                return $self->end($w);
        }
        elsif ($w eq 'query' && $self->strand($w) == 0){
                return $self->end($w);
        }

        else {
                return $self->start($w);
        }
}
#-------------------------------------------------------------------------------
sub cigar_string {
    my $self = shift;
    
    return($self->{_CIGAR}) if($self->{_CIGAR});

    my $cigar = '';
    
    my @nt = $self->query_string() =~ /(.)/g;
    my @aa = $self->hit_string() =~ /(.)/g;
    
    die "ERROR: query and hit string lengths do not match correctly\n".
	"in Bio::Search:HSP::PhatHSP::protein2genome\n".
	"for hit ".$self->name()."\n"  if(@nt != @aa);
    
    my $type = ''; # M, I, R, F, R
    my $value = 0;
    
    my $q_buf ='';
    my $h_buf = '';
    
    for(my $i = 0; $i < @aa; $i++){
	#starts with a capital letter so process string buffer
	if($aa[$i] =~ /A-Z/ && $q_buf){
	    my $q_gaps = @{[$q_buf =~ /([^A-Za-z])/g]};
	    my $h_gaps = @{[$h_buf =~ /([^A-Za-z])/g]};
	    
	    my $found = ''; #what the buffer says should be done
	    my $step = 0; # what the step value is
	    
	    #M
	    if($q_gaps == 0 && $h_gaps == 0){
		$found = 'M';
		$step = 1;
	    }
	    #I
	    elsif($q_gaps == 3){
		$found = 'I';
		$step = 1;
	    }
	    #D
	    elsif($h_gaps == 3){
		$found = 'D';
		$step = 1;
	    }
	    #R
	    elsif($q_gaps){
		$found = 'R';
		$step = $q_gaps;
	    }
	    #F
	    elsif($h_gaps){
		$found = 'F';
		$step = $h_gaps;
	    }
	    
	    #cigar buffer was reverse frame shift which simultaneously itterates the match
	    if($type eq 'R'){
		$cigar .= $type.$value if($value && $type);
		$type = 'M';
		$value = 1;
	    }
	    
	    #cigar type different than current so process cigar buffer
	    if($found ne $type || $type eq 'F'){
		$cigar .= $type.$value if($value && $type);
		$type = $found;
		$value = 0;
	    }
	    
	    $value += $step;
	    
	    $q_buf ='';
	    $h_buf = '';
	}
	
	#add character to buffer
	$h_buf .= $aa[$i];		
	$q_buf .= $nt[$i];		
	
	#cooresponds to one amino acid or is last character so process string buffer
	if(length($h_buf) == 3 || $i == @aa - 1){
	    my $q_gaps = @{[$q_buf =~ /([^A-Za-z])/g]};
	    my $h_gaps = @{[$h_buf =~ /([^A-Za-z])/g]};
	    
	    my $found = ''; #what the buffer says should be done
	    my $step = 0; # what the step value is
	    
	    #M
	    if($q_gaps == 0 && $h_gaps == 0){
		$found = 'M';
		$step = 1;
	    }
	    #I
	    elsif($q_gaps == 3){
		$found = 'I';
		$step = 1;
	    }
	    #D
	    elsif($h_gaps == 3){
		$found = 'D';
		$step = 1;
	    }
	    #R
	    elsif($q_gaps){
		$found = 'R';
		$step = $q_gaps;
	    }
	    #F
	    elsif($h_gaps){
		$found = 'F';
		$step = $h_gaps;
	    }
	    
	    #cigar buffer was reverse frame shift which simultaneously itterates the match
	    if($type eq 'R'){
		$cigar .= $type.$value if($value && $type);
		$type = 'M';
		$value = 1;
	    }
	    
	    #cigar type different than current so process cigar buffer
	    if($found ne $type || $type eq 'F'){
		$cigar .= $type.$value if($value && $type);
		$type = $found;
		$value = 0;
	    }
	    
	    $value += $step;
	    
	    $q_buf ='';
	    $h_buf = '';
	    
	    #this was the last character so process cigar buffer
	    if($i == @aa - 1){
		$cigar .= $type.$value if($value && $type);
	    }
	}
    }
    
    if($self->strand('query') eq '-1'){
       my @set = $cigar =~ /([A-Z]\d+)/g;
       @set = reverse(@set);
       $cigar = join('', @set);
    }

    return ($self->{_CIGAR} = $cigar);
}
#-------------------------------------------------------------------------------
sub hasRun {

        print " method hasRun is not available for ";
        print " Bio::Search:HSP::PhatHSP::protein2genome\n";

        die;

}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Bio::Search:HSP::PhatHSP::protein2genome::AutoLoader called for: ",
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


