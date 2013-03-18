#-------------------------------------------------------------------------------
#------                               Iterator::GFF3                     -------
#-------------------------------------------------------------------------------
package Iterator::GFF3;
use strict;
use Iterator;
use Iterator::Fasta;
use vars qw(@ISA @EXPORT $VERSION $AUTOLOAD);
use Exporter;
use FileHandle;
use Bio::Tools::GFF;
use Fasta;
use File::Temp qw(tempfile);
use Scalar::Util qw(openhandle);

@ISA = qw(Iterator::Fasta Iterator);

#-------------------------------------------------------------------------------
#------------------------------- FUNCTIONS -------------------------------------
#-------------------------------------------------------------------------------
sub new {
        my $class      = shift; #evaluator::gff3_to_phatHit::Maker
	my $gff_file   = shift; #gff file, can be for whole genome
	my $fasta_file = shift; #multi-fasta file for gff3

        my $self = {};
        bless ($self, $class);

	$fasta_file = $gff_file if (! $fasta_file);
        $self->fileName($fasta_file);
	$self->fileHandle($fasta_file);
	$self->{gff_file} = $gff_file;

	my $fh = $self->fileHandle();
	if (! openhandle($fh)){ #checks to see if file handle is open
	    die "ERROR: No open filehandle Iterator::GFF3\n"; 
	}

	if($gff_file eq $fasta_file){
	    my $line;
	    my $last = $fh->getpos();
	    while($line = <$fh>){
		if ($line =~ /^\#\#FASTA/){
		    $self->startPos($fh->getpos());
		    last;
		}
		elsif($line =~ /^\>/){
		    $self->startPos($last);
		    $fh->setpos($self->startPos);
		    last;
		}
		$last = $fh->getpos();
	    }
	}

	return $self;
}
#-------------------------------------------------------------------------------
sub find{
   my $self = shift;
   my $id = shift;

   my $fasta = $self->SUPER::find($id);
   my $features = $self->getFeatures($id);

   return ( $fasta, $features);
}
#-------------------------------------------------------------------------------
sub nextEntry{
   my $self = shift;

   my $fasta = $self->SUPER::nextEntry();
   my $id = Fasta::getSeqID(\$fasta);
   my $features = $self->getFeatures($id);

   return ($fasta, $features);
}
#-------------------------------------------------------------------------------
sub nextFasta{
   my $self = shift;

   return $self->SUPER::nextEntry();
}
#----------------------------------------------------------------------
sub getFeatures {
   my $self = shift;
   my $id  = shift;
   
   my $parser = new Bio::Tools::GFF(-gff_version => 3,
				    -file        => $self->{gff_file},
				   );
   
   my @features;
   while (my $f = $parser->next_feature()) {      
      push(@features, $f) if(! defined($id) || $f->seq_id eq $id);
   }
   
   return \@features;
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

        #print STDERR "Iterator::GFF3::AutoLoader called for: ",
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
