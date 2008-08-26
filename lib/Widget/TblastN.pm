#------------------------------------------------------------------------
#----                        Widget::tblastn                         ---- 
#------------------------------------------------------------------------
package Widget::tblastn;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Bio::DB::Fasta;
@ISA = qw(
	Widget
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
#------------------------------------------------------------------------------
sub run {
        my $self    = shift;
        my $command = shift;

        my $exe       = '/usr/local/bdgp/wublast/tblastn';
        my $mat       = '/usr/local/bdgp/wublast/matrix/aa/BLOSUM80';
        my $db        = $self->db();
        my $projDir   = $self->projectDir();
        my $fastaFile = $self->queryFastaFile();

        my $out = "$projDir/blastReport";

        if (defined($command)){
                $self->print_command($command);
                system("$command");
        }
        else {

		system("$exe $db  $fastaFile  -V10 -B10 -topComboN=1 wordmask=seg > $out");
                $self->blastReport($out);
        }
}
#------------------------------------------------------------------------------
sub parse {
	my $self   = shift;
	my $params = shift;

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

        print STDERR "Widget::BlastX::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


