#------------------------------------------------------------------------
#----                        Widget::geneid                          ---- 
#------------------------------------------------------------------------
package Widget::geneid;
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

	my $exe       = '/users/myandell/bin/geneid/geneid.Oct_09_2002';
	my $projDir   = $self->projectDir();
	my $fastaFile = $self->queryFastaFile();

	my $out = "$projDir/geneidReport";

	if (defined($command)){
		$self->print_command($command);
		system("$command");
	}
	else {

		$self->print_command();
 		system("$exe -G  $fastaFile > $out");
		$self->report($out);
	}
}
#-------------------------------------------------------------------------------
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

        print STDERR "Widget::geneid::AutoLoader called for: ",
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


