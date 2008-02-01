#------------------------------------------------------------------------
#----                        Widget::qrna                            ---- 
#------------------------------------------------------------------------
package Widget::qrna;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
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
sub get_result {
	my $self = shift;
	my @args = @_;

	my $str ="\$self->{results}";
	foreach my $arg (@args){
		if ($arg =~ /\d+/){
			$str .= "->[$arg]";
		}
		else {
			$str .= "->{$arg}"; 	
		}
	}
	print STDERR "$str\n";
	return eval $str;
}
#------------------------------------------------------------------------------
sub run {
	my $self    = shift;
	my $command = shift;

	my $exe       = '/home/myandell/work/qrna/qrna-2.0.0/src/qrna';
	my $projDir   = $self->projectDir();
	my $fastaFile = $self->queryFastaFile();

	my $out = "$projDir/qrnaReport" if defined($projDir);

	if (defined($command)){
		$self->print_command($command);
		system("$command");
	}
	else {

		die;
		$self->print_command();
 		#system();
		$self->blastReport($out);
	}
}
#-------------------------------------------------------------------------------
sub parse {
	my $self = shift;
	my $file = shift;



	if (! -e $file){
		print STDERR "$file does not exist!\n";
		return;
	}
	my $fh = new FileHandle();
	   $fh->open($file);

	my %data;
	my $i = 0;
	while (my $line = <$fh>){

		chomp($line);
		if    ($line =~ s/#\s+ qrna\s+//){
			$data{version} = $line;
		}
		elsif ($line =~ s/#\s+PAM\s+model\s+=\s+//){
			$data{prot_matrix} = $line;
		}
		elsif ($line =~ s/#\s+RNA\s+model\s+=\s+//){
			$data{rna_model} = $line;	
		}
		elsif ($line =~ s/#\s+RIBOPROB\s+matrix\s+=\s+//){
			$data{ribo_matrix} = $line;
		}
		elsif ($line =~ s/#\s+seq\s+file\s+=\s+//){
			$data{seq_file} = $line;
		}
		elsif ($line =~ s/#\s+#seqs\:\s+//){
			#$data{seqs} = $line
		}
		elsif ($line =~ s/^>//){
			my ($name, $length) = $line =~ /(\S+)\s+\((\d+)\)/;
			push(@{$data{seqs}},{name => $name, length => $length});
			
		}
		elsif ($line =~ s/Divergence\s+time\s+\(variable\)\:\s+//){
			$data{divergence_time} = $line;
		}
		elsif ($line =~ s/^length alignment\:\s+//){
			my ($l, $i, $m, $g) = 
			$line =~ 
			/^(\d+)\s+\(id=(\S+)\)\s+\(mut=(\S+)\)\s+\(gap=(\S+)\)/; 

			$data{alignment_stats} = {length => $l, 
			                              id => $i, 
			                             mut => $m, 
			                             gap => $g,
			                         };

		}
		elsif ($line =~ /^\w{3}\s+ends/){
			my $w = $line =~ /\*/ ? 1: 0;

			my ($m) = $line =~ /^(\w{3})\s+ends\s+/;

			my ($s) = $line =~ /\(([+-])\)/;

			my ($o, $l) = $line =~ /=\s+\((\d+)\.\.\[(\d+)\]/;	

			push(@{$data{ends}{$m}}, {strand => $s,
			                          offset => $o,
			                          length => $l,
						  winner => $w,
			                          });
		}
		elsif ($line =~ s/^winner\s+=\s+//){
			$line =~ s/\s+//g;
			$data{winner} = $line; 
		}
		elsif ($line =~ /^\s+OTH\s+=\s+\S+\s+COD\s+=\s+\S+\s+RNA\s+=\s+\S+/){
			my ($o, $c, $r) = 
			$line=~/\s+OTH\s+=\s+(\S+)\s+COD\s+=\s+(\S+)\s+RNA\s+=\s+(\S+)/; 

			$data{scores}{raw}{OTH} = $o;
			$data{scores}{raw}{COD} = $c;
			$data{scores}{raw}{RNA} = $r;
		}
                elsif ($line =~ /^\s+logoddspostOTH/){

                        my ($o) = $line=~/\s+logoddspostOTH\s+=\s+(\S+)/;

			my ($c) = $line=~/\s+logoddspostCOD\s+=\s+(\S+)/;

			my ($r) = $line=~/\s+logoddspostRNA\s+=\s+(\S+)/;

                        $data{scores}{lod_post}{OTH} = $o;
                        $data{scores}{lod_post}{COD} = $c;
                        $data{scores}{lod_post}{RNA} = $r;
                }
                elsif ($line =~ /^\s+sigmoidalOTH/){

                        my ($o) = $line=~/\s+sigmoidalOTH\s+=\s+(\S+)/;

                        my ($c) = $line=~/\s+sigmoidalCOD\s+=\s+(\S+)/;

                        my ($r) = $line=~/\s+sigmoidalRNA\s+=\s+(\S+)/;

                        $data{scores}{sigmoidal}{OTH} = $o;
                        $data{scores}{sigmoidal}{COD} = $c;
                        $data{scores}{sigmoidal}{RNA} = $r;
                }


	}

	$fh->close();

	$self->{results} = \%data;

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

        print STDERR "Widget::blastp::AutoLoader called for: ",
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


