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
sub run {
	my $self    = shift;
	my $command = shift;

	if (defined($command)){
		$self->print_command($command);
		system("$command");
	}
	else {

		die "No command given in Widget::qrna::run\n";
	}
}
#-------------------------------------------------------------------------------
sub parse {
	my $self = shift;
	my $file = shift;

	$/ = 'length alignment';

	if (! -e $file){
		print STDERR "$file does not exist!\n";
		return;
	}
	my $fh = new FileHandle();
	   $fh->open($file);

	my $i = 0;

	my $header = <$fh>;

	my %data;
	$data{run_data}{out_file} = $file;

	parse_header(\%data, $header);

	while (my $line = <$fh>){

		chomp($line);

		my $segment = parse_segments($line);

		push(@{$data{segments}}, $segment);
	}


	$/ = "\n";

	$fh->close();

	$self->{results} = \%data;
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse_header {
	my $data   = shift;
	my $header = shift;

	my @lines = split(/\n/, $header);

	while (defined (my $line = shift(@lines))){
		chomp($line);
		$line =~ s/^\#//;
		$line =~ s/^\s+//;
		if    ($line =~ /^qrna/){
			$data->{run_data}->{version} = $line;
		}
		elsif ($line =~ /Rate-generating PAM model/){
			($data->{run_data}->{PAM}) = 
			$line =~ /Rate-generating PAM model =  (\S+)/;
		}
		elsif ($line =~ /Rate-generating RIBOPROB matrix/){
			($data->{run_data}->{RIBOPROB}) = 
			$line =~ /Rate-generating RIBOPROB matrix =  (\S+)/;
                }
                elsif ($line =~ /seq file/){
			($data->{run_data}->{seq_file}) = 
			$line =~ /seq file  =  (\S+)/;
                }
                elsif ($line =~ /\#seqsi\:/){
			$data->{run_data}->{num_seqs} = $line;
                }
                elsif ($line =~ /window version/){
			($data->{run_data}->{window_size}) = 
			$line =~ /window version:\s+(.*)$/;
                }
                elsif ($line =~ />/){
			$line =~ s/>//;
			push(@{$data->{run_data}->{seq_ids}}, $line);
                }
                elsif ($line =~ /length of whole alignment/){
			($data->{run_data}->{tot_align_length}) =
			$line =~ /length of whole alignment after removing common gaps:\s+(\d+)/;
                }
                elsif ($line =~ /Divergence time/){
			($data->{run_data}->{divergence_times}) =
			$line =~ /Divergence time \(variable\)\:\s+(.*)$/;
                }
                elsif ($line =~ /\[alignment/){
			my ($id, $mt, $gp) = 
			$line =~ /alignment\s+ID\s+=\s+(\S+)\s+MUT\s+=\s+(\S+)\s+GAP\s+=\s+(\S+)$/;
			$data->{run_data}->{alignment_stats}->{id} = $id;
			$data->{run_data}->{alignment_stats}->{mt} = $mt;
			$data->{run_data}->{alignment_stats}->{gp} = $gp;
                }



	}

}
#-------------------------------------------------------------------------------
sub parse_segments {
	my $segment = shift;

        my @lines = split(/\n/, $segment);

	my %seg;
        while (defined (my $line = shift(@lines))){
		chomp($line);
		$line =~ s/^\s+//;

		if    ($line =~ /gap\=\d+/){
			my ($l, $i, $m, $g) = $line =~/:\s+(\d+)\s+\(id=(\S+)\)\s+\(mut=(\S+)\)\s+\(gap=(\S+)\)/;
			$seg{align_stats}{length} = $l;
			$seg{align_stats}{id}     = $i;	
			$seg{align_stats}{mt}     = $m;
			$seg{align_stats}{gp}     = $g;
		}
		elsif ($line =~ /^pos[XY]/){
			my ($p, $b, $e) = $line =~ /pos([XY])\:.*\[(\d+)\-(\d+)\]\(.*/;
			$seg{pos}{$p} = [$b, $e];
		}
                elsif ($line =~ /^SS\s+[\.\1]/){
			chomp($line);
			$line =~ s/^\s+//;
			
			my $a = $line;
			   $a =~ s/^\s+//;
			   $a =~ s/SS//;
			my $b = shift(@lines);
			   $b =~ s/SS//;
			   $b =~ s/^\s+//;
	
			my $x = shift(@lines);
		        my ($x_id, $x_str) = $x =~ /(\S+)\s+(\S+)/;

			my $y = shift(@lines);
			my ($y_id, $y_str) = $y =~ /(\S+)\s+(\S+)/;

			$seg{alignment}{SS1} .= $a;
			$seg{alignment}{SS2} .= $b;
			$seg{alignment}{X}   .= $x_str;
			$seg{alignment}{Y}   .= $y_str;
			$seg{alignment}{X_ID} = $x_id;
			$seg{alignment}{Y_ID} = $y_id;
                }
                elsif ($line =~ /\w{3}\s+ends\s+\*/){
			my ($mod) =  $line =~ /^(\w{3}).*/;
			my ($s)   =  $line =~ /\(([\+\-])\).*/;
			my ($off) =  $line =~ /\=\s+\((\d+)\.+/;
			my ($len) =  $line =~ /\[(\d+)\]/;

			$seg{MODELS}{$mod}{strand} = $s;
			$seg{MODELS}{$mod}{offset} = $off;
			$seg{MODELS}{$mod}{length} = $len;
                }
                elsif ($line =~ /^OTH\s\=\s+[\d\.]+\s+COD.+/){
			my ($oth) = $line =~ /OTH\s\=\s+(\-?[\d\.]+)/;
			my ($cod) = $line =~ /COD\s\=\s+(\-?[\d\.]+)/;
			my ($rna) = $line =~ /RNA\s\=\s+(\-?[\d\.]+)/;

			$seg{MODELS}{OTH}{raw_score} = $oth;
			$seg{MODELS}{COD}{raw_score} = $cod;
			$seg{MODELS}{RNA}{raw_score} = $rna;
                }
                elsif ($line =~ /^logoddspost/){
                        my ($oth) = $line =~ /logoddspostOTH\s\=\s+(\-?[\d\.]+)/;
                        my ($cod) = $line =~ /logoddspostCOD\s\=\s+(\-?[\d\.]+)/;
                        my ($rna) = $line =~ /logoddspostRNA\s\=\s+(\-?[\d\.]+)/;

                        $seg{MODELS}{OTH}{lod_post_p} = $oth;
                        $seg{MODELS}{COD}{lod_post_p} = $cod;
                        $seg{MODELS}{RNA}{lod_post_p} = $rna;

                }
                elsif ($line =~ /^sigmoidal/){
                        my ($oth) = $line =~ /sigmoidalOTH\s\=\s+(\-?[\d\.]+)/;
                        my ($cod) = $line =~ /sigmoidalCOD\s\=\s+(\-?[\d\.]+)/;
                        my ($rna) = $line =~ /sigmoidalRNA\s\=\s+(\-?[\d\.]+)/;

                        $seg{MODELS}{OTH}{sigmoidal_score} = $oth;
                        $seg{MODELS}{COD}{sigmoidal_score} = $cod;
                        $seg{MODELS}{RNA}{sigmoidal_score} = $rna;

                }
                elsif ($line =~ /^winner/){
			($seg{winning_model}) = $line =~ /^winner\s+=\s+(\w{3})/;
                }


	}
	return \%seg;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::qrna::AutoLoader called for: ",
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


