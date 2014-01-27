#------------------------------------------------------------------------
#----                        Widget::primer3                         ---- 
#------------------------------------------------------------------------
package Widget::primer3;
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
sub write_in_file {
	my $self     = shift;
	my $location = shift;

	my $fh = new FileHandle();
	   $fh->open(">$location");

	foreach my $d (@{$self->directives}){
		my ($key, $value) =  each %{$d};
		print $fh "$key\=$value\n";
	}
	print $fh "\=\n";
	$fh->close();
}
#------------------------------------------------------------------------------
sub run {
	my $self     = shift;
	my $in_file  = shift;
	my $out_file = shift;

	my $exe = '/usr/local/bin/primer3';

	my $command = "cat $in_file | $exe  > $out_file";

	$self->print_command($command);
	system("$command");

	$self->parse($out_file);

	system("rm $in_file");
	system("rm $out_file");
}
#-------------------------------------------------------------------------------
sub add_directive {
	my $self  = shift;
	my $key   = shift;
	my $value = shift;

	push(@{$self->{directives}}, {$key => $value});
}
#-------------------------------------------------------------------------------
sub directives {
	my $self = shift;

	return $self->{directives} || [];
}
#-------------------------------------------------------------------------------
sub _rearrange_results {
	my $self = shift;

	my $id_base;
	my %primers;
	while (my $datum = shift(@{$self->results})){
		my ($k, $v) = each %{$datum};
		if ($k =~ /PRIMER_(LEFT|RIGHT)/){
			if ($k =~ /PRIMER_(LEFT|RIGHT)_(\d+)_(\w+)/){
				$primers{$1}[$2]{$3} = $v;
			}
			elsif ($k =~ /PRIMER_(LEFT|RIGHT)_(\d+)$/) {
				my $one = $1;
				my $two = $2; 

				my ($o, $l) = $v =~ /(\d+)\,(\d+)/;
				$primers{$one}[$two]{'OFFSET'} = $o;
				$primers{$one}[$two]{'LENGTH'} = $l;
			}
                        elsif ($k =~ /PRIMER_(LEFT|RIGHT)$/){
				my $one = $1;
				my ($o, $l) = $v =~ /(\d+)\,(\d+)/;
				$primers{$one}[0]{'OFFSET'} = $o;
				$primers{$one}[0]{'LENGTH'} = $l;
                        }
			elsif ($k =~ /PRIMER_(LEFT|RIGHT)_(\w+)$/){
				$primers{$1}[0]{$2} = $v;
			}
		}
		else {
			$id_base = $v if $k eq 'PRIMER_SEQUENCE_ID';
		}
	}
	foreach my $key (keys %primers){
		for (my $i = 0; $i < @{$primers{$key}}; $i++){
			$primers{$key}[$i]{ID} = "$id_base\.$i";
		}
	}
	$self->primers(\%primers);
}
#-------------------------------------------------------------------------------
sub primer {
	my $self  = shift;
	my $up_dn = shift;
	my $i     = shift;

	return $self->primers->{$up_dn}->[$i];
}
#-------------------------------------------------------------------------------
sub primers {
	my $self    = shift;
	my $primers = shift;

	if (defined($primers)){
		$self->{primers} = $primers;
	}
	return $self->{primers} || [];
}
#-------------------------------------------------------------------------------
sub results {
	my $self = shift;

	return $self->{results} || [];
}
#-------------------------------------------------------------------------------
sub parse {
	my $self = shift;
	my $file = shift;

	my $fh = new FileHandle();
	   $fh->open("$file");

	my %results;
	while (my $line = <$fh>){
		chomp($line);
		last if $line eq '=';
		my @data = split(/\=/, $line);
		if ($data[0] eq 'PRIMER_ERROR'){
			push(@{$self->{errors}},  {$data[0] => $data[1]});
		}
		else {
			push(@{$self->{results}}, {$data[0] => $data[1]});
		}
	}
	$self->_rearrange_results();
	$fh->close();

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

        print STDERR "Widget::primer3::AutoLoader called for: ",
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


