#-------------------------------------------------------------------------------
#------                            Fasta                               ---------
#-------------------------------------------------------------------------------
package Fasta;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
	my $class = shift;
	my $fasta = shift;

	my $self = {};

	bless $self, $class;

	$self->fasta($fasta);

	return $self;
}
#-------------------------------------------------------------------------------
sub identifier {
	my $self = shift;
	my $id   = shift;

        if (defined($id)){
                $self->{identifier} = $id;
        }
        elsif (defined($self->{identifier})){
                return $self->{identifier};
        }
	else {
		my ($id) = $self->def() =~ /^>(\S+)/;
		$self->{identifier} = $id;
		return $self->{identifier};
	}
}
#-------------------------------------------------------------------------------
sub fasta {
	my $self  = shift;
	my $fasta = shift;

	if (defined($fasta)){
		$self->{fasta} = $fasta;
	}
	else {
		return $self->{fasta};
	}
}
#-------------------------------------------------------------------------------
sub seq {
        my $self = shift;
        my $seq  = shift;

        if  (defined($seq) && $seq eq 'revComp' && !defined($self->{seq})){
                my @a = split(/\n/, $self->fasta());
                shift(@a);
                $self->{seq} = join('', @a);;
                return revComp($self->{seq});
        }
        if (defined($seq) && $seq eq 'revComp'){
                return revComp($self->{seq});
        }
        elsif (defined($seq)){
                $self->{seq} = $seq;
        }
        elsif (defined($self->{seq})){
                return $self->{seq};
        }
        else {
                my @a = split(/\n/, $self->fasta());
		shift(@a);
                $self->{seq} = join('', @a);;
                return $self->{seq};
        }
}
#-------------------------------------------------------------------------------
sub def {
	my $self = shift;
	my $def  = shift;

        if (defined($def)){
                $self->{def} = $def;
        }
	elsif (defined($self->{def})){
		return $self->{def};
	}
        else {
		my @a = split(/\n/, $self->fasta());
		$self->{def} = shift(@a);
		return $self->{def};
        }
}
#-------------------------------------------------------------------------------
#------------------------------- SUBS ------------------------------------------
#-------------------------------------------------------------------------------
sub  getWantedFromMulti {
        my $multiFasta = shift;
        my $wanted     = shift;

        $/ = "\n>";

        my $fh = new FileHandle;
           $fh->open("$multiFasta");

	my @fastas;
        while(my $line = <$fh>){
                $line =~ s/>//;
                $line = ">".$line;
		if (!defined($wanted) || $line =~ /$wanted/){
			push(@fastas, $line);
                }
        }
        $/ = "\n";
        $fh->close;

	return \@fastas;
}
#-----------------------------------------------------------------------------
sub getDef {
        my $fasta = shift;

	my @fasta = split(/\n/, $fasta);
        my $seq = '';
        while(my $l = shift(@fasta)){
                chomp($l);
                return $l if $l =~ /^>/;
        }

}
#-------------------------------------------------------------------------------
sub getSeq {
        my $fasta = shift;

	my @fasta = split(/\n/, $fasta);

        my $seq = '';
        while(my $l = shift(@fasta)){
                chomp($l);
                next if $l =~ /^>/;
                $seq .= $l;
        }

        return \$seq;

}
#-------------------------------------------------------------------------------
sub revComp {
	my $seq = shift;

        my($rc);

        $rc = $seq;
        $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX
                 /tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
        $rc = reverse scalar ($rc);

        return($rc);

}
#-------------------------------------------------------------------------------
sub toFasta {
        my $def = shift;
        my $seq = shift;

        my $fasta = $def."\n";

	$fasta .=  formatSeq($seq, 60);
        return \$fasta;

}
#-------------------------------------------------------------------------------
sub formatSeq {
        my $seq = shift;
	my $l   = shift;

	my $fasta = '';
        for (my $i=0; $i< length($$seq);$i+=$l){
                $fasta .= substr($$seq, $i, $l)."\n";
        }
        return $fasta;

}
#-------------------------------------------------------------------------------
sub toBpos {
        my $def = shift;
        my $seq = shift;

        my $fasta = $def."\n";

	my $bpos = "0 " x length($$seq);

        for (my $i=0; $i< length($bpos);$i+=60){
                $fasta .= substr($bpos, $i, 60)."\n";
        }
        return \$fasta;

}

#-------------------------------------------------------------------------------
sub toQual {
	my $def = shift;
	my $seq = shift;

	$$seq  =~ s/^\s+//;

	my @values = split(/\s+/, $$seq);
        my $fasta  = $def."\n";

	my $j = 0;
        for (my $i=0; $i< @values ;$i++){
		if ($j < 20){
			my $v = $values[$i];

			$fasta .= length($values[$i]) == 1 ? $values[$i]."  " 
			                                   : $values[$i]." ";
			$j++;
		}
		elsif ($j == 20){
			$fasta .= $values[$i]."\n";
			$j = 0;
		}
		else {
			die "dead in toQual\n";
		}
        }
	$fasta .= "\n" unless $fasta =~ /\n$/; 

        return \$fasta;

}
#-------------------------------------------------------------------------------
sub shift2 (\@) {
	return splice(@{$_[0]}, 0, 2);
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Fasta::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#----------------------------------------------------------------------------

1;
