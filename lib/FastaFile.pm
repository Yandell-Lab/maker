#-------------------------------------------------------------------------------
#------                            FastaFile                           ---------
#-------------------------------------------------------------------------------
package FastaFile;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;

@ISA = qw(
          );

#-------------------------------------------------------------------------------
#------------------------------- SUBS ------------------------------------------
#-------------------------------------------------------------------------------
sub  getWantedFromMulti {
        my $multiFasta = shift;
        my $wanted     = shift;

        $/ = "\n>";

        my $fh = new FileHandle;
           $fh->open("$multiFasta") || die "couldn't open $multiFasta\n$!";

	my @fastas;
        while(my $line = <$fh>){
                $line =~ s/>//;
                $line = ">".$line;
		if (!defined($wanted) || $line =~ /$wanted/){
			push(@fastas, $line);
                }
        }
        $/ = "\n";
        $fh->close || die "couldn't close $multiFasta\n$!";

	return \@fastas;
}
#-----------------------------------------------------------------------------
sub getFasta {
	my $fastaFile = shift;

	my $seq = getSeq($fastaFile);
	my $def = getDef($fastaFile);	

	return toFasta($def, $seq);
}
#-----------------------------------------------------------------------------
sub writeFile {
	my $f     = shift;
	my $loc   = shift;

	my $fasta = (ref($f) eq '') ? \$f : $f; 

        my $fh = new FileHandle();
           $fh->open(">$loc") || die "couldn't open $loc\n$!";

	print $fh $$fasta;

	$fh->close() || die "couldn't close $loc\n$!";

	return $loc;
}
#-----------------------------------------------------------------------------
sub getDef {
        my $fastaFile = shift;

        my $fh = new FileHandle();
           $fh->open($fastaFile) || die "couldn't open $fastaFile\n$!";

        my $seq = '';
        while(my $l = <$fh>){
                chomp($l);
                return $l if $l =~ /^>/;
	    }
        $fh->close() || die "couldn't close $fastaFile\n$!";

}
#-------------------------------------------------------------------------------
sub getName {
        my $fastaFile = shift;

        my $fh = new FileHandle();
           $fh->open($fastaFile) || die "couldn't open $fastaFile\n$!";

        my $seq = '';
        while(my $l = <$fh>){
                chomp($l);
		return $1 if $l =~ /^>(\S+)/;
        }
        $fh->close() || die "couldn't close $fastaFile\n$!";

}
#-------------------------------------------------------------------------------
sub getSeq {
	my $fastaFile = shift;

	my $fh = new FileHandle();
	   $fh->open($fastaFile) || die "couldn't open $fastaFile\n$!";

	my $seq = '';
	while(my $l = <$fh>){
		chomp($l);
		next if $l =~ /^>/;
		$seq .= $l;
	}
	$fh->close() || die "couldn't close $fastaFile\n$!";

	return \$seq;
	
}
#-------------------------------------------------------------------------------
sub revComp {
	my $seq = shift;

	$seq =~ tr/ACGTYRKMB/TGCARYMKV/;

	die "fix this! dead in FastaFile::revComp\n";
	return reverse($seq);
}
#-------------------------------------------------------------------------------
sub toFasta {
        my $def = shift;
        my $seq = shift;

        my $fasta = $def."\n";

        for (my $i=0; $i< length($$seq);$i+=60){
                $fasta .= substr($$seq, $i, 60)."\n";
        }
        return \$fasta;

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
1;









