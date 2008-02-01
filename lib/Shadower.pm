#-------------------------------------------------------------------------------
#------                            Shadower                            ---------
#-------------------------------------------------------------------------------
package Shadower;
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
sub shadowSequence {
	my $sequence = shift;
	my $features = shift;
	my $flank    = shift;

	$flank = 0 unless defined($flank);

	$sequence = \lc($$sequence);
        foreach my $p (@{$features}){
                my $b = $p->[0];
                my $e = $p->[1];

                ($b, $e) = ($e, $b) if $e < $b;

		my $f = $b - $flank;

                $b = $f > 0 ? $b - $flank : 1;
                $e = $e + $flank;

                my $l = $e - $b + 1;
		
                substr($$sequence, $b -1 , $l, uc(substr($$sequence, $b -1 , $l)));
        }

	return $sequence;
}
#-------------------------------------------------------------------------------
sub reverseMaskSequence {
        my $sequence = shift;
        my $features = shift;
        my $flank    = shift;
        my $char     = shift;

        $char = 'N' unless defined($char);

        $flank = 0 unless defined($flank);

        my $sSeq = shadowSequence($sequence, $features, $flank);

        $$sSeq =~ s/[a-z]/$char/g;

        return \uc($$sSeq);

}
#-------------------------------------------------------------------------------
sub softMaskSequence {
        my $sequence = shift;
        my $features = shift;
        my $flank    = shift;

	print STDERR "softmasking can only be done once AND DONE LAST!\n";
	print STDERR "this function is alpha!\n";

	$flank = 0 unless defined($flank);

	my $sSeq = shadowSequence($sequence, $features, $flank);

	$$sSeq =~ tr/a-zA-Z/A-Za-z/;

	return $$sSeq;

}
#-------------------------------------------------------------------------------
sub maskSequence {
        my $sequence = shift;
        my $features = shift;
        my $flank    = shift;
	my $char     = shift;

	$char = 'N' unless defined($char);

	$flank = 0 unless defined($flank);

	my $sSeq = shadowSequence($sequence, $features, $flank);

	$$sSeq =~ s/[A-Z]/$char/g;

	return \uc($$sSeq);
}
#-------------------------------------------------------------------------------
sub getUpperCasedSegments {
        my $sequence = shift;
        my $flank    = shift;

        $flank = 0 unless defined($flank);

        my @pieces;
        while ($$sequence =~ m/([A-Z]+)/g ) {
                my $e = pos($$sequence);
                my $b = $e - length($1);
                push(@pieces, {b => $b , e => $e, piece => $1});
        }
        return \@pieces;

}
#-------------------------------------------------------------------------------
sub getPieces {
        my $sequence = shift;
        my $features = shift;
        my $flank    = shift;

	$flank = 0 unless defined($flank);

        my $sSeq = shadowSequence($sequence, $features, $flank);

        my @pieces;
        while ($$sSeq =~ m/([A-Z]+)/g ) {
                my $e = pos($$sSeq);
                my $b = $e - length($1);
                push(@pieces, {b => $b , e => $e, piece => $1});
        }
	return \@pieces;
}
#-------------------------------------------------------------------------------
sub getNegativePieces {
        my $sequence = shift;
        my $features = shift;
        my $flank    = shift;

        $flank = 0 unless defined($flank);

        my $sSeq = shadowSequence($sequence, $features, $flank);


        my @pieces;
        while ($$sSeq =~ m/([a-z]+)/g ) {
                my $e = pos($$sSeq);
                my $b = $e - length($1);
                push(@pieces, {b => $b , e => $e, piece => $1});
        }
        return \@pieces;
}
#-------------------------------------------------------------------------------
sub shift2 (\@) {
	return splice(@{$_[0]}, 0, 2);
}
#-------------------------------------------------------------------------------
1;
