#-------------------------------------------------------------------------------
#------                            Shadower                            ---------
#-------------------------------------------------------------------------------
package Shadower;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Bit::Vector;
use Carp;

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
sub shadowVector {
	my $max      = shift;
	my $features = shift;
	my $flank    = shift;

	$max++; #I pretend coordinate 0 doesn't exist

	$flank = 0 unless defined($flank);

	my $vector = Bit::Vector->new($max);
        foreach my $p (@{$features}){
                my $b = $p->[0];
                my $e = $p->[1];

                ($b, $e) = ($e, $b) if $e < $b;

		my $f = $b - $flank;

                $b = $f > 0 ? $b - $flank : 1;
                $e = $e + $flank;
		
                $vector->Interval_Fill($b,$e);
        }

	return $vector;
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
                my $b = $e - length($1) + 1;
                push(@pieces, {b => $b , e => $e, piece => $1});
        }
        return \@pieces;

}
#-------------------------------------------------------------------------------
sub getPieces {
        my $sequence = shift;
        my $features = shift; #coordinates not objects
        my $flank    = shift;

	$flank = 0 unless defined($flank);

	my $pieces = getVectorPieces($features, $flank);

	if(ref($sequence) eq 'SCALAR'){
	    foreach my $p (@$pieces){
		my $seq = substr($$sequence, $p->{b}-1, ($p->{e}-$p->{b})+1);
		$p->{piece} = $seq;
	    }
	}
	else{
	    foreach my $p (@$pieces){
		my $seq = $sequence->subseq($p->{b}, $p->{e});
		$p->{piece} = $seq;
	    }
	}

	return $pieces;
}
#-------------------------------------------------------------------------------
sub getVectorPieces {
        my $features = shift; #coordinates not objects
        my $flank    = shift;

	if(ref($features) ne 'ARRAY'){
	    confess "ERROR: Not an array reference - ref: ".ref($features)."\n";
	}
	return [] if(!@$features);
	$flank = 0 unless defined($flank);

	foreach my $f (@$features){
	    ($f->[0], $f->[1]) = ($f->[1], $f->[0]) if($f->[0] > $f->[1]);
	}

	#initialize pieces array
	my @pieces;
	my @fs = sort {$a->[0] <=> $b->[0]} @$features;
	push(@pieces, {b => $fs[0][0]-$flank, e => $fs[0][1]+$flank});
	$pieces[-1]{b} = 1 if($pieces[-1]{b} < 1);

	#now grow pieces
	for(my $i = 1; $i < @fs; $i++){
	    my $p = $pieces[-1];
	    my $f = $fs[$i];

	    #FYI this is always true --> $p->{b} <= $f->[0]
	    if($f->[0] - $flank <= $p->{e} + 1){ #plus 1 because being neigbors is sufficient
		$p->{e} = $f->[1]+$flank if($f->[1]+$flank > $p->{e});
	    }
	    else{
		push(@pieces, {b => $f->[0]-$flank, e => $f->[1]+$flank});
		$pieces[-1]{b} = 1 if($pieces[-1]{b} < 1);
	    }
	}

	#my $top = 0;
	#foreach my $f (@$features){
	#    $top = $f->[1]+$flank if($f->[1]+$flank > $top);
	#    $top = $f->[0]+$flank if($f->[0]+$flank > $top);
	#}

        #my $sVec = shadowVector($top, $features, $flank);

        #my @pieces;
	#my $s = 0;
	#while (($s < $sVec->Size()) && (my ($min,$max) = $sVec->Interval_Scan_inc($s))){
	#    $s = $max + 2;
        #
	#    push(@pieces, {b => $min , e => $max});
	#}
        #

	return \@pieces;
}
sub getVectorPieces_old {
        my $features = shift; #coordinates not objects
        my $flank    = shift;

	if(ref($features) ne 'ARRAY'){
	    confess "ERROR: Not an array reference - ref: ".ref($features)."\n";
	}
	return [] if(!@$features);
	$flank = 0 unless defined($flank);


	my $top = 0;
	foreach my $f (@$features){
	    $top = $f->[1]+$flank if($f->[1]+$flank > $top);
	    $top = $f->[0]+$flank if($f->[0]+$flank > $top);
	}

        my $sVec = shadowVector($top, $features, $flank);

        my @pieces;
	my $s = 0;
	while (($s < $sVec->Size()) && (my ($min,$max) = $sVec->Interval_Scan_inc($s))){
	    $s = $max + 2;
        
	    push(@pieces, {b => $min , e => $max});
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
                my $b = $e - length($1) + 1;
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
