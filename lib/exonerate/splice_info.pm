#------------------------------------------------------------------------
#----                       exonerate::splice_info                   ---- 
#------------------------------------------------------------------------
package exonerate::splice_info;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use CGL::TranslationMachine;
use FastaSeq;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub splice_code {
	my $d = shift;
	my $a = shift;

	my $str = '';
        if    ($d =~ /^gt$/i && $a =~ /^ag$/i) {
       		$str = '+';
        }
        elsif ($d =~ /^ct$/i && $a =~ /^ac$/i) {
                $str = '-';
        }
        elsif ($d =~ /^at$/i && $a =~ /^ac$/i){
               $str = '+';
        }
        elsif ($d =~ /^gt$/i && $a =~ /^at$/i){
               $str = '-';
        }

        else {
              $str = '0';
        }

	return $str;
}
#------------------------------------------------------------------------
sub get_numbers {
	my $splice_str = shift;

	my $plus  = $splice_str =~ tr/\+/\+/;
	my $minus = $splice_str =~ tr/\-/\-/;
	my $zero  = $splice_str =~ tr/0/0/;

	my $length = length($splice_str);

	my %data;

	$data{p} = $plus;
	$data{m} = $minus;
	$data{z} = $zero;
	$data{l} = $length;

	return \%data;
}
#------------------------------------------------------------------------
sub needs_to_be_revcomped {
	my $hit = shift;
	my $seq_o = shift;

	if ($hit->num_hsps == 1){
	    my $seq;
	    my $r_seq;
	    
	    if($seq_o){
		my $B = [$hit->hsps()]->[0]->nB('query') - 1; #make array space value
		my $E = [$hit->hsps()]->[0]->nE('query') - 1; #make array space value
		my $strand = $hit->strand('query');
		my $length = abs($E - $B) +1; #substr length
		$seq = FastaSeq::substr_o($seq_o, $B, $length);
		$r_seq = Fasta::revComp($seq);
		($seq, $r_seq) = ($r_seq, $seq) if($strand == -1);
	    }
	    else{
		$seq   = [$hit->hsps()]->[0]->seq('query')->seq();
		$r_seq = Fasta::revComp(\$seq);
	    }
	    
	    $seq =~ s/[\-\_\+\=\|]//g;
	    $r_seq =~ s/[\-\_\+\=\|]//g;
	    
	    my $tM = new CGL::TranslationMachine();
	    (my $p_seq , undef) = $tM->longest_translation($seq);
	    (my $r_p_seq , undef) = $tM->longest_translation($r_seq);
	    
	    my $p_len = length($p_seq);
	    my $r_len = length($r_p_seq);
	    if( $p_len >= $r_len){
		return 0;
	    }
	    else{
		return 1;
	    }
	}

	my $str = get_splice_str($hit);

	my $num = get_numbers($str);

	my $ratio = $num->{m}/$num->{l};

	if ($ratio >= .75){
		return 1;
	}
	else {
		return 0;
	}

}
#------------------------------------------------------------------------
sub get_splice_str {
    my $hit = shift;
    
    my $sorted = PhatHit_utils::sort_hits($hit);
    
    my $splice_str = '';
    for (my $i = 0; $i < @{$sorted} - 1; $i++){
	my $pre_hsp = $sorted->[$i];
	my $pos_hsp = $sorted->[$i+1];
		
	my $code = splice_code($pre_hsp->donor(),
			       $pos_hsp->acceptor(),
			       );
	$splice_str .= $code;
    }
    
    return $splice_str;
}
#------------------------------------------------------------------------
sub set_donors_acceptors {
    my $hit = shift;
    my $seq = shift;

    die "ERROR: No seq given in splice_info::set_donors_acceptors\n" if(! $seq);

    my $sorted = PhatHit_utils::sort_hits($hit);
    for (my $i = 1; $i < @{$sorted}; $i++){
        my $pre_hsp = $sorted->[$i-1];
        my $pos_hsp = $sorted->[$i];

	if(! exists $pre_hsp->{donor}){
	    my $strand = $pre_hsp->strand('query');
	    my $E = $pre_hsp->nE('query') - 1; #make array space value
	    my $length = 2; #substr length
	    my $p = ($strand == 1) ? $E + 1 : $E - $length; #substr start position
	    my $donor = FastaSeq::substr_o($seq, $p, $length);

	    $donor = Fasta::revComp($donor) if($strand == -1);
	    $pre_hsp->{donor} = $donor;
	}
	if(! exists $pos_hsp->{acceptor}){
	    my $strand = $pos_hsp->strand('query');
	    my $B = $pos_hsp->nB('query') - 1; #make array space value
	    my $length = 2; #substr length
	    my $p = ($strand == 1) ? $B - $length : $B + 1; #substr start position
	    my $acceptor = FastaSeq::substr_o($seq, $p, $length);
	    $acceptor = Fasta::revComp($acceptor) if($strand == -1);
	    $pos_hsp->{acceptor} = $acceptor;
	}
    }
} 
#------------------------------------------------------------------------
1;


