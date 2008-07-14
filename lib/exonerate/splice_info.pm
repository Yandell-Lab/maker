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
use SimpleCluster;
use compare;
use cluster;
use CGL::TranslationMachine;

@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub splice_code {
	my $d = shift;
	my $a = shift;

	my $str = '';
        if    ($d eq 'gt' && $a eq 'ag') {
       		$str = '+';
        }
        elsif ($d eq 'ct' && $a eq 'ac') {
                $str = '-';
        }
        elsif ($d eq 'at' && $a eq 'ac'){
               $str = '+';
        }
        elsif ($d eq 'gt' && $a eq 'at'){
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

	if ($hit->num_hsps == 1){
	        my $seq = [$hit->hsps()]->[0]->seq('query')->seq();
		$seq =~ s/[\-\_\+\=\|]//g;

		my $r_seq = [$hit->hsps()]->[0]->seq('query')->revcom()->seq;
		$r_seq =~ s/[\-\_\+\=\|]//g;
		
		my $tM = new CGL::TranslationMachine();
		(my $p_seq , undef) = $tM->longest_translation($seq);
		(my $r_p_seq , undef) = $tM->longest_translation($r_seq);

		if( length($p_seq) >= length($r_p_seq) ){
		        return 0 if $hit->strand('hit') ==  1;
		        return 1 if $hit->strand('hit') == -1;
		}
		else{
		        return 1 if $hit->strand('hit') ==  1;
		        return 0 if $hit->strand('hit') == -1;
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
	for (my $i = 1; $i < @{$sorted}; $i++){
		my $pre_hsp = $sorted->[$i-1];
		my $pos_hsp = $sorted->[$i];
 
		my $code = splice_code($pre_hsp->donor(),
		                       $pos_hsp->acceptor(),
				      );
		$splice_str .= $code;
	}

	return $splice_str;
}
#------------------------------------------------------------------------
1;


