#------------------------------------------------------------------------
#----              evaluator::scoring.pm	               ---- 
#------------------------------------------------------------------------
package evaluator::scoring;
use strict;
use vars qw(@ISA @EXPORT);
use Exporter;
use FileHandle;

@ISA = qw(
       );
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------


#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub valuecal {

#This function calculates the C-value of a single single. The parameter that it 
#gets is a hash which has the following keys: QI (which is an array reference),
#type, junk, alt_splicing, and SO_CODE(which is an array reference).
	my %arg= @_;
        my $gene;
	$gene->{QI}= $arg{QI};
	$gene->{SO_CODE}=$arg{SO_CODE};
	$gene->{type}=$arg{type};
	$gene->{junk}=$arg{junk};
	$gene->{alt_splicing}=$arg{alt_splicing};
	$gene->{completion}=$arg{completion};
	
	if ( !defined($gene->{type}) || ($gene->{type} eq 'unknown')) 
		{ $gene->{type}='mRNA';}
	if ($gene->{type} eq 'protein-coding' || $gene->{type} eq 'Coding_transcript') 
		{ $gene->{type}='mRNA';}

        my ($c1_1,$c1_2,$c1_3,$c1_4,$c1_5,$c1_6,$c1_7,$c1_8,$c1_9);
        my ($alt,$junk,$bonus);

                if    (($gene->{type} eq "mRNA")
				&& ($gene->{QI}->[6]==1) ) {
                        $c1_1 =0; $c1_2=0  ; $c1_3=50; $c1_4=25; $c1_5= 0;
                        $c1_6 =30;$c1_7=0;   $c1_8=0;  $c1_9=0;  $alt = 30;
                        $junk =10;$bonus=60;
                }

                elsif (($gene->{type} eq "mRNA")
				&& ($gene->{QI}->[6]!=1) ) {
                        $c1_1 =0; $c1_2=100; $c1_3=50; $c1_4=25; $c1_5= 70;
                        $c1_6 =30;$c1_7=0;   $c1_8=0;  $c1_9=0;  $alt = 30;
                        $junk =10;$bonus=0;
                }

                elsif (($gene->{type} ne "mRNA")
				&& ($gene->{QI}->[6]==1) ) {
                        $c1_1 =0; $c1_2=0;   $c1_3=50; $c1_4=0 ; $c1_5=  0;
                        $c1_6 = 0;$c1_7=0;   $c1_8=0;  $c1_9=0;  $alt = 30;
                        $junk = 0;$bonus=120; 
                }

                elsif (($gene->{type} ne "mRNA")
				&& ($gene->{QI}->[6]!=1) ) {
                        $c1_1 =0; $c1_2=100; $c1_3=50; $c1_4=0;  $c1_5= 70;
                        $c1_6 = 0;$c1_7=0;   $c1_8=0;  $c1_9=0;  $alt = 30;
                        $junk = 0;$bonus=0;
                }


                my ($utr5,$utr3, $so_a, $so_c);
                if ($gene->{QI}->[0]==0) {$utr5=0;}
                else {$utr5=10;}
                if ($gene->{QI}->[7]==0) {$utr3=0;}
                else {$utr3=10;}
                if ($gene->{SO_CODE}->[0]!=0) {$so_a=-50;}
                else {$so_a=0;}
                if ($gene->{SO_CODE}->[2]!=0) {$so_c=50;}
                else {$so_c=0;}

                my $value   = $utr5          + $gene->{QI}->[1]*$c1_2 +
                      $gene->{QI}->[2]*$c1_3 + $gene->{QI}->[3]*$c1_4 +
                      $gene->{QI}->[4]*$c1_5 + $gene->{QI}->[5]*$c1_6 +
                      $gene->{QI}->[6]*$c1_7 + $utr3                  +
                      $gene->{QI}->[8]*$c1_9 + $so_a + $so_c +
                      $gene->{alt_splicing}*$alt+
                      $gene->{junk}*$junk +$bonus;
	return $value;
}

#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
1;


