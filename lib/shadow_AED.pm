#------------------------------------------------------------------------
#----                          shadow_AED.pm                         ---- 
#------------------------------------------------------------------------
package shadow_AED;
use strict;

sub get_abAED{
    my $hits = shift;
    my $tran = shift;

    my $sum;
    foreach my $h (@$hits){
	$sum += get_AED([$h], $tran);
    }
    return 1 if(! @$hits);
    return $sum/@$hits;
}

sub get_eAED {
   my $hits = shift;
   my $tran = shift;

   return 1 if(! @{$hits} || ! $tran);
   
   my @ests = grep {$_->algorithm =~ /est2genome|^est_gff|^altest_gff|^blastn|^tblastx/} @$hits;
   my @prots = grep {$_->algorithm =~ /protein2genome|^protein_gff|^blastx/} @$hits;

   #build splice site index for trans
   my %splices;
   my @hsps = $tran->hsps;
   for( my $i = 0; $i < @hsps; $i++){
       my $hsp = $hsps[$i];
       $splices{start}{$hsp->start}++ unless($i == 0);
       $splices{end}{$hsp->end}++ unless($i == @hsps - 1);
   }

   #filter out hsps that don't have a matching splice site
   my @keepers;
   foreach my $h (@ests){
       my @hsps = $h->hsps;
       if(@hsps <= 1){
	   push(@keepers, @hsps);
	   next;
       }

       foreach my $hsp (@hsps){
	   
       }
   }

}

sub get_AED {
   my $hits = shift;
   my $tran = shift;

   return 1 if(! @{$hits} || ! $tran);

   my ($start, $end) = ($tran->start('query'), $tran->end('query'));

   foreach my $hit (@{$hits}){
      my ($hs, $he) = ($hit->start('query'), $hit->end('query'));

      $start = $hs if($hs < $start);
      $end = $he if($he > $end);
   }

   my $length = $end - $start + 1;
   my $offset = $start; # do not - 1 so as to give array space coors
   my @b_seq = map {0} (1..$length);
   
   #map out hit space
   foreach my $hit (@{$hits}){
      my @hsps = $hit->hsps() if defined $hit->hsps();      
      
      foreach my $hsp (@hsps){
	 my $s = $hsp->start('query') - $offset;
	 my $e = $hsp->end('query') - $offset;

	 #array space coors
	 die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
	 die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

	 @b_seq[$s..$e] = map {1} ($s..$e);
      }
   }

   #map out transcript space
   my @hsps = $tran->hsps() if defined $tran->hsps();      
   
   foreach my $hsp (@hsps){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);
      
      foreach my $i ($s..$e){
	 $b_seq[$i] += 2;
      }
   }
   
   #calculate AED
   my %index = (0 => 0,
		1 => 0,
		2 => 0,
		3 => 0,
	       );

   foreach my $i (@b_seq){
      $index{$i}++;
   }

   my $spec = $index{3}/($index{2} + $index{3}); #specificity
   my $sens = $index{3}/($index{1} + $index{3}); #sensitivity
   my $AED = 1 - ($spec + $sens)/2;

   return $AED;
}

1;
