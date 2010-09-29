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

   #get coordinates to build array in memory
   my ($start, $end) = ($tran->start('query'), $tran->end('query'));
   foreach my $hit (@{$hits}){
       my ($hs, $he) = ($hit->start('query'), $hit->end('query'));

       $start = $hs if($hs < $start);
       $end = $he if($he > $end);
   }

   #build array in memory
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

      @b_seq[$s..$e] = map {2} ($s..$e); #replaces hit
   }
   
   #seperate out hit types
   my @ests = grep {$_->algorithm =~ /est2genome|^est_gff|^altest_gff|^blastn|^tblastx/} @$hits;
   my @prots = grep {$_->algorithm =~ /protein2genome|^protein_gff|^blastx/} @$hits;

   #==do EST filtering
   my @keepers;

   #build splice site index for transcript
   my %splices;
   for( my $i = 0; $i < @hsps; $i++){
       my $hsp = $hsps[$i];
       $splices{start}{$hsp->start}++ unless($i == 0);
       $splices{end}{$hsp->end}++ unless($i == @hsps - 1);
   }      

   #filter out hsps that don't have a matching splice site
   foreach my $e (@ests){
       my @ehsps = $e->hsps;

       if(@ehsps == 1){ #how to handel single exon ESTs
	   my $ehsp = $ehsps[0];

	   #are ESTs fully contained or hanging off ends of gene
	   my $aB = $ehsp->start;
	   my $aE = $ehsp->end;
	   foreach(my $i = 0; $i < @hsps, $i++){
	       my $hsp = $hsps[$i];
	       my $bB = $hsp->start;
	       my $bE = $hsp->end;

	       my $class = compare::compare($aB, $aE, $bB, $bE);
	       if($class =~ /[IbB1]/ ){ #single HSP EST fully contained in exon of gene
		   push(@keepers, $ehsp);
	       }
	       elsif(@hsps == 1 && $class =~ /[i]/){ #single exon gene fully contained in single HSP EST
		   push(@keepers, $ehsp);
	       }
	       elsif(@hsps == 1 && $class =~ /[aAzZ]/){ #single exon gene overlapped by single HSP EST
		   push(@keepers, $ehsp);
	       }
	       #elsif($i == 0 && $class =~ /[AZ]/){ #single HPS EST hanging off first exon of multi-exon transcript
	       #    push(@keepers, $ehsp);
	       #}
	       #elsif($i == @hsps - 1 && $class =~ /[az]/){ #single HPS EST hanging off last exon of multi-exon transcript
	       #    push(@keepers, $ehsp);
	       #}
	   }
       }
       elsif(@hsps == 1){ #how to handle multi-exon ESTs and single exon gene model
	   my $hsp = $hsps[0];

	   #is the model fully contained in EST HSP or hanging off ends of EST
	   my $bB = $ehsp->start;
	   my $bE = $ehsp->end;
	   foreach(my $i = 0; $i < @ehsps, $i++){
	       my $ehsp = $ehsps[$i];
	       my $aB = $hsp->start;
	       my $aE = $hsp->end;

	       my $class = compare::compare($aB, $aE, $bB, $bE);
	       if($class =~ /[iaA1]/ ){ #single exon gene fully contained in EST HSP
		   push(@keepers, $ehsp);
	       }
	       #elsif($i == 0 && $class =~ /[Bz]/){ #single exon gene hanging off first HSP of EST
	       #    push(@keepers, $ehsp);
	       #}
	       #elsif($i == @hsps - 1 && $class =~ /[bZ]/){ #single exon gene hanging off last HSP of EST
	       #    push(@keepers, $ehsp);
	       #}
	   }
       }
       else{ #handle multi-exon ESTs and multi-exon gene model
	   for(my $i = 0; $i < @ehsps; $i++){
	       if($i != 0 && splices{start}{$ehsp->start}){ #a first splice site in EST also first splice site in gene
		   push(@keepers, $ehsp);
	       }
	       elsif($i != @hsps - 1 && splices{end}{$ehsp->end}){ #a last splice site in EST also last splice site in gene
		   push(@keepers, $ehsp);
	       }
	   }
       }
   }


   #map keeper ESTs
   foreach my $hsp (@keepers){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

      #only add to transcript overlap
      foreach my $i ($s..$e){
	  $b_seq[$i] += 2 if($b_seq[$i] == 1);
      }
   }

   #==filter proteins hits



   #==calculate AED
   my %index = (0 => 0,
		1 => 0,
		2 => 0,
		3 => 0,
	       );

   foreach my $i (@b_seq){
      $index{$i}++;
   }

   #catch error caused by bad GFF3 input
   die "ERROR: The feature being compared appears to be missing\n".
       "some of it's structure.  This can happen when you use\n".
       "a malformed GFF3 file as input to one of MAKER's evidence\n".
       "passthrough options. Failed on ". $tran->name." (from shadow_AED)\n"
       if($index{2} + $index{3} == 0 || $index{1} + $index{3} == 0);

   my $spec = $index{3}/($index{2} + $index{3}); #specificity
   my $sens = $index{3}/($index{1} + $index{3}); #sensitivity
   my $AED = 1 - ($spec + $sens)/2;

   return $AED;
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

   #catch error caused by bad GFF3 input
   die "ERROR: The feature being compared appears to be missing\n".
       "some of it's structure.  This can happen when you use\n".
       "a malformed GFF3 file as input to one of MAKER's evidence\n".
       "passthrough options. Failed on ". $tran->name." (from shadow_AED)\n"
       if($index{2} + $index{3} == 0 || $index{1} + $index{3} == 0);

   my $spec = $index{3}/($index{2} + $index{3}); #specificity
   my $sens = $index{3}/($index{1} + $index{3}); #sensitivity
   my $AED = 1 - ($spec + $sens)/2;

   return $AED;
}

1;
