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

   #==calculate bp in evidence
   my %index = (0 => 0,
		1 => 0,
		2 => 0,
		3 => 0,
	       );

   foreach my $i (@b_seq){
       $index{$i}++ if($i == 1);
   }

   #map out transcript space
   my @hsps = sort {$a->start <=> $b->start} $tran->hsps() if defined $tran->hsps();
   my $buf = 0;

   foreach my $hsp (@hsps){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

      @b_seq[$s..$e] = map {2} ($s..$e); #replaces hit
   }

   #==calculate bp in hit
   foreach my $i (@b_seq){
       $index{$i}++ if($i == 2);
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
       my @ehsps = sort {$a->start <=> $b->start} $e->hsps;
       
       EHSP: for (my $i = 0; $i < @ehsps; $i++){
	   my $ehsp = $ehsps[$i];
	   my $aB = $ehsp->start;
	   my $aE = $ehsp->end;
	   my $aL = abs($aE - $aB) + 1;
	   
	   #splice site matches so keep EST HSP
	   if($i != 0 && $splices{start}{$aB}){ #a first splice site in EST also first splice site in gene
	       push(@keepers, $ehsp);
	       next EHSP;
	   }
	   elsif($i != @hsps - 1 && $splices{end}{$aE}){ #a last splice site in EST also last splice site in gene
	       push(@keepers, $ehsp);
	       next EHSP;
	   }
	   
	   #no splice site match so check overlap of EST HSP with exon
	   for(my $j = 0; $j < @hsps; $j++){
	       my $hsp = $hsps[$j]; 
	       my $bB = $hsp->start;
	       my $bE = $hsp->end;
	       my $bL = abs($bE - $bB) + 1;
	       
	       my $class = compare::compare($aB, $aE, $bB, $bE);
	       next if($class eq '0');
	       
	       my $oB = ($aB > $bB) ? $aB : $bB; #overlap begin
	       my $oE = ($aE < $bE) ? $aE : $bE; #overlap end
	       my $oL = abs($oE - $oB) + 1;
	       
	       if(@hsps == 1 && @ehsps == 1 && $oL/$bL >= .9){ #at least 90% of single exon gene overlapped by single HSP EST
		   push(@keepers, $ehsp);
		   next EHSP;
	       }
	       
	       if($oL/$bL >= .9 && $oL/$aL >= .9){ #at least 90% of exon and 90% of single HSP EST overlap
		   push(@keepers, $ehsp);
		   next EHSP;
	       }
	   }
       }
   }

   #map keeper EST HSPs
   foreach my $hsp (@keepers){
      my $s = $hsp->start('query') - $offset;
      my $e = $hsp->end('query') - $offset;
      
      #array space coors
      die "ERROR: Start value not permited!!\n" if($s >= $length || $s < 0);
      die "ERROR: End value not permited!!\n" if($e < 0 || $e >= $length);

      #only add to transcript overlap
      foreach my $i ($s..$e){
	  $b_seq[$i] += 1 if($b_seq[$i] == 2);
      }
   }

   #==map protein hits by phase
   my @ok_frames;
   foreach my $p (@prots){
       foreach my $phsp ($p->hsps){
	   my $start = $phsp->start('query') - $offset;
	   my $end = $phsp->end('query') - $offset;

	   my $pos = $start; #array position
	   my $cigar = $phsp->cigar_string();
	   if(! $cigar){ #if no gap attribute than we assume translation begins at first bp
	       my $length = abs($end - $start) + 1;
	       $cigar .= 'M'.int($length/3);
	   }

	   my @gap = $cigar =~ /([A-Z]\d+)/g;
	   foreach my $g (@gap){
	       $g =~ /([A-Z])(\d+)/;
	       if($1 eq 'F'){
		   $pos += $2;
	       }
	       elsif($1 eq 'R'){
		   $pos -= $2;
	       }
	       elsif($1 eq 'D'){
		   $pos += ($2 * 3);
	       }
	       elsif($1 eq 'M'){
		   my $go = $2;
		   while($go--){
		       if($p->strand('query') == 1){
			   $ok_frames[$pos+0]->{0}++;
			   $ok_frames[$pos+1]->{1}++;
			   $ok_frames[$pos+2]->{2}++;
		       }
		       else{
			   $ok_frames[$pos+2]->{0}++;
			   $ok_frames[$pos+1]->{1}++;
			   $ok_frames[$pos+0]->{2}++;
		       }
		       $pos += 3;
		   }
	       }
	   }
       }
   }

   #compare phase to transcript
   my $aB = $tran->{_TSTART}{query} - $offset;
   my $aE = $tran->{_TEND}{query} - $offset;
   my $phase = 0;
   ($aB, $aE) = ($aE, $aB) if($aB > $aE);
   foreach my $hsp (@{PhatHit_utils::sort_hits($tran, 'query')}){
       my $bB = $hsp->start('query') - $offset;
       my $bE = $hsp->end('query') - $offset;
       my $class = compare::compare($aB, $aE, $bB, $bE);
       next if ($class eq '0');

       $bB = $aB if($aB > $bB);
       $bE = $aE if($aE < $bE);
       

       my @select =  ($tran->strand('query') == 1) ? ($bB..$bE) : reverse($bB..$bE);
       foreach my $i (@select){
	   $b_seq[$i] += 1 if($b_seq[$i] == 2 && $ok_frames[$i]->{$phase});
	   $phase = ($phase + 1) % 3
       }
   }

   #==calculate overlap
   foreach my $i (@b_seq){
       $index{$i}++ if($i == 3);
   }

   #catch error caused by bad GFF3 input (i.e. hits with no HSPs)
   die "ERROR: The feature being compared appears to be missing\n".
       "some of it's structure.  This can happen when you use\n".
       "a malformed GFF3 file as input to one of MAKER's evidence\n".
       "passthrough options. Failed on ". $tran->name." (from shadow_AED)\n"
       if($index{2} + $index{3} == 0 || $index{1} + $index{3} == 0);

   my $spec = $index{3}/($index{2}); #specificity
   my $sens = $index{3}/($index{1}); #sensitivity
   my $eAED = 1 - ($spec + $sens)/2;

   return $eAED;
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
