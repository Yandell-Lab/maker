#------------------------------------------------------------------------
#----                          evaluator::quality_seq                ----
#------------------------------------------------------------------------
package evaluator::quality_seq;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use Exporter;
use Shadower;
@ISA = qw(
       );

our %MATRIX;
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub prep {
	my $eat     = shift;
	my $pol_est = shift;
	my $pol_pro = shift;
	my $blastx  = shift;
	my $snaps   = shift;
	my $seq     = shift;

	our %MATRIX = get_matrix();

#	my $eat_t_bs   = get_binary_string([$eat], $seq,   'hsps');
#	my $eat_c_bs   = get_binary_string([$eat], $seq,   'cdss');
#	my $blastx_bs  = get_binary_string($blastx, $seq,  'hsps');
#	my $pol_est_bs = get_binary_string($pol_est, $seq, 'hsps');
#	my $pol_pro_bs = get_binary_string($pol_pro, $seq, 'hsps');
#	my $snaps_bs   = get_binary_string($snaps, $seq,   'hsps');
	
	#---mod
	my $eat_t_bs   = get_binary_array([$eat], $eat,   'hsps');
	my $eat_c_bs   = get_binary_array([$eat], $eat,   'cdss');
	my $blastx_bs  = get_binary_array($blastx, $eat,  'hsps');
	my $pol_est_bs = get_binary_array($pol_est, $eat, 'hsps');
	my $pol_pro_bs = get_binary_array($pol_pro, $eat, 'hsps');
	my $snaps_bs   = get_binary_array($snaps, $eat,   'hsps');
	#---

	my $q_data = 
	get_q_seq($eat_t_bs, $eat_c_bs, $blastx_bs, $pol_est_bs, $pol_pro_bs , $snaps_bs);

	
	my $q_seq =  join("|", @{$q_data});

	return \$q_seq;
}
#------------------------------------------------------------------------
sub get_matrix {

	my %codes;
        $codes{'snap blastx pol_pro pol_est'}{c} = 10;
        $codes{'snap blastx pol_pro pol_est'}{u} = -9;
        $codes{'snap blastx pol_pro pol_est'}{i} = -9;
        $codes{'blastx pol_pro pol_est'}{c} = 8;
        $codes{'blastx pol_pro pol_est'}{u} = -9;
        $codes{'blastx pol_pro pol_est'}{i} = -8;
        $codes{'snap pol_pro pol_est'}{c} = 10;
        $codes{'snap pol_pro pol_est'}{u} = -9;
        $codes{'snap pol_pro pol_est'}{i} = -9;
        $codes{'pol_pro pol_est'}{c} = 8;
        $codes{'pol_pro pol_est'}{u} = -9;
        $codes{'pol_pro pol_est'}{i} = -8;
        $codes{'snap blastx pol_est'}{c} = 8;
        $codes{'snap blastx pol_est'}{u} = -9;
        $codes{'snap blastx pol_est'}{i} = -9;
        $codes{'blastx pol_est'}{c} = 7;
        $codes{'blastx pol_est'}{u} = -2;
        $codes{'blastx pol_est'}{i} = -8;
        $codes{'snap pol_est'}{c} = 6;
        $codes{'snap pol_est'}{u} = -5;
        $codes{'snap pol_est'}{i} = -8;
        $codes{'pol_est'}{c} = 5;
        $codes{'pol_est'}{u} = 10;
        $codes{'pol_est'}{i} = -8;
        $codes{'snap blastx pol_pro'}{c} = 4;
        $codes{'snap blastx pol_pro'}{u} = -9;
        $codes{'snap blastx pol_pro'}{i} = -9;
        $codes{'blastx pol_pro'}{c} = 4;
        $codes{'blastx pol_pro'}{u} = -8;
        $codes{'blastx pol_pro'}{i} = -7;
        $codes{'snap pol_pro'}{c} = 4;
        $codes{'snap pol_pro'}{u} = -8;
        $codes{'snap pol_pro'}{i} = -8;
        $codes{'pol_pro'}{c} = 3;
        $codes{'pol_pro'}{u} = -8;
        $codes{'pol_pro'}{i} = -7;
        $codes{'snap blastx'}{c} = 2;
        $codes{'snap blastx'}{u} = -7;
        $codes{'snap blastx'}{i} = -5;
        $codes{'blastx'}{c} = 1;
        $codes{'blastx'}{u} = -5;
        $codes{'blastx'}{i} = -3;
        $codes{'snap'}{c} = 0;
        $codes{'snap'}{u} = -4;
        $codes{'snap'}{i} = -3;
        $codes{'none'}{c} = -9;
        $codes{'none'}{u} = -9;
        $codes{'none'}{i} = 10;

	return %codes;
}
#------------------------------------------------------------------------
sub get_q_seq {
   my $eat_t_bs   = shift;
   my $eat_c_bs   = shift;
   my $blastx_bs  = shift;
   my $pol_est_bs = shift;
   my $pol_pro_bs = shift;
   my $snaps_bs   = shift;
   
   my @scores;
   my $scale = 1;
   for (my $i = 0; $i < @$eat_t_bs; $i++){
      my $area;
      if ($eat_c_bs->[$i] == 1){
	 # in CDS
	 $area = 'c';
      }
      elsif ($eat_t_bs->[$i] == 1 && $eat_c_bs->[$i] == 0){
	 # in UTR
	 $area = 'u';
      }
      elsif ($eat_t_bs->[$i] == 0){
	 # intergenic or intron
	 $area = 'i';
      }
   
      my $key = '';
      $key .= 'snap '    if $snaps_bs->[$i];
      $key .= 'blastx '  if $blastx_bs->[$i];
      $key .= 'pol_pro ' if $pol_pro_bs->[$i];
      $key .= 'pol_est'  if $pol_est_bs->[$i];
      $key = 'none' unless($key);

      my $score = $MATRIX{$key}{$area};
      push(@scores, $score);
   } 
   
   my $averages = average_scores(\@scores, 3);
   
   my $q_seq = prune_to_transcript_only($averages, $eat_t_bs);
   
   return $q_seq;
}
#------------------------------------------------------------------------
#old 12/15/2008
# sub get_q_seq {
#         my $eat_t_bs   = shift;
#         my $eat_c_bs   = shift;
#         my $blastx_bs  = shift;
#         my $pol_est_bs = shift;
#         my $pol_pro_bs = shift;
#         my $snaps_bs   = shift;
#
# 	my @scores;
# 	my $scale = 1;
# 	for (my $i = 0; $i < length($$eat_t_bs); $i++){
# 		my $eat_t   = substr($$eat_t_bs,   $i, $scale);
# 		my $eat_c   = substr($$eat_c_bs,   $i, $scale);
# 		my $blastx  = substr($$blastx_bs,  $i, $scale);
# 		my $pol_est = substr($$pol_est_bs, $i, $scale);
# 		my $pol_pro = substr($$pol_pro_bs, $i, $scale);
# 		my $snaps   = substr($$snaps_bs,   $i, $scale);
#	
# 		my $column = $eat_t.$eat_c.$snaps.$blastx.$pol_pro.$pol_est;
#
# 		my $score = score_column($column);
# 		push(@scores, $score);
# 	} 
#
# 	my $averages = average_scores(\@scores, 3);
#
#
# 	my $q_seq = prune_to_transcript_only($averages, $eat_t_bs);
#
# 	return $q_seq;
# }
#------------------------------------------------------------------------
sub prune_to_transcript_only {
	my $scores   = shift;
	my $eat_t_bs = shift;

	my $size = @{$scores};

	die "Bad data in quality_seq::prune_to_transcript_only\n"
	if @{$eat_t_bs} != $size;

	my @final;
	for (my $i = 0; $i < @{$eat_t_bs}; $i++){
		push(@final, $scores->[$i]) 
		if ($eat_t_bs->[$i] == 1);
	}
	return \@final;
}
#------------------------------------------------------------------------
sub average_scores {
	my $scores = shift;
	my $scale  = shift;

	my @averages;
	for (my $i = 0; $i < @{$scores}; $i++){
		my $data = shift_multi($scores, $i, $scale);
		my $sum  = 0;

		my $size = 0;
		foreach my $datum (@{$data}){
			next unless defined $datum;
			$sum += $datum;
			$size++;
		}

		my $ave =  $size == 0 ? 0 : substr($sum/$size, 0, 5);
		push(@averages, $ave);
	}

	return \@averages;
}
#------------------------------------------------------------------------
sub shift_multi {
	my $array = shift;
	my $x     = shift;
	my $num   = shift;

	my @data;
	for (my $i = 0; $i < $num; $i++){
		push(@data, $array->[$x + $i]);	
	}
	return \@data;
}
#------------------------------------------------------------------------
sub score_column {
	my $column = shift;

	my $area;
	if     (substr($column, 1, 1) == 1){
		# in CDS
		$area = 'c';
	}
	elsif (substr($column, 0, 1) == 1 && substr($column, 1, 1) == 0){
		# in UTR
		$area = 'u';
	}
	elsif (substr($column, 0, 1) == 0){
		# intergenic or intron
		$area = 'i';
	}

	my $snap    = substr($column, 2, 1);
	my $blastx  = substr($column, 3, 1);
	my $pol_pro = substr($column, 4, 1);
	my $pol_est = substr($column, 5, 1);

	my $key = '';
	   $key .= 'snap '    if $snap;
	   $key .= 'blastx '  if $blastx;
	   $key .= 'pol_pro ' if $pol_pro;
	   $key .= 'pol_est'  if $pol_est;

	  $key = 'none' unless($key);

	my $score = $MATRIX{$key}{$area};
	
	return $score;
}
#------------------------------------------------------------------------
sub get_binary_string {
	my $hits = shift;
	my $seq  = shift;
	my $flag = shift;

	my @coors;
	foreach my $hit (@{$hits}){
		my @hsps;
		if ($flag eq 'hsps'){
			@hsps = $hit->hsps() if defined $hit->hsps();
		}
		elsif ($flag eq 'cdss'){
			@hsps = @{$hit->cdss()} if defined $hit->cdss();
		}
		else {
		}
		foreach my $hsp (@hsps){
			push(@coors, [$hsp->nB('query'), $hsp->nE('query')]);
		}
	}
	
	my $masked_seq = Shadower::maskSequence($seq, \@coors, 0, '1');
	  $$masked_seq =~ s/[^1]/0/g;
	
	return $masked_seq;
}
#------------------------------------------------------------------------
sub get_binary_array {
	my $hits = shift;
	my $tran = shift;
	my $flag = shift;

	my ($start, $end) = ($tran->start('query'), $tran->end('query'));
	my $length = $end - $start + 1;
	my $offset = $start; # do not - 1 so as to give array coors
	my @b_seq = map {0} (1..$length);

	my @coors;
	foreach my $hit (@{$hits}){
		my @hsps;
		if ($flag eq 'hsps'){
			@hsps = $hit->hsps() if defined $hit->hsps();
		}
		elsif ($flag eq 'cdss'){
			@hsps = @{$hit->cdss()} if defined $hit->cdss();
		}
		else {
		}

		foreach my $hsp (@hsps){
		   my $s = $hsp->start('query') - $offset;
		   my $e = $hsp->end('query') - $offset;
		   next if($s >= $length); #array space
		   next if($e < 0); #array space
		   $s = 0 if($s < 0); #array space
		   $e = $length - 1 if($e >= $length); #array space
		   @b_seq[$s..$e] = map {1} ($s..$e);
		}
	}
		
	return \@b_seq;
}	
#------------------------------------------------------------------------
1;


