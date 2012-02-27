#------------------------------------------------------------------------
#----          evalautor::fathom_utils.pm   			             ---- 
#------------------------------------------------------------------------
package evaluator::fathom_utils;
use strict;
use FileHandle;
use evaluator::pseudo_hit;
use Fasta;
use FastaFile;

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------

sub snap_backwards {
	my $box = shift;
	my $CTL = shift;
	my $the_void = shift;

	my $eat = $box->{transcript};
	my $seq = $box->{seq};
	my $name = $box->{t_name};
	$name = Fasta::seqID2SafeID($name);

	my $pseudo_hit = evaluator::pseudo_hit::convert_to_pseudo_hit
				($eat);

	my ($chop_hit, $chop_seq) = chophit($pseudo_hit, $seq, $CTL->{pred_flank});
	my $noUTR_hit = remove_utr($chop_hit, $box->{'5_len'}, $box->{'3_len'});

	my $fasta_file = create_fasta_file($chop_seq, $name, $the_void);
	my $anno_file =  create_anno_file ($noUTR_hit, $name, $the_void);
	
	my $fathom_report = run_fathom($fasta_file, $anno_file, $name, $the_void, $CTL);

	my ($snap_score, $all_scores) = parse_fathom_output($fathom_report, $name);
	
	return {'overall_score'		=>	$snap_score,
		'all_scores'		=>	$all_scores,
		};
}	

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub remove_utr {
	my $hit = shift;
	my $five = shift;
	my $three = shift;

	($five, $three) = ($three, $five) if $hit->{strand} == -1;

	my @hsp_1;
	foreach my $hsp (@{$hit->{hsps}}) {
		my $hsp_length = $hsp->[1]-$hsp->[0] +1;
		if ($five == 0) {
			push @hsp_1, [$hsp->[0], $hsp->[1]];
		}
		elsif ($hsp_length <= $five) {
			$five -= $hsp_length;
		}
		else {
			push @hsp_1, [(($hsp->[0]) + $five), $hsp->[1]];
			$five = 0;
		}
	}

	my @hsp_2;
	for(my $i = $#hsp_1; $i>=0; $i--) {
		my $hsp = $hsp_1[$i];
		my $hsp_length = $hsp->[1] - $hsp->[0] + 1;
		if ($three == 0) {
			unshift @hsp_2, [$hsp->[0], $hsp->[1]];
		}
		elsif ($hsp_length <= $three) {
			$three -= $hsp_length;
		}
		else {
			unshift @hsp_2, [$hsp->[0], (($hsp->[1]) - $three)];
			$three = 0;
		}
	}

	my $new_hit = { 'b' => $hsp_2[0]->[0],
			'e' => $hsp_2[$#hsp_2]->[1],
			'strand' => $hit->{strand},
			'hsps' => \@hsp_2,
			};

	return $new_hit;
}		
#-------------------------------------------------------------------------------
sub parse_fathom_output {
	my $file = shift;
	my $name = shift;

	my $fh = new FileHandle;
	$fh->open($file);
	
	my @items;
	my $status = 1;

	while (my $line = <$fh>) {
		chomp $line;

		my @i = split /\t/, $line;
		next unless scalar @i == 9;
		next unless $i[8] eq $name;

		push @items, [$i[1], $i[4]];
	}

	my @sorted_items = sort {$a->[0] <=> $b->[0]} @items;
	
	my @scores;
	my $score = 0;
	foreach my $item (@sorted_items) {
		push @scores, $item->[1];

		if ($item->[1] =~ /^[-\.\d]+$/ 
			&& $item->[1] ne '.') {
			$score += $item->[1];
		}

		else { 
			$status = -1;
		}

	}

	return ($score, \@scores, $status);
}
		
#-------------------------------------------------------------------------------
sub run_fathom {
	my $fasta_file = shift;
	my $anno_file = shift;
	my $name = shift;
	my $the_void = shift;
	my $CTL = shift;

        my $file_name = $the_void."/".$name;
        $file_name .= '.fathom.output';

	my $command = '';
	$command .= $CTL->{fathom};
	$command .= " $anno_file $fasta_file";
	$command .= " -score-genes ";
	$command .= $CTL->{snaphmm};
	#$command .= " -errors-ok";
	$command .= " >$file_name";

	print STDERR "Running fathom over $name ...\n";

	`$command`;

	return $file_name;
}

#-------------------------------------------------------------------------------
sub create_anno_file {
	my $hit = shift;
	my $name = shift;
	my $the_void = shift;

	my $anno = '';
	$anno .= ">$name\n";

	foreach my $hsp (@{$hit->{hsps}}) {
		my $left = $hsp->[0];
		my $right = $hsp->[1];

		my $label;
		if (scalar @{$hit->{hsps}} == 1) {
			$label = 'Esngl';
		}
		elsif ( $hit->{b} == $left ) {
			if ($hit->{strand} == 1) { $label = 'Einit';}
			else { $label = 'Eterm'; }
		}
		elsif ( $hit->{e} == $right) {
			if ($hit->{strand} == 1) { $label = 'Eterm';}
			else { $label = 'Einit'; }
		}
		else { $label = 'Exon'; }

		$anno .= $label."\t";
		if ($hit->{strand} == 1) {
			$anno .= $left."\t".$right."\t";
		}
		else {
			$anno .= $right."\t".$left."\t";
		}

		$anno .= $name."\n";
	}

	my $file_name = $the_void."/".$name;
        $file_name .= '.fathom.hint';

	my $fh = new FileHandle;
	$fh->open(">$file_name");
	print $fh $anno;
	$fh->close;

	return $file_name;

}
	
#-------------------------------------------------------------------------------
sub chophit {
	my $pseudo_hit = shift;
	my $seq = shift;
	my $flank = shift;

	my ($left , $right);
	$left = ($pseudo_hit->{b} - $flank >= 1)? ($pseudo_hit->{b} - $flank) :
			1;
	$right = ($pseudo_hit->{e} + $flank <= length($$seq))? ($pseudo_hit
			->{e} + $flank ): length($$seq);

	my $chop_seq = substr $$seq, $left-1, ($right - $left + 1);

	my @new_hsps;
	foreach my $hsp (@{$pseudo_hit->{hsps}}) {
		push @new_hsps, [$hsp->[0]-$left+1, $hsp->[1]-$left+1];
	}


	my $chop_hit = { 'b'	=> $pseudo_hit->{b} - $left + 1,
			 'e'	=> $pseudo_hit->{e} - $left + 1,
			 'strand' => $pseudo_hit->{strand},
			 'hsps' => \@new_hsps,
		       };

	return ($chop_hit, \$chop_seq);
	
}

#-------------------------------------------------------------------------------
sub create_fasta_file {
	my $seq = shift;
	my $name = shift;
	my $the_void = shift;

	my $file_name = $the_void."/".$name;
	$file_name .= '.fathom.fasta';
	
	my $fasta = Fasta::toFasta('>'.$name, $seq);
	FastaFile::writeFile(\$fasta, $file_name);
	
	return $file_name;
}
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
1;


