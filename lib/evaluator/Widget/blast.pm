#------------------------------------------------------------------------
#---- 	            	blast.pm			             ---- 
#------------------------------------------------------------------------
package evaluator::Widget::blast;
use strict;
use FileHandle;
use lib '/home/hao/me/lib/';
use Fasta;
use FastaFile;
use Widget::blastn;
use PhatHit_utils;
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub short_est_support {
	my %arg = @_;

	my $splice_sites= $arg{splice_sites};
	my $side_thre   = $arg{side_thre};


	my $hits_locations = get_blastn_hits(%arg);

	my %support_times ;
	foreach my $site (@$splice_sites) {
		$support_times{$site} = 0;
		foreach my $est (@$hits_locations) {
			my ($nB, $nE) = @$est;
			next unless ($site>=($nB+$side_thre) &&
					($site+$side_thre-1)<=$nE);
			$support_times{$site} ++;
		}
	}
	
	return \%support_times;
			
}

#------------------------------------------------------------------------------
sub get_blastn_hits {
	my %arg = @_;
	
        my $path        = $arg{path};
        my $cpus        = $arg{cpus};
        my $seq         = $arg{seq};
        my $db          = $arg{db};   ## Database should  have been built!
        my $t_name      = $arg{t_name};
        my $window_size = $arg{window_size};
        my $expection   = $arg{expection};
        my $blastn      = $arg{blastn};
	my $hspmax	= $arg{hspmax};
	my $gspmax	= $arg{gspmax};
	my $split_hit   = $arg{split_hit};
	my $percov	= $arg{percov};
	my $percid	= $arg{percid};
	my $bit_min	= $arg{hsp_bit_min};

	my %params = (  'significance' => $expection,
			'hsp_bit_min'  => $bit_min,
			'percov'       => $percov,
			'percid'       => $percid,
			'split_hit'    => $split_hit,
			);

	

	my $tmp_seq_file = $path.'/'.$t_name."_transcript_seq.fasta";
	
	my $fasta = Fasta::toFasta(">$t_name", $seq);
	FastaFile::writeFile($fasta, $tmp_seq_file);

	my $out_name = $path."/".$t_name."_blastn.shortreads";
	
	runBlastn($tmp_seq_file, $db, $out_name, $blastn, 
		 $expection, $cpus, $hspmax, $gspmax, $split_hit,
		 $percov, $percid, $bit_min);

	my $keepers = Widget::blastn::parse($out_name, \%params);
	PhatHit_utils::add_offset($keepers, 0);
	
	my @shattered_hits;
	foreach my $hit (@$keepers) {
		my $hits  = PhatHit_utils::shatter_hit($hit);
		push @shattered_hits, @$hits;
	}


	my @coors;
	foreach my $hit (@shattered_hits) {
		my ($nB, $nE) = PhatHit_utils::get_span_of_hit($hit, 'query');
		($nB, $nE) = ($nE, $nB) if ($nB>$nE);
		push @coors, [$nB, $nE];
	}

	return \@coors;	

}

#------------------------------------------------------------------------------
sub runBlastn {
	my $tmp_seq_file = shift;
	my $db		 = shift;
	my $out_name	 = shift;
	my $blastn	 = shift;
	my $expection	 = shift;
	my $cpus	 = shift;
	my $hspmax	 = shift;
	my $gspmax 	 = shift;
	my $split_hit	 = shift;
	my $percov	 = shift;
	my $percid	 = shift;
	my $bit_min	 = shift;

	return if -e $out_name;
	
	my $command = $blastn;
	$command .= " $db $tmp_seq_file B=100000 V=100000 E=$expection";
	$command .= " R=3";
	$command .= " W=15";
	$command .= " M=1";
	$command .= " N=-3";
	$command .= " Q=3";
	$command .= " Z=1000";
	$command .= " Y=500000000";
	$command .= " cpus=$cpus";
	$command .= " topcomboN=1";
	$command .= " hspmax=$hspmax";
	$command .= " gspmax=$gspmax";
	$command .= " hspsepqmax=$split_hit";
	$command .= " hspsepsmax=$split_hit";
	$command .= " lcmask";
	$command .= " wordmask=seg";
	$command .= " gi";
	$command .= " > $out_name";

	`$command`;
}

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
1;


