#------------------------------------------------------------------------
#----             		maker_gff3.pm	             ---- 
#------------------------------------------------------------------------
package evaluator::maker_gff3;
use strict;

use PostData;
use FileHandle;
use Tie::File;

use Fasta;

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------

sub new {
   my ($class, @args) = @_;

   my $self = {};

   bless ($self, $class);

   my $file = shift @ARGS;
   my @lines;

	tie @array, 'Tie::File', filename or die ...;

        my $fh = new FileHandle;
        die "Can't open $file.\n" unless $fh->open("<$file");
	
	#parse annotation portion of file
        while (defined(my $line=<$fh>)) {
                last if $line=~/^\#\#FASTA/;
		next if $line=~/^\#/;
                chomp $line;
                my @item=split /\t/, $line;
                if ($#item==8) { push @lines,  \@item; }
        }

	#parse fasta portion of file (can be multifasta)
	my $fasta;
	while (my $line = <$fh>) {
		next if $line=~/^>/;
		if ( $line=~ /^([A-Za-z]+)\n/) {
			$fasta .= $1;
		}
	}

	$fh->close;

	my $info = {'file'=>$file, 'gff'=>\@lines, 'fasta'=>\$fasta};

	my $rm_hits 	    = get_phat_hits($info,'repeatmasker');

	my $snap_hits 	    = get_phat_hits($info,'snap');
	my $augustus_hits   = get_phat_hits($info, 'augustus');
	my $fgenesh_hits    = get_phat_hits($info, 'fgenesh');
	
	my $blastn_hits     = get_phat_hits($info,'blastn');
	my $blastx_hits     = get_phat_hits($info,'blastx');
	my $pro2genome_hits = get_phat_hits($info,'protein2genome');
	my $est2genome_hits = get_phat_hits($info,'est2genome');

	##Then return all the hits to evaluator!!

	return { 'rm' 	    => $rm_hits,
		 'snap'     => $snap_hits,
		 'augustus' => $augustus_hits,
		 'fgenesh'  => $fgenesh_hits,
		 'blastn'   => $blastn_hits,
		 'blastx'   => $blastx_hits,
		 'pro2ge'   => $pro2genome_hits,
		 'est2ge'   => $est2genome_hits,
	       };
}


#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub creating_temp_fasta {
	my $gff_file = shift;
	
	my $fh = new FileHandle;
	$fh->open($gff_file);

	my $content;

	while (my $line = <$fh>) {
		last if $line =~ /^##FASTA/;
	}

	while (my $line = <$fh>) {
		$content .= $line;
	}
	$fh->close;

	my ($file_name) = $content =~ /^>(.+?)[\n|\s]/;
	$file_name .= '.fastatemp';

	$fh->open(">$file_name");
	
	print $fh $content;
	$fh->close;

	return $file_name;
	
}
#-------------------------------------------------------------------------------

sub get_phat_hits {
	my $info = shift;
	my $hit_type = shift;


	my $lines = $info->{gff};
	my $pre_hits  = get_selected_hits($lines, $hit_type);

	return undef unless $pre_hits;
	
	my $q_len = length(${$info->{fasta}});
	my $q_name = $info->{gff}->[0]->[0];
	
	my @hits;
	while (my ($pre_hit_id, $pre_hit) = each %$pre_hits) {
		
		my $hit_len = abs($pre_hit->{e}-$pre_hit->{s}) + 1;

		my $hit = new Hits::PhatHit('-name'	    => $pre_hit->{id},
					    '-description'  => $pre_hit->{description},
					    '-algorithm'    => $hit_type,
					    '-length'       => $q_len,
					    '-score'        => get_pre_hit_score($pre_hit),
					   );

		$hit->{hit_id} = $pre_hit->{id};
		$hit->{id}     = $pre_hit->{id};
		$hit->{pre}= $pre_hit;

		$hit->queryLength($q_len);
	
		my $hit_start = 1;
		foreach my $part (@{$pre_hit->{parts}}) {
			my @args;
			my $hsp_seq = get_hsp_seq($part->{s}, $part->{e},
					$part->{str}, $info->{fasta});		
				## The hsp seq is revcomped!
			my $hit_end = abs($part->{e}-$part->{s}) + $hit_start;

			push @args, '-query_start';
			push @args, $part->{s};

			push @args, '-query_seq';
			push @args, $hsp_seq;

			push @args, '-score';
			if ($part->{score} eq '.') {push @args, 10000000;}
			else {push @args, $part->{score};}
		
			push @args, '-homology_seq';
			push @args, $hsp_seq;  ## Later maybe we can add real
						# homology sequence.

			push @args, '-hit_start';
			if ($part->{t_s}) {push @args, $part->{t_s};}
			else {push @args, $hit_start;}

			push @args, '-hit_seq';
			push @args, $hsp_seq;
		
			push @args, '-hsp_length';
			push @args, abs($part->{e}-$part->{s}) +1;
	
			push @args, '-identical';
			push @args, abs($part->{e}-$part->{s}) +1;

			push @args, '-hit_length';
			push @args, abs($part->{e}-$part->{s}) +1;

			push @args, '-query_name';
			push @args, $q_name;

			push @args, '-algorithm';
			if ($hit_type eq 'snap') {push @args, 'SNAP';}
			else {push @args, $hit_type; }
			
			push @args, '-bits';
			push @args, undef;

			push @args, '-evalue';
			push @args, undef;

                        push @args, '-pvalue';
                        push @args, undef;

			push @args, '-query_length';
			push @args, $q_len;
	
			push @args, '-query_end';
			push @args, $part->{e};
	
			push @args, '-conserved';
			push @args, abs($part->{e} - $part->{s}) + 1;
			
			push @args, '-hit_name';
			push @args, $part->{id};

			push @args, '-hit_end';
			if ($part->{t_e}) {push @args, $part->{t_e};}
			else {push @args, $hit_end;}
	
			push @args, '-query_gaps';
			push @args, 0;
	
			push @args, '-hit_gaps';
			push @args, 0;

			my $hsp = new Hits::PhatHsp(@args);
			   $hsp->queryName($q_name);
			   $hsp->{pre} = $part;

                        if ($part->{str} eq '+' ){
                                $hsp->{_strand_hack}->{query} = 1;
                        }
                        else {
                                $hsp->{_strand_hack}->{query} = -1;
                        }
			if ($part->{t_strand}) {
				if ($part->{t_strand} eq '+') {
					$hsp->{_strand_hack}->{hit}   = 1;
				}
				else {
					$hsp->{_strand_hack}->{hit}   = -1;
				}
			}
			else {$hsp->{_strand_hack}->{hit} =1;}
                        $hit->add_hsp($hsp);

			$hit_start += abs($part->{e} - $part->{s}) + 1;
		}
		push @hits, $hit;
	}

	## Is there any other thing needs to be done? (like set minimum score)

	return \@hits;
}
	

#-------------------------------------------------------------------------------
sub get_hsp_seq {
	my $s = shift;
	my $e = shift;
	my $str = shift;
	my $fasta = shift;

	my $length = abs($e-$s) +1;
	my $hsp_seq = substr($$fasta, $s-1, $length);
	$hsp_seq = Fasta::revComp($hsp_seq) if $str eq '-';

	return $hsp_seq;
}

#-------------------------------------------------------------------------------
sub get_selected_hits {
	my $lines = shift;
	my $hit_type = shift;

	my @unordered;
	foreach my $items (@$lines) {
		next unless $items->[1] eq $hit_type;

		my $data ={};
		$data->{seq_id}      = $items->[0];
		$data->{class}       = $items->[1];
		$data->{type}  	     = $items->[2];
		$data->{s}	     = $items->[3];
		$data->{e}	     = $items->[4];
		$data->{score}	     = $items->[5];
		$data->{str}	     = $items->[6];
		$data->{description} = $items->[8];

		($data->{id}, $data->{target}, $data->{parent}, $data->{name})
			= parse_description($data->{description});

		($data->{t_id},$data->{t_s}, $data->{t_e}, $data->{t_strand})
			= parse_target($data->{target}) 
				if $data->{target};

		push @unordered, $data;
	}

	return undef unless $#unordered >0;

	my $hits = {};
	foreach my $one (@unordered) {
		$hits->{$one->{id}} = $one if $one->{type} ne 'match_part';
	}

	foreach my $one (@unordered) {
		push (@{$hits->{$one->{parent}}->{parts}}, $one) 
			if $one->{type} eq 'match_part';
	}

	return $hits;
}		
#------------------------------------------------------------------------
sub parse_target {
	my $target = shift;

	my @items = split / /, $target;
	return @items;

}
#------------------------------------------------------------------------
sub parse_description {
	my $info = shift;

	my @items = split /;/, $info;

	my $things = {};
	foreach my $item (@items) {
		$things->{$1} = $2 if $item=~ /^(.+?)=(.*)$/;
	}
	
	return ($things->{ID},$things->{Target},
		$things->{Parent}, $things->{Name});

}
#------------------------------------------------------------------------
sub get_pre_hit_score {
	my $pre_hit = shift;
	
	return $pre_hit->{score} if $pre_hit->{score} ne '.';

	my $total;
	foreach my $part (@{$pre_hit->{parts}}) {

		return '.' if $part->{score} eq '.';
		$total += $part->{score} 
			if $part->{score} =~ /(-|\d|\.)+/;
	}
	return $total;
}
#------------------------------------------------------------------------
1;


