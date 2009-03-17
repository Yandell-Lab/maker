#------------------------------------------------------------------------
#----                          maker::join                           ----
#------------------------------------------------------------------------
package maker::join;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use FileHandle;
use PostData;
use Exporter;
use PhatHit_utils;
use compare;
use cluster;
use clean;
@ISA = qw(
       );
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub get_exon_lengths {
	my $f    = shift;
	my $n    = shift;
	my $type = shift;

	my $length = 0;

	my $i = 0;
	my $sorted = PhatHit_utils::sort_hits($f, 'query');
	foreach my $hsp (@{$sorted}){	
		if ($type == 3 && $i < $n) { $i++; next};

		$length += 
		abs($hsp->nB('query') - $hsp->nE('query')) + 1;	

		last if ($type == 5 && $i ==$n);

		$i++;
	}
	return $length;
}
#------------------------------------------------------------------------
sub bing {
	my $thing = shift;
	my $what  = shift;

	if ($what eq 'pre' or $what eq 'pos'){
		if (defined($thing->{$what})){
			return length($thing->{$what});
		}
		else {
			return -1;
		}
	}
	else {
		return $thing->{$what};	
	}
}
#------------------------------------------------------------------------
sub best_five_prime_extension {
           bing($b, 'pre') <=> bing($a, 'pre') 
                           ||
        bing($b, 'five_l') <=> bing($a, 'five_l');
}

#------------------------------------------------------------------------
sub best_three_prime_extension {
	    bing($b, 'pos') <=> bing($a, 'pos')
	    		    ||
	bing($b, 'three_l') <=> bing($a, 'three_l');

}
#------------------------------------------------------------------------
sub find_best_one {
        my $g    = shift;
        my $ests = shift;

        my $strs = get_strings($ests, $g);

        my @candidates;
        foreach my $datum (@{$strs}){
                my $g_to_e_str = $datum->[0];
                my $e_to_g_str = $datum->[1];
                my $est        = $datum->[2];

                #print " PRIOR find_best_one e_to_g_str:$e_to_g_str g_to_e_str:$g_to_e_str\n";
                #sleep 3;

		next unless $g_to_e_str =~ /^[^0]$/;

		my ($pre) = $e_to_g_str =~ /(0*[^0]{1}).*/; 
		my ($pos) = $e_to_g_str =~ /.*([^0]{1}.*)/;

                #print " POST find_best_one pre:$pre pos:$pos\n";
                #sleep 3;

                push(@candidates, load_candidate($e_to_g_str, $est, $pre, $pos));
        }
        my @sorted_five = sort best_five_prime_extension @candidates;
        my $best_five = shift(@sorted_five);
        my @sorted_three = sort best_three_prime_extension @candidates;
        my $best_three = shift(@sorted_three);

        return ($best_five, $best_three);
}
#------------------------------------------------------------------------
sub find_best_five {
        my $g    = shift;
        my $ests = shift;

        my $strs = get_strings($ests, $g);

        my @candidates;
        foreach my $datum (@{$strs}){
                my $g_to_e_str = $datum->[0];
                my $e_to_g_str = $datum->[1];
                my $est        = $datum->[2];

		#print STDERR " PRIOR find_best_five e_to_g_str:$e_to_g_str g_to_e_str:$g_to_e_str\n";
                #sleep 3;

		# always bad, but rules out single exon genes too!
                #next if $g_to_e_str =~ /I/;

		# no extension to be had
                next if $e_to_g_str =~ /^1+.*$/;

                # can't be a 5 prime ext unless first exon of g is B, z or 1.
                next unless $g_to_e_str =~ /^[Bz1].*$/;
        
                # cant be same splice form if has an internal 0.
		#next if     $e_to_g_str =~ /[AZ]1*[^1]+[^0]/;

		# cant be same splice form if has an internal 0.
                #next if     $g_to_e_str =~ /[Bz]1*[^1]+[^0]/;

                my ($pre, $pos) = $e_to_g_str =~ /(0*[ZbA1])(.*)$/;
                
                #print STDERR " POST find_best_five e_to_g_str:$e_to_g_str pre:$pre pos:$pos\n";
                #sleep 3;

		push(@candidates, load_candidate($e_to_g_str, $est, $pre, $pos));
        }
        my @sorted_five = sort best_five_prime_extension @candidates;
        my $best_five = shift(@sorted_five);
        return $best_five;
}
#------------------------------------------------------------------------
sub find_best_three {
	my $g    = shift;
	my $ests = shift;

	my $strs = get_strings($ests, $g);

	my @candidates;
	foreach my $datum (@{$strs}){
		my $g_to_e_str = $datum->[0];
		my $e_to_g_str = $datum->[1];
		my $est        = $datum->[2];

		#print STDERR " PRIOR find_best_three e_to_g_str:$e_to_g_str g_to_e_str:$g_to_e_str\n";
                #sleep 3;

		# no extension to be had
		next if $e_to_g_str =~ /^.*1+$/;

		#print "AAAAAA\n";
		# always bad, but rules out single exon genes too!
                #next if $g_to_e_str =~ /I/;

		#print "BBBBBB\n";

		# can't be a 3 prime ext unless last exon of g is b or 1.
		next unless $g_to_e_str =~ /^.*[b1]$/;
	
		#print "CCCCCC\n";

		# cant be same splice form if has an internal not 1.
		#next if  $g_to_e_str =~ /^.*[^0][^1]+.*[b1]0*$/;	

		#print "DDDDD\n";

		my ($pre, $pos) = $e_to_g_str =~ /(.*)([a1]0*)$/;
		
		#print STDERR " POST find_best_three e_to_g_str:$e_to_g_str pre:$pre pos:$pos\n";
		#sleep 3;

	       push(@candidates, load_candidate($e_to_g_str, $est, $pre, $pos));

	}
        my @sorted_three = sort best_three_prime_extension @candidates;
        my $best_three = shift(@sorted_three);

	#PostData($best_three);
	#die;
	return $best_three;
}
#------------------------------------------------------------------------
sub load_candidate {
	my $e_to_g_str = shift;
	my $est        = shift;
	my $pre        = shift;
	my $pos        = shift;

        my $l_pre = defined($pre) ? length($pre) : -1;
        my $l_pos = defined($pos) ? length($pos) : -1;

        my $l_str = length($e_to_g_str);

        my $five_join  = defined($pre) ? $l_pre - 1 : -1;
      
        my $three_join = defined($pos) ? $l_str - $l_pos: -1;
 
        my $pre_exon_length = 
        defined($pre) ? get_exon_lengths($est, $five_join,  5) : - 1;

        my $pos_exon_length = 
        defined($pos) ? get_exon_lengths($est, $three_join, 3) : -1;

        my $f_j_class = substr($e_to_g_str, $five_join, 1);
        my $t_j_class = substr($e_to_g_str, $three_join, 1);

        my $num_hsps = $est->hsps();

        my $candidate = {'three_j'   => $three_join, 
                         'five_j'    => $five_join,
                         'three_l'   => $pos_exon_length,
                         'five_l'    => $pre_exon_length,
                         'str'       => $e_to_g_str,
                         'pre'       => $pre,
                         'pos'       => $pos,
                         'f_j_class' => $f_j_class,
                         't_j_class' => $t_j_class, 
                         'num_hsps'  => $num_hsps, 
                         'f'         => $est,
                        };

	return $candidate;
}
#------------------------------------------------------------------------
sub get_strings {
	my $ests  = shift;
	my $g     = shift;

	my @strs;
	foreach my $e (@{$ests}){
		next unless $e->strand('query') == $g->strand('query');

		my $e_to_g_str = compare::compare_phat_hits($e, $g, 'query', 3);

		my $g_to_e_str = compare::compare_phat_hits($g, $e, 'query', 3);

		push(@strs, [$g_to_e_str, $e_to_g_str, $e]);
	}
	return \@strs;
}
#------------------------------------------------------------------------
sub join_f {
	my $b_5         = shift;
	my $g           = shift;
	my $b_3         = shift;
	my $q_seq       = shift;

	#$b_5->{f} = $b_5->{f}->name;
	#$b_3->{f} = $b_3->{f}->name;

	#PostData($b_5);
	#print "LLLL\n";
	#PostData($b_3);
	#print "GGGGGGG\n";
	#$g->show();
	#print "PPPPPPPPPPPPP\n";

	#return $g unless (defined($b_5) || defined($b_3));

	my $sorted_g = PhatHit_utils::sort_hits($g, 'query');

	my $sorted_5 = defined($b_5) ? PhatHit_utils::sort_hits($b_5->{f}, 'query')
	                             :[];

	my $sorted_3 = defined($b_3) ? PhatHit_utils::sort_hits($b_3->{f}, 'query')
	                             :[];

	my $join_offset_5 = $b_5->{five_j};
	my $join_offset_3 = $b_3->{five_j};

	my @anno_hsps;

	#-- push on b_5's upsteam hsps
	my $i = 0;
	foreach my $hsp (@{$sorted_5}){
		last if $i ==  $b_5->{five_j};
		push(@anno_hsps, clone_hsp($hsp));	
		$i++;
	}
	#-- merge the 5-prime overlaping features
	if ( !defined($join_offset_5) || $join_offset_5 == -1){
		push(@anno_hsps, clone_hsp($sorted_g->[0]));
	}
	else {
		my $merged_hsp = merge_hsp($b_5,
		                           $sorted_g->[0],
			                   $q_seq,
		                           5,
		                          );

		push(@anno_hsps, $merged_hsp);
	}

	#-- push on the gene pred's interior hsps
	for (my $i = 1; $i < @{$sorted_g} -1; $i++){
		push(@anno_hsps, clone_hsp($sorted_g->[$i]));
	}

        #-- merge the 3-prime overlaping features
	my $omega = $#{$sorted_g};
        if (!defined($join_offset_3) || $join_offset_3 == -1){
                push(@anno_hsps, clone_hsp($sorted_g->[$omega]));
        }
        else {
		# sepcial case for single exon genes
		my $omega_exon;
		if ($g->num_hsps == 1 && defined($anno_hsps[0])){
			$omega_exon = pop(@anno_hsps);
		}
		else {
			$omega_exon = $sorted_g->[$omega];	
		}
		my $merged_hsp = merge_hsp($b_3,
                                           $omega_exon,
                                           $q_seq,
                                           3,
                                           );

                push(@anno_hsps, $merged_hsp);
        }

        #-- push on b_3's downsteam hsps added 52006
	if (defined($b_3->{three_j})){
		for (my $i = $b_3->{three_j} + 1; $i < @{$sorted_3}; $i++){
			push(@anno_hsps, clone_hsp($sorted_3->[$i]));
		} 
	}

        #my $n = @anno_hsps;
        #print "LLLLLLLLLLLLLL:$n\n";
        #die;

	my $new_total_score = 0;
	my $length = 0;
	foreach my $hsp (@anno_hsps){
	    my $score = $hsp->score();
	    $score = 0 if( ! $score || $score eq '.' || $score eq 'NA');
	    $new_total_score += $score;
	    $length += $hsp->length();
        }

	my $hit_class = ref($g);

	my $new_f = new $hit_class('-name'         => $g->name,
                                   '-description'  => $g->description,
                                   '-algorithm'    => $g->algorithm,
                                   '-length'       => $length,
			           '-score'        => $new_total_score, 
                                 );

	my @evidence;
	push(@evidence, $g->name);
	push(@evidence, $b_5->{f}->name) if defined($b_5->{f});
	push(@evidence, $b_3->{f}->name) if defined($b_3->{f});

	$new_f->evidence(\@evidence);
	$new_f->queryLength(length($$q_seq));

	foreach my $hsp (@anno_hsps){
		$new_f->add_hsp($hsp);
	}

	return $new_f;
}
#------------------------------------------------------------------------
sub clone_hsp{
    my $hsp = shift;

    my @args;

    push(@args, '-query_start');
    push(@args, $hsp->start('query'));

    push(@args, '-query_seq');
    push(@args, $hsp->query_string);

    push(@args, '-score');
    push(@args, $hsp->score);

    push(@args, '-homology_seq');
    push(@args, $hsp->homology_string);

    push(@args, '-hit_start');
    push(@args, $hsp->start('hit'));

    push(@args, '-hit_seq');
    push(@args, $hsp->hit_string);

    push(@args, '-hsp_length');
    push(@args, $hsp->length('total'));

    push(@args, '-identical');
    push(@args, $hsp->{IDENTICAL});

    push(@args, '-hit_length');
    push(@args, $hsp->length('hit'));

    push(@args, '-query_name');
    push(@args, $hsp->{QUERY_NAME});

    push(@args, '-algorithm');
    push(@args, $hsp->algorithm);

    push(@args, '-bits');
    push(@args, $hsp->bits);

    push(@args, '-evalue');
    push(@args, $hsp->evalue);

    push(@args, '-pvalue');
    push(@args, $hsp->pvalue);

    push(@args, '-query_length');
    push(@args, $hsp->length('query'));

    push(@args, '-query_end');
    push(@args, $hsp->end('query'));

    push(@args, '-conserved');
    push(@args, $hsp->{CONSERVED});

    push(@args, '-hit_name');
    push(@args, $hsp->name);

    push(@args, '-hit_end');
    push(@args, $hsp->end('hit'));

    push(@args, '-query_gaps');
    push(@args, $hsp->{QUERY_GAPS});

    push(@args, '-hit_gaps');
    push(@args, $hsp->{HIT_GAPS});

    my $REF = ref($hsp);
    my $clone = new $REF(@args);

    $clone->{queryName} = $hsp->{queryName};
    #-------------------------------------------------
    # setting strand because bioperl is all messed up!
    #------------------------------------------------

    $clone->{_strand_hack}->{query} = $hsp->{_strand_hack}->{query};
    $clone->{_strand_hack}->{hit}   = $hsp->{_strand_hack}->{query};

    return $clone;
}
#------------------------------------------------------------------------
sub merge_hsp {
	my $ext   = shift;
	my $g_hsp = shift;
	my $q_seq = shift;
	my $type  = shift;

	my $offset;
	my $class;
	if ($type == 5){
		$offset = $ext->{five_j};	
		$class  = $ext->{f_j_class};
	}
	else {
		$offset = $ext->{three_j};
		$class  = $ext->{t_j_class};
	}

	my $sorted_ext = PhatHit_utils::sort_hits($ext->{f}, 'query');
	my $ext_hsp = $sorted_ext->[$offset];

	$g_hsp = clone_hsp($g_hsp);

	if     ($class eq '1'){
		return $g_hsp;
	}
	elsif ($type == 5 && $class eq 'A'){
		if ($g_hsp->strand('query') == 1 && $ext_hsp->strand('query') == 1){
			$g_hsp->query->location->start($ext_hsp->start);
		}
		elsif($g_hsp->strand('query') == - 1 && $ext_hsp->strand('query') == -1) {
			$g_hsp->query->location->end($ext_hsp->end);
		}
		else {
			die "dead in auto_annotator::merge_hsp\n";
		}
	}
        elsif ($type == 3 && $class eq 'a'){
		if ($g_hsp->strand('query') == 1 && $ext_hsp->strand('query') == 1){
                        $g_hsp->query->location->end($ext_hsp->end);
                }
                elsif ($g_hsp->strand('query') == - 1 && $ext_hsp->strand('query') == -1) {
                        $g_hsp->query->location->start($ext_hsp->start);
                }
		else {
			die "dead in auto_annotator::merge_hsp\n";
		}
        }

	my $exon = {'b'      => $g_hsp->nB('query'),
	            'e'      => $g_hsp->nE('query') ,
		    'strand' => $g_hsp->strand('query'),
		   };

	my $new_exon_seq = Widget::snap::get_exon_seq($exon, $q_seq);

	$g_hsp->{_query_string}    = $new_exon_seq;
	$g_hsp->{_hit_string}      = $new_exon_seq;
	$g_hsp->{_homology_string} = $new_exon_seq;

	$g_hsp->{'_sequenceschanged'} = 1;
	return $g_hsp;
}
#------------------------------------------------------------------------
1;


