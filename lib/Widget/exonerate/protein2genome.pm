#------------------------------------------------------------------------
#----                   Widget::exonerate::protein2genome            ---- 
#------------------------------------------------------------------------
package Widget::exonerate::protein2genome;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget::exonerate;
use Bio::Search::Hit::PhatHit::protein2genome;
use Bio::Search::HSP::PhatHSP::protein2genome;
@ISA = qw(
	Widget::exonerate
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class  = shift;
        my @args   = @_;

        my $self = $class->SUPER::new(@args);

	bless ($self, $class);
        return $self;
}
#------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_blocks {
	my $file = shift;

      my $fh = new FileHandle();
         $fh->open($file);

	local $/ = 'C4 Alignment:';

	my @chunks;
	while (my $line = <$fh>){
		push(@chunks, $line);
	} 
	$fh->close();

	local $/ = "\n";

	#checks if exonerate realy finished
	unless(grep {/completed exonerate analysis/} $chunks[-1]){
            unlink($file);
            die"ERROR: The file $file appears to be incomplete\n".
                "MAKER will need to delete the file, before trying again\n\n";
        }

	return \@chunks;
}
#-------------------------------------------------------------------------------
sub get_exon_coors {
	my $v    = shift;
	my $type = shift;

	my $pos_q = $v->{q_b};
	my $pos_t = $v->{t_b};

	my $exon = 0;

	my @data;
        foreach my $o (@{$v->{operations}}){
		if    ($o->{state} eq 'M'){

		        #print "ZBBBBBBBBBBBB:$pos_q exon:$exon\n";
		        if ($v->{q_b} < $v->{q_e}){
			    $data[$exon]{q}{b} = $pos_q + 1
				unless defined($data[$exon]{q}{b});
		        }
		        else {
			    $data[$exon]{q}{b} = $pos_q
				unless defined($data[$exon]{q}{b});
		        }

                        if ($v->{t_b} < $v->{t_e}){
			    $data[$exon]{t}{b} = $pos_t + 1
				unless defined($data[$exon]{t}{b});
                        }
                        else {
			    $data[$exon]{t}{b} = $pos_t
				unless defined($data[$exon]{t}{b});
                        }

                     	if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

			if ($v->{t_strand} == 1){
				$pos_t += $o->{t};
			}
			else {
				$pos_t -= $o->{t};
			}

			#print "ZCCCCCCCCCC:$pos_q exon:$exon\n";
		}
		elsif ($o->{state} eq 'G'){

                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }
		}
                elsif ($o->{state} eq 'N'){
			die "dead in est2genomic::get_hsp_coors:N\n";
                }
                elsif ($o->{state} eq '5'){
			if ($type eq '5I3'){
			        if ($v->{q_b} < $v->{q_e}){
				    $data[$exon]{q}{e} = $pos_q;
				}
				else{
				    $data[$exon]{q}{e} = $pos_q + 1;
				}
				
				if ($v->{t_b} < $v->{t_e}){
				    $data[$exon]{t}{e} = $pos_t;
				}
				else{
				    $data[$exon]{t}{e} = $pos_t + 1;
				}
				$data[$exon]{q}{strand} = $v->{q_strand};
				$data[$exon]{t}{strand} = $v->{t_strand};

				#print "ZDDDDDDDDD:$pos_q exon:$exon\n";

                                #fix 0 length exons in report
                                if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
                                    $data[$exon] = undef;
                                    $exon--; #undo iteration
                                }

				$exon++;
			}
                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }

                }
                elsif ($o->{state} eq '3'){
			if ($type eq '3I5'){
			        if ($v->{q_b} < $v->{q_e}){
				    $data[$exon]{q}{e} = $pos_q;
				}
				else{
				    $data[$exon]{q}{e} = $pos_q + 1;
				}
				
				if ($v->{t_b} < $v->{t_e}){
				    $data[$exon]{t}{e} = $pos_t;
				}
				else{
				    $data[$exon]{t}{e} = $pos_t + 1;
				}
				$data[$exon]{q}{strand} = $v->{q_strand};
				$data[$exon]{t}{strand} = $v->{t_strand};
                						
				#print "ZEEEEEEEEE:$pos_q exon:$exon\n";

                                #fix 0 length exons in report
                                if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
                                    $data[$exon] = undef;
                                    $exon--; #undo iteration
                                }

				$exon++;
			}
                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }

                }
                elsif ($o->{state} eq 'I'){

                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }

                }
                elsif ($o->{state} eq 'S'){
		        if ($v->{q_b} < $v->{q_e}){
			    $data[$exon]{q}{b} = $pos_q + 1
				unless defined($data[$exon]{q}{b});
		        }
		        else {
			    $data[$exon]{q}{b} = $pos_q
				unless defined($data[$exon]{q}{b});
		        }

                        if ($v->{t_b} < $v->{t_e}){
			    $data[$exon]{t}{b} = $pos_t + 1
				unless defined($data[$exon]{t}{b});
                        }
                        else {
			    $data[$exon]{t}{b} = $pos_t
				unless defined($data[$exon]{t}{b});
                        }
		   
                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }

                }
                elsif ($o->{state} eq 'F'){

                        if ($v->{q_strand} == 1){
                                $pos_q += $o->{q};
                        }
                        else {
                                $pos_q -= $o->{q};
                        }

                        if ($v->{t_strand} == 1){
                                $pos_t += $o->{t};
                        }
                        else {
                                $pos_t -= $o->{t};
                        }

                }
		else {
			die "unknown state in Widget::exonerate::est2genome::get_hsp_coors!\n";
		}

        }
	
        if ($v->{q_b} < $v->{q_e}){
            $data[$exon]{q}{e} = $pos_q;
        }
        else{
            $data[$exon]{q}{e} = $pos_q + 1;
        }

        if ($v->{t_b} < $v->{t_e}){
            $data[$exon]{t}{e} = $pos_t;
        }
        else{
            $data[$exon]{t}{e} = $pos_t + 1;
        }
        $data[$exon]{q}{strand} = $v->{q_strand};
        $data[$exon]{t}{strand} = $v->{t_strand};

	#fix 0 length exons in report
	if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
	   delete($data[$exon]);
	   $exon--; #undo iteration
	}
	
	#my $new_data = fix_exon_coors(\@data);
	#return $new_data;
	return \@data;
}

#-------------------------------------------------------------------------------
sub fix_exon_coors {#no longer needed 7-26-2008
        my $data = shift;
        foreach my $exon (@{$data}){

		$exon->{q}{b}++ if $exon->{q}{b} == 0;
		$exon->{t}{b}++ if $exon->{t}{b} == 0;

                $exon->{q}{e}++ if $exon->{q}{e} == 0;
                $exon->{t}{e}++ if $exon->{t}{e} == 0;

                if ($exon->{t}->{strand} == 1){
                        $exon->{t}->{b}++;
                }  
                else {
                        $exon->{t}->{b}++;
                }
        }
}
#-------------------------------------------------------------------------------
sub get_model_order {
	my $v = shift;

	my $str = '';
	foreach my $o (@{$v->{operations}}){
		$str .= $o->{state};
	}

	my $type;
	if ($str =~ /3I5/ && $str =~ /5I3/){
		$type = 'mixed';
		warn "MIXED MODEL in Widget/est2genome!\n";
		warn "TELL MARK Y!\n";
		#sleep 5;
	}
	elsif ($str =~ /5I3/){
		$type = '5I3';
	}
	elsif ($str =~ /3I5/){
		$type = '3I5';
	}

	return $type;
}
#-------------------------------------------------------------------------------
sub assemble {
	my $bhd       = shift;
	my $bad       = shift;
	my $v         = shift;
	my $q_seq_len = shift;
	my $t_seq_len = shift;

	my $type = get_model_order($v);

	my $exons = get_exon_coors($v, $type);

	add_align_strs($bad, $exons);
	add_align_attr($exons, $q_seq_len, $t_seq_len);

        my $phat_hit = 
	new Bio::Search::Hit::PhatHit::protein2genome('-name'        => $v->{q_id},
						      '-description' => $v->{q_id},
						      '-algorithm'   => 'exonerate::protein2genome',
						      '-length'     => $q_seq_len,
						      );

        $phat_hit->queryLength($t_seq_len);

	my $i = -1;
	foreach my $exon (@{$exons}){
		$i++;
		#print "XXXXXXXXX $i XXXXXXXXXXXX\n";
		#next unless $i ==3;
		my $args = load_args($exon, $v);	
		my $hsp = new Bio::Search::HSP::PhatHSP::protein2genome(@{$args});
                   $hsp->queryName($v->{q_id});
		#-------------------------------------------------
                # setting strand because bioperl is all fucked up!
                #------------------------------------------------
                $hsp->{_strand_hack}->{query} = $exon->{t}->{strand};
                $hsp->{_strand_hack}->{hit}   = 0; #proteins are not stranded
		$hsp->{_indentical_hack}      = $exon->{identical};

		#print substr($hsp->query_string(), 0, 100)."\n";
		#print substr($hsp->homology_string(), 0, 100)."\n";
		#print substr($hsp->hit_string(), 0, 100)."\n";

		#PostData($exon);
		#$hsp->show();
		
		my $q_pos = $hsp->nE('query') -1 ;
;
		#my $h_pos = $hsp->equivalent_pos_in_alignment_partner('query', $q_pos);
		#my $q_char = $hsp->whatIsThere('query', $q_pos);
		#my $h_char  = $hsp->whatIsThere('hit', $h_pos);
		#my $q_test  =  $hsp->equivalent_pos_in_alignment_partner('hit', $h_pos);
		#print "q_pos:$q_pos q_char:$q_char h_pos:$h_pos h_char:$h_char q_test:$q_test\n";
		#die;

		$phat_hit->add_hsp($hsp);	
	}
	return $phat_hit;
}
#-------------------------------------------------------------------------------
sub split_aa_str {
        my $q_aa_str = shift;

        my $reg_ex = 'Target\s+Intron\s+\d+';

        my @q_aa_strs = split(/$reg_ex/, $q_aa_str);

        foreach my $str (@q_aa_strs){
                $str =~ s/\s*[<>]+\s*$//;
                $str =~ s/^\s*[<>]+\s*//;
        }

        return \@q_aa_strs;

}
#-------------------------------------------------------------------------------
sub add_align_strs {
	my $bad   = shift;
	my $exons = shift;

	my $q_aa_str = $bad->{q_aa_str};

	if ($q_aa_str =~ /Target Intron/){
	        #build array of exons for string
		my $q_aa_strs = split_aa_str($q_aa_str);

		#now build arrays for other string types
		my $m_strs    = [];
		my $t_aa_strs = [];
		my $t_nc_strs = [];

		my $i = 0; #report count
		my $o = 0; #offset
		foreach my $q_aa_part (@{$q_aa_strs}){
		    my $L = length($q_aa_part);

		    my $m_str    = substr($bad->{m_str},
					  $o,
					  $L,
					  );

		    my $t_aa_str = substr($bad->{t_aa_str},
					  $o,
					  $L,
					  );

		    my $t_nc_str = substr($bad->{t_nc_str},
					  $o,
					  $L,
					  );

		    if($L != 0) { #don't add 0 length exons
			push(@$m_strs, $m_str);
			push(@$t_aa_strs, $t_aa_str);
			push(@$t_nc_strs, $t_nc_str);
		    }

		    $o += $L + 28 + length($i + 1);
		    $i++;
		}

		@{$q_aa_strs} = grep {$_} @{$q_aa_strs}; #remove 0 length seqs

		#now correct splice site crossing features to conform to restraints
		#cooresponding to HSP objects and the GFF3 Gap attribute. 
		for(my $i = 0; $i < @{$q_aa_strs}; $i++){
		    #first move the amino acid that crosses the splice site to the next
		    #exon and replace it with gap characters, then add spaces to the next
		    #exon's nucleotide sequence to account for the moved amino acid
		    $q_aa_strs->[$i] =~ /\{([A-Za-z]+)\}$/;
		    if((my $rep = $1) && $i != @{$q_aa_strs} - 1){
			$q_aa_strs->[$i + 1] = $rep . $q_aa_strs->[$i + 1];
			$rep =~ s/./-/g;
			$q_aa_strs->[$i] =~ s/\{([A-Za-z]+)\}$/$rep/;
			$t_nc_strs->[$i + 1]= $rep . $t_nc_strs->[$i + 1];
		    }

		    #also move amino acid for translated target string
		    $t_aa_strs->[$i] =~ /\{([A-Za-z]+)\}$/;
                    if((my $rep = $1) && $i != @{$q_aa_strs} - 1){
                        $t_aa_strs->[$i + 1] = $rep . $t_aa_strs->[$i + 1];
                        $rep =~ s/./-/g;
                        $t_aa_strs->[$i] =~ s/\{([A-Za-z]+)\}$/$rep/;
		    }

		    #now duplicate homology string crossing the splice site to the next exon
		    #to account for the moved amino acid
		    $m_strs->[$i] =~ /\{(\|+)\}$/;
		    if((my $rep = $1) && $i != @{$q_aa_strs} - 1){
			$m_strs->[$i + 1] = $rep . $m_strs->[$i + 1];
		    }

		    #now safely replace '{' and '}' characters
		    $q_aa_strs->[$i] =~ s/[\{\}]//g;
		    $t_aa_strs->[$i] =~ s/[\{\}]//g;
		    $t_nc_strs->[$i] =~ s/[\{\}]//g;
		    $m_strs->[$i]    =~ s/[\{\}]//g;

		    #finally add string to exons
		    $exons->[$i]->{q_aa_str} = $q_aa_strs->[$i];
		    $exons->[$i]->{t_aa_str} = $t_aa_strs->[$i];
		    $exons->[$i]->{t_nc_str} = $t_nc_strs->[$i];
		    $exons->[$i]->{m_str}    = $m_strs->[$i];
		}
	}
	else {
		$exons->[0]->{q_aa_str} = $bad->{q_aa_str};
		$exons->[0]->{m_str}    = $bad->{m_str};
		$exons->[0]->{t_aa_str} = $bad->{t_aa_str};
		$exons->[0]->{t_nc_str} = $bad->{t_nc_str};
	}

	return $exons;
}
#-------------------------------------------------------------------------------
sub add_align_attr {
	my $exons     = shift;
	my $q_seq_len = shift;
	my $t_seq_len = shift;

	foreach my $e (@{$exons}){

		my $num_pipes = $e->{m_str} =~ tr/\|\+/\|\+/;
		my $num_con   = $e->{m_str} =~ tr/\|\+\!\:\./\|\+\!\:\./;
		my $m_str_len = length($e->{m_str});

		my $identical = $num_pipes/3;
		my $conserved = $num_con/ 3;

		my ($q_gaps) = $e->{q_aa_str} =~ tr/\<\-\>/\<\-\>/;
		my ($t_gaps) = $e->{t_nc_str} =~ tr/\-/\-/;

		my $q_hit_len = $e->{q_aa_str} =~ tr/A-Z/A-Z/;
		my $t_hit_len = $e->{t_nc_str} =~ tr/A-Z/A-Z/;

		$e->{identical} = $identical;
		$e->{conserved} = $conserved;
		$e->{q_gaps}    = $q_gaps;
		$e->{t_gaps}    = $t_gaps;
		$e->{q_hit_len} = $q_hit_len;
		$e->{t_hit_len} = $t_hit_len;	

		$e->{q_seq_len} = $q_seq_len;
        	$e->{t_seq_len} = $t_seq_len;

		
	}
}
#-------------------------------------------------------------------------------
sub load_args {
	my $exon = shift;
	my $v    = shift;

	my ($t_b, $t_e);
	if ($exon->{t}->{b} < $exon->{t}->{e}){
		$t_b = $exon->{t}->{b};
		$t_e = $exon->{t}->{e};
	}
	else {
		$t_b = $exon->{t}->{e};
                $t_e = $exon->{t}->{b};
	}

	my $q_b = $exon->{q}->{b};
	my $q_e = $exon->{q}->{e};

        my @args;

        push(@args, '-query_start');
        push(@args, $t_b);

	#reverse hit and target because exonerate reverses them from what is expected
        push(@args, '-query_seq');
        push(@args, $exon->{t_nc_str});

        push(@args, '-score');
        push(@args, $v->{score});

        push(@args, '-homology_seq');
        push(@args, $exon->{m_str});

        push(@args, '-hit_start');
        push(@args, $q_b);

	#reverse hit and target because exonerate reverses them from what is expected
        push(@args, '-hit_seq');
        push(@args, $exon->{q_aa_str});

        push(@args, '-hsp_length');
        push(@args, length($exon->{t_aa_str}));

        push(@args, '-identical');
	push(@args, $exon->{identical});

        push(@args, '-hit_length');
        push(@args, $exon->{q_hit_len});

        push(@args, '-query_name');
        push(@args, $v->{t_id});

        push(@args, '-algorithm');
        push(@args, 'exonerate::protein2genome');

        push(@args, '-bits');
        push(@args, $v->{score}); # bioperl hack!

        push(@args, '-evalue');
        push(@args, 'NA');

        push(@args, '-pvalue');
        push(@args, 'NA');

        push(@args, '-query_length');
        push(@args, $exon->{t_seq_len});

        push(@args, '-query_end');
        push(@args, $t_e);

        push(@args, '-conserved');
        push(@args, $exon->{conserved});

        push(@args, '-hit_name');
        push(@args, $v->{q_id});

        push(@args, '-hit_end');
        push(@args, $q_e);

        push(@args, '-query_gaps');
        push(@args, $exon->{t_gaps});

        push(@args, '-hit_gaps');
        push(@args, $exon->{q_gaps});

        push(@args, '-stranded');
        push(@args, 'NONE');



	return \@args;

}
#-------------------------------------------------------------------------------
sub get_gaps {
	my $type = shift;
	my $v    = shift;

	my $gaps = 0;
	foreach my $o (@{$v->{operations}}){
		next unless $o->{state} eq 'G';
		$gaps += $o->{$type};	
	}
	return $gaps;
}
#-------------------------------------------------------------------------------
sub parse {
	my $file = shift;
	my $q_seq_length = shift;
	my $t_seq_length = shift;	

	my $blocks = get_blocks($file);

	my $command;
	my $hostname;
	
	my @hits;
	my $BID = -1;
	foreach my $block  (@{$blocks}){
		$BID++;
		#next unless $BID ==7;
		my ($block_head_data, $block_align_data, $vulgarity);
		my @b = split(/\n/, $block);
		if ($b[0] =~ /Command line/){
			($command, $hostname) = parse_header(\@b);
		}
		elsif ($b[2] =~ /Query\:/ && $b[3] =~ /Target:/){

			($block_head_data, $block_align_data, $vulgarity) 
			= parse_block(\@b);

			my $phat_hit = assemble($block_head_data,
			                         $block_align_data,
			                         $vulgarity,
						 $q_seq_length,
			                         $t_seq_length,
			                         );
 
			push(@hits, $phat_hit);
		}
	}
	return \@hits;
}
#-------------------------------------------------------------------------------
sub parse_block_head {
	my $b = shift;

	
	my %data;
	foreach my $l (@{$b}){
        	chomp($l);
                if    ($l =~ /Query\:/){
                        ($data{q_id}) = $l =~ /\s+Query\:\s+(.*)$/;
                }
                elsif ($l =~ /Target\:/){
                        ($data{t_id}) = $l =~ /\s+Target\:\s+(.*)$/;
                }
                elsif ($l =~ /Model\:/){
                       ($data{model}) = $l =~ /\s+Model\:\s+(.*)$/;
                }
               elsif ($l =~ /Raw\s+score\:\s+/){
                       ($data{r_score}) = $l =~ /\s+Raw\s+score\:\s+(.*)$/;
                }
               elsif ($l =~ /Query\s+range\:/){
                       ($data{q_b}, $data{q_e}) = 
			$l =~ /\s+Query\s+range\:\s+(\d+)\s+\-\>\s+(\d+)$/;
                }
                elsif ($l =~ /Target\s+range\:/){
                       ($data{t_b}, $data{t_e}) = 
			$l =~ /\s+Target\s+range\:\s+(\d+)\s+\-\>\s+(\d+)$/;
                }
	}
	return \%data;
}
#-------------------------------------------------------------------------------
sub parse_block_tail {
	my $b = shift;

	my ($cigar, $vulgarity);
        foreach my $l (@{$b}){
                chomp($l);
                if ($l =~ /cigar\:/){
                        #$cigar = parse_cigar($l);
                }
                elsif ($l =~ /vulgar\:/){
                        $vulgarity = parse_vulgar($l);
                }


        }
	return $vulgarity;
}
#-------------------------------------------------------------------------------
sub parse_vulgar {
	my $l = shift;

	#print "VULGAR:$l\n";

	my @f = split(/\s+/, $l);

	my @left = splice(@f, 0, 10);

	my %vulgarity;
	while (defined (my $s = shift(@f))){
		my $q = shift @f;
		my $t = shift @f; 

		die "dead in parse_vulgar!\n"
                unless defined($s) 
		&&     defined($q) 
	        &&     defined($t);
		
		push(@{$vulgarity{operations}}, { 'state' => $s, 'q' => $q, 't' => $t});
	}

	
	$vulgarity{score}    = $left[9];
	$vulgarity{q_id}     = $left[1];
	$vulgarity{q_b}      = $left[2];
	$vulgarity{q_e}      = $left[3];
	$vulgarity{q_strand} = 1;

	$vulgarity{t_id}     = $left[5];
	$vulgarity{t_b}      = $left[6];
        $vulgarity{t_e}      = $left[7];
        $vulgarity{t_strand} = $left[8] eq '+' ? 1 : -1;
	
	
	 return \%vulgarity;
}
#-------------------------------------------------------------------------------
sub parse_block {
        my $b = shift;

	my @block_head = splice(@{$b},2, 7);
	my @block_tail;
	push(@block_tail, pop(@{$b}));
	push(@block_tail, pop(@{$b}));
	push(@block_tail, pop(@{$b}));	
	push(@block_tail, pop(@{$b}));

	my $block_head_data = parse_block_head(\@block_head);

	my $vulgarity       = parse_block_tail(\@block_tail);

	shift @{$b};
	shift @{$b};

	my $block_align_data = parse_block_align($b);

	return ($block_head_data, $block_align_data, $vulgarity); 
}
#-------------------------------------------------------------------------------
sub parse_block_align {
	my $b = shift;

	my %data;
	$data{q_aa_str} = '';
	$data{m_str}    = '';
	$data{t_aa_str} = '';
	$data{t_nc_str} = '';

	my ($w, $x, $y, $z) = '';
	while (defined(my $l = shift(@{$b}))){
		next unless $l =~ /\S+/;
		my ($lead, $q_aa) = $l =~ /(\s+\d+\s+\:\s)(\s*\S.*\S\s*)\s+\:/;
		   ($lead, $q_aa) = $l =~ /(\s+\d+\s+\:\s)(\s*\S\s*)\s+\:/
                   if !defined($lead) || !defined($q_aa);;

		my $o = length($lead);

		my $m    = substr(shift(@{$b}), $o);
		my $t_aa = substr(shift(@{$b}), $o); 
		my $t_nc = substr(shift(@{$b}), $o, length($q_aa));
 
		die "no q_aa\n" unless defined($q_aa);
		die "no m\n" unless defined($m);
		die "no t_aa\n" unless defined($t_aa);
		die "no t_nc\n" unless defined($t_nc);

		$data{q_aa_str} .= $q_aa;
		$data{m_str}    .= $m;
		$data{t_aa_str} .= $t_aa;
		$data{t_nc_str} .= $t_nc;

		#print $q_aa."\n";
		#print $m."\n";
		#print $t_aa."\n";
		#print $t_nc."\n";
		#print "\n";
	}
	
	return \%data;	
}
#-------------------------------------------------------------------------------
sub parse_header {
	my $b = shift;

	my ($command)  = $b->[0] =~ /Command line: \[(.+)\]/;
	my ($hostname) = $b->[1] =~ /Hostname\: \[(.+)\]/;

	return ($command, $hostname);
}
#-------------------------------------------------------------------------------
sub keepers {
        my $sio    = shift;
        my $params = shift;


	#print "XXXXXXXXXXX\n";

        my $result = $sio->next_result();

	#PostData($result);

        my @keepers;
        my $start = $result->hits();
        while(my $hit = $result->next_hit) {
		$hit->show();
=head;
                my $significance = $hit->significance();
                $significance = "1".$significance if  $significance =~ /^e/;
                $hit->queryLength($result->query_length);
                $hit->queryName($result->query_name);
                #next unless $significance < $params->{significance};
=cut;
                my @hsps;
                while(my $hsp = $hit->next_hsp) {
			print "start q:".$hsp->start('hit')."\n";
			print "end q:".$hsp->end('hit')."\n";
                        #$hsp->query_name($result->query_name);

                        #push(@hsps, $hsp) if $hsp->bits > $params->{hsp_bit_min};
                }
                $hit->hsps(\@hsps);
                push(@keepers, $hit) if $hit->hsps();
        }
        my $end     = @keepers;
        my $deleted = $start - $end;
        print STDERR "deleted:$deleted hits\n" unless $main::quiet;

        return \@keepers;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::RepeatMasker::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


