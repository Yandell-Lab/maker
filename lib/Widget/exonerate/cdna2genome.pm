#------------------------------------------------------------------------
#----                   Widget::exonerate::cdna2genome                ---- 
#------------------------------------------------------------------------
package Widget::exonerate::cdna2genome;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget::exonerate;
use Bio::Search::Hit::PhatHit::cdna2genome;
use Bio::Search::HSP::PhatHSP::cdna2genome;
use PhatHit_utils;
use exonerate::splice_info;
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

	#checks if exonerate really finished
	unless(@chunks && grep {/completed exonerate analysis/} $chunks[-1]){
	    unlink($file);
	    die "ERROR: The file $file appears to be incomplete\n".
		"MAKER will need to delete the file, before trying again\n\n";
	}

	return \@chunks;
}
#-------------------------------------------------------------------------------
sub get_exon_coors_r {
        my $v      = shift;
        my $type  = shift;

        my $pos_q = $v->{q_b};
        my $pos_t = $v->{t_b};

        my $exon = 0;

        my @data;
        foreach my $o (@{$v->{operations}}){
                if    ($o->{state} eq 'M' || $o->{state} eq 'C'){

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
                elsif ($o->{state} eq 'N'){
                        die "dead in cdna2genomic::get_hsp_coors:N\n";
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

				#fix 0 length exons in report
				if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
				    $data[$exon] = undef;
				    $exon--; #undo iteration
				}

				$exon++;
			}
                        #print "ZEEEEEEEEE:$pos_q exon:$exon\n";

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
                        die "unknown state in Widget::exonerate::cdna2genome::get_hsp_coors!\n";
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
sub get_exon_coors_f {
	my $v    = shift;
	my $type = shift;

	my $pos_q = $v->{q_b};
	my $pos_t = $v->{t_b};

	my $exon = 0;

	my @data;
        foreach my $o (@{$v->{operations}}){
		if    ($o->{state} eq 'M' || $o->{state} eq 'C'){

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
                elsif ($o->{state} eq 'N'){
			die "dead in cdna2genomic::get_hsp_coors:N\n";
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
			die "unknown state in Widget::exonerate::cdna2genome::get_hsp_coors!\n";
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
sub fix_exon_coors { #no longer needed 7-26-2008
	my $data = shift;

	my @new_data;
	foreach my $exon (@{$data}){
		next unless defined($exon->{q}->{strand}); # deperate hack on 52006
		$exon->{q}{b}++ if $exon->{q}{b} == 0;
		$exon->{t}{b}++ if $exon->{t}{b} == 0;

                $exon->{q}{e}++ if $exon->{q}{e} == 0;
                $exon->{t}{e}++ if $exon->{t}{e} == 0;

		if ($exon->{q}->{strand} == 1){
			$exon->{t}{b}++;
		}
		else {
			$exon->{t}{b}++;
		}	
		push(@new_data, $exon);
	}
	return \@new_data;
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
		warn "MIXED MODEL in Widget/cdna2genome!\n";
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

	my $exons;
	if ($v->{q_strand} == 1){
		$exons = get_exon_coors_f($v, $type);
	}
	else {
		$exons = get_exon_coors_r($v, $type);
	}

	add_align_strs($bad, $exons);
	add_align_attr($exons, $q_seq_len, $t_seq_len);

        my $phat_hit = 
	new Bio::Search::Hit::PhatHit::cdna2genome('-name'        => $v->{q_id},
						   '-description' => $v->{q_id},
						   '-algorithm'   => 'exonerate::cdna2genome',
						   '-length'     => $q_seq_len,);

        $phat_hit->queryLength($t_seq_len);

	my $i = -1;
	foreach my $exon (@{$exons}){
		$i++;
		#print "XXXXXXXXX $i XXXXXXXXXXXX\n";
		#next unless $i ==2;
		my $args = load_args($exon, $v);	
		my $hsp = new Bio::Search::HSP::PhatHSP::cdna2genome(@{$args});
                   $hsp->queryName($v->{q_id});
		#-------------------------------------------------
                # setting strand because bioperl is all fucked up!
                #------------------------------------------------
                $hsp->{_strand_hack}->{query} = $exon->{t}->{strand};
                $hsp->{_strand_hack}->{hit}   = $exon->{q}->{strand};
		$hsp->{_indentical_hack}      = $exon->{identical};

		#print substr($hsp->query_string(), 0, 100)."\n";
		#print substr($hsp->homology_string(), 0, 100)."\n";
		#print substr($hsp->hit_string(), 0, 100)."\n";

		#my $q_pos = $hsp->strand('query') == 1  ? $hsp->nE('query') - 1 
                #                                        : $hsp->nE('query') + 1;

		#my $h_pos  = $hsp->equivalent_pos_in_alignment_partner('query', $q_pos);
		#my $q_char = $hsp->whatIsThere('query', $q_pos);
		#my $h_char = $hsp->whatIsThere('hit', $h_pos);
		#my $q_test =  $hsp->equivalent_pos_in_alignment_partner('hit', $h_pos);

		#print "q_pos:$q_pos q_char:$q_char h_pos:$h_pos h_char:$h_char q_test:$q_test\n";

		#PostData($exon);
		#$hsp->show();

		$hsp->donor($exon->{donor});
		$hsp->acceptor($exon->{acceptor});

		$phat_hit->add_hsp($hsp);	
	}

	if (exonerate::splice_info::needs_to_be_revcomped($phat_hit)){
	    $phat_hit = PhatHit_utils::copy($phat_hit, 'both');
	}


	return $phat_hit;
}
#-------------------------------------------------------------------------------
sub split_nc_str {
	my $q_nc_str = shift;

	my $reg_ex = 'Target\s+Intron\s+\d+';

	my @q_nc_strs = split(/$reg_ex/, $q_nc_str);

	foreach my $str (@q_nc_strs){
	    $str =~ s/\s*[<>]+\s*$//;
	    $str =~ s/^\s*[<>]+\s*//;
	}

	return \@q_nc_strs;

}
#-------------------------------------------------------------------------------
sub get_donor {
	my $t_nuc_str = shift;
	my $offset    = shift;
	my $q_nuc_str = shift;

	my $d;
	if (length($t_nuc_str) == $offset + length($q_nuc_str)){
		$d = undef;
	}
	else {
		$d = substr($t_nuc_str, $offset + length($q_nuc_str) , 2 );
	}

	die "Messed up donor:$d in Widget::exonerate::cdna2genome::get_donor!\n"
	if (defined($d) &&  lc($d) ne $d);

	return $d
}
#-------------------------------------------------------------------------------
sub get_acceptor {
        my $t_nuc_str = shift;
        my $offset    = shift;
        my $q_nuc_str = shift;

	my $a;
	if ($offset == 0){
		$a = undef;
	}
	else {
		$a = substr($t_nuc_str, $offset -2, 2 );
	}

        die "Messed up acceptor:$a in Widget::exonerate::cdna2genome::get_acceptor!\n"
        if defined($a) && lc($a) ne $a;

        return $a 
}
#-------------------------------------------------------------------------------
sub add_align_strs {
	my $bad   = shift;
	my $exons = shift;

	my $q_nc_str = $bad->{q_nc_str};

	if ($q_nc_str =~ /Target Intron/){
		my $q_nc_strs = split_nc_str($q_nc_str);

		my $i = 0; #report count
		my $j = 0; #practical count
		my $o = 0; #offset
		foreach my $q_nc_part (@{$q_nc_strs}){
			my $L = length($q_nc_part);

			my $m_str    = substr($bad->{m_str},    
				              $o, 
			                      $L,
			                      );

			my $t_nc_str = substr($bad->{t_nc_str}, 
			                      $o,
			                      $L,
			                      );

			my $L2 = 1;
			if(@{$q_nc_strs} != @$exons){
			    $L2 = length($t_nc_str) - scalar($t_nc_str =~ tr/\-/\-/);
			}

			my $don = get_donor($bad->{t_nc_str}, $o, $q_nc_part);
			my $acc = get_acceptor($bad->{t_nc_str}, $o, $q_nc_part);

			if($L != 0 && $L2 != 0) { #don't add 0 length exons
			    $exons->[$j]->{donor}    = $don;
			    $exons->[$j]->{acceptor} = $acc;
			    $exons->[$j]->{q_nc_str} = $q_nc_part;
			    $exons->[$j]->{m_str}    = $m_str;
			    $exons->[$j]->{t_nc_str} = $t_nc_str;
			    $j++;
			}

			$o += $L + 28 + length($i + 1);
			$i++;
		}
	}
	else {
		$exons->[0]->{q_nc_str} = $bad->{q_nc_str};
		$exons->[0]->{m_str}    = $bad->{m_str};
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

		my $identical = $num_pipes;
		my $conserved = $num_con;

		my ($q_gaps) = $e->{q_nc_str} =~ tr/\-/\-/;
		my ($t_gaps) = $e->{t_nc_str} =~ tr/\-/\-/;

		my $q_hit_len = $e->{q_nc_str} =~ tr/A-Z/A-Z/;
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

        push(@args, '-query_seq');
        push(@args, $exon->{t_nc_str});

        push(@args, '-score');
        push(@args, $v->{score});

        push(@args, '-homology_seq');
        push(@args, $exon->{m_str});

        push(@args, '-hit_start');
        push(@args, $q_b);

        push(@args, '-hit_seq');
        push(@args, $exon->{q_nc_str});

        push(@args, '-hsp_length');
        push(@args, length($exon->{t_nc_str}));

        push(@args, '-identical');
        #push(@args, $exon->{identical});
	push(@args, $exon->{identical});

        push(@args, '-hit_length');
        push(@args, $exon->{q_hit_len});

        push(@args, '-query_name');
        push(@args, $v->{t_id});

        push(@args, '-algorithm');
        push(@args, 'exonerate::cdna2genome');

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
		#next unless $BID ==1;
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

	#filter for best match as they can be mapped to both strands	
	my @keepers;
	@hits = sort {$b->score <=> $a->score} @hits;
	foreach my $h (@hits){
	    my $bad;
	    foreach my $k (@keepers){
		if(compare::compare($h->start('query'), $h->end('query'),
				    $k->start('query'), $k->end('query'))
		  ){
		    $bad = 1;
		    last;
		}
	    }
	    push(@keepers, $h) unless($bad);
	}

	return \@keepers;
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
	$vulgarity{q_strand} = $left[4] eq '+' ? 1 : -1;

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
		my ($lead, $q_nc) = $l =~ /(\s+\d+\s+\:\s)(\s*\S.*\S\s*)\s+\:/;
		   ($lead, $q_nc) = $l =~ /(\s+\d+\s+\:\s)(\s*\S\s*)\s+\:/
		   if !defined($lead) || !defined($q_nc);

		#print "$lead\n";
		#print "$q_nc\n";

		my $o = length($lead);

		my @buf;
		while(my $l2 = shift(@{$b})){
		    push(@buf, $l2);
		    last if($l2 =~ /(\s+\d+\s+\:\s)(\s*\S.*\S\s*)\s+\:/ ||
			    $l2 =~ /(\s+\d+\s+\:\s)(\s*\S\s*)\s+\:/);
		}
		
		my $m;
		my $t_nc;
		
		if(@buf == 2){
		    $m    = substr($buf[0], $o);
		    $t_nc = substr($buf[1], $o, length($q_nc));
		}
		else{
		    unshift(@buf, ' 'x(length($q_nc)+$o)) if(@buf == 3); #missing extra query line
		    my @q_nc1 = split('', $q_nc);
		    my @q_nc2 = split('', substr($buf[0], $o));
		    $m = substr($buf[1], $o);
		    my @t_nc1 = split('', substr($buf[2], $o));
		    my @t_nc2 = split('', substr($buf[3], $o));

		    #fix q_nc to remove protein codon
		    undef $q_nc;
		    for(my $i = 0; $i < @q_nc1; $i++){
			if(defined($q_nc1[$i]) && $q_nc1[$i] ne ' '){
			    $q_nc .= $q_nc1[$i];
			}
			else{
			    $q_nc .= $q_nc2[$i];
			}
		    }

		    #fix t_nc1 to remove protein codon
		    undef $t_nc;
		    for(my $i = 0; $i < @q_nc1; $i++){ #restricted by q_nc
			if(defined($t_nc2[$i]) && $t_nc2[$i] ne ' '){
			    $t_nc .= $t_nc2[$i];
			}
			else{
			    $t_nc .= $t_nc1[$i];
			}
		    }
		}

		#print "$m\n";
		#print "$t_nc\n";
		#print "XXXXXXXXXXXXXXXXXX\n";

		die "no q_nc\n" unless defined($q_nc);
		die "no m\n"    unless defined($m);
		die "no t_nc\n" unless defined($t_nc);

		$q_nc =~ s/(\{|\})//g;
		$m =~ s/(\{|\})//g;
		$t_nc =~ s/(\{|\})//g;
		$data{q_nc_str} .= $q_nc;
		$data{m_str}    .= $m;
		$data{t_nc_str} .= $t_nc;

		#print $q_nc."\n";
		#print $m."\n";
		#print $t_nc."\n";
		#print "\n";
	}

	$data{q_nc_str} =~ s/\<\-\>/\-\-\-/g;
	$data{t_nc_str} =~ s/\<\-\>/\-\-\-/g;
	
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
		#$hit->show();
=head;
                my $significance = $hit->significance();
                $significance = "1".$significance if  $significance =~ /^e/;
                $hit->queryLength($result->query_length);
                $hit->queryName($result->query_name);
                #next unless $significance < $params->{significance};
=cut;
                my @hsps;
                while(my $hsp = $hit->next_hsp) {
			#print "start q:".$hsp->start('hit')."\n";
			#print "end q:".$hsp->end('hit')."\n";
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


