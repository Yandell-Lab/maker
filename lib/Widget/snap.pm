#------------------------------------------------------------------------
#----                        Widget::snap                            ---- 
#------------------------------------------------------------------------
package Widget::snap;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use snap::PhatHit;
use snap::PhatHsp;
@ISA = qw(
	Widget
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
sub run {
	my $self    = shift;
	my $command = shift;

	if (defined($command)){
		$self->print_command($command);
		system("$command");
	}
	else {
		die "you must give Widget::snap a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
        my $report = shift;
        my $params = shift;
	my $q_file = shift;

	my $iterator = new Iterator::Fasta($q_file);
        my $fasta = $iterator->nextEntry();

        my $def     = Fasta::getDef($fasta);
        my $q_seq   = Fasta::getSeq($fasta);

        my ($q_name)  = $def =~ /^>(.+)/;

	my $fh = new FileHandle();
	   $fh->open($report);

	<$fh>; # shift off the header
	my %g;
	my $i = -1;
	while (my $line = <$fh>){
		chomp($line);
		my @stuff = split(/\s+/, $line);

		$i = 0  if  !defined($g{$stuff[8]});

		my $b = $stuff[1];
		my $e = $stuff[2];
		#($b, $e) = ($e, $b) if $stuff[3] eq '-';

		$g{$stuff[8]}[$i]{type}     = $stuff[0];
		$g{$stuff[8]}[$i]{b}        = $stuff[1];
		$g{$stuff[8]}[$i]{e}        = $stuff[2];
		$g{$stuff[8]}[$i]{strand}   = $stuff[3] eq '+' ? 1 : -1;
		$g{$stuff[8]}[$i]{score}    = $stuff[4];
		$g{$stuff[8]}[$i]{five}     = $stuff[5];
		$g{$stuff[8]}[$i]{three}    = $stuff[6];
		$g{$stuff[8]}[$i]{frame}    = $stuff[7];	
		$i++;
	}
	$fh->close();
	my $phat_hits = load_phat_hits($q_name, $q_seq, \%g, $params);
	return $phat_hits;
}
#-------------------------------------------------------------------------------
sub get_xdef {
	my $seq     = shift;
	my $p_coors  = shift;
	my $n_coors = shift;
	my $s       = shift;
	my $c_u     = shift;
	my $i_u     = shift;
	my $offset  = shift;
	my $i_flank = shift;	

        my $p_pieces = Shadower::getPieces($seq, $p_coors, 0);

        my @xdef;
	push(@xdef, ">xdef file offset: $offset i_flank:$i_flank");
        for (my $i = 0; $i < @{$p_pieces}; $i++){
                my $p = $p_pieces->[$i];
                my $c_b = $p->{b} - $offset;
                my $c_e = $p->{e} - $offset;
		my $l = "Coding\t".$c_b."\t".$c_e;
		$l .= "\t$s\t$c_u\t\.\t\.\t\.\tADJ";
                push(@xdef, $l);
        }

	my $n_pieces = Shadower::getPieces($seq, $n_coors, 0);

	return \@xdef if !defined($n_pieces->[1]);

        my @intron;

	my $num = @{$n_pieces};
        for (my $i = @{$n_pieces} - 1; $i > 0;  $i--){
                my $p_r = $n_pieces->[$i];
                my $p_l = $n_pieces->[$i - 1];
                my $i_b = ($p_l->{e} - $offset) + $i_flank;
                my $i_e = ($p_r->{b} - $offset) - $i_flank;
	
		next if abs($i_b - $i_e) < 2*$i_flank;
		next if abs($i_b - $i_e) < 25;

		my $l = "Intron\t".$i_b."\t".$i_e;
                $l .= "\t$s\t$i_u\t\.\t\.\t\.\tADJ";

		push(@xdef, $l);

		my $c = "Coding\t".$i_b."\t".$i_e;
		$c .= "\t$s\t-100\t\.\t\.\t\.\tADJ";
		push(@xdef, $c);


        }
	return \@xdef;
}
#-------------------------------------------------------------------------------
sub adjust_xdef {
	my $xdef       = shift;
	my $offset     = shift;
	my $chunk_size = shift;
	my @new_xdef;
	my $i = 0;

	my $end = $offset + $chunk_size;

	foreach my $l (@{$xdef}){
		if ($i == 0){
			push(@new_xdef, $l);
		}
		else {
			my @fields = split (/\s+/, $l);
			if ($fields[1] >=  $offset && $fields[2] <= $end){
				push(@new_xdef, $l);
			}
			else {
			} 
		}
			
		$i++;
	}
	return \@new_xdef;
}
#-------------------------------------------------------------------------------
sub write_xdef_file {
	my $xdef = shift;
	my $loc  = shift;

	my $fh = new FileHandle();
	   $fh->open(">$loc") || die "couldn't open $loc";	

	foreach my $l (@{$xdef}){
		print $fh $l."\n";
	}
	$fh->close();

}
#-------------------------------------------------------------------------------
sub get_exon_seq {
	my $exon  = shift;
	my $q_seq = shift;

       my $e_b = $exon->{b};
       my $e_e = $exon->{e};
        
        my $length = abs($e_e - $e_b) + 1;
        
        ($e_b, $e_e) = ($e_e, $e_b) if $exon->{strand} == -1;

        my $e_seq = substr($$q_seq, $e_b - 1, $length);

        $e_seq = Fasta::revComp($e_seq) if $exon->{strand} == -1;

	return $e_seq;
}
#-------------------------------------------------------------------------------
sub total_score {
	my $g = shift;

	my $total = 0;
	foreach my $exon (@{$g}){
		$total += $exon->{score};
	}
	return $total;
}
#-------------------------------------------------------------------------------
sub load_phat_hits {
	my $q_name  = shift;
	my $q_seq   = shift;
	my $g       = shift;
	my $params  = shift;

	my $q_len = length($$q_seq);

	my @keepers;
	foreach my $gene (keys %{$g}){
	
		my $total_score = total_score($g->{$gene});

		my %hsps;
		my $i = 0;
		my $f = new snap::PhatHit('-name'         => $gene,
                                          '-description'  => 'NA',
                                          '-algorithm'    => 'snap',
                                          '-length'       => $q_len,
					  '-score'        => $total_score,
                                          );

                $f->queryLength($q_len);

		my $hit_start = 1;
		foreach my $exon (@{$g->{$gene}}){
			my @args;
			my $exon_seq = get_exon_seq($exon, $q_seq); # revcomped!
			my $hit_end = abs($exon->{e} - $exon->{b}) + $hit_start;

 			push(@args, '-query_start');
        		push(@args, $exon->{b});

       			push(@args, '-query_seq');
       			push(@args, $exon_seq);

        		push(@args, '-score');
        		push(@args, $exon->{score});
       			push(@args, '-homology_seq');
       			push(@args, $exon_seq);

        		push(@args, '-hit_start');
        		push(@args, $hit_start);

       			push(@args, '-hit_seq');
       			push(@args, $exon_seq);

        		push(@args, '-hsp_length');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-identical');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_length');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-query_name');
        		push(@args, $q_name);

        		push(@args, '-algorithm');
        		push(@args, 'SNAP');

        		push(@args, '-bits');
        		push(@args, 'NA');

        		push(@args, '-evalue');
        		push(@args, 'NA');

        		push(@args, '-pvalue');
        		push(@args, 'NA');

        		push(@args, '-query_length');
        		push(@args, $q_len);

        		push(@args, '-query_end');
        		push(@args, $exon->{e});

        		push(@args, '-conserved');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_name');
        		push(@args, $gene.":".$i.":".$exon->{type});

        		push(@args, '-hit_end');
        		push(@args, $hit_end);

        		push(@args, '-query_gaps');
        		push(@args, 0);

        		push(@args, '-hit_gaps');
        		push(@args, 0);

        		my $hsp = new snap::PhatHsp(@args);
        		   $hsp->queryName($q_name);
        		#-------------------------------------------------
        		# setting strand because bioperl is all messed up!
        		#------------------------------------------------
			if ($exon->{strand} == 1 ){
				$hsp->{_strand_hack}->{query} = 1;
				$hsp->{_strand_hack}->{hit}   = 1;
			}
			else {
				$hsp->{_strand_hack}->{query} = -1;
				$hsp->{_strand_hack}->{hit}   =  1;
			}
			$f->add_hsp($hsp);

			$hit_start += abs($exon->{e} - $exon->{b}) + 1; 
			
			$i++;
		}
              	push(@keepers, $f);

        }
        return keepers(\@keepers, $params);

}
#-------------------------------------------------------------------------------
sub keepers {
        my $genes   = shift;
	my $params = shift; 


	my $start = @{$genes};
        my @keepers;
	foreach my $phat_hit (@{$genes}){
               	my @hsps;
		my $total_score = 0;
              while(my $hsp = $phat_hit->next_hsp) {
                        push(@hsps, $hsp)
			if $hsp->score() > $params->{min_exon_score};
			
			$total_score += $hsp->score();
             }
             $phat_hit->hsps(\@hsps);
             push(@keepers, $phat_hit) 
	     if ($phat_hit->hsps() && $total_score > $params->{min_gene_score});	

	}
        my $end     = @keepers;
        my $deleted = $start - $end;
        print STDERR "deleted:$deleted genes\n";

        return \@keepers;
}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::snap::AutoLoader called for: ",
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


