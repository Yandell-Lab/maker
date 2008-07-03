#------------------------------------------------------------------------
#----                        Widget::augustus                        ---- 
#------------------------------------------------------------------------
package Widget::augustus;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use augustus::PhatHit;
use augustus::PhatHsp;
use PhatHit_utils;

@ISA = qw(
	Widget
       );

my $OPT_F; # global option -f to force cleanup
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
#------------------------------------------------------------------------
sub get_aug_shot {
        my $seq           = shift;
        my $def           = shift;
        my $id            = shift;
        my $the_void      = shift;
        my $set           = shift;
        my $pred_flank    = shift;
        my $pred_command  = shift;
        my $q_id          = shift;
	   $OPT_F         = shift;

        my ($shadow_seq, $strand, $offset, $xdef) =
            prep_for_genefinder($seq, $set, $pred_flank, $id, $q_id);

        my $shadow_fasta = Fasta::toFasta($def." $id offset:$offset",
                                          $shadow_seq,
                                         );


        my $alt_aug_command;
        if ($strand == 1){
                $alt_aug_command = $pred_command.' --strand=forward';
        }
        else {
                 $alt_aug_command = $pred_command.' --strand=backward';
        }

        my $gene_preds = augustus($$shadow_fasta,
                                  $the_void,
                                  $id,
                                  $strand,
                                  $offset,
                                  $xdef,
                                  $alt_aug_command,
                                 );


        return ($gene_preds, $strand);
}
#------------------------------------------------------------------------
sub prep_for_genefinder {
        my $seq    = shift;
        my $set    = shift;
        my $flank  = shift;
        my $seq_id = shift;
        my $q_id   = shift;

        my $gomiph = $set->{gomiph};
        my $mia    = $set->{mia};
        my $augs  = $set->{preds};
        my $ests   = $set->{ests};
        my @t_data;

        push(@t_data, @{$gomiph})  if defined($gomiph);
        push(@t_data, @{$augs})   if defined($augs);
        push(@t_data, $mia)        if defined($mia);
        push(@t_data, @{$ests})    if defined($ests);

        my $p_set_coors = PhatHit_utils::get_hsp_coors($gomiph, 'query');

        my $n_set_coors =
        defined($ests) ? PhatHit_utils::get_hsp_coors($ests, 'query')
                      : [];

        my @coors;
        my $plus  = 0;
        my $minus = 0;


        my $least;
        my $most;
        foreach my $hit (@t_data){
                foreach my $hsp ($hit->hsps()){
                        my $s = $hsp->start('query');
                        my $e = $hsp->end('query');

                        $least = $s if !defined($least) || $s < $least;
                        $most  = $e if !defined($most)  || $e > $most;

                        if ($hsp->strand('query') == 1) {
                                $plus++;
                        }
                        else {
                                $minus++;
                        }
                }
        }

        my @span_coors = [$least, $most];

        my $p = Shadower::getPieces($seq, \@span_coors, $flank);

        my $final_seq = $p->[0]->{piece};

        my $offset    = $p->[0]->{b} - 1;

        my $strand = $plus > $minus ? 1 : -1;

        my $i_flank = 2;

        my $xdef = get_xdef($seq,
                            $p_set_coors,
                            $n_set_coors,
                            $strand == 1 ? '+' : '-',
                            $offset,
                            $i_flank,
                            $q_id,
                           );


        return (\$final_seq, $strand, $offset, $xdef);
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
#------------------------------------------------------------------------
sub augustus {
        my $fasta      = shift;
        my $the_void   = shift;
        my $seq_id     = shift;
        my $strand     = shift;
        my $offset     = shift;
        my $xdef       = shift;
        my $command    = shift;

        my $aug_keepers = [];

        my $file_name = "$the_void/$seq_id\.auto_annotator\.$offset\.augustus.fasta";

        my $o_file    = "$the_void/$seq_id\.$offset\.auto_annotator\.augustus";

        my $xdef_file = "$the_void/$seq_id\.$offset\.auto_annotator\.xdef\.augustus";

        write_xdef_file($xdef, $xdef_file) if defined $xdef;

        FastaFile::writeFile(\$fasta, $file_name);

        $command .= ' --hintsfile='.$xdef_file if -e $xdef_file;

            $command .= " $file_name";
            $command .= " > $o_file";


        if (-e $o_file && ! $OPT_F){
                print STDERR "re reading augustus report.\n";
                print STDERR "$o_file\n";
        }
        else {
                print STDERR "running  augustus.\n";
		my $w = new Widget::augustus();
                $w->run($command);
        }


        my %params;
           $params{min_exon_score}  = -100;
           $params{min_gene_score}  = -20;

        my $keepers = parse($o_file,
                           \%params,
                            $file_name,
                           );

        PhatHit_utils::add_offset($keepers,
                                  $offset,
                                  );

        return $keepers;
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
		die "you must give Widget::augustus a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_xdef {
       my $seq       = shift;
        my $p_coors  = shift;
        my $n_coors = shift;
        my $s       = shift;
        my $offset  = shift;
        my $i_flank = shift;
	my $q_id    = shift;

	my $group = 'gb|bogus;';

        my $p_pieces = Shadower::getPieces($seq, $p_coors, 0);

        my @xdef;
        for (my $i = 0; $i < @{$p_pieces}; $i++){
                my $p = $p_pieces->[$i];
                my $c_b = $p->{b} - $offset;
                my $c_e = $p->{e} - $offset;

		my $l  = "$q_id\tPROTEIN\tCDSpart";
		   $l .= "\t".$c_b."\t".$c_e."\t"."1e-1000"."\t".$s;
                   $l .= "\t".'.'."\t".'group=p_pieces;'."source=P";

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

                my $l  = "$q_id\tEST-INTRON\tintronpart";
                   $l .= "\t".$i_b."\t".$i_e."\t"."1e-1000"."\t".$s;
                   $l .= "\t".'.'."\t".'group=n_pieces;'."source=E";

                push(@xdef, $l);
        }

        for (my $i = 0; $i < @{$n_pieces}; $i++){
                my $p = $n_pieces->[$i];
                my $e_b = $p->{b} - $offset;
                my $e_e = $p->{e} - $offset;

                my $l  = "$q_id\tEST-EXON\texonpart";
                   $l .= "\t".$e_b."\t".$e_e."\t"."1e-1000"."\t".$s;
                   $l .= "\t".'.'."\t".'group=n_peices;'."source=E";

                push(@xdef, $l);
        }

        return \@xdef;

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
sub parse_gene {
	my $stuff = shift;



	shift(@{$stuff});

	my %this_gene;
	my $i = 0;
	my $gene_name;
	my $score;
	foreach my $datum (@{$stuff}){

		next if $datum =~ /### end gene/;

		my @fields = split(/\s+/, $datum);


		if ($fields[1] eq 'AUGUSTUS' && $fields[2] eq 'gene'){
			$gene_name = $fields[0].'-'.$fields[8];
			$score     = $fields[5];

			$this_gene{name}  = $gene_name;
			$this_gene{score} = $score;
		}
		if ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'CDS' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

			push(@{$this_gene{CDS}}, {'b'      => $fields[3],
			                          'e'      => $fields[4],
					          'strand' => $strand,
					          'frame'  => $fields[7],
						  'score'  => $score,
						  'type'   => 'CDS',
				                  });

		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'start_codon' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

			$this_gene{start_codon} = {'b'      => $fields[3],
						   'e'      => $fields[4], 
		                                  'strand' => $strand,
			                          'frame'  => $fields[7], 
		                                  'score'  => $score,
						  'type'   => 'start_codon',
                                                  };

		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'stop_codon' ){
			my $strand = $fields[6] eq '+' ? 1 : -1;

                        $this_gene{stop_codon} = {'b'      => $fields[3],
			                          'e'      => $fields[4],
                                                  'strand' => $strand,
                                                  'frame'  => $fields[7],
                                                  'score'  => $score,
						  'type'   => 'stop_codon',
                                                  };

		}
		elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'gene'){
		}
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'transcript'){
                }
                elsif ($fields[1] eq 'AUGUSTUS' && $fields[2]  eq 'intron'){
                }
                elsif ($fields[1] eq 'protein' && $fields[2]  eq 'sequence'){
                }
                elsif ($fields[0] eq '#'){
                }
		else {
			PostData(\@fields);	
			die "dead in Widget::augustus::parse_gene unknown feature type:".$fields[2]."\n";
		}
		
	}

	return \%this_gene; 
}
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

	$/ = "### gene";

	my $fh = new FileHandle();
	   $fh->open($report);

	<$fh>; # shift off the header

	my @genes;
	while (my $line = <$fh>){
		chomp($line);

		my @stuff = split(/\n/, $line);

		my $gene = parse_gene(\@stuff);
		push(@genes, $gene);
	}

	$/= "\n";

	$fh->close();

	refactor(\@genes);
	my $phat_hits = load_phat_hits($q_name, $q_seq, \@genes, $params);
	return $phat_hits;
}
#-------------------------------------------------------------------------------
sub refactor {
	my $genes = shift;

	foreach my $g (@{$genes}){
		my $CDS = $g->{CDS};
		my $start_codon = $g->{start_codon} || undef;
		my $stop_codon  = $g->{stop_codon}  || undef;
		foreach my $cds (reverse sort sort_cdss @{$CDS}){
			
	                if ($cds->{strand} == 1){
                        	$cds->{e} = $stop_codon->{e};
                	}
                	else {
                        	$cds->{b} = $stop_codon->{b};
                	}

			last;
		}
	}
}
#-------------------------------------------------------------------------------
sub sort_cdss {
	if ($a->{strand} == 1){
		$a->{b} <=> $b->{b};
	}
	else {
		$b->{b} <=> $a->{b}
	}
}
#-------------------------------------------------------------------------------
sub load_phat_hits {
	my $q_name  = shift;
	my $q_seq   = shift;
	my $genes   = shift;
	my $params  = shift;

	my $q_len = length($$q_seq);

	my @keepers;
	foreach my $g (@{$genes}){
	
		my $total_score = $g->{score};

		my %hsps;
		my $i = 0;
		my $f = new augustus::PhatHit('-name'         => $g->{name},
                                              '-description'  => 'NA',
                                              '-algorithm'    => 'augustus',
                                              '-length'       => $q_len,
					      '-score'        => $total_score,
                                             );

                $f->queryLength($q_len);

		my $hit_start = 1;
		foreach my $exon (@{$g->{CDS}}){
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
        		push(@args, 'AUGUSTUS');

        		push(@args, '-bits');
        		push(@args, 'NA');

        		push(@args, '-evalue');
        		push(@args, 'NA');

        		push(@args, '-pvalue');
        		push(@args, $exon->{score});

        		push(@args, '-query_length');
        		push(@args, $q_len);

        		push(@args, '-query_end');
        		push(@args, $exon->{e});

        		push(@args, '-conserved');
        		push(@args, abs($exon->{e} - $exon->{b}) + 1);

        		push(@args, '-hit_name');
        		push(@args, $g->{name}.":".$i.":".$exon->{type});

        		push(@args, '-hit_end');
        		push(@args, $hit_end);

        		push(@args, '-query_gaps');
        		push(@args, 0);

        		push(@args, '-hit_gaps');
        		push(@args, 0);

        		my $hsp = new augustus::PhatHsp(@args);
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

        print STDERR "Widget::augustus::AutoLoader called for: ",
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


