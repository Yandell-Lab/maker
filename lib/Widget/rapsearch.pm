#------------------------------------------------------------------------
#----                       Widget::rapsearch                        ---- 
#------------------------------------------------------------------------
package Widget::rapsearch;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::Search::Hit::PhatHit::rapsearch;
use Bio::Search::HSP::PhatHSP::rapsearch;
use Exporter;
use PostData;
use FileHandle;
use Widget;
use PhatHit_utils;
use IPC::Open3;
use Symbol;

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
      my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
      my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);

      my $all_err;
      {
	  local $/ = \1;
	  while (my $line = <$CHLD_ERR>){
	      $all_err .= $line;
	      print STDERR $line unless($main::quiet);
	  }
      }
      waitpid $pid, 0;
      if ($? != 0 && ! ignorable($all_err, $command)){
	 #run a second time on failure
	 sleep 15;
	 $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	 $all_err = '';
	 { 
	     local $/ = \1;
	     while (my $line = <$CHLD_ERR>){
		 $all_err .= $line;
		 print STDERR $line unless($main::quiet);
	     }
	 }
	 waitpid $pid, 0;
	 
	 if ($? != 0 && ! ignorable($all_err, $command)){
	    die "ERROR: RapSearch failed\n";
	 }
      }
   }
   else {
      die "you must give Widget::rapsearch a command to run!\n";
   }
}

sub ignorable{
    my $err = shift;
    my $command = shift;

    my $ignorable = 0;
    #if($err =~ /There are no valid contexts/){
    #	$ignorable = 1;
    #}
    #else{
    #	$ignorable = 0;
    #}

    return $ignorable;
}

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
    my $report       = shift;
    my $params       = shift;

    #these should be in report?
    my $q_seq_length = $params->{query_length};
    my $t_seq_length = 'NA';
    my $depth        = $params->{depth};
    my $t_index      = GI::build_fasta_index($params->{db});
    
    my $hitType = 'Bio::Search::Hit::PhatHit::rapsearch';
    my $hspType = 'Bio::Search::HSP::PhatHSP::rapsearch';

    my $count = 0;
    my $hit;
    my $keepers = [];
    my $fh = new FileHandle();
    $fh->open($report);
    while (my $first = <$fh>){
	chomp($first);
	next if(!$first); #just in case
	
	#Full Alignment Example
	#chr10 vs ARATH bits=23.1 log(E-value)=0.98 identity=46.1% aln-len=13 mismatch=7 gap-openings=0 nFrame=5
	#Query: 79399 CHRCRGFGHVIRD 79361
	#             C  C+G GH+I++
	#Sbjct:   221 CQTCKGTGHIIKE 233
	#
	#.........
	
	#parse IDs (Example -->)
	#chr10 vs ARATH
	my ($q_name, $t_name, $s_line);
	if($first =~ /^(\S+)\s+vs\s+(\S+)\s+(.*)/){
	    $q_name = $1;
	    $t_name = $2;
	    $s_line = $3;
	}
	else{
	    die "ERROR: File is not a rapSearch report or is corrupt: $report\n";
	}
	
	#make new hit if necessary
	if(!$hit || $hit->name ne $t_name){
	    #parse/cluster/filter at same time
	    if($hit){
		push(@$keepers, @{keepers([$hit], $params)});
		$keepers = cluster::shadow_cluster($depth, $keepers)
		    if($depth && @$keepers >= $depth*20);
	    }

	    $t_seq_length = $t_index->get_length_by_id($t_name);
	    if(! defined($t_seq_length) && !$t_index->{__REBUILT}){
		my $db = $params->{db};
		$db =~ s/\.\d+$//;
		my $list = [grep {/\.\d+$/} <$db.*>];
		if(@$list > 1){
		    $t_index = GI::build_fasta_index($list);
		    $t_seq_length = $t_index->get_length_by_id($t_name);
		}
		$t_index->{__REBUILT} = 1; #set flag so I don't rebuild again
	    }

	    $hit = $hitType->new('-name'        => $t_name,
				 '-description' => '',
				 '-algorithm'   => 'rapsearch',
				 '-length'      => $t_seq_length);
	    $count++; #one new hit
	}
	
	#parse stats (Example -->)
	#bits=23.1 log(E-value)=0.98 identity=46.1% aln-len=13 mismatch=7 gap-openings=0 nFrame=5
	my %stats = map {split(/=/, )} split(/\s+/, $s_line);
	
	#parse alignments (Example -->)
	#Query: 79399 CHRCRGFGHVIRD 79361
	#             C  C+G GH+I++
	#Sbjct:   221 CQTCKGTGHIIKE 233
	my $query  = <$fh>; #query
	my $match  = <$fh>; #match
	my $target = <$fh>; #target
	my $empty  = <$fh>; #empty line
	
	my ($q_start, $q_str, $q_end) = $query  =~ /Query\:[\t\s]+(\d+)\s+(\S+)\s+(\d+)/;
	my ($t_start, $t_str, $t_end) = $target =~ /Sbjct\:[\t\s]+(\d+)\s+(\S+)\s+(\d+)/;
	my $q_strand = ($q_start < $q_end) ? 1 : -1;
	my $t_strand = 0; #proteins are not stranded 
	
	my $pos = index($query, $q_str, 5);
	my $m_str = substr($match, $pos, length($q_str));
	
	#add additional stats
	$stats{query_gaps} = scalar($q_str =~ tr/\-/\-/);
	$stats{hit_gaps}   = scalar($t_str =~ tr/\-/\-/);
	$stats{identical}  = scalar($m_str =~ tr/a-zA-Z/a-zA-Z/);
	$stats{conserved}  = $stats{identical} + scalar($m_str =~ tr/\+/\+/);	

	#make HSP arguments
	my @args;
	push(@args, '-query_name');
	push(@args, $q_name);

	push(@args, '-hit_name');
	push(@args, $t_name);

	push(@args, '-algorithm');
	push(@args, 'rapsearch');	

	push(@args, '-evalue');
	push(@args, 10**$stats{'log(E-value)'});

	push(@args, '-pvalue');
	push(@args, $stats{'log(E-value)'}); #hack

	push(@args, '-bits');
	push(@args, $stats{bits});

	push(@args, '-score');
	push(@args, 'NA');
	
	push(@args, '-identical');
	push(@args, $stats{identical});

	push(@args, '-conserved');
	push(@args, $stats{conserved});
	
	push(@args, '-query_start');
	push(@args, $q_start);

	push(@args, '-query_end');
	push(@args, $q_end);	

	push(@args, '-query_length');
	push(@args, $stats{'aln-len'} - $stats{query_gaps});
	
	push(@args, '-query_gaps');
	push(@args, $stats{query_gaps});

	push(@args, '-query_seq');
	push(@args, $q_str);

	push(@args, '-hit_start');
	push(@args, $t_start);

	push(@args, '-hit_end');
	push(@args, $t_end);
	
	push(@args, '-hit_length');
	push(@args, $stats{'aln-len'} - $stats{hit_gaps});

	push(@args, '-hit_gaps');
	push(@args, $stats{hit_gaps});
	
	push(@args, '-hit_seq');
	push(@args, $t_str);

	push(@args, '-hsp_length');
	push(@args, $stats{'aln-len'});
	
	push(@args, '-homology_seq');
	push(@args, $m_str);
	
	my $hsp = $hspType->new(@args);
	$hsp->queryName($q_name);
	$hsp->query_name($q_name);
	$hit->queryLength($q_seq_length);

	#-------------------------------------------------
	# setting strand because bioperl is all messed up!
	#-------------------------------------------------
	$hsp->{_strand_hack}->{query} = $q_strand;
	$hsp->{_strand_hack}->{hit}   = $t_strand;
	$hsp->{_indentical_hack}      = $stats{identical};
	
	$hit->add_hsp($hsp);
    }
    $fh->close();

    #add final hit
    if($hit){
	push(@$keepers, @{keepers([$hit], $params)});
    }

    #final clustering
    if($depth){
	$keepers = cluster::shadow_cluster($depth, $keepers);
	$keepers = GI::flatten($keepers);
    }

    #report
    my $end     = @$keepers;
    my $deleted = $count - $end;
    print STDERR "deleted:$deleted hits\n" unless $main::quiet;

    return $keepers;
}
#-------------------------------------------------------------------------------
sub keepers {
    my $hits   = shift;
    my $params = shift; 

    my @keepers;

    my $split_hit   = $params->{split_hit} || 0;
    my $hsp_bit_min = $params->{hsp_bit_min} || 0;

    #iterate hits
    foreach my $hit (@$hits){
	my $q_length = $hit->queryLength;
	my $cutoff = ($params->{is_last}) ? $q_length+1 : $q_length - $split_hit;
	my $scutoff = ($params->{is_first}) ? 0 : 1 + $split_hit;
	
	#iterate HSPs
	my @hsps;
	foreach my $hsp ($hit->hsps) {
	    #clear a little memory (will I need this later?)
	    $hsp->cigar_string; #holds alignment
	    $hsp->{QUERY_SEQ} = '';
	    $hsp->{HIT_SEQ} = '';
	    $hsp->{HOMOLOGY_SEQ} = '';

	    push(@hsps, $hsp) if($hsp->bits > $hsp_bit_min);
	}
	$hit->hsps(\@hsps);

	#split hits
	my $hits = PhatHit_utils::split_hit_by_strand([$hit]);
	$hits = PhatHit_utils::split_hit_on_intron($hits, $split_hit) if($split_hit);
	
	#filter split hits
	foreach my $h (@{$hits}) {
	    next unless($h->hsps());
	    
	    #just in case fix strand for messed up hit
	    $h = PhatHit_utils::copy($h, 'both') if ($h->strand('hit') < 0);
	    
	    my $s = $h->start('query');
	    my $e = $h->end('query');
	    if($split_hit && ($s <=  $scutoff || $e >= $cutoff)){
		$h->start(); #force create of stored start and end
		$h->strand(); #force creation
		$h->getLengths(); #force creation of stored lengths
		if(!_keep_filt($h, $params)){ #memory optimization
		    $h->{_hsps} = [];
		    $h->{_remove} = 1;
		}
		push(@keepers, $h);
		next;
	    }

	    #filtering
	    next unless(_keep_filt($h, $params));	    
	}
    }
    
    return \@keepers;
}
#-----------------------------------------------------------------------------
sub _keep_filt {
    my $h = shift;
    my $params = shift || {};

    #filtering
    #next unless PhatHit_utils::is_contigous($h);
    if(defined($params->{significance})){
	my $sig_thr = $params->{significance};
	my $significance = $h->significance();
	$significance = "1".$significance if  $significance =~ /^e/;
	$significance = 0                 if  $significance =~ /0\.$/;
	return 0 if($significance > $sig_thr);
    }
    if(defined($params->{percov})){
	my $pcov = $params->{percov} || 0;
	return 0 if($h->pAh < $pcov/2);
    }
    if(defined($params->{percid})){
	my $pid = $params->{percid} || 0;
	return 0 if($h->hsp('best')->frac_identical() < $pid &&
		    $h->frac_identical < $pid);
    }

    return 1;
}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::rapsearch::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


