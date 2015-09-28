#------------------------------------------------------------------------
#----                        Widget::snoscan                            ---- 
#------------------------------------------------------------------------
package Widget::snoscan;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Fasta;
use FastaFile;
use Iterator::Fasta;
use Bio::Search::Hit::PhatHit::snoscan;
use Bio::Search::HSP::PhatHSP::snoscan;
use PhatHit_utils;
use IPC::Open3;
use Symbol;
use FastaSeq;
use Error qw(:try);
use Error::Simple;
use GI;
@ISA = qw(
	Widget
       );

my $OPT_F; # global varible for maker clean up
my $LOG; # global varible for maker runlog

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
		{		
		    local $/ = \1;
		    while (my $line = <$CHLD_ERR>){
			print STDERR $line unless($main::quiet);
		    }
		}
		waitpid $pid, 0;
		die "ERROR: Snoscan failed\n" if($? != 0);
	}
	else {
		die "you must give Widget::snoscan a command to run!\n";
	}
}
##-------------------------------------------------------------------------------
##------------------------------ FUNCTIONS --------------------------------------
##-------------------------------------------------------------------------------
sub parse {
    shift @_ if($_[0] =~ /^Widget\:\:/); #ignore
    my $report = shift;
    my $params = shift;
    my $q_file = shift;
    
    my $def;
    my $q_seq;

    if(ref($q_file) eq 'FastaChunk'){ #object not scalar
	$q_seq = $q_file->seq_w_flank;
	$def = $q_file->def_w_flank;
    }
    elsif($q_file =~ /^>/){
	$def     = Fasta::getDef(\$q_file);
	$q_seq   = Fasta::getSeqRef(\$q_file);
    }
    else{
	my $index = GI::build_fasta_index($q_file);
	my ($id) = $index->get_all_ids;
	$def = ">$id";
	$q_seq   = $index->get_Seq_by_id($id);

    }
    
    my $q_name  = Fasta::def2SeqID($def);
    my $fh = new FileHandle();
    $fh->open($report);

    my %d;
    while (my $line = <$fh>){
        chomp($line);
	next unless $line =~ /^>>/;

        my @stuff = split(/\s+/, $line);
	my $score = $stuff[2];
        my $key = $stuff[3];
        my $match = $stuff[5];
        my ($b, $e) = $key =~ /\((\S+?)-(\S+?)\)/;
        my $strand = ($b < $e) ? 1 : -1;

	if(!defined($d{$key})){
	    $d{$key}{'b'} = $b;
	    $d{$key}{'e'} = $e;
	    $d{$key}{'strand'} = $strand;
	    $d{$key}{'score'} = $score;
	    $d{$key}{'best'}  = $match;
	}
	elsif($score > $d{$key}{'score'}){
	    $d{$key}{'score'} = $score;
	    $d{$key}{'best'}  = $match;
	}
	push(@{$d{$key}{'match'}}, $match);
    }
    $fh->close();
    my %g;
    my $i = 0;
    foreach my $key (keys %d){
	my $best = $d{$key}{'best'};
	my $id = "snoscan_$i\__$best";
	$g{$id}[0]{strand}   = $d{$key}{'strand'};
	$g{$id}[0]{type}     = 'exon';
	$g{$id}[0]{b}        = $d{$key}{'b'}; 
	$g{$id}[0]{e}        = $d{$key}{'e'};
	$g{$id}[0]{score}    = $d{$key}{'score'};
	$g{$id}[0]{match}    = join(',', @{$d{$key}{'match'}});
	$i++;
    }

    my $phat_hits = load_phat_hits($q_name, $q_seq, \%g, $params);
    return $phat_hits;
}
##-------------------------------------------------------------------------------
sub get_exon_seq {
    my $exon  = shift;
    my $q_seq = shift;
    
    my $e_b = $exon->{b};
    my $e_e = $exon->{e};
    
    my $length = abs($e_e - $e_b) + 1;
    
    ($e_b, $e_e) = ($e_e, $e_b) if ($e_b > $e_e);
    
    my $e_seq = substr_o($q_seq, $e_b - 1, $length);
    
    $e_seq = Fasta::revComp($e_seq) if $exon->{strand} == -1;
    
    return $e_seq;
}
##-------------------------------------------------------------------------------
sub total_score {
    my $g = shift;
    
    my $total = 0;
    foreach my $exon (@{$g}){
	$total += $exon->{score};
    }
    return $total;
}
##-------------------------------------------------------------------------------
sub load_phat_hits {
    my $q_name  = shift;
    my $q_seq   = shift;
    my $g       = shift;
    my $params  = shift;
    
    my $q_len = length_o($q_seq);
    my @keepers;
    foreach my $gene (keys %{$g}){	
	#build hit
	my %hsps;
	my $i = 0;
	my $f = new Bio::Search::Hit::PhatHit::snoscan('-name'         => $gene,
						       '-description'  => 'NA',
						       '-algorithm'    => 'snoscan',
						       '-length'       => $q_len,
						       '-score'        => '',
						       );
	
	$f->queryLength($q_len);
	#build hsps
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
	    push(@args, 'snoscan');
	    
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
	    
	    my $hsp = new Bio::Search::HSP::PhatHSP::snoscan(@args);
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

    my $final = \@keepers;

    return $final;
}
##-----------------------------------------------------------------------------
sub AUTOLOAD {
    my ($self, $arg) = @_;
    
    my $caller = caller();
    use vars qw($AUTOLOAD);
    my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
    $call =~/DESTROY/ && return;
    
    print STDERR "Widget::snoscan::AutoLoader called for: ",
    "\$self->$call","()\n";
    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
    
    if (defined($arg)){
	$self->{$call} = $arg;
    }
    else {
	return $self->{$call};
    }
}
##------------------------------------------------------------------------
1;


