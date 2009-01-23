#!/usr/bin/perl -w 
use strict;
use Getopt::Std;

##### Initialize Threshhold  ####
my @thresh = ();
use vars qw($opt_h $opt_c $opt_e $opt_o $opt_a $opt_t $opt_l);

push @thresh, 0.5;
push @thresh, 0.5;
push @thresh, 0.5;
push @thresh, 0;
push @thresh, 0;
push @thresh, 75;

getopts("hc:e:o:a:t:");
my $usage = "maker2zff.pl directory name [options]

directory  - the location of the gff files to use for the hmm
name       - a name for the outfile, will produce name.dna and name.ann

OPTIONS
For determining which genes are High Confidence for Retraining, there are 6 criteria.
-c fraction  The fraction of splice sites confirmed by an EST alignment, default 0.5
-e fraction  The fraction of exons that overlap an EST alignment, default 0.5
-o fraction  The faction of exons that overlap any evidence (EST or Protein), default 0.5
-a fraction  The fraction of splice sites confirmed by an ab-initio SNAP prediction, default 0
-t fraction  The fraction of exons the overlap an ab-initio SNAP prediction, default 0
-l number    The min length of the protein sequence produced by the mRNA
";

if ($opt_c) {$thresh[0] = $opt_c}
if ($opt_e) {$thresh[1] = $opt_e}
if ($opt_o) {$thresh[2] = $opt_o}
if ($opt_a) {$thresh[3] = $opt_a}
if ($opt_t) {$thresh[4] = $opt_t}
if ($opt_l) {$thresh[5] = $opt_l}

my %id2name = ();
my %status = ();
my %parent = ();
my %exons = ();
my %hc = ();
my %seq = ();

my $dir = shift @ARGV;
my $outfile = shift @ARGV;

if ($opt_h) {die $usage;}
if ((not $dir) || (not $outfile)) {
    die $usage;
}

#### Get all of the GFF file in the directory  ####
opendir DIR, "$dir" or die $!;
my @files = grep /^.+\.gff$/, readdir DIR;
close DIR;

#### Go through the GFF file and determine which genes are HC  ####

foreach my $file (@files) {
    open GFF, "<$dir/$file" or die $!;
    while (my $line = <GFF>) {
	chomp($line);
	if ($line =~ m/^\s*#/) {
	    next;
        } elsif ($line =~ m/^\s*$/) {
	    next;
        } elsif ($line =~ m/^>(\S+)/) {
	    my $header = $1;
	    $seq{$header} = "";
	    while (<GFF>) {
		$seq{$header} = join("", $seq{$header}, $_)
	    }
        } else {
	    my ($seqid, $source, $tag, $start, $end, $score, $strand, $phase, $annot) = split(/\t/, $line);
	    my %annotation = split(/[;=]/, $annot);
	    if ($tag eq "mRNA" && $source eq "maker") {
		my $id = $annotation{'ID'};
		my $lname = $annotation{'Name'};
		my $parent = $annotation{'Parent'};
		my ($name, $qi) = split(/\s+/, $lname);
		my $ishc = is_hc($qi);
		if ($ishc == 1 ) {
		    $hc{$id} = 1;
		}
	    }elsif ($tag eq "CDS" && $source eq "maker") {
		my $parent = $annotation{'Parent'};
		if ($strand eq "-") {
		    my $temp = $start;
		    $start = $end;
		    $end = $temp;
		}
		push @{$exons{$seqid}}, [$start, $end, $parent]
	    }
	}
    }
}

#### Print out the exon locations of the HC genes ####

open(ZFF, ">$outfile\.ann") or die;
open DNA, ">$outfile\.dna" or die $!;

foreach my $seqid (sort {$a cmp $b} keys %exons) {
    my @exons = @{$exons{$seqid}};
    print ZFF ">",$seqid,"\n";
    print DNA ">",$seqid,"\n",$seq{$seqid},"\n";
    foreach my $e (@exons) {
        my ($start, $end, $pname) = @{$e};
	if ($hc{$pname}) {
	    print ZFF join(" ", "Exon ",$start, $end, $pname),"\n";
	}
    }
}

################## Subroutines ###################

sub is_hc {
    my $qi = shift @_;
    my @q = split(/\|/, $qi);
    my @qual = (@q[1..5],$q[8]);
    my $hc = 1;
    foreach my $i (0..5) {
	if ($thresh[$i] < $qual[$i]) {
	    $hc = 0;
	}
    }
    return $hc;
}
