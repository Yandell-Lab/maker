#!/usr/bin/perl -w 
use strict;
use Getopt::Std;

##### Initialize Threshhold  ####
my @thresh = ();
my $thrAED = 0.5;
use vars qw($opt_h $opt_c $opt_e $opt_o $opt_a $opt_t $opt_l $opt_x);

push @thresh, 0.5;
push @thresh, 0.5;
push @thresh, 0.5;
push @thresh, 0;
push @thresh, 0;
push @thresh, 75;


getopts("hc:e:o:a:t:x");
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
-x number    Max AED to allow 0.5 is default
";

if ($opt_c) {$thresh[0] = $opt_c}
if ($opt_e) {$thresh[1] = $opt_e}
if ($opt_o) {$thresh[2] = $opt_o}
if ($opt_a) {$thresh[3] = $opt_a}
if ($opt_t) {$thresh[4] = $opt_t}
if ($opt_l) {$thresh[5] = $opt_l}
if ($opt_x) {$thrAED = $opt_x}

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
        if ($line =~ m/\#\#FASTA/) {
	    my $header;
	    while (my $d = <GFF>) {
		if($d =~ /^>(\S+)/){
		    $header = $1;
		    $seq{$header} = "";
		}
		else{
		    $seq{$header} .= $d;
		}
	    }
        }
	elsif($line =~ /^\s*\#|^\n$|^\s*$/){
	    next;
	}
	else {
	    my ($seqid, $source, $tag, $start, $end, $score, $strand, $phase, $annot) = split(/\t/, $line);
	    my %annotation = split(/[;=]/, $annot);
	    if ($tag eq "mRNA" && $source eq "maker") {
		my $id = $annotation{'ID'};
		my $lname = $annotation{'Name'};
		my $parent = $annotation{'Parent'};
		my ($name, $qi) = split(/\s+/, $lname);
		if(! $qi){
		    ($qi) = $line =~ /_QI\=([^\;\n]+)/;
		}
		my ($AED) = $line =~ /_AED\=([^\;\n]+)/;
		my $ishc = is_hc($qi, $AED);
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
    my $AED = shift @_;

    my @q = split(/\|/, $qi);
    my @qual = (@q[1..5],$q[8]);
    my $hc = 1;
    foreach my $i (0..5) {
	if ($qual[$i] <= $thresh[$i]) {
	    $hc = 0;
	}
    }

    $hc = 0 if($AED > $thrAED);

    return $hc;
}
