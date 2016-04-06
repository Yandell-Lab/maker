#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Location::Split;
use Bio::Location::Simple;

#This script was adapted from Jason Stajich's zff2augustus_gbk.pl script
#-->https://github.com/hyphaltip/genome-scripts/blob/master/gene_prediction/

#usage
my $usage = "
Usage:
     zff2genbank.pl <export.ann> <export.dna>

Description:
     Converts SNAP ZFF format to GenBank format for Augustus training.
     Works with files produced by the fathom script during SNAP training
     (i.e. export.ann and export.dna).

";

#get command line options
my $zff = shift @ARGV;
my $dna = shift @ARGV;

#print usage
if(!defined($zff) || !defined($dna)){
    print $usage;
    exit(0);
}

#validate options
my $err;
foreach my $file ($zff, $dna){
    $err .= "ERROR: File does not exist: $file\n" if(!-f $file);
}
die $err if($err);

#open files
my $dbh = Bio::DB::Fasta->new($dna);
my $out = Bio::SeqIO->new(-format => 'genbank',-fh => \*STDOUT);
open(my $fh => $zff) or die $!;

#process ZFF
my $seq;
my @location;
while(1) {
    my $line = <$fh>;

    #sequence line
    if(!defined($line) || $line =~ /^>(\S+)/){
	# process the previous sequence first
	if($seq){
	    my $loc = (@location == 1) ? shift(@location) : Bio::Location::Split->new();
	    $loc->add_sub_Location($_) for(sort {$a->start <=> $b->start} @location);
	    my $gene = Bio::SeqFeature::Generic->new(-primary_tag => 'CDS',
						     -location => $loc);
	    $seq->add_SeqFeature($gene);
	    $out->write_seq($seq);
	    @location = ();
	}
	last if(!defined($line));

	my $seqid = $1;
	my $seqstr = $dbh->seq($seqid);
	die("ERROR: cannot find $seqid in the input file $dna\n") if(!defined $seqstr);
	
	$seq = Bio::Seq->new(-seq => $seqstr, -id => $seqid);
	$seq->add_SeqFeature(Bio::SeqFeature::Generic->new(-primary_tag => 'source',
							   -start => 1,
							   -end   => $dbh->length($seqid)));
	next;
    }

    #feature location line
    my @F = split(/\t/, $line);
    die "ERROR: Input does not appear to be ZFF\n" if(@F != 4 && @F != 9);
    
    my ($B, $E, $s) = ($F[1] >= $F[2]) ? ($F[1], $F[2], 1) : ($F[2], $F[1], -1);
    push @location, Bio::Location::Simple->new(-start  => $B,
					       -end    => $E,
					       -strand => $s);
}
close($fh);

#finished
exit(0);
