#!/usr/bin/perl -w
#make_better_wbannot.pl

use Bio::DB::GFF;
use Bio::Tools::GFF;
use Bio::FeatureIO::gff;

my $dir = $ARGV[0];
if (not $dir) {
    die "add_utr_gff.pl <directory>\n
where directory is the location of your gff files.\n"
}

opendir DIR, "$dir" or die $!;
my @files = grep/.+\.gff$/, readdir DIR;
close DIR;

foreach my $file (@files) {
    chomp($file);
    my $outfile = "$dir/$file";
    $outfile =~ s/\.gff/\.wutr\.gff3/;
    open NEW, ">$outfile" or die $!;
    print NEW "##gff-version 3\n";
    my %exon = ();
    my %parent = ();
    my %id2name = ();
    my %cds = ();
    my %strand = ();
    
    my $parser = new Bio::Tools::GFF(-gff_version => 3,
				     -file        => "$dir/$file");
    
    while (my $feature = $parser->next_feature()) {
	my $tag = $feature->primary_tag;
	my $start = $feature->start;
	my $end = $feature->end;
	if ($feature->source_tag eq "maker") {
	    if ($tag eq "mRNA") {
		my @id = $feature->get_tag_values('ID');
		my $id = $id[0];
		my @name = $feature->get_tag_values('Name');
		my $name = $name[0];
		$id2name{$id} = $name;
		$strand{$name} = $feature->strand;
	    }elsif ($tag eq "CDS") {
		my @parents = $feature->get_tag_values('Parent');
		foreach my $p (@parents) {
		    my $pid = $p;
		    my $pname = $id2name{$pid};
		    $strand{$pname} = $feature->strand;
		    push @{$cds{$pname}}, ($start, $end);
		}
	    }elsif ($tag eq "exon") {
		my @parents = $feature->get_tag_values('Parent');
		foreach my $p (@parents) {
		    my $pid = $p;
		    my $pname = $id2name{$pid};
		    push @{$exons{$pname}}, ($start, $end);
		}
	    }
	}
    }
    my %utr = ();
 F1:foreach my $mrna (keys %exons) {
	my @exons = sort {$a <=> $b} @{$exons{$mrna}};
	
	next F1 unless $cds{$mrna};
	my @cds = sort {$a <=> $b} @{$cds{$mrna}};
	
	$cds_start = $cds[0];
	$cds_end = $cds[$#cds];
	$num_exons = scalar(@exons)/2;
	
	my $i = 0;
	my $j = 1;
	foreach (1..$num_exons) {
	    my $exon_start = $exons[$i];
	    my $exon_end = $exons[$j];
	    if (not $exon_start) {
		print "no exon start\n";
	    }
	    if ($exon_start < $cds_start) {
		if ($exon_end < $cds_start) {
		    push @{$utr{$mrna}}, ["lower_utr",$exon_start, $exon_end];
		}else {
		    push @{$utr{$mrna}}, ["lower_utr",$exon_start, $cds_start - 1];
		}
	    }elsif ($exon_end > $cds_end) {
		if ($exon_start > $cds_end) {
		    push @{$utr{$mrna}}, ["upper_utr",$exon_start, $exon_end];
		}else {
		    push @{$utr{$mrna}}, ["upper_utr",$cds_end+1, $exon_end];
		}
	    }
	    $i += 2;
	    $j += 2;
	}
    }
    open GFF, "<$dir/$file" or die $!;
    my $line = <GFF>;
 W1:while ($line = <GFF>) {
	chomp($line);
	print NEW $line,"\n";
	my @features = split(/\t/, $line);
	my $seqid = $features[0];
	my $source = $features[1];
	next W1 unless $source;
	if ($source eq "maker") {
	    my $tag = $features[2];
	    my $start = $features[3];
	    my $end = $features[4];
	    my $strand = $features[6];
	    if ($tag eq "mRNA") {
		%annotations = split(/[=;]/, $features[8]);
		my $name = (split(/\s+/, $annotations{'Name'}))[0];
		$strand = $strand{$name};
		if ($utr{$name}) {
		    my $utr_id = 1;
		    foreach $utr (@{$utr{$name}}) {
			($type, $start, $end) = @{$utr};
			my $utr_start = $start;
			my $utr_end =  $end;
			if ($strand == 1 && $type eq "lower_utr") {
			    $utrtag = "five_prime_UTR";
			    $gstrand = "+";
			}elsif ($strand == 1 && $type eq "upper_utr") {
			    $utrtag = "three_prime_UTR";
			    $gstrand = "+";
			}elsif ($strand == -1 &&  $type eq "lower_utr") {
			    $utrtag = "three_prime_UTR";
			    $gstrand = "-";
			}else {
			    $utrtag = "five_prime_UTR";
			    $gstrand = "-";
			}
			my $id = "ID=$name\:utr\:$utr_id";
			my $parent = "Parent=$name";
			$annot = join(";", $id, $parent);
			print NEW join("\t",$seqid,"maker",$utrtag,$utr_start,$utr_end,".", $gstrand,".",$annot),"\n";
			$utr_id ++;
		    }
		}
	    }
	}
    }
}
