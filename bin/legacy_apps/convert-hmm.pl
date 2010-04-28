#!/usr/bin/perl
use strict; use warnings;
use Getopt::Std;
our ($opt_u, $opt_d);
getopts('u:d:');

my $UPSTREAM   = 500;
my $DOWNSTREAM = 500;

die "
usage: $0 [options] <snap hmm>
options
  -u <int>   [default $UPSTREAM]
  -d <int>   [default $DOWNSTREAM]
" unless @ARGV == 1;

$DOWNSTREAM = $opt_d if $opt_d;
$UPSTREAM   = $opt_u if $opt_u;

MAINLOOP: while (<>) {
	if (/^zoeHMM/) {
		my ($z, $n, $s, $t, $d, $m) = split;
		$s++; # states increase by 1
		$d++; # durations increase by 1
		$m++; # models increase by 1
		print "zoeHMM PRO-SNAP:$n $s $t $d $m\n\n";
	} elsif (/<STATES>/) {
		print; print "\n";
		my $blank = <>;
		while (<>) {
			last if /^\s/;
			if (/Inter/) {
				print "UTR5\t1.0\t0.0\t0\t0\tgeometric\n";
				print "UTR3\t0.0\t1.0\t0\t0\tgeometric\n";
			} elsif (/Intron/) {
				print "Intron\t0.0\t0.0\t0\t0\tgeometric\n";
			} else {
				print;
			}
		}
	} elsif (/<STATE_TRANSITIONS>/) {
		print "\n"; print; print "\n";
		my $blank = <>;
		while (<>) {
			last if /^\s/;
			my ($s1, $s2, $p) = split;
			if    ($s1 eq 'Inter') {$s1 = 'UTR5'}
			elsif ($s2 eq 'Inter') {$s2 = 'UTR3'}
			print join("\t", $s1, $s2, $p), "\n";
		}
	} elsif (/<PHASE_PREFERENCES>/) {
		print "\n"; print; print "\n";
		my $blank = <>;
		while (<>) {
			last if /^\s/;
			print;
		}
	} elsif (/<STATE_DURATIONS>/) {
		print "\n"; print; print "\n";
		my $blank = <>;
		while (<>) {
			last MAINLOOP if /<SEQUENCE_MODELS>/;
			if (/^Inter (\d+)/) {
				# skip over Inter duration
				while (<>) {last if $_ !~ /\S/}
				# replace with UTR durations
				print "UTR5 1\n\tGEOMETRIC 0 -1\n\t\t$UPSTREAM\n\n";
				print "UTR3 1\n\tGEOMETRIC 0 -1\n\t\t$DOWNSTREAM\n\n";
			} else {
				print;
			}
		}
	}
}

print "<SEQUENCE_MODELS>\n";
while (<>) {
	if (/^Inter(.+)/) {
		my (@utr5, @utr3);
		push @utr5, "UTR5$1\n";
		push @utr3, "UTR3$1\n";
		while (<>) {
			last if $_ !~ /\S/;
			push @utr5, $_;
			push @utr3, $_;
		}
		print @utr5;
		print "\n";
		print @utr3;
	}
	print;
}
