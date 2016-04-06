#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(ceil);
use Getopt::Long;
use File::Temp;
use Cwd;
use File::Which;
use File::Spec;
use File::Copy;

my $usage = "
Usage:
     train_augustus.pl <genes.gb> <species>

Options:
     cpus      <INT>    Number of CPUs to use for training
     max       <INT>    Maximum models to use for optimization
     help|?             Prints this usage statement

";

#get options
my $cpus = 1;
my $max;
GetOptions ("cpus=i" => \$cpus, 
	    "max=i"  => \$max, 
	    "help|?" => sub { print $usage; exit(0); })
    or die("ERROR: Failed command line arguments\n");

my $gb = shift;
my $species = shift;

#validate options
if(!$gb || !-f $gb){
    print $usage;
    die "ERROR: Failure to supply GenBank input file\n";
}
if(!defined($species)){
    print $usage;
    die "ERROR: Failure to supply species\n";
}
$gb = Cwd::abs_path($gb);
$species =~ s/[\s\t]+/_/g;

my $count = `grep -c "LOCUS" $gb`;
$max = $count if(!$max || $max > $count);
die "ERROR: Training ZFF has no annotations\n" if(!$count);

#get scripts
my $aug_exe = File::Which::which('augustus');
die "ERROR: Could not find Augustus executable\n" if(!$aug_exe);

my $cwd = Cwd::cwd();
(my $base = Cwd::abs_path($aug_exe)) =~ s/\/bin\/augustus$//;
$base = $cwd if(!-d $base);

my $acp = $ENV{AUGUSTUS_CONFIG_PATH};
$acp = "$base/config" if(!$acp || !-d $acp);
die "ERROR: Could not find AUGUSTUS_CONFIG_PATH\n" if(!-d $acp);

my $scripts = "$base/scripts";
my ($nsp_exe) = grep {$_ && -f $_} map {("$scripts/$_", File::Which::which($_))} ('new_species.pl');
my ($trn_exe) = grep {$_ && -f $_} map {("$scripts/$_", File::Which::which($_))} ('etraining');
my ($opt_exe) = grep {$_ && -f $_} map {("$scripts/$_", File::Which::which($_))} ('optimize_augustus.pl');
my ($ran_exe) = grep {$_ && -f $_} map {("$scripts/$_", File::Which::which($_))} ('randomSplit.pl');
die "ERROR: Could not find Augustus accessory scripts\n" if(!$nsp_exe || !$trn_exe || !$opt_exe || !$ran_exe);

#prepare
my $tdir = File::Temp::tempdir(CLEANUP => 1, TEMPLATE => 'augtrainXXXXX', DIR => File::Spec->tmpdir);
my $tacp = "$tdir/config";
my $tgb  = "$tdir/genes.gb";
File::Copy::copy($gb, $tgb);
mkdir("$tacp");
mkdir("$tacp/species");
symlink("$acp/species/generic", "$tacp/species/generic");
symlink("$acp/extrinsic", "$tacp/extrinsic");
symlink("$acp/profile", "$tacp/profile");
symlink("$acp/model", "$tacp/model");
$ENV{AUGUSTUS_CONFIG_PATH} = $tacp;

#train
print STDERR "STATUS: Training ...\n";
my $k = ($cpus < $max) ? $cpus : $max;
$k = 8 if($k < 8);
chdir($tdir); #work in temp directory
do_system("$ran_exe $tgb $max")  or die "ERROR: Failed to create optimization set\n";
do_system("$nsp_exe --species=$species") or die "ERROR: Failed to create new Augustus species\n";
do_system("$trn_exe --species=$species $gb") or die "ERROR: Failed to initialize training for Augustus species\n";
do_system("$opt_exe --species=$species --onlytrain=$tgb.train --kfold=$k --cpus=$cpus $tgb.test") or die "ERROR: Failed to optimize Augustus training\n";
do_system("$trn_exe --species=$species $gb") or die "ERROR: Failed final Augustus training\n";
chdir($cwd); #change back

#evaluate
#print STDERR "STATUS: Evaluating final parameter performance ...\n";
#do_system("$ran_exe $tgb 100")  or die "ERROR: Failed to create final evaluation set\n";
#do_system("$aug_exe --species=$species $tgb.test | grep -A 22 Evaluation")  or die "ERROR: Failed to evaluate final training performance\n";

#result
do_system("tar -C $tacp/species -zcvf $cwd/$species.tgz $species") or die "ERROR: Could not package Augustus parameters\n";
print STDERR "STATUS: Final species parameters --> $cwd/$species.tgz\n";

exit(0);

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
sub do_system {
    system(@_);

    return 1 if($? == 0);
}
