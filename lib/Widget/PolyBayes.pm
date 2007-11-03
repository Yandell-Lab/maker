#------------------------------------------------------------------------
#----                        Widget::PolyBayes                       ---- 
#------------------------------------------------------------------------
package Widget::PolyBayes;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Bio::DB::Fasta;
use Snp;
@ISA = qw(
	Widget
       );

BEGIN {
	$ENV{'POLYBAYES_LIB'} = '/users/myandell/polybayes/PolyBayes_3_0/lib';
}
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
sub runCrossMatch {
        my $wDir = shift;

        my $ests   = "$wDir/edit_dir/CLUSTER.members.fasta";
        my $anchor = "$wDir/edit_dir/CLUSTER.anchor.fasta";

        my $out    = "$wDir/edit_dir/CLUSTER.cm";
        my $cm   = '/usr/local/bdgp/bin/cross_match';
        my $args = "-discrep_lists -tags -masklevel 5 $ests $anchor";

        system("$cm $args > $out");


}
#-----------------------------------------------------------------------------
sub runAnchoredAlignment {
        my $wDir = shift;

        my $polyB = "/users/myandell/polybayes/PolyBayes_3_0/bin/polybayes.pl";

        my $p = "$wDir/edit_dir/CLUSTER";

        my $aDna  = "$wDir/edit_dir/CLUSTER.anchor.fasta";

        my $mDna    = "$wDir/edit_dir/CLUSTER.members.fasta";
        my $mQual   = "$wDir/edit_dir/CLUSTER.members.fasta.qual";
        my $mBpos   = "$wDir/edit_dir/CLUSTER.members.fasta.bpos";

        my $cm      = "$wDir/edit_dir/CLUSTER.cm";
        my $aceOut  = "$wDir/edit_dir/CLUSTER.polybayes.alignment.ace";

        my $rOut    = "$wDir/edit_dir/CLUSTER.polybayes.alignment.out";
        my $phdPath = "$wDir/phd_dir";


        my $command  = "$polyB -writeAlignment -project $p -analysis polybayes ";
           $command .= "-cluster $p  -inputFormat map -anchorDna $aDna ";
           $command .= "-memberDna $mDna  -memberBaseQuality $mQual ";
           $command .= "-memberBasePosition $mBpos -crossMatch $cm ";
           $command .= "-anchorBaseQualityDefault 40 -nofilterParalogs ";
           $command .= "-paralogFilterMinimumBaseQuality 12 ";
           $command .= "-noscreenSnps -writeAce -aceOut $aceOut ";
           $command .= "-writeReport -reportOut $rOut -debug -monitor ";
           $command .= "-writePhdFiles -phdFilePathOut  $phdPath";

        system("$command");

}
#-----------------------------------------------------------------------------
sub runScreenSNPS {

        my $wDir = shift;

        my $polyB = "/users/myandell/polybayes/PolyBayes_3_0/bin/polybayes.pl";

        my $p = "$wDir/edit_dir/CLUSTER";

        my $aDna  = "$wDir/edit_dir/CLUSTER.anchor.fasta";

        my $mDna    = "$wDir/edit_dir/CLUSTER.members.fasta";
        my $mQual   = "$wDir/edit_dir/CLUSTER.members.fasta.qual";
        my $mBpos   = "$wDir/edit_dir/CLUSTER.members.fasta.bpos";

        my $cm      = "$wDir/edit_dir/CLUSTER.cm";
        my $aceOut  = "$wDir/edit_dir/CLUSTER.polybayes.snp.ace";

        my $rOut    = "$wDir/edit_dir/CLUSTER.polybayes.snp.out";
        my $phdPath = "$wDir/phd_dir";

        my $aceIn = "$wDir/edit_dir/CLUSTER.polybayes.alignment.ace";
        my $phdIn = "$wDir/phd_dir";


        my $command  = "$polyB -writeAlignment -project $p -analysis polybayes ";
           $command .= "-cluster $p  -inputFormat ace -aceIn $aceIn ";
           $command .= " -readPhdFiles -phdFilePathIn $phdIn ";
           $command .= "-anchorBaseQualityDefault 40 -filterParalogs ";
           $command .= "-thresholdNative 0.65 -paralogFilterMinimumBaseQuality 12 ";
           $command .= "-screenSnps -considerAnchor -considerTemplateConsensus ";
           $command .= "-priorPoly 0.003 -thresholdSnp 0.5 -writeAce -aceOut $aceOut ";
           $command .= "-writeReport -reportOut $rOut -debug -monitor -writePhdFiles ";
           $command .= "-phdFilePathOut $phdPath";

        system("$command");

}
#-------------------------------------------------------------------------------
sub parse {
	my $self = shift;
        my $file = shift;

	$self->file($file);

        my $fh = new FileHandle();
           $fh->open("$file");

	my %alignedLengths;
        while (my $line = <$fh>){


		if    ($line =~ /ALIGNED_LENGTH/){
			my ($name)   = $line =~ /.*PRIMARY_SEQUENCE\:\s+(\S+).*/;
			#my ($length) = $line =~ /.*ALIGNED_LENGTH\:\s+(\S+).*/;
			my ($length) = $line =~ /.*PADDED_LENGTH\:\s+(\S+).*/;
			$alignedLengths{$name} = $length;

		}
                next unless $line =~ /SNP_ID/;
		my $snp = new Snp(); 
		my %snp;
                my ($depth) = $line =~ /SEQUENCE_DEPTH:\s+(\d+)/;
                my ($var)   = $line =~ /VAR:\s+(\w+\-*)/;
                my ($aVar, $rVar)  = $var =~ /([\w\-])([\w\-])/;
                my ($aName) = $line=~ /ALIGNMENT:\s+(\S+)\s+/;
                my ($a, $c, $g, $t, $d) =
                $line =~ /SEQUENCE_ALLELES:\s+a=(\d+),c=(\d+),g=(\d+),t=(\d+)-=(\d+)/;

		my ($aSnpPos) = $line =~ /.*UNPADDED_POS:\s+(\d+).*/; 
		#my ($aLength) = $line =~ /.*UNPADDED_LENGTH:\s+(\d+).*/;
		my ($aLength) = $line =~ /.*UNPADDED_LENGTH:\s+(\d+).*/;

		$alignedLengths{$aName} = $aLength;

		my ($pVar) = $line =~ /P_VAR\:\s+(\d?\.?\d+e?\-?\d+)\s+/;
		my ($p1)   = $line =~ /P_1\:\s+(\d?\.?\d+e?\-?\d+)\s+/;
                my ($p2)   = $line =~ /P_2\:\s+(\d?\.?\d+e?\-?\d+)\s+/;
                my ($p3)   = $line =~ /P_3\:\s+(\d?\.?\d+e?\-?\d+)\s+/;
                my ($p4)   = $line =~ /P_4\:\s+(\d?\.?\d+e?\-?\d+)\s+/;

		$snp->aLength($aLength);
		$snp->aSnpPos($aSnpPos);

		$snp->pVar($pVar);
		$snp->p1($p1);
                $snp->p2($p2);
                $snp->p3($p3);
                $snp->p4($p4);
		
		$snp->depth($depth);
		$snp->var($var);
		$snp->aName($aName);
		$snp->count('a', $a);
		$snp->count('t', $t);		
                $snp->count('g', $g);
                $snp->count('c', $c);
                $snp->count('-', $d);

                my ($data) = $line =~ /COLUMN:\s+(.*)$/;


		my $column = parseColumn($data);
		my $col = new Column();
		   $col->addColumn($column);
		   $col->aName($snp->aName);
		   $col->addAlignedLengths(\%alignedLengths);

		$snp->addColumn($col, 'original');

		my $unique = Column::makeUnique($col);
		   $unique->aName($snp->aName);
		$snp->addColumn($unique, 'unique');

		$self->addSnp($snp);

        }

        $fh->close();
}
#-------------------------------------------------------------------------------
sub hasSnp {
	my $self = shift;

	if (defined(@{$self->{snps}})){
		return 1;
	}
	else {
		return 0;
	}
}
#-------------------------------------------------------------------------------
sub addSnp {
	my $self = shift;
	my $snp  = shift;

	 push(@{$self->{snps}}, $snp) if defined($snp);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parseColumn {
        my $data = shift;

        my @fields = split("TEMPLATE", $data);

        shift(@fields);

        my %column;
        foreach my $f (@fields){
                my @stuff = split(",", $f);

                my $sName       = $stuff[4];
                my $qV          = $stuff[8];
                my $sB          = $stuff[7];
		my $orientation = $stuff[10];
		my $snpPos      = $stuff[9];
		my $status      = $stuff[11];
		   $status      =~ s/\|//;
		   $status      =~ s/\;SEQUENCE//;
		$column{$sName}{sO} =  $orientation;
		$column{$sName}{sP} =  $snpPos;
		$column{$sName}{aS} =  $status;
		$column{$sName}{qV} =  $qV;
                $column{$sName}{qV} =  $qV;
                $column{$sName}{sB} =  $sB;
        }


        return \%column;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::RepeatMasker::AutoLoader called for: ",
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


