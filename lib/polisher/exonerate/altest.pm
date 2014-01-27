#------------------------------------------------------------------------
#----                     polisher::exonerate::altest                   ---- 
#------------------------------------------------------------------------
package polisher::exonerate::altest;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Widget::exonerate::cdna2genome;
use polisher;
use polisher::exonerate;
use Exporter;
use PostData;
use FileHandle;
use PostData;
use Exporter;
use Fasta;
use FastaFile;
use URI::Escape;
use PhatHit_utils;

@ISA = qw(
	polisher::exonerate
	polisher
       );

#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub polish {
	my $g_file     = shift;
	my $e_file     = shift;
	my $o_file     = shift;
	my $g_len      = shift;
	my $e_len      = shift;	
	my $the_void   = shift;
	my $offset     = shift || 0;
	my $exe        = shift;
	my $percent    = shift;
	my $min_intron = shift;
	my $max_intron = shift;
	my $matrix     = shift;

	$e_len = polisher::prep($e_file) if(!$e_len);
	$g_len = polisher::prep($g_file) if(!$g_len);

	my $hits = a_exonerate($g_file, 
			       $e_file, 
			       $o_file,
			       $the_void, 
			       $e_len, 
			       $g_len,
			       $exe,
			       $percent,
			       $min_intron,
			       $max_intron,
			       $matrix,
			      );


	return undef unless $hits->[0];

	PhatHit_utils::add_offset($hits, $offset);

	return $hits;
}
#------------------------------------------------------------------------
sub a_exonerate {
        my $g_file   = shift;
        my $e_file   = shift;
	my $o_file     = shift;
        my $the_void = shift;
        my $e_len    = shift;
        my $g_len    = shift;
	my $exe      = shift;
	my $percent  = shift;
	my $min_intron = shift;
	my $max_intron = shift;
	my $matrix   = shift;


        runExonerate($g_file, $e_file, $o_file, $exe, $percent, $min_intron, $max_intron, $matrix);

        return Widget::exonerate::cdna2genome::parse($o_file, $e_len, $g_len);
}
#-----------------------------------------------------------------------------
sub runExonerate {
        my $t_file  = shift;
        my $q_file  = shift;
        my $o_file  = shift;
	my $exe     = shift;
	my $percent = shift;
	my $min_intron = shift || 20;
	my $max_intron = shift || 200000;
	my $matrix  = shift;


        my $command  = "$exe  -q $q_file -t $t_file -Q dna -T dna";
	   $command .= " --model cdna2genome ";
	   $command .= " --minintron $min_intron --maxintron $max_intron --showcigar";
	   $command .= " --percent $percent";
	if ($matrix) {
	    $command .= " --dnasubmat $matrix";	
	}
	$command .= " > $o_file";
	
        my $w = new Widget::exonerate::cdna2genome();
	
        if (-e $o_file){
	    print STDERR "re reading exonerate report.\n" unless($main::quiet);
	    print STDERR "$o_file\n"unless($main::quiet);
        }
        else {
	    print STDERR "running  cdna2genome search.\n"unless($main::quiet);
	    $w->run($command);
        }
	
}
#-----------------------------------------------------------------------------
1;


