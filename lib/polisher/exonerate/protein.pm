#------------------------------------------------------------------------
#----                   polisher::exonerate::protein                 ---- 
#------------------------------------------------------------------------
package polisher::exonerate::protein;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Widget::exonerate::protein2genome;
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
	my $p_file     = shift;
	my $o_file     = shift;
	my $g_len      = shift;
        my $p_len      = shift;
	my $the_void   = shift;
	my $offset     = shift || 0;
	my $exe        = shift;
	my $percent    = shift;
	my $min_intron = shift;
	my $max_intron = shift;
	my $matrix     = shift;

	$p_len = polisher::prep($p_file) if(!$p_len);
        $g_len = polisher::prep($g_file) if(!$g_len);

	my $hits = p_exonerate($g_file, 
			       $p_file, 
			       $o_file,
			       $the_void, 
			       $p_len, 
			       $g_len,
			       $exe,
			       $percent,
			       $min_intron,
			       $max_intron,
			       $matrix
			       );


	return undef unless $hits->[0];

	PhatHit_utils::add_offset($hits, $offset);

	return $hits;
}
#------------------------------------------------------------------------
sub p_exonerate {
        my $g_file   = shift;
        my $p_file   = shift;
	my $o_file   = shift;
        my $the_void = shift;
        my $p_len    = shift;
        my $g_len    = shift;
	my $exe      = shift;
	my $percent  = shift;
	my $min_intron = shift;
	my $max_intron = shift;
	my $matrix   = shift;

        runExonerate($g_file, $p_file, $o_file, $exe, $percent, $min_intron, $max_intron, $matrix);

        return Widget::exonerate::protein2genome::parse($o_file, $p_len, $g_len);
	
    }
#-----------------------------------------------------------------------------
sub runExonerate {
        my $t_file = shift;
        my $q_file = shift;
        my $o_file = shift;
	my $exe    = shift;
	my $percent = shift;
	my $min_intron = shift || 20;
	my $max_intron = shift || 200000;
	my $matrix = shift;

        my $command  = "$exe -q $q_file -t $t_file -Q protein -T dna ";
	$command .= "-m protein2genome --softmasktarget ";
	$command .= " --percent $percent";
	if ($matrix) {
	    $command .= " --proteinsubmat $matrix";
	}
	$command .= " --showcigar ";
	$command .= " > $o_file";
	
        my $w = new Widget::exonerate::protein2genome();
	
        if (-e $o_file){
	    print STDERR "re reading exonerate report.\n" unless($main::quiet);
	    print STDERR "$o_file\n" unless($main::quiet);
        }
        else {
	    print STDERR "running  exonerate search.\n" unless($main::quiet);
	    $w->run($command, $o_file);
        }
    }
#-----------------------------------------------------------------------------
1;


