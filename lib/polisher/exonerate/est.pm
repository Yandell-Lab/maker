#------------------------------------------------------------------------
#----                     polisher::exonerate::est                   ---- 
#------------------------------------------------------------------------
package polisher::exonerate::est;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Widget::exonerate::est2genome;
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

@ISA = qw(
	polisher::exonerate
	polisher
       );

my $OPT_F; #GLOBAL VALUE
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub polish {
	my $g_file     = shift;
	my $e_file     = shift;
	my $the_void   = shift;
	my $offset     = shift || 0;
	my $ext        = shift || 0;
	my $exe        = shift;
	my $percent    = shift;
	my $matrix     = shift;
	$OPT_F         = shift;

	my ($g_id, $e_id, $e_len, $g_len) = 
	polisher::prep($g_file, $e_file);

	my $hits = e_exonerate($g_file, 
			       $e_file, 
			       $the_void, 
			       $g_id, 	
			       $e_id, 
			       $e_len, 
			       $g_len,
			       $ext,
			       $exe,
			       $percent,
			       $matrix,
			      );


	return undef unless $hits->[0];

	foreach my $f (@{$hits}){
		polisher::add_offset($offset, $f);
	}
	return $hits;
}
#------------------------------------------------------------------------
sub e_exonerate {
        my $g_file   = shift;
        my $e_file   = shift;
        my $the_void = shift;
        my $g_id     = shift;
        my $e_id     = shift;
        my $e_len    = shift;
        my $g_len    = shift;
	my $ext      = shift;
	my $exe      = shift;
	my $percent  = shift; 
	my $matrix   = shift;

	my $safe_g_id = uri_escape($g_id,  #build a safe name for file names from the sequence identifier
				   '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				   );

	my $safe_e_id = uri_escape($e_id,  #build a safe name for file names from the sequence identifier
				   '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
				   );

        my $o_file    = "$the_void/$safe_g_id\.$safe_e_id\.$ext\.est_exonerate";

        runExonerate($g_file, $e_file, $o_file, $exe, $percent, $matrix);

        return Widget::exonerate::est2genome::parse($o_file, $e_len, $g_len);

}
#-----------------------------------------------------------------------------
sub runExonerate {
        my $t_file  = shift;
        my $q_file  = shift;
        my $o_file  = shift;
	my $exe     = shift;
	my $percent = shift;
	my $matrix  = shift;

        my $command  = "$exe  -q $q_file -t $t_file -Q dna -T dna";
	   $command .= " --model est2genome ";
	   $command .= " --minintron 20 --showcigar";
	   $command .= " --percent $percent";
	if ($matrix) {
	    $command .= " --dnasubmat $matrix";	
	}
	$command .= " > $o_file";
	
        my $w = new Widget::exonerate::est2genome();
	
        if (-e $o_file && ! $OPT_F){
	    print STDERR "re reading exonerate report.\n";
	    print STDERR "$o_file\n";
        }
        else {
	    print STDERR "running  est2genome search.\n";
	    $w->run($command);
        }
	
}
#-----------------------------------------------------------------------------
1;


