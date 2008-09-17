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

@ISA = qw(
	polisher::exonerate
	polisher
       );

my $OPT_F; #GLOBAL VALUE
my $LOG;
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub polish {
	my $g_file     = shift;
	my $p_file     = shift;
	my $the_void   = shift;
	my $offset     = shift || 0;
	my $ext        = shift || 0;
	my $exe        = shift;
	my $percent    = shift;
	my $matrix     = shift;
	$OPT_F         = shift;
	$LOG           = shift;

	my ($g_id, $p_id, $p_len, $g_len) = 
	polisher::prep($g_file, $p_file);

	my $hits = p_exonerate($g_file, 
			    $p_file, 
		            $the_void, 
		            $g_id, 	
		            $p_id, 
		            $p_len, 
		            $g_len,
			    $ext,
			    $exe,
			    $percent,
			    $matrix
		);


	return undef unless $hits->[0];

	foreach my $f (@{$hits}){
		polisher::add_offset($offset, $f);
	}
	return $hits;
}
#------------------------------------------------------------------------
sub p_exonerate {
        my $g_file   = shift;
        my $p_file   = shift;
        my $the_void = shift;
        my $g_id     = shift;
        my $p_id     = shift;
        my $p_len    = shift;
        my $g_len    = shift;
	my $ext      = shift;
	my $exe      = shift;
	my $percent  = shift;
	my $matrix   = shift;

        my $safe_g_id = uri_escape($g_id,  #build a safe name for file names from the sequence identifier
                                   '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
                                   );

        my $safe_p_id = uri_escape($p_id,  #build a safe name for file names from the sequence identifier
                                   '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
                                   );

        my $o_file    = "$the_void/$safe_g_id\.$safe_p_id\.$ext\.p_exonerate";
	
	$LOG->add_entry("STARTED", $o_file, "") if(defined $LOG);
        runExonerate($g_file, $p_file, $o_file, $exe, $percent, $matrix);
	$LOG->add_entry("FINISHED", $o_file, "") if(defined $LOG);

        return Widget::exonerate::protein2genome::parse($o_file, $p_len, $g_len);
	
    }
#-----------------------------------------------------------------------------
sub runExonerate {
        my $t_file = shift;
        my $q_file = shift;
        my $o_file = shift;
	my $exe    = shift;
	my $percent = shift;
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
	
        if (-e $o_file && ! $OPT_F){
	    print STDERR "re reading exonerate report.\n" unless($main::quiet);
	    print STDERR "$o_file\n" unless($main::quiet);
        }
        else {
	    print STDERR "running  exonerate search.\n" unless($main::quiet);
	    $w->run($command);
        }
	
	
    }
#-----------------------------------------------------------------------------
1;


