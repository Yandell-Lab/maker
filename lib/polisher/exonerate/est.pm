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
my $LOG;
#------------------------------------------------------------------------
#--------------------------- FUNCTIONS ----------------------------------
#------------------------------------------------------------------------
sub polish {
	my $g_file     = shift;
	my $e_file     = shift;
	my $the_void   = shift;
	my $offset     = shift || 0;
	my $exe        = shift;
	my $percent    = shift;
	my $matrix     = shift;
	$OPT_F         = shift;
	$LOG           = shift;

	my ($e_len, $g_len) = polisher::prep($g_file, $e_file);

	my $hits = e_exonerate($g_file, 
			       $e_file, 
			       $the_void, 
			       $e_len, 
			       $g_len,
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
        my $e_len    = shift;
        my $g_len    = shift;
	my $exe      = shift;
	my $percent  = shift; 
	my $matrix   = shift;

        my ($g_name) = $g_file =~ /([^\/]+)$/;
        $g_name =~ s/\.fasta$//;
        my ($e_name) = $e_file =~ /([^\/]+)$/;
        $e_name =~ s/\.fasta$//;

        my $o_file    = "$the_void/$g_name.$e_name.est_exonerate";

	$LOG->add_entry("STARTED", $o_file, "") if(defined $LOG);   
        runExonerate($g_file, $e_file, $o_file, $exe, $percent, $matrix);
	$LOG->add_entry("FINISHED", $o_file, "") if(defined $LOG);   

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
	    print STDERR "re reading exonerate report.\n" unless($main::quiet);
	    print STDERR "$o_file\n"unless($main::quiet);
        }
        else {
	    print STDERR "running  est2genome search.\n"unless($main::quiet);
	    $w->run($command);
        }
	
}
#-----------------------------------------------------------------------------
1;


