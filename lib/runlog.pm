#-------------------------------------------------------------------------------
#------                           runlog                               ---------
#-------------------------------------------------------------------------------
package runlog;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use File::Find::Rule;

@ISA = qw();
$VERSION = 0.1;

#===make list of internal variables to log
my @ctl_to_log = ('est',
		  'protein',
		  'repeat_protein',
		  'rmlib',
		  'rm_gff',
		  'predictor',
		  'snaphmm',
		  'augustus_species',
		  'model_org',
		  'max_dna_len',
		  'split_hit',
		  'snap_flank',
		  'te_remove',
		  'single_exon',
		  'alt_peptide',
		  'percov_blastn',
		  'percid_blastn',
		  'eval_blastn',
		  'bit_blastn',
		  'percov_blastx',
		  'percid_blastx',
		  'eval_blastx',
		  'bit_blastx',
		  'e_perc_cov',
		  'ep_score_limit',
		  'en_score_limit',
		 );
   
my @opt_to_log = ('R',
		  'GFF',
		  'PREDS'
		 );

my %SEEN;

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $self = {};
    my $class = shift;

    bless ($self, $class);

    $self->{CTL_OPTIONS} = shift;
    $self->{OPT} = shift;
    $self->{the_void} = shift;
    $self->{file_name} = shift || "run.log";
    $self->{file} = $self->{the_void} . "/" . $self->{file_name};

    $self->{file} =~ s/\/\//\//g;

    $self->_load_old_log();
    $self->_clean_files();
    $self->_write_new_log();

    return $self;
}
#-------------------------------------------------------------------------------
sub _load_old_log {
   my $self = shift;

   $self->{die_count} = 0;

   my $log_file = $self->{file};
   my %logged_vals;

   if (-e $log_file){
      open (IN, "< $log_file");      
      while( defined (my $line = <IN>)){
	 chomp $line;
      
	 my ($type, $key, $value) = split ("\t", $line);
	 $logged_vals{$type}{$key} = defined($value) ? $value : '';
	 
	 $self->{die_count} = $value if($type eq 'DIED' && $key eq 'COUNT');
      }
      close(IN);
   }

   $self->{old_log} = \%logged_vals;
}
#-------------------------------------------------------------------------------
sub _clean_files{
   my $self = shift;

   my $the_void = $self->{the_void};
   my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};
   my %OPT = %{$self->{OPT}};
   
   #get gff file name
   my $log_file = $self->{file};
   my $gff_file = $the_void;
   $gff_file =~ s/theVoid\.([^\/]+$)/$1.gff/;

   #===id files types that can not be re-used
   my $continue_flag = 1; #signal whether to continue with this seq
   
   my %logged_vals = %{$self->{old_log}}; #values from existing log
   my %rm_key; #key to what type of old files to remove
   my @files; #list of files to remove
   my @dirs; #list of directories to remove

   if (-e $log_file) {
      if (exists $logged_vals{DIED}) {
	 if($OPT{f} && ! $SEEN{$log_file}){
	    $self->{die_count} = 0; #reset the die count
	    $continue_flag = 1; #re-run
	    $rm_key{force}++;
	    $SEEN{$log_file}++;
	 }
	 elsif ($OPT{retry}){
	    $continue_flag = 2;	#rerun died
	    $continue_flag = -3 if($self->{die_count} > $OPT{retry}); #only let die up to count
	    $rm_key{retry}++ if ($continue_flag == 2);
	 }
	 else{
	    $continue_flag = -1; #skip died
	 }
      }
      elsif ($OPT{f}) {
	 $rm_key{force}++;
	 $continue_flag = 1;	#run/re-run
      }
      else {
	 $continue_flag = 0 if (-e $gff_file); #don't re-run finished
      
	 #CHECK COMMAND LINE OPTIONS FOR CHANGES
	 while (my $key = each %{$logged_vals{OPT}}) {
	    my $log_val = $logged_vals{OPT}{$key};
	    my $opt_val = '';
	    $opt_val = $OPT{$key} unless (! defined $OPT{$key});
	 
	    #if previous log options are not the same as current the options
	    if ($log_val ne $opt_val) {
	       print STDERR "WARNING: Change in command line flag \'$key\' since last run\n";
	       $continue_flag = 1; #re-run because opts changed
	    
	       if (-e $gff_file) {
		  $rm_key{gff}++; #always rebuild gff when some option has changed
	       }       
	    
	       if ($key eq 'R' ||
		   $key eq 'GFF'
		  ) {
		  $rm_key{all}++;
	       }
	    }
	 }
	 
	 #CHECK CONTROL FILE OPTIONS FOR CHANGES
	 while (my $key = each %{$logged_vals{CTL_OPTIONS}}) {
	    #needed only for MPI
	    if (defined $CTL_OPTIONS{"old_$key"}){
	       $CTL_OPTIONS{$key} = $CTL_OPTIONS{"old_$key"};
	    }

	    my $log_val = $logged_vals{CTL_OPTIONS}{$key};
	    my $ctl_val = '';
	    $ctl_val = $CTL_OPTIONS{$key} if(defined $CTL_OPTIONS{$key});

	    #if previous log options are not the same as current control file options
	    if ($log_val ne $ctl_val) {
	       print STDERR "WARNING: Control file option \'$key\' has changed\n".
	       "Old:$log_val\tNew:$ctl_val\n\n";
	    
	       $continue_flag = 1; #re-run because ctlopts changed
	    
	       if (-e $gff_file) {
		  $rm_key{gff}++; #always rebuild gff when some option has changed
	       }
	    
	       if ($key eq 'max_dna_len') {
		  $rm_key{all}++;
	       }
	    
	       if ( ($key eq 'repeat_protein' ||
		     $key eq 'rm_lib' ||
		     $key eq 'model_org' ||
		     $key eq 'te_remove') &&
		    ! $OPT{R} &&
		    ! $OPT{GFF}
		  ) {
		  $rm_key{all}++;
	       }
	    
	       if ( $key eq 'rm_gff' &&  
		    $OPT{GFF} &&
		    ! $OPT{R}
		  ) {
		  $rm_key{all}++;
	       }
	    
	       if ($key eq 'snaphmm' ||
		   $key eq 'snap_flank' ||
		   $key eq 'single_exon'
		  ) {
		  $rm_key{snap}++;
	       }
	    
	       if ($key eq 'augustus_species' ||
		   $key eq 'snap_flank' ||
		   $key eq 'single_exon'
		  ) {
		  $rm_key{augustus}++;
	       }
	    
	       if ($key eq 'split_hit' ||
		   $key eq'ep_score_limit' ||
		   $key eq'en_score_limit'
		  ) {
		  $rm_key{exonerate}++;
	       }
	    
	       if ($key eq 'protein' ||
		   $key eq 'alt_peptide' ||
		   $key eq 'eval_blastx'
		  ) {
		  $rm_key{blastx}++;
		  $rm_key{exonerate}++;
	       }
	    
	       if ($key eq 'est' ||
		   $key eq 'eval_blastn'
		  ) {
		  $rm_key{blastn}++;
		  $rm_key{exonerate}++;
	       }
	    }
	 }
      
	 #CHECK FOR FILES THAT DID NOT FINISH
	 while (my $key = each %{$logged_vals{STARTED}}) {
	    if (! exists $logged_vals{FINISHED}{$key}) {
	       print STDERR "WARNING: The file $key from the last run did not\n".
	       "finish and must be erased\n";
	       push(@files, $key);
	    }
	 }
      }
   
      #===Remove files that can no longer be re-used
   
      #-print file type specific warnings and delete files
      if (exists $rm_key{retry}) {
	 print STDERR "WARNING: Old data must be removed before re-running this sequence\n";
      
	 #delete everything in the void
	 File::Path::rmtree($the_void);
	 File::Path::mkpath($the_void);
	 unlink($gff_file) if(-e $gff_file);
      }
      elsif (exists $rm_key{force}) {
	 print STDERR "WARNING: All old files will be erased before continuing\n";
      
	 #delete everything in the void
	 File::Path::rmtree($the_void);
	 File::Path::mkpath($the_void);
	 unlink($gff_file) if(-e $gff_file);
      }
      elsif (exists $rm_key{all}) {
	 print STDERR "WARNING: Changes in control files make re-use of all old data impossible\n".
	 "All old files will be erased before continuing\n";
      
	 #delete everything in the void
	 File::Path::rmtree($the_void);
	 File::Path::mkpath($the_void);
	 unlink($gff_file) if(-e $gff_file);
      }
      else {
	 if (exists $rm_key{snap}) {
	    print STDERR "WARNING: Changes in control files make re-use of old Snap data impossible\n".
	    "Old Snap files will be erased before continuing\n";
	 
	    my @f = File::Find::Rule->file()->name( '*snap' )->in($the_void);
	    push (@files, @f);
	 }
	 if (exists $rm_key{augustus} && defined($logged_vals{CTL_OPTIONS}{augustus})) {
	    print STDERR "WARNING: Changes in control files make re-use of old Augustus data impossible\n".
	    "Old Augustus files will be erased before continuing\n";
	 
	    my @f = File::Find::Rule->file()->name( '*augustus' )->in($the_void);
	    push (@files, @f);
	 }
	 if (exists $rm_key{blastn}) {
	    print STDERR "WARNING: Changes in control files make re-use of old Blastn data impossible\n".
	    "Old Blastn files will be erased before continuing\n";
	 
	    my @f = File::Find::Rule->file()->name( '*blastn*' )->in($the_void);
	    push (@files, @f);
	 
	    my @d = File::Find::Rule->directory()->name( '*blastn*' )->in($the_void);
	    push (@dirs, @d);
	 }
	 if (exists $rm_key{blastx}) {
	    print STDERR "WARNING: Changes in control files make re-use of old Blastx data impossible\n".
	    "Old Blastx files will be erased before continuing\n";
	 
	    my @f = File::Find::Rule->file()->name( '*blastx' )->in($the_void);
	    push (@files, @f);
	 
	    my @d = File::Find::Rule->directory()->name( '*blastx*' )->in($the_void);
	    push (@dirs, @d);
	 }
	 if (exists $rm_key{exonerate}) {
	    print STDERR "WARNING: Changes in control files make re-use of old Exonerate data impossible\n".
	    "Old Exonerate files will be erased before continuing\n";
	 
	    my @f = File::Find::Rule->file()->name( '*exonerate' )->in($the_void);
	    push (@files, @f);
	 }
	 if (exists $rm_key{gff}) {
	    print STDERR "WARNING: The gff file $gff_file must now be removed.\n";
	    push (@files, $gff_file);
	 }
      
	 #delete files in the void
	 foreach my $file (@files) {
	    unlink($file);
	 }
      
	 #delete directories in the void
	 foreach my $dir (@dirs) {
	    File::Path::rmtree($dir);
	 }
      }
   }

   $self->{continue_flag} = $continue_flag;
}
#-------------------------------------------------------------------------------
#==build log file for current settings
sub _write_new_log {
   my $self = shift;

   my $log_file = $self->{file};
   
   my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};
   my %OPT = %{$self->{OPT}};

   return if ($self->{continue_flag} <= 0);

   open (LOG, "> $log_file");

   #log command line options
   foreach my $key (@opt_to_log) {
      my $opt_val = '';
      $opt_val = $OPT{$key} unless (! defined $OPT{$key});
      
      print LOG "OPT\t$key\t$opt_val\n";
   }
   
   #log control file options
   foreach my $key (@ctl_to_log) {
      #needed only for MPI
      if (defined $CTL_OPTIONS{"old_$key"}){
	 $CTL_OPTIONS{$key} = $CTL_OPTIONS{"old_$key"};
      }

      my $ctl_val = '';
      $ctl_val = $CTL_OPTIONS{$key} if(defined $CTL_OPTIONS{$key});
				  
      print LOG "CTL_OPTIONS\t$key\t$ctl_val\n";
   }
   close(LOG);
}
#-------------------------------------------------------------------------------
sub add_entry {
   my $self = shift;

   my $type  = shift;
   my $key   = shift;
   my $value = shift;

   my $log_file = $self->{file};

   open(LOG, ">> $log_file");
   print LOG "$type\t$key\t$value\n";
   close(LOG);
}
#-------------------------------------------------------------------------------
sub get_die_count {
   my $self = shift;

   return $self->{die_count};
}
#-------------------------------------------------------------------------------
sub get_continue_flag {
   my $self = shift;

   #-----------------
   #  0 => already ran don't re-run
   #  1 => run or re-run with current settings
   #  2 => re-run previously died
   # -1 => don't re-run, it died
   # -2 => not finished but skip because $OPT{retry} is in force
   # -3 => skip when $OPT{retry} is in force because die_count is too high
   #-----------------

   return $self->{continue_flag};
}
#-------------------------------------------------------------------------------
1;
