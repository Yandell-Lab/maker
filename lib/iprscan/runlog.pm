#-------------------------------------------------------------------------------
#------                       iprscan::runlog                          ---------
#-------------------------------------------------------------------------------
package iprscan::runlog;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Fasta;
use File::NFSLock;

@ISA = qw();
$VERSION = 0.1;

#===make list of internal variables to log
my @ctl_to_log = ('infile',
		  'outfile',
		  'appl',
		  'nocrc',
		  'seqtype',
		  'trtable',
		  'trlen',
		  'goterms',
		  'iprlookup',
		  'format'
		 );

my %SEEN;

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $self = {};
    my $class = shift;
    my @args = @_;

    bless ($self, $class);    

    if($self->_initialize(@args)){
       $self->_load_old_log();
       $self->_clean_files();
       $self->_write_new_log();
    }

    $self->report_status();

    return $self;
}
#-------------------------------------------------------------------------------
sub _initialize {
   my $self = shift;
   $self->{CTL_OPTIONS} = shift;
   $self->{params} = shift;
   $self->{file_name} = shift || "run.log";

   print STDERR "\n\n\n--Next Contig--\n\n" unless($main::quiet);

   $self->{CWD} = $self->{CTL_OPTIONS}->{CWD} || Cwd::cwd();

   #lock must be persitent in object or it is detroyed outside of block
   unless($self->{LOCK} = new File::NFSLock($self->{file_name}, 'NB', undef, 40)){
       $self->{continue_flag} = -3;
       return 0;
   }

   return 1;
}
#-------------------------------------------------------------------------------
sub _load_old_log {
   my $self = shift;

   $self->{die_count} = 0;

   my $log_file = $self->{file_name};
   my %logged_vals;

   if (-e $log_file){#load log file if available
      print STDERR "Processing run.log file...\n" unless($main::quiet);
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

    my $the_void = $self->{params}->{the_void};
    my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};
    
    #get ipr file name
    my $log_file = $self->{file_name};
    my $ipr_file = $the_void;
    my $out_base = $the_void;
    $ipr_file =~ s/theVoid\.([^\/]+)$/$1.ipr/;
    $out_base =~ s/theVoid\.[^\/]+$//;
    
    #===id files types that can not be re-used
    my $continue_flag = 1; #signal whether to continue with this seq
    
    my %logged_vals = %{$self->{old_log}}; #values from existing log
    my %rm_key; #key to what type of old files to remove
    my @files; #list of files to remove
    my @dirs; #list of directories to remove
    
    if (-e $log_file) {
	if (exists $logged_vals{DIED}) {
	    if($CTL_OPTIONS{force} && ! $SEEN{$log_file}){
		$self->{die_count} = 0; #reset the die count
		$continue_flag = 1; #re-run
		$rm_key{force}++;
		$SEEN{$log_file}++;
	    }
            elsif($CTL_OPTIONS{always_try} && ! $SEEN{$log_file}){
                $self->{die_count} = 0; #reset the die count
                $continue_flag = 1; #re-run
                $SEEN{$log_file}++;
	    }
	    else{
		$continue_flag = ($CTL_OPTIONS{clean_try}) ? 2 : 3;	#rerun died
		$continue_flag = -1 if($self->{die_count} > $CTL_OPTIONS{retry}); #only let die up to count
		$rm_key{retry}++ if ($continue_flag == 2);
		$SEEN{$log_file}++;
	    }
	}
	elsif ($CTL_OPTIONS{force} && ! $SEEN{$log_file}) {
	    $rm_key{force}++;
	    $continue_flag = 1;	#run/re-run
	    $SEEN{$log_file}++;
	}
	elsif ($CTL_OPTIONS{again} && ! $SEEN{$log_file}){
	    $continue_flag = 1; #run/re-run
	    $rm_key{ipr}++;
	    $SEEN{$log_file}++;
	}
	else {
	    $continue_flag = 0 if (-e $ipr_file); #don't re-run finished
	    $SEEN{$log_file}++;
	}
	
	if($continue_flag >= 0 || $continue_flag == -1){
	    #CHECK CONTROL FILE OPTIONS FOR CHANGES
	    my $cwd = ($self->{CWD}) ? $self->{CWD} : Cwd::cwd();
	    while (my $key = each %{$logged_vals{CTL_OPTIONS}}) {
		my $log_val = '';
		if(defined $logged_vals{CTL_OPTIONS}{$key}){	       
		    $log_val = $logged_vals{CTL_OPTIONS}{$key};
		    if($key eq 'appl'){
			#don't care about order
			my @set = split(',', $logged_vals{CTL_OPTIONS}{appl});
			@set = sort @set;
			$log_val = join(',', @set);
		    }
		}
		
		my $ctl_val = '';
		if(defined $CTL_OPTIONS{$key}){
		    $ctl_val = $CTL_OPTIONS{$key};
		    $ctl_val =~ s/^$cwd\/*//;
		    if($key eq 'appl'){
			#don't care about order
			my @set = sort @{$CTL_OPTIONS{appl}};
			$ctl_val = join(',', @set);
		    }
		}
		
		#if previous log options are not the same as current control file options
		if ($log_val ne $ctl_val) {
		    
		    print STDERR "MAKER WARNING: Control file option \'$key\' has changed\n".
			"Old:$log_val\tNew:$ctl_val\n\n" unless($main::quiet);
		    
		    $continue_flag = 1; #re-run because ctlopts changed
		    
		    $rm_key{ipr}++; #always rebuild IPR when some option has changed
		    
		    #certain keys that affect everything
		    if ($key eq 'nocrc' ||
			$key eq 'seqtype' ||
			$key eq 'trlen' ||
			$key eq 'trtable' ||
			$key eq 'goterms' ||
			$key eq 'iprlookup' ||
			$key eq 'format'
			) {
			$rm_key{all}++;
		    }
		}
	    }
	    
	    #CHECK FOR FILES THAT DID NOT FINISH
	    while (my $key = each %{$logged_vals{STARTED}}) {
		if (! exists $logged_vals{FINISHED}{$key}) {
		    print STDERR "MAKER WARNING: The file $key\n".
			"did not finish on the last run and must be erased\n" unless($main::quiet);
		    push(@files, $key);
		    $key =~ /([^\/]+)$/;
		    my $rm_f_name = $1;
		    $rm_f_name =~ s/\.fasta$//;
		    my @d = <$the_void/*$rm_f_name*>;
		    foreach my $d (@d){
			push (@dirs, $d) if (-d $d);
		    }
		}
	    }
	}

	#===Remove files that can no longer be re-used
	
	#-print file type specific warnings and delete files
	if (exists $rm_key{force}) {
	    print STDERR "WARNING: All old files will be erased before continuing\n" unless($main::quiet);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    	    
	    #remove all other files
	    unlink($ipr_file) if(-e $ipr_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	}
	elsif (exists $rm_key{retry}) {
	    print STDERR "WARNING: Old data must be removed before re-running this sequence\n" unless($main::quiet);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    unlink($ipr_file) if(-e $ipr_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	}
	elsif (exists $rm_key{all}) {
	    print STDERR "WARNING: Changes in control files make re-use of all old data impossible\n".
		"All old files will be erased before continuing\n" unless($main::quiet);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    unlink($ipr_file) if(-e $ipr_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	}
	elsif (exists $rm_key{ipr}) {
	    print STDERR "WARNING: Any preexisting IPR files for this contig must now be removed.\n" unless($main::quiet);
	    push (@files, $ipr_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	}
	    
	#delete files in the void
	foreach my $file (@files) {
	    unlink($file) if(-e $file);
	}
	
	#delete directories in the void
	foreach my $dir (@dirs) {
	    File::Path::rmtree($dir);
	}
	
	#just in case, this will help remove temp_dirs
	my @d = <$the_void/*.temp_dir>;
	foreach my $d (@d){
	    File::Path::rmtree($d) if (-d $d);
	}
    }
    
    $self->{continue_flag} = $continue_flag;
}
#-------------------------------------------------------------------------------
#==build log file for current settings
sub _write_new_log {
   my $self = shift;

   my $log_file = $self->{file_name};
   
   my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};

   return if ($self->{continue_flag} <= 0);

   open (LOG, "> $log_file");

   #log control file options
   my $cwd = ($self->{CWD}) ? $self->{CWD} : Cwd::cwd();
   foreach my $key (@ctl_to_log) {
      my $ctl_val = '';
      if(defined $CTL_OPTIONS{$key}){
	 $ctl_val = $CTL_OPTIONS{$key} ;
	 $ctl_val =~ s/^$cwd\/*//;
	 if($key eq 'appl'){
	     #don't care about order
	     my @set = sort @{$CTL_OPTIONS{appl}};
	     $ctl_val = join(',', @set);
	 }
      }	  
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

   my $log_file = $self->{file_name};
   my $cwd = ($self->{CWD}) ? $self->{CWD} : Cwd::cwd();

   #this line hides unnecessarilly deep directory details
   #this is important for maker webserver security
   if(defined $key && defined $type){
       $key =~ s/^$cwd\/*// if($type =~ /^STARTED$|^FINISHED$/);
   }

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
sub report_status {
   my $self = shift;
   my $flag = $self->{continue_flag};
   my $die_count = $self->{die_count};
   my $seq_id = $self->{params}->{seq_id};
   my $seq_out_name = Fasta::seqID2SafeID($seq_id);
   my $out_dir = $self->{params}->{out_dir};
   my $fasta_ref = $self->{params}->{fasta_ref};
   #my $length = $self->{params}->{seq_length};

   #maintain lock only if status id positive (possible race condition?)
   $self->{LOCK}->maintain(30) if($flag > 0);

   if($flag == 0){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig has already been processed!!\n",
		   "mpi_iprscan will now skip to the next contig.\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",		   
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now starting the contig!!\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 2){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
                   "All contig related data will be erased before continuing!!\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == -1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig failed $die_count time!!\n",
		   "mpi_iprscan will not try again!!\n",
		   "The contig will be stored in a fasta file that you can use for debugging.\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",
		   "FASTA: $out_dir/$seq_out_name.died.fasta\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);

      open (my $DFAS, "> $out_dir/$seq_out_name.died.fasta");
      print $DFAS $$fasta_ref;
      close ($DFAS);
   }
   elsif($flag == -2){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Skipping the contig because it is too short!!\n",
		   "SeqID: $seq_id\n",
                   #"Length: $length\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == -3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Another instance of is processing this contig!!\n",
                   "SeqID: $seq_id\n",
                   #"Length: $length\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
  }
   else{
      die "ERROR: No valid continue flag\n";
   }
}
#-------------------------------------------------------------------------------
sub get_continue_flag {
   my $self = shift;
   my $flag = $self->{continue_flag};
   my $message;

   if($flag == 0){
      $message = 'FINISHED'; #already ran don't re-run
   }
   elsif($flag == 1){
      $message = 'STARTED'; #run with current settings
   }
   elsif($flag == 2){
      $message = 'RETRY_CLEAN'; #re-run previously died
   }
   elsif($flag == 3){
       $message = 'RETRY'; #re-run previously died
   }
   elsif($flag == -1){
      $message = 'DIED_SKIPPED_PERMANENT'; #don't re-run, it died
   }
   elsif($flag == -2){
      $message = 'SKIPPED_SMALL'; #not finished but skipped
   }
   elsif($flag == -3){
       $message = ''; #no short message, as contig is running elsewhere
   }
   else{
      die "ERROR: No valid continue flag\n";
   }

   return $flag, $message;
}
#-------------------------------------------------------------------------------
#used to pull lock off before serialization of runlog object
sub strip_off_lock {
    my $self = shift;
    my $lock = $self->{LOCK};

    $self->{LOCK} = undef;

    return $lock;
}
#-------------------------------------------------------------------------------
sub unlock {
    my $self = shift;

    $self->{LOCK}->unlock if(defined $self->{LOCK});
}
#-------------------------------------------------------------------------------
sub DESTROY {
    my $self = shift;

    $self->unlock;
}
#-------------------------------------------------------------------------------
1;
