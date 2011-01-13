#-------------------------------------------------------------------------------
#------                       iprscan::runlog                          ---------
#-------------------------------------------------------------------------------
package iprscan::runlog;

BEGIN {
    @AnyDBM_File::ISA = qw(DB_File GDBM_File NDBM_File SDBM_File);
}

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use Fasta;
use File::NFSLock;
use AnyDBM_File;
use URI::Escape;

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
       $self->_compare_and_clean();
       $self->_write_new_log();
    }

    $self->report_status();

    return $self;
}
#-------------------------------------------------------------------------------
sub _initialize {
   my $self = shift;
   $self->{CTL_OPT} = shift;
   $self->{params} = shift;
   $self->{file_name} = shift || "run.log";
   $self->{CWD} = $self->{CTL_OPT}->{CWD} || Cwd::cwd();

   print STDERR "\n\n\n--Next Contig--\n\n" unless($main::quiet);

   #lock must be persitent in object or it is detroyed outside of block
   if(my $lock = new File::NFSLock($self->{file_name}, 'NB', undef, 90)){
       $self->{continue_flag} = 1;
       $self->{LOCK} = $lock;
       return 1;
   }
   else{
       $self->{continue_flag} = -3;
       return 0;
   }
}
#-------------------------------------------------------------------------------
sub _load_old_log {
   my $self = shift;

   $self->{die_count} = 0;
   $self->{old_shared_id} = '';

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
         $self->{old_shared_id} = $key if($type eq 'SHARED_ID');
      }
      close(IN);
   }

   $self->{old_log} = \%logged_vals;
}
#-------------------------------------------------------------------------------
sub _compare_and_clean {
    my $self = shift;

    #make sure lock is till mine
    if(! $self->{LOCK}->still_mine){
        $self->{continue_flag} = -3;
        return;
    }

    $self->{CWD} = Cwd::getcwd() if(!$self->{CWD});
    my $CWD = $self->{CWD};
    my $the_void = $self->{params}->{the_void};
    my %CTL_OPT = %{$self->{CTL_OPT}};
    
    #get ipr file name
    my $log_file = $self->{file_name};
    my $ipr_file = $the_void;
    my $out_base = $the_void;
    $out_base =~ s/theVoid\.[^\/]+$//;
    $ipr_file =~ s/theVoid\.([^\/]+)$//;
    $ipr_file .= uri_escape($1, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.') .".ipr";
    my $name = $1;

    #===id files types that can not be re-used
    my $continue_flag = $self->{continue_flag}; #signal whether to continue with this seq
    my %rm_key; #key to what type of old files to remove
    my @files; #list of files to remove
    my @dirs; #list of directories to remove

    #evaluate existing runlog for contig
    if (-e $log_file) {
        my %logged_vals = %{$self->_load_old_log}; #values from existing log                                                                  
        my $shared = $self->{old_shared_id} eq $CTL_OPT{_shared_id};

        if($self->{params}->{seq_length} < $self->{CTL_OPT}->{min_contig} && $shared){
            $self->{continue_flag} = -4; #handled elsewhere, don't run, don't report
            return;
        }
        elsif(-e $ipr_file && $shared){
            $self->{continue_flag} = -4; #handled elsewhere, don't run, don't report
            return;
        }
        elsif($self->{params}->{seq_length} < $self->{CTL_OPT}->{min_contig} && !$shared){
            $continue_flag = -2; #skip short
            $self->{die_count} = 0; #reset the die count
            $rm_key{force}++; #will destroy ipr and fastas for this contig
        }
        elsif ($CTL_OPT{force} && !$shared) {
            $rm_key{force}++;
            $self->{die_count} = 0; #reset the die count
            $continue_flag = 1; #run/re-run
        }
        elsif ($CTL_OPT{again} && !$shared){
            $rm_key{ipr}++;
            $self->{die_count} = 0; #reset the die count
            $continue_flag = 1; #run/re-run
        }
        elsif(-e $ipr_file){
            $continue_flag = 0; #don't re-run finished
        }
        elsif ($CTL_OPT{always_try}){
            $self->{die_count} = 0; #reset the die count
            $continue_flag = 1; #re-run
        }
        elsif ($logged_vals{DIED}) {
            $continue_flag = ($CTL_OPT{clean_try}) ? 2 : 3;     #rerun died
            $continue_flag = -1 if($self->{die_count} > $CTL_OPT{retry}); #only let die up to count
            $continue_flag = -4 if($shared && $continue_flag == -1);
            $rm_key{force}++ if($continue_flag == 2);
        }

        my %skip;
        if(!$rm_key{force}){
            #CHECK CONTROL FILE OPTIONS FOR CHANGES
            foreach my $key (@ctl_to_log) {
                my $log_val = '';
                if(defined $logged_vals{CTL_OPTIONS}{$key}){
                    $log_val = $logged_vals{CTL_OPTIONS}{$key};
                    $log_val =~ s/(^|\,)$CWD\/*/$1/g;
                    $log_val = join(',', sort split(',', $log_val));
		}
		
		my $ctl_val = '';
		if(defined $CTL_OPT{$key}){
		    $ctl_val = $CTL_OPT{$key};
		    $ctl_val =~ s/(^|\,)$CWD\/*/$1/g;
		    $ctl_val = join(',', sort split(',', $ctl_val));
		}
		
		#if previous log options are not the same as current control file options
		if ($log_val ne $ctl_val) {	    
		    print STDERR "WARNING: Control file option \'$key\' has changed\n".
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
		    print STDERR "WARNING: The file $key\n".
			"did not finish on the last run and must be erased\n" unless($main::quiet);
		    push(@files, $key);

                    #this next step will get both temp directories and files that
                    #have the incorrect location in the log but are in theVoid
		    $key =~ /([^\/]+)$/;
		    my $rm_f_name = $1;

                    if(! -e $key && -e "$the_void/$rm_f_name"){
                        push(@files, $rm_f_name);
                    }


		    $rm_f_name =~ s/\.fasta$//;
		    my @d = <$the_void/*$rm_f_name*>;
		    foreach my $d (@d){
			push (@dirs, $d) if (-d $d);
		    }
		}
	    }
	}

	#===Remove files that can no longer be re-used
	if(! $self->{LOCK}->still_mine){
	    warn "ERROR: Lock broken in runlog\n";
	    $self->{continue_flag} = -3;
	    return;
	}

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
	elsif (exists $rm_key{all}) {
	    print STDERR "WARNING: Changes in control files make re-use of all old data impossible\n".
		"All old files will be erased before continuing\n" unless($main::quiet);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);

	    #remove all other files
	    unlink($ipr_file) if(-e $ipr_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	}
	else{
	    if (exists $rm_key{ipr}) {
		print STDERR "WARNING: Any preexisting IPR files for this contig must now be removed.\n" unless($main::quiet);
		push (@files, $ipr_file);
		push (@files, @{[<$out_base/*.fasta>]});
	    }
	    else{
		#see if died fasta exists
		my $d_fasta = $ipr_file;
		$d_fasta =~ s/\.ipr$/.died.fasta/;
		if(-e $d_fasta){
		    push (@files, $d_fasta);
		}
	    }
	}
    
	#delete files in the void
	foreach my $file (@files) {
	    unlink($file) if(-e $file);
	}
	
	#delete directories in the void
	foreach my $dir (@dirs) {
	    File::Path::rmtree($dir);
	}
    }
    
    $self->{continue_flag} = $continue_flag;
}
#-------------------------------------------------------------------------------
#==build log file for current settings
sub _write_new_log {
   my $self = shift;

   return if($self->{continue_flag} <= -3);

   $self->{CWD} = Cwd::getcwd() if(!$self->{CWD});
   my $CWD = $self->{CWD};
   my $log_file = $self->{file_name};
   my %CTL_OPT = %{$self->{CTL_OPT}};

   open (LOG, "> $log_file");

   #log control file options
   print LOG "SHARED_ID\t$CTL_OPT{_shared_id}\t\n";
   foreach my $key (@ctl_to_log) {
       my $ctl_val = '';
       if(defined $CTL_OPT{$key}){
	   $ctl_val = $CTL_OPT{$key};
	   $ctl_val =~ s/(^|\,)$CWD\/*/$1/g;
	   $ctl_val = join(',', sort split(',', $ctl_val));
       }
       print LOG "CTL_OPTIONS\t$key\t$ctl_val\n";
   }
   $self->add_entry("DIED","COUNT",$self->get_die_count()) if($self->{continue_flag} == -1);
   close(LOG);
}
#-------------------------------------------------------------------------------
sub add_entry {
   my $self = shift;

   my $type  = shift;
   my $key   = shift;
   my $value = shift;

   $self->{CWD} = Cwd::getcwd() if(!$self->{CWD});
   my $CWD = $self->{CWD};

   my $log_file = $self->{file_name};

   #this line hides unnecessarilly deep directory details
   #this is important for maker webserver security
   #also important when moving work directory around
   if(defined $key && defined $type){
       $key =~ s/^$CWD\/*// if($type =~ /^STARTED$|^FINISHED$/);
   }

   open(my $LOG, ">> $log_file");
   print $LOG "$type\t$key\t$value\n";
   close($LOG);
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
   my $length = $self->{params}->{seq_length};

   #maintain lock only if status id positive (possible race condition?)
   $self->{LOCK}->maintain(30) if($flag > 0);

   if($flag == 0){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig has already been processed!!\n",
		   "iprscan will now skip to the next contig.\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",		   
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now starting the contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 2){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
                   "All contig related data will be erased before continuing!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == 3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
   }
   elsif($flag == -1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig failed $die_count time!!\n",
		   "iprscan will not try again!!\n",
		   "The contig will be stored in a fasta file that you can use for debugging.\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
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
                   "Length: $length\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::quiet);
  }
   elsif($flag == -4){
       #do nothing
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
   elsif($flag == -4){
       $message = ''; #no short message, as contig was handled elsewhere
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
sub compare_comma_lists{
    my $val1 = shift;
    my $val2 = shift;

    my @set1 = map {$_ =~ /^([^\:]+)/; $1} split(',', $val1);
    my @set2 = map {$_ =~ /^([^\:]+)/; $1} split(',', $val2);

    my %count;

    foreach my $v (@set1, @set2){
        $count{$v}++;
    }

    my @change = map {/([^\/]+)$/; $1} grep {$count{$_} == 1} keys %count;
    my @keep = map {/([^\/]+)$/; $1} grep {$count{$_} != 1} keys %count;

    return (\@change, \@keep);
}
#-------------------------------------------------------------------------------
sub DESTROY {
    my $self = shift;

    $self->unlock;
}
#-------------------------------------------------------------------------------
1;
