#-------------------------------------------------------------------------------
#------                           runlog                               ---------
#-------------------------------------------------------------------------------
package runlog;

BEGIN {
    @AnyDBM_File::ISA = qw(DB_File GDBM_File NDBM_File SDBM_File);
}

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use IO::File;
use Fasta;
use File::NFSLock;
use AnyDBM_File;
use GI;

@ISA = qw();
$VERSION = 0.1;

#===make list of internal variables to log
my @ctl_to_log = ('genome_gff',
		  'other_gff',
		  'est',
		  'est_reads',
		  'altest',
		  'est_gff',
		  'altest_gff',
		  'protein',
		  'protein_gff',
		  'model_org',
		  'repeat_protein',
		  'rmlib',
		  'rm_gff',
		  'organism_type',
		  'predictor',
		  'snaphmm',
		  'gmhmm',
		  'augustus_species',
		  'fgenesh_par_file',
		  'model_gff',
		  'pred_gff',
		  'max_dna_len',
		  'split_hit',
		  'pred_flank',
		  'min_protein',
		  'AED_threshold',
		  'single_exon',
		  'single_length',
		  'keep_preds',
		  'map_forward',
		  'alt_peptide',
		  'evaluate',
		  'blast_type',
		  'softmask',
		  'pcov_blastn',
		  'pid_blastn',
		  'eval_blastn',
		  'bit_blastn',
		  'pcov_rm_blastx',
		  'pid_rm_blastx',
		  'eval_rm_blastx',
		  'bit_rm_blastx',
		  'pcov_blastx',
		  'pid_blastx',
		  'eval_blastx',
		  'bit_blastx',
		  'pcov_tblastx',
		  'pid_tblastx',
		  'eval_tblastx',
		  'bit_tblastx',
		  'ep_score_limit',
		  'en_score_limit',
		  'enable_fathom',
		  'unmask',
		  'model_pass',
		  'est_pass',
		  'altest_pass',
		  'protein_pass',
		  'rm_pass',
		  'other_pass',
		  'pred_pass',
		  'run'
		 );

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
       $self->_compare_and_clean();
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
   $self->{CWD} = $self->{CTL_OPTIONS}->{CWD};

   print STDERR "\n\n\n--Next Contig--\n\n" unless($main::quiet);

   #make sure the log is only accesed by this process
   if(my $lock = new File::NFSLock($self->{file_name}, 'NB', undef, 40)){
       my $min_contig = $self->{CTL_OPTIONS}->{min_contig};
       my $length = $self->{params}->{seq_length};
   
       #log can be processed and compared
       if($length >= $min_contig){
	   $self->{continue_flag} = 1;
	   $self->{LOCK} = $lock;

	   return 1;
       }
       #skip if this is a short contig
       else{
	   $self->{continue_flag} = -2; # -2 is skip signal
	   $lock->unlock;
	       
	   return 0;
       }
   }
   #belongs to another process so skip
   else{
       $self->{continue_flag} = -3;

       return 0;
   }
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
sub _compare_and_clean {
    my $self = shift;

    #make sure lock is till mine
    if(! $self->{LOCK}->still_mine){
	$self->{continue_flag} = -3;
	return;
    }

    my $CWD = $self->{CWD};
    my $the_void = $self->{params}->{the_void};
    my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};

    #get gff file name
    my $log_file = $self->{file_name};
    my $gff_file = $the_void;
    my $out_base = $the_void;
    $out_base =~ s/theVoid\.[^\/]+$//;
    $gff_file =~ s/theVoid\.([^\/]+)$/$1.gff/;
    my $name = $1;

    #===id files types that can not be re-used
    my $continue_flag = $self->{continue_flag}; #signals whether to continue with this seq
    
    my %logged_vals = %{$self->{old_log}}; #values from existing log
    my %rm_key; #key to what type of old files to remove
    my @files; #list of files to remove
    my @dirs; #list of directories to remove
    
    if ($continue_flag > 0 && -e $log_file) {

	die "ERROR: Database timed out in runlog::_clean_files\n\n"
	    unless (my $lock = new File::NFSLock($CTL_OPTIONS{SEEN_file}, 'EX', 120, 30));

	#seen is set in ds_utility.pm
	my %SEEN;
	tie (%SEEN, 'AnyDBM_File', $CTL_OPTIONS{SEEN_file});

	if (exists $logged_vals{DIED}) {
	    if($CTL_OPTIONS{force} && ! defined $SEEN{$name}){
		$self->{die_count} = 0; #reset the die count
		$continue_flag = 1; #re-run
		$rm_key{force}++;
	    }
	    elsif($CTL_OPTIONS{always_try} && ! defined $SEEN{$name}){
		$self->{die_count} = 0; #reset the die count
		$continue_flag = 1; #re-run
	    }
	    else{
		$continue_flag = ($CTL_OPTIONS{clean_try}) ? 2 : 3;	#rerun died
		$continue_flag = -1 if($self->{die_count} > $CTL_OPTIONS{retry}); #only let die up to count
		$rm_key{retry}++ if ($continue_flag == 2);
	    }
	}
	elsif ($CTL_OPTIONS{force} && ! defined $SEEN{$name}) {
	    $rm_key{force}++;
	    $continue_flag = 1;	#run/re-run
	}
	elsif ($CTL_OPTIONS{again} && ! defined $SEEN{$name}){
	    $continue_flag = 1; #run/re-run
	    $rm_key{gff}++;
	}
	else {
	    $continue_flag = 0 if (-e $gff_file); #don't re-run finished
	}

	$SEEN{$name} = $continue_flag;
	untie %SEEN;
	$lock->unlock;

	if($continue_flag >= -1){
	    #CHECK CONTROL FILE OPTIONS FOR CHANGES
	    my $cwd = ($CWD) ?$CWD : Cwd::getcwd();

	    foreach my $key (@ctl_to_log) {
		#changes to hidden key only important when nothing else changed
		if($key eq 'run'){
		    next if(exists $rm_key{gff});
		}

		#these are only sometimes important
		if($key =~ /^est_pass$|^altest_pass$|^protein_pass$|^rm_pass$/ ||
		   $key =~ /^pred_pass$|^model_pass$|^other_pass$/
		   ){
		    next unless($CTL_OPTIONS{genome_gff});
		    my $old = (exists $logged_vals{CTL_OPTIONS}{genome_gff}) ? 
			$logged_vals{CTL_OPTIONS}{genome_gff} : '';
		    $old =~ s/^$cwd\/*//;
		    my $new = $CTL_OPTIONS{genome_gff};
		    $new =~ s/^$cwd\/*//;

		    #only continue if change not already happening
		    #because of a new gff3 file
		    next unless($old eq $new);
		}

                #these are only sometimes important
                if($key =~ /^map_forward$/){
                    next unless($CTL_OPTIONS{genome_gff} || $CTL_OPTIONS{model_gff});
                }
		
		my $log_val = '';
		if(defined $logged_vals{CTL_OPTIONS}{$key}){	       
		    $log_val = $logged_vals{CTL_OPTIONS}{$key};
		    if($key eq 'repeat_protein'){
			#don't care about absolute location
			$log_val =~ s/.*\/(te_proteins.fasta)$/$1/;
		    }
		    elsif($key eq 'run'){
			#don't care about order
			my @set = split(',', $logged_vals{CTL_OPTIONS}{run});
			@set = sort @set;
			$log_val = join(',', @set);
		    }
		}
		
		my $ctl_val = '';
		if(defined $CTL_OPTIONS{$key}){
		    $ctl_val = $CTL_OPTIONS{$key};
		    $ctl_val =~ s/^$cwd\/*//;
		    if($key eq 'repeat_protein'){
			#don't care about absolute location
			$ctl_val =~ s/.*\/(te_proteins.fasta)$/$1/;
		    }
		    elsif($key eq 'run'){
			#don't care about order
			my @set = sort @{$CTL_OPTIONS{_run}};
			$ctl_val = join(',', @set);
		    }
		}

		#organism_type was always eukaryotic before and not logged
		if($key eq 'organism_type' && ! $log_val){
		    $log_val = 'eukaryotic';
		}

		#softmask was always true before and not logged
		if($key eq 'softmask' && $log_val eq ''){
		    $log_val = 1;
		}

		#AED_thrshold was always 1 before and not logged
		if($key eq 'AED_threshold' && $log_val eq ''){
		    $log_val = 1;
		}

		#if previous log options are not the same as current control file options
		if ($log_val ne $ctl_val) {
		    
		    print STDERR "MAKER WARNING: Control file option \'$key\' has changed\n".
			"Old:$log_val\tNew:$ctl_val\n\n" unless($main::qq);
		    
		    $continue_flag = 1; #re-run because ctlopts changed
		    
		    $rm_key{gff}++; #always rebuild gff when some option has changed
		    
		    #certain keys that don't affect preds
		    if($key ne 'evaluate' &&
		       $key ne 'enable_fathom' &&
		       $key ne 'keep_preds' &&
		       $key ne 'other_pass' &&
		       $key ne 'other_gff' &&
		       $key ne 'map_forward'
		      ){
			$rm_key{preds}++; #almost all changes affect final predictions
		    }
		    
		    if ($key eq 'organism_type') {
			$rm_key{all}++;
		    }
		    
		    if ($key eq 'max_dna_len') {
			$rm_key{all}++;
		    }
		    
		    if ($key eq 'rm_gff' ||
			$key eq 'model_org' ||
			$key eq 'rmlib'		   
			) {
			$rm_key{all}++;
		    }
		    
		    if ($key eq 'repeat_protein' ||
			$key eq 'pcov_rm_blastx' ||
			$key eq 'pcid_rm_blastx' ||
			$key eq 'eval_rm_blastx' ||
			$key eq 'bit_rm_blastx' ||
			$key eq 'blast_type'
			) {
			$rm_key{all_but}++;
		    }
		    
		    if ($key eq 'snaphmm') {
			$rm_key{snap}++;
		    }
		    
		    if ($key eq 'augustus_species') {
			$rm_key{augustus}++;
		    }
		    
		    if ($key eq 'fgenesh_par_file') {
			$rm_key{fgenesh}++;
		    }
		    
		    if ($key eq 'gmhmm') {
			$rm_key{genemark}++;
		    }
		    
		    if ($key eq 'split_hit' ||
			$key eq'ep_score_limit' ||
			$key eq'en_score_limit'
			) {
			$rm_key{e_exonerate}++;
			$rm_key{p_exonerate}++;
		    }
		    
		    if ($key eq 'protein' ||
			$key eq 'alt_peptide' ||
			$key eq 'eval_blastx' ||
			$key eq 'softmask'
			) {
			$rm_key{blastx}++;
			$rm_key{p_exonerate}++;
		    }
		    
		    if ($key eq 'est') {
			$rm_key{est_blastn}++;
			$rm_key{e_exonerate}++;
		    }
		    
		    if ($key eq 'est_reads') {
			$rm_key{read_blastn}++;
		    }
		    
		    if ($key eq 'eval_blastn' ||
			$key eq 'split_hit'
			) {
			$rm_key{blastn}++;
			$rm_key{e_exonerate}++;
		    }
		    
		    if ($key eq 'altest' ||
			$key eq 'eval_tblastx' ||
			$key eq 'split_hit' ||
			$key eq 'softmask'
			) {
			$rm_key{tblastx}++;
		    }
		}
	    }
	    
	    #CHECK FOR FILES THAT DID NOT FINISH
	    while (my $key = each %{$logged_vals{STARTED}}) {
		if (! exists $logged_vals{FINISHED}{$key}) {
		    print STDERR "MAKER WARNING: The file $key\n".
			"did not finish on the last run and must be erased\n" unless($main::qq);

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
	
	#-print file type specific warnings and delete files
	if (exists $rm_key{force}) {
	    print STDERR "MAKER WARNING: All old files will be erased before continuing\n" unless($main::qq);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    
	    #remove evaluator output
	    File::Path::rmtree("$out_base/evaluator");
	    
	    #remove all other files
	    unlink($gff_file) if(-e $gff_file);
	    my @f = <$out_base/*.fasta>;
	    push (@files, @f);
	}
	elsif (exists $rm_key{retry}) {
	    print STDERR "MAKER WARNING: Old data must be removed before re-running this sequence\n" unless($main::qq);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    unlink($gff_file) if(-e $gff_file);
	    
	    #remove evaluator output
	    File::Path::rmtree("$out_base/evaluator");
	}
	elsif (exists $rm_key{all}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of all old data impossible\n".
		"All old files will be erased before continuing\n" unless($main::qq);
	    
	    #delete everything in the void
	    File::Path::rmtree($the_void);
	    File::Path::mkpath($the_void);
	    unlink($gff_file) if(-e $gff_file);
	    
	    #remove evaluator output
	    File::Path::rmtree("$out_base/evaluator");
	}
	elsif (exists $rm_key{all_but}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of all but RepeatMasker data impossible\n".
		"All old non-RepeatMasker files will be erased before continuing\n" unless($main::qq);
	    
	    #delete everything in the void
	    my @f = <$the_void/*>;
	    @f = grep(!/(\.out|\.cat|\.tbl)$/, @f);
	    
	    #delete files in the void
	    foreach my $f (@f) {
		unlink($f) if(-f $f);
		File::Path::rmtree($f) if(-d $f);
	    }
	}
	else {
	    if (exists $rm_key{preds}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of hint based predictions impossible\n".
		    "Old hint based prediction files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*auto_annotator*>;
		push (@files, @f);
	    }
	    if (exists $rm_key{snap}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old SNAP data impossible\n".
		    "Old SNAP files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*snap*>;
		push (@files, @f);
	    }
	    if (exists $rm_key{augustus}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old Augustus data impossible\n".
		    "Old Augustus files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*augustus*>;
		push (@files, @f);
	    }
	    if (exists $rm_key{fgenesh}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old FGENESH data impossible\n".
		    "Old FGENESH files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*fgenesh*>;
		push (@files, @f);
	    }
	    if (exists $rm_key{genemark}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old GeneMark data impossible\n".
		    "Old GeneMark files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*genemark*>;
		push (@files, @f);
	    }
	    if (exists $rm_key{blastn}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of all old EST Blastn data impossible\n".
		    "Old EST Blastn files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*blastn*>;
		foreach my $f (@f){
		    push (@files, $f) if (-f $f);
		    push (@dirs, $f) if (-d $f);
		}
	    }
	    else{
		if (exists $rm_key{est_blastn}) {
		    print STDERR "MAKER WARNING: Changes in control files make re-use of assembled EST Blastn data impossible\n".
			"Old EST Blastn files will be erased before continuing\n" unless($main::qq);
		    
		    my @f = <$the_void/*est_blastn*>;
		    foreach my $f (@f){
			push (@files, $f) if (-f $f);
			push (@dirs, $f) if (-d $f);
		    }
		}
		elsif (exists $rm_key{read_blastn}) {
		    print STDERR "MAKER WARNING: Changes in control files make re-use of unassembled EST Blastn data impossible\n".
			"Old EST reads Blastn files will be erased before continuing\n" unless($main::qq);
		    
		    my @f = <$the_void/*read_blastn*>;
		    foreach my $f (@f){
			push (@files, $f) if (-f $f);
			push (@dirs, $f) if (-d $f);
		    }
		}
	    }
	    
	    if (exists $rm_key{tblastx}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old tBlastx data impossible\n".
		    "Old tBlastx files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*tblastx*>;
		foreach my $f (@f){
		    push (@files, $f) if (-f $f);
		    push (@dirs, $f) if (-d $f);
		}
	    }
	    
	    if (exists $rm_key{blastx}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old Blastx data impossible\n".
		    "Old Blastx files will be erased before continuing\n" unless($main::qq);
	 
		my @f = <$the_void/*blastx*>;

		my ($te) = $CTL_OPTIONS{repeat_protein} =~ /([^\/]+)$/;

		if($te){
		    $te =~ s/\.fasta$//;    
		    @f = grep { ! /$te\.blastx$/} @f;
		}
		    
		foreach my $f (@f){
		    push (@files, $f) if (-f $f);
		    push (@dirs, $f) if (-d $f);
		}

	    }
	    
	    if (exists $rm_key{e_exonerate}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old EST Exonerate data impossible\n".
		    "Old EST Exonerate files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*est_exonerate*>;
		push (@files, @f);
	    }
	    
	    if (exists $rm_key{p_exonerate}) {
		print STDERR "MAKER WARNING: Changes in control files make re-use of old protein Exonerate data impossible\n".
		    "Old protein Exonerate files will be erased before continuing\n" unless($main::qq);
		
		my @f = <$the_void/*p_exonerate*>;
		push (@files, @f);
	    }
	    
	    if (exists $rm_key{gff}) {
		print STDERR "MAKER WARNING: Any preexisting GFF3 and fasta files for this contig must now be removed.\n"
		    unless($main::qq);
		push (@files, $gff_file);
		push (@files, @{[<$out_base/evaluator/*.eva>]});
		push (@files, @{[<$out_base/*maker*.fasta>]});
	    }

	    #see if died fasta exists
	    my $d_fasta = $gff_file;
	    $d_fasta =~ s/\.gff$/.died.fasta/;
	    if(-e $d_fasta){
                push (@files, $d_fasta);
            }

	    #delete files in the void
	    foreach my $file (@files) {
		unlink($file);
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
    }
    
    $self->{continue_flag} = $continue_flag;
}
#-------------------------------------------------------------------------------
#==build log file for current settings
sub _write_new_log {
   my $self = shift;

   my $CWD = $self->{CWD};
   my $log_file = $self->{file_name};
   
   my %CTL_OPTIONS = %{$self->{CTL_OPTIONS}};

   return if ($self->{continue_flag} <= 0);

   open (LOG, "> $log_file");

   #log control file options
   my $cwd = ($CWD) ? $CWD : Cwd::getcwd();
 
   foreach my $key (@ctl_to_log) {
      my $ctl_val = '';
      if(defined $CTL_OPTIONS{$key}){
	 $ctl_val = $CTL_OPTIONS{$key} ;
	 $ctl_val =~ s/^$cwd\/*//;
	 if($key eq 'repeat_protein'){
	    #don't care about absolute location
	    $ctl_val =~ s/.*\/(te_proteins.fasta)$/$1/;
	 }
	 elsif($key eq 'run'){
	     #don't care about order
	     my @set = sort @{$CTL_OPTIONS{_run}};
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

   my $CWD = $self->{CWD};
   my $log_file = $self->{file_name};
   my $cwd = ($CWD) ?$CWD : Cwd::getcwd();

   #this line hides unnecessarilly deep directory details
   #this is important for maker webserver security
   #also important when moving work directory around
   if(defined $key && defined $type){
       $key =~ s/^$cwd\/*// if($type =~ /^STARTED$|^FINISHED$/);
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
		   "Maker will now skip to the next contig.\n",
		   "Run maker with the -f flag to force Maker to recompute all contig data.\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",		   
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == 1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now starting the contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == 2){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
                   "All contig related data will be erased before continuing!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == 3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "Retry: $die_count!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == -1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig failed $die_count time!!\n",
		   "Maker will not try again!!\n",
		   "The contig will be stored in a fasta file that you can use for debugging.\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "FASTA: $out_dir/$seq_out_name.died.fasta\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);

      open (my $DFAS, "> $out_dir/$seq_out_name.died.fasta");
      print $DFAS $$fasta_ref;
      close ($DFAS);
   }
   elsif($flag == -2){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Skipping the contig because it is too short!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == -3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Another instance of maker is processing this contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
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
sub are_same_opts {
    shift @_ if(ref($_[0]) eq 'runlog'); #can be method or function
    my %CTL_OPT = %{shift @_};
    my @files = @_;

    my %LOG_OPT;
    if(ref($files[0]) eq 'HASH'){
	%LOG_OPT = %{$files[0]};
    }
    elsif(ref($files[0]) eq 'ARRAY'){
	%LOG_OPT = GI::load_control_files($files[0], {parse => 1})
    }
    else{
	%LOG_OPT = GI::load_control_files(\@files, {parse => 1})
    }

    foreach my $key (@ctl_to_log){
	my $ctl = (defined $CTL_OPT{$key}) ? $CTL_OPT{$key} : '';
	my $log = (defined $LOG_OPT{$key}) ? $LOG_OPT{$key}: '';

	$ctl = join(',', sort split(',', $ctl));
	$log = join(',', sort split(',', $log));

	if($ctl ne $log){
	    return 0;
	}
    }

    return 1;
}
#-------------------------------------------------------------------------------
sub DESTROY {
    my $self = shift;

    $self->unlock;
}
1;
