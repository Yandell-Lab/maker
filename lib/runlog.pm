#-------------------------------------------------------------------------------
#------                           runlog                               ---------
#-------------------------------------------------------------------------------
package runlog;

use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use IO::File;
use Fasta;
use File::NFSLock;
use GI;
use URI::Escape;
use Carp;

@ISA = qw();
$VERSION = 0.1;

#===make list of internal variables to log
my @ctl_to_log = ('maker_gff',
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
		  'est2genome',
		  'altest2genome',
		  'snaphmm',
		  'gmhmm',
		  'augustus_species',
		  'fgenesh_par_file',
		  'run_evm',
		  'model_gff',
		  'pred_gff',
		  'max_dna_len',
		  'split_hit',
		  'min_intron',
		  'pred_flank',
		  'pred_stats',
		  'min_protein',
		  'AED_threshold',
		  'single_exon',
		  'single_length',
		  'keep_preds',
		  'map_forward',
		  'est_forward',
		  'correct_est_fusion',
		  'alt_splice',
		  'always_complete',
		  'alt_peptide',
		  'evaluate',
		  'blast_type',
		  'use_rapsearch',
		  'softmask',
		  'pcov_blastn',
		  'pid_blastn',
		  'eval_blastn',
		  'bit_blastn',
		  'depth_blastn',
		  'pcov_rm_blastx',
		  'pid_rm_blastx',
		  'eval_rm_blastx',
		  'bit_rm_blastx',
		  'pcov_blastx',
		  'pid_blastx',
		  'depth_blastx',
		  'eval_blastx',
		  'bit_blastx',
		  'pcov_tblastx',
		  'pid_tblastx',
		  'eval_tblastx',
		  'bit_tblastx',
		  'depth_tblastx',
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

   #make sure the log is only accesed by this process
   if(my $lock = new File::NFSLock($self->{file_name}, 'NB', undef, 40)){
	   $self->{continue_flag} = 1;
	   $self->{LOCK} = $lock;
	   return 1;
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
   $self->{old_shared_id} = '';

   print STDERR "Processing run.log file...\n"
       if(-f $self->{file_name} && !$main::quiet);

   my %logged_vals;
   my @files = ($self->{file_name});
   foreach my $log_file (@files){       
       my %stat;
       if (-f $log_file){ #load log file if available
	   open (IN, "< $log_file");      
	   while( defined (my $line = <IN>)){
	       chomp $line;
	       
	       my ($type, $key, $value) = split ("\t", $line);

	       if($type eq 'FINISHED' && defined $stat{STARTED}{$key}){
		   #delete or hash can become very large
		   delete($stat{STARTED}{$key});
		   delete($stat{FINISHED}{$key});
		   next;
	       }
	       elsif($type eq 'STARTED' && defined $stat{FINISHED}{$key}){
		   #delete or hash can become very large
		   delete($stat{STARTED}{$key});
		   delete($stat{FINISHED}{$key});
		   next;
	       }
	       elsif($type eq '###'){
		   #progress checkpoint reached, all analyses completed up to this point
		   $stat{STARTED} = {};
		   $stat{FINISHED} = {};
		   next;
	       }
	       
	       if($type eq 'FINISHED' || $type eq 'STARTED' || $type eq 'LOGCHILD'){
		   $stat{$type}{$key} = defined($value) ? $value : '';
	       }
	       else{
		   $key = 'maker_gff' if($key eq 'genome_gff'); #for backwards compatability
		   $logged_vals{$type}{$key} = defined($value) ? $value : '';
	       }
	       
	       if($type eq 'DIED' && $key eq 'COUNT'){
		   $self->{die_count} = $value;
	       }
	       elsif($type eq 'SHARED_ID'){
		   $self->{old_shared_id} = $key;
	       }
	   }
	   close(IN);
   
	   while(my $type = each %stat){
	       while(my $key = each %{$stat{$type}}){
		   if($type eq 'LOGCHILD'){
		       push(@files, $key);
		       next;
		   }
		   $logged_vals{$type}{$key} = $stat{$type}{$key};
	       }
	   }
       }
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

    #get gff file name
    my $log_file = $self->{file_name};
    my $gff_file = $the_void;
    my $out_base = $the_void;
    $out_base =~ s/theVoid\.[^\/]+$//;
    $gff_file =~ s/theVoid\.([^\/]+)$//;
    $gff_file .= uri_escape($1, '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:\.') .".gff";
    my $name = $1;

    #===id files types that can not be re-used
    my $continue_flag = $self->{continue_flag}; #signals whether to continue with this seq
    my %rm_key; #key to what type of old files to remove
    my %skip;
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
	elsif(-e $gff_file && $shared){
	    $self->{continue_flag} = -4; #handled elsewhere, don't run, don't report
	    return;
	}
	elsif($self->{params}->{seq_length} < $self->{CTL_OPT}->{min_contig} && !$shared){
	    $continue_flag = -2; #skip short
	    $self->{die_count} = 0; #reset the die count
	    $rm_key{force}++; #will destroy gff and fastas for this contig
	}
	elsif ($CTL_OPT{force} && !$shared) {
	    $rm_key{force}++;
	    $self->{die_count} = 0; #reset the die count
	    $continue_flag = 1;	#run/re-run
	}
	elsif ($CTL_OPT{again} && !$shared){
	    $rm_key{gff}++;
	    $self->{die_count} = 0; #reset the die count
	    $continue_flag = 1; #run/re-run
	}
	elsif(-e $gff_file){
	    $continue_flag = 0; #don't re-run finished
	}
	elsif ($CTL_OPT{always_try}){
	    $self->{die_count} = 0; #reset the die count
	    $continue_flag = 1; #re-run	    
	}
	elsif($logged_vals{DIED}) {
	    $continue_flag = ($CTL_OPT{clean_try}) ? 2 : 3;	#rerun died
	    $continue_flag = -1 if($self->{die_count} >= $CTL_OPT{tries}); #only let die up to count
	    $continue_flag = -4 if($shared && $continue_flag == -1);
	    $rm_key{force}++ if($continue_flag == 2);	    
	}

        #already belongs to another process
        if($continue_flag == -3 || $continue_flag == 4){
              $self->{continue_flag} = $continue_flag;
              return;
        }

        #set to never try, so just let it init, cleanup, and then stop
        if($CTL_OPT{tries} == 0 && $continue_flag != 0){
             $continue_flag = -1;
        }

	if(!$rm_key{force}){
	    #CHECK CONTROL FILE OPTIONS FOR CHANGES
	    foreach my $key (@ctl_to_log) {
		#ignore hidden key if something else already changed
		if($key eq 'run'){
		    next if(exists $rm_key{gff});
		}

		#these are only sometimes important
		if($key =~ /^est_pass$|^altest_pass$|^protein_pass$|^rm_pass$/ ||
		   $key =~ /^pred_pass$|^model_pass$|^other_pass$/
		   ){
		    next unless($CTL_OPT{maker_gff});
		    my $old = $logged_vals{CTL_OPTIONS}{maker_gff} || '';
		    $old =~ s/(^|\,)$CWD\/*/$1/g;
		    $old = join(',', sort split(',', $old));
		    my $new = $CTL_OPT{maker_gff};
		    $new =~ s/(^|\,)$CWD\/*/$1/g;
		    $new = join(',', sort split(',', $new));

		    #only continue if change not already happening
		    #because of a new gff3 file
		    next unless($old eq $new);
		}

                #these are only sometimes important
                if($key =~ /^map_forward$/){
                    next unless($CTL_OPT{maker_gff} || $CTL_OPT{model_gff});
                }
		
		my $log_val = '';
		if(defined $logged_vals{CTL_OPTIONS}{$key}){	       
		    $log_val = $logged_vals{CTL_OPTIONS}{$key};
		    $log_val =~ s/(^|\,)$CWD\/*/$1/g;
		    $log_val = join(',', sort split(',', $log_val));
		    if($key eq 'repeat_protein'){
			#don't care about absolute location
			$log_val =~ s/[^\,]+\/data\/(te_proteins.fasta)(\,|\:|$)/$1$2/;
		    }			
		}
	    
		my $ctl_val = '';
		if(defined $CTL_OPT{$key}){
		    $ctl_val = $CTL_OPT{$key};
		    $ctl_val =~ s/(^|\,)$CWD\/*/$1/g;
		    $ctl_val = join(',', sort split(',', $ctl_val));
		    if($key eq 'repeat_protein'){
			#don't care about absolute location
			$ctl_val =~ s/[^\,]+\/data\/(te_proteins.fasta)(\,|\:|$)/$1$2/;
		    }
		}

	        #always_complete was off before and not logged
		if($key eq 'always_complete' && ! $log_val){
		    $log_val = 0;
		}

	        #alt_splice was off before and not logged
		if($key eq 'alt_splice' && ! $log_val){
		    $log_val = 0;
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

		#min_intron was always 20 before and not logged
		if($key eq 'min_intron' && $log_val eq ''){
		    $log_val = 20;
		}

		#est_forward was always 0 before and not logged (is is also sometimes hidden)
		if($key eq 'est_forward'){
		    $ctl_val = 0 if($ctl_val eq '');
		    $log_val = 0 if($log_val eq '');
		}

		#correct_est_fusion was always 0 before and not logged
		if($key eq 'correct_est_fusion'){
		    $ctl_val = 0 if($ctl_val eq '');
		    $log_val = 0 if($log_val eq '');
		}

	        #run_evm was off before and not logged
		if($key eq 'run_evm'){
		    $ctl_val = 0 if($ctl_val eq '');
		    $log_val = 0 if($log_val eq '');
		}

	        #est2genome was previously part of predictor
		if($key eq 'est2genome' && $log_val eq ''){
		    $log_val = (grep {!/altest/} grep {/est2genome/} $logged_vals{predictor}) ? 1 : 0;
		}

	        #altest2genome was previously part of predictor
		if($key eq 'altest2genome' && $log_val eq ''){
		    $log_val = (grep {/altest2genome/} $logged_vals{predictor}) ? 1 : 0;
		}

	        #protein2genome was previously part of predictor
		if($key eq 'protein2genome' && $log_val eq ''){
		    $log_val = (grep {/protein2genome/} $logged_vals{predictor}) ? 1 : 0;
		}

		#if previous log options are not the same as current control file options
		if ($log_val ne $ctl_val) {
		    print STDERR "MAKER WARNING: Control file option \'$key\' has changed\n".
			"Old:$log_val\tNew:$ctl_val\n\n" unless($main::qq);
		    
		    $continue_flag = 1; #re-run because ctlopts changed
		    
		    $rm_key{gff}++; #always rebuild gff when some option has changed

		    #find changed list members
		    my ($change, $keep) = compare_comma_lists($ctl_val, $log_val);

		    #certain keys that don't affect preds
		    if($key ne 'evaluate' &&
		       $key ne 'enable_fathom' &&
		       $key ne 'keep_preds' &&
		       $key ne 'other_pass' &&
		       $key ne 'other_gff' &&
		       $key ne 'map_forward' &&
		       $key ne 'est2genome' &&
		       $key ne 'altest2genome' &&
		       $key ne 'protein2genome' &&
		       @$change
		      ){
			$rm_key{preds}++; #almost all changes affect final predictions
		    }
		    
		    if ($key eq 'organism_type') {
			$rm_key{all}++;
		    }
		    
		    if ($key eq 'max_dna_len') {
			$rm_key{all}++;
		    }
		    
		    if ($key eq 'model_org') {
			$rm_key{all}++;
		    }
		    
		    if(	$key eq 'pcov_rm_blastx' ||
			$key eq 'pcid_rm_blastx' ||
			$key eq 'eval_rm_blastx' ||
			$key eq 'bit_rm_blastx' ||
			$key eq 'blast_type' ||
			$key eq 'use_rapsearch'
			) {
			$rm_key{all_but}++;
		    }

		    if ($key eq 'split_hit' ||
			$key eq 'min_intron' ||
			$key eq'ep_score_limit' ||
			$key eq'en_score_limit'
			) {
			$rm_key{e_exonerate}++;
			$rm_key{a_exonerate}++;
			$rm_key{p_exonerate}++;
		    }

		    if ($key eq 'alt_peptide' ||
			$key eq 'eval_blastx' ||
			$key eq 'softmask'
			) {
			$rm_key{blastx}++;
			$rm_key{p_exonerate}++;
		    }

		    if ($key eq 'eval_tblastx' ||
			$key eq 'split_hit' ||
			$key eq 'softmask'
			) {
			$rm_key{tblastx}++;
			$rm_key{a_exonerate}++;
		    }

		    if ($key eq 'eval_blastn' ||
			$key eq 'split_hit'
			) {
			$rm_key{blastn}++;
			$rm_key{e_exonerate}++;
		    }

		    if ($key eq 'rm_gff') {
			$rm_key{all}++;
		    }

		    if ($key eq 'rmlib') {
			$rm_key{all_but_rb}++;
			$skip{'specific.out'} = $keep;
		    }

		    if ($key eq 'repeat_protein') {
			$rm_key{all_but}++;
			$skip{repeatrunner} = $keep;
		    }
		    
		    if ($key eq 'snaphmm') {
			$rm_key{snap}++;
			$skip{snap} = $keep;
		    }
		    
		    if ($key eq 'augustus_species') {
			$rm_key{augustus}++;
			$skip{augustus} = $keep;
		    }
		    
		    if ($key eq 'fgenesh_par_file') {
			$rm_key{fgenesh}++;
			$skip{fgenesh} = $keep;
		    }
		    
		    if ($key eq 'gmhmm') {
			$rm_key{genemark}++;
			$skip{genemark} = $keep;
		    }

		    if ($key eq 'run_evm') {
			$rm_key{evm}++;
			#$skip{evm} = $keep;
		    }

		    if ($key eq 'protein') {
			$rm_key{blastx}++;
			$rm_key{p_exonerate}++ if(@$change);
			$skip{blastx} = $keep;
		    }		    		   

		    if ($key eq 'est') {
			$rm_key{blastn}++;
			$rm_key{e_exonerate}++ if(@$change);
			$skip{blastn} = $keep;
		    }
		    
		    if ($key eq 'altest') {
			$rm_key{tblastx}++;
			$rm_key{a_exonerate}++ if(@$change);
			$skip{tblastx} = $keep;
		    }
		}
	    }
	     
	    #CHECK FOR FILES THAT DID NOT 
	    while (my $key = each %{$logged_vals{STARTED}}) {
		if (! defined $logged_vals{FINISHED}{$key}) {
		    print STDERR "MAKER WARNING: The file $key\n".
			"did not finish on the last run and must be erased\n" unless($main::qq);

		    push(@files, $key);

		    #this next step will get both temp directories and files that
		    #have the incorrect location in the log but are in theVoid
		    $key =~ /([^\/]+)$/;
		    my $rm_f_name = $1;

		    if(! -e $key && -e "$the_void/$rm_f_name"){
			push(@files, "$the_void/$rm_f_name");
		    }
		    elsif(! -e $key){
			my @f = <$the_void/*/$rm_f_name>;			
			push(@files, @f);
		    }

		    $rm_f_name =~ s/\.fasta$//;
		    my @d = (<$the_void/*$rm_f_name*>, <$the_void/*/*$rm_f_name*>);
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
    }
    elsif($self->{params}->{seq_length} < $self->{CTL_OPT}->{min_contig}){
	$continue_flag = -2; #skip short
	$self->{die_count} = 0; #reset the die count
	$rm_key{force}++; #will destroy gff and fastas for this contig
    }

    #-print file type specific warnings and delete files
    if (exists $rm_key{force} || $rm_key{cleantry}) {
	print STDERR "MAKER WARNING: All old files will be erased before continuing\n" unless($main::qq);
	
	#delete everything in the void
	File::Path::rmtree($the_void);
	File::Path::mkpath($the_void);
	
	#remove all other files
	push (@files, $gff_file) if(-e $gff_file);
	push (@files, @{[<$out_base/*.fasta>]});
	push (@dirs, "$out_base/evaluator");
    }
    elsif (exists $rm_key{all}) {
	print STDERR "MAKER WARNING: Changes in control files make re-use of all old data impossible\n".
	    "All old files will be erased before continuing\n" unless($main::qq);
	
	#delete everything in the void
	File::Path::rmtree($the_void);
	File::Path::mkpath($the_void);
	
	#remove all other files	    
	push (@files, $gff_file) if(-e $gff_file);
	push (@files, @{[<$out_base/*.fasta>]});
	push (@dirs, "$out_base/evaluator");
    }
    elsif (exists $rm_key{all_but_rb}) { #all but repbase reports
	print STDERR "MAKER WARNING: Changes in control files make re-use of all but some RepeatMasker data impossible\n".
	    "All old non-RepeatMasker and some RepeatMasker files will be erased before continuing\n" unless($main::qq);
	
	#delete everything in the void
	my @f = <$the_void/*>;
	grep {-f $_ || !/\/\d+$/} @f;
	push(@f, <$the_void/*/*>);
	@f = grep {!/(\.rb\.out|\.rb\.cat|\.rb\.tbl)$/} @f;
	
	foreach my $e (@{$skip{'specific.out'}}){
	    @f = grep {!/\.$e\.specific\.(out|cat|tbl)$/} @f;
	}
	
	#delete files in the void
	foreach my $f (@f){
	    push (@files, $f) if (-f $f);
	    push (@dirs, $f) if (-d $f);
	}
	
	#remove all other files
	push (@files, $gff_file) if(-e $gff_file);
	push (@files, @{[<$out_base/*.fasta>]});
	push (@dirs, "$out_base/evaluator");
    }
    elsif (exists $rm_key{all_but}) {
	print STDERR "MAKER WARNING: Changes in control files make re-use of all but RepeatMasker data impossible\n".
	    "All old non-RepeatMasker files will be erased before continuing\n" unless($main::qq);
	
	#delete everything in the void
	my @f = <$the_void/*>;
	grep {-f $_ || !/\/\d+$/} @f;
	push(@f, <$the_void/*/*>);
	@f = grep(!/(\.out|\.cat|\.tbl)$/, @f);
	
	#delete files in the void
	foreach my $f (@f){
	    push (@files, $f) if (-f $f);
	    push (@dirs, $f) if (-d $f);
	}
	
	#remove all other files	    
	push (@files, $gff_file) if(-e $gff_file);
	push (@files, @{[<$out_base/*.fasta>]});
	push (@dirs, "$out_base/evaluator");
    }
    else {
	if (exists $rm_key{preds}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of hint based predictions impossible\n".
		"Old hint based prediction files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.auto_annotator.*>, <$the_void/*/*.auto_annotator.*>);
	    push (@files, @f);
	}
	if (exists $rm_key{snap}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old SNAP data impossible\n".
		"Old SNAP files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.snap>, <$the_void/*/*.snap>);
	    foreach my $e (@{$skip{snap}}){
		@f = grep {!/\.$e\.snap$/} @f;
	    }
	    push (@files, @f);
	}
	if (exists $rm_key{augustus}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old Augustus data impossible\n".
		"Old Augustus files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.augustus>, <$the_void/*/*.augustus>);
	    foreach my $e (@{$skip{augustus}}){
		@f = grep {!/\.$e\.augustus$/} @f;
	    }
	    push (@files, @f); 
	}
	if (exists $rm_key{fgenesh}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old FGENESH data impossible\n".
		"Old FGENESH files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.fgenesh>, <$the_void/*/*.fgenesh>);
	    foreach my $e (@{$skip{fgenesh}}){
		@f = grep {!/\.$e\.fgenesh$/} @f;
	    }
	    push (@files, @f);
	}
	if (exists $rm_key{genemark}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old GeneMark data impossible\n".
		"Old GeneMark files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.genemark>, <$the_void/*/*.genemark>);
	    foreach my $e (@{$skip{genemark}}){
		@f = grep {!/\.$e\.genemark$/} @f;
	    }
	    push (@files, @f);
	}
	if (exists $rm_key{evm}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old EVM data impossible\n".
		"Old EVM files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.evm>, <$the_void/*/*.evm>);
	    #foreach my $e (@{$skip{evm}}){
	    #	@f = grep {!/\.$e\.evm$/} @f;
	    #}
	    push (@files, @f);
	}
	if (exists $rm_key{blastn}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of all old EST Blastn data impossible\n".
		"Old EST Blastn files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.blastn>, <$the_void/*/*.blastn>);
	    foreach my $e (@{$skip{blastn}}){
		@f = grep {!/\.$e\.blastn$/} @f;
	    }
	    foreach my $f (@f){
		push (@files, $f) if (-f $f);
		push (@dirs, $f) if (-d $f);
	    }
	}	    
	if (exists $rm_key{tblastx}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old tBlastx data impossible\n".
		"Old tBlastx files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.tblastx>, <$the_void/*/*.tblastx>);
	    foreach my $e (@{$skip{tblastx}}){
		@f = grep {!/\.$e\.tblastx$/} @f;
	    }
	    foreach my $f (@f){
		push (@files, $f) if (-f $f);
		push (@dirs, $f) if (-d $f);
	    }
	}
	if (exists $rm_key{blastx}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old Blastx data impossible\n".
		"Old Blastx files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.blastx>, <$the_void/*/*.blastx>);
	    
	    foreach my $e (@{$skip{blastx}}){
		@f = grep {!/\.$e\.blastx$/} @f;
	    }
	    foreach my $f (@f){
		push (@files, $f) if (-f $f);
		push (@dirs, $f) if (-d $f);
	    }
	    
	}
	
	if (exists $rm_key{repeatrunner}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old repeatrunner data impossible\n".
		"Old repeatrunner files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.repeatrunner>, <$the_void/*/*.repeatrunner>);
	    
	    foreach my $e (@{$skip{repeatrunner}}){
		@f = grep {!/\.$e\.repeatrunner$/} @f;
	    }
	    foreach my $f (@f){
		push (@files, $f) if (-f $f);
		push (@dirs, $f) if (-d $f);
	    }
	    
	}
	
	if (exists $rm_key{e_exonerate}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old EST Exonerate data impossible\n".
		"Old EST Exonerate files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.est_exonerate>, <$the_void/*/*.est_exonerate>);
	    push (@files, @f);
	}
	
	if (exists $rm_key{a_exonerate}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old altEST Exonerate data impossible\n".
		"Old altEST Exonerate files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.alt_exonerate>, <$the_void/*/*.alt_exonerate>);
	    push (@files, @f);
	}
	
	if (exists $rm_key{p_exonerate}) {
	    print STDERR "MAKER WARNING: Changes in control files make re-use of old protein Exonerate data impossible\n".
		"Old protein Exonerate files will be erased before continuing\n" unless($main::qq);
	    
	    my @f = (<$the_void/*.p_exonerate>, <$the_void/*/*.p_exonerate>);
	    push (@files, @f);
	}
	
	if (exists $rm_key{gff}) {
	    print STDERR "MAKER WARNING: Any preexisting GFF3 and fasta files for this contig must now be removed.\n"
		unless($main::qq);
	    push (@files, $gff_file);
	    push (@files, @{[<$out_base/*.fasta>]});
	    push (@files, @{[<$the_void/*.holdover>]});
	    push (@files, @{[<$the_void/*/*.holdover>]});
	    push (@files, @{[<$the_void/*.section>]});
	    push (@files, @{[<$the_void/*/*.section>]});
	    push (@files, @{[<$the_void/evidence_*.gff>]});
	    push (@files, @{[<$the_void/*/evidence_*.gff>]});
	    push (@files, "$the_void/query.masked.fasta");
	    push (@files, "$the_void/*/query.masked.fasta");
	    push (@files, "$the_void/query.masked.gff");
	    push (@files, "$the_void/*/query.masked.gff");
	    push (@dirs, "$out_base/evaluator");
	}
	else{
	    #see if died fasta exists
	    my $d_fasta = $gff_file;
	    $d_fasta =~ s/\.gff$/.died.fasta/;
	    if(-e $d_fasta){
		push (@files, $d_fasta);
	    }
	}
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
    my @d = (<$the_void/*.temp_dir>, <$the_void/*/*.temp_dir>);
    foreach my $d (@d){
	File::Path::rmtree($d) if (-d $d);
    }    

    #just let it run init/cleanup and then stop
    if($CTL_OPT{tries} == 0 && $continue_flag != 0){
        $continue_flag = -1;
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
   
   open (my $LOG, "> $log_file");
   
   #log control file options
   print $LOG "SHARED_ID\t$CTL_OPT{_shared_id}\t\n";
   foreach my $key (@ctl_to_log) {
       my $ctl_val = '';
       if(defined $CTL_OPT{$key}){
	   $ctl_val = $CTL_OPT{$key};
	   $ctl_val =~ s/(^|\,)$CWD\/*/$1/g;
	   $ctl_val = join(',', sort split(',', $ctl_val));
	   if($key eq 'repeat_protein'){
	       #don't care about absolute location
	       $ctl_val =~ s/[^\,]+\/data\/(te_proteins.fasta)(\,|\:|$)/$1$2/;
	   }
       }	  
       print $LOG "CTL_OPTIONS\t$key\t$ctl_val\n";
   }
   print $LOG join("\t", "DIED","COUNT",$self->get_die_count)."\n"
       if($self->{continue_flag} == -1 && $self->get_die_count);
   close($LOG);      
}
#-------------------------------------------------------------------------------
#sets a log child to hold values seperate from the whole log
sub set_child {
    my $self = shift;
    my $tag = shift;

    if(! defined $tag){
	$self->{logchild} = undef;
	return;
    }

    my $log_file = $self->{file_name};
    my $file = $log_file;
    my $the_void = $self->{params}->{the_void};
    $file =~ s/(.*\/)?([^\/]+)$/$the_void\/$2.child.$tag/;
    
    $self->{logchild} = $file;

    #make note of log childs existance
    open(my $LOG, ">> $log_file");
    #add more than once just incase (important when multiple nodes write here)
    print $LOG "LOGCHILD\t$file\t\n".
	       "LOGCHILD\t$file\t\n".
	       "LOGCHILD\t$file\t\n";
    close($LOG);

    #instantiate as empty
    open(my $CHILD, ">$file");
    close($CHILD);
}
#-------------------------------------------------------------------------------
sub add_entry {
   my $self = shift;

   my $type  = shift;
   my $key   = shift;
   my $value = shift;

   $self->{CWD} = Cwd::getcwd() if(!$self->{CWD});
   my $CWD = $self->{CWD};

   my $log_file = ($self->{logchild}) ? $self->{logchild} : $self->{file_name};

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

   return $self->{die_count} || 0;
}
#-------------------------------------------------------------------------------
sub report_status {
   my $self = shift;
   my $flag = $self->{continue_flag};
   my $die_count = $self->get_die_count;
   my $seq_id = $self->{params}->{seq_id};
   my $seq_out_name = Fasta::seqID2SafeID($seq_id);
   my $out_dir = $self->{params}->{out_dir};
   my $the_void = $self->{params}->{the_void};
   my $length = $self->{params}->{seq_length};
   my $CTL_OPT = $self->{CTL_OPT};
   my $tries = $die_count + 1;

   #maintain lock only if status id positive (possible race condition?)
   if($flag > 0){
       $flag = -3 if(!$self->{LOCK}->maintain(30));
   }

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
		   "Tries: $tries!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == 3){
      print STDERR "#---------------------------------------------------------------------\n",
                   "Now retrying the contig!!\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "Tries: $tries!!\n",
                   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);
   }
   elsif($flag == -1){
      print STDERR "#---------------------------------------------------------------------\n",
                   "The contig failed $die_count times!!\n",
		   "Maker will not try again!!\n",
		   "The contig will be stored in a fasta file that you can use for debugging.\n",
		   "SeqID: $seq_id\n",
                   "Length: $length\n",
		   "FASTA: $out_dir/$seq_out_name.died.fasta\n",
		   "#---------------------------------------------------------------------\n\n\n"
		       unless($main::qq);

      my $g_index = GI::build_fasta_index($CTL_OPT->{_g_db});
      my $q_seq_obj = $g_index->get_Seq_by_id($seq_id);

      #still no sequence? try rebuilding the index and try again
      if (! $q_seq_obj) {
	  print STDERR "WARNING: Cannot find >$seq_id, trying to re-index the fasta.\n";
	  $g_index->reindex();
	  $q_seq_obj = $g_index->get_Seq_by_id($seq_id);
	  if (! $q_seq_obj) {
	      print STDERR "stop here: $seq_id\n";
	      confess "ERROR: Fasta index error\n";
	  }
      }
      
      open (my $DFAS, "> $out_dir/$seq_out_name.died.fasta");
      print $DFAS ${&Fasta::seq2fastaRef($seq_id, \ ($q_seq_obj->seq))};
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
   elsif($flag == -4){
       #do nothing
   }
   else{
      confess "ERROR: No valid continue flag\n";
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
      confess "ERROR: No valid continue flag\n";
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
sub fasta_file {
   my $self = shift;
   my $arg = shift;

   if($arg){
       $self->{fasta_file} = $arg;
   }

   return $self->{fasta_file};
}
#-------------------------------------------------------------------------------
#this method is depricated
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
1;
