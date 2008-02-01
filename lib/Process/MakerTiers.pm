#! /usr/bin/perl -w

package Process::MakerTiers;

use FindBin;
use lib "$FindBin::Bin/../..";

use strict;
use Process::MakerChunk;
use File::Path;
use URI::Escape;

#-----------------------------------------------------------------------------
#------------------------------GLOBAL VARIABLES-------------------------------
#-----------------------------------------------------------------------------
my $TIER_ID = 0;
my %INDEX;
my %SEEN;

#-----------------------------------------------------------------------------
#-----------------------------------METHODS-----------------------------------
#-----------------------------------------------------------------------------
sub new {
    my ($class, @args) = @_;

    my $self = {};

    bless ($self, $class);

    $self->{TIER_ID} = $TIER_ID;
    $TIER_ID++;

    if (@args) {
	my $arg = shift @args;
	if (ref $arg eq 'Process::MakerChunk') {
	    $self = $arg->clone();
	}
	else {
	    $self->{VARS}{fasta} = $arg;
	    $self->{VARS}{CTL_OPTIONS} = shift @args;
	    $self->{VARS}{OPT} = shift @args;
	    $self->{VARS}{DS_FH} = shift @args;
	    $self->_initialize();
	}
    }

    $INDEX{$self->{TIER_ID}} = $self;

    return $self;
}

#--------------------------------------------------------------
sub _initialize {
    my $self        = shift;
    my $fasta       = $self->{VARS}{fasta};
    my %CTL_OPTIONS = %{$self->{VARS}{CTL_OPTIONS}};
    my %OPT         = %{$self->{VARS}{OPT}};

    #-------------------------CHUNK
    $self->{VARS}{query_def} = Fasta::getDef($fasta); #Get fasta header
    $self->{VARS}{query_seq} = Fasta::getSeq($fasta); #Get fasta sequence
    ($self->{VARS}{seq_id})  = $self->{VARS}{query_def} =~ /^>(\S+)/; #Get sequence identifier off of fasta header
    
    #--error checking for no-unique ids in multi-fasta
    if (exists $SEEN{$self->{VARS}{seq_id}}){
	warn "ERROR:  The multi-fasta file contains non-unique sequence ids.\n",
	     "The id " . $self->{VARS}{seq_id} . " occurs more than once.\n";

	$SEEN{$self->{VARS}{seq_id}}++;
	my $new_id = $self->{VARS}{seq_id}.".V".$SEEN{$self->{VARS}{seq_id}};

	warn "\n\n\nThe id ".$self->{VARS}{seq_id}." will now be changed to $new_id \n";

	$self->{VARS}{seq_id} = $new_id;
	$SEEN{$new_id}++;
    }
    else{
	$SEEN{$self->{VARS}{seq_id}}++;
    }

    #--build a safe name for file names from the sequence identifier  
    $self->{VARS}{seq_out_name} = uri_escape($self->{VARS}{seq_id},
					    '\*\?\|\\\/\'\"\{\}\<\>\;\,\^\(\)\$\~\:'
					   );

    #--set up void directory where analysis is stored
    $self->{VARS}{out_dir} = $CTL_OPTIONS{'out_base'};

    if (exists $CTL_OPTIONS{'datastore'}) {
	$self->{VARS}{out_dir} = $CTL_OPTIONS{'datastore'}->id_to_dir($self->{VARS}{seq_out_name});
	
	$CTL_OPTIONS{'datastore'}->mkdir($self->{VARS}{seq_out_name}) ||
	    die "ERROR: could not make directory $self->{VARS}{out_dir}\n";
	
	my $fh = $self->{VARS}{DS_FH};
	print $fh "$self->{VARS}{seq_id}\t$self->{VARS}{out_dir}\n";
    }

    $self->{VARS}{the_void}  = build_the_void($self->{VARS}{seq_out_name},
					     $OPT{a}, $self->{VARS}{out_dir}
					    );
    #-------------------------CHUNK
    $self->{TERMINATE} = 0;
    $self->{LEVEL}{CURRENT} = -1;
    $self->next_level();
}

#--------------------------------------------------------------
sub id {
    my $self = shift;
   
    return $self->{TIER_ID};
}

#--------------------------------------------------------------
sub next_chunk {
    my $self = shift;

    my $current = $self->{LEVEL}{CURRENT};

    if (my $chunk = shift @{$self->{LEVEL}{$current}{CHUNKS}}) {
	return $chunk;
    }

    return undef;
}

#--------------------------------------------------------------
sub _build_chunk {
    my $self = shift;
    my $level = shift;
    my $args = shift;

    my $chunk_id = $self->id().":".$level.":". $self->{LEVEL}{$level}{CHUNK_COUNT};
    
    my $chunk = Process::MakerChunk->new($level, $args, $chunk_id);
    push (@{$self->{LEVEL}{$level}{CHUNKS}}, $chunk);
    $self->{LEVEL}{$level}{CHUNK_COUNT}++;
}
#--------------------------------------------------------------
sub polish_results{
    my $self = shift;
    my $level = $self->{LEVEL}{CURRENT};

    unless ($level >= 0 && $self->level_finished()){
	return 0;
    }

    if (not $self->{LEVEL}{$level}{RESULTS}){
	return 0;
    }

    foreach my $result (@{$self->{LEVEL}{$level}{RESULTS}}){
	my @results = @{$result};

	if ($level == 0) {
	    #------------------------RESULTS
	    $self->{VARS}{temp_masked_fasta} = shift @results;
	    $self->{VARS}{rma_keepers}       = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 1) {
	    #------------------------RESULTS
	    push (@{$self->{VARS}{repeat_blastx_keepers}}, @{shift @results});
	    #------------------------RESULTS
	}
	elsif ($level == 2) {
	    #------------------------RESULTS
	    $self->{VARS}{rm_keepers}   = shift @results;
	    $self->{VARS}{masked_fasta} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 3) {
	    #------------------------RESULTS
	    $self->{VARS}{IOX}           = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 4) {
	    #------------------------RESULTS
	    push (@{$self->{VARS}{blastx_keepers}}, @{shift @results});
	    #------------------------RESULTS
	}
	elsif ($level == 5) {
	    #------------------------RESULTS
	    push (@{$self->{VARS}{tblastx_keepers}}, @{shift @results});
	    #------------------------RESULTS
	}
	elsif ($level == 6) {
	    #------------------------RESULTS
	    push (@{$self->{VARS}{blastn_keepers}}, @{shift @results});
	    #------------------------RESULTS
	}
	elsif ($level == 7) {
	    #------------------------RESULTS
	    $self->{VARS}{IOX}             = shift @results;
	    $self->{VARS}{blastx_clusters} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 8) {
	    #------------------------RESULTS
	    $self->{VARS}{exonerate_p_clusters} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 9) {
	    #------------------------RESULTS
	    $self->{VARS}{IOX}              = shift @results;
	    $self->{VARS}{tblastx_clusters} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 10) {
	    #------------------------RESULTS
	    $self->{VARS}{blastn_clusters} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 11) {
	    #------------------------RESULTS
	    $self->{VARS}{exonerate_e_clusters} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 12) {
	    #------------------------RESULTS
	    $self->{VARS}{IOX} = shift @results;
	    $self->{VARS}{snaps} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 13) {
	    #------------------------RESULTS
	    $self->{VARS}{annotations} = shift @results;
	    #------------------------RESULTS
	}
	elsif ($level == 14) {
	    #------------------------RESULTS
	    $self->{TERMINATE} = 1;
	    $INDEX{$self->id()} = undef;
	    #------------------------RESULTS
	}
    }

    return 1;
}
#--------------------------------------------------------------
sub next_level {
    my $self  = shift;
    my $level = $self->{LEVEL}{CURRENT};
    my %OPT   = %{$self->{VARS}{OPT}};
    my @args;
    
    if ($level  >= 0) {
	if ($self->level_finished() && not $self->terminated()){

	    #--get results for current level
	    $self->polish_results();
	}
	else {
	    return 0;
	}
    }	    
    
    #--now go up one level
    $level++;
    $self->{LEVEL}{CURRENT} = $level;
    
    #--initiate level variables
    $self->{LEVEL}{$level}{CHUNK_COUNT} = 0;
    $self->{LEVEL}{$level}{RESULT_COUNT} = 0;
    $self->{LEVEL}{$level}{CHUNKS} = undef;
    $self->{LEVEL}{$level}{RESULTS} = undef;
    $self->{LEVEL}{$level}{ERROR} = undef;
    
    #--select variables to send to Process::MakerChunk object
    if ($level == 0) {
	$self->{VARS}{masked_fasta} = \$self->{VARS}{fasta};
	$self->{VARS}{rm_keepers} = [];

	return $self->next_level() if ($OPT{R});
	
	#------------------------ARGS_IN
	@args =( $self->{VARS}{fasta},
		 $self->{VARS}{the_void},
		 $self->{VARS}{seq_out_name},
		 $self->{VARS}{query_seq},
		 $self->{VARS}{query_def},
		 $self->{VARS}{CTL_OPTIONS},
		 $self->{VARS}{OPT}{f}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 1){
	return $self->next_level() if ($OPT{R});

	foreach my $repeat_protein (@{$self->{VARS}{CTL_OPTIONS}{repeat_protein}}){
	    #------------------------ARGS_IN
	    @args =( $self->{VARS}{temp_masked_fasta},
		     $repeat_protein,
		     $self->{VARS}{the_void},
		     $self->{VARS}{query_seq},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	    #------------------------ARGS_IN
		
	    #-------------------------CHUNK
	    $self->_build_chunk($level,\@args);
	    #-------------------------CHUNK
	}
	
	return 1;
    }
    elsif ($level == 2){
	return $self->next_level() if ($OPT{R});

	#------------------------ARGS_IN
	@args =( $self->{VARS}{query_seq},
		 $self->{VARS}{rma_keepers},
		 $self->{VARS}{repeat_blastx_keepers},
		 $self->{VARS}{query_def},
		 $self->{VARS}{the_void}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 3){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{masked_fasta},
		 $self->{VARS}{the_void},
		 $self->{VARS}{out_dir},
		 $self->{VARS}{seq_out_name},
		 $self->{VARS}{rm_keepers},
		 $self->{VARS}{CTL_OPTIONS}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 4){
	foreach my $protein (@{$self->{VARS}{CTL_OPTIONS}{protein}}){
	    #------------------------ARGS_IN
	    @args =( $self->{VARS}{masked_fasta},
		     $protein,
		     $self->{VARS}{the_void},
		     $self->{VARS}{query_seq},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	    #------------------------ARGS_IN

	    #-------------------------CHUNK
	    $self->_build_chunk($level,\@args);
	    #-------------------------CHUNK
	}
	
	return 1;
    }
    elsif($level == 5){
	foreach my $alt_est (@{$self->{VARS}{CTL_OPTIONS}{alt_est}}){
	    #------------------------ARGS_IN
	    @args =( $self->{VARS}{masked_fasta},
		     $alt_est,
		     $self->{VARS}{the_void},
		     $self->{VARS}{query_seq},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	    #------------------------ARGS_IN

	    #-------------------------CHUNK
	    $self->_build_chunk($level,\@args);
	    #-------------------------CHUNK
	}

	return 1;
    }
    elsif($level == 6){
	foreach my $est (@{$self->{VARS}{CTL_OPTIONS}{est}}){
	    #------------------------ARGS_IN
	    @args =( $self->{VARS}{masked_fasta},
		     $est,
		     $self->{VARS}{the_void},
		     $self->{VARS}{query_seq},
		     $self->{VARS}{seq_out_name},
		     $self->{VARS}{CTL_OPTIONS},
		     $self->{VARS}{OPT}{f}
		   );
	    #------------------------ARGS_IN

	    #-------------------------CHUNK
	    $self->_build_chunk($level,\@args);
	    #-------------------------CHUNK
	}

	return 1;
    }
    elsif($level == 7){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{blastx_keepers},
		 $self->{VARS}{tblastx_keepers},
		 $self->{VARS}{blastn_keepers},
		 $self->{VARS}{IOX},
		 $self->{VARS}{query_seq},
		 $self->{VARS}{the_void}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 8){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{fasta},
		 $self->{VARS}{blastx_clusters},
		 $self->{VARS}{CTL_OPTIONS}{old_protein},
		 $self->{VARS}{the_void},
		 $self->{VARS}{CTL_OPTIONS},
		 $self->{VARS}{OPT}{f}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 9){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{exonerate_p_clusters},
		 $self->{VARS}{IOX},
		 $self->{VARS}{tblastx_keepers},
		 $self->{VARS}{query_seq}
	       );
	#------------------------ARGS_IN

    }
    elsif($level == 10){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{blastn_keepers},
		 $self->{VARS}{query_seq}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 11){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{fasta},
		 $self->{VARS}{blastn_clusters},
		 $self->{VARS}{CTL_OPTIONS}{old_est},
		 $self->{VARS}{the_void},
		 $self->{VARS}{CTL_OPTIONS},
		 $self->{VARS}{OPT}{f}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 12){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{exonerate_e_clusters},
		 $self->{VARS}{IOX},
		 $self->{VARS}{masked_fasta},
		 $self->{VARS}{the_void},
		 $self->{VARS}{seq_out_name},
		 $self->{VARS}{CTL_OPTIONS},
		 $self->{VARS}{OPT}{f}
	       );
      	#------------------------ARGS_IN
    }
    elsif($level == 13){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{fasta},
		 $self->{VARS}{masked_fasta},
		 $self->{VARS}{exonerate_p_clusters},
		 $self->{VARS}{exonerate_e_clusters},
		 $self->{VARS}{blastx_clusters},
		 $self->{VARS}{snaps},
		 $self->{VARS}{the_void},
		 $self->{VARS}{CTL_OPTIONS},
		 $self->{VARS}{OPT}{f}
	       );
	#------------------------ARGS_IN
    }
    elsif($level == 14){
	#------------------------ARGS_IN
	@args =( $self->{VARS}{blastx_clusters},
		 $self->{VARS}{blastn_clusters},
		 $self->{VARS}{tblastx_clusters},
		 $self->{VARS}{exonerate_p_clusters},
		 $self->{VARS}{exonerate_e_clusters},
		 $self->{VARS}{query_seq},
		 $self->{VARS}{seq_id},
		 $self->{VARS}{annotations},
		 $self->{VARS}{rm_keepers},
		 $self->{VARS}{snaps},
		 $self->{VARS}{out_dir},
		 $self->{VARS}{seq_out_name},
		 $self->{VARS}{IOX},
		 $self->{VARS}{the_void},
		 $self->{VARS}{CTL_OPTIONS}
	       );
	#------------------------ARGS_IN
    }
    else{
	return 0;
    }

    #-------------------------CHUNK
    $self->_build_chunk($level,\@args);
    #-------------------------CHUNK

    return 1;
}

#--------------------------------------------------------------
sub clone {
    my $self = shift;
    
    my $clone = dclone($self);

    return $clone;
}

#--------------------------------------------------------------
sub level_finished {
    my $self = shift;
    my $level = $self->{LEVEL}{CURRENT};

    if ($self->{LEVEL}{$level}{CHUNK_COUNT} == $self->{LEVEL}{$level}{RESULT_COUNT}){
	return 1;
    }

    return 0;
}

#--------------------------------------------------------------
sub terminated {
    my $self = shift;

    return $self->{TERMINATE};
}

#--------------------------------------------------------------
sub update_chunk {
    my $self = shift;
    my $chunk = shift;

    my $id = $chunk->id();
    my @result = $chunk->result();

    my ($tier_id, $level_num, $chunk_num) = split (":", $id);
    
    push (@{$INDEX{$tier_id}->{LEVEL}{$level_num}{RESULTS}}, \@result);
    $INDEX{$tier_id}->{LEVEL}{$level_num}{RESULT_COUNT}++;
    
    $INDEX{$tier_id}->{LEVEL}{$level_num}{ERROR} .= $chunk->error();
    print STDERR $chunk->error();
}

#-----------------------------------------------------------------------------
#------------------------------------SUBS-------------------------------------
#-----------------------------------------------------------------------------
sub build_the_void {
    my $seq_id  = shift;
    my $run_id  = shift;
    my $out_dir = shift;

    $out_dir =~ s/\/$//;

    my $vid = "theVoid\.$seq_id\.$run_id";   
    my $the_void = "$out_dir/$vid";
    mkpath ($the_void);

    return $the_void;
}

#-----------------------------------------------------------------------------
1;
