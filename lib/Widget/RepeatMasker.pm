#------------------------------------------------------------------------
#----                        Widget::RepeatMasker                    ---- 
#------------------------------------------------------------------------
package Widget::RepeatMasker;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
use Bio::Search::Hit::PhatHit::repeatmasker;
use Bio::Search::HSP::PhatHSP::repeatmasker;
use IPC::Open3;
use Symbol;

@ISA = qw(
	Widget
       );

#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $class  = shift;
        my @args   = @_;

        my $self = $class->SUPER::new(@args);

	bless ($self, $class);
        return $self;
}
#------------------------------------------------------------------------------
sub run {
	my $self    = shift;
	my $command = shift;

	my $exe       = '/usr/local/bdgp/RepeatMasker/RepeatMasker';
	my $lib       = 'bdgp.dros.lib';
	my $projDir   = $self->projectDir();
	my $fastaFile = $self->queryFastaFile();

	if (defined($command)){
		$self->print_command($command);
		my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
		my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);

		my $err;
		{
		    local $/ = \1;
		    while (my $line = <$CHLD_ERR>){
			print STDERR $line unless($main::quiet);
			$err .= $line;
			kill (9, $pid) if ($err =~ /refinelib\) does not exist/); #kill rather than wait for failure
		    }
		}
		waitpid $pid, 0;
		
		#try again because of error caused by RM algorithm
		#I submitted a bug fix, they ignored it, so this is how I get arround it
		if($? != 0 && $err =~ /refinelib\) does not exist/){
		    print STDERR "Reconfiguring command and trying again.\n" unless($main::quiet);
		    $command =~ s/\-species\s+[^\s]+/-species mammalia/;
		    $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
		    $err = ();
		    while (my $line = <$CHLD_ERR>){
			print STDERR $line unless($main::quiet);
			$err.= $line;
		    }
		    waitpid $pid, 0;
		    die "ERROR: RepeatMasker failed\n" if $? != 0;
		}
		elsif($? != 0){
		    die "ERROR: RepeatMasker failed\n";
		}
	}
	else {
		$self->print_command();
		system("$exe $fastaFile -lib $lib -dir $projDir");
		my $name = $self->queryName();
		$self->maskedFastaFile("$projDir/$name\.masked");
	}

}
#-------------------------------------------------------------------------------
sub maskedFastaFile {
	my $self = shift;
	my $name = shift;

	if    (defined($name)){
		$self->{maskedFastaFile} = $name;
	}
	elsif (defined($self->{maskedFastaFile})){
		return $self->{maskedFastaFile};
	}
	else {
		my $name = $self->queryName();
		my $projDir = $self->projectDir();
		my $masked = "$projDir/$name\.masked";

		$self->{maskedFastaFile} = $masked;

		return $self->{maskedFastaFile};
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_fields {

	my $line = shift;

	my @fields = split(/\s+/, $line);


        if ($fields[10] =~ /[A-Za-z\-]/){
		return @fields;
        }
        elsif ($fields[10] =~ /\(?\d+\)?$/) {
                # file has missing field!!
                my @pre = splice(@fields, 0, 10);

                my @new_fields = (@pre, 'missing_from_file', @fields);

                return @new_fields;

        }
	else {
		die " unknown event in widget::repeatmasker::get_fields line:$line\n";
	}

}
#-------------------------------------------------------------------------------
sub parse {
	my $file     = shift;
	my $q_name   = shift;
	my $q_length = shift;

	my $fh = new FileHandle();
	   $fh->open("$file") || die "ERROR: Could not open \'$file\'\n";

	my %hsps;
	my $count; #used to see if file finished
	while (my $line = <$fh>){
	        $count++;
		next if $line =~ /There were no repetitive sequences detected in/;
		next unless $line =~ /^\s*\d/;

		$line =~ s/^\s+//;

		#my @fields = split(/\s+/, $line);
		my @fields = get_fields($line);

		my $q_begin = $fields[5];
		my $q_end   = $fields[6];

		my ($q_remainder) = $fields[7] =~/\((\d+)\)/;
		my $q_seq_length  = $fields[6] + $q_remainder; 	


		my $h_begin;
		my $h_end;
		my $h_remainder;
		my $h_seq_length;
		if ($fields[8] eq 'C'){
			$h_end   = abs $fields[12];
			$h_begin = abs $fields[13];
			
			($h_remainder) = $fields[11] =~/\((\d+)\)/;
		}
		else {
			$h_begin = $fields[11];
			$h_end   = $fields[12];
			($h_remainder) = $fields[13] =~/\((\d+)\)/;
		}
		$h_begin = 1 if($h_begin == 0); #only happens when RepeatMasker configured with with rmblast

		$h_seq_length = $h_end + $h_remainder;		
		my $identical = (1 - $fields[1]/100) * abs($q_end - $q_begin);
		my $q_gaps    = (1 - $fields[2]/100) * abs($q_end - $q_begin);
		my $h_gaps    = (1 - $fields[3]/100) * abs($h_end - $h_begin);
		my $h_name = "species:$fields[9]|genus:$fields[10]";
                my @args;

                push(@args, '-query_start');
                push(@args, $q_begin);

                #push(@args, '-query_seq');
                #push(@args, $fields[4]);

                push(@args, '-score');
                push(@args, $fields[0]);

                #push(@args, '-homology_seq');
                #push(@args, $m);

                push(@args, '-hit_start');
                push(@args, $h_begin);

                #push(@args, '-hit_seq');
                #push(@args, $fields[10]);

                push(@args, '-hsp_length');
                push(@args, $q_end - $q_begin);

                push(@args, '-identical');
                push(@args, $identical);

                push(@args, '-hit_length');
                push(@args, $h_seq_length);

                push(@args, '-query_name');
                push(@args, $q_name);

                push(@args, '-algorithm');
                push(@args, 'repeatmasker');

                push(@args, '-bits');
                push(@args, 'NA');

                push(@args, '-evalue');
                push(@args, 'NA');

                push(@args, '-pvalue');
                push(@args, 'NA');

                push(@args, '-query_length');
                push(@args, $q_seq_length);

                push(@args, '-query_end');
                push(@args, $q_end);

                push(@args, '-conserved');
                push(@args, $identical);

                push(@args, '-hit_name');
                push(@args, $h_name);

                push(@args, '-hit_end');
                push(@args, $h_end);

                push(@args, '-query_gaps');
                push(@args, $q_gaps);

                push(@args, '-hit_gaps');
                push(@args, $h_gaps);

		my $hsp = Bio::Search::HSP::PhatHSP::repeatmasker->new(@args);
                   $hsp->queryName($q_name);
                #-------------------------------------------------
                # setting strand because bioperl is all fucked up!
                #------------------------------------------------
                if ($q_begin < $q_end){
                        $hsp->{_strand_hack}->{query} = 1;
                }
                else {
                        $hsp->{_strand_hack}->{query} = -1;
                }
                if ($h_begin < $h_end){
                        $hsp->{_strand_hack}->{hit} = 1;
                }
                else {
                        $hsp->{_strand_hack}->{hit} = -1;
                }
                #-------------------------------------------------
                # I hate bioperl!
                #-------------------------------------------------

		push(@{$hsps{$h_name}}, $hsp);

	}
	$fh->close();
	
	#checks if RepeatMasker really finished
	unless($count){
            #unlink($file);
            #die "ERROR: The file $file appears to be incomplete\n".
            #    "MAKER will need to delete the file, before trying again\n\n";
        }

	my @keepers;
	foreach my $key (keys %hsps){
		 my $f =
		     Bio::Search::Hit::PhatHit::repeatmasker->new('-name' => $key,
								  '-description'  => 'NA',
								  '-algorithm'    => 'repeatmasker',
								  '-length'       => $q_length,
								  );

		$f->queryLength($q_length);
		foreach my $hsp (@{$hsps{$key}}){
			$f->add_hsp($hsp);	
		}
		push(@keepers, $f);
	}
	return \@keepers;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::RepeatMasker::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


