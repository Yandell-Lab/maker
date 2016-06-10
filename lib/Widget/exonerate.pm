#------------------------------------------------------------------------
#----                        Widget::exonerate                       ---- 
#------------------------------------------------------------------------
package Widget::exonerate;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;
use Widget;
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
    my $o_file  = shift;
    
    if (defined($command)){
	$self->print_command($command);
	my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
	my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	{
	    local $/ = \1;
	    while (my $line = <$CHLD_ERR>){
		print STDERR $line unless($main::quiet);
	    }
	}
	waitpid $pid, 0;
	
	#try a second time
	if($? != 0){
	    $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	    {
		local $/ = \1;
		while (my $line = <$CHLD_ERR>){
		    print STDERR $line unless($main::quiet);
		}
	    }
	    waitpid $pid, 0;
	    
	    if($? != 0){
		unlink($o_file) if($o_file && -f $o_file);
		die "ERROR: Exonerate failed\n";
	    }
	}
    }
    else {
	die " Widget::exonerate::run needs a command!\n";
	
    }   
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub get_model_order {
    my $v = shift;

    my $str = '';
    foreach my $o (@{$v->{operations}}){
	$str .= $o->{state};
    }

    my $type;
    if ($str =~ /3I5/ && $str =~ /5I3/){
	die "ERROR: Mixed model in Widget::exonerate!\n";
    }
    elsif ($str =~ /5I3/){
	$type = '5I3';
    }
    elsif ($str =~ /3I5/){
	$type = '3I5';
    }
    
    if(!$type){
	if($str =~ /[A-HJ-Z]53|53[A-HJ-Z]/ && $str =~ /[A-HJ-Z]35|35[A-HJ-Z]/){
	    die "ERROR: Mixed zero intron model in Widget::exonerate!\n";
	}
	elsif($str =~ /[A-HJ-Z]53|53[A-HJ-Z]/){
	    $type = '5I3';
	}
	elsif($str =~ /[A-HJ-Z]35|35[A-HJ-Z]/){
	    $type = '3I5';
	}
	elsif($str =~ /[35]/){
	    die "ERROR: Could not determine model type in Widget::exonerate!\n";
	}
    }
    
    return $type;
}
#------------------------------------------------------------------------
sub get_exon_coors {
    my $v    = shift;
    my $type = shift;
    
    my $pos_q = $v->{q_b};
    my $pos_t = $v->{t_b};
    my $exon = 0;
    
    my @data;
    foreach my $o (@{$v->{operations}}){
        my $st = $o->{state};
        if($st eq 'M' || $st eq 'G' || $st eq 'F' || $st eq 'S' || $st eq 'C'){ #these are part of exon
            #set exon query start
            if ($v->{q_b} < $v->{q_e}){
                $data[$exon]{q}{b} = $pos_q + 1
                    unless defined($data[$exon]{q}{b});
            }
            else {
                $data[$exon]{q}{b} = $pos_q
                    unless defined($data[$exon]{q}{b});
            }
	    
            #set target query start
            if ($v->{t_b} < $v->{t_e}){
                $data[$exon]{t}{b} = $pos_t + 1
                    unless defined($data[$exon]{t}{b});
            }
            else {
                $data[$exon]{t}{b} = $pos_t
                    unless defined($data[$exon]{t}{b});
            }
        }
        elsif ($st eq '5'|| $st eq '3'){ #splice site
	    if(substr($type, 0, 1) eq $st){ #5 goes with 5I3 & 3 with 3I5
		#set exon query end
		if ($v->{q_b} < $v->{q_e}){
		    $data[$exon]{q}{e} = $pos_q;
		}
		else{
		    $data[$exon]{q}{e} = $pos_q + 1;
		}
		
		#set target query end
		if ($v->{t_b} < $v->{t_e}){
		    $data[$exon]{t}{e} = $pos_t;
		}
		else{
		    $data[$exon]{t}{e} = $pos_t + 1;
		}
		$data[$exon]{q}{strand} = $v->{q_strand};
		$data[$exon]{t}{strand} = $v->{t_strand};
		
		#remove 0 length exons from report (exonerate bug)
		if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
		    $data[$exon] = undef;
		    $exon--; #undo iteration
		}
		
		$exon++;
	    }
        }
        elsif($st eq 'I'){ #intron
            #do nothing
        }
        elsif ($st eq 'N'){ #Non-equivalenced region
            die "dead in Widget::exonerate::get_exon_coors:N\n";
        }
        else {
            die "unknown state in Widget::exonerate::get_exon_coors!\n";
        }
	
        #move position for every type
        if ($v->{q_strand} == 1){
            $pos_q += $o->{q};
        }
        else {
            $pos_q -= $o->{q};
        }
	
        if ($v->{t_strand} == 1){
            $pos_t += $o->{t};
        }
        else {
            $pos_t -= $o->{t};
        }
    }
    
    #finish last exon
    if ($v->{q_b} < $v->{q_e}){
        $data[$exon]{q}{e} = $pos_q;
    }
    else{
        $data[$exon]{q}{e} = $pos_q + 1;
    }
    if ($v->{t_b} < $v->{t_e}){
        $data[$exon]{t}{e} = $pos_t;
    }
    else{
        $data[$exon]{t}{e} = $pos_t + 1;
    }
    $data[$exon]{q}{strand} = $v->{q_strand};
    $data[$exon]{t}{strand} = $v->{t_strand};
    
    #remove 0 length exons in report (exonerate bug)
    if(! defined($data[$exon]{q}{b}) || ! defined($data[$exon]{t}{b})){
        delete($data[$exon]);
        $exon--; #undo iteration
    }
    
    return \@data;
}
#------------------------------------------------------------------------
1;


