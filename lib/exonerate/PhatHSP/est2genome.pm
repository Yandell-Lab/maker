package exonerate::PhatHSP::est2genome;

use strict;
use warnings;

use Bio::Search::HSP::PhatHSP::Base;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
	    Bio::Search::HSP::PhatHSP::Base
	   );
}

sub new
 {
        my $class  = shift;
	my @args   = @_;

	my $self = $class->SUPER::new(@args);

	bless ($self, $class);


	return $self;
}
#------------------------------------------------------------------------------
sub strand {
        my $self = shift;
        my $w    = shift;

        $w = 'hit' if $w eq 'sbjct';

        return $self->{_strand_hack}->{$w};
}
#------------------------------------------------------------------------------
sub donor {
	my $self = shift;
	my $d    = shift;

	if (defined($d)){
		$self->{donor} = $d;
	}
	else {
		return $self->{donor};
	}
}
#------------------------------------------------------------------------------
sub acceptor {
        my $self = shift;
        my $a    = shift;

        if (defined($a)){
                $self->{acceptor} = $a;
        }
        else {
                return $self->{acceptor};
        }
}
#------------------------------------------------------------------------------
sub whatIsThere  {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

	unless (defined($pos)){
		#print STDERR "no postition defined on $w\n";
		#print STDERR "in sub whatIsThere\n";
		return undef;
	}
	my @a = $w eq 'query' ?
	  split('',$self->query_string()) :  split('',$self->hit_string());
 
	my $i = $self->nB($w);

	my $delta = $self->strand($w) == 1 ? 1 : -1; 

        while (my $nuc = shift(@a)){
                if ($nuc ne '-'){
                        return $nuc if $i == $pos;
                        $i+= $delta;
                }
                else {
                        #
                }
        }

        return undef;
}

sub show
 {
	my $self = shift;

	my $d = defined($self->donor) ? $self->donor : 'UNDEF';
	my $a = defined($self->acceptor) ? $self->acceptor : 'UNDEF';

        print "--------------------------------------\n";
	print ref($self), "\n";
	print "   ".$self->hit->seq_id()."\n";
        print "--------------------------------------\n";
	print $self->query->seq_id()."\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDt:".$self->frac_identical('total')."\n";
	print "FRAME:".$self->frame()."\n";
	print "DONOR:".$d."\n";
	print "ACCEP:".$a."\n";
}

sub name
 {
	my $self = shift;
	
	return $self->hit->seq_id();
}
sub nB
 {
        my $self = shift;
        my $w    = shift;

        die unless defined($w);

        if ($self->strand($w) == 1){
                return $self->start($w);
        }
        elsif ($w eq 'query' && $self->strand($w) == 0){
                return $self->start($w);
        }
        else {
                return $self->end($w);
        }
}

sub nE
 {
        my $self = shift;
        my $w    = shift;

        if ($self->strand($w) == 1){
                return $self->end($w);
        }
        elsif ($w eq 'query' && $self->strand($w) == 0){
                return $self->end($w);
        }

        else {
                return $self->start($w);
        }
}

sub _check
 {
	my $self = shift;
	my $q_i  = shift;
	my $h_i  = shift;
	my $w    = shift;
	my $pos  = shift;
	my $m    = shift;
	my $q    = shift;
	my $h    = shift;

	if (!defined($pos)){
		#print STDERR "position not defined on $w\n";
		return undef;
	}
	elsif($w eq 'query'){
		return ($q_i, $h_i, $q, $h, $m) if $q_i == $pos;
	}
	elsif ($w eq 'hit') {
		return ($q_i, $h_i, $q, $h, $m) if $h_i == $pos;
	}
	else {
	}

	return undef;

}
sub _set_i
 {
	my $self = shift;
	my $ms   = shift;
	my $qs   = shift;
	my $hs   = shift;


	my ($m, $m_i, $q, $q_i, $h, $h_i);

	$m   = $ms->[0];
        $m_i = $ms->[1];

       	$q   = $qs->[0];
       	$q_i = $qs->[1];

       	$h   = $hs->[0];
       	$h_i = $hs->[1];

        unless (defined($q_i)){
                $q_i = $self->nB('query') unless defined($q_i);
                $h_i = $self->nB('hit')   unless defined($h_i);
                $m_i = 0;

                return ($m_i, $q_i, $h_i);
        }

	if ($q ne '-' && $h ne '-'){
               if ($self->strand('hit') == -1) {$h_i-=1} else{$h_i+=1};
               if ($self->strand('query') == -1){$q_i-=1} else{$q_i+=1};
               $m_i++;
	}
        elsif ($q eq '-' && $h ne '-'){
		if ($self->strand('hit') == -1){$h_i-=1} else{$h_i+=1};
                $m_i++;
        }
        elsif ($q ne '-' && $h eq '-'){
		if ($self->strand('query') == -1){$q_i-=1} else{$q_i+=1};
                $m_i++;
        }


	return $self->_set_exit($m_i, $q_i, $h_i);
}

sub _set_exit
 {
	my $self = shift;
	my $m_i  = shift;
	my $q_i  = shift;
	my $h_i  = shift;

	if ($self->nE('query') == $q_i){
		if ($self->strand('hit') == 0){
			#$h_i -=2;
		}
		else {
			#$h_i +=2;
		}
	}
	return ($m_i, $q_i, $h_i);
}
################################################## subroutine header end ##

sub AUTOLOAD
 {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "exonerate::PhatHsp::est2genome::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

