package exonerate::PhatHSP::protein2genome;

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
sub whatIsThere
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

        if ($w eq 'query'){
              my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
                return undef unless defined($q_i);
                if ($q eq '---'){
                        return undef;
                }
                else {
                        return $q;
                }

        }
        else {
                my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
                return undef unless defined($h_i);
                if ($h eq '---'){
                        return undef;
                }
                else {
                        return $h;
                }
        }

}

#------------------------------------------------------------------------------
sub show
 {
	my $self = shift;

        print "--------------------------------------\n";
	print ref($self), "\n";
	print "   ".$self->hit->seq_id()."\n";
        print "--------------------------------------\n";
	print "hit:".$self->query->seq_id()."\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDq:".$self->frac_identical('query')."\n";
	print "FRAME:".$self->frame()."\n";
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

	if ($w eq 'query'){
		return ($q_i, $h_i, $q, $h, $m) if $q_i == $pos;
	}
        elsif ($w eq 'hit') {
                #print "q_i:$q_i h_i:$h_i w:$w pos:$pos m:$m q:$q h:$h\n";
                return ($q_i, $h_i, $q, $h, $m) if $h_i == $pos;

                if ($pos == $self->nE('hit') && $h_i == $pos){
                        if ($self->strand('hit') == 1){
                                return ($q_i, $h_i - 2, $q, $h, $m);
                        }
                        else {
                                return ($q_i, $h_i, $q, $h, $m);
                        }
                }
                elsif ($self->strand('hit') == 1){
                        my $distance = $pos - $h_i;

                        return ($q_i, $h_i, $q, $h, $m)
                        if ($distance > 0 && $distance < 3);
                }
                else {

                        my $distance =  $h_i - $pos;
                        return ($q_i, $h_i, $q, $h, $m)
                          if ($distance > 0 && $distance < 3);
                }
        }
        else {
                die "dead in PhatHsp::_check err:2\n";
        }
	return undef;

}
#------------------------------------------------------------------------------
sub prep {
	my $str = shift;

	my @str = split('', $str);
	
	my @data;
	while (@str){
		push(@data, shift(@str).shift(@str).shift(@str));
	}

	return @data;
}
#------------------------------------------------------------------------------
sub build_scan_inputs {
	my $self = shift;

        my $q_str = $self->query_string();
        my $h_str = $self->hit_string();
	my $m_str = $self->homology_string();

	my @pieces;
        while ($m_str =~ m/([#]+)/g ) {
                my $e = pos($m_str);
                my $b = $e - length($1);
                push(@pieces, {b => $b ,e => $e,l  => length($1), piece => $1});
        }


	my @sorted = sort {$a->{b} <=> $b->{b}} @pieces;
	my @m;
	my @q;
	my @h;

	my $o = 0;
	while (my $p = shift(@sorted)){
		my $l = $p->{b} - $o;
		my $m_substr = substr($m_str, $o, $l);
		my $q_substr = substr($q_str, $o, $l);
		my $h_substr = substr($h_str, $o, $l);

		push(@m, prep($m_substr));
		push(@q, prep($q_substr));
		push(@h, prep($h_substr));

		push(@m, $p->{piece});
		push(@q, $p->{piece});
		push(@h, $p->{piece});

		$o = $p->{b} + $p->{l};	
	}
	push(@m, prep(substr($m_str, $o)));
	push(@q, prep(substr($q_str, $o)));
	push(@h, prep(substr($h_str, $o)));

	die "q, h, and match strings not equal!\n"
        unless ($#q == $#h && $#h == $#m);

	return (\@m, \@q, \@h);
}
#------------------------------------------------------------------------------
sub _scan
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

	my ($m_array, $q_array, $h_array) = $self->build_scan_inputs();

        my ($m_i, $q_i, $h_i);
        while  (my $m = shift(@{$m_array})){

		my $q = shift(@{$q_array});
		my $h = shift(@{$h_array});

                ($m_i, $q_i, $h_i) = $self->_set_i([$m, $m_i],
                                                   [$q, $q_i],
                                                   [$h, $h_i],
                                                  );

                #print "m:$m m_i:$m_i q:$q q_i:$q_i h:$h h_i:$h_i\n";
                #print "q_i:$q_i h_i:$h_i w:$w pos:$pos m:$m q:$q h:$h\n";
                my @data = $self->_check($q_i, $h_i, $w, $pos, $m, $q, $h);
                return @data if defined($data[0]);
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

	if      ($m =~ /#/){
		my $len = length($m);
		if ($self->strand('hit') == -1){$h_i-=$len} else{$h_i+=$len};
		$m_i+= $len;
		
	}
	elsif   ($m =~ /[\|\:\.\!]/){
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
		$q_i++;
		$m_i++;
	}
	elsif ($m =~ /\s+/  && $q ne '<->'  && $h ne '---'){
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
		$q_i++;
		$m_i++;
	}
        elsif ($m =~ /\s+/  && $q eq '<->'  && $h ne '---'){
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
                $m_i++;
        }
        elsif ($m =~ /\s+/ && $q ne '<->' && $h eq '---'){
                $m_i++;
		$q_i++;
        }
	else {
		die  "in _set_i unknown char :$m\n";
	}

	return $self->_set_exit($m_i, $q_i, $h_i);
}

sub _set_exit
 {
	my $self = shift;
	my $m_i  = shift;
	my $q_i  = shift;
	my $h_i  = shift;

	if ($self->nE('query') == $q_i +1  ){
		if ($self->strand('hit') == 1){
			#$h_i -=3;
		}
		else {
			#$h_i +=3;
		}
	}
	return ($m_i, $q_i, $h_i);
}


=head2 AUTOLOAD

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub AUTOLOAD
 {
        my $self = shift;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "blastn::PhatHsp::AutoLoader called for: ",
        #      "\$self->$call","()\n";
        #print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

