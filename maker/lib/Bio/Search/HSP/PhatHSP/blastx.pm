###################################################### main header begin ##

=head1 Bio::Search::HSP::PhatHSP::blastx

Bio::Search::HSP::PhatHSP::blastx - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastx", "check type.");

=head1 DESCRIPTION

Stub documentation for this module was created by
ExtUtils::ModuleMaker.  And, then it was poked, prodded, and otherwise
massaged into it's current form by George.

Hopefully the module author wasn't negligent enough to leave the stub
unedited.

Blah blah blah.

=head1 USAGE

Expand on the examples from the SYNOPSIS.

=head1 BUGS

Not yet.

=head1 AUTHOR

 Mark Yandell
 myandell@fruitfly.org
 http://www.yandell-lab.org

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

List other relevant resources.

=head1 FUNCTION/METHOD DOCUMENTATION

The rest of this document describes the {class,package}'s methods and
subroutines.  Private methods are usually preceded with an underscore
(_) and should not be considered part of the supported interface (they
may change without warning).

=cut

###################################################### main header end   ##

package Bio::Search::HSP::PhatHSP::blastx;

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

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastx", "check type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub new
 {
        my $class  = shift;
	my @args   = @_;

	my $self = $class->SUPER::new(@args);

	bless ($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 isContained

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_before = $hsp->isContained('query', 8000);
 my $q_in = $hsp->isContained('query', 8400);
 my $q_after = $hsp->isContained('query', 9125);

 my $h_before = $hsp->isContained('hit', 1300);
 my $h_in = $hsp->isContained('hit', 1350);
 my $h_after = $hsp->isContained('hit', 1610);

=for example end

=for example_testing
  is($q_before, 0, "Check for position before the query sequence.");
  is($q_in, 1, "Check for position in the query sequence.");
  is($q_after, 0, "Check for position after the query sequence.");
  is($h_before, 0, "Check for position before the hit sequence.");
  is($h_in, 1, "Check for position in the hit sequence.");
  is($h_after, 0, "Check for position after the hit sequence.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub isContained {
  my $self = shift;
  my $w    = shift;
  my $pos  = shift;

  if ($w eq 'query'){
    # query strand == 1
    if ($self->strand('query') == 1){
      return 1 if ($pos >= $self->nB('query') && $pos <= $self->nE('query'));
    }
    # query strand == 0 or -1
    else {
      return 1 if ($pos <= $self->nB('query') && $pos >= $self->nE('query'));
    }
  }
  else {
    # hit strand == 0
    if ($self->strand('hit') == 0){
      return 1 if ($pos >= $self->nB('hit') && $pos <= $self->nE('hit'));
    }
    # hit strand == 1 or -1 falls out the bottom...
  }
  return 0;
}

################################################ subroutine header begin ##

=head2 whatIsThere

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_char = $hsp->whatIsThere('query', 9118); # the first match starts here
 my $h_char = $hsp->whatIsThere('hit', 1328); # ditto in subject

 my $q_gap = $hsp->whatIsThere('query', 8940); # there's a query gap here
 my $h_gap = $hsp->whatIsThere('hit', 1394);

=for example end

=for example_testing
  is($q_char, "M", "check base at postition 9118 in query.");
  is($h_char, "M", "check base at postition 1328 in hit.");
  is($q_gap, "K", "check for the gap at pos 8937 in query.");
  is($h_gap, "N", "check base at postition 1394 in hit.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub whatIsThere
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

        if ($w eq 'hit'){
                my @h = split('',$self->hit_string());

                my $i = $self->nB($w);

                while (my $h_aa = shift(@h)){
                        if ($h_aa ne '-'){
                                return $h_aa if $i == $pos;
                                $i++;
                        }
                        else {
                        }
                }
		return undef;
        }
        else {
                my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
		#print " XXX w:$w pos:$pos q_i:$q_i h_i:$h_i q:$q m:$m h:$h \n";
                return undef unless defined($q_i);
                if ($q eq '-'){
                        return undef;
                }
                else {
                        return $q;
                }
        }

}

################################################ subroutine header begin ##

=head2 show

 Usage     : How to use this function/method

=for example begin

  use PROTO;
  my $foo = new PROTO;

=for example end

=for example_testing
  isa_ok($foo, "PROTO", "Check if it's the right type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub show
 {
	my $self = shift;

        print "--------------------------------------\n";
	print ref($self), "\n";
# XXXX change me to self->name()...	
	print "   ".$self->hit->seqname()."\n";
        print "--------------------------------------\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDt:".$self->frac_identical('total')."\n";
	print "FRAME:".$self->frame()."\n";
}

################################################ subroutine header begin ##

=head2 name

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $name = $hsp->name();

=for example end

=for example_testing
  is($name, "roo|orf1", "check type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub name
 {
	my $self = shift;
	
	return $self->hit->seq_id();
}

################################################ subroutine header begin ##

=head2 nB

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_begin = $hsp->nB('query');
 my $h_natural_begin = $hsp->nB('hit');

=for example end

=for example_testing
  is($q_natural_begin, 9118, "Check the query's natural begin.");
  is($h_natural_begin, 1328, "Check the hit's natural begin.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nB
 {
	my $self = shift;
	my $w    = shift;

	die unless defined($w);

	if ($self->strand($w) == 1){
		return $self->start($w);
	}
        if ($self->strand($w) == -1){
                return $self->end($w);
        }
	elsif ($w eq 'hit' && $self->strand($w) == 0){
		return $self->start($w);
	}
	else {
		die "dead in blastx::PhatHsp::nB\n";
	}
}

################################################ subroutine header begin ##

=head2 nE

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastx',
              'sample_data/blastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_end = $hsp->nE('query');
 my $h_natural_end = $hsp->nE('hit');

=for example end

=for example_testing
  is($q_natural_end, 8354, "Check the query's natural end.");
  is($h_natural_end, 1583, "Check the hit's natural end.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub nE
 {
        my $self = shift;
        my $w    = shift;

        if ($self->strand($w) == 1){
                return $self->end($w);
        }
        if ($self->strand($w) == -1){
                return $self->start($w);
        }
        elsif ($w eq 'hit' && $self->strand($w) == 0){
                return $self->end($w);
        }

        else {
		die "dead in blastx::PhatHsp::nB\n";
        }
}

################################################ subroutine header begin ##

=head2 _check

 Usage     : cigar_string

 Purpose   : produces gff3 compatible cigar sring
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

sub cigar_string {
    my $self = shift;

    return($self->{_CIGAR}) if($self->{_CIGAR});

    my $cigar = '';
    my $q_str = ($self->strand("hit") == -1) ? reverse($self->query_string()) : $self->query_string();
    my $h_str = ($self->strand("hit") == -1) ? reverse($self->hit_string()) : $self->hit_string();

    my @q = $q_str =~ /(.)/g;
    my @h = $h_str =~ /(.)/g;

    #hack to fix a bioperl bug where trailing '-' is not added
    if(@q - @h == 1){
	push(@h, '-');
	$self->hit_string($self->hit_string().'-');
    }
    elsif(@h - @q == 1){
	push(@q, '-');
	$self->query_string($self->query_string().'-');
    }

    die "ERROR: query and hit string lengths do not match correctly\n".
        "in Bio::Search:HSP::PhatHSP::blastx\n".
	"for hit ".$self->name()."\n" if(@q != @h);

    my $cigr = '';
    my $type = ''; # M, I, D
    my $value = 0;

    for(my $i = 0; $i < @h; $i++){
        my $found = '';

        #M
        if($q[$i] ne '-' && $h[$i] ne '-'){
            $found = 'M';
        }
        #I
        elsif($q[$i] eq '-'){
            $found = 'I';
        }
        #D
        elsif($h[$i] eq '-'){
            $found = 'D';
        }

        if($found eq $type){
            $value++;
        }
        else{
            $cigar .= "$type$value" if($type);
            $type = $found;
            $value = 1;
        }
    }

    if($value){
       $cigar .= "$type$value" if($type);
    }

    return ($self->{_CIGAR} = $cigar);
}

################################################ subroutine header begin ##

=head2 _check

 Usage     : *private

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :
           :
 See Also  :

=cut

################################################## subroutine header end ##

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

        if ($w eq 'hit'){
		#print "w:$w pos:$pos q_i:$q_i h_i:$h_i h:$h m:$m q:$q \n";
                return ($q_i, $h_i, $q, $h, $m) if $h_i == $pos;
        }
        elsif ($w eq 'query') {
                #print "ZZZ w:$w pos:$pos q_i:$q_i h_i:$h_i q:$q m:$m h:$h \n";
                if ($pos == $self->nE('query') && $q_i == $pos){
                        if ($self->strand('query') == 1){
                                return ($q_i - 2 , $h_i, $q, $h, $m);
                        }
                        else {
                                return ($q_i , $h_i, $q, $h, $m);
                        }
                }
		elsif ($q_i == $pos){
			return ($q_i, $h_i, $q, $h, $m);
		}
                elsif ($self->strand('query') == -1){
			return undef if $pos < $self->nE('query');
                        my $distance = $q_i - $pos;

			#print "distance:$distance\n";
                        return ($q_i, $h_i, $q, $h, $m)
                        if ($distance > 0 && $distance < 3);
                }
                else {
# XXXX
			return undef if $pos < $self->nB('query');
                        my $distance = $pos - $q_i;
                        return ($q_i, $h_i, $q, $h, $m)
                        if ($distance > 0 && $distance < 3);
                }

        }
        else {
                die "dead in PhatHsp::_check err:2\n";
        }

        return undef;
}


################################################ subroutine header begin ##

=head2 _set_i

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

	if   ($m eq $q && $q eq $h){
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$h_i++;
		$m_i++;
	}
	elsif ($m eq '+'){
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$h_i++;
		$m_i++;
	}
	elsif ($m eq ' ' && $q ne '-' && $h ne '-'){
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$h_i++;
		$m_i++;
	}
        elsif ($m eq ' ' && $q eq '-' && $h ne '-'){
		$h_i++;
                $m_i++;
        }
        elsif ($m eq ' ' && $q ne '-' && $h eq '-'){
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
                $m_i++;
        }


	return $self->_set_exit($m_i, $q_i, $h_i);
}

################################################ subroutine header begin ##

=head2 _set_exit

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

sub _set_exit
 {
	my $self = shift;
	my $m_i  = shift;
	my $q_i  = shift;
	my $h_i  = shift;

	#print "FFF q_i:$q_i h_i:$h_i\n";	
	if ($self->nE('query') == $q_i){
		if ($self->strand('query') == -1){
			$q_i +=2;
		}
		else {
			$q_i -=2;
		}
	}
	return ($m_i, $q_i, $h_i);
}

################################################ subroutine header begin ##

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

        if($ENV{CGL_CHATTER}) {
	    print STDERR "blastx::PhatHsp::AutoLoader called for: ",
	    "\$self->$call","()\n";
	    print STDERR "call to AutoLoader issued from: ", $caller, "\n";
	}

        if (@_){
                $self->{$call} = shift;
        }

	return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

