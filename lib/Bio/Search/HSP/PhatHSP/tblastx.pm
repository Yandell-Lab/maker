###################################################### main header begin ##

=head1 Bio::Search::HSP::PhatHSP::tblastx

Bio::Search::HSP::PhatHSP::tblastx - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::tblastx", "check type.");

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

package Bio::Search::HSP::PhatHSP::tblastx;

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::tblastx", "check type.");

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

=head2 whatIsThere

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_char = $hsp->whatIsThere('query', 869); # the first match starts here
 my $h_char = $hsp->whatIsThere('hit', 70542); # ditto in subject

 my $h_gap = $hsp->whatIsThere('hit', 70559); # ditto in subject

=for example end

=for example_testing
  is($q_char, "S", "check amino acid at postition 869 in query.");
  is($h_char, "S", "check amino acid at postition 70542 in hit.");
  is($h_gap, "E", "check amino acid at postition 70559 in hit.");

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

	my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
        if ($w eq 'query'){
		my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
		return undef unless defined($q_i);
                if ($q eq '-'){
                        return undef;
                }
                else {
                        return $q;
                }

        }
        else {
                my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
                return undef unless defined($h_i);
                if ($h eq '-'){
                        return undef;
                }
                else {
                        return $h;
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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $name = $hsp->name();

=for example end

=for example_testing
  is($name, "Contig190.4:0", "Check the hit's name.");

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
	
	return $self->hit->seqname();
}

################################################ subroutine header begin ##

=head2 nB

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_begin = $hsp->nB('query');
 my $h_natural_begin = $hsp->nB('hit');

=for example end

=for example_testing
  is($q_natural_begin, 869, "Check the beginning of the query.");
  is($h_natural_begin, 70542, "Check the beginning of the hit.");

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
	elsif ($w eq 'query' && $self->strand($w) == 0){
		return $self->start($w);
	}
	else {
		return $self->end($w);
	}
}

################################################ subroutine header begin ##

=head2 nE

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('tblastx',
              'sample_data/tblastx.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_end = $hsp->nE('query');
 my $h_natural_end = $hsp->nE('hit');

=for example end

=for example_testing
  is($q_natural_end, 3412, "Check the end of the query.");
  is($h_natural_end, 73097, "Check the end of the hit.");

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
        elsif ($w eq 'query' && $self->strand($w) == 0){
                return $self->end($w);
        }

        else {
                return $self->start($w);
        }
}

################################################ subroutine header begin ##

=head2 _check

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

	#print "q_i:$q_i h_i:$h_i w:$w pos:$pos m:$m q:$q h:$h\n";
	if ($w eq 'query'){
		if    ($q_i == $pos){
	              return ($q_i, $h_i, $q, $h, $m);
		}
		elsif ($self->strand('query') == 1){
			return undef if $pos < $self->nB('query');
			return undef if $pos > $self->nE('query');
                        if ($pos <=  $self->nE('query')
                        &&  $pos - $q_i < 3
                        &&  $pos - $q_i >= 0){
				 return ($q_i, $h_i, $q, $h, $m);
			}

		}
		elsif ($self->strand('query') == -1){
			return undef if $pos < $self->nE('query');
			return undef if $pos > $self->nB('query');
			if ($self->nE('query') - $pos  < 3
			&&  $h_i - $pos < 3){
				return ($q_i, $h_i, $q, $h, $m);
			}
		}
		else {
			die "dead in PhatHsp::_check err:1\n";	
		}
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

		        my $distance = $h_i - $pos;
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
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$m_i++;
	}
	elsif ($m eq '+'){
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$m_i++;
	}
	elsif ($m eq ' ' && $q ne '-' && $h ne '-'){
		if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
		if ($self->strand('query') == -1){$q_i-=3} else{$q_i+=3};
		$m_i++;
	}
        elsif ($m eq ' ' && $q eq '-' && $h ne '-'){
                if ($self->strand('hit') == -1){$h_i-=3} else{$h_i+=3};
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
 Comments  : handle falling off of the end of HSP (both ends, name is
           : historical (hysterical?)
 See Also  :

=cut

################################################## subroutine header end ##

# XXXX This is the same as tblastn but probably shouldn't be, looks like
# XXXX a cut and pasto.

sub _set_exit {
	my $self = shift;
	my $m_i  = shift;	# position in middle string
	my $q_i  = shift;	# position in query string
	my $h_i  = shift;	# position in hit string

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
	    print STDERR "tblastx::PhatHsp::AutoLoader called for: ",
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

