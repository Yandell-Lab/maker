###################################################### main header begin ##

=head1 Bio::Search::HSP::PhatHSP::blastp

Bio::Search::HSP::PhatHSP::blastp - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastp", "check type.");

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

package Bio::Search::HSP::PhatHSP::blastp;

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastp", "check type.");

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

=head2 followingPos

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastp", "check type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

# XXXX DEAD CODE?
sub followingPos
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

}

################################################ subroutine header begin ##

=head2 isContained

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $q_before = $hsp->isContained('query', 0);
 my $q_in = $hsp->isContained('query', 25);
 my $q_after = $hsp->isContained('query', 1000);

 my $h_before = $hsp->isContained('hit', 0);
 my $h_in = $hsp->isContained('hit', 53);
 my $h_after = $hsp->isContained('hit', 1000);

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

  # in the query.
  if ($w eq 'query'){
    return 1 if ($pos >= $self->nB('query') && $pos <= $self->nE('query'));
  }
  # in the hit.
  else {
    return 1 if ($pos >= $self->nB('hit') && $pos <= $self->nE('hit'));
  }
  return 0;
}


################################################ subroutine header begin ##

=head2 whatIsThere

 Usage     : How to use this function/method

=for example
  use Bio::Search::HSP::PhatHSP::Base;
  my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs
             ('blastp', 'sample_data/blastp.sample.report');

=for example begin

  my $hsp = $hsps->[0];
  my $q_char = $hsp->whatIsThere('query', 1);
  my $h_char = $hsp->whatIsThere('hit', 1);

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::blastp");
  # these are at the start of the hsp.
  is($hsp->whatIsThere('query', 1), "M",
     "blastp: ask what's at position 1 in the query.");
  is($hsp->whatIsThere('hit', 1), "M",
     "blastp: ask what's at position 1 in the hit.");


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

        if ($w eq 'query'){
                my @q = split('',$self->query_string());

                my $i = $self->nB($w);

                while (my $q_aa = shift(@q)){
                        if ($q_aa ne '-'){
                                return $q_aa if $i == $pos;
                                $i++;
                        }
                        else {
                        }
                }
		return undef;
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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $name = $hsp->name();

=for example end

=for example_testing
  is($name, "CG9522-PA.3", "Check the hit's name.");

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_begin = $hsp->nB('query');
 my $h_natural_begin = $hsp->nB('hit');

=for example end

=for example_testing
  is($q_natural_begin, 1, "Check the query's natural begin.");
  is($h_natural_begin, 1, "Check the hit's natural begin.");

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

	if ($self->strand($w) == 0){
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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastp',
              'sample_data/blastp.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_end = $hsp->nE('query');
 my $h_natural_end = $hsp->nE('hit');

=for example end

=for example_testing
  is($q_natural_end, 623, "Check the query's natural end.");
  is($h_natural_end, 615, "Check the hit's natural end.");

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

        if ($self->strand($w) == 0){
                return $self->end($w);
        }
        else {
		die "dead in blastx::PhatHsp::nB\n";
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

	if ($w eq 'query'){
		return ($q_i, $h_i, $q, $h, $m) if $q_i == $pos;
	}
	elsif ($w eq 'hit') {
		return ($q_i, $h_i, $q, $h, $m) if $h_i == $pos;
	}
	else {
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
		$h_i++;
		$q_i++;
		$m_i++;
	}
	elsif ($m eq '+'){
		$h_i++;
		$q_i++;
		$m_i++;
	}
	elsif ($m eq ' ' && $q ne '-' && $h ne '-'){
		$h_i++;
		$q_i++;
		$m_i++;
	}
        elsif ($m eq ' ' && $q eq '-' && $h ne '-'){
		$h_i++;
                $m_i++;
        }
        elsif ($m eq ' ' && $q ne '-' && $h eq '-'){
                $m_i++;
		$q_i++;
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
	    print STDERR "blastp::PhatHsp::AutoLoader called for: ",
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

