###################################################### main header begin ##

=head1 Bio::Search::HSP::PhatHSP::qrna

Bio::Search::HSP::PhatHSP::qrna - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::qrna", "check type.");

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

package Bio::Search::HSP::PhatHSP::qrna;

use strict;
use warnings;
use PostData;
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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness

=for example end

=for example_testing
  isa_ok($hsp, "Bio::Search::HSP::PhatHSP::qrna", "check type.");

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_char = $hsp->whatIsThere('query', 23); # the first match starts here
 my $h_char = $hsp->whatIsThere('hit', 16827); # ditto in subject

 my $q_gap = $hsp->whatIsThere('query', 38); # there's a query gap here
 my $h_gap = $hsp->whatIsThere('hit', 16847);

=for example end

=for example_testing
  is($q_char, "A", "check base at postition 23 in query.");
  is($h_char, "A", "check base at postition 16827 in hit.");
  is($q_gap, "T", "check for the gap at pos 38 in query.");
  is($h_gap, "G", "check base at postition 16847 in hit.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub whatIsThere  {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

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
################################################## subroutine header end ##

sub add_qrna_segments  {
        my $self = shift;
        my $segs = shift;

	foreach my $s (@{$segs}){
		push(@{$self->{qrna_segments}}, $s);
	}
}
################################################# subroutine header end ##

sub add_post_model_p  {
        my $self = shift;

        foreach my $s (@{$self->{qrna_segments}}){
		foreach my $m (keys %{$s->{MODELS}}){
			my $sig_s = $s->{MODELS}->{$m}->{sigmoidal_score};
			my $p     = exp($sig_s)/(1 + exp($sig_s));
			$self->{MODELS}->{$m}->{post_model_p} = $p;

		}		

        }
}

################################################## subroutine header end ##
sub qrna_string  {
        my $self = shift;


	return $self->{qrna_string} if defined($self->{qrna_string});

	my $segs = $self->qrna_segments();
	
	warn " no qrna_segments for hsp in PhatHSP::qrna::qrna_string\n"
	unless $segs;
	warn " you may need to add them see PhatHSP::qrna::add_qrna_segments\n"
	unless $segs;

	my $native_seq = $self->query_string();
	foreach my $s (@{$segs}){
		my $winner  = $s->{winning_model};
                my $winner_score = $s->{MODELS}->{$winner}->{sigmoidal_score};

                my $c;
                if ($winner eq 'OTH'){
                    $c = 'O';
                }
                elsif ($winner eq 'COD'){
                                        $c = 'P';
                } 
                elsif ($winner eq 'RNA'){
                                $c = 'R';
                }
                my $offset  = $s->{MODELS}->{$winner}->{offset};
                my $length  = $s->{MODELS}->{$winner}->{length};
                my $pos_x   = $s->{pos}->{X}->[0];

		my $start = $offset;
		my $end   = $offset + $length;

		make_qrna_str(\$native_seq, 
		              $start, 
		              $end, 
		              $c,
		              );
	}
	
	$self->{qrna_string} = $native_seq;
}

################################################## subroutine header end ##
sub make_qrna_str {
	my $hsp_seq = shift;
	my $start   = shift;
	my $end     = shift;
	my $code    = shift;

	for (my $i = $start; $i < $end; $i++){
		substr($$hsp_seq, $i, 1) = $code;
	} 

}
################################################## subroutine header end ##
sub best_segment  {
        my $self = shift;

	my $segs = $self->qrna_segments();
	

	my @foo;
	foreach my $s (@{$segs}){
		my $win = $s->{winning_model};
		my $score = $s->{MODELS}->{$win}->{sigmoidal_score};

		push(@foo, [$score, $s]);

	}

	my @sorted = sort { $b->[0] <=> $a->[0] } @foo;

	return $sorted[0]->[1];
}

################################################ subroutine header begin ##

=head2 debug_show

 Usage     : How to use this function/method

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub debug_show
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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $name = $hsp->name();

=for example end

=for example_testing
  is($name, "3197985", "Check the name of the sequence that was hit.");

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_begin = $hsp->nB('query');
 my $h_natural_begin = $hsp->nB('hit');

=for example end

=for example_testing
  is($q_natural_begin, 23, "check query's natural begin.");
  is($h_natural_begin, 16827, "check hit's natural begin.");

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
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('qrna',
              'sample_data/qrna.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hsps is filled in by test harness
 my $q_natural_end = $hsp->nE('query');
 my $h_natural_end = $hsp->nE('hit');

=for example end

=for example_testing
  is($q_natural_end, 698, "check query's natural begin.");
  is($h_natural_end, 17492, "check hit's natural begin.");

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
	elsif ($m eq '|'){	# XXXX george added this
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

# XXXX NoOp...
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
	    print STDERR "qrna::PhatHsp::AutoLoader called for: ",
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

