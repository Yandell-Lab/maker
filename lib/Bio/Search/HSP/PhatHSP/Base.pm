###################################################### main header begin ##

=head1 Bio::Search::HSP::PhatHSP::Base

Bio::Search::HSP::PhatHSP::Base - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  # kinky type check to search up inheritance tree (man UNIVERSAL)
  is($hsp->isa("Bio::Search::HSP::PhatHSP::Base"), 1, "Check type.");

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

package Bio::Search::HSP::PhatHSP::Base;

use strict;
use warnings;

use Bio::Search::HSP::GenericHSP;

BEGIN {
  use vars qw( $VERSION @ISA );

  $VERSION     = 0.01;
  @ISA = qw(
	    Bio::Search::HSP::GenericHSP
	   );
}

################################################ subroutine header begin ##

=head2 new

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  # kinky type check to search up inheritance tree (man UNIVERSAL)
  is($hsp->isa("Bio::Search::HSP::PhatHSP::Base"), 1, "Check type.");

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

	#die;
	my $self = $class->SUPER::new(@args);

	bless ($self, $class);

	return $self;
}

################################################ subroutine header begin ##

=head2 whatIsInTheMiddle

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $q_mid1 = $hsp->whatIsInTheMiddle('query', 23);
 my $q_mid2 = $hsp->whatIsInTheMiddle('query', 25);
 my $h_mid1 = $hsp->whatIsInTheMiddle('hit', 16827);
 my $h_mid2 = $hsp->whatIsInTheMiddle('hit', 16829);

=for example end

=for example_testing
  is($q_mid1, '|', "Check the middle character using query coordinates.");
  is($h_mid1, '|', "Check the middle character using hit coordinates.");
  is($q_mid2, ' ', "Check the middle character using query coordinates.");
  is($h_mid2, ' ', "Check the middle character using hit coordinates.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
           : XXXX correct behaviour, can't specify a location in 
           : query that's a gap in the alignment.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub whatIsInTheMiddle
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

# XXXX rewrite to factor out common code.
        if ($w eq 'query'){
                my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
                return undef unless defined($q_i);
                if ($q eq '-'){
                        return undef;
                }
                else {
                        return $m;
                }
        }
        else {
                my ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);
                return undef unless defined($h_i);
                if ($h eq '-'){
                        return undef;
                }
                else {
                        return $m;
                }
        }

}


################################################ subroutine header begin ##

=head2 nB

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

sub nB {
    my $self = shift;
    my $w    = shift || 'query';
    
    return ($self->strand($w) eq '-1') ? $self->end($w) : $self->start($w);
}

################################################ subroutine header begin ##

=head2 nE

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

sub nE {
    my $self = shift;
    my $w    = shift || 'query';
    
    return ($self->strand($w) eq '-1') ? $self->start($w) : $self->end($w);
}

################################################ subroutine header begin ##

=head2 equivalent_pos_in_alignment_partner

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $h_via_q1 = $hsp->equivalent_pos_in_alignment_partner('query', 23);
 my $q_via_h1 = $hsp->equivalent_pos_in_alignment_partner('hit', 16827);
 # does it handle a gappy position.
 my $h_via_q2 = $hsp->equivalent_pos_in_alignment_partner('query', 43);
 my $q_via_h2 = $hsp->equivalent_pos_in_alignment_partner('hit', 16842);

=for example end

=for example_testing
  is($h_via_q1, 16827, "Find position in hit using query position 23.");
  is($q_via_h1, 23, "Find position in query using hit position 16827.");
  is($h_via_q2, undef, "Find position in hit using query position 43.");
  is($q_via_h2, undef, "Find position in query using hit position 16842.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub equivalent_pos_in_alignment_partner {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

        my ($q_i, $h_i, $q, $h, $m);
        my $equiv_pos = undef;
        my $equiv_char = undef;

	# look down the alignment for the relevant info
        ($q_i, $h_i, $q, $h, $m) = $self->_scan($w, $pos);

	# choose the appropriate postition/character depending on context.
	$equiv_pos = ($w eq 'query') ? $h_i : $q_i;
	$equiv_char = ($w eq 'query') ? $h : $q;

	return undef if (! defined($equiv_pos)); # no equiv. pos. in partner
	return undef if ($equiv_char eq '-'); # equiv. pos. is a gap.
	return $equiv_pos;	# phew, actually got an answer.
}

################################################ subroutine header begin ##

=head2 id

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $id1 = $hsp->id("Testing321");
 my $id2 = $hsp->id();
 my $id3 = $hsp->id(undef);

=for example end

=for example_testing
  # kinky type check to search up inheritance tree (man UNIVERSAL)
  is($id1, 'Testing321', "Check the id getter/setter (setting).");
  is($id2, 'Testing321', "Check the id getter/setter (getting).");
  is($id3, undef, "Check the id getter/setter (clearing).");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub id
 {
        my $self = shift;

        if (@_){
                $self->{id} = shift;
        }

	return $self->{id};
}

################################################ subroutine header begin ##

=head2 isContained

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $q_before = $hsp->isContained('query', 20);
 my $q_in = $hsp->isContained('query', 25);
 my $q_after = $hsp->isContained('query', 1000);

 my $h_before = $hsp->isContained('hit', 1000);
 my $h_in = $hsp->isContained('hit', 16900);
 my $h_after = $hsp->isContained('hit', 17500);

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
    # query strand == 1 or query strand == 0
    if ($self->strand('query') == 1 || $self->strand('query') == 0){
      # XXXX covers blastn and blastp bioperl strand ambiguity.
      return 1 if ($pos >= $self->nB('query') && $pos <= $self->nE('query'));
    }
    # query strand == -1
    else {
      return 1 if ($pos <= $self->nB('query') && $pos >= $self->nE('query'));
    }
  }
  else {
    # hit strand == 1
    if ($self->strand('hit') == 1) {
      return 1 if ($pos >= $self->nB('hit') && $pos <= $self->nE('hit'));
    }
    # hit strand == 0 or query strand == -1
    else {
      return 1 if ($pos <= $self->nB('hit') && $pos >= $self->nE('hit'));
    }
  }
  return 0;
}

################################################ subroutine header begin ##

=head2 _scan

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

sub _scan
 {
        my $self = shift;
        my $w    = shift;
        my $pos  = shift;

        my @q = split('',$self->query_string());
        my @h = split('',$self->hit_string());
        my @m = split('',$self->homology_string());

        die "q, h, and match strings not equal!\n"
        unless ($#q == $#h && $#h == $#m);

        my ($m_i, $q_i, $h_i);
        while  (my $m = shift(@m)){
                my $q = shift(@q);
                my $h = shift(@h);

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

################################################ subroutine header begin ##

=head2 show

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

sub show
 {
	my $self = shift;

        print "--------------------------------------\n";
        print "-- ", ref($self), " --\n";
	print "   ".$self->hit->seqname()."\n";
        print "--------------------------------------\n";
	print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
	print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
        print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
	print "fIDt:".$self->frac_identical('total')."\n";
	print "bits:".$self->bits."\n";
	print "significance:".$self->significance."\n";
}

################################################ subroutine header begin ##

=head2 name

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $name = $hsp->name();

=for example end

=for example_testing
  is($name, "3197985", "Check the name of the hit.");

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

=head2 strand

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $strand = $hsp->strand();

=for example end

=for example_testing
  is($name, "3197985", "Check the name of the hit.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub strand
 {
	my $self = shift;
	my $what = shift;

	if ($what =~ /subject|sbjct/i){
	    $what = 'hit';
	}

	if($what eq 'hit'){
	    return $self->{_strand_hack}{hit} if(exists $self->{_strand_hack});
	    return $self->{STRAND_HIT} if(exists $self->{STRAND_HIT});
	    return $self->SUPER::strand('hit');
	}
	elsif($what eq 'query'){
	    return $self->{_strand_hack}{query} if(exists $self->{_strand_hack});
	    return $self->{STRAND_QUERY} if(exists $self->{STRAND_QUERY});
	    return $self->SUPER::strand('query');
	}
	else{
	    return $self->SUPER::strand($what);
	}
}

################################################ subroutine header begin ##

=head2 _get_args

 Usage     : *private*

 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  : Create an 'args array' for 'newing' a generic Bioperl
           : HSP object.
 See Also  : Bio::yada-yada::HSP XXXX

=cut

################################################## subroutine header end ##

sub _get_args
 {
	my $self = shift;


	my @args;
	push(@args, '-query_start');
	push(@args, $self->nB('query'));

	push(@args, '-score');
	push(@args, $self->score());

	push(@args, '-homology_seq');
	push(@args, $self->homology_string());

        push(@args, '-hit_start');
        push(@args, $self->nB('hit'));

        push(@args, '-hit_seq');
        push(@args, $self->seq_str('hit'));

        push(@args, '-hsp_length');
        push(@args, $self->length);

        push(@args, '-identical');
        push(@args, $self->num_identical);

        push(@args, '-hit_length');
        push(@args, $self->length('hit'));

        push(@args, '-query_name');
        push(@args, $self->query->seqname);

        push(@args, '-algorithm');
        push(@args, $self->algorithm);

        push(@args, '-bits');
        push(@args, $self->bits);

        push(@args, '-evalue');
        push(@args, $self->evalue);

        push(@args, '-pvalue');
        push(@args, $self->pvalue);

        push(@args, '-query_length');
        push(@args, $self->length('query'));

        push(@args, '-query_end');
        push(@args, $self->nE('query'));

        push(@args, '-conserved');
        push(@args, $self->num_conserved);

        push(@args, '-hit_name');
        push(@args, $self->hit->seqname());

        push(@args, '-hit_end');
        push(@args, $self->nE('hit'));

	my $q_frame_2 =  $self->query->frame() + 1;
	my $s_frame_2 =  $self->hit->frame()   + 1;

	$q_frame_2 = -1*$q_frame_2 if $self->strand('query') == -1;
	$s_frame_2 = -1*$s_frame_2 if $self->strand('hit') == -1;	

        push(@args, '-query_frame');
        push(@args, $q_frame_2);

        push(@args, '-hit_frame');
        push(@args, $s_frame_2);

        push(@args, '-query_seq');
        push(@args, $self->seq_str('query'));

	return \@args;
}

################################################ subroutine header begin ##

=head2 hasRun

 Usage     : How to use this function/method

=for example
 use Bio::Search::HSP::PhatHSP::Base;
 my $hsps = Bio::Search::HSP::PhatHSP::Base::_getTestHSPs('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hsp = $hsps->[0];		# $hits is filled in by test harness
 my $yes = $hsp->hasRun('|', 6);
 my $no = $hsp->hasRun('|', 10);

=for example end

=for example_testing
  is($yes, 1, "Does the hsp has a run of 6 matches (it should)?");
  is($no, 0, "Does the hsp has a run of 10 matches (it shouldn't)?");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub hasRun
 {
	my $self    = shift;
	my $m       = shift;
	my $l       = shift;
	my $offset  = shift;
	my $lRegion = shift;

	my $homology_string;
	if (defined($offset)){
		return 0 unless defined(substr($self->homology_string, $offset, $lRegion));

		$homology_string =
		substr($self->homology_string, $offset, $lRegion);

	}
	else {
		$homology_string = $self->homology_string();
	}

	my $test;
	if ($m eq 'A-Z+'){
		$test = '[A-Z\+]'."{".$l."}"
	}
	else {
		$test =  $m =~ /[\|\+]/ ?  "[".quotemeta($m)."]"."{".$l."}" : "[".$m."]"."{".$l."}";
	}

	if ($homology_string =~ /$test/){
		return 1;
	}
	else {
		return 0;
	}
}

################################################ subroutine header begin ##

=head2 _getTestHSPs

 Usage     : *private*

 Purpose   : Load a blast report so that the inline tests can have
           : something to play with.
 Returns   : A reference to an array of PhatHSPs.
 Argument  : A type ('blastn', 'blastx', etc...) and a report filename.
 Throws    :
 Comments  :
           :
 See Also  : Bio::SearchIO

=cut

################################################## subroutine header end ##


use Bio::SearchIO;
use Bio::Search::Hit::HitFactory;
use Bio::Search::HSP::HSPFactory;

sub _getTestHSPs {
  my($type, $report) = @_;

  my $sio;			# the search IO object

  my $hitFactory;		# hit and hsp factory objects for search
  my $hspFactory;
  my $result;			# a blast result
  my $hit;			# a blast hit
  my $hsp;			# a blast hsp
  my $hsp_aref;			# reference to an array of hsp's


  $sio = new Bio::SearchIO(-format => 'blast', -file   => $report);

  $hitFactory = new Bio::Search::Hit::HitFactory(-type =>
                'Bio::Search::Hit::PhatHit::' . $type);
  $hspFactory = new Bio::Search::HSP::HSPFactory(-type =>
                'Bio::Search::HSP::PhatHSP::' . $type);

  $sio->_eventHandler->register_factory('hit', $hitFactory);
  $sio->_eventHandler->register_factory('hsp', $hspFactory);

  while($result = $sio->next_result) {
    while($hit = $result->next_hit) {
      while(my $hsp = $hit->next_hsp) {
	push @{$hsp_aref}, $hsp;
      }
    }
  }

  return($hsp_aref);
}

################################################ subroutine header begin ##

=head2 _getTestHSPs

 Usage     : 

 Purpose   : produces a cigar string for use in GFF3 etc.
 Returns   : A string
 Argument  : None
 Throws    :
 Comments  :
           :
 See Also  : 

=cut

################################################## subroutine header end ##


sub cigar_string {
    my $self = shift;

    return($self->{_CIGAR}) if($self->{_CIGAR});

    my $string = $self->SUPER::cigar_string();
    my $cigar = '';
    my $type = '';
    my $value = 0;
    while($string =~ /(\d+)([A-Z])/g){
	if($type eq $2){
	    $value += $1;
	}
	else{
	    $cigar .= "$type$value" if($value && $type);
	    $type = $2;
	    $value = $1;
	}
    }

    $cigar .=  "$type$value" if($value && $type);

    return ($self->{_CIGAR} = $cigar);
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
	    print STDERR "PhatHsp::AutoLoader called for: ",
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

