###################################################### main header begin ##

=head1 Bio::Search::Hit::PhatHit::Base

Bio::Search::Hit::PhatHit::Base - The platonic ideal of a CGL module.

=head1 SYNOPSIS

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  # kinky type check to search up inheritance tree (man UNIVERSAL)
  is($hit->isa("Bio::Search::Hit::PhatHit::Base"), 1, "check type.");

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

package Bio::Search::Hit::PhatHit::Base;

use strict;
use warnings;

use Bio::Search::Hit::GenericHit;
use Bit::Vector;

BEGIN {
   use vars qw( $VERSION @ISA );

   $VERSION     = 0.01;
   @ISA = qw(
	     Bio::Search::Hit::GenericHit
	    );
}

################################################ subroutine header begin ##

=head1 new

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness

=for example end

=for example_testing
  # kinky Test::More type check to search up inheritance tree (man UNIVERSAL)
  is($hit->isa("Bio::Search::Hit::PhatHit::Base"), 1, "check type.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub new {
   my $class  = shift;
   my @args   = @_;


   my $self = $class->SUPER::new(@args);

   bless ($self, $class);

   return $self;
}

################################################ subroutine header begin ##

=head1 nB

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_begin = $hit->nB('query');
 my $h_natural_begin = $hit->nB('hit');
=for example end

=for example_testing
  is($q_natural_begin, 23, "Check the query's natural begin.");
  is($h_natural_begin, 16596, "Check the hit's natural begin.");

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

   if ($w eq 'sbjct' || $w eq 'subject'){
       $w = 'hit';
   }

   #set both nB for query and hit on first time through
   if(! defined $self->{nB}){
       $self->{nB} = {};
       $self->nB('query');
       $self->nB('hit');
   }

   if(! defined $self->{nB}{$w}){
       my $low;
       my $high;
       foreach my $hsp ($self->hsps){
	   if(!$hsp->{HIT_START} || !$hsp->{HIT_END}){
	       die;
	   }
	   $low  = $hsp->start($w) if(!$low || $hsp->start($w) < $low);
	   $high = $hsp->end($w) if(!$high || $hsp->end($w) > $high);
       }

       if ($self->strand($w) == 1 || $self->strand($w) == 0) {
	   $self->{nB}{$w} =  $low;
       }
       elsif ($self->strand($w) == -1) {
	   $self->{nB}{$w} =  $high;
       }
       else {
	   die "dead in PhatHit::Base::nB\n";
       }
   }

   return $self->{nB}{$w};
}

################################################ subroutine header begin ##

=head1 nE

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_end = $hit->nE('query');
 my $h_natural_end = $hit->nE('hit');
=for example end

=for example_testing
  is($q_natural_end, 742, "Check the query's natural end.");
  is($h_natural_end, 17492, "Check the hit's natural end.");

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

   if ($w eq 'sbjct' || $w eq 'subject'){
       $w = 'hit';
   }

   #set both nB for query and hit on first time through
   if(! defined $self->{nE}){
       $self->{nE} = {};
       $self->nE('query');
       $self->nE('hit');
   }

   if(! defined $self->{nE}{$w}){
       my $low;
       my $high;
       foreach my $hsp ($self->hsps){
	   $low  = $hsp->start($w) if(!$low || $hsp->start($w) < $low);
	   $high = $hsp->end($w) if(!$high || $hsp->end($w) > $high);
       }

       if ($self->strand($w) == -1) {
	   $self->{nE}{$w} = $low;
       }
       elsif ($self->strand($w) == 1 || $self->strand($w) == 0) {
	   $self->{nE}{$w} = $high;
       }
       else {
	   die "dead in blastn::PhatHit::nE\n";
       }
   }

   return $self->{nE}{$w};
}

################################################ subroutine header begin ##

=head1 nB

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_begin = $hit->nB('query');
 my $h_natural_begin = $hit->nB('hit');
=for example end

=for example_testing
  is($q_natural_begin, 23, "Check the query's natural begin.");
  is($h_natural_begin, 16596, "Check the hit's natural begin.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub start {
   my $self = shift;
   my $w    = shift;
     
   my $b = $self->nB($w);
   my $e = $self->nE($w);

   ($b, $e) = ($e, $b) if ($b > $e);

   return $b;
}

################################################ subroutine header begin ##

=head1 end

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_natural_begin = $hit->end('query');
 my $h_natural_begin = $hit->end('hit');
=for example end

=for example_testing
  is($q_natural_begin, 23, "Check the query's natural begin.");
  is($h_natural_begin, 16596, "Check the hit's natural begin.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub end {
   my $self = shift;
   my $w    = shift;
     
   my $b = $self->nB($w);
   my $e = $self->nE($w);

   ($b, $e) = ($e, $b) if ($b > $e);

   return $e;
}


################################################ subroutine header begin ##

=head1 id

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $id1 = $hit->id('Testing321');
 my $id2 = $hit->id();
 my $id3 = $hit->id(undef);

=for example end

=for example_testing
  # kinky type check to search up inheritance tree (man UNIVERSAL)
  is($id1, 'Testing321', "Check id getter/setter (setting).");
  is($id2, 'Testing321', "Check id getter/setter (getting).");
  is($id3, undef, "Check id getter/setter (unsetting).");

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

   if (@_) {
      $self->{id} = shift;
   }

   return $self->{id};
}

################################################ subroutine header begin ##

=head1 strand

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_strand = $hit->strand('query');
 my $h_strand = $hit->strand('hit');
 my @strand = $hit->strand();

=for example end

=for example_testing
  is($q_strand, 1, "Check query strand.");
  is($h_strand, 1, "Check hit strand.");
  is(scalar(@strand), 2, "Check the size of the strand array");
  is($strand[0], 1, "Check query strand from the array.");
  is($strand[1], 1, "Check hit strand from the array.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub strand {
   my ($self, $seqType) = @_;

   $seqType ||= (wantarray ? 'list' : 'query');
   $seqType = 'sbjct' if $seqType eq 'hit';

   # no hsps because of memory optimization
   if(! defined($self->{strand}{query})){
       my (%qstr, %hstr);
       foreach my $hsp ( $self->hsps ) {
	   my $q = $hsp->strand('query');
	   my $h = $hsp->strand('hit');
	   $qstr{$q}++;
	   $hstr{$h}++;
       }
       
       #changed 2/12/2009 to act like current version of Bioperl
       my $topq;
       my $qstr;
       while(my $key = each %qstr){
	   if(! $topq || $qstr{$key} > $topq){
	       $qstr = $key;
	       $topq = $qstr{$key};
	   }
       }
       
       my $toph;
       my $hstr;
       while(my $key = each %hstr){
	   if(! $toph || $hstr{$key} > $toph){
	       $hstr = $key;
	       $toph = $hstr{$key};
	   }
       }
       
       $self->{strand}{query} = $qstr;
       $self->{strand}{hit} = $hstr;
   }
   
   if($seqType =~ /list|array/i){
       return ($self->{strand}{query}, $self->{strand}{hit});
   }
   elsif($seqType eq 'query'){
       return $self->{strand}{query};
   }
   else{
       return $self->{strand}{hit};
   }
}

################################################ subroutine header begin ##

=head1 hsps

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my @hsps = $hit->hsps();
 my $hsp_count = $hit->hsps();

=for example end

=for example_testing
  isa_ok($hsps[0], "Bio::Search::HSP::PhatHSP::blastn", "check type.");
  is($hsp_count, 2, "Check the number of hsps in the hit.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

# XXXX
# think about changing this to return a reference, and to be in
# cjm-normal-form.
sub hsps {
   my $self = shift;
   my $hsps = shift;

   if (not ref $self->{'_hsps'}) {
      # XXXX throw!
      $self->throw("Can't get HSPs: data not collected.");
   }
   elsif (defined($hsps)) {
      $self->{'_hsps'} = [@$hsps]; #creates copy of the array reference
      delete $self->{nB};
      delete $self->{nE};
   }
   return wantarray ? @{$self->{'_hsps'}} : scalar(@{$self->{'_hsps'}});

}

################################################ subroutine header begin ##

# XXXX percent aligned query

=head1 pAq

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $percent_aligned = $hit->pAq(775); # query length from blast report.

=for example end

=for example_testing
  is($percent_aligned, 153.8, "Check percent aligned in query.");

 Purpose   : Calculate the "percent aligned query" statistic....
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub pAq {
   my $self    = shift;
   my $ql      = shift;

   #XXXX danger will robinson
   $ql = $self->queryLength() unless defined($ql);

   die "\$self->queryLength undefined; define or provide as arg!\n"
       unless defined($ql);

   my ($laq, $lah) = $self->getLengths();
   my $perAl = ($laq/$ql);

   return sprintf '%.4f', $perAl;
}

################################################ subroutine header begin ##

# XXXX percent aligned hit
=head1 pAh

Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $percent_aligned = $hit->pAh();

=for example end

=for example_testing
  is($percent_aligned, 1.613, "Check percent aligned in hit.");

 Purpose   : Computer the "percent aligned hit" statistic....
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub pAh {
   my $self = shift;
   my $hl   = shift;

   $hl = $self->length() unless defined($hl);

   die "\$self->length undefined; define or provide as arg!\n"
       unless defined($hl);

   my ($laq, $lah) = $self->getLengths();
   my $perAl = ($lah/$hl);

   return sprintf '%.4f', $perAl;
}

################################################ subroutine header begin ##

=head1 getLengths

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my ($q_length, $h_length) = $hit->getLengths();

=for example end

=for example_testing
  is($q_length, 1192, "Check the query length.");
  is($h_length, 1179, "Check the hit length.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub getLengths {
   my $self = shift;

   if (defined ($self->{LAlnQ}) && defined ($self->{LAlnH})){
         return $self->{LAlnQ}, $self->{LAlnH};
   }

   my $q_e = $self->nE('query');
   my $q_b = $self->nB('query');
   my $h_e = $self->nE('hit');
   my $h_b = $self->nB('hit');

   ($q_b, $q_e) = ($q_e, $q_b) if ($q_b > $q_e);
   ($h_b, $h_e) = ($h_e, $h_b) if ($h_b > $h_e);

   my $qLen = abs($q_e - $q_b) + 1;
   my $hLen = abs($h_e - $h_b) + 1;

   my $qOffset = $q_b - 1;
   my $hOffset = $h_b - 1;

   my $q_vec = Bit::Vector->new($qLen + 1); #pretend 0 coor doesn't exist
   my $h_vec = Bit::Vector->new($hLen + 1); #pretend 0 coor doesn't exist

   foreach my $s ($self->hsps) {
      my $q_hspEnd = $s->end('query') - $qOffset;
      my $h_hspEnd = $s->end('hit') - $hOffset;

      my $q_hspStart = $s->start('query') - $qOffset;
      my $h_hspStart = $s->start('hit') - $hOffset;

      $q_vec->Interval_Fill($q_hspStart, $q_hspEnd);
      $h_vec->Interval_Fill($h_hspStart, $h_hspEnd);
   }
   
   my $laq = $q_vec->Norm();
   my $lah = $h_vec->Norm();

   $self->{LAlnQ} = $laq;
   $self->{LAlnH} = $lah;

   return ($laq, $lah);
}

################################################ subroutine header begin ##

=head1 show

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

   print "---------------------------------------------------------\n";
   print "-- ", ref($self), " -------------\n";
   print "---------------------------------------------------------\n";
   print $self->name()."\n";
   print "is_split:".$self->is_split()."\n";
   print "sQ:".$self->strand('query')." sH:". $self->strand('hit')."\n";
   print "nBq:".$self->nB('query')." nEq:".$self->nE('query')."\n";
   print "nBh:".$self->nB('hit')." nEh:".$self->nE('hit')."\n";
   print "pAq:". $self->pAq."\n";
   print "pAH:". $self->pAh."\n";
   print "E/P:". $self->significance()."\n";
   print "queryLength:".$self->queryLength()."\n";
   print "sbjctLength:".$self->length()."\n";
   my $i = 0;
   foreach my $hsp ($self->hsps) {
      print "-------------- HSP:$i -----------------\n";
      $hsp->show();
      $i++;
   }
   print "---------------------------------------------------------\n";
}

################################################ subroutine header begin ##

=head1 frac_identical

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_hsps = $hit->frac_identical();
 my $h_hsps = $hit->frac_identical();

=for example end

=for example_testing
 is($q_hsps->[0]->frac_identical(), .90, "Check fraction identical.");

 Purpose   : Give fraction identical without HSP tiling.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub frac_identical {
   my $self = shift;

   return $self->{_f_id} if (defined ($self->{_f_id}));
   
   my $t_f;
   my $t_l;
   
   foreach my $s ($self->hsps){
       $t_l += $s->length;
       $t_f += $s->num_identical;
   }

   $self->{_f_id} = sprintf("%.3f", $t_f/$t_l);
  
   return $self->{_f_id};
}

################################################ subroutine header begin ##
=head1 revSortFeatures

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];          # $hits is filled in by test harness
 my $q_hsps = $hit->revSortFeatures('query');
 my $h_hsps = $hit->revSortFeatures('hit');

=for example end

=for example_testing
 is($q_hsps->[0]->nE('query'), 742, "Check the end of the first hsp (query).");
 is($h_hsps->[0]->nE('hit'), 17492, "Check the end of the first hsp (hit).");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.
                                                                                                                                                                 
=cut

################################################## subroutine header end ##

sub revSortFeatures
{
    my $self = shift;
    my $who  = shift;

    my @rev = reverse sort {$a->end($who) <=> $b->end($who)} $self->hsps;

    return \@rev;
}


################################################ subroutine header begin ##

=head1 sortFeatures

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $q_hsps = $hit->sortFeatures('query');
 my $h_hsps = $hit->sortFeatures('hit');

=for example end

=for example_testing
 is($q_hsps->[0]->nB('query'), 23, "Check the end of the first hsp (query).");
 is($h_hsps->[0]->nB('hit'), 16596, "Check the end of the first hsp (hit).");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub sortFeatures
{
   my $self = shift;
   my $who  = shift;

   my @subfeatures = sort {$a->start($who) <=> $b->start($who)} $self->hsps;

   return \@subfeatures;
}

################################################ subroutine header begin ##

=head1 sortedHSPs

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $hsps = $hit->sortedHSPs();

=for example end

=for example_testing

 Purpose   : Sorts HSPs to be in order od transcription/translation.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub sortedHSPs
{
   my $self = shift;
   my $who  = 'query';

   my @subfeatures = ($self->strand == 1) ? sort {$a->start($who) <=> $b->start($who)} $self->hsps :
       sort {$b->end($who) <=> $a->end($who)} $self->hsps;

   return \@subfeatures;
}

################################################ subroutine header begin ##

=head1 queryLength

 Usage     : How to use this function/method

=for example
 use Bio::Search::Hit::PhatHit::Base;
 my $hits = Bio::Search::Hit::PhatHit::Base::_getTestHits('blastn',
              'sample_data/blastn.sample.report');

=for example begin

 my $hit = $hits->[0];		# $hits is filled in by test harness
 my $query_length = $hit->queryLength(775);

=for example end

=for example_testing
  is($query_length, 775, "Check the query length.");

 Purpose   : What the subroutine does.
 Returns   : The types and values it returns.
 Argument  : Required and optional input.
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.
 See Also  : Other things that might be useful.

=cut

################################################## subroutine header end ##

sub queryLength
{
   my $self = shift;

   if (@_) {
      $self->{queryLength} = shift;
   }

   return $self->{queryLength};
}

################################################ subroutine header begin ##

=head1 _getTestHits

 Usage     : *private*

 Purpose   : Load a blast report so that the inline tests can have
           : something to play with.
 Returns   : A reference to an array of PhatHits.
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

sub _getTestHits {
   my($type, $report) = @_;

   my $sio;			# the search IO object

   my $hitFactory;	      # hit and hsp factory objects for search
   my $hspFactory;
   my $result;			# a blast result
   my $hit;			# a blast hit
   my $hit_aref;		# reference to an array of hits


   $sio = new Bio::SearchIO(-format => 'blast', -file   => $report);

   $hitFactory = new Bio::Search::Hit::HitFactory(-type =>
						  'Bio::Search::Hit::PhatHit::' . $type);
   $hspFactory = new Bio::Search::HSP::HSPFactory(-type =>
						  'Bio::Search::HSP::PhatHSP::' . $type);

   $sio->_eventHandler->register_factory('hit', $hitFactory);
   $sio->_eventHandler->register_factory('hsp', $hspFactory);

   while ($result = $sio->next_result) {
      while ($hit = $result->next_hit) {
	 push @{$hit_aref}, $hit;
      }
   }

   return($hit_aref);
}

################################################ subroutine header begin ##

=head1 AUTOLOAD

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

   if ($ENV{CGL_CHATTER}) {
      print STDERR "PhatHit::AutoLoader called for: ",
      "\$self->$call","()\n";
      print STDERR "call to AutoLoader issued from: ", $caller, "\n";
   }

   if (@_) {
      $self->{$call} = shift;
   }

   return $self->{$call};
}

1; #this line is important and will help the module return a true value
__END__

