#------------------------------------------------------------------------
#----                        Widget::blastn                          ---- 
#------------------------------------------------------------------------
package Widget::blastn;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::Search::Hit::PhatHit::blastn;
use Exporter;
use PostData;
use FileHandle;
use Widget;
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

	if (defined($command)){
		$self->print_command($command);
		system("$command");
	}
	else {
		die "you must give Widget::blastn a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
        my $report = shift;
        my $params = shift;

        my $hitType = 'Bio::Search::Hit::PhatHit::blastn';
        my $hspType = 'Bio::Search::HSP::PhatHSP::blastn';

        my $sio = new Bio::SearchIO(-format => 'blast',
                                    -file   => $report,
                                   );

        my $hspFactory = new Bio::Search::HSP::HSPFactory(-type =>$hspType);
        my $hitFactory = new Bio::Search::Hit::HitFactory(-type =>$hitType);


        $sio->_eventHandler->register_factory('hsp', $hspFactory);
        $sio->_eventHandler->register_factory('hit', $hitFactory);

        return keepers($sio, $params);

}
#-------------------------------------------------------------------------------
sub keepers {
        my $sio    = shift;
	my $params = shift; 

        my @keepers;
	my $start = 0;
	
        while (my $result = $sio->next_result()){
	    $start += $result->hits();
	    while(my $hit = $result->next_hit) {
                my $significance = $hit->significance();
                $significance = "1".$significance if  $significance =~ /^e/;
		$significance = 0                 if  $significance =~ /0\./;
                $hit->queryLength($result->query_length);
                $hit->queryName($result->query_name);
                next unless $significance < $params->{significance};
                my @hsps;
                while(my $hsp = $hit->next_hsp) {
		    $hsp->query_name($result->query_name);
		    
		    push(@hsps, $hsp) if $hsp->bits > $params->{hsp_bit_min};
		    
		}
                $hit->hsps(\@hsps);
		if ( $hit->strand('query') eq '-1/1' || $hit->strand('query') eq '1/-1'){
		    my ($p_hit, $m_hit) = clean_mixed_strands($hit);
		    push(@keepers, $p_hit) unless $p_hit eq 'no_hsps';
		    push(@keepers, $m_hit) unless $m_hit eq 'no_hsps';
		}
                else {
		    push(@keepers, $hit) if $hit->hsps();
		}
		
	    }
	}
        my $end     = @keepers;
        my $deleted = $start - $end;
        print STDERR "deleted:$deleted hits\n";
	
        return \@keepers;
    }
#-----------------------------------------------------------------------------
sub clean_mixed_strands {
        my $hit = shift;

        my @plus;
        my @minus;
        while(my $hsp = $hit->next_hsp) {
                push(@plus, $hsp) if $hsp->strand('query') ==  1;
                push(@minus,$hsp) if $hsp->strand('query') == -1;
        }

      my $p_hit = new Bio::Search::Hit::PhatHit::blastn(
	                              '-name'        => $hit->name()." plus HSPs only",
                                      '-description' => $hit->name()." cleaned_mixed_strands",
                                      '-algorithm'   => 'blastx_',
                                      '-length'      => $hit->length(),
                                     );
      my $m_hit = new Bio::Search::Hit::PhatHit::blastn(
	                              '-name'        => $hit->name()." minus HSPs only",
                                      '-description' => $hit->name()." cleaned_mixed_strands",
                                      '-algorithm'   => 'blastn',
                                      '-length'      => $hit->length(),
                                     );

	if (defined $plus[0]){
		$p_hit->hsps(\@plus);
		$p_hit->queryName($hit->queryName);
	}
	else {
		$p_hit = 'no_hsps';
	}

        if (defined $minus[0]){
                $m_hit->hsps(\@minus);
		$m_hit->queryName($hit->queryName);
        }
        else {
        }
                $m_hit = 'no_hsps';

        return ($p_hit, $m_hit);

}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::blastn::AutoLoader called for: ",
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


