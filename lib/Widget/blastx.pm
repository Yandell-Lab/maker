#------------------------------------------------------------------------
#----                        Widget::blastx                          ---- 
#------------------------------------------------------------------------
package Widget::blastx;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::Search::Hit::PhatHit::blastx;
use Exporter;
use PostData;
use FileHandle;
use Widget;
use PhatHit_utils;
use IPC::Open3;

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
		my $pid = open3(\*CHLD_IN, \*CHLD_OUT, \*CHLD_ERR, $command);
		local $/ = \1;
		while (my $line = <CHLD_ERR>){
		   print STDERR $line unless($main::quiet);
		}
		waitpid $pid, 0;
	}
	else {
		die "you must give Widget::blastx a command to run!\n";
	}
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
        my $report = shift;
        my $params = shift;

        my $hitType = 'Bio::Search::Hit::PhatHit::blastx';
        my $hspType = 'Bio::Search::HSP::PhatHSP::blastx';

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
      my $hits = [$result->hits()];
      $start += @{$hits};
      
      foreach my $hit (@{$hits}){
	 $hit->queryLength($result->query_length);
	 $hit->queryName($result->query_name);
	 
	 my @hsps;
	 while(my $hsp = $hit->next_hsp) {
		 $hsp->query_name($result->query_name);
		 
		 push(@hsps, $hsp) if $hsp->bits > $params->{hsp_bit_min};
	      }
	 $hit->hsps(\@hsps);
      }

      $hits = PhatHit_utils::split_hit_by_strand($hits);
      if (exists $params->{split_hit}){
	 $hits = PhatHit_utils::split_hit_on_intron($hits, $params->{split_hit});
      }

      foreach my $hit (@{$hits}) {
	 my $significance = $hit->significance();
	 $significance = "1".$significance if  $significance =~ /^e/;
	 $significance = 0                 if  $significance =~ /0\./;
	 next unless $significance < $params->{significance};

	 #next unless $hit->pAh > $params->{percov};
	 #next unless $hit->hsp('best')->frac_identical() > $params->{percid};
	 #next unless PhatHit_utils::is_contigous($hit);
	 
	 push(@keepers, $hit) if $hit->hsps();
      }
   }
   my $end     = @keepers;
   my $deleted = $start - $end;
   print STDERR "deleted:$deleted hits\n" unless $main::quiet;
   
   return \@keepers;
}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "Widget::blastx::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if (defined($arg)){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------

1;


