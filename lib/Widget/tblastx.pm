#------------------------------------------------------------------------
#----                        Widget::tblastx                          ---- 
#------------------------------------------------------------------------
package Widget::tblastx;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Bio::Search::Hit::PhatHit::tblastx;
use Exporter;
use PostData;
use FileHandle;
use Widget;
use PhatHit_utils;
use IPC::Open3;
use Symbol;

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
      my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym);
      my $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);

      my $all_err;
      {
	  local $/ = \1;
	  while (my $line = <$CHLD_ERR>){
	      $all_err .= $line;
	      print STDERR $line unless($main::quiet);
	  }
      }
      waitpid $pid, 0;
      if ($? != 0 && ! ignoreable($all_err, $command)){
	 #try again on filure
	 sleep 15;
	 $pid = open3($CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command);
	 $all_err = '';
	 {
	     local $/ = \1;
	     while (my $line = <$CHLD_ERR>){
		 $all_err .= $line;
		 print STDERR $line unless($main::quiet);
	     }
	 }
	 waitpid $pid, 0;
	 
	 if ($? != 0 && ! ignoreable($all_err, $command)){
	    die "ERROR: TBLASTX failed\n";
	 }
      }
   }
   else {
      die "you must give Widget::tblastx a command to run!\n";
   }
}

sub ignorable{
    my $err = shift;
    my $command = shift;

    my $ignorable = 0;
    if($err =~ /There are no valid contexts/){
        $ignorable = 1;
    }
    elsif($err =~ /invalid query sequence/){
        $ignorable = 1;
    }
    elsif($err =~ /calculate ungapped Karlin-Altschul/){
        $ignorable = 1;
    }
    elsif($err =~ /sequence has zero length/){
        $ignorable = 1;
    }
    elsif($err =~ /Sequence contains no data/){
        $ignorable = 1;
    }
    else{
        $ignorable = 0;
    }

    #concatenate error into file (fixes BLAST+ 2.29 issue)
    if($ignorable && $command && $command =~ /\s\-+(out|output_file)\s+([^\s\t\n]+)/){
        my $file = $2;
        return $ignorable if(! -f $file);

        open(OUT, ">>$file");
        print OUT $err;
        close(OUT);
    }

    return $ignorable;
}

#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub parse {
        my $report = shift;
        my $params = shift;

        my $hitType = 'Bio::Search::Hit::PhatHit::tblastx';
        my $hspType = 'Bio::Search::HSP::PhatHSP::tblastx';

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
   my $clusters = [];
   my $start = 0;
   my $ok; #context override (check if report is finished)

   my $depth = $params->{depth};
   my $split_hit = $params->{split_hit} || 0;
   my $hsp_bit_min = $params->{hsp_bit_min} || 0;

   #iterate through results
   my $i = 0;
   while (my $result = $sio->next_result()){
      if(!$ok && !$result->get_statistic('posted_date')){
	  my $ok;
	  open(my $IN, '<', $sio->file);
	  while(my $line = <$IN>){
	      if(ignorable($line)){
		  $ok = 1;
		  last;
	      }
	  }
	  close($IN);
	  die "ERROR: TBLASTX does not appear to be finished in Widget::tblastx::keepers\n" if(!$ok);
      }

      $start += $result->num_hits;

      #iterate hits
      my $q_length = $result->query_length;
      my $q_name = $result->query_name;
      my $db_name = $result->database_name;
      my $cutoff = ($params->{is_last}) ? $q_length+1 : $q_length - $split_hit;
      my $scutoff = ($params->{is_first}) ? 0 : 1 + $split_hit;

      while (my $hit = $result->next_hit){
          $hit->queryLength($q_length);
          $hit->queryName($q_name);
          $hit->database_name($db_name);

          #iterate HSP
	  my @hsps;
          while(my $hsp = $hit->next_hsp) {
	      #clear a little memory (will I need this later?)
	      $hsp->cigar_string; #holds alignment
	      $hsp->{QUERY_SEQ} = '';
	      $hsp->{HIT_SEQ} = '';
	      $hsp->{HOMOLOGY_SEQ} = '';	      

              $hsp->query_name($q_name);
              push(@hsps, $hsp) if($hsp->bits > $hsp_bit_min);
          }
          $hit->hsps(\@hsps);

          #split hits
	  my $hits = PhatHit_utils::split_hit_by_strand([$hit]);
          $hits = PhatHit_utils::split_hit_on_intron($hits, $split_hit) if($split_hit);

	  #fix strand
	  foreach my $h (@{$hits}) {
              next unless($h->hsps());	      
              #fix strand for messed up ncbi blast
	      $hit = PhatHit_utils::copy($hit, 'both') if ($hit->strand('hit') < 0);
	  }

	  #fix duplication across strands (common in tblastx)
	  my $fixdup;
	  @$hits = sort {$b->score <=> $a->score} @$hits;
	  foreach my $h (@$hits){
              next unless($h->hsps());
	      my $bad;
	      foreach my $k (@$fixdup){
		  if(compare::compare($h->start('query'), $h->end('query'),
				      $k->start('query'), $k->end('query'))
		      ){
		      $bad = 1;
		      last;
		  }
	      }
	      push(@$fixdup, $h) unless($bad);
	  }
	  $hits = $fixdup;

          #filter split hits	      
	  foreach my $h (@{$hits}) {
              next unless($h->hsps());
              my $s = $h->start('query');
              my $e = $h->end('query');
              if($split_hit && ($s <=  $scutoff || $e >= $cutoff)){
		  $h->start(); #force create of stored start and end
		  $h->strand(); #force creation
		  $h->getLengths(); #force creation of stored lengths
		  if(!_keep_filt($h, $params)){ #memory optimization
		      $h->{_hsps} = [];
		      $h->{_remove} = 1;
		  }
                  push(@keepers, $h);
                  next;
	      }
	      
	      #filtering
	      next unless(_keep_filt($h, $params));

              #keep
	      if(!$depth){
                  push(@keepers, $h);
                  next;
	      }
              else{
                  push(@$clusters, $h);
                  $i++;
              }
          }

          #parse, cluster, and filter at same time
	  next if(! $depth || $i < $depth*20);
          $i = 0;

          $clusters = cluster::shadow_cluster($depth, $clusters);
      }
  }

   #final filter
   if($depth && @$clusters){
       $clusters = cluster::shadow_cluster($depth, $clusters);
       push(@keepers, @{GI::flatten($clusters)});
   }
   
   my $end     = @keepers;
   my $deleted = $start - $end;
   print STDERR "deleted:$deleted hits\n" unless $main::quiet;
   
   return \@keepers;
}
#-----------------------------------------------------------------------------
sub _keep_filt {
    my $h = shift;
    my $params = shift;

    #filtering
    #next unless PhatHit_utils::is_contigous($h);
    if(defined($params->{significance})){
	my $sig_thr = $params->{significance};
	my $significance = $h->significance();
	$significance = "1".$significance if  $significance =~ /^e/;
	$significance = 0                 if  $significance =~ /0\.$/;
	return 0 if($significance > $sig_thr);
    }
    if(defined($params->{percov})){
	my $pcov = $params->{percov} || 0;
	return 0 if($h->pAh < $pcov/2);
    }
    if(defined($params->{percid})){
	my $pid = $params->{percid} || 0;
	return 0 if($h->hsp('best')->frac_identical() < $pid &&
		    $h->frac_identical < $pid);
    }

    return 1;
}
#-----------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::tblastx::AutoLoader called for: ",
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


