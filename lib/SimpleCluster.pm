#------------------------------------------------------------------------------
#----                            SimpleCluster.pm                          ---- 
#------------------------------------------------------------------------------
package SimpleCluster;
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
use PostData;
use FileHandle;

@ISA = qw(
       );
#------------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------------
#------------------------------------------------------------------------------
sub new {
        my $self  = shift;
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub pairs {
        my $matrix    = shift;
        my $threshold = shift;

	print STDERR "getting Pairs\n";
        my @pairs;
         while (my $iUid = each %{$matrix}){
                 while (my $jUid = each %{$matrix->{$iUid}}) {
                        next if $iUid eq $jUid;
                        my $Sij = $matrix->{$iUid}->{$jUid};
                        next unless defined($Sij);
                        next if (defined($threshold) && $Sij < $threshold);
                        next if $Sij == 0;
                        push(@pairs, [$iUid, $jUid]);
                }
        }
        return \@pairs;
}
#---------------------------------------------------------------------------
sub singleLinkageClusters {
        my $pairs     = shift;

        my $cId = -1;
        my %cIds;
        my %rMap;
        print STDERR "doing single linkage clustering\n";
        foreach my $pair (@{$pairs}){
                my ($mUidI, $mUidJ) = @{$pair};
                if    (!defined($cIds{$mUidI}) && !defined($cIds{$mUidJ})){
                        $cId++;
                        $cIds{$mUidI} = $cId;
                        $cIds{$mUidJ} = $cId;
                        push(@{$rMap{$cIds{$mUidI}}}, $mUidI, $mUidJ);
                }
                elsif (defined($cIds{$mUidI}) && !defined($cIds{$mUidJ})){
                        my $cId = $cIds{$mUidI};
                        $cIds{$mUidJ} = $cId;
                        push(@{$rMap{$cId}}, $mUidJ);
                }
                elsif (!defined($cIds{$mUidI}) && defined($cIds{$mUidJ})){
                        my $cId = $cIds{$mUidJ};
                        $cIds{$mUidI} = $cId;
                        push(@{$rMap{$cId}}, $mUidI);
                }
                elsif (defined($cIds{$mUidI}) && defined($cIds{$mUidJ})){
                        next if ($cIds{$mUidI} == $cIds{$mUidJ} 
			     && $cIds{$mUidI}  eq $cIds{$mUidJ});

                        my @copy;
                        my $cIdI = $cIds{$mUidI};
                        my $cIdJ = $cIds{$mUidJ};
                        foreach my $mUid (@{$rMap{$cIdJ}}){
				$cIds{$mUid}= $cIdI;
                                push(@copy, $mUid);
                        }
                        foreach my $mUid (@{$rMap{$cIdI}}){
                                push(@copy, $mUid);
                        }
                        $rMap{$cIdI} = \@copy;
			undef $rMap{$cIdJ};
                }
        }
        return \%rMap;
}
#----------------------------------------------------------------------------
sub joinAlternates {
        my $a          = shift;
        my $b          = shift;
        my $alternates = shift;
        unless (defined(@{$alternates})){        #-- start
                push(@{$alternates->[0]}, $a);
                push(@{$alternates->[0]}, $b);
                return $alternates;
        }
        my $swallow = 0;
        for (my $i = 0; $i < @{$alternates}; $i++){
                my %hash;
                my @newGroup;
                foreach my $f (@{$alternates->[$i]}){
                        $hash{$f->name} = 1;
                         push(@newGroup, $f);
                }
                if (!defined($hash{$a->name}) && defined($hash{$b->name})){
                        $swallow = 1;
                        push(@newGroup, $a);
                }
                elsif (defined($hash{$a->name}) && !defined($hash{$b->name})){
                        $swallow = 1;
                        push(@newGroup, $b);
                }
                elsif (!defined($hash{$a->name}) && !defined($hash{$b->name})){
                        #-- does not go in this group
                }
                elsif (defined($hash{$a->name}) && defined($hash{$b->name})){
                        $swallow = 1;
                        #-- just comming at it from the opposite direction
                }
                $alternates->[$i] = \@newGroup;
        }
        if ($swallow == 0){
                my $counts = @{$alternates};
                push(@{$alternates->[$counts]}, $a);
                push(@{$alternates->[$counts]}, $b);
        }
        return $alternates;
}
#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        print STDERR "SimpleCluster::AutoLoader called for: ",
              "\$self->$call","()\n";
        print STDERR "call to AutoLoader issued from: ", $caller, "\n";

        if ($arg){
                $self->{$call} = $arg;
        }
        else {
                return $self->{$call};
        }
}
#------------------------------------------------------------------------
1;


