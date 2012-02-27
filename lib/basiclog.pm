#-------------------------------------------------------------------------------
#------                           runlog                               ---------
#-------------------------------------------------------------------------------
package basiclog;
use strict;
use vars qw(@ISA @EXPORT $VERSION %SEEN);
use Exporter;

@ISA = qw();

#-------------------------------------------------------------------------------
#------------------------------- Methods ---------------------------------------
#-------------------------------------------------------------------------------
sub new {
    my $self = {};
    my $class = shift;

    bless ($self, $class);

    my $file = shift || "run.log";
    $self->{file} = $file;

    open(my $OUT, "> $file");
    close($OUT);

    return $self;
}
#-------------------------------------------------------------------------------
sub add_entry {
   my $self = shift;

   my @vals = @_;
   
   foreach my $val (@vals){
       chomp($val);
   }

   my $log_file = $self->{file};

   open(LOG, ">> $log_file");
   print LOG join("\t", @vals)."\n";
   close(LOG);
}
#-------------------------------------------------------------------------------
sub continue_flag {
   my $self = shift;
   my $id = shift;

   $SEEN{$id}++;

   my $flag = 0;
   $flag = 1 if($SEEN{$id} <= 2); #retry only once

   return $flag;
}
#-------------------------------------------------------------------------------
1;
