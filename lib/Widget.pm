#----------------------------------------------------------------------------
#----                             Widget                                 ---- 
#----------------------------------------------------------------------------
package Widget;
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
        my $class = shift;
	my @args   = @_;

        my $self = {};
        bless $self;

	return $self;
}
#------------------------------------------------------------------------------
sub run {
	my $self  = shift;
	my $params = shift;
	die "run is an abstract method\n";
}
#-------------------------------------------------------------------------------
sub print_command {
	my $self    = shift;
	my $command = shift;

	print STDERR "#--------- command -------------#\n" unless ($main::quiet);
	print STDERR ref($self).":\n" unless ($main::quiet);
	if (defined($command)){
		print STDERR $command."\n" unless ($main::quiet);
	}
	else {
		print STDERR "executing default command...\n" unless ($main::quiet);
	}
	print STDERR "#-------------------------------#\n" unless ($main::quiet);
}
#-------------------------------------------------------------------------------
sub build_redirect {
   my $self = shift;
   my $command = shift;

   my $current = '';
   my $last = '';
   my $next2last = '';
   
   my $ignore = 0;
   my $quote = '';
   my $record = 0;
   my $redirect = 1;
   my $line;
   
   while ($$command =~ /(.)/g){
      $next2last = $last;
      $last = $current;
      $current = $1;
      
      if($current eq "\'" && $last ne "\\"){
	 if ($quote eq $current){
	    $ignore = 0;
	    $quote = '';
	 }
	 elsif($quote eq "\""){
	    #do nothing
	 }
	 else{
	    $ignore = 1;
	    $quote = $current;
	 }
      }
      elsif($current eq "\"" && $last ne "\\"){
	 if ($quote eq $current){
	    $ignore = 0;
	    $quote = '';
	 }
	 elsif($quote eq "\'"){
	    #do nothing
	 }
	 else{
	    $ignore = 1;
	    $quote = $current;
	 }
      }
      elsif($current eq "\>" && !$ignore && $last ne "\\" ){
	 if ($last =~ /\d/ && $next2last =~ /[\s\t]/){
	    if ($last eq '1'){#alternate STDOUT redirect
	       $record = 1;
	       $line .= $last if (!defined($line))
	    }
	    elsif($last eq '2'){#can't redirect STDERR
	       $redirect = 0;
	    }
	 }
	 elsif($last eq "\&" && $next2last ne "\\"){#can't redirect STDERR
	    $redirect = 0;
	 }
	 else{
	    $record = 1;
	 }
      }
      
      if($record){
	 $line .= $current;
      }
   }
   
   if ($redirect){
      if (defined $line){
	 $$command =~ s/$line$//;
      }
      
      $$command .= " 2>&1";
      $$command .= " " . $line if (defined($line));
   }
}
#-------------------------------------------------------------------------------
sub queryName {
	my $self = shift;
	my $name = shift;

	if    (defined($name)){
		$self->{queryName} = $name;
	}
	elsif (defined($self->{queryName})){
		return $self->{queryName};
	}
	else {
		my $file = $self->queryFastaFile();
		my ($name) = $file =~ /.*\/(\S+)$/;
		$self->{queryName} = $name;
		return $self->{queryName};
	}
}
#-------------------------------------------------------------------------------
sub queryFasta {
	my $self = shift;

	my $fh = new FileHandle();
	   $fh->open($self->queryFastaFile);

	my $fasta = '';
	while (my $line .= <$fh>){ $fasta .= $line};
	$fh->close;
	return $fasta;	
}
#-------------------------------------------------------------------------------
sub parse {
	my $self   = shift;
	my $params = shift;

	die "parse is an abstract method\n";
}
#-------------------------------------------------------------------------------
#------------------------------ FUNCTIONS --------------------------------------
#-------------------------------------------------------------------------------
sub toFasta {
        my $def = shift;
        my $seq = shift;

        my $fasta = $def."\n";
        for (my $i=0; $i< length($seq);$i+=60){
                $fasta .= substr($seq, $i, 60). "\n";
        }
        return $fasta;
}

#-------------------------------------------------------------------------------
sub AUTOLOAD {
        my ($self, $arg) = @_;

        my $caller = caller();
        use vars qw($AUTOLOAD);
        my ($call) = $AUTOLOAD =~/.*\:\:(\w+)$/;
        $call =~/DESTROY/ && return;

        #print STDERR "Widget::AutoLoader called for: ",
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


