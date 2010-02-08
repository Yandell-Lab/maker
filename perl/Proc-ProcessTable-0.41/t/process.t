# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

use strict;
use Test;
BEGIN { plan tests => 3 }

# check wether ProcProcessTable is there
use Proc::ProcessTable; 

# Test code

$SIG{CHLD} = sub{wait;};
my ( $got, $field );

my $t = new Proc::ProcessTable;

# Is there a process called cron
foreach $got ( @{$t->table} )
{
#  print STDERR $got->pid, "  ", $got->fname, "\n";
  print STDERR "--------------------------------\n";
  foreach $field ($t->fields){
    my $v = $got->{$field};
    if (ref($v) eq "ARRAY")
    {
      $v = "\"" . join ("\",\"", @$v) . "\"";
    }
    print STDERR $field, ":  ", (defined $v ? $v : "<undef>"), "\n";
  }

}

# fork a child process
my $child_pid = fork;

if ( $child_pid )
{ 
  # parent, fork returned pod of the child process
  foreach $got ( @{$t->table} )
  {
    if( $got->pid == $child_pid )
    {
      ok(1);    # pid of the child process found

      if( $got->kill(9) )
      {
        ok(1);
      }
      else
      {
        ok(0);
        kill 9, $child_pid;
        exit -1;
      }

      sleep 2;

      # the child process should be dead now
      foreach $got ( @{$t->table} )
      {
        if( $got->pid == $child_pid )
        {
          ok(0);
          kill 9, $child_pid;
          exit -1;
        }
      }
      ok(1);
      exit 0;
    }
  }

  # pid of child was never found
  ok(0);
  exit -1;

}
else 
{ 
  # child, fork returned 0
  # child process will be killed soon
  sleep 10000;
}
