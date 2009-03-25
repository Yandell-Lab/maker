use Test;
BEGIN { plan tests => 7 }
use DBI;
my $dbh = DBI->connect("dbi:SQLite:dbname=foo", "", "", { RaiseError => 1 });
ok($dbh);
my $sth = $dbh->prepare("SELECT * FROM f");
ok($sth);
ok($sth->execute);
my $row = $sth->fetch;
ok($row);
ok(@$row, 3);
print join(", ", @$row), "\n";
my $rows = $sth->execute;
ok($rows);
ok($sth->fetch);
$sth->finish;
$dbh->disconnect;
