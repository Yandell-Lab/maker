use Test;
BEGIN { plan tests => 9 }
use DBI;
my $dbh = DBI->connect("dbi:SQLite:dbname=foo", "", "");
ok($dbh);
my $sth = $dbh->prepare("INSERT INTO f VALUES (?, ?, ?)");
ok($sth);
ok(my $rows = $sth->execute("Fred", "Bloggs", "fred\@bloggs.com"));
ok($rows == 1);
ok($dbh->func('last_insert_rowid'));
ok($sth->execute("test", "test", "1"));
ok($sth->execute("test", "test", "2"));
ok($sth->execute("test", "test", "3"));
ok($dbh->do("delete from f where f1='test'") == 3);
$dbh->disconnect;
