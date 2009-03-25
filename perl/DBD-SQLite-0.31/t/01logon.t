use Test;
BEGIN { plan tests => 3 }
use DBI;
my $dbh = DBI->connect("dbi:SQLite:dbname=foo", "", "");
ok($dbh);
ok($dbh->{sqlite_version});
ok($dbh->{sqlite_encoding});
print "sqlite_version=$dbh->{sqlite_version}, sqlite_encoding=$dbh->{sqlite_encoding}\n";
$dbh->disconnect;

