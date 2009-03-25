eval {
    require DBI;
    require DBI::DBD;
    die "Too old" unless $DBI::VERSION >= 1.03;
};
use ExtUtils::MakeMaker;
use Config;
use strict;
use FindBin;

WriteMakefile(
    'NAME'          => 'DBD::SQLite',
    'VERSION_FROM'  => 'lib/DBD/SQLite.pm', # finds $VERSION
    'PREREQ_PM'     => {DBI => 1.21}, # e.g., Module::Name => 1.1
    'OBJECT'        => '$(O_FILES)',
    'INC'           => '-I$(DBI_INSTARCH_DIR)',
    'OPTIMIZE'      => "-O6 -DNDEBUG=1 -DSQLITE_PTR_SZ=$Config{ptrsize}",
    'clean'         => { FILES => 'SQLite.xsi config.h' },
    'DESTDIR' => "$FindBin::Bin/../",
    'INSTALLDIRS' => "perl",
    'INSTALLMAN1DIR' => "none",
    'INSTALLMAN3DIR' => "none",
    'INSTALLARCHLIB' => "lib/",
    'INSTALLPRIVLIB' => "lib/",
    'INSTALLBIN'     => "bin/",
    'INSTALLSCRIPT'  => "bin/",
);

package MY;
sub postamble { 
  DBI::DBD::dbd_postamble(@_);
}
sub libscan {
    my ($self, $path) = @_;
    ($path =~ m/\~$/) ? undef : $path;
}
