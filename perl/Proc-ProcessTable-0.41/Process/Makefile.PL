use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'Proc::ProcessTable::Process',
    'VERSION_FROM' => 'Process.pm', # finds $VERSION
    'DESTDIR' => "$FindBin::Bin/../",
    'INSTALLDIRS' => "perl",
    'INSTALLMAN1DIR' => "none",
    'INSTALLMAN3DIR' => "none",
    'INSTALLARCHLIB' => "lib/",
    'INSTALLPRIVLIB' => "lib/",
    'INSTALLBIN'     => "bin/",
    'INSTALLSCRIPT'  => "bin/",
);
