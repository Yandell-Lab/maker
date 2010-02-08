use ExtUtils::MakeMaker;
use FindBin;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

sub MY::c_o {
  package MY; # so that "SUPER" works right
  my $inherited = shift->SUPER::c_o(@_);
  $inherited =~ s/\$\*\.c/\$\(C_FILES\)/; 
  $inherited;
}

WriteMakefile(
    'NAME'	=> 'Proc::ProcessTable',
    'VERSION_FROM' => 'ProcessTable.pm', # finds $VERSION
    'LIBS'	=> [''],   # e.g., '-lm' 
    'DEFINE'	=> '',     # e.g., '-DHAVE_SOMETHING' 
    'INC'	=> '',     # e.g., '-I/usr/include/other' 
    'LDFROM'    => '$(O_FILES)',
    'OBJECT'    => 'ProcessTable.o OS.o',
    'PREREQ_PM'    => { 'File::Find' => 0, 'Storable' => 0 },
    'clean'     => { FILES => 'OS.c' },
    'dist' => { COMPRESS=>"gzip", SUFFIX=>"gz" },
    'DESTDIR' => "$FindBin::Bin/../",
    'INSTALLDIRS' => "perl",
    'INSTALLMAN1DIR' => "none",
    'INSTALLMAN3DIR' => "none",
    'INSTALLARCHLIB' => "lib/",
    'INSTALLPRIVLIB' => "lib/",
    'INSTALLBIN'     => "bin/",
    'INSTALLSCRIPT'  => "bin/",
);
