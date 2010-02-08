symlink "os/SunOS.c", "OS.c" || die "Could not link os/SunOS.c to os/OS.c\n";
$self->{LIBS} = "-lkvm";
