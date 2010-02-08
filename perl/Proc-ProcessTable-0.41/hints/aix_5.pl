symlink "os/aix_getprocs.c", "OS.c" || die "Could not link os/aix_getprocs.c to os/OS.c\n";

$self->{LIBS} = "-lodm -lcfg";
