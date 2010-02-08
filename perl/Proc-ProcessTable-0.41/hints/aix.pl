symlink "os/aix.c", "OS.c" || die "Could not link os/aix.c to os/OS.c\n";

$self->{LIBS} = "-lodm -lcfg";
