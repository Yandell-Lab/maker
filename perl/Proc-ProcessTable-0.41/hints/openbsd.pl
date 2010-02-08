symlink "os/OpenBSD.c", "OS.c" || die "Could not link os/OpenBSD.c to os/OS.c\n";
$self->{LIBS} = ['-lkvm'];
