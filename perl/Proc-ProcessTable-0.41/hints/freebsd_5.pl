symlink "os/FreeBSD-kvm.c", "OS.c" || die "Could not link os/FreeBSD-kvm.c to os/OS.c\n";
$self->{LIBS} = ['-lkvm'];
