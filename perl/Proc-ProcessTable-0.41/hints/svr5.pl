# Check this is UnixWare 7 or above.
`uname` =~ /UnixWare/ || die "Not a supported SVR5\n";

# If this OS version supports the new /proc filesystem, use it.
`uname -v` =~ /^([\d\.]+)$/;
unless( $1 >= 7 ){
    die "Unsupported OS Version: $1\n";
}

symlink "os/UnixWare.c", "OS.c" || die "Could not link os/UnixWare.c to os/OS.c\n";
$self->{DEFINE} = "-D_KMEMUSER";
