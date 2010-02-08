symlink "os/IRIX.c", "OS.c" || die "Could not link os/IRIX.c to OS.c\n";

# If this OS version supports the new /proc filesystem, use it; 
# otherwise default to ioctl-proc
#`uname -r` =~ /^(\d+\.\d+)/;
#if( $1 > 5.5 ){
#    $self->{DEFINE} = "-DPROC_FS";
#}

# I really hope we won't go beyond IRIX 999.999
`uname -r` =~ /^(\d+)\.(\d+)/;
$self->{DEFINE} = sprintf "-DIRIX_VERSION=0x%03d%03d", $1,$2;
