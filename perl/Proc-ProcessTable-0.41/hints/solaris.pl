symlink "os/Solaris.c", "OS.c" || die "Could not link os/Solaris.c to os/OS.c\n";

# Historical Note (please correct me if this is wrong!):
#
# A linux-like /proc filesystem was introduced with Solaris 2.6, with
# process information available in pseudo-files under pseudo-directories
# named for process ids in /proc. Previous to 2.6 Solaris used a
# half-baked ioctl-based /proc scheme. Previous to this (SunOS, maybe
# much earlier versions of Solaris?) process information was only
# available through /dev/kmem, which required root access (ps is setuid
# on SunOS-- yuck). All three schemes are supported and available in 2.6,
# but at least the ioctl-based /proc support is supposed to be going away.

# If this OS version supports the new /proc filesystem, use it; 
# otherwise default to ioctl-proc
`uname -r` =~ /^(\d+\.\d+)/;
if( $1 > 5.5 ){
    $self->{DEFINE} = $self->{DEFINE} . " -DPROC_FS";
}

# For reasons I don't understand, we have to turn off the large file
# environment flags in order to compile in the large file environment
#10/28/2002:
# Should not set CCLFAGS to empty.  Should only remove these two symbol while preserving others.
if( `/usr/bin/getconf LFS_CFLAGS` =~ /-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64/){
    $self->{CCFLAGS} = $Config{ccflags};
    $self->{CCFLAGS} =~ s/-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64//;
}
