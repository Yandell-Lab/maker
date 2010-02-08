#include <sys/types.h>
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* We need to pass in a cap for ignore, lower for store on object */
/* We can just lc these! */
static char Defaultformat[] = "iiiiiiisssssssissiiiiiii";

/* Mapping of field to type */
static char* Fields[] = {
    "pid",
    "ppid",
    "uid",
    "euid",
    "gid",
    "pgrp",
    "sess",

    "flags",
    "sflags",
    "start",
    "time",
    "wchan",
    "state",

    "ttydev",
    "ttynum",
    "fname",
    "cmndline",
    "priority",
    "nice",
    "vmsize",
    "rssize",
    "tsize",
    "dsize",
    "ssize"
};

/*

     pid_t   ki_pid;                 // Process identifier
     pid_t   ki_ppid;                // parent process id 
     pid_t   ki_pgid;                // process group id  

     udev_t  ki_tdev;                // controlling tty dev
     pid_t   ki_sid;                 // Process session ID  

     uid_t   ki_uid;                 // effective user id
     uid_t   ki_ruid;                // Real user id
     gid_t   ki_rgid;                // Real group id

     vm_size_t ki_size;              // virtual size
     segsz_t ki_rssize;              // current resident set size in pages
     segsz_t ki_tsize;               // text size (pages) XXX
     segsz_t ki_dsize;               // data size (pages) XXX
     segsz_t ki_ssize;               // stack size (pages)   

     u_int64_t ki_runtime;           // Real time in microsec  
     struct  timeval ki_start;       // starting time

     long    ki_sflag;               // PS_* flags
     long    ki_flag;                // P_* flags
     char    ki_stat;                // S* process status  

     char    ki_wmesg[WMESGLEN+1];   // wchan message
     char    ki_comm[COMMLEN+1];     // command name 

     char    ki_nice;                // Process "nice" value
     struct  priority ki_pri;        // process priority

*/

