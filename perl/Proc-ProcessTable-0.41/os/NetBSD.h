#include <sys/types.h>
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/stat.h>
#include <ctype.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct procstat {
  char comm[MAXCOMLEN+1];
  int pid;
  int ppid;
  int pgid;
  int sid;
  int tdev_maj;
  int tdev_min;
  char flags[256]; /* XXX */
  int start;
  int start_mic;
  int utime;
  int utime_mic;
  int stime;
  int stime_mic;
  char wchan[256]; /* XXX */
  int euid;
  int ruid;
  int rgid;
  int egid;
  char groups[256]; /* XXX */
};

/* We need to pass in a cap for ignore, lower for store on object */
/* We can just lc these! */
static char Defaultformat[] = "iiiiiissssssiisssis";

/* Mapping of field to type */
static char* Fields[] = {
  "uid",
#define F_UID 0 

  "gid",
#define F_GID 1

  "pid",
#define F_PID 2

  "ppid",
#define F_PPID 3

  "pgrp",
#define F_PGRP 4

  "sess",
#define F_SESS 5

  "flags",
#define F_FLAGS 6

  "utime",
#define F_UTIME 7

  "stime",
#define F_STIME 8

  "time",
#define F_TIME 9

  "wchan",
#define F_WCHAN 10

  "start",
#define F_START 11

  "euid",
#define F_EUID 12

  "egid",
#define F_EGID 13

  "fname",
#define F_FNAME 14

  "state",
#define F_STATE 15

  "ttydev",
#define F_TTYDEV 16

  "ttynum",
#define F_TTYNUM 17

  "cmndline"
#define F_CMNDLINE 18

#define F_LASTFIELD 18
};
