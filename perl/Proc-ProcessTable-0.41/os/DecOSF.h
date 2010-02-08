#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/proc.h> /* needed for process state constants */
#include <sys/table.h> /* needed for process state constants */
#include <sys/statvfs.h>
#include <sys/types.h>
#include <sys/procfs.h>

#ifdef i386
#undef SP
#define SP sp
#endif

/* How many processes to get at once from table() */
#define PROCCNT 64

/****************************************/
/* Process state strings that we return */
/****************************************/
#define SLEEP  "sleep"
#define WAIT   "wait"       
#define RUN    "run"
#define IDLE   "idle"
#define ZOMBIE "defunct"
#define STOP   "stop"
#define ONPROC "onprocessor" 

/* Digital Unix is an all-or-nothing deal, all this stuff comes out of 
   one structure, so we don't need to dick around with the format much */
static char Format[] = "iiiiiiViiiilllllllssss";

/* Mapping of field to type */
static char* Fields[] = {
  "uid",
  "gid",
  "euid",
  "egid",
  "suid",
  "sgid",
  "groups",
  "pid",
  "ppid",
  "pgrp",
  "sess",
  "priority",
  "ttynum",
  "flags",
  "time",
  "size",
  "rss",
  "start",
  "fname",
  "pctcpu",
  "state",
  "cmndline"
};

