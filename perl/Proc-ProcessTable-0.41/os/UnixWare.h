#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/procfs.h>
#include <sys/statvfs.h>
#include <sys/types.h>
#include <sys/fault.h>
#include <sys/syscall.h>
#include <sys/param.h>

/****************************************/
/* Process state strings that we return */
/****************************************/
#define SLEEP  "sleep"
#define WAIT   "wait"       
#define RUN    "run"
#define IDLE   "idle"
#define STOP   "stop"
#define ONPROC "onprocessor" 

/* Solaris is an all-or-nothing deal, all this stuff comes out of 
   one structure, so we don't need to dick around with the format much */
static char Format[] = "iiiiiiiilllllllslsis";

/* Mapping of field to type */
static char* Fields[] = {

  "uid",
  "gid",
  "pid",
  "ppid",
  "pgrp",
  "sess",

  "priority",
  "nice",
  "ttynum",
  "flags",
  "time",

  "timensec",

  "size",
  "rss",
  "wchan",

  "fname",

  "start",

  "state",
  "onpro",
  "cmndline"
};









