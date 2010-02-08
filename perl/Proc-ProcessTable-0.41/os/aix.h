/***************************************************************************/
/* Includes                                                                */
/***************************************************************************/

#include <stdlib.h>
#include <procinfo.h>
#include <unistd.h>
#include <string.h>
#include <utmp.h>
#include <time.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/proc.h>
#include <sys/time.h>

#include <odmi.h>
#include <sys/cfgodm.h>
#include <sys/cfgdb.h>


/***************************************************************************/
/* Defines                                                                 */
/***************************************************************************/

/*
 * Process states
 *
 */

#define SLEEP  "sleep"
#define WAIT   "wait"       
#define RUN    "run"
#define IDLE   "idle"
#define ZOMBIE "defunct"
#define STOP   "stop"
#define UWAIT  "uwait"
#define ACTIVE "active"    
/*****************************************************************************/
/* Copyright (c) 1998, David Paquet. All rights reserved.                    */
/* This file is free software; you can redistribute it and/or modify it      */
/* under the same terms as Perl itself.                                      */
/*****************************************************************************/

/*
 * Arbitrary constants
 * 
 */

/* Grab the maximum argument length */
#include <sys/limits.h>

#define MAX_PROCS	1024	/* Pretty overloaded isn't it ? */
#define MAXARGLN	ARG_MAX


/*
 * Some landmarks ...
 *
 */

#define F_STAT       9
#define F_TTY       27
#define F_PRM       32
#define F_COMM      33
#define F_FLAST     33


/***************************************************************************/
/* Globals                                                                 */
/***************************************************************************/

static unsigned long long Sysmem;
static int  PageSize;
static int  ProcessNumber;

static char Fullformat[]    = "llllllllsslsllllllllllllllllllllsss";
static char Zombformat[]    = "lllllllllslslllllll";


static char* ZombFields[] = {
  "pid",
  "ppid",
  "sess",
  "pgrp",
  "uid",
  "suid",
  "priority",
  "nice",
  "pctcpu",
  "stat",
  "flags",
  "wchan",
  "wtype",
  "adspace",
  "majflt",
  "minflt",
  "utime",
  "stime",
  "size" };


static char* FullFields[] = {
  "pid",
  "ppid",
  "sess",
  "pgrp",
  "uid",
  "suid",
  "priority",
  "nice",
  "pctcpu",
  "stat",
  "flags",
  "wchan",
  "wtype",
  "adspace",
  "majflt",
  "minflt",
  /*    "utime", */    /* field valid for zombies only, see <procinfo.h> */
  /*    "stime", */    /* field valid for zombies only, see <procinfo.h> */
  "size",
  "luid",
  "euid",
  "gid",
  "start",
  "utime",
  "stime",
  "cutime",
  "cstime",
  "tsize",
  "ttyp",
  "ttynum",
  "ttympx",
  "drss",
  "trss",
  "dvm",
  "pctmem",
  "comm",
  "cmndline"
};


