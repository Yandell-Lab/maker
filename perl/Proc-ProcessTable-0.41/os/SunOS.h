/*
 * Copyright (c) 2001 by Shawn A. Clifford <shawn.a.clifford@lmco.com>
 * This file may be distributed under the same terms as Perl.
 *
 * Modification History:
 *
 * Who	When		Description
 * ---	----------	--------------------------------------------
 * SAC	30July2001	Original code
 */

#include <stdio.h>
#include <kvm.h>
#include <fcntl.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/proc.h>
#include <sys/user.h>
#include <sys/session.h>
#include <sys/sysmacros.h>
#include <sys/limits.h>
#include <errno.h>
#include <malloc.h>

#ifndef perror
   void perror(char *);
#endif
#ifndef printf
   int printf(char *, ...);
#endif
#ifndef bzero
   void bzero(char *, int);
#endif
#ifndef strncat
   char *strncat(char *, char *, int);
#endif
#ifndef strcat
   char *strcat(char *, char*);
#endif

extern void bless_into_proc(char* format, char** fields, ...);

static char *Format = "iiiiiiiiiiiiissllsiii";

static char *Fields[] = {
    "uid",    /* real uid */
#             define F_UID		0
    "gid",    /* real gid */
#             define F_GID		1        
    "euid",   /* effective uid */
#             define F_EID		2
    "egid",   /* effective gid */
#             define F_EGID		3
    "pid",    /* process id */
#             define F_PID		4
    "ppid",   /* parent pid */
#             define F_PPID		5
    "pgrp",   /* process group leader (pid) */
#             define F_PGRP		6
    "priority",/* priority, negative is high */
#             define F_PRIORITY		7
    "flags",  /* process flags */
#             define F_FLAGS		8
    "size",   /* data + stack (in KBytes) */
#             define F_SIZE		9
    "rss",    /* resident set size (in KBytes) */
#             define F_RSS		10
    "nice",   /* nice for cpu usage */
#             define F_NICE		11
    "time",   /* seconds resident */
#             define F_TIME		12
    "fname",  /* file name of running image */
#             define F_FNAME		13
    "cmndline",/* entire command line */
#             define F_CMNDLINE		14
    "cpticks",/* ticks of cpu time, for pctcpu */
#             define F_CPTICKS		15
    "pctcpu", /* (decayed) %cpu for this process */
#             define F_PCTCPU		16
    "state",  /* current run state (eg. sleep, wait, ..) */
#             define F_STATE		17
    "sess",   /* aka sid */
#             define F_SESS		18
    "sid",    /* session id */
#             define F_SID		19
    "ttynum"  /* minor device number for the tty */
#             define F_TTYNUM		20
#             define F_LASTFIELD	20
};

static char *States[] = { 
   "", "sleep", "wait", "run", "idle", "zombie", "stop"
};
