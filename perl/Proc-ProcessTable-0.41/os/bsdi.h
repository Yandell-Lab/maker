#include <sys/types.h>
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/stat.h>
#include <sys/proc.h>
#include <sys/sysctl.h>
#include <fcntl.h>
#include <kvm.h>
#include <limits.h>
#include <unistd.h>
#include <sys/user.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Grab the maximum argument length */
#include <sys/syslimits.h>

#define MAXARGLN	ARG_MAX

/**********************************************************/
/* Copyright (c) 1999, Magic Software Development, Inc.   */ 
/* Author:  Sean Ray Eskins. <sean@gilasoft.com>          */
/* This file is free software; it can be modified and/or  */
/* redistributed under the same terms as Perl itself.     */
/**********************************************************/

/**********************************/
/* Process state names we return  */
/**********************************/
#define SLEEP   "SLEEP"
#define RUN     "RUN"
#define IDLE    "IDLE"
#define STOP    "STOP"
#define ZOMBIE  "ZOMBIE"
#define UNKNOWN "UNKNOWN"

/* We are using kinfo_proc structure and kvm_getprocs to access */
/* processes from the kernel, so we don't need the procstat structure */
/* present in perl modules for other operating systems (ie: we */
/* are not accessing the processes from a proc directory) */
/* All processes are accessed straight from the kernel, so we don't */
/* need to mess around with the format much. */

/*******************************************************/
/* Some scripts needed to round off to top two digits  */
/* of "miliseconds", so that only the first two digits */
/* are printed to the m:ss.[miliseconds] time field of */
/* Fields[] below.                                     */
/*******************************************************/


/* power(int base, int pow) -- takes integer base to power pow. . . */ 
/*                             for use with roundit(int integer). */

int power(int base, int pow) {
  int count, result;
  count=1;
  result=1;
  while (count <= pow) {
    result=base*result;
    count++;
  }
  return(result);
} 


/* roundit(int) -- rounds a 4-, 5-, or 6-digit number to its highest */
/*                 two places. */


int roundit(int integer) {
  char string[10];
  int length;
  int newnum;
  int cutoff;
  int pow10;
  sprintf(string, "%d", integer);
  length=strlen(string);
  cutoff=length - 2;
  pow10=power(10, cutoff);
  newnum=integer/pow10;
  return(newnum); 
}
  


/* We need to pass in a cap for ignore, lower for store on object */
/* We can just lc these! */
static char Defaultformat[] = "iiiiiisissssssls";

/* Mapping of field to type */
static char* Fields[] = {
  "uid",                   /* ruid from kernel */
#define F_UID 0 

  "gid",                   /* rgid from kernel */
#define F_GID 1

  "pid",
#define F_PID 2

  "ppid",
#define F_PPID 3

  "pgrp",
#define F_PGRP 4

  "priority", 
#define F_PRIORITY 5

  "sess",
#define F_SESSION 6

  "leader",
#define F_LEADER 7

  "time",
#define F_TIME 8              

  "wchan",
#define F_WCHAN 9

  "fname",                /* command name (without args) */
#define F_NAME  10

  "state",
#define F_STATE 11

  "started",
#define F_STARTED 12

  "ttydev",
#define F_TTYDEV 13

  "ttynum",
#define F_TTYNUM 14
  "cmndline"
#define F_CMNDLINE 15

#define F_LASTFIELD 15
};



