/*
 * Copyright (c) 2002, Target Corporation.  All Rights Reserved.
 * This file is free software; you can redistribute it and/or modify
 * it under the same terms as Perl itself.
 *
 * Author: James FitzGibbon <james.fitzgibbon@target.com>
 *
 * based on aix.h distributed with Proc::ProcessTable v0.35, which
 * is Copyright (c) 1998, David Paquet.
 *
 */

 
/* Descriptive Process States */
#define     SLEEP       "sleep"
#define     WAIT        "wait"
#define     RUN         "run"
#define     IDLE        "idle"
#define     ZOMBIE      "defunct"
#define     STOP        "stop"
#define     UWAIT       "uwait"
#define     ACTIVE      "active"


/* How many processes to grab at a time
 *
 * internally, this module will need
 * sizeof(struct procsinfo64) * PROCS_TO_FETCH
 * memory while in ->table(), so you can change this constant
 * if you are tight for memory
 *
 */
#define     PROCS_TO_FETCH      1000


/* Length of the various string fields */
#define     PCT_LENGTH       7
#define     STATE_LENGTH        10


/* Format string */
static char Defaultformat[] = "iiiiiiiiiiiiisiiijjlllljjjiiijjjssss";


/* Mapping of field to position in format string */
static char* Fields[] = {
    "pid",      /* int */
#define F_PID 0

    "ppid",     /* int */
#define F_PPID 1

    "sess",     /* int */
#define F_SESS 2

    "pgrp",     /* int */
#define F_PGRP 3

    "uid",      /* int */
#define F_UID 4

    "suid",     /* int */
#define F_SUID 5

    "luid",     /* int */
#define F_LUID 6

    "euid",     /* int */
#define F_EUID 7

    "gid",      /* int */
#define F_GID 8

    "egid",     /* int */
#define F_EGID 9

    "priority", /* int */
#define F_PRIORITY 10

    "nice",     /* int */
#define F_NICE 11

    "thcount",  /* int */
#define F_THCOUNT 12

    "stat",     /* int -> string */
#define F_STAT 13

    "flags",    /* int */
#define F_FLAGS 14

    "flags2",   /* int */
#define F_FLAGS2 15

    "adspace",  /* long */
#define F_ADSPACE 16

    "majflt",   /* long */
#define F_MAJFLT 17

    "minflt",   /* long */
#define F_MINFLT 18

    "utime",    /* long */
#define F_UTIME 19

    "stime",    /* long */
#define F_STIME 20

    "cutime",   /* long */
#define F_CUTIME 21

    "cstime",   /* long */
#define F_CSTIME 22

    "start",    /* long */
#define F_START 23

    "size",     /* long */
#define F_SIZE 24

    "tsize",    /* long */
#define F_TSIZE 25

    "ttyp",     /* int */
#define F_TTYP 26

    "ttynum",   /* int */
#define F_TTYNUM 27

    "ttympx",   /* int */
#define F_TTYMPX 28

    "drss",     /* long */
#define F_DRSS 29

    "trss",     /* long */
#define F_TRSS 30

    "dvm",      /* long */
#define F_DVM 31

    "pctmem",   /* float -> string */
#define F_PCTMEM 32

    "pctcpu",   /* float -> string */
#define F_PCTCPU 33

    "comm",     /* string */
#define F_COMM 34

    "cmndline", /* string */
#define F_CMNDLINE 35
#define F_LASTFIELD 35

};

/*
 * EOF
 */
