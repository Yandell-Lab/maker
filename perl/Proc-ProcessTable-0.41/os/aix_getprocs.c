/*
 * Copyright (c) 2002, Target Corporation.  All Rights Reserved.
 * This file is free software; you can redistribute it and/or modify
 * it under the same terms as Perl itself.
 *
 * Author: James FitzGibbon <james.fitzgibbon@target.com>
 *
 * based on aix.c distributed with Proc::ProcessTable v0.35, which
 * is Copyright (c) 1998, David Paquet.
 *
 */

#include <odmi.h>
#include <procinfo.h>
#include <stdlib.h>
#include <sys/cfgodm.h>
#include <sys/types.h>

#include "EXTERN.h"
#include "perl.h"
#include "os/aix_getprocs.h"


/* convert a struct timeval to seconds *
 * TVALU_* presumes that tv_usec is microseconds
 * TVALN_* presumes that tv_usec is nanoseconds (struct procsinfo64)
 */
#define	TVALU_TO_SEC(x)		(x.tv_sec + ((double)x.tv_usec / 1000000.0))
#define	TVALN_TO_SEC(x)		(x.tv_sec + ((double)x.tv_usec / 1000000000.0))


void bless_procs(struct procsinfo64 *, int);


static unsigned long long memory;
static int pagesize = 0;
static int ncpus = 0;
static double now_time = 0.0;
static char format[F_LASTFIELD+1];


char *OS_initialize() {

	struct CuAt* odm_object;
	int num_fetched;

	/* get the number of processors online */
	ncpus = sysconf(_SC_NPROCESSORS_ONLN);
	if( ncpus == -1 ) {		/* sysconf error */
		ncpus = 1;
	}

	/* get the page size in bytes */
	pagesize = getpagesize();

	/* get the amount of physical memory */
	if( 0 != odm_initialize() ) {
		/* fprintf(stderr, "cannot initialize ODM in Proc::ProcessTable::OS_initialize (AIX)!\n"); */
		ppt_warn("cannot initialize ODM in Proc::ProcessTable::OS_initialize (AIX)!");
	} else {
		odm_object = (struct CuAt*)getattr("sys0", "realmem", 0, &num_fetched);
		memory = strtoull(odm_object->value, 0, 10);
		odm_terminate();
	}

    memory = memory * 1024;

	return NULL;

}


void OS_get_table() {

	struct procsinfo64 *procs = NULL;
	int index = 0;
	int fetched = 0;
	struct timeval now_tval;

	/* allocate memory to hold the procs */
/*	procs = New(0, procs, PROCS_TO_FETCH, struct procsinfo64); */
    procs = (struct procsinfo64 *)malloc(sizeof(struct procsinfo64) * PROCS_TO_FETCH);
	if(NULL == procs) {
		/* fprintf(stderr, "cannot allocate memory in Proc::ProcessTable::OS_get_table!\n"); */
		ppt_warn("cannot allocate memory in Proc::ProcessTable::OS_get_table!");
		return;
	}

	/* get current time of day */
	gettimeofday(&now_tval, 0);
	now_time = TVALU_TO_SEC(now_tval);

	/* keep on grabbing chunks of processes until getprocs returns a smaller
       block than we asked for */
	while( (fetched = getprocs(procs, sizeof(struct procsinfo64),
                               NULL, 0, &index, PROCS_TO_FETCH))
           >= PROCS_TO_FETCH) {

		bless_procs(procs, fetched);
	}

	/* bless the last block of procs */
	bless_procs(procs, fetched);

	/* release the memory */
	Safefree(procs);

	return;

}


void bless_procs(struct procsinfo64 *procs, int count)
{

	int index;
	char cmndline[ARG_MAX];
	char comm[ARG_MAX];
	char pctmem[PCT_LENGTH];
	char pctcpu[PCT_LENGTH];
	char state[STATE_LENGTH];
	char *c = NULL;
	struct procsinfo64 *curproc = NULL;
    struct procsinfo curproc_for_getargs;
	int zombie;
	int done;
	long utime, stime, cutime, cstime;

    /* initialize */
    Zero(&curproc_for_getargs, 1, struct procsinfo);

	for( index = 0, curproc = procs; index < count; index++, curproc = procs+index ) {

		/* reset char fields */
		memset(cmndline, 0, ARG_MAX);
		memset(comm, 0, ARG_MAX);
		memset(pctmem, 0, PCT_LENGTH);
		memset(pctcpu, 0, PCT_LENGTH);
		memset(state, 0, STATE_LENGTH);

		/* skip procs with a state of NONE */
		if( curproc->pi_state == SNONE ) {
			continue;
		}

		/* copy the format into something we can use */
		strcpy(format, Defaultformat);

		/* presume we are not a zombie */
		zombie = 0;

		/* set a descriptive name for the process state */
		switch( curproc->pi_state ) {

			case SSLEEP:
					strcpy(state, SLEEP);
					break;

			case SRUN:
					strcpy(state, RUN);
					break;

			case SIDL:
					strcpy(state, IDLE);
					break;

			case SZOMB:
					strcpy(state, ZOMBIE);
					zombie = 1;
					break;

			case SSTOP:
					strcpy(state, STOP);
					break;

			case SACTIVE:
					strcpy(state, ACTIVE);
					break;

			default:
					format[F_STAT] = 'S';
					break;

		}

		/* calc utime, stime, cutime, cstime */
		if( zombie ) {
			utime = curproc->pi_utime;
			stime = curproc->pi_stime;
		} else {
			utime = TVALN_TO_SEC(curproc->pi_ru.ru_utime);
			stime = TVALN_TO_SEC(curproc->pi_ru.ru_stime);
			cutime = TVALN_TO_SEC(curproc->pi_cru.ru_utime);
			cstime = TVALN_TO_SEC(curproc->pi_cru.ru_stime);
		}

		/* calc pctcpu */
		sprintf(pctcpu, "%3.2f", ((utime + stime) * 100 /
                                  ( now_time - curproc->pi_start )) / ncpus );

		/* calc pctmem */
		if( memory == 0 ) {
			format[F_PCTMEM] = 'S';
		} else {
			sprintf(pctmem, "%3.2f", (curproc->pi_drss + curproc->pi_trss) * pagesize * 100 / (float)memory);
		}

		/* determine comm & cmndline */
		if( (curproc->pi_flags & SKPROC) == SKPROC ) {
			if( curproc->pi_pid == 0 ) {
				snprintf(comm, ARG_MAX, "kproc (swapper)");
				snprintf(cmndline, ARG_MAX, "kproc (swapper)");
			} else {
				snprintf(comm, ARG_MAX, "kproc (%s)", curproc->pi_comm);
				snprintf(cmndline, ARG_MAX, "kproc (%s)", curproc->pi_comm);
			}
		} else {
			snprintf(comm, ARG_MAX, "%s", curproc->pi_comm);
            curproc_for_getargs.pi_pid = curproc->pi_pid;
			if( getargs(&curproc_for_getargs, sizeof(struct procsinfo), cmndline, ARG_MAX) < 0 ) {
				snprintf(cmndline, ARG_MAX, "%s", curproc->pi_comm);
			} else {
				/* replace NUL characters in command line with spaces */
				c = cmndline;
				done = 0;
				while( ! done ) {
					if( *c == '\0' ) {
						if( *(c+1) == '\0' ) {
							done = 1;
						} else {
							*c++ = ' ';
						}
					} else {
						c++;
					}
				}
			}
		}

		/* sanity check: bless_into_procs doesn't like negative ints */
		if( curproc->pi_cred.cr_uid == -1 ) {
			format[F_EUID] = 'I';
			curproc->pi_cred.cr_uid = 0;
		}
		if( curproc->pi_ttyd == -1 ) {
			format[F_TTYNUM] = 'I';
			curproc->pi_ttyd = 0;
		}
		if( curproc->pi_ttympx == -1 ) {
			format[F_TTYMPX] = 'I';
			curproc->pi_ttympx = 0;
		}
		
		/* give it to Perl */
		bless_into_proc(format, Fields,
			curproc->pi_pid,		/* pid */
			curproc->pi_ppid,		/* ppid */
			curproc->pi_sid,		/* sess */
			curproc->pi_pgrp,		/* pgrp */
			curproc->pi_uid,		/* uid */
			curproc->pi_suid,		/* suid */
			curproc->pi_cred.cr_luid,	/* luid */
			curproc->pi_cred.cr_uid,	/* euid */
			curproc->pi_cred.cr_rgid,	/* gid */
			curproc->pi_cred.cr_gid,	/* egid */
			curproc->pi_pri,		/* priority */
			curproc->pi_nice,		/* nice */
			curproc->pi_thcount,	/* thcount */
			state,
			curproc->pi_flags,		/* flags */
			curproc->pi_flags2,		/* flags2 */
			curproc->pi_adspace,	/* adspace */
			curproc->pi_majflt,		/* majflt */
			curproc->pi_minflt,		/* minflt */
			(long)(utime * 100),	/* utime */
			(long)(stime * 100),	/* stime */
			(long)(cutime * 100),	/* cutime */
			(long)(cstime * 100),	/* cstime */
			curproc->pi_start,		/* start */
			curproc->pi_size,		/* size */
			curproc->pi_tsize,		/* tsize */
			curproc->pi_ttyp,		/* ttyp */
			curproc->pi_ttyd,		/* ttynum */
			curproc->pi_ttympx,		/* ttympx */
			curproc->pi_drss,		/* drss */
			curproc->pi_trss,		/* trss */
			curproc->pi_dvm,		/* dvm */
			pctmem,					/* pctmem */
			pctcpu,					/* pctcpu */
			comm,					/* comm */
			cmndline				/* cmndline */
		);
	}

}


/*
 * EOF
 */
