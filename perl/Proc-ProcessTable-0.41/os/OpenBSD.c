/*
 * Copyright (c) 2006, William Yodlowsky <bsd@openbsd.rutgers.edu>
 * This file is free software; it can be modified and/or
 * redistributed under the same terms as Perl itself
 */

#include <sys/types.h>
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/proc.h>
#include <sys/stat.h>
#include <sys/sysctl.h>
#include <sys/syslimits.h>
#include <sys/user.h>

#include <ctype.h>
#include <fcntl.h>
#include <kvm.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* No need for the procstat structure since we don't use /proc */

/* We need to pass in a cap for ignore, lower for store on object */
/* We can just lc these! */
static char Defaultformat[] = "liiiiiiiiiiiissss";

/* Mapping of field to type */
static char* Fields[] = {
	"ttynum",
	"uid",
	"gid",
	"euid",
	"egid",
	"pid",
	"ppid",
	"pgrp",
	"sess",
	"time",
	"utime",
	"stime",
	"start",
	"fname",
	"state",
	"ttydev",
	"cmndline"
};
#define F_LASTFIELD 16

/* Set up simple bounds checking */
#define STRLCPY(num,targ,src) if (strlcpy(targ,src,sizeof(targ)) >= sizeof(targ)) \
	ppt_croak("call:%d source string is too big to copy into buffer", num);

#define STRLCAT(num,targ,src) if (strlcat(targ,src,sizeof(targ)) >= sizeof(targ)) \
	ppt_croak("call:%d source string is too big to append to buffer", num);

extern void bless_into_proc(char* format, char** fields, ...);

char* OS_initialize() {
	return(NULL);
}

void OS_get_table() {
	kvm_t *kd;
	char errbuf[_POSIX2_LINE_MAX];
	struct kinfo_proc2 *procs;
	int count;
	int i, argcount;
	int ttynum;
	long start;
	char *ttydev;
	char cmndline[ARG_MAX+1];
	char **pargv;

	/* for bless_into_proc */
	static char format[F_LASTFIELD + 2];
	char state[20];

	/* open the kvm interface */
	if ((kd = kvm_open(NULL, NULL, NULL, KVM_NO_FILES, errbuf)) == NULL) {
		ppt_croak("kvm_open: %s", errbuf);
	}

	/* get processes */
	if ((procs = kvm_getproc2(kd, KERN_PROC_ALL, 0, sizeof(*procs), &count)) == NULL) {
		kvm_close(kd);
		ppt_croak("kvm_getproc2: %s", kvm_geterr(kd));
	}

	/* bless_into_proc each process's information */
	for (i=0; i < count; i++) {
		STRLCPY(1,format,Defaultformat);

		/* get ttydev */
		ttynum=procs[i].p_tdev;
		ttydev=devname(ttynum, S_IFCHR);
		if (ttydev == NULL) ttydev = "??";

		/* process state */
		switch (procs[i].p_stat) {
			case SIDL:
				STRLCPY(2,state,"IDLE");
				break;
			case SRUN:
				STRLCPY(3,state,"RUN");
				break;
			case SSLEEP:
				STRLCPY(4,state,"SLEEP");
				break;
			case SSTOP:
				STRLCPY(5,state,"STOP");
				break;
			case SZOMB:
				STRLCPY(6,state,"ZOMBIE");
				break;
			default:
				STRLCPY(7,state,"UNKNOWN");
				break;
		}

		/* arguments */
		cmndline[0] = NULL;
		pargv = kvm_getargv2(kd, (const struct kinfo_proc2 *) &(procs[i]), 0);
		if (pargv) {
			argcount = 0;
			while (pargv[argcount] && strlen(cmndline) <= ARG_MAX) {
				STRLCAT(1,cmndline,pargv[argcount]);
				STRLCAT(2,cmndline," ");
				argcount++;
			}
		}

		/* everything else is straightforward, bless the lot */
		bless_into_proc( format,
			Fields,
			ttynum,
			procs[i].p_ruid,
			procs[i].p_rgid,
			procs[i].p_uid,
			procs[i].p_gid,
			procs[i].p_pid,
			procs[i].p_ppid,
			procs[i].p__pgid,
			procs[i].p_sid,
			procs[i].p_rtime_sec,
			procs[i].p_uutime_sec,
			procs[i].p_ustime_sec,
			procs[i].p_ustart_sec,
			procs[i].p_comm,
			state,
			ttydev,
			cmndline
		);
	}
	if (kd) {
		kvm_close(kd);
	}
}
