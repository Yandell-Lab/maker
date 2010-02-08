/*-
 * This code relies heavily on the Darwin "ps" command, which is available
 * from Apple in the adv_cmds portion of the Darwin distribution. The portions
 * of this code which were included from that source are:
 *
 * Copyright (c) 1990, 1993, 1994
 *	The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * The portions of this code which were necessary to tie into the Perl
 * Proc::ProcessTable module are:
 *
 * Copyright (c) 2003, Thomas R. Wyant, III
 *
 * and may be reused under the same terms as Perl itself.
 */

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <mach/mach_init.h>
#include <mach/mach_port.h>
#include <mach/mach_traps.h>
#include <mach/mach_types.h>
#include <mach/shared_memory_server.h>
#include <mach/task.h>
#include <mach/thread_act.h>
#include <mach/time_value.h>
#include <mach/vm_map.h>
#include <stdlib.h>
#include <sys/proc.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/sysctl.h>
#include <kvm.h>
#include <unistd.h>


/* these are in sys/sysctl.h */

#define KI_PROC(ki) (&(ki)->ki_p->kp_proc)
#define KI_EPROC(ki) (&(ki)->ki_p->kp_eproc)
#define STATE_MAX       7

typedef struct thread_values {
	struct thread_basic_info tb;
	/* struct policy_infos	schedinfo; */
	union {
		struct policy_timeshare_info tshare;
		struct policy_rr_info rr;
		struct policy_fifo_info fifo;
	} schedinfo;
} thread_values_t;

struct usave {
	struct	timeval u_start;
	struct	rusage u_ru;
	struct	rusage u_cru;
	char	u_acflag;
	char	u_valid;
};

typedef struct kinfo {
	struct kinfo_proc *ki_p;	/* proc structure */
	struct usave ki_u;	/* interesting parts of user */
	char *ki_args;		/* exec args */
	char *ki_env;		/* environment */
	task_port_t task;
	int state;
	int cpu_usage;
	int curpri;
	int basepri;
	int swapped;
	struct task_basic_info tasks_info;
	struct task_thread_times_info times;
	/* struct policy_infos	schedinfo; */
	union {
		struct policy_timeshare_info tshare;
		struct policy_rr_info rr;
		struct policy_fifo_info fifo;
	} schedinfo;
	int	invalid_tinfo;
	int	thread_count;
	thread_port_array_t thread_list;
	thread_values_t *thval;
	int	invalid_thinfo;
} KINFO;


#define TIME_IN_MICROSECONDS 1

#ifdef TESTING

/* Note that FREE_BUFFERS is to be defined appropriately before
 * DIE_HORRIBLY is used.
 */
#define DIE_HORRIBLY(s) {FREE_BUFFERS perror (s); exit (1);}

#else /* def TESTING */

#define DIE_HORRIBLY(s) {FREE_BUFFERS perror (s); return;}

extern void bless_into_proc(char* format, char** fields, ...);

#ifdef TIME_IN_MICROSECONDS
static char *Format = "iiiiiiiiiiiiiiijjjllslsssss";
#else /* def TIME_IN_MICROSECONDS */
static char *Format = "iiiiiiiiiiiiiiilllllslsssss";
#endif /* def TIME_IN_MICROSECONDS */

static char *Fields[] = {
	"pid",		/* Process ID */
#				define F_PID	0
	"ppid",		/* Parent pid */
#				define F_PPID	1
	"pgrp",		/* Process group leader */
#				define F_PGRP	2
	"uid",		/* real uid */
#				define F_UID	3
	"gid",		/* real gid */
#				define F_GID	4
	"euid",		/* effective uid */
#				define F_EID	5
	"egid",		/* effective gid */
#				define F_EGID	6
	"suid",		/* saved uid */
#				define F_SUID	7
	"sgid",		/* saved gid */
#				define F_SGID	8
	"priority",	/* proiority */
#				define F_PRIORITY 9
	"size",		/* virtual size (Kbytes) */
#				define F_SIZE	10
	"rss",		/* resident set size (Kbytes) */
#				define F_RSS	11
	"flags",	/* process flags */
#				define F_FLAGS	12
	"nice",		/* nice for CPU usage */
#				define F_NICE	13
	"sess",		/* session pointer */
#				define F_SESS	14
	"time",		/* total time (system + user, centi- or microseconds, depending) */
#				define F_CPTICKS	15
	"stime",	/* system time (centi- or micrseconds, depending) */
#				define F_STIME		16
	"utime",	/* user time (centi- or microseconds, depending) */
#				define F_UTIME		17
	"start",	/* Start time in seconds since epoch */
#				define F_START		18
	"wchan",	/* Wait channel (addr of system call) */
#				define F_WCHAN		19
	"ttydev",
#				define F_TTYDEV		20
	"ttynum",	/* device number for tty, or -1 */
#				define F_TTYNUM		21
	"pctcpu",	/* Percent cpu as string representing float */
#				define F_PCTCPU		22
	"pctmem",	/* Percent memory as string representing float */
#				define F_PCTMEM		23
	"state",	/* current run state (e.g. sleep, wait ...) */
#				define F_STATE		24
    "cmndline",	/* entire command line */
#				define F_CMNDLINE	25
	"fname",	/* filename */
#				define F_FNAME		26
#				define F_LASTFIELD	26
};

#endif /* def TESTING */


static char *States[] = {
	"", "idle", "run", "sleep", "stop", "zombie"
};
