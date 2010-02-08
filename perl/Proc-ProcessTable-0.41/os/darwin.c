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
 * Copyright (c) 2003, 2004 by Thomas R. Wyant, III
 *
 * and may be reused under the same terms as Perl itself.
 */
 
#include "os/darwin.h"

/*
 * static void getproclline (KINFO *k, char ** command_name, int *cmdlen, int eflg);
 */
static void getproclline (KINFO *k, char ** command_name, int *cmdlen,
	int eflg, int show_args);
static int get_task_info (KINFO *ki);

int mempages = 0;

#ifdef TESTING

void OS_get_table (void);
char *OS_initialize (void);

int main (int argc, char **argv) {
OS_initialize ();
OS_get_table ();
exit (0);
}

#endif


char* OS_initialize(void) {

	size_t oldlen;
	int mib[2];

	oldlen = sizeof(mempages);
	mib[0] = CTL_HW;
	mib[1] = HW_PHYSMEM;

	sysctl(mib, 2, &mempages, &oldlen, NULL, 0);

	return NULL;
}

/* The appropriate FREE_BUFFERS definition for OS_get_table */

#define FREE_BUFFERS \
	{	if (kprocbuf != NULL) free (kprocbuf); \
		if (kinfo != NULL) free (kinfo); \
	}

void OS_get_table(void) {

	size_t bufSize = 0;
	char *command_name;
	int cmdlen;
	int i;
	KINFO *kinfo = NULL;
	struct kinfo_proc *kp;
	struct kinfo_proc *kprocbuf = NULL;
	int local_error=0;
	int mib[4] = { CTL_KERN, KERN_PROC, KERN_PROC_ALL, 0 };
	int nentries;
	size_t orig_bufSize = 0;
	char pctcpu[6];
	char pctmem[6];
	int retry_count = 0;
	int state;

	if (sysctl(mib, 4, NULL, &bufSize, NULL, 0) < 0)
		DIE_HORRIBLY ("Failure calling sysctl")

	if ((kprocbuf= kp = (struct kinfo_proc *)malloc(bufSize)) == NULL)
		DIE_HORRIBLY ("Memory allocation failure")

	retry_count = 0;
	orig_bufSize = bufSize;
	for(retry_count=0; ; retry_count++) {
	/* retry for transient errors due to load in the system */
		local_error = 0;
		bufSize = orig_bufSize;
		if ((local_error = sysctl(mib, 4, kp, &bufSize, NULL, 0)) < 0) {
			if (retry_count < 1000) {
			/* 1 sec back off */
				sleep(1);
				continue;
			}
			DIE_HORRIBLY ("Failure calling sysctl")
		} else if (local_error == 0) {
			break;
		}
	/* 1 sec back off */
	sleep(1);
	}

	/* This has to be after the second sysctl since the bufSize
		may have changed.  */
	nentries = bufSize/ sizeof(struct kinfo_proc);

	if ((kinfo = malloc(nentries * sizeof(KINFO))) == NULL)
		DIE_HORRIBLY ("Memory allocation failure")
	memset(kinfo, 0, (nentries * sizeof(*kinfo)));

	/* the loop was stolen from ps - but it backs through the data.
	 * We're going through forward, which means we need to play
	 * slightly different games.
	 */
#if 1
	kp += nentries - 1;
	for (i = 1; i <= nentries; i++, --kp) {
#else
	for (i = nentries; --i >= 0; ++kp) {
#endif
		struct extern_proc *p;
		struct eproc *e;
		KINFO *ki;
#ifdef TESTING
		char *ttname = NULL;
#endif /* def TESTING */

		time_value_t total_time, system_time, user_time;

		ki = &kinfo[i];	
		ki->ki_p = kp;

		get_task_info(ki);

		p = KI_PROC (ki);
		e = KI_EPROC (ki);
		state = p->p_stat == SZOMB ? SZOMB : ki->state;
		getproclline (ki, &command_name, &cmdlen, 0, 1);

		user_time = ki->tasks_info.user_time;
		time_value_add (&user_time, &ki->times.user_time);
		system_time = ki->tasks_info.system_time;
		time_value_add (&system_time, &ki->times.system_time);
		total_time = user_time;
		time_value_add (&total_time, &system_time);

#ifndef	TH_USAGE_SCALE
#define	TH_USAGE_SCALE	1000
#endif	TH_USAGE_SCALE
#define usage_to_percent(u)	((u*100)/TH_USAGE_SCALE)
#define usage_to_tenths(u)	(((u*1000)/TH_USAGE_SCALE) % 10)
		sprintf (pctcpu, "%d.%01d", usage_to_percent (ki->cpu_usage),
			usage_to_tenths (ki->cpu_usage));
		sprintf (pctmem, "%.1f", ((float) ki->tasks_info.resident_size)
			* 100 / mempages);

#ifdef TIME_IN_MICROSECONDS
#define time_value_to_ticks(x) 1000000LL * (x).seconds + (x).microseconds
#else /* def TIME_IN_MICROSECONDS */
#define time_value_to_ticks(x) ((long) (x).seconds * 100 + \
	((long) (x).microseconds + 5000) / 10000)
#endif /* def TIME_IN_MICROSECONDS */

#ifdef TESTING

		/* Since we're testing stand-alone, we'll provide a ttydev value
		 * for the look of the thing. We don't have to do this when live,
		 * since the .xs module does this automagically when it sees the
		 * ttynum value. This means ttydev must be passed BEFORE
		 * ttynum.
		 */

		if (e->e_tdev != NODEV) ttname = devname (e->e_tdev, S_IFCHR);
		if (ttname == NULL) ttname = "";

#ifdef DEBUGGING
		printf ("\nPid: %d; i:%d; kp: %p\n", p->p_pid, i, kp);
#else /* def DEBUGGING */
		printf ("\nPid: %d\n", p->p_pid);
#endif /* def DEBUGGING */
		printf ("    ppid: %d\n", e->e_ppid);
		printf ("    pgid: %d\n", e->e_pgid);
		printf ("     uid: %d\n", e->e_pcred.p_ruid);
		printf ("     gid: %d\n", e->e_pcred.p_rgid);
		printf ("    euid: %d\n", e->e_ucred.cr_uid);
		printf ("    egid: %d\n", e->e_ucred.cr_gid);
		printf ("    suid: %d\n", e->e_pcred.p_svuid);
		printf ("    sgid: %d\n", e->e_pcred.p_svgid);
		printf ("priority: %u\n", ki->curpri);
		printf ("    size: %lu Kb\n", (u_long)ki->tasks_info.virtual_size/1024);
		printf ("     rss: %lu Kb\n", (u_long)ki->tasks_info.resident_size/1024);
		printf ("   flags: %#0x\n", p->p_flag);
		printf ("    nice: %d\n", p->p_nice);
		printf (" session: %p\n", e->e_sess);
#ifdef TIME_IN_MICROSECONDS
		printf ("    time: %lld microseconds\n", time_value_to_ticks (total_time));
		printf ("   stime: %lld microseconds\n", time_value_to_ticks (system_time));
		printf ("   utime: %lld microseconds\n", time_value_to_ticks (user_time));
#else /* def TIME_IN_MICROSECONDS */
		printf ("    time: %ld centiseconds\n", time_value_to_ticks (total_time));
		printf ("   stime: %ld centiseconds\n", time_value_to_ticks (system_time));
		printf ("   utime: %ld centiseconds\n", time_value_to_ticks (user_time));
#endif /* def TIME_IN_MICROSECONDS */
		printf ("   start: %ld\n", (unsigned long) p->p_starttime.tv_sec);
		printf ("   wchan: %p\n", p->p_wchan);
		printf ("  ttydev: %s\n", ttname);
		printf ("  ttynum: %ld\n", (long) e->e_tdev);
		printf ("    %%cpu: %s\n", pctcpu);
		printf ("    %%mem: %s\n", pctmem);
		printf ("   state: %d (%s)\n", state, States[state]);
		printf ("     cmd: %s\n", cmdlen ? command_name : " (none available)");
		printf ("   fname: %s\n", p->p_comm);

#else /* def TESTING */

		/* Send if off to Perl */

		bless_into_proc (Format, Fields,
			p->p_pid,
			e->e_ppid,
			e->e_pgid,
			e->e_pcred.p_ruid,
			e->e_pcred.p_rgid,
			e->e_ucred.cr_uid,
			e->e_ucred.cr_gid,
			e->e_pcred.p_svuid,
			e->e_pcred.p_svgid,
			ki->curpri,
			ki->tasks_info.virtual_size/1024,
			ki->tasks_info.resident_size/1024,
			p->p_flag,
			p->p_nice,
			e->e_sess,
			time_value_to_ticks (total_time),
			time_value_to_ticks (system_time),
			time_value_to_ticks (user_time),
			p->p_starttime.tv_sec,
			p->p_wchan,
			"",		/* The .xs code provides ttydev automagically */
			e->e_tdev,
			pctcpu,
			pctmem,
			States[state],
			cmdlen ? command_name : "",
			p->p_comm
			);

#endif /* def TESTING */

		if (command_name != NULL) free (command_name);
	}

	FREE_BUFFERS

}

/* We're done with this definition of FREE_BUFFERS. Get rid of it in
 * the hope that using it below (as a consequence of an inappropriate
 * use of DIE_HORRIBLY, for example) will generate a more obvious
 * error message.
 */

#undef FREE_BUFFERS


/*
 * The interface used to fetch the command arguments changed
 * drastically between Jaguar (MacOS 10.2, Darwin 6.something) and
 * Panther (MacOS 10.3, Darwin version unknown to me at this time).
 * The corresponding module of the ps command changed likewise.
 * Unfortunately, we have to either keep both interfaces around,
 * or abandon Panther (stupid!), or abandon Jaguar (which I'm reluctant
 * to do, since I have a not-so-sneaking sympathy for those who upgrade
 * no software before its time). -- TRW
 */

#ifdef KERN_PROCARGS2

/*
 * The following is pretty much verbatim from module print.c of the
 * Panther version of the ps command. Specifically, it's from
 * adv_cmds-63. But the calling sequence has been modified for our
 * convenience. Specifically, the global variable eflg has been made
 * into an argument. -- TRW
 */

/*
 * Get command and arguments.
 *
 * If the global variable eflg is non-zero and the user has permission to view
 * the process's environment, the environment is included.
 */
static void
getproclline(KINFO *k, char **command_name, int *cmdlen, int eflg, int show_args)
{
	int		mib[3], argmax, nargs, c = 0;
	size_t		size;
	char		*procargs, *sp, *np, *cp;
/*  Made into a command argument. -- TRW
 *	extern int	eflg;
 */

	/* Get the maximum process arguments size. */
	mib[0] = CTL_KERN;
	mib[1] = KERN_ARGMAX;

	size = sizeof(argmax);
	if (sysctl(mib, 2, &argmax, &size, NULL, 0) == -1) {
		goto ERROR_A;
	}

	/* Allocate space for the arguments. */
	procargs = (char *)malloc(argmax);
	if (procargs == NULL) {
		goto ERROR_A;
	}

	/*
	 * Make a sysctl() call to get the raw argument space of the process.
	 * The layout is documented in start.s, which is part of the Csu
	 * project.  In summary, it looks like:
	 *
	 * /---------------\ 0x00000000
	 * :               :
	 * :               :
	 * |---------------|
	 * | argc          |
	 * |---------------|
	 * | arg[0]        |
	 * |---------------|
	 * :               :
	 * :               :
	 * |---------------|
	 * | arg[argc - 1] |
	 * |---------------|
	 * | 0             |
	 * |---------------|
	 * | env[0]        |
	 * |---------------|
	 * :               :
	 * :               :
	 * |---------------|
	 * | env[n]        |
	 * |---------------|
	 * | 0             |
	 * |---------------| <-- Beginning of data returned by sysctl() is here.
	 * | argc          |
	 * |---------------|
	 * | exec_path     |
	 * |:::::::::::::::|
	 * |               |
	 * | String area.  |
	 * |               |
	 * |---------------| <-- Top of stack.
	 * :               :
	 * :               :
	 * \---------------/ 0xffffffff
	 */
	mib[0] = CTL_KERN;
	mib[1] = KERN_PROCARGS2;
	mib[2] = KI_PROC(k)->p_pid;

	size = (size_t)argmax;
	if (sysctl(mib, 3, procargs, &size, NULL, 0) == -1) {
		goto ERROR_B;
	}

	memcpy(&nargs, procargs, sizeof(nargs));
	cp = procargs + sizeof(nargs);

	/* Skip the saved exec_path. */
	for (; cp < &procargs[size]; cp++) {
		if (*cp == '\0') {
			/* End of exec_path reached. */
			break;
		}
	}
	if (cp == &procargs[size]) {
		goto ERROR_B;
	}

	/* Skip trailing '\0' characters. */
	for (; cp < &procargs[size]; cp++) {
		if (*cp != '\0') {
			/* Beginning of first argument reached. */
			break;
		}
	}
	if (cp == &procargs[size]) {
		goto ERROR_B;
	}
	/* Save where the argv[0] string starts. */
	sp = cp;

	/*
	 * Iterate through the '\0'-terminated strings and convert '\0' to ' '
	 * until a string is found that has a '=' character in it (or there are
	 * no more strings in procargs).  There is no way to deterministically
	 * know where the command arguments end and the environment strings
	 * start, which is why the '=' character is searched for as a heuristic.
	 */
	for (np = NULL; c < nargs && cp < &procargs[size]; cp++) {
		if (*cp == '\0') {
			c++;
			if (np != NULL) {
				/* Convert previous '\0'. */
				*np = ' ';
			}
			/* Note location of current '\0'. */
			np = cp;

			if (!show_args) {
				/*
				 * Don't convert '\0' characters to ' '.
				 * However, we needed to know that the
				 * command name was terminated, which we
				 * now know.
				 */
				break;
			}
		}
	}

	/*
	 * If eflg is non-zero, continue converting '\0' characters to ' '
	 * characters until no more strings that look like environment settings
	 * follow.
	 */
	if ( (eflg != 0) && ( (getuid() == 0) || (KI_EPROC(k)->e_pcred.p_ruid == getuid()) ) ) {
		for (; cp < &procargs[size]; cp++) {
			if (*cp == '\0') {
				if (np != NULL) {
					if (&np[1] == cp) {
						/*
						 * Two '\0' characters in a row.
						 * This should normally only
						 * happen after all the strings
						 * have been seen, but in any
						 * case, stop parsing.
						 */
						break;
					}
					/* Convert previous '\0'. */
					*np = ' ';
				}
				/* Note location of current '\0'. */
				np = cp;
			}
		}
	}

	/*
	 * sp points to the beginning of the arguments/environment string, and
	 * np should point to the '\0' terminator for the string.
	 */
	if (np == NULL || np == sp) {
		/* Empty or unterminated string. */
		goto ERROR_B;
	}

	/* Make a copy of the string. */
	*cmdlen = asprintf(command_name, "%s", sp);

	/* Clean up. */
	free(procargs);
	return;

	ERROR_B:
	free(procargs);
	ERROR_A:
	*cmdlen = asprintf(command_name, "(%s)", KI_PROC(k)->p_comm);
}


#else	/* #ifdef KERN_PROCARGS2 */

/*
 * The following code is pretty much verbatim from module print.c of
 * the ps command. The version of print.c is unknown. It is identical
 * to, but certainly earlier than, the version in adv_cmds-46, which
 * is the most recent version that goes with Jaguar. The calling
 * sequence has been modified to pass eflg as an argument (rather than
 * a global variable), and to be the same as the Panther version. Also,
 * some band-aid code has been inserted late in the module to cover an
 * apparent bug. The Panther version of this subroutine was completely
 * rewritten, and uses a different sysctl funtion. -- TRW.
 */
static void getproclline (KINFO *k, char ** command_name, int *cmdlen,
	int eflg, int show_args)
/* show_args = 1, display environment; 0 = don't. */
{

	/*
	 *	Get command and arguments.
	 */
	int		command_length;
	char * cmdpath;
	volatile int 	*ip, *savedip;
	volatile char	*cp;
	int		nbad;
	char		c;
	char		*end_argc;
	int 		mib[4];
	char *		arguments;
	size_t		arguments_size = 4096;
	int			len=0;
	volatile unsigned int *valuep;
	unsigned int value;
	int blahlen=0, skiplen=0;

	/* A sysctl() is made to find out the full path that the command
	   was called with.
	*/

	*command_name = NULL;
	*cmdlen = 0;

	mib[0] = CTL_KERN;
	mib[1] = KERN_PROCARGS;
	mib[2] = KI_PROC(k)->p_pid;
	mib[3] = 0;

	arguments = (char *) malloc(arguments_size);
	if (sysctl(mib, 3, arguments, &arguments_size, NULL, 0) < 0) {
		goto retucomm;
	}
	end_argc = &arguments[arguments_size];

	ip = (int *)end_argc;
	ip -= 2;		/* last arg word and .long 0 */
	while (*--ip)
		if (ip == (int *)arguments)
		goto retucomm;

		savedip = ip;
		savedip++;
		cp = (char *)savedip;
	while (*--ip)
	    if (ip == (int *)arguments)
		goto retucomm;
        ip++;
        
        valuep = (unsigned int *)ip;
        value = *valuep;
        if ((value & 0xbfff0000) == 0xbfff0000) {
                ip++;ip++;
		valuep = ip;
               blahlen = strlen((char *) ip);
                skiplen = (blahlen +3 ) /4 ;
                valuep += skiplen;
                cp = (char *)valuep;
                while (!*cp) {
                    cp++;
                }
                savedip = (int *) cp;
        }

	nbad = 0;

	for (cp = (char *)savedip; cp < (end_argc-1); cp++) {
	    c = *cp & 0177;
	    if (c == 0)
		*cp = ' ';
	    else if (c < ' ' || c > 0176) {
		if (++nbad >= 5*(eflg+1)) {
		    *cp++ = ' ';
		    break;
		}
		*cp = '?';
	    }
		else if (eflg == 0 && c == '=') {
		while (*--cp != ' ')
			if (cp <= (char *)ip)
			break;
		break;
		}
	}

	*cp = 0;
#if 0
	while (*--cp == ' ')
	    *cp = 0;
#endif
	cp = (char *)savedip;
	command_length = end_argc - cp;	/* <= MAX_COMMAND_SIZE */

	if (cp[0] == '-' || cp[0] == '?' || cp[0] <= ' ') {
		/*
		 *	Not enough information - add short command name
		 */
		/*
		 *	I have a report of this section of the code failing under
		 *	Panther because command_length < 0. When Jaguar hits this
		 *  code it typically has a largeish value of command_length
		 *  (200 bytes plus). Both are clearly bogus, though in the
		 *  case of Panther, the problem is benign. The problem is
		 *  that Jaguar uses a completely different sysctl call to
		 *  get the data, so I can't simply use that code. The
		 *  print.c module of the ps command from adv_cmds_43 (the
		 *  latest Jaguar version) is identical to the one I based
		 *  this code on originally. So no help there. Until I can
		 *  figure something better, we'll just have to rely on this
		 *  band-aid. A possible way to proceed, once I upgrade to
		 *  Panther myself, is to conditionalize on the existence of
		 *  KERN_PROCARGS2 to decide which version of this subroutine
		 *  to use. Sigh. -- TRW
		 */
#ifdef DEBUGGING
		fprintf (stdout, "Debug - getproclline found short cmd; pid %d command_length %d\n",
			KI_PROC(k)->p_pid, command_length);
#endif
		if (command_length > 0) {
			len = ((unsigned)command_length + MAXCOMLEN + 5);
			cmdpath = (char *)malloc(len);
			(void) strncpy(cmdpath, (const char *) cp, command_length);
			(void) strcat(cmdpath, " (");
			(void) strncat(cmdpath, KI_PROC(k)->p_comm,
				MAXCOMLEN+1);
			(void) strcat(cmdpath, ")");
			*command_name = cmdpath;
			*cmdlen = len;
			}
		  else {
		  	cmdpath = (char *)malloc(2);
		  	strncpy (cmdpath, "", 2);
		  	*command_name = cmdpath;
			*cmdlen = 0;
			}
		free(arguments);
		return;
	}
	else {
		cmdpath = (char *)malloc((unsigned)command_length + 1);
		(void) strncpy(cmdpath, (const char *) cp, command_length);
		cmdpath[command_length] = '\0';
		*command_name = cmdpath;
		*cmdlen = command_length;
		free(arguments);
		return;
	}

retucomm:
	len = (MAXCOMLEN + 5);
	cmdpath = (char *)malloc(len);
	(void) strcpy(cmdpath, " (");
	(void) strncat(cmdpath, KI_PROC(k)->p_comm,
			MAXCOMLEN+1);
	(void) strcat(cmdpath, ")");
	*cmdlen = len;
	*command_name = cmdpath;
	free(arguments);
	return;
}

#endif	/* #ifdef KERN_PROCARGS2 */

static int mach_state_order (int s, long sleep_time);
static int thread_schedinfo (KINFO *ki, thread_port_t thread,
	policy_t pol, void * buf);

static int get_task_info (KINFO *ki) 
{
	kern_return_t   	error;
	unsigned int		info_count = TASK_BASIC_INFO_COUNT;
	unsigned int 		thread_info_count = THREAD_BASIC_INFO_COUNT;
	pid_t				pid;
	int j, err = 0;

	pid = KI_PROC(ki)->p_pid;
	if (task_for_pid(mach_task_self(), pid, &ki->task) != KERN_SUCCESS) {
		return(1);
	}
	info_count = TASK_BASIC_INFO_COUNT;
	error = task_info(ki->task, TASK_BASIC_INFO, (task_info_t) &ki->tasks_info, &info_count);
	if (error != KERN_SUCCESS) {
		ki->invalid_tinfo=1;
#ifdef DEBUG
		mach_error("Error calling task_info()", error);
#endif
		return(1);
	}
	{
		vm_region_basic_info_data_64_t	b_info;
		vm_address_t					address = GLOBAL_SHARED_TEXT_SEGMENT;
		vm_size_t					size;
		mach_port_t					object_name;

		/*
		 * try to determine if this task has the split libraries
		 * mapped in... if so, adjust its virtual size down by
		 * the 2 segments that are used for split libraries
		 */
		info_count = VM_REGION_BASIC_INFO_COUNT_64;
		error = vm_region_64(ki->task, &address, &size, VM_REGION_BASIC_INFO,
					(vm_region_info_t)&b_info, &info_count, &object_name);
		if (error == KERN_SUCCESS) {
			if (b_info.reserved && size == (SHARED_TEXT_REGION_SIZE) &&
				ki->tasks_info.virtual_size > (SHARED_TEXT_REGION_SIZE + SHARED_DATA_REGION_SIZE))
					ki->tasks_info.virtual_size -= (SHARED_TEXT_REGION_SIZE + SHARED_DATA_REGION_SIZE);
		}
	}
	info_count = TASK_THREAD_TIMES_INFO_COUNT;
	error = task_info(ki->task, TASK_THREAD_TIMES_INFO, (task_info_t) &ki->times, &info_count);
	if (error != KERN_SUCCESS) {
		ki->invalid_tinfo=1;
#ifdef DEBUG
		mach_error("Error calling task_info()", error);
#endif
		return(1);
	}

	switch(ki->tasks_info.policy) {

		case POLICY_TIMESHARE :
			info_count = POLICY_TIMESHARE_INFO_COUNT;
			error = task_info(ki->task, TASK_SCHED_TIMESHARE_INFO, (task_info_t) &ki->schedinfo.tshare, &info_count);
			if (error != KERN_SUCCESS) {
				ki->invalid_tinfo=1;
#ifdef DEBUG
				mach_error("Error calling task_info()", error);
#endif
				return(1);
			}

			ki->curpri = ki->schedinfo.tshare.cur_priority;
			ki->basepri = ki->schedinfo.tshare.base_priority;
			break;

		case POLICY_RR :
	 		info_count = POLICY_RR_INFO_COUNT;
			error = task_info(ki->task, TASK_SCHED_RR_INFO, (task_info_t) &ki->schedinfo.rr, &info_count);
			if (error != KERN_SUCCESS) {
				ki->invalid_tinfo=1;
#ifdef DEBUG
				mach_error("Error calling task_info()", error);
#endif
				return(1);
			}

			ki->curpri = ki->schedinfo.rr.base_priority;
			ki->basepri = ki->schedinfo.rr.base_priority;
			break;

		case POLICY_FIFO :
			info_count = POLICY_FIFO_INFO_COUNT;
			error = task_info(ki->task, TASK_SCHED_FIFO_INFO, (task_info_t) &ki->schedinfo.fifo, &info_count);
			if (error != KERN_SUCCESS) {
				ki->invalid_tinfo=1;
#ifdef DEBUG
				mach_error("Error calling task_info()", error);
#endif
				return(1);
			}

			ki->curpri = ki->schedinfo.fifo.base_priority;
			ki->basepri = ki->schedinfo.fifo.base_priority;
			break;
	}

	ki->invalid_tinfo=0;

	ki->cpu_usage=0;
	error = task_threads(ki->task, &ki->thread_list, &ki->thread_count);
	if (error != KERN_SUCCESS) {
		mach_port_deallocate(mach_task_self(),ki->task);
#ifdef DEBUG
		mach_error("Call to task_threads() failed", error);
#endif
		return(1);
	}
	err=0;
	ki->state = STATE_MAX;
	//ki->curpri = 255;
	//ki->basepri = 255;
	ki->swapped = 1;
	ki->thval = malloc(ki->thread_count * sizeof(struct thread_values));
	if (ki->thval != NULL) {
		for (j = 0; j < ki->thread_count; j++) {
			int tstate;
			thread_info_count = THREAD_BASIC_INFO_COUNT;
			error = thread_info(ki->thread_list[j], THREAD_BASIC_INFO,
				(thread_info_t)&ki->thval[j].tb,
				&thread_info_count);
			if (error != KERN_SUCCESS) {
#ifdef DEBUG
				mach_error("Call to thread_info() failed", error);
#endif
				err=1;
			}
			error = thread_schedinfo(ki, ki->thread_list[j],
				ki->thval[j].tb.policy, &ki->thval[j].schedinfo);
			if (error != KERN_SUCCESS) {
#ifdef DEBUG
				mach_error("Call to thread_info() failed", error);
#endif
				err=1;
			}
			ki->cpu_usage += ki->thval[j].tb.cpu_usage;
			tstate = mach_state_order(ki->thval[j].tb.run_state,
					ki->thval[j].tb.sleep_time);
			if (tstate < ki->state)
				ki->state = tstate;
			if ((ki->thval[j].tb.flags & TH_FLAGS_SWAPPED ) == 0)
				ki->swapped = 0;
			mach_port_deallocate(mach_task_self(),
				ki->thread_list[j]);
		}
		free (ki->thval);
		ki->thval = NULL;
	}
	ki->invalid_thinfo = err;
	/* Deallocate the list of threads. */
	error = vm_deallocate(mach_task_self(), 
		(vm_address_t)(ki->thread_list),
		 sizeof(thread_port_array_t) * ki->thread_count);
	if (error != KERN_SUCCESS) {
#ifdef DEBUG
		mach_error("Trouble freeing thread_list", error);
#endif
	}

	mach_port_deallocate(mach_task_self(),ki->task);
	return(0);
}

static int mach_state_order (int s, long sleep_time)
{      
	switch (s) {
	case TH_STATE_RUNNING:		return(1);
	case TH_STATE_UNINTERRUPTIBLE:
								return(2);
	case TH_STATE_WAITING:		return((sleep_time > 20) ? 4 : 3);
	case TH_STATE_STOPPED:		return(5);
	case TH_STATE_HALTED:		return(6);  
	default:					return(7); 
	}
}

static int thread_schedinfo (KINFO *ki, thread_port_t thread,
	policy_t pol, void * buf)
{
	unsigned int		count;
	int ret = KERN_FAILURE;

	switch (pol) {

	case POLICY_TIMESHARE:
		count = POLICY_TIMESHARE_INFO_COUNT;
		ret = thread_info(thread, THREAD_SCHED_TIMESHARE_INFO,
					(thread_info_t)buf, &count);
		if((ret == KERN_SUCCESS) && (ki->curpri < (((struct policy_timeshare_info *)buf)->cur_priority)))
			ki->curpri  = ((struct policy_timeshare_info *)buf)->cur_priority;
		break;

	case POLICY_FIFO:
		count = POLICY_FIFO_INFO_COUNT;
		ret = thread_info(thread, THREAD_SCHED_FIFO_INFO,
					buf, &count);
		if((ret == KERN_SUCCESS) && (ki->curpri < (((struct policy_fifo_info *)buf)->base_priority)))
			ki->curpri  = ((struct policy_fifo_info *)buf)->base_priority;
		break;

	case POLICY_RR:
		count = POLICY_RR_INFO_COUNT;
		ret = thread_info(thread, THREAD_SCHED_RR_INFO,
					buf, &count);
		if((ret == KERN_SUCCESS) && (ki->curpri < (((struct policy_rr_info *)buf)->base_priority)))
			ki->curpri  = ((struct policy_rr_info *)buf)->base_priority;
		break;
	}
	return(ret);
}