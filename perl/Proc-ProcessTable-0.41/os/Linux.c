
#include "os/Linux.h"

unsigned long Hertz;
/* NOTE: Before this was actually milliseconds even though it said microseconds, now it is correct. */
#define JIFFIES_TO_MICROSECONDS(x) (((x)*1e6)/Hertz)
static int init_Hertz_value(void);

/* Given a path to a /proc/XXX/stat file and a pointer to a procstat
   struct, fill the struct */
struct procstat* get_procstat( char* path, struct procstat* prs)
{
	FILE* fp;
	int result;

	/* hmmm...  no file pointer???  Then Return NULL */
	if( (fp = fopen( path, "r" )) == NULL )
		return NULL;

	result = fscanf(fp, 
			"%d %s %c %d %d %d %d %d %u %u %u %u %u %Ld %Ld %Ld %Ld %d %d %u %u %lu %u %u %u %u %u %u %u %u %d %d %d %d %u",
			&prs->pid, 
			prs->comm,		/* char comm[FILENAME_MAX]; */
			&prs->state, 
			&prs->ppid, 
			&prs->pgrp, 
			&prs->session, 
			&prs->tty, 
			&prs->tpgid,
			&prs->flags,
			&prs->minflt,
			&prs->cminflt,
			&prs->majflt, 
			&prs->cmajflt, 
			&prs->utime,
			&prs->stime,
			&prs->cutime,
			&prs->cstime,
			&prs->counter,
			&prs->priority,
			&prs->timeout,
			&prs->itrealvalue,
			&prs->starttime,
			&prs->vsize,
			&prs->rss,
			&prs->rlim,
			&prs->startcode,
			&prs->endcode,
			&prs->startstack,
			&prs->kstkesp,
			&prs->kstkeip,
			&prs->signal,
			&prs->blocked,
			&prs->sigignore,
			&prs->sigcatch,
			&prs->wchan);

	fclose(fp);

	/* 35 items in scanf's list... It's all or nothing baby */
	if (result != 35)
		return NULL;

	prs->utime = JIFFIES_TO_MICROSECONDS(prs->utime); 
	prs->stime = JIFFIES_TO_MICROSECONDS(prs->stime); 
	prs->cutime = JIFFIES_TO_MICROSECONDS(prs->cutime); 
	prs->cstime = JIFFIES_TO_MICROSECONDS(prs->cstime); 
	prs->starttime /= Hertz;
	prs->timeout = JIFFIES_TO_MICROSECONDS(prs->timeout); 

	return prs;
}

/* Make sure /proc is mounted and initialize some global vars */
char* OS_initialize()
{
	char cbuf[1024];
	FILE* fp;
	struct statfs sfs;

	static char* no_proc = "/proc unavailable";
	if( statfs("/proc", &sfs) == -1 )
		return no_proc;

	/* Get boottime from /proc/stat */
	Btime = 0;
	if( (fp = fopen( "/proc/stat", "r" )) != NULL ){ 
		while(!feof(fp)) {
			if(fscanf(fp,"btime %ld", &Btime) == 1)
				break;
			if (fgets(cbuf, 1024, fp) == NULL)
				break;
				
		}
		fclose(fp);
	}

	/* Get total system memory from /proc/meminfo and convert to pages */
	Sysmem = 0;
	if( (fp = fopen( "/proc/meminfo", "r" )) != NULL ) { 
		while(!feof(fp)) {
			if(fscanf(fp,"MemTotal: %u", &Sysmem) == 1){
			        Sysmem *= 1024; /* Convert Sysmem from kB to bytes */ 
				Sysmem /= getpagesize(); /* Convert Sysmem from bytes to pages */
				break;
			}
			if (fgets(cbuf, 1024, fp) == NULL)
				break;
		}
		fclose(fp);
	}
	
	init_Hertz_value();
	return NULL;
}



void OS_get_table()
{
	DIR *procdir;
	struct dirent *procdirp;
	char pathbuf[PATH_MAX];
	struct stat filestat;
	FILE* fp;

	/* for bless_into_proc, we make 2 to guarantee that strings terminate */
	struct procstat prs[1];
	char fname[NAME_MAX];
	unsigned long start;
	char state[32];
	char cmndline[ARG_MAX]; 
	char pctmem[32];
	char pctcpu[32];
	char cbuf[1024];
	static char format[F_LASTFIELD + 2];

	size_t pagesize = getpagesize();
  
	/* used by exec, added by scip */
	int size;
	char exec[ARG_MAX];

	/* used by cwd, added by scip */
	char curdir[ARG_MAX];

	/* used by the euid/egid stuff, added by scip */
	int dummyid, euid, suid, fuid, egid, sgid, fgid, i;

	if( (procdir = opendir("/proc")) == NULL )
		return;

	/* Iterate through all the process entries (numeric) under /proc */
	while( (procdirp = readdir(procdir)) != NULL ) {

		/* Only look at this file if it's a proc id; that is, all numbers */
		if( strtok(procdirp->d_name, "0123456789") != NULL )
			continue;

		/* zero our format */
		strncpy(format, Defaultformat, sizeof(format));

		/* get puid and pgid from proc file */
		sprintf(pathbuf, "%s%s", "/proc/", procdirp->d_name);
		if( stat( pathbuf, &filestat ) != -1 ) {
			format[F_UID] = tolower(format[F_UID]); /* uid */
			format[F_GID] = tolower(format[F_GID]); /* gid */
		}

		/* We'll initialize ALL strings here...  why? 
			cause if any of the file read code should fail, 
			then the perl bless call will NOT try to read garbage.
			In theory, it doesn't matter.  However just in case...
			 */
		pctcpu[0] = pctmem[0] = state[0] = '\0';
		cmndline[0] = fname[0] = exec[0] = '\0';
		curdir[0] = '\0';

		/* get stuff out of /proc/PROC_ID/stat */
		memset(prs, 0, sizeof(struct procstat));
		strcat( pathbuf, "/stat" );

		/* I know that we we want to collect everything we can... 
			BUT, if we don't get a procstat struct, everything is pretty useless 
			I mean, we don't even get a pid.  So how can this be useful?
			So... we skip this one if no prs exists */
		if( get_procstat(pathbuf, prs) == NULL )
			continue;

		/* make all these fields valid */
		for( i = F_PID; i <= F_WCHAN; i++ )
			format[i] = tolower(format[i]);

		/* Get rid of the parens. There's probably a better way... */
		strcpy(fname, strtok(prs->comm, "()"));
		format[F_FNAME] = tolower(format[F_FNAME]); /* fname */

		/* starttime is seconds since boot; convert to unix time */
		if( Btime != 0 ) {
			start = prs->starttime + Btime;
			format[F_START] = tolower(format[F_START]); /* start */
		}

		/* calculate pctcpu - NOTE: This assumes the cpu time is in microsecond units! */
		sprintf( pctcpu, "%3.2f", (float) 100 * ((prs->utime + prs->stime)/1e6) / (time(NULL) - start) );
		format[F_PCTCPU] =  tolower(format[F_PCTCPU]); /* pctcpu */

		/* convert character state into one of our "official" readable
			states */
		switch(prs->state)
		{
			case 'R':
				strcpy(state, RUN);
				format[F_STATE] = tolower(format[F_STATE]); /* state */
				break;
			case 'S':
				strcpy(state, SLEEP);
				format[F_STATE] = tolower(format[F_STATE]); /* state */
				break;
			case 'D':
				strcpy(state, UWAIT);
				format[F_STATE] = tolower(format[F_STATE]); /* state */
				break;
			case 'Z':
				strcpy(state, ZOMBIE);
				format[F_STATE] = tolower(format[F_STATE]); /* state */
				break;
			case 'T':
				strcpy(state, STOP);
				format[F_STATE] = tolower(format[F_STATE]); /* state */
				break;
		}

		/* calculate pctmem */
		if( Sysmem != 0 ){  
			sprintf( pctmem, "%3.2f", (float) (prs->rss * 100 / Sysmem));
			format[F_PCTMEM] = tolower(format[F_PCTMEM]); /* pctmem */ 
		}

		/* 
		 * get the executable info, added by scip
		 * we do _not_ check if the readlink call failed, because
		 * we could simply not have permissions to read the link
		 * which then simply results in an empty "exec" field.
		 *
		 * Should have been ARG_MAX not PATH_MAX
		 *
		 * Also... we check the return value of readlink.
		 */
		sprintf(pathbuf, "%s%s%s", "/proc/", procdirp->d_name, "/exe");
		if ((size = readlink(pathbuf, exec, ARG_MAX - 1)) >= 0) {
			exec[size] = '\0';
			format[F_EXEC] = tolower(format[F_EXEC]);
		}

		/*
		 * get the euid, egid and so on out of /proc/$$/status
		 * where the 2 lines in which we are interested in are:
		 * [5] Uid:    500     500     500     500
		 * [6] Gid:    500     500     500     500
		 * added by scip
		 *
		 * modified to use a simpler parsing routine.  since fscanf, and fgets are
		 * probably implemented in assembly (or parts of it) it's safe to assume
		 * pretty good performance.
		 */
		sprintf(pathbuf, "%s%s%s", "/proc/", procdirp->d_name, "/status");
		if ( (fp = fopen( pathbuf, "r" )) != NULL) {
			i = 0;
			while(!feof(fp)) {
				if(fscanf(fp, "Uid: %d %d %d %d", &dummyid, &euid, &suid, &fuid) == 4) {
					format[F_EUID] = tolower(format[F_EUID]);
					format[F_SUID] = tolower(format[F_SUID]);
					format[F_FUID] = tolower(format[F_FUID]);
					i++;
					continue;
				}
				if(fscanf(fp, "Gid: %d %d %d %d", &dummyid, &egid, &sgid, &fgid) == 4) {
					format[F_EGID] = tolower(format[F_EGID]);
					format[F_SGID] = tolower(format[F_SGID]);
					format[F_FGID] = tolower(format[F_FGID]);
					i++;
					continue;
				}

				/* we'll short curcuit this if we have found what we wanted */
				if (i >= 2)
					break;

				/* this removes the whole line if it doesn't match at all */		
				if (fgets(cbuf, 1024, fp) == NULL)
					break;
			}
			fclose(fp);
		}

		/*
		 * get the cwd info, added by scip
		 */
		sprintf(pathbuf, "%s%s%s", "/proc/", procdirp->d_name, "/cwd");
		if ((size = readlink(pathbuf, curdir, ARG_MAX - 1)) >= 0) {
			curdir[size] = '\0';
			format[F_CWD] = tolower(format[F_CWD]);
		}

		/* get stuff out of /proc/PROC_ID/cmdline */
		sprintf(pathbuf, "%s%s%s", "/proc/", procdirp->d_name, "/cmdline");
		if( (fp = fopen( pathbuf, "r" )) != NULL ) { 
			size_t got;
			if( (got = fread(cmndline, sizeof(char), ARG_MAX, fp)) < 1 ) {
				strncpy(cmndline, fname, ARG_MAX);
				cmndline[ARG_MAX - 1] = '\0'; 
			} else {
				size_t i;
				for(i = 0; i < got; i++)
					if( cmndline[i] == '\0' )
						cmndline[i] = ' ';

				cmndline[got] = '\0'; 
			}
			format[F_CMNDLINE] = tolower(format[F_CMNDLINE]);
			fclose(fp);
		}

		/* Make sure our format is not all x's, which happens
			when a process is killed before we can read its info 
			and we get a blank object 

			-- unecessary because we don't return an object unless it has 
				at least a pid --
			*/
		/*if( strpbrk(format,"sil" ) == NULL)
			continue; */

		/* Go ahead and bless into a perl object */
		bless_into_proc( format,           
			Fields,
			filestat.st_uid,
			filestat.st_gid,
		 
			prs->pid,
			prs->ppid,
			prs->pgrp,
			prs->session,
			prs->priority,
			prs->tty,
			prs->flags,
			prs->minflt,
			prs->cminflt,
			prs->majflt,
			prs->cmajflt,
			prs->utime,
			prs->stime,
			prs->cutime,
			prs->cstime,
			prs->utime + prs->stime, /* FIXME check units w/solaris for consistency */
			prs->cutime + prs->cstime,
			prs->vsize,
			prs->rss * pagesize, 
			prs->wchan,
	 
			fname,
			start,
			pctcpu,
			state,
			pctmem,
			cmndline,
			exec,
			euid, suid, fuid,
			egid, sgid, fgid,
			curdir
		);
	}
	closedir(procdir);
}



  /* Get Hertz to convert jiffies to seconds */
  /* LIFTED FROM procps 2.0.2 */
#define BAD_OPEN_MESSAGE                                        \
"Error: /proc must be mounted\n"                                \
"  To mount /proc at boot you need an /etc/fstab line like:\n"  \
"      /proc   /proc   proc    defaults\n"                      \
"  In the meantime, mount /proc /proc -t proc\n"

#define STAT_FILE    "/proc/stat"
static int stat_fd = -1;
#define UPTIME_FILE  "/proc/uptime"
static int uptime_fd = -1;
#define LOADAVG_FILE "/proc/loadavg"
static int loadavg_fd = -1;
#define MEMINFO_FILE "/proc/meminfo"
static int meminfo_fd = -1;

static char buf[1024];

/* This macro opens filename only if necessary and seeks to 0 so
 * that successive calls to the functions are more efficient.
 * It also reads the current contents of the file into the global buf.
 */
#define FILE_TO_BUF(filename, fd) do{                           \
    static int n;                                               \
    if (fd == -1 && (fd = open(filename, O_RDONLY)) == -1) {    \
        ppt_warn(BAD_OPEN_MESSAGE);                   		\
        close(fd);                                              \
        return 0;                                               \
    }                                                           \
    lseek(fd, 0L, SEEK_SET);                                    \
    if ((n = read(fd, buf, sizeof buf - 1)) < 0) {              \
        perror(filename);                                       \
        close(fd);                                              \
        fd = -1;                                                \
        return 0;                                               \
    }                                                           \
    buf[n] = '\0';                                              \
}while(0)
                                                                                       
/***********************************************************************
 * Some values in /proc are expressed in units of 1/HZ seconds, where HZ
 * is the kernel clock tick rate. One of these units is called a jiffy.
 * The HZ value used in the kernel may vary according to hacker desire.
 * According to Linus Torvalds, this is not true. He considers the values
 * in /proc as being in architecture-dependant units that have no relation
 * to the kernel clock tick rate. Examination of the kernel source code
 * reveals that opinion as wishful thinking.
 *
 * In any case, we need the HZ constant as used in /proc. (the real HZ value
 * may differ, but we don't care) There are several ways we could get HZ:
 *
 * 1. Include the kernel header file. If it changes, recompile this library.
 * 2. Use the sysconf() function. When HZ changes, recompile the C library!
 * 3. Ask the kernel. This is obviously correct...
 *
 * Linus Torvalds won't let us ask the kernel, because he thinks we should
 * not know the HZ value. Oh well, we don't have to listen to him.
 * Someone smuggled out the HZ value. :-)
 *
 * This code should work fine, even if Linus fixes the kernel to match his
 * stated behavior. The code only fails in case of a partial conversion.
 *
 */

static int init_Hertz_value(void) 
{
#ifdef HZ
	Hertz = (unsigned long)HZ;    /* <asm/param.h> */
#else

	unsigned long user_j, nice_j, sys_j, other_j;  /* jiffies (clock ticks) */
	double up_1, up_2, seconds;
	unsigned long jiffies, h;

	do{
		FILE_TO_BUF(UPTIME_FILE,uptime_fd);  sscanf(buf, "%lf", &up_1);
		/* uptime(&up_1, NULL); */
		FILE_TO_BUF(STAT_FILE,stat_fd);
		sscanf(buf, "cpu %lu %lu %lu %lu", &user_j, &nice_j, &sys_j, &other_j);
		FILE_TO_BUF(UPTIME_FILE,uptime_fd);  sscanf(buf, "%lf", &up_2);
		/* uptime(&up_2, NULL); */
	} while((long)( (up_2-up_1)*1000.0/up_1 )); /* want under 0.1% error */

	jiffies = user_j + nice_j + sys_j + other_j;
	seconds = (up_1 + up_2) / 2;
	h = (unsigned long)( (double)jiffies/seconds );
	switch(h)
	{
		case   48 ...   52 :  Hertz =   50; break;
		case   58 ...   62 :  Hertz =   60; break;
		case   95 ...  105 :  Hertz =  100; break; /* normal Linux */
		case  124 ...  132 :  Hertz =  128; break;
		case  195 ...  204 :  Hertz =  200; break; /* normal << 1 */
		case  253 ...  260 :  Hertz =  256; break;
		case  393 ...  408 :  Hertz =  400; break; /* normal << 2 */
		case  790 ...  808 :  Hertz =  800; break; /* normal << 3 */
		case  990 ... 1010 :  Hertz = 1000; break;
		case 1015 ... 1035 :  Hertz = 1024; break; /* Alpha */
		default:
		Hertz = (sizeof(long)==sizeof(int)) ? 100UL : 1024UL;
		ppt_warn("Unknown HZ value! (%ld) Assume %ld.\n", h, Hertz);
	}
#endif /* HZ */
	return 0; /* useless, but FILE_TO_BUF has a return in it */
}

#include <fcntl.h>


