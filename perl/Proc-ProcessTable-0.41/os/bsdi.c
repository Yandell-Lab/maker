/**************************************************************/
/* Copyright (c) 1999, Magic Software Development, Inc.       */ 
/* Author:  Sean Ray Eskins. <sean@gilasoft.com>              */
/* This file is free software; it can be modified and/or      */
/* redistributed under the same terms as Perl itself.         */
/**************************************************************/

#include "bsdi.h"

extern void bless_into_proc(char* format, char** fields, ...);

char* OS_initialize() {
   return(NULL);
}

void OS_get_table() {
  kvm_t *kd;
  char errbuf[_POSIX2_LINE_MAX];
  struct kinfo_proc *procs;          /* array of processes */
  int count;                         /* returns number of processes */
  int i, j;
  int seconds, minutes, secleft;     /* for time[20] */
  int milsec, shortmsec;             /* for miliseconds of time[20] */      
  int ttynum;
  long start;
  char *ttydev;
  static time_t now;                 /* for started[20] */
  struct tm *tp;                     /* for month/day/hour/min/AM/PM fields of started[20] */
  int SECSPERHOUR = 3600;
  int SECSPERDAY = 24 * 3600;
  pid_t sesid;
  int length;
  char cmndline[MAXARGLN+1];
  char ** argv;

  /* for bless_into_proc */
  static char format[F_LASTFIELD + 1];
  
  /* variables to hold some values for bless_into_proc */
  char state[20];
  char cputime[20];
  char started[20];
  char session[20];
  char shortsess[20];

  /* Open the kvm interface, get a descriptor */
  if ((kd = kvm_open(NULL, NULL, NULL, 0, errbuf)) == NULL) {
    /* fprintf(stderr, "kvm_open: %s\n", errbuf); */
    ppt_croak("kvm_open: ", errbuf);
  }  
 
  /* Get the list of processes. */
  if ((procs = kvm_getprocs(kd, KERN_PROC_ALL, 0, &count)) == NULL) {
	kvm_close(kd);
    /* fprintf(stderr, "kvm_getprocs: %s\n", kvm_geterr(kd)); */
    ppt_croak("kvm_getprocs: ", kvm_geterr(kd));
  }

  /* Iterate through the processes in kinfo_proc, sending proc info */
  /* to bless_into_proc for each proc */
  for (i=0; i < count; i++) {
    static struct pstats ps;
    static struct session seslead;
    strcpy(format, Defaultformat);
	
    /* get ttydev */
    ttynum=procs[i].kp_eproc.e_tdev;
    ttydev=devname(ttynum, S_IFCHR);
    if (ttydev == NULL) ttydev = "??";

    /* get the state of processes */
    switch (procs[i].kp_proc.p_stat) 
      {
      case SIDL:
        strcpy(state, IDLE);
        break;
      case SRUN:
        strcpy(state, RUN);
        break;
      case SSLEEP:
        strcpy(state, SLEEP);
        break;
      case SSTOP:
        strcpy(state, STOP);
        break;
      case SZOMB: 
        strcpy(state, ZOMBIE);
        break;
      default:
        strcpy(state, UNKNOWN);
        break;
      }

    /* get the cpu time of processes */
    seconds=procs[i].kp_proc.p_rtime.tv_sec;
    milsec=procs[i].kp_proc.p_rtime.tv_usec;
    shortmsec=roundit(milsec);  

    if (seconds < 60) {
      if (seconds < 10)
        sprintf(cputime, "0:0%d.%d", seconds, shortmsec);
      else
        sprintf(cputime, "0:%d.%d", seconds, shortmsec);
    }
    else {
      minutes=seconds/60;
      secleft=seconds-(minutes * 60);
      if (secleft < 10)
        sprintf(cputime, "%d:0%d.%d", minutes, secleft, shortmsec);
      else
        sprintf(cputime, "%d:%d.%d", minutes, secleft, shortmsec);
    }

    /* get the start time of process (when started) */
    /* fill the pstats struct using kvm_read */
    if (kvm_read(kd, (u_long)procs[i].kp_proc.p_stats, (char *)&ps, sizeof(ps)) == sizeof(ps)) {
      start=ps.p_start.tv_sec;
      tp=localtime(&start);
      if (!now)
        (void)time(&now);
      if (now - ps.p_start.tv_sec < 24 * SECSPERHOUR) {
        static char fmt[] = __CONCAT("%l:%", "M%p");
        (void)strftime(started, sizeof(started) - 1, fmt, tp);
      } 
      else if (now - ps.p_start.tv_sec < 7 * SECSPERDAY) {
        static char fmt[] = __CONCAT("%a%", "I%p");
        (void)strftime(started, sizeof(started) - 1, fmt, tp);
      } 
      else
        (void)strftime(started, sizeof(started) - 1, "%e%b%y", tp);
    }

    /* get the session ID (ie: session pointer ID) */
    sprintf(session, "%x", (u_long)procs[i].kp_eproc.e_sess);   
    length=strlen(session);
    for (j=0; j < length; j++) {
      if (session[1] == '1') 
        shortsess[j]=session[j+1];
      else
        shortsess[j]=session[j+2];
    } 
    
    /* fill the session leader and proc struct using kvm_read */
    if (kvm_read(kd, (u_long)procs[i].kp_eproc.e_sess, (char *)&seslead, sizeof(seslead)) == sizeof(seslead)) {
      static struct proc leader; 
      if (kvm_read(kd, (u_long)seslead.s_leader, (char *)&leader, sizeof(leader)) == sizeof(leader)) {
        sesid=leader.p_pid;     
      }
    }

	/* retrieve the arguments */
	cmndline[0] = '\0';
	argv = kvm_getargv(kd, (const struct kinfo_proc *) &(procs[i]) , 0);
	if (argv) {
	  int j = 0;
	  while (argv[j] && strlen(cmndline) <= MAXARGLN) {
		strcat(cmndline, argv[j]);
		strcat(cmndline, " ");
		j++;
	  }
	}

    /* access everything else directly from the kernel, send it */
    /* into bless_into_proc */
    bless_into_proc( format,
                     Fields,
                     procs[i].kp_eproc.e_pcred.p_ruid,
                     procs[i].kp_eproc.e_pcred.p_rgid,
                     procs[i].kp_proc.p_pid,
                     procs[i].kp_eproc.e_ppid,
                     procs[i].kp_eproc.e_pgid,
                     procs[i].kp_proc.p_priority - PZERO,
                     shortsess,
                     sesid,
                     cputime,
                     procs[i].kp_eproc.e_wmesg,
                     procs[i].kp_proc.p_comm,
                     state,
                     started,
                     ttydev,
                     ttynum,
					 cmndline
                     );
  }
  if (kd) {
	kvm_close(kd);
  }
}
    




