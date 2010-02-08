#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "os/DecOSF.h"

/* Make sure /proc is mounted */
char* OS_initialize()
{
  struct statvfs svfs;

  static char* no_proc = "/proc unavailable";
  if( statvfs("/proc", &svfs) == -1 )
  {
    return no_proc;
  }

  return NULL;
}


/* FIXME we should get minimum info like process ID and ownership from
   file stat-- does this work for IOCTL-proc? Does it for FS-proc? It
   does on linux... */

void OS_get_table()
{
  int maxproc, procnum, procbase;
  int psdata;
  char pathbuf[MAXPATHLEN];

  /* defined in <sys/table.h> */
  struct tbl_procinfo pstbl[PROCCNT];

  /* defined in <sys/procfs.h> */
  struct prpsinfo psbuf;
  struct prcred pscredbuf;

  /* variables to hold some values for bless_into_proc */
  char state[20]; 
  char pctcpu[7];
  char pctmem[7];
  long pr_age, now;
  uid_t *psgroups;
  int i;
  AV *group_array;
  SV *group_ref;

  if( (maxproc = table( TBL_PROCINFO, 0, NULL, (unsigned int) -1, 0 )) == -1 )
    return;

  /* loop over all processes */
  procbase = -1;
  for( procnum = 0; procnum < maxproc; procnum ++)
  {
    /* Get the table entries */
    if( (procbase < 0) || (procnum >= (procbase + PROCCNT)) )
    {
      procbase = procnum;
      if( table( TBL_PROCINFO, procbase, pstbl, PROCCNT,
	  sizeof (struct tbl_procinfo)) == -1 )
	continue;
    }

    /* See if the process exists */
    if( pstbl[procnum % PROCCNT].pi_status == PI_EMPTY )
      continue;

    /* Construct path of the form /proc/proc_number */
    sprintf( pathbuf, "/proc/%d", pstbl[procnum % PROCCNT].pi_pid );

    if( (psdata = open( pathbuf, O_RDONLY )) == -1 ) continue;

    if( ioctl(psdata, PIOCPSINFO, &psbuf) == -1 ) continue; 
    if( ioctl(psdata, PIOCCRED, &pscredbuf) == -1 ) continue; 
    if( (psgroups = malloc(pscredbuf.pr_ngroups * sizeof (u_int))) == NULL )
      continue;
    if( ioctl(psdata, PIOCGROUPS, psgroups) == -1 )
    {
      free(psgroups);
      continue;
    }

    close(psdata);

    /* translate process state, makros defined in <sys/proc.h> */
    switch( psbuf.pr_state)
    {
      case SSLEEP: 
        strcpy(state, SLEEP);
        break;
      case SWAIT:
        strcpy(state, WAIT);
        break;
      case SRUN:
        strcpy(state, RUN);
        break;
      case SZOMB:
        strcpy(state, ZOMBIE);
        break;
      case SSTOP:
        strcpy(state, STOP);
        break;
      case SIDL:
        strcpy(state, IDLE);
        break;
      }

    now = time(NULL);
    pr_age = now - psbuf.pr_start.tv_sec;	
    if (pr_age < 1) {
	sprintf( pctcpu, "%3.2f", (float) 0.0); /* pr_time.tv_sec is 0 */
    } else {
       sprintf( pctcpu, "%3.2f", (float) ((psbuf.pr_time.tv_sec * 100) / pr_age));
    }

    /* create groups array */
    group_array = newAV();
    for (i = 0; i < pscredbuf.pr_ngroups; i++)
    {
      av_push(group_array, newSViv(psgroups[i]));
    }
    free(psgroups);
    group_ref = newRV_noinc((SV *) group_array);

    bless_into_proc( Format,           
                     Fields,

                     psbuf.pr_uid,           /* uid, uid_t is int */
                     psbuf.pr_gid,           /* gid, gid_t is int */
		     pscredbuf.pr_euid,      /* euid, uid_t is int */
		     pscredbuf.pr_egid,      /* egid, uid_t is int */
		     pscredbuf.pr_suid,      /* suid, uid_t is int */
		     pscredbuf.pr_sgid,      /* sgid, uid_t is int */
		     group_ref,              /* groups, perl ref to array */
                     psbuf.pr_pid,           /* pid, pid_t is int */
                     psbuf.pr_ppid,          /* ppid, pid_t is int */
                     psbuf.pr_pgrp,          /* pgrp, pid_t is int */ 
                     psbuf.pr_sid,           /* sess, pid_t is int */
                     psbuf.pr_pri,           /* priority, long, */   
                     psbuf.pr_ttydev,        /* ttynum, dev_t is int */
                     psbuf.pr_flag,          /* flags, u_long */
                     psbuf.pr_time.tv_sec,   /* time, long */
                     psbuf.pr_size * getpagesize(),   /* size (bytes) */
                     psbuf.pr_rssize * getpagesize(), /* rss (bytes)  */
                     psbuf.pr_start.tv_sec,  /* start */
                     psbuf.pr_fname,         /* fname */
                     pctcpu,                 /* pctcpu */
                     state,                  /* state */
                     psbuf.pr_psargs         /* cmndline */ 
                   );

  }
}
