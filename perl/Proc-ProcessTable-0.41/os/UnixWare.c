
#include "os/UnixWare.h"

/* Make sure /proc is mounted */
char* OS_initialize(){
  struct statvfs svfs;

  static char* no_proc = "/proc unavailable";
  if( statvfs("/proc", &svfs) == -1 ){
    return no_proc;
  }

  return NULL;
}


/* FIXME we should get minimum info like process ID and ownership from
   file stat-- does this work for IOCTL-proc? Does it for FS-proc? It
   does on linux... */

void OS_get_table(){
  DIR *procdir;
  struct dirent *procdirp;
  int psdata;
  char pathbuf[MAXPATHLEN];
  struct psinfo psbuf;
  lwpstat_t pr_state;
  
  /* variables to hold some values for bless_into_proc */
  char state[20]; 
  
  if( (procdir = opendir( "/proc" )) == NULL ) return;
  
  while( (procdirp = readdir(procdir)) != NULL ){
    
    /* Only look at this file if it's a proc id; that is, all numbers */
    if( strtok(procdirp->d_name, "0123456789") != NULL ){ continue; }
      
    /* Construct path of the form /proc/proc_number */
    strcpy( pathbuf, "/proc/"); 
    strcat( pathbuf, procdirp->d_name );
    strcat( pathbuf, "/psinfo" );
      
    if( (psdata = open( pathbuf, O_RDONLY )) == -1 ) continue;
	
    read(psdata, (void *) &psbuf, sizeof(struct psinfo) );
    close(psdata);

    /* translate process state */
    pr_state = psbuf.pr_lwp.pr_state;
    switch( pr_state )
      {
      case SSLEEP: 
	strcpy(state, SLEEP);
	break;
      case SRUN:
	strcpy(state, RUN);
	break;
      case SSTOP:
	strcpy(state, STOP);
	break;
      case SIDL:
	strcpy(state, IDLE);
	break;
      case SONPROC:
	strcpy(state, ONPROC);
	break;
      }

    bless_into_proc( Format,           
		     Fields,
		     
		     psbuf.pr_uid,           /* uid */
		     psbuf.pr_gid,           /* gid */
		     psbuf.pr_pid,           /* pid */
		     psbuf.pr_ppid,          /* ppid */
		     psbuf.pr_pgid,          /* pgrp */ 
		     psbuf.pr_sid,           /* sess */
		     psbuf.pr_lwp.pr_pri,    /* priority */   
		     psbuf.pr_lwp.pr_nice,   /* nice */
		     psbuf.pr_ttydev,        /* ttynum */
		     psbuf.pr_flag,          /* flags */
		     psbuf.pr_time.tv_sec,   /* time */
		     psbuf.pr_time.tv_nsec,   /* time nanosec */
		     psbuf.pr_size * getpagesize(),    /* size (bytes) */
		     psbuf.pr_rssize * getpagesize(),  /* rss (bytes)  */
		     psbuf.pr_lwp.pr_wchan,   /* wchan */ 
		     psbuf.pr_fname,         /* fname */
		     psbuf.pr_start.tv_sec,  /* start */
		     state,                  /* state */
		     psbuf.pr_lwp.pr_onpro,  /* on which processor */
		     psbuf.pr_psargs         /* cmndline */ 
		    );
    
  }
  closedir(procdir);
}
