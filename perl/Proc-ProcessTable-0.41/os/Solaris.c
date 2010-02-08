
#include "os/Solaris.h"

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
#include <sys/types.h>
void OS_get_table(){
  DIR *procdir;
  char p[sizeof(struct dirent) + 1024] ;
  struct dirent *procdirp = (struct dirent *)p;
  int psdata;
  char pathbuf[MAXPATHLEN];
  DIR *procdirlwp;
  struct dirent *procdirlwpp;
  char pathbuflwp[MAXPATHLEN];


#if defined(PROC_FS)
  struct psinfo psbuf;
#else
  struct prpsinfo psbuf;
#endif
  
  /* variables to hold some values for bless_into_proc */
  char state[20]; 
  char pctcpu[7];
  int numthr;
  char pctmem[7];

/* MR: quick hack for different readdir_r versions */ 
#if defined(_POSIX_PTHREAD_SEMANTICS)
 struct dirent *procdirp_r;
#endif
 
  if( (procdir = opendir( "/proc" )) == NULL ) return;
  
#if defined(_POSIX_PTHREAD_SEMANTICS)
  while( readdir_r(procdir, procdirp, &procdirp_r ) == 0 && procdirp_r != NULL ){
#else
  while( readdir_r(procdir, procdirp ) != NULL ){
#endif
    
    /* Only look at this file if it's a proc id; that is, all numbers */
    if( strtok(procdirp->d_name, "0123456789") != NULL ){ continue; }
      
    /* Construct path of the form /proc/proc_number */
    strcpy( pathbuf, "/proc/"); 
    strcat( pathbuf, procdirp->d_name );
    strcpy( pathbuflwp, pathbuf);
    strcat( pathbuflwp, "/lwp/");
    if( (procdirlwp = opendir( pathbuflwp )) != NULL ){
      numthr = 0;
      while( (procdirlwpp = readdir(procdirlwp)) != NULL ){
      numthr++;
      }
      closedir(procdirlwp);
    }
    numthr = numthr - 2;

      
#if defined(PROC_FS)
    strcat( pathbuf, "/psinfo" ); /* Solaris 2.6 has process info here */
#endif
      
    if( (psdata = open( pathbuf, O_RDONLY )) == -1 ) continue;
	
#if defined(PROC_FS)
    read(psdata, (void *) &psbuf, sizeof(struct psinfo) );
#else
    if( ioctl(psdata, PIOCPSINFO, &psbuf) == -1 ) continue; 
#endif
    close(psdata);

    /* translate process state */
#if defined(PROC_FS)
    switch( psbuf.pr_lwp.pr_state )
#else
    switch( psbuf.pr_state)
#endif
      {
      case SSLEEP: 
	strcpy(state, SLEEP);
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
      case SONPROC:
	strcpy(state, ONPROC);
	break;
      }

    /* These seem to be 2 bytes, 1st byte int part, 2nd byte fractional */
    /* Perl can handle stringy numbers of the form 1.5 */
    sprintf( pctcpu,  "%5.2f%", ((double)psbuf.pr_pctcpu)/0x8000*100  );
    sprintf( pctmem,  "%5.2f%", ((double)psbuf.pr_pctmem)/0x8000*100 );
    
    bless_into_proc( Format,           
		     Fields,
		     
		     psbuf.pr_uid,           /* uid */
		     psbuf.pr_gid,           /* gid */
		     psbuf.pr_euid,          /* euid */
		     psbuf.pr_egid,          /* egid */
		     psbuf.pr_pid,           /* pid */
		     psbuf.pr_ppid,          /* ppid */
#if defined(PROC_FS)
		     psbuf.pr_pgid,          /* pgrp */ 
#else
		     psbuf.pr_pgrp,          /* pgrp */ 
#endif
		     psbuf.pr_sid,           /* sess */
#if defined(PROC_FS)
		     psbuf.pr_lwp.pr_pri,    /* priority */   
		     psbuf.pr_lwp.pr_nice,   /* nice */
#else
		     psbuf.pr_pri,           /* priority */   
		     psbuf.pr_nice,          /* nice */
#endif
		     psbuf.pr_ttydev,        /* ttynum */
		     psbuf.pr_flag,          /* flags */
		     psbuf.pr_time.tv_sec,   /* time */
		     psbuf.pr_ctime.tv_sec,  /* ctime */
#if defined(PROC_FS)
		     psbuf.pr_time.tv_nsec,   /* time nanosec */
		     psbuf.pr_ctime.tv_nsec,  /* ctime nanosec */
		     psbuf.pr_size * 1024,    /* size (bytes) */
		     psbuf.pr_rssize * 1024,  /* rss (bytes)  */
		     psbuf.pr_lwp.pr_wchan,   /* wchan */ 
#else
		     psbuf.pr_bysize,        /* size (bytes) */
		     psbuf.pr_byrssize,      /* rss (bytes)  */
		     psbuf.pr_wchan,         /* wchan */
#endif

		     psbuf.pr_fname,         /* fname */
		     psbuf.pr_start.tv_sec,  /* start */
		     pctcpu,                 /* pctcpu */
		     state,                  /* state */
#if defined(PROC_FS)
		     psbuf.pr_lwp.pr_onpro,  /* on which processor */
#endif
		     pctmem,                 /* pctmem */
		     psbuf.pr_psargs,        /* cmndline */ 
		     numthr                  /* numthr */ 
		    );
    
  }
  closedir(procdir);
}
