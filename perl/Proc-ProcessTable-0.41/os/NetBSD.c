
#include "os/NetBSD.h"

/* Given a path to a /proc/XXX/status file and a pointer to a procstat
   struct, fill the struct */
struct procstat* get_procstat( char* path, struct procstat* prs){
  FILE* fp;

  if( (fp = fopen( path, "r" )) != NULL ){ 
    fscanf(fp, 
	   "%s %d %d %d %d %d,%d %s %d,%d %d,%d %d,%d %s %d %d %d,%d,%s",
	   &prs->comm, 
	   &prs->pid, 
	   &prs->ppid,
	   &prs->pgid,
	   &prs->sid, 
	   &prs->tdev_maj, 
	   &prs->tdev_min, 
	   &prs->flags, 
	   &prs->start,
	   &prs->start_mic,
	   &prs->utime,
	   &prs->utime_mic,
	   &prs->stime, 
	   &prs->stime_mic, 
	   &prs->wchan,
	   &prs->euid,
	   &prs->ruid,
	   &prs->rgid,
	   &prs->egid,
	   &prs->groups
	   );
    fclose(fp);
    return(prs);
  }
  else{
    return(NULL);
  }
}

/* Make sure /proc is mounted and initialize some global vars */
char* OS_initialize(){
  char cbuf[1024]; 
  char bbuf[32]; 
  struct statfs sfs;

  static char* no_proc = "/proc unavailable";
  if( statfs("/proc", &sfs) == -1 ){
    return no_proc;
  }

  return NULL;
}

void OS_get_table(){
  DIR *procdir;
  struct dirent *procdirp;
  FILE *fp;
  char pathbuf[PATH_MAX];

  /* for bless_into_proc */
  struct procstat prs; 
  static char format[F_LASTFIELD + 1];
  char cmndline[ARG_MAX];

  if( (procdir = opendir("/proc")) == NULL ){
    return;
  }

  /* Iterate through all the process entries (numeric) under /proc */
  while( (procdirp = readdir(procdir)) != NULL ){

    /* Only look at this file if it's a proc id; that is, all numbers */
    if( strtok(procdirp->d_name, "0123456789") == NULL ){
      /* zero our format */
      strcpy(format, Defaultformat);

      sprintf(pathbuf, "%s%s", "/proc/", procdirp->d_name);

      /* get stuff out of /proc/PROC_ID/status */
      memset(&prs, 0, sizeof(prs));
      strcat( pathbuf, "/status" );
      if (get_procstat(pathbuf, &prs) != NULL) {
	int i;
	char utime[20];
	double utime_f;
	char stime[20];
	double stime_f;
	char time[20];
	double time_f;
	char start[20];
	double start_f;
	char *ttydev;
	int ttynum;

	utime_f = prs.utime + prs.utime_mic/1000000;
	stime_f = prs.stime + prs.stime_mic/1000000;
	time_f  = utime_f + stime_f;
	start_f = prs.start + prs.start_mic/1000000;

	sprintf(utime, "%f", utime_f);
	sprintf(stime, "%f", stime_f);
	sprintf(time, "%f", time_f);
	sprintf(start, "%f", start_f);

	ttynum = makedev(prs.tdev_maj, prs.tdev_min);
	ttydev = devname(ttynum, S_IFCHR);
	if (ttydev == NULL) ttydev = "??";

	/* get stuff out of /proc/PROC_ID/cmdline */
	sprintf(pathbuf, "%s%s%s", "/proc/", procdirp->d_name, "/cmdline");
	if( (fp = fopen( pathbuf, "r" )) != NULL ){ 
	  size_t got;
	  if( (got = fread(cmndline, sizeof(char), ARG_MAX, fp)) > 0 ){
	    size_t i;
	    for(i = 0; i < got; i++){
	      if( cmndline[i] == '\0' ) cmndline[i] = ' ';
	    }
	    cmndline[got] = '\0'; /* necessary? */
	    
	    format[F_CMNDLINE] = tolower(format[F_CMNDLINE]);
	  }
	  fclose(fp);
	}

	bless_into_proc( format,
			 Fields,
			 prs.ruid,
			 prs.rgid,
			 prs.pid,
			 prs.ppid,
			 prs.pgid,
			 prs.sid,
			 prs.flags,
			 utime,
			 stime,
			 time,
			 prs.wchan,
			 start,
			 prs.euid,
			 prs.egid,
			 prs.comm,
			 prs.wchan,
			 ttydev,
			 ttynum,
			 cmndline
			 );
      }
    }
  }
  closedir(procdir);
}

