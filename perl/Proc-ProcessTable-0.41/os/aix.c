/*****************************************************************************/
/* Copyright (c) 1998, David Paquet. All rights reserved.                    */
/* This file is free software; you can redistribute it and/or modify it      */
/* under the same terms as Perl itself.                                      */
/*****************************************************************************/

#include "os/aix.h"


char* OS_initialize() {
  
  struct CuAt*      obj;
  int               how_many;
  

  Sysmem = 0;

  /*
   * Get the real memory size via ODM
   *
   */


  if (odm_initialize() == 0) {
    obj = (struct CuAt*)getattr ("sys0", "realmem", 0, &how_many);
    Sysmem = strtoull(obj->value, 0, 10);
    odm_terminate();
  }

  else {
    printf("BIG PROLEM !\n");
  }
  
  Sysmem = Sysmem * 1024;
  
  /*
   * Get The number of processors
   *
   */
   ProcessNumber = sysconf(_SC_NPROCESSORS_ONLN);
   if ( ProcessNumber == -1 ) {
     ProcessNumber = 1;
   }


  /*
   * Get the page size in bytes
   *
   */
  
  PageSize = getpagesize();
  
  return NULL;
}


void OS_get_table() {
  int                  i, proc_nb;
  struct procinfo      pr_buff[MAX_PROCS];
  struct userinfo      uinfo;
  char                 format[F_FLAST + 1];
  char                 wchan[15], pctcpu[7], pctmem[7], state[10];
  char                 Args[MAXARGLN+1], Arglist[MAXARGLN+1], Comm[MAXARGLN+1];
  int                  argcount;
  struct timeval       now_tval;
  double               utime, stime, cutime, cstime, now;
  
  
  strcpy(format, Fullformat);
  
  proc_nb = getproc(pr_buff, MAX_PROCS, sizeof(struct procinfo));
  
  gettimeofday(&now_tval, (void *)NULL);
  now = (double)now_tval.tv_sec + (double) now_tval.tv_usec / 1000000.0;
  
  for(i=0; i<proc_nb; i++) {
    
    if ( pr_buff[i].pi_wchan != NULL ) {
      sprintf(wchan, "%p", pr_buff[i].pi_wchan);
    }
    
    if ( pr_buff[i].pi_stat == SNONE ) {
      continue;
    }
    
    
    switch (pr_buff[i].pi_stat){
    case SSLEEP :
      strcpy(state, SLEEP);
      break;
    case SRUN :
      strcpy(state, RUN);
      break;
    case SIDL :
      strcpy(state, IDLE);
      break;
    case SZOMB :
      strcpy(state, ZOMBIE);
      break;
    case SSTOP :
      strcpy(state, STOP);
      break;
    case SACTIVE :
      strcpy(state, ACTIVE);
      break;
    default:
      format[F_STAT] = 'S';
    }

    
    if ( state == ZOMBIE 
	 || getuser(&pr_buff[i], sizeof(struct procinfo),
		    &uinfo,      sizeof(struct userinfo) ) < 0 ) {
      
      bless_into_proc( Zombformat,
		       ZombFields,
		       pr_buff[i].pi_pid,
		       pr_buff[i].pi_ppid,
		       pr_buff[i].pi_sid,
		       pr_buff[i].pi_pgrp,
		       pr_buff[i].pi_uid,
		       pr_buff[i].pi_suid,
		       pr_buff[i].pi_pri,
		       pr_buff[i].pi_nice,
		       pr_buff[i].pi_cpu,
		       /* pr_buff[i].pi_stat, */
		       state,
		       pr_buff[i].pi_flag,
		       wchan,
		       pr_buff[i].pi_wtype,
		       pr_buff[i].pi_adspace,
		       pr_buff[i].pi_majflt,
		       pr_buff[i].pi_minflt,
		       pr_buff[i].pi_utime,
		       pr_buff[i].pi_stime,
		       pr_buff[i].pi_size * PageSize );
      
      continue;
    } 
    
    /*
     * Command line args processing
     * 
     */

    Args[0] = '\0';

    if ( pr_buff[i].pi_flag & SKPROC ) {
      if ( pr_buff[i].pi_pid == 0 ) {
	strcpy(Args, "kproc (swapper)");
	strcpy(Comm, "kproc (swapper)");
      }
      else {
	sprintf(Args, "kproc (%s)", uinfo.ui_comm);
	sprintf(Comm, "kproc (%s)", uinfo.ui_comm);
      }
    }
    else {
      strncpy(Comm, uinfo.ui_comm, MAXARGLN);
      Comm[MAXARGLN] = '\0';
      if (getargs(&pr_buff[i], sizeof(struct procinfo), Arglist, MAXARGLN) < 0) {
	sprintf(Args, "%s", uinfo.ui_comm);
      }
      else {
	/* Returns a succession of strings seperated by a null
           characters, 2 nulls for the end (see getargs info/man page) */
	argcount = -1;
	while (++argcount < MAXARGLN) {
	  /* Copy everything but replace the end of the arg with a
             space */
	  if (Arglist[argcount] != '\0') {
	    Args[argcount] = Arglist[argcount];
	  }
	  else {
	    /* Is this the last arg, then hop out of the loop */
	    if (Arglist[argcount+1] == '\0') {
	      /* Terminate the arguments */
	      Args[argcount] = '\0';
	      break;
	    }
	    /* Seperate arguments with a space */
	    Args[argcount] = ' ';
	  }
	}
      }
    }
    
    
    /*
     * Convert time values into seconds
     * 
     */

    
    utime = uinfo.ui_ru.ru_utime.tv_sec + 
            (double) uinfo.ui_ru.ru_utime.tv_usec / 1000000.0;
    
    stime = uinfo.ui_ru.ru_stime.tv_sec + 
            (double) uinfo.ui_ru.ru_stime.tv_usec / 1000000.0;

    cutime = uinfo.ui_cru.ru_utime.tv_sec + 
            (double) uinfo.ui_cru.ru_utime.tv_usec / 1000000.0;
    
    cstime = uinfo.ui_cru.ru_stime.tv_sec + 
            (double) uinfo.ui_cru.ru_stime.tv_usec / 1000000.0;

    /*
     * percentage calculation
     * 
     */
        
    pctcpu[0] = pctmem[0]= '\0';      
/* compute %CPU in SMP environment */    
    sprintf( pctcpu, 
	     "%3.2f",
	     ((utime + stime ) * 100 / ( now - uinfo.ui_start )) / ProcessNumber );
    
    if ( Sysmem == 0 ) {
      format[F_PRM] = 'S';
    }
    else {
      sprintf( pctmem,
	       "%3.2f",
	       (uinfo.ui_drss + uinfo.ui_trss) * PageSize * 100 / 
	       (float)Sysmem );
    }


    /*
     * Give it all to Perl
     * (we convert time values to hundredth of a second
     *  to keep in sync with the Linux version )
     */
    
    bless_into_proc( format,
		     FullFields,
		     pr_buff[i].pi_pid,
		     pr_buff[i].pi_ppid,
		     pr_buff[i].pi_sid,
		     pr_buff[i].pi_pgrp,
		     pr_buff[i].pi_uid,
		     pr_buff[i].pi_suid,
		     pr_buff[i].pi_pri,
		     pr_buff[i].pi_nice,
		     /* pr_buff[i].pi_cpu,  */
		     pctcpu,
		     /* pr_buff[i].pi_stat, */
		     state,
		     pr_buff[i].pi_flag,
		     wchan,
		     pr_buff[i].pi_wtype,
		     pr_buff[i].pi_adspace,
		     pr_buff[i].pi_majflt,
		     pr_buff[i].pi_minflt,
		     ((pr_buff[i].pi_size * PageSize) - uinfo.ui_tsize)/1024,
		     uinfo.ui_luid,
		     uinfo.ui_uid,
		     uinfo.ui_gid,
		     uinfo.ui_start,
		     (long) (utime  * 100),
		     (long) (stime  * 100),
		     (long) (cutime * 100),
		     (long) (cstime * 100),
		     uinfo.ui_tsize/1024,
		     uinfo.ui_ttyp,
		     uinfo.ui_ttyd,
		     uinfo.ui_ttympx,
		     ((uinfo.ui_drss + uinfo.ui_trss) * PageSize)/1024,
		     (uinfo.ui_trss * PageSize)/1024,
		     uinfo.ui_dvm,
		     /* uinfo.ui_prm,*/
		     pctmem,
		     Comm,
                     Args ); 
    
  } 
  
}

