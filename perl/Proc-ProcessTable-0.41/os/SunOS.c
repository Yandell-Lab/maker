/*
 * Copyright (c) 2001 by Shawn A. Clifford <shawn.a.clifford@lmco.com>
 * This file may be distributed under the same terms as Perl.
 *
 * Modification History:
 *
 * Who	When		Description
 * ---	----------	--------------------------------------------
 * SAC	30July2001	Original code
 */

#include "os/SunOS.h"

char* OS_initialize(void) {
    return NULL;
}

void OS_get_table(void) {
   struct proc *p;
   struct user *u;
   struct ucred cr;
   struct sess sess;
   char **arg;
   char **env;
   char cmdline[_POSIX_ARG_MAX];
   char fname[_POSIX_PATH_MAX];
   kvm_t *kd;
   int i, count;
   int ttynum;
        
   /* Open the kernel for reading */
   if ((kd = kvm_open(NULL, NULL, NULL, O_RDONLY, NULL)) == NULL) {
      ppt_croak("kvm_open: can't open kernel");
   }
   
   /*
    *  Loop over all processes
    */ 
   while ((p = kvm_nextproc(kd)) != NULL) {
     
      /* Get the u-area for this process or skip this process */
      if ((u = kvm_getu(kd, p)) == NULL) continue;

      /* Get the command line arguments for this process or skip */
      bzero(fname, sizeof(fname));
      if (kvm_getcmd(kd, p, u, &arg, &env) < 0) {
         sprintf(fname, "%s", u->u_comm);
         sprintf(cmdline, "%s", u->u_comm);
      } else {
         sprintf(fname, "%s", arg[0]);
         bzero(cmdline, sizeof(cmdline));
         i=0;
         while (arg[i] != NULL) {
            count = sizeof(cmdline) - strlen(cmdline) - 1;
            strncat(cmdline, arg[i++], count);
            if (arg[i] != NULL) strcat(cmdline, " ");
         }
      }
      
      /* Get the process credentials */
      if (kvm_read(kd, p->p_cred, &cr, sizeof(struct ucred)) < 0) continue;

      /* Get the session info */
      kvm_read(kd, p->p_sessp, &sess, sizeof(struct sess));
      ttynum = minor(sess.s_ttyd);
      if (major(sess.s_ttyd) == 0) ttynum = -1;
      
      /* Send if off to Perl */
      bless_into_proc(  Format,
			Fields,
			cr.cr_ruid,
			cr.cr_rgid,
			cr.cr_uid,
			cr.cr_gid,
			p->p_pid,
			p->p_ppid,
			p->p_pgrp,
			p->p_pri,
			p->p_flag,
			(p->p_dsize+p->p_ssize)*4,
			p->p_rssize*4,
			p->p_nice,
			p->p_time,
			fname,
			cmdline,
			p->p_cpticks,
			p->p_pctcpu,
			States[(int)p->p_stat],
			sess.s_sid,
			sess.s_sid,
			ttynum );
   }

   /* Close the kernel and exit */
   kvm_close(kd);
}
