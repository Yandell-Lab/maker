#include "os/FreeBSD-kvm.h"
#include <unistd.h>
#include <sys/sysctl.h>
#include <sys/proc.h>
#include <sys/user.h>                                                          
#include <kvm.h>



/* Make sure /proc is mounted and initialize some global vars */
char* OS_initialize(){
  return NULL;
}

void OS_get_table(){

  kvm_t *kd;
  char errbuf[2048];
  struct kinfo_proc *procs;          /* array of processes */
  int count;                         /* returns number of processes */
  int i;
  char ** argv;

  char state[20];
  char start[20];
  char time[20];
  char flag[20];
  char sflag[20];

  static char format[128];
  char cmndline[ARG_MAX];
  int priority;

  /* Open the kvm interface, get a descriptor */
  if ((kd = kvm_open(NULL, NULL, NULL, 0, errbuf)) == NULL) {
     fprintf(stderr, "kvm_open: %s\n", errbuf);
     ppt_croak("kvm_open: ", errbuf);
  }  

  /* Get the list of processes. */
  if ((procs = kvm_getprocs(kd, KERN_PROC_ALL, 0, &count)) == NULL) {
     kvm_close(kd);
     fprintf(stderr, "kvm_getprocs: %s\n", kvm_geterr(kd));
     ppt_croak("kvm_getprocs: ", kvm_geterr(kd));
  }

  /* Iterate through the processes in kinfo_proc, sending proc info */
  /* to bless_into_proc for each proc */

  for (i=0; i < count; i++) {
     static struct pstats ps;
     static struct session seslead;
     strcpy(format, Defaultformat);

     // retrieve the arguments
     cmndline[0] = '\0';
     argv = kvm_getargv(kd, (const struct kinfo_proc *) &(procs[i]) , 0);
     if (argv) {
       int j = 0;
       while (argv[j] && strlen(cmndline) <= ARG_MAX) {
         strcat(cmndline, argv[j]);
         strcat(cmndline, " ");
         j++;
       }
     }

     switch (procs[i].ki_stat) {
        case SSTOP:
            strcpy(state, "stop");
            break;
        case SSLEEP:
            strcpy(state, "sleep");
            break;
        case SRUN:
            strcpy(state, "run");
            break;
        case SIDL:
            strcpy(state, "idle");
            break;
        case SWAIT:
            strcpy(state, "wait");
            break;
        case SLOCK:
            strcpy(state, "lock");
            break;
        case SZOMB:
            strcpy(state, "zombie");
            break;
        default:
            strcpy(state, "???");
            break;
     }


     sprintf(start, "%d.%d", procs[i].ki_start.tv_sec, procs[i].ki_start.tv_usec);
     sprintf(time, "%.6f", procs[i].ki_runtime/1000000.0);
     sprintf(flag, "0x%04x", procs[i].ki_flag);
     sprintf(sflag, "0x%04x", procs[i].ki_sflag);

     bless_into_proc( format,
                      Fields,

                      procs[i].ki_pid,
                      procs[i].ki_ppid,
                      procs[i].ki_ruid,
                      procs[i].ki_uid,
                      procs[i].ki_rgid,
                      procs[i].ki_pgid,
                      procs[i].ki_sid,

                      flag,
                      sflag,

                      start, 
                      time,

                      procs[i].ki_wmesg,
                      state,

                      "??", // will be resolved automatically
                      procs[i].ki_tdev,

                      procs[i].ki_comm,
                      cmndline,

                      procs[i].ki_pri.pri_user,
                      procs[i].ki_nice,

                      procs[i].ki_size,              // virtual size
                      procs[i].ki_rssize,              // current resident set size in pages
                      procs[i].ki_tsize,               // text size (pages) XXX
                      procs[i].ki_dsize,               // data size (pages) XXX
                      procs[i].ki_ssize               // stack size (pages)   
              );


   }
   if (kd) kvm_close(kd);
}

