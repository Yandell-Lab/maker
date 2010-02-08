/*alan.martin@oracle.com
 * Copyright (c) 1998 by Mike Romberg ( romberg@fsl.noaa.gov )
 * This file may be distributed under the same terms as Perl.
 *
 * This will probably only work under HPUX-10 or later.
 *
 * 8/26/99 Added "fname" field for consistency with other OS's - D. Urist
 *
 */

#define _PSTAT64 /* For David Good's 64-bit HPUX 11.0 patch */

#include <stdlib.h>
#include <sys/param.h>
#include <sys/pstat.h>

#define BURST 30 /* How many pstat structs to get per syscall */

extern void bless_into_proc(char* format, char** fields, ...);

static char *Format = 
"llllllllllllllllsllssllllsllslllllllllllllllllllllSSSllllSSSll";

static char *Fields[] = {
"uid",
"pid",
"ppid",
"dsize",
"tsize",
"ssize",
"nice",
"ttynum", /* */
"pgrp",
"pri",
"addr",
"cpu",
"utime",
"stime",
"start",
"flag",
"state",
"wchan",
"procnum",
"cmndline",
"fname", 
"time",
"cpticks",
"cptickstotal",
"fss",
"pctcpu",
"rssize",
"suid",
"ucomm", /* char */
"shmsize",
"mmsize",
"usize",
"iosize",
"vtsize",
"vdsize",
"vssize",
"vshmsize",
"vmmsize",
"vusize",
"viosize",
"minorfaults", /* ulong */
"majorfaults", /* ulong */
"nswap", /* ulong */
"nsignals", /* ulong */
"msgrcv", /* ulong */
"msgsnd", /* ulong */
"maxrss",
"sid",
"schedpolicy",
"ticksleft",
"rdir", /* */
"cdir", /* */
"text", /* */
"highestfd",
"euid",
"egid",
"ioch",
"usercycles", /* */
"systemcycles", /* */
"interruptcycles", /* */
"gid",
"lwpid"
};

static char *States[] = { 
"", "sleep", "run", "stop", "zombie", "uwait", "other" 
};

char* OS_initialize()
    {
    return NULL;
    }

void OS_get_table()
    {
    struct pst_status pst[BURST];

    int i, count;
    int idx = 0;
    char buff[256]; /* used to format %cpu which is the only float. */

    while ((count = pstat_getproc(pst, sizeof(pst[0]), BURST, idx)) > 0)
        {
        for (i = 0; i < count; i++)
            {
            sprintf(buff, "%f", pst[i].pst_pctcpu * 100);
            bless_into_proc(Format, Fields,
              (long) pst[i].pst_uid,
              (long) pst[i].pst_pid,
              (long) pst[i].pst_ppid,
              (long) pst[i].pst_dsize,
              (long) pst[i].pst_tsize,
              (long) pst[i].pst_ssize,
              (long) pst[i].pst_nice,
              (long) makedev(pst[i].pst_term.psd_major, pst[i].pst_term.psd_minor),
              (long) pst[i].pst_pgrp,
              (long) pst[i].pst_pri,
              (long) pst[i].pst_addr,
              (long) pst[i].pst_cpu,
              (long) pst[i].pst_utime,
              (long) pst[i].pst_stime,
              (long) pst[i].pst_start,
              (long) pst[i].pst_flag,
              States[pst[i].pst_stat],
              (long) pst[i].pst_wchan,
              (long) pst[i].pst_procnum,
              pst[i].pst_cmd,
              pst[i].pst_ucomm,
              (long) pst[i].pst_cptickstotal/100,
              (long) pst[i].pst_cpticks,
              (long) pst[i].pst_cptickstotal,
              (long) pst[i].pst_fss,
              buff,
              (long) pst[i].pst_rssize,
              (long) pst[i].pst_suid,
              pst[i].pst_ucomm,
              (long) pst[i].pst_shmsize,
              (long) pst[i].pst_mmsize,
              (long) pst[i].pst_usize,
              (long) pst[i].pst_iosize,
              (long) pst[i].pst_vtsize,
              (long) pst[i].pst_vdsize,
              (long) pst[i].pst_vssize,
              (long) pst[i].pst_vshmsize,
              (long) pst[i].pst_vmmsize,
              (long) pst[i].pst_vusize,
              (long) pst[i].pst_viosize,
              (long) pst[i].pst_minorfaults,
              (long) pst[i].pst_majorfaults,
              (long) pst[i].pst_nswap,
              (long) pst[i].pst_nsignals,
              (long) pst[i].pst_msgrcv,
              (long) pst[i].pst_msgsnd,
              (long) pst[i].pst_maxrss,
              (long) pst[i].pst_sid,
              (long) pst[i].pst_schedpolicy,
              (long) pst[i].pst_ticksleft,
              "",
              "",
              "",
              (long) pst[i].pst_highestfd,
              (long) pst[i].pst_euid,
              (long) pst[i].pst_egid,
              (long) pst[i].pst_ioch,
              "",
              "",
              "",
              (long) pst[i].pst_gid,
              (long) pst[i].pst_lwpid);
            }
        idx = pst[count-1].pst_idx + 1;
        }
    }


