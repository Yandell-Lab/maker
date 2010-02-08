/* 

	Adapted from ps.cc by J Robert Ray <jrray@jrray.org>


   ps.cc

   Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002 Red Hat, Inc.

This file is part of Cygwin.

This software is a copyrighted work licensed under the terms of the
Cygwin license.  Please consult the file "CYGWIN_LICENSE" for
details. */

#include <stdio.h>
#include <windows.h>
#include <time.h>
#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <pwd.h>
#include <sys/cygwin.h>
#include <tlhelp32.h>
#include <psapi.h>

#include "os/cygwin.h"

typedef BOOL (WINAPI *ENUMPROCESSMODULES)(
  HANDLE hProcess,      // handle to the process
  HMODULE * lphModule,  // array to receive the module handles
  DWORD cb,             // size of the array
  LPDWORD lpcbNeeded    // receives the number of bytes returned
);

typedef DWORD (WINAPI *GETMODULEFILENAME)(
  HANDLE hProcess,
  HMODULE hModule,
  LPTSTR lpstrFileName,
  DWORD nSize
);

typedef HANDLE (WINAPI *CREATESNAPSHOT)(
    DWORD dwFlags,
    DWORD th32ProcessID
);

// Win95 functions
typedef BOOL (WINAPI *PROCESSWALK)(
    HANDLE hSnapshot,
    LPPROCESSENTRY32 lppe
);

typedef struct external_pinfo external_pinfo;

ENUMPROCESSMODULES myEnumProcessModules;
GETMODULEFILENAME myGetModuleFileNameEx;
CREATESNAPSHOT myCreateToolhelp32Snapshot;
PROCESSWALK myProcess32First;
PROCESSWALK myProcess32Next;

static int init_win_result = FALSE;

static BOOL WINAPI dummyprocessmodules (
  HANDLE hProcess,      // handle to the process
  HMODULE * lphModule,  // array to receive the module handles
  DWORD cb,             // size of the array
  LPDWORD lpcbNeeded    // receives the number of bytes returned
)
{
  lphModule[0] = (HMODULE) *lpcbNeeded;
  *lpcbNeeded = 1;
  return 1;
}

static DWORD WINAPI GetModuleFileNameEx95 (
  HANDLE hProcess,
  HMODULE hModule,
  LPTSTR lpstrFileName,
  DWORD n
)
{
  HANDLE h;
  DWORD pid = (DWORD) hModule;
  PROCESSENTRY32 proc;

  h = myCreateToolhelp32Snapshot (TH32CS_SNAPPROCESS, 0);
  if (!h)
    return 0;

  proc.dwSize = sizeof (proc);
  if (myProcess32First(h, &proc))
    do
      if (proc.th32ProcessID == pid)
	{
	  CloseHandle (h);
	  strcpy (lpstrFileName, proc.szExeFile);
	  return 1;
	}
    while (myProcess32Next (h, &proc));
  CloseHandle (h);
  return 0;
}

int
init_win ()
{
  OSVERSIONINFO os_version_info;
  HMODULE h;

  memset (&os_version_info, 0, sizeof os_version_info);
  os_version_info.dwOSVersionInfoSize = sizeof (OSVERSIONINFO);
  GetVersionEx (&os_version_info);

  if (os_version_info.dwPlatformId == VER_PLATFORM_WIN32_NT)
    {
      h = LoadLibrary ("psapi.dll");
      if (!h)
	return 0;
      myEnumProcessModules = (ENUMPROCESSMODULES) GetProcAddress (h, "EnumProcessModules");
      myGetModuleFileNameEx = (GETMODULEFILENAME) GetProcAddress (h, "GetModuleFileNameExA");
      if (!myEnumProcessModules || !myGetModuleFileNameEx)
	return 0;
      return 1;
    }

  h = GetModuleHandle("KERNEL32.DLL");
  myCreateToolhelp32Snapshot = (CREATESNAPSHOT)GetProcAddress (h, "CreateToolhelp32Snapshot");
  myProcess32First = (PROCESSWALK)GetProcAddress (h, "Process32First");
  myProcess32Next  = (PROCESSWALK)GetProcAddress (h, "Process32Next");
  if (!myCreateToolhelp32Snapshot || !myProcess32First || !myProcess32Next)
    return 0;

  myEnumProcessModules = dummyprocessmodules;
  myGetModuleFileNameEx = GetModuleFileNameEx95;
  return 1;
}

static char *
start_time (external_pinfo *child)
{
  time_t st = child->start_time;
  time_t t = time (NULL);
  static char stime[40] = {'\0'};
  char now[40];

  strncpy (stime, ctime (&st) + 4, 15);
  strcpy (now, ctime (&t) + 4);

  if ((t - st) < (24 * 3600))
    return (stime + 7);

  stime[6] = '\0';

  return stime;
}

#define FACTOR (0x19db1ded53ea710LL)
#define NSPERSEC 10000000LL

/* Convert a Win32 time to "UNIX" format. */
long __stdcall
to_time_t (FILETIME *ptr)
{
  /* A file time is the number of 100ns since jan 1 1601
     stuffed into two long words.
     A time_t is the number of seconds since jan 1 1970.  */

  long rem;
  long long x = ((long long) ptr->dwHighDateTime << 32) + ((unsigned)ptr->dwLowDateTime);
  x -= FACTOR;                  /* number of 100ns between 1601 and 1970 */
  rem = x % ((long long)NSPERSEC);
  rem += (NSPERSEC / 2);
  x /= (long long) NSPERSEC;            /* number of 100ns in a second */
  x += (long long) (rem / NSPERSEC);
  return x;
}

void
OS_get_table()
{
  external_pinfo *p;
  int uid;
  cygwin_getinfo_types query = CW_GETPINFO;
  char ch;
  int pid;
  char *pstate;
  char pname[MAX_PATH];
  HMODULE hm[1000];
  char uname[128];
  char *fields;

  uid = getuid ();

  (void) cygwin_internal (CW_LOCK_PINFO, 1000);

  if (query == CW_GETPINFO && !init_win_result)
    query = CW_GETPINFO;

  for (pid = 0;
       (p = (external_pinfo *) cygwin_internal (query, pid | CW_NEXTPID));
       pid = p->pid)
    {
      pstate = " ";
      if (p->process_state & PID_STOPPED)
	pstate = "stopped";
      else if (p->process_state & PID_TTYIN)
	pstate = "ttyin";
      else if (p->process_state & PID_TTYOU)
	pstate = "ttyout";

      if (p->process_state & (PID_ORPHANED | PID_EXITED))
        strcpy (pname, "<defunct>");
      else if (p->ppid)
	{
	  char *s;
	  pname[0] = '\0';
	  cygwin_conv_to_posix_path (p->progname, pname);
	  s = strchr (pname, '\0') - 4;
	  if (s > pname && strcasecmp (s, ".exe") == 0)
	    *s = '\0';
	}
      else if (query == CW_GETPINFO)
	{
	  DWORD n;
	  FILETIME ct, et, kt, ut;
	  HANDLE h = OpenProcess (PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
	  			  FALSE, p->dwProcessId);
	  if (!h)
	    continue;
	  n = p->dwProcessId;
	  if (!myEnumProcessModules (h, hm, sizeof (hm), &n))
	    n = 0;
	  if (!n || !myGetModuleFileNameEx (h, hm[0], pname, MAX_PATH))
	    strcpy (pname, "*** unknown ***");
	  if (GetProcessTimes (h, &ct, &et, &kt, &ut))
	    p->start_time = to_time_t (&ct);
	  CloseHandle (h);
	}

        {
          struct passwd *pw;

          if ((pw = getpwuid (p->version >= EXTERNAL_PINFO_VERSION_32_BIT ?
	                      p->uid32 : p->uid)))
            strcpy (uname, pw->pw_name);
          else
            sprintf (uname, "%u", (unsigned)
	             (p->version >= EXTERNAL_PINFO_VERSION_32_BIT ?
		      p->uid32 : p->uid));
        }

	if (query == CW_GETPINFO) {
		fields = "iiiiisiis";
	} else {
		fields = "iiiiisIis";
	}

	bless_into_proc(fields, Fields, 

		p->version >= EXTERNAL_PINFO_VERSION_32_BIT ? p->uid32 : p->uid,
		p->pid,
		p->ppid,
		p->pgid,
		p->dwProcessId,
		pname,
		p->start_time,
		p->ctty,
		pstate
	);

    }
  (void) cygwin_internal (CW_UNLOCK_PINFO);
}

char*
OS_initialize()
{
	init_win_result = init_win();

	return NULL;
}

