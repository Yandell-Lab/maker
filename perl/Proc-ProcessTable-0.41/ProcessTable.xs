#ifdef	PROCESSTABLE_THREAD
#define __REENTRANT
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#ifdef	PROCESSTABLE_THREAD
#include <pthread.h>
#endif

/*  As of version 5.005_something it seems sv_undef has been
supplanted by PL_sv_undef. */
#ifdef sv_undef
#define PL_sv_undef sv_undef
#endif

/* dTHX was used in perl 5.005 */
#ifndef dTHX
#define dTHX dTHR
#endif

/********************/
/* General includes */
/********************/
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <stdarg.h>

/* prototypes to make the compiler shut up */
void ppt_warn(const char*, ...);
void ppt_die(const char*, ...);
void store_ttydev(HV*, unsigned long);
void bless_into_proc(char* , char**, ...);
void OS_get_table();
char* OS_initialize();

char** Fields = NULL; 
int Numfields;

/* Cache a pointer the TTY device number hash for quick lookups */
HV* Ttydevs;

/* This holds a pointer to the list of process objects we will build */
AV* Proclist;

/* Our local varargs warn which can be called as extern by code
 * that doesn't know Perl internals (and thus doesn't have a
 * warn() defined).
 *
 * I think vwarn() and vcroak() have been changed to warn() and 
 * croak() in perl 5.8?? warn and croak exist in 5.6, but don't 
 * seem to accept format args.
 */
void ppt_warn(const char *pat, ...) {
    dTHX;
    va_list args;
    va_start(args, pat);
    vwarn(pat, &args);
    va_end(args);
}

/* same with croak */
void ppt_croak(const char *pat, ...) {
    dTHX;
    va_list args;
    va_start(args, pat);
    vcroak(pat, &args);
    va_end(args);
}

/* Look up the tty device, given the ttynum and store it */
void store_ttydev( HV* myhash, unsigned long ttynum ){
  SV** ttydev;
  char ttynumbuf[1024]; 
  
  sprintf(ttynumbuf, "%lu", ttynum);
  if( 
     Ttydevs != NULL &&
     (ttydev = hv_fetch(Ttydevs, ttynumbuf, strlen(ttynumbuf), 0)) != NULL 
     ){
    hv_store(myhash, "ttydev", strlen("ttydev"), newSVsv(*ttydev), 0); 
  }
  else{
    /* hv_store(myhash, "ttydev", strlen("ttydev"), &PL_sv_undef, 0); */ 

/*     Stuff an empty string into the hash if there is no tty; this */
/*     way the ttydev method won't return undef for nonexistent ttys. I'm */
/*     not sure if this is really the right behavior... */

    hv_store(myhash, "ttydev", strlen("ttydev"), newSVpv("",0), 0); 

  }
}


/**********************************************************************/
/* This gets called by OS-specific get_table                          */
/* format specifies what types are being passed in, in a string       */
/* containing these specifiers:                                       */
/*   S    ignore this string                                          */
/*   s    string                                                      */
/*   I    ignore this int                                             */
/*   i    int                                                         */
/*   L    ignore this long                                            */
/*   l    long                                                        */
/*   J    ignore this long-long                                       */
/*   j    long-long                                                   */
/*   U 	  ignore this unsigned                                        */
/*   u 	  unsigned                                                    */
/*   V    perl scalar value                                           */
/* fields is an array of pointers to field names                      */
/* following that is a var args list of field values                  */
/**********************************************************************/
void bless_into_proc(char* format, char** fields, ...){
  va_list args;
  char* key;
  char* s_val;
  SV *SV_val;
  int i_val;
  unsigned u_val;
  long l_val;
  long long ll_val;

  HV* myhash;
  SV* ref;
  HV* mystash;
  SV* blessed;

  /* Blech */
  if(Fields == NULL){
    Fields = fields; 
    Numfields = strlen(format);
  }

  myhash = newHV(); /* create a perl hash */

  va_start(args, fields);
  while( *format ){
    key = *fields; 
    switch(*format)
      {
      case 'S': /* ignore; creates an undef value for this key in the hash */
	va_arg(args, char *);
	hv_store(myhash, key, strlen(key), &PL_sv_undef, 0);
	break;
      case 's':  /* string */
	s_val = va_arg(args, char *);
	hv_store(myhash, key, strlen(key), newSVpv(s_val, strlen(s_val)), 0);
	break;

      case 'I':  /* ignore; creates an undef value for this key in the hash */
	va_arg(args, int);
	hv_store(myhash, key, strlen(key), &PL_sv_undef, 0);
	break;
      case 'i':  /* int */
	i_val = va_arg(args, int);
	hv_store(myhash, key, strlen(key), newSViv(i_val), 0);

	/* Look up and store the tty if this is ttynum */
	if( !strcmp(key, "ttynum") ) store_ttydev( myhash, i_val );
	break;

      case 'U':  /* ignore; creates an undef value for this key in the hash */
	va_arg(args, unsigned );
	hv_store(myhash, key, strlen(key), &PL_sv_undef, 0);
	break;
      case 'u':  /* int */
	u_val = va_arg(args, unsigned);
	hv_store(myhash, key, strlen(key), newSVuv(u_val), 0);
	break;

      case 'L':  /* ignore; creates an undef value for this key in the hash */
	va_arg(args, long);
	hv_store(myhash, key, strlen(key), &PL_sv_undef, 0);
	break;
      case 'l':  /* long */
	l_val = va_arg(args, long);
	hv_store(myhash, key, strlen(key), newSVnv(l_val), 0);

	/* Look up and store the tty if this is ttynum */
	if( !strcmp(key, "ttynum") ) store_ttydev( myhash, l_val );
	break;

      case 'J':  /* ignore; creates an undef value for this key in the hash */
	va_arg(args, long long);
	hv_store(myhash, key, strlen(key), &PL_sv_undef, 0);
	break;
      case 'j':  /* long long */
	ll_val = va_arg(args, long long);
	hv_store(myhash, key, strlen(key), newSVnv(ll_val), 0);
	break;

      case 'V':  /* perl scalar value */
	SV_val = va_arg(args, SV *);
	hv_store(myhash, key, strlen(key), SV_val, 0);
	break;

      default:
	croak("Unknown data format type `%c' returned from OS_get_table", *format);
	va_end(args); 
      }
    
    format++;
    fields++;
  }

  /* objectify the hash */
  ref = newRV_noinc((SV*) myhash);                        /* create ref from hash pointer */
  mystash = gv_stashpv("Proc::ProcessTable::Process", 1); /* create symbol table for this obj */
  blessed = sv_bless(ref, mystash);                       /* bless it */
  /* push it onto the array */
  av_push(Proclist, blessed);

  va_end(args);
}

/**********************************************************************/
/* Generic funcs generated by h2xs                                    */
/**********************************************************************/

static int
not_here(s)
char *s;
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(name, arg)
char *name;
int arg;
{
    errno = 0;
    switch (*name) {
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}

#ifdef	PROCESSTABLE_THREAD
pthread_mutex_t _mutex_table;
pthread_mutex_t _mutex_new;

void
mutex_op(int lock, pthread_mutex_t *mutex)
{
	if (lock == 0) {	/*unlock*/
		pthread_mutex_unlock(mutex);
	} else {		/*lock*/
		pthread_mutex_lock(mutex);
	}
}
#endif

void
mutex_new(int lock)
{
#ifdef	PROCESSTABLE_THREAD
	mutex_op(lock, &_mutex_new);
#endif
}

void
mutex_table(int lock)
{
#ifdef	PROCESSTABLE_THREAD
	mutex_op(lock, &_mutex_table);
#endif
}

MODULE = Proc::ProcessTable		PACKAGE = Proc::ProcessTable		
PROTOTYPES: DISABLE

BOOT:
#ifdef	PROCESSTABLE_THREAD
	pthread_mutex_init(&_mutex_table, NULL);
	pthread_mutex_init(&_mutex_new, NULL);
#endif

void
mutex_new(lock)
	int lock

void
mutex_table(lock)
	int lock

double
constant(name,arg)
	char *		name
	int		arg

SV*
table(obj)
     SV*  obj
     CODE:

     HV* hash;
     SV** fetched;


     mutex_table(1);
     /* Cache a pointer to the tty device hash */ 
     Ttydevs = perl_get_hv("Proc::ProcessTable::TTYDEVS", FALSE);

     /* dereference our object to a hash */
     hash = (HV*) SvRV(obj);

     /* If the Table array already exists on our object we clear it
        and store a pointer to it in Proclist */
     if( hv_exists(hash, "Table", 5) ){
       /* fetch the hash entry */
       fetched = hv_fetch(hash, "Table", 5, 0);
       /* what's stored in the hash is a ref to the array, so we need
          to dereference it */
       Proclist = (AV*) SvRV(*fetched);
       av_clear(Proclist);
     }
     else{
       /* Otherwise we create the array and store it on the object. */
       Proclist = newAV();
       hv_store(hash, "Table", 5, newRV_noinc((SV*)Proclist), 0);
     }

     /* Call get_table to build the process objects and push them onto
        the Proclist */
     OS_get_table();

     /* Return a ref to our process list */
     RETVAL = newRV_inc((SV*) Proclist);

     mutex_table(0);
     
     OUTPUT:
     RETVAL

void
fields(obj)
     SV*  obj
     PPCODE:

     int i;
     SV* my_sv;

     if( Fields == NULL ){
       PUSHMARK(SP);
       XPUSHs(obj);
       PUTBACK;
       perl_call_method("table", G_DISCARD);
     }

     EXTEND(SP,Numfields);
     for (i=0; i < Numfields; i++ ){
       my_sv = newSVpv(Fields[i],0);
       PUSHs(sv_2mortal(my_sv));
     }

void 
_initialize_os(obj)
     SV*  obj
     CODE:
     char* error;

     if( (error = OS_initialize()) != NULL ){
       croak(error);
     }
