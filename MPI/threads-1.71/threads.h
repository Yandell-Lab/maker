#ifndef _THREADS_H_
#define _THREADS_H_

/* Needed for 5.8.0 */
#ifndef CLONEf_JOIN_IN
#  define CLONEf_JOIN_IN        8
#endif
#ifndef SAVEBOOL
#  define SAVEBOOL(a)
#endif

/* Added in 5.11.x */
#ifndef G_WANT
#  define G_WANT                (128|1)
#endif

#endif
