/*
** 2001 September 15
**
** The author disclaims copyright to this source code.  In place of
** a legal notice, here is a blessing:
**
**    May you do good and not evil.
**    May you find forgiveness for yourself and forgive others.
**    May you share freely, never taking more than you give.
**
*************************************************************************
** This file contains code to implement a pseudo-random number
** generator (PRNG) for SQLite.
**
** Random numbers are used by some of the database backends in order
** to generate random integer keys for tables or random filenames.
**
** $Id: random.c,v 1.22 2004/02/14 19:14:50 matt Exp $
*/
#include "sqliteInt.h"
#include "os.h"


/*
** Get a single 8-bit random value from the RC4 PRNG.  The Mutex
** must be held while executing this routine.
**
** Why not just use a library random generator like lrand48() for this?
** Because the OP_NewRecno opcode in the VDBE depends on having a very
** good source of random numbers.  The lrand48() library function may
** well be good enough.  But maybe not.  Or maybe lrand48() has some
** subtle problems on some systems that could cause problems.  It is hard
** to know.  To minimize the risk of problems due to bad lrand48()
** implementations, SQLite uses this random number generator based
** on RC4, which we know works very well.
*/
static int randomByte(){
  int t;

  /* All threads share a single random number generator.
  ** This structure is the current state of the generator.
  */
  static struct {
    int isInit;          /* True if initialized */
    int i, j;            /* State variables */
    int s[256];          /* State variables */
  } prng;

  /* Initialize the state of the random number generator once,
  ** the first time this routine is called.  The seed value does
  ** not need to contain a lot of randomness since we are not
  ** trying to do secure encryption or anything like that...
  **
  ** Nothing in this file or anywhere else in SQLite does any kind of
  ** encryption.  The RC4 algorithm is being used as a PRNG (pseudo-random
  ** number generator) not as an encryption device.
  */
  if( !prng.isInit ){
    int i;
    char k[256];
    prng.j = 0;
    prng.i = 0;
    sqliteOsRandomSeed(k);
    for(i=0; i<256; i++){
      prng.s[i] = i;
    }
    for(i=0; i<256; i++){
      int t;
      prng.j = (prng.j + prng.s[i] + k[i]) & 0xff;
      t = prng.s[prng.j];
      prng.s[prng.j] = prng.s[i];
      prng.s[i] = t;
    }
    prng.isInit = 1;
  }

  /* Generate and return single random byte
  */
  prng.i = (prng.i + 1) & 0xff;
  prng.j = (prng.j + prng.s[prng.i]) & 0xff;
  t = prng.s[prng.i];
  prng.s[prng.i] = prng.s[prng.j];
  prng.s[prng.j] = t;
  t = prng.s[prng.i] + prng.s[prng.j];
  return prng.s[t & 0xff];
}

/*
** Return an random 8-bit integer.
*/
int sqliteRandomByte(){
  int r;
  sqliteOsEnterMutex();
  r = randomByte();
  sqliteOsLeaveMutex();
  return r;
}

/*
** Return a random 32-bit integer.  The integer is generated by making
** 4 calls to sqliteRandomByte().
*/
int sqliteRandomInteger(){
  int r;
  int i;
  sqliteOsEnterMutex();
  r = randomByte();
  for(i=1; i<4; i++){
    r = (r<<8) + randomByte();
  }
  sqliteOsLeaveMutex();
  return r;
}
