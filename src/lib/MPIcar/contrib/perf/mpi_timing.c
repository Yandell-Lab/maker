/******************************************************************************
* MPI Timing Program - C Version
* FILE: mpi_timing.c
* OTHER FILES:  make.mpi_timing.c
* DESCRIPTION:  MPI timing example code.  C version.
*   In this example code, a MPI communication timing test is performed.
*   The processor with taskid = 0 will send "reps" number of messages to
*   the processor with taskid = 1, waiting for a reply between each rep.
*   Before and after timings are made for each rep and an average 
*   calculated when completed.
* AUTHOR: Blaise Barney - adapted from pvm C version
* CONVERTED TO MPI: George L. Gusciora (1/25/95)
* Slight modifications by: Josh Wilmes
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#define NUMBER_REPS     200
#define MESSAGE_SIZE    128

int main(argc,argv)
int argc;
char *argv[];
{
   int reps;                            /* number of samples per test */
   int dt1, dt2;                        /* time for one iter */
   int at1, at2;                        /* accum. time */
   int n;
   char *inmsg, *outmsg;                /* Buffer containing message */
   int type;
   int numtasks, taskid;
   int rc, dest, source, nbytes;
   double starttime, endtime;
   MPI_Status status;

   rc = MPI_Init(&argc,&argv);
   rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   if (rc != 0)
      printf ("error initializing MPI and obtaining task ID information\n");
   else
      printf ("mpi_pi_mm MPI task ID = %d\n", taskid);
   if (numtasks != 2)
   {
      fprintf(stderr, "Error; Set environment variable MP_PROCS to 2\n");
      exit(1);
   }

   at1 = 0;
   inmsg = (char *) malloc(MESSAGE_SIZE);
   outmsg = (char *) malloc(MESSAGE_SIZE);
   type = 1;
   reps = NUMBER_REPS;

   if (taskid == 0)
   {
      /* round-trip timing test */
      printf("Doing round trip test, minimal message size, %d reps.\n",reps);
      dest = 1;
      source = 1;
      for (n = 1; n <= reps; n++)
      {
         MPI_Barrier(MPI_COMM_WORLD);
	 starttime = MPI_Wtime();          /* before time */
         /* send message to worker - message type set to 1.  If  */
         /* return code is less than zero quit */
         rc = MPI_Send(outmsg, MESSAGE_SIZE, MPI_CHAR, dest, type, MPI_COMM_WORLD);
         if (rc < 0)
         {
            fprintf(stderr, "Send error in processor 1\n");
            exit(1);
         }
         /* Now wait to receive the echo reply from the worker  */
         /* Quit if return code */
         /* is less than zero */
         rc = MPI_Recv(inmsg, MESSAGE_SIZE, MPI_CHAR, source, type, MPI_COMM_WORLD, &status);
         if (rc < 0)
         {
            fprintf(stderr, "Recieve error in processor 0\n");
            exit(1);
         }

         endtime = MPI_Wtime();  /* after time */

         /* calculate round trip time and print */
	 dt1 = (endtime - starttime) * 1000000;
         printf("round trip# %2d   uSec = %8d\n", n, dt1);

         at1 += dt1;
      }
      printf("\n*** Round Trip Avg uSec = %d\n", at1 / reps);
   } else if (taskid == 1)
   {
      dest = 0;
      source = 0;
      for (n = 1; n <= reps; n++)
      {
         MPI_Barrier(MPI_COMM_WORLD);
         rc = MPI_Recv(inmsg, MESSAGE_SIZE, MPI_CHAR, source, type, MPI_COMM_WORLD, &status);
         if (rc < 0)
         {
            fprintf(stderr, "Recieve error in processor 0\n");
            exit(1);
         }
         rc = MPI_Send(outmsg, MESSAGE_SIZE, MPI_CHAR, dest, type, MPI_COMM_WORLD);
         if (rc < 0)
         {
            fprintf(stderr, "Send error in processor 1\n");
            exit(1);
         }
      }
   }
   MPI_Finalize();
   exit(0);
}
