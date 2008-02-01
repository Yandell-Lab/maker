#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include <mpi.h>

#define MPI_STRING ((MPI_Datatype)34)
/*#define  ARGV_DEBUG*/

#ifdef FLOAT_HACK
# undef MPI_FLOAT
# define MPI_FLOAT MPI_DOUBLE
#endif

#include "utils.c"

/* don't ask- mpicc hoses up these definitions */
#undef VERSION
#undef XS_VERSION
#define VERSION "0.03"
#define XS_VERSION "0.03"

MODULE = Parallel::MPIcar		PACKAGE = Parallel::MPIcar
PROTOTYPES: DISABLE

double
constant(name,arg)
	char *		name
	int		arg


int
MPI_Comm_size(comm)
	MPI_Comm	comm
      PREINIT:
        int d;
      CODE:
        MPIpm_errhandler("MPI_Comm_size",  MPI_Comm_size(comm,&d));
	RETVAL = d;
      OUTPUT:
        RETVAL


int
MPI_Comm_rank(comm)
	MPI_Comm	comm
      PREINIT:
        int d;
      CODE:
        MPIpm_errhandler("MPI_Comm_rank", MPI_Comm_rank(comm,&d));
	RETVAL = d;
      OUTPUT:
        RETVAL

void
MPI_Init()
      PREINIT:
	AV *args_av;
	SV **sv_tmp;
	SV *sv_tmp2;
	SV *argv0;
        int argc;
	char **argv;
	int i;
      CODE:
        /* Get @ARGV */
        args_av = perl_get_av("main::ARGV", TRUE);
	if(args_av == NULL)
	    croak("Parallel::MPIcar: $ARGV not set. Oops");

        /* Get $0 */
        argv0 = perl_get_sv("main::0", FALSE);
	if(argv0 == NULL)
	    croak("Parallel::MPIcar: $0 not set. Oops");

        /* We can't run MPI from a -e or a - script. */
	if(!strncmp(SvPV(argv0, PL_na), "-", 1) ||
           !strncmp(SvPV(argv0, PL_na), "-e", 2))
	    croak("Parallel::MPIcar: Cannot use MPI with command line script");
 
	/* debug */
#ifdef ARGV_DEBUG
	printf("[%d] av_len=%d\n",getpid(),av_len(args_av));
 	for (i=0 ; i <= av_len(args_av) ; i++) {
	    sv_tmp = av_fetch(args_av,i,0);
	    printf("[%d] $ARGV[%d]=%s\n",getpid(),i,SvPV(*sv_tmp, PL_na));
	}
#endif
 
        /* argc = $#ARGV+1  +1 for argv[0] */
	argc = av_len(args_av)+2;
	if (argc == 1) {
	    croak("MPI_Init: no arguments found in the argv");
	} else {
            /* build up argv by setting argv[0] from $0
	       and the rest from @ARGV

               We add an extra NULL to prevent a coredump
               during global destruction
	    */
	    argv = (char**) malloc((argc+1) * sizeof(char*));
            argv[0] = strdup(SvPV(argv0, PL_na));
	    for(i=1; i<argc; i++) {
	        sv_tmp = av_fetch(args_av,i-1,0);
	        if (sv_tmp == NULL) {
		    argv[i] = NULL;
	        } else {
		    argv[i] = strdup(SvPV(*sv_tmp, PL_na));
	        }
            }
            argv[argc] = NULL; /* prevents coredumps */
	}

	/* debug */
#ifdef ARGV_DEBUG
	printf("[%d] argc=%d\n",getpid(),argc);
 	for (i=0;i<argc;i++)
	    printf("[%d] argv[%d]=%s\n",getpid(),i,argv[i]);
#endif
        
        /* Call the actual function */
	MPIpm_errhandler("MPI_Init",MPI_Init(&argc, &argv));

        /* Allow later MPI funcs to return to our error handler */
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	/* debug */
#ifdef ARGV_DEBUG
	printf("[%d] argc=%d\n",getpid(),argc);
 	for (i=0;i<argc;i++)
	    printf("[%d] argv[%d]=%s\n",getpid(),i,argv[i]);
#endif

        /* Now copy argv back to @ARGV by clearing out @ARGV and pushing
	   each arg back onto @ARGV. */
        if(argc > 1) {
	    av_extend(args_av, argc-1);
            av_clear(args_av);
	    for(i=1;i<argc;i++) {
                sv_tmp2 = newSVpv(argv[i], 0);
                sv_tmp2 = SvREFCNT_inc(sv_tmp2);
		av_push(args_av, sv_tmp2);
	    }
	} else {
            /* No args; clear @ARGV */
            av_clear(args_av);
	}

	/* debug */
#ifdef ARGV_DEBUG
	printf("[%d] av_len=%d\n",getpid(),av_len(args_av));
 	for (i=0 ; i <= av_len(args_av) ; i++) {
	    sv_tmp = av_fetch(args_av,i,0);
	    if (sv_tmp == NULL)
		printf("[%d] $ARGV[%d]=undef\n",getpid(),i);
	    else {
		printf("[%d] $ARGV[%d]=%s\n",getpid(),i,SvPV(*sv_tmp, PL_na));
	    }
	}
#endif


void
MPI_Finalize()
    PREINIT:
	int rc;
    CODE:
	rc = MPI_Finalize();
        MPIpm_errhandler("MPI_Finalize",rc);


void
MPI_Send(ref, count, datatype, dest, tag, comm)
    SV* ref
    int count
    MPI_Datatype datatype
    int	dest
    int	tag
    MPI_Comm comm
  PREINIT:
    int len;
    void* buf;
    int ret;
  CODE:     
    if (! SvROK(ref)) 
	croak("MPI_Send: First argument is not a reference!");

    if (SvTYPE(SvRV(ref)) == SVt_PVHV) {
	croak("MPI_Send: Hashes are not supported yet");
    } else if (SvTYPE(SvRV(ref)) == SVt_PVAV) {
#ifdef SEND_DEBUG
	int i;
#endif /* SEND_DEBUG */
	AV *stuff = (AV*) SvRV(ref);
        if(count > (av_len(stuff)+1)) {
            printf("MPI_Send: count param is larger than given array.  Using "
                 "array length.\n");
            count = av_len(stuff)+1;
        }
#ifdef SEND_DEBUG
	printf("[%d] av_len=%d\n",getpid(),av_len(stuff));
 	for (i=0 ; i <= av_len(stuff) ; i++) {
	    SV **sv_tmp = av_fetch(stuff,i,0);
	    if (sv_tmp == NULL)
		printf("[%d] $stuff[%d]=undef\n",getpid(),i);
	    else {
		printf("[%d] $stuff[%d]=%s\n",getpid(),i,SvPV(*sv_tmp, PL_na));
	    }
	}
#endif /* SEND_DEBUG */
        len = MPIpm_packarray(&buf, stuff, datatype, count);
#ifdef SEND_DEBUG
        printf("[%d] len=%d\n[%d] ", getpid(), len, getpid());
	for(i=0;i<len;i++) {
	    printf("%02x ", (unsigned char)((char*)buf)[i]);
	    if((i!=0) && (i%16) == 0) printf("\n[%d] ", getpid());
	}
	printf("\n");
#endif /* SEND_DEBUG */
        if(datatype == MPI_STRING)
	    MPIpm_errhandler("MPI_Send",
		             MPI_Send(&len, 1, MPI_INT, dest, tag, comm));
	MPIpm_errhandler("MPI_Send",
	                 MPI_Send(buf, len, MPI_CHAR, dest, tag, comm));
    } else {
	buf = (void*) malloc(MPIpm_bufsize(datatype,NULL,count));
	MPIpm_packscalar(buf,SvRV(ref),datatype);

	ret = MPI_Send(buf,count,datatype,dest,tag,comm);

	free(buf);
	MPIpm_errhandler("MPI_Send",ret);
    }


void
MPI_Recv(ref, count, datatype, source, tag, comm)
	SV *	ref
	int     count
	MPI_Datatype	datatype
	int	source
	int	tag
	MPI_Comm comm
      PREINIT:
        void* buf;
        int ret;
	MPI_Status status;
#ifdef SEND_DEBUG
        int i;
#endif
      PPCODE:
	if (! SvROK(ref)) 
            croak("MPI_Recv: First argument is not a reference!");

	if (SvTYPE(SvRV(ref)) == SVt_PVHV) {
            croak("MPI_Recv: Hashes are not supported yet.");
	} else if (SvTYPE(SvRV(ref)) == SVt_PVAV) {
	    int len;
	    AV *stuff = (AV*) SvRV(ref);
            switch(datatype) {
	      case MPI_STRING:
	        MPI_Recv(&len, 1, MPI_INT, source, tag, comm, &status);
		break;
	      case MPI_INT:
                len = count * sizeof(int);
                break;
#ifndef FLOAT_HACK
	      case MPI_FLOAT:
                len = count * sizeof(float);
                break;
#endif
	      case MPI_DOUBLE:
                len = count * sizeof(double);
                break;
            }
#ifdef SEND_DEBUG
	    printf("[%d] len=%d\n", getpid(), len);
#endif
            buf = (char *) malloc(len);
            ret = MPI_Recv(buf, len, MPI_CHAR, source, tag, comm, &status);
#ifdef SEND_DEBUG
	    printf("[%d] len=%d\n[%d] ", getpid(), len, getpid());
	    for(i=0;i<len;i++) {
		printf("%02x ", (unsigned char)((char*)buf)[i]);
		if((i!=0) && (i%16) == 0) printf("\n[%d] ", getpid());
	    }
	    printf("\n");
#endif /* SEND_DEBUG */
            MPIpm_unpackarray(buf, &stuff, datatype, count);
	    MPIpm_errhandler("MPI_Recv",ret);
	    /* return the status as a 4 element array:
	     * (count,MPI_SOURCE,MPI_TAG,MPI_ERROR) */
	    XPUSHs(sv_2mortal(newSViv(status.count)));
	    XPUSHs(sv_2mortal(newSViv(status.MPI_SOURCE)));
	    XPUSHs(sv_2mortal(newSViv(status.MPI_TAG)));
	    XPUSHs(sv_2mortal(newSViv(status.MPI_ERROR)));
	} else {
	  buf = (void*) malloc(MPIpm_bufsize(datatype,NULL,count));

	  ret = MPI_Recv(buf,count,datatype,source,tag,comm,&status);

	  if (datatype == MPI_CHAR){
	      sv_setpvn(SvRV(ref), buf, count);
	  }	  
	  else{
	      MPIpm_unpackscalar(buf,SvRV(ref),datatype);
	  }

	  free(buf);
	  MPIpm_errhandler("MPI_Recv",ret);

	  /* return the status as a 4 element array:
	   * (count,MPI_SOURCE,MPI_TAG,MPI_ERROR) */
	  XPUSHs(sv_2mortal(newSViv(status.count)));
	  XPUSHs(sv_2mortal(newSViv(status.MPI_SOURCE)));
	  XPUSHs(sv_2mortal(newSViv(status.MPI_TAG)));
	  XPUSHs(sv_2mortal(newSViv(status.MPI_ERROR)));
	}


int
MPI_Barrier(comm)
	MPI_Comm	comm
     CODE:
        MPIpm_errhandler("MPI_Barrier",MPI_Barrier(comm));


int
MPI_Bcast(ref, count, datatype, root, comm)
	SV *    ref
	int     count
        MPI_Datatype 	datatype
	int	root
	MPI_Comm	comm
      PREINIT:
        void* buf;
        int ret;
      CODE:     
	if (! SvROK(ref)) 
            croak("MPI_Bcast: First argument is not a reference!");

	if (SvTYPE(SvRV(ref)) == SVt_PVAV) {
	    AV* array = (AV*) SvRV(ref);
            int len;
            int rank;
	    MPI_Comm_rank(comm, &rank);
	    if(rank == root)
		MPIpm_packarray(&buf, array, datatype, count);
	    else
		buf = (void*) malloc(MPIpm_bufsize(datatype,SvRV(ref),count));

	    ret = MPI_Bcast(buf,count,datatype,root,comm);
	    if(rank != root)
		MPIpm_unpackarray(buf,&array,datatype, count);
	    MPIpm_errhandler("MPI_Bcast",ret);
#if 0
	    croak("MPI_Bcast: Arrays are not implemented yet.\n");
#endif
	} else {
	    buf = (void*) malloc(MPIpm_bufsize(datatype,SvRV(ref),count));
	    MPIpm_packscalar(buf,SvRV(ref),datatype);

	    ret = MPI_Bcast(buf,count,datatype,root,comm);

	    MPIpm_unpackscalar(buf,SvRV(ref),datatype);
	    free(buf);
	    MPIpm_errhandler("MPI_Bcast",ret);
	}


double
MPI_Wtime()


double
MPI_Wtick()


int
MPI_Initialized()
      PREINIT:
        int flag;
      CODE:
        MPIpm_errhandler("MPI_Initialized",  MPI_Initialized(&flag));
	RETVAL = flag;
      OUTPUT:
        RETVAL


void
MPI_Abort(comm, errorcode)
	MPI_Comm	comm
	int	errorcode
      CODE:
	MPIpm_errhandler("MPI_Abort", MPI_Abort(comm,errorcode));


int
MPI_Reduce(sendref, recvref, count, datatype, op, root, comm)
	SV *	sendref
	SV *	recvref
	int	count
	MPI_Datatype	datatype
	MPI_Op	op
	int	root
	MPI_Comm	comm
      PREINIT:
        void* sendbuf, *recvbuf;
        int ret;
      CODE:     
	if (! SvROK(sendref) || ! SvROK(recvref))
            croak("MPI_Reduce: First two arguments must be references!");

	if (SvTYPE(SvRV(sendref)) == SVt_PVAV ||
            SvTYPE(SvRV(recvref)) == SVt_PVAV) {
	  croak("MPI_Reduce: Arrays are not yet implemented");
	} else {
	  sendbuf = (void*)malloc(MPIpm_bufsize(datatype,SvRV(sendref),count));
	  recvbuf = (void*)malloc(MPIpm_bufsize(datatype,SvRV(recvref),count));
	  MPIpm_packscalar(sendbuf,SvRV(sendref),datatype);

	  ret = MPI_Reduce(sendbuf,recvbuf,count,datatype,op,root,comm);

	  MPIpm_unpackscalar(recvbuf,SvRV(recvref),datatype);
	  free(sendbuf);
	  free(recvbuf);
	  MPIpm_errhandler("MPI_Reduce",ret);
	}


int
MPI_Allreduce(sendref, recvref, count, datatype, op, comm)
	SV *	sendref
	SV *	recvref
	int	count
	MPI_Datatype	datatype
	MPI_Op	op
	MPI_Comm	comm
      PREINIT:
        void* sendbuf, *recvbuf;
        int ret;
      CODE:     
	if (! SvROK(sendref) || ! SvROK(recvref))
            croak("MPI_Allreduce: First two arguments must be references!");

	if (SvTYPE(SvRV(sendref)) == SVt_PVAV ||
            SvTYPE(SvRV(recvref)) == SVt_PVAV) {
	  croak("MPI_Allreduce: Arrays are not yet implemented");
	} else {
	  sendbuf = (void*)malloc(MPIpm_bufsize(datatype,SvRV(sendref),count));
	  recvbuf = (void*)malloc(MPIpm_bufsize(datatype,SvRV(recvref),count));
	  MPIpm_packscalar(sendbuf,SvRV(sendref),datatype);

	  ret = MPI_Allreduce(sendbuf,recvbuf,count,datatype,op,comm);

	  MPIpm_unpackscalar(recvbuf,SvRV(recvref),datatype);
	  free(sendbuf);
	  free(recvbuf);
	  MPIpm_errhandler("MPI_Allreduce",ret);
	}


int
MPI_Scatter(sendref, sendcnt, sendtype, recvref, recvcnt, recvtype, root, comm)
	SV *    sendref
        int     sendcnt
        MPI_Datatype 	sendtype
	SV *    recvref
        int     recvcnt
        MPI_Datatype 	recvtype
	int	root
	MPI_Comm	comm
      PREINIT:
        void* sendbuf, *recvbuf;
        int ret;
      CODE:     
	if (! SvROK(sendref) || ! SvROK(recvref))
            croak("MPI_Scatter: First and Fourth arguments must be references!");

	if (SvTYPE(SvRV(sendref)) == SVt_PVAV)
        {
            int rank,len;
            MPI_Comm_rank(comm, &rank);
            if(rank == root)
		len = MPIpm_packarray(&sendbuf, (AV*)SvRV(sendref), sendtype,0);
	    recvbuf = (void*) malloc(MPIpm_bufsize(recvtype,(SV*)SvRV(sendref),recvcnt));
	    ret = MPI_Scatter(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm); 
	    if(SvTYPE(SvRV(recvref)) == SVt_PVAV) {
		AV *recv = (AV*) SvRV(recvref);
		MPIpm_unpackarray(recvbuf,&recv,recvtype,recvcnt);
            } else {
		SV *recv = (SV*) SvRV(recvref);
		MPIpm_unpackscalar(recvbuf,recv,recvtype);
            }
	    MPIpm_errhandler("MPI_Scatter",ret);
	} else {
	  sendbuf = (void*)calloc(MPIpm_bufsize(sendtype,SvRV(sendref),sendcnt)+1,1);
	  recvbuf = (void*)calloc(MPIpm_bufsize(recvtype,SvRV(recvref),recvcnt)+1,1);
	  MPIpm_packscalar(sendbuf,SvRV(sendref),sendtype);

	  ret = MPI_Scatter(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm); 
	  MPIpm_unpackscalar(recvbuf,SvRV(recvref),recvtype);
#ifdef WHY_DOES_THIS_MAKE_IT_SEGFAULT
	  free(sendbuf);
	  free(recvbuf);
#endif
	  MPIpm_errhandler("MPI_Scatter",ret);
	}


int
MPI_Gather(sendref, sendcnt, sendtype, recvref, recvcnt, recvtype, root, comm)
	SV *    sendref
        int     sendcnt
        MPI_Datatype 	sendtype
	SV *    recvref
        int     recvcnt
        MPI_Datatype 	recvtype
	int	root
	MPI_Comm	comm
      PREINIT:
        void* sendbuf, *recvbuf;
        int ret;
      CODE:     
	if (! SvROK(sendref) || ! SvROK(recvref))
            croak("MPI_Gather: First and Fourth arguments must be references!");

	if (SvTYPE(SvRV(sendref)) == SVt_PVAV ||
            SvTYPE(SvRV(recvref)) == SVt_PVAV)
	{
	    croak("MPI_Gather: Arrays are not implemented yet.");
	} else {
	  sendbuf = (void*)malloc(MPIpm_bufsize(sendtype,SvRV(sendref),sendcnt));
	  recvbuf = (void*)malloc(MPIpm_bufsize(recvtype,SvRV(recvref),recvcnt));
	  MPIpm_packscalar(sendbuf,SvRV(sendref),sendtype);

	  ret = MPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,recvcnt,recvtype,root,comm);

	  MPIpm_unpackscalar(recvbuf,SvRV(recvref),recvtype);
	  free(sendbuf);
	  free(recvbuf);
	  MPIpm_errhandler("MPI_Gather",ret);
	}

int
MPI_Sendrecv(sendref, sendcount, sendtype, dest, sendtag, recvref, recvcount, recvtype, source, recvtag, comm)
	SV *	sendref
	int	sendcount
	MPI_Datatype	sendtype
	int	dest
	int	sendtag
	SV *	recvref
	int	recvcount
	MPI_Datatype	recvtype
	int	source
	int	recvtag
	MPI_Comm	comm
      PREINIT:
        void* sendbuf, *recvbuf;
        int ret;
	MPI_Status status;
      PPCODE:     
	if (! SvROK(sendref) || ! SvROK(recvref))
            croak("MPI_Sendrecv: First and Fourth arguments must be references!");

	if (SvTYPE(SvRV(sendref)) == SVt_PVAV &&
            SvTYPE(SvRV(recvref)) == SVt_PVAV)
	{
            AV* array = (AV*) SvRV(recvref);
	    int len;
	    recvbuf = malloc(MPIpm_bufsize(recvtype, SvRV(recvref), recvcount));
	    len = MPIpm_packarray(&sendbuf, (AV*)SvRV(sendref), sendtype, sendcount);
	    ret = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                               recvbuf, recvcount, recvtype, source, recvtag,
                               comm, &status);

	    MPIpm_unpackarray(recvbuf,&array,recvtype,recvcount);
	    MPIpm_errhandler("MPI_Sendrecv",ret);
	} else {
	  sendbuf = (void*)malloc(MPIpm_bufsize(sendtype,SvRV(sendref),sendcount));
	  recvbuf = (void*)malloc(MPIpm_bufsize(recvtype,SvRV(recvref),recvcount));
	  MPIpm_packscalar(sendbuf,SvRV(sendref),sendtype);
          MPIpm_packscalar(recvbuf,SvRV(recvref),recvtype);
	  
	  ret = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest,
                             sendtag, recvbuf, recvcount, recvtype,
			     source, recvtag, comm, &status);

	  MPIpm_unpackscalar(recvbuf,SvRV(recvref),recvtype);
	  free(sendbuf);
	  free(recvbuf);
	  MPIpm_errhandler("MPI_Sendrecv",ret);
	}

	/* return the status as a 4 element array:
	 * (count,MPI_SOURCE,MPI_TAG,MPI_ERROR) */
	XPUSHs(sv_2mortal(newSViv(status.count)));
	XPUSHs(sv_2mortal(newSViv(status.MPI_SOURCE)));
	XPUSHs(sv_2mortal(newSViv(status.MPI_TAG)));
	XPUSHs(sv_2mortal(newSViv(status.MPI_ERROR)));
