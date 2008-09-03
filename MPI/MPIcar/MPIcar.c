/*
 * This file was generated automatically by xsubpp version 1.9508 from the
 * contents of MPIcar.xs. Do not edit this file, edit MPIcar.xs instead.
 *
 *	ANY CHANGES MADE HERE WILL BE LOST!
 *
 */

#line 1 "MPIcar.xs"
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

#line 32 "MPIcar.c"
XS(XS_Parallel__MPIcar_constant); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_constant)
{
    dXSARGS;
    if (items != 2)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::constant(name, arg)");
    {
	char *	name = (char *)SvPV_nolen(ST(0));
	int	arg = (int)SvIV(ST(1));
	double	RETVAL;
	dXSTARG;

	RETVAL = constant(name, arg);
	XSprePUSH; PUSHn((double)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Comm_size); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Comm_size)
{
    dXSARGS;
    if (items != 1)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Comm_size(comm)");
    {
	MPI_Comm	comm;
#line 35 "MPIcar.xs"
        int d;
#line 61 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 37 "MPIcar.xs"
        MPIpm_errhandler("MPI_Comm_size",  MPI_Comm_size(comm,&d));
	RETVAL = d;
#line 74 "MPIcar.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Comm_rank); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Comm_rank)
{
    dXSARGS;
    if (items != 1)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Comm_rank(comm)");
    {
	MPI_Comm	comm;
#line 47 "MPIcar.xs"
        int d;
#line 90 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 49 "MPIcar.xs"
        MPIpm_errhandler("MPI_Comm_rank", MPI_Comm_rank(comm,&d));
	RETVAL = d;
#line 103 "MPIcar.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Init); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Init)
{
    dXSARGS;
    if (items != 0)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Init()");
    {
#line 57 "MPIcar.xs"
	AV *args_av;
	SV **sv_tmp;
	SV *sv_tmp2;
	SV *argv0;
        int argc;
	char **argv;
	int i;
#line 124 "MPIcar.c"
#line 65 "MPIcar.xs"
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
#line 221 "MPIcar.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_Parallel__MPIcar_MPI_Finalize); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Finalize)
{
    dXSARGS;
    if (items != 0)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Finalize()");
    {
#line 165 "MPIcar.xs"
	int rc;
#line 235 "MPIcar.c"
#line 167 "MPIcar.xs"
	rc = MPI_Finalize();
        MPIpm_errhandler("MPI_Finalize",rc);
#line 239 "MPIcar.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_Parallel__MPIcar_MPI_Send); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Send)
{
    dXSARGS;
    if (items != 6)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Send(ref, count, datatype, dest, tag, comm)");
    {
	SV*	ref = ST(0);
	int	count = (int)SvIV(ST(1));
	MPI_Datatype	datatype;
	int	dest = (int)SvIV(ST(3));
	int	tag = (int)SvIV(ST(4));
	MPI_Comm	comm;
#line 180 "MPIcar.xs"
    int len;
    void* buf;
    int ret;
#line 261 "MPIcar.c"

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    datatype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "datatype is not of type MPI_Datatype");

	if (sv_derived_from(ST(5), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(5)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 184 "MPIcar.xs"
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
#line 326 "MPIcar.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_Parallel__MPIcar_MPI_Recv); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Recv)
{
    dXSARGS;
    if (items != 6)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Recv(ref, count, datatype, source, tag, comm)");
    SP -= items;
    {
	SV *	ref = ST(0);
	int	count = (int)SvIV(ST(1));
	MPI_Datatype	datatype;
	int	source = (int)SvIV(ST(3));
	int	tag = (int)SvIV(ST(4));
	MPI_Comm	comm;
#line 244 "MPIcar.xs"
        void* buf;
        int ret;
	MPI_Status status;
#ifdef SEND_DEBUG
        int i;
#endif
#line 352 "MPIcar.c"

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    datatype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "datatype is not of type MPI_Datatype");

	if (sv_derived_from(ST(5), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(5)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 251 "MPIcar.xs"
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
#line 435 "MPIcar.c"
	PUTBACK;
	return;
    }
}

XS(XS_Parallel__MPIcar_MPI_Barrier); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Barrier)
{
    dXSARGS;
    if (items != 1)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Barrier(comm)");
    {
	MPI_Comm	comm;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 324 "MPIcar.xs"
        MPIpm_errhandler("MPI_Barrier",MPI_Barrier(comm));
#line 460 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Bcast); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Bcast)
{
    dXSARGS;
    if (items != 5)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Bcast(ref, count, datatype, root, comm)");
    {
	SV *	ref = ST(0);
	int	count = (int)SvIV(ST(1));
	MPI_Datatype	datatype;
	int	root = (int)SvIV(ST(3));
	MPI_Comm	comm;
#line 335 "MPIcar.xs"
        void* buf;
        int ret;
#line 480 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    datatype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "datatype is not of type MPI_Datatype");

	if (sv_derived_from(ST(4), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(4)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 338 "MPIcar.xs"
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
#line 528 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Wtime); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Wtime)
{
    dXSARGS;
    if (items != 0)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Wtime()");
    {
	double	RETVAL;
	dXSTARG;

	RETVAL = MPI_Wtime();
	XSprePUSH; PUSHn((double)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Wtick); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Wtick)
{
    dXSARGS;
    if (items != 0)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Wtick()");
    {
	double	RETVAL;
	dXSTARG;

	RETVAL = MPI_Wtick();
	XSprePUSH; PUSHn((double)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Initialized); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Initialized)
{
    dXSARGS;
    if (items != 0)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Initialized()");
    {
#line 381 "MPIcar.xs"
        int flag;
#line 574 "MPIcar.c"
	int	RETVAL;
	dXSTARG;
#line 383 "MPIcar.xs"
        MPIpm_errhandler("MPI_Initialized",  MPI_Initialized(&flag));
	RETVAL = flag;
#line 580 "MPIcar.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Abort); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Abort)
{
    dXSARGS;
    if (items != 2)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Abort(comm, errorcode)");
    {
	MPI_Comm	comm;
	int	errorcode = (int)SvIV(ST(1));

	if (sv_derived_from(ST(0), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 394 "MPIcar.xs"
	MPIpm_errhandler("MPI_Abort", MPI_Abort(comm,errorcode));
#line 604 "MPIcar.c"
    }
    XSRETURN_EMPTY;
}

XS(XS_Parallel__MPIcar_MPI_Reduce); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Reduce)
{
    dXSARGS;
    if (items != 7)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Reduce(sendref, recvref, count, datatype, op, root, comm)");
    {
	SV *	sendref = ST(0);
	SV *	recvref = ST(1);
	int	count = (int)SvIV(ST(2));
	MPI_Datatype	datatype;
	MPI_Op	op;
	int	root = (int)SvIV(ST(5));
	MPI_Comm	comm;
#line 407 "MPIcar.xs"
        void* sendbuf, *recvbuf;
        int ret;
#line 626 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(3), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(3)));
	    datatype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "datatype is not of type MPI_Datatype");

	if (sv_derived_from(ST(4), "MPI_Op")) {
	    IV tmp = SvIV((SV*)SvRV(ST(4)));
	    op = INT2PTR(MPI_Op,tmp);
	}
	else
	    Perl_croak(aTHX_ "op is not of type MPI_Op");

	if (sv_derived_from(ST(6), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(6)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 410 "MPIcar.xs"
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
#line 669 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Allreduce); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Allreduce)
{
    dXSARGS;
    if (items != 6)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Allreduce(sendref, recvref, count, datatype, op, comm)");
    {
	SV *	sendref = ST(0);
	SV *	recvref = ST(1);
	int	count = (int)SvIV(ST(2));
	MPI_Datatype	datatype;
	MPI_Op	op;
	MPI_Comm	comm;
#line 439 "MPIcar.xs"
        void* sendbuf, *recvbuf;
        int ret;
#line 690 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(3), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(3)));
	    datatype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "datatype is not of type MPI_Datatype");

	if (sv_derived_from(ST(4), "MPI_Op")) {
	    IV tmp = SvIV((SV*)SvRV(ST(4)));
	    op = INT2PTR(MPI_Op,tmp);
	}
	else
	    Perl_croak(aTHX_ "op is not of type MPI_Op");

	if (sv_derived_from(ST(5), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(5)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 442 "MPIcar.xs"
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
#line 733 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Scatter); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Scatter)
{
    dXSARGS;
    if (items != 8)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Scatter(sendref, sendcnt, sendtype, recvref, recvcnt, recvtype, root, comm)");
    {
	SV *	sendref = ST(0);
	int	sendcnt = (int)SvIV(ST(1));
	MPI_Datatype	sendtype;
	SV *	recvref = ST(3);
	int	recvcnt = (int)SvIV(ST(4));
	MPI_Datatype	recvtype;
	int	root = (int)SvIV(ST(6));
	MPI_Comm	comm;
#line 473 "MPIcar.xs"
        void* sendbuf, *recvbuf;
        int ret;
#line 756 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    sendtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "sendtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(5), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(5)));
	    recvtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "recvtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(7), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(7)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 476 "MPIcar.xs"
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
#line 813 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Gather); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Gather)
{
    dXSARGS;
    if (items != 8)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Gather(sendref, sendcnt, sendtype, recvref, recvcnt, recvtype, root, comm)");
    {
	SV *	sendref = ST(0);
	int	sendcnt = (int)SvIV(ST(1));
	MPI_Datatype	sendtype;
	SV *	recvref = ST(3);
	int	recvcnt = (int)SvIV(ST(4));
	MPI_Datatype	recvtype;
	int	root = (int)SvIV(ST(6));
	MPI_Comm	comm;
#line 521 "MPIcar.xs"
        void* sendbuf, *recvbuf;
        int ret;
#line 836 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    sendtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "sendtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(5), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(5)));
	    recvtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "recvtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(7), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(7)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 524 "MPIcar.xs"
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
#line 880 "MPIcar.c"
    }
    XSRETURN(1);
}

XS(XS_Parallel__MPIcar_MPI_Sendrecv); /* prototype to pass -Wmissing-prototypes */
XS(XS_Parallel__MPIcar_MPI_Sendrecv)
{
    dXSARGS;
    if (items != 11)
	Perl_croak(aTHX_ "Usage: Parallel::MPIcar::MPI_Sendrecv(sendref, sendcount, sendtype, dest, sendtag, recvref, recvcount, recvtype, source, recvtag, comm)");
    SP -= items;
    {
	SV *	sendref = ST(0);
	int	sendcount = (int)SvIV(ST(1));
	MPI_Datatype	sendtype;
	int	dest = (int)SvIV(ST(3));
	int	sendtag = (int)SvIV(ST(4));
	SV *	recvref = ST(5);
	int	recvcount = (int)SvIV(ST(6));
	MPI_Datatype	recvtype;
	int	source = (int)SvIV(ST(8));
	int	recvtag = (int)SvIV(ST(9));
	MPI_Comm	comm;
#line 558 "MPIcar.xs"
        void* sendbuf, *recvbuf;
        int ret;
	MPI_Status status;
#line 908 "MPIcar.c"
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(2), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(2)));
	    sendtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "sendtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(7), "MPI_Datatype")) {
	    IV tmp = SvIV((SV*)SvRV(ST(7)));
	    recvtype = INT2PTR(MPI_Datatype,tmp);
	}
	else
	    Perl_croak(aTHX_ "recvtype is not of type MPI_Datatype");

	if (sv_derived_from(ST(10), "MPI_Comm")) {
	    IV tmp = SvIV((SV*)SvRV(ST(10)));
	    comm = INT2PTR(MPI_Comm,tmp);
	}
	else
	    Perl_croak(aTHX_ "comm is not of type MPI_Comm");
#line 562 "MPIcar.xs"
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
#line 971 "MPIcar.c"
	PUTBACK;
	return;
    }
}

#ifdef __cplusplus
extern "C"
#endif
XS(boot_Parallel__MPIcar); /* prototype to pass -Wmissing-prototypes */
XS(boot_Parallel__MPIcar)
{
    dXSARGS;
    char* file = __FILE__;

    XS_VERSION_BOOTCHECK ;

        newXS("Parallel::MPIcar::constant", XS_Parallel__MPIcar_constant, file);
        newXS("Parallel::MPIcar::MPI_Comm_size", XS_Parallel__MPIcar_MPI_Comm_size, file);
        newXS("Parallel::MPIcar::MPI_Comm_rank", XS_Parallel__MPIcar_MPI_Comm_rank, file);
        newXS("Parallel::MPIcar::MPI_Init", XS_Parallel__MPIcar_MPI_Init, file);
        newXS("Parallel::MPIcar::MPI_Finalize", XS_Parallel__MPIcar_MPI_Finalize, file);
        newXS("Parallel::MPIcar::MPI_Send", XS_Parallel__MPIcar_MPI_Send, file);
        newXS("Parallel::MPIcar::MPI_Recv", XS_Parallel__MPIcar_MPI_Recv, file);
        newXS("Parallel::MPIcar::MPI_Barrier", XS_Parallel__MPIcar_MPI_Barrier, file);
        newXS("Parallel::MPIcar::MPI_Bcast", XS_Parallel__MPIcar_MPI_Bcast, file);
        newXS("Parallel::MPIcar::MPI_Wtime", XS_Parallel__MPIcar_MPI_Wtime, file);
        newXS("Parallel::MPIcar::MPI_Wtick", XS_Parallel__MPIcar_MPI_Wtick, file);
        newXS("Parallel::MPIcar::MPI_Initialized", XS_Parallel__MPIcar_MPI_Initialized, file);
        newXS("Parallel::MPIcar::MPI_Abort", XS_Parallel__MPIcar_MPI_Abort, file);
        newXS("Parallel::MPIcar::MPI_Reduce", XS_Parallel__MPIcar_MPI_Reduce, file);
        newXS("Parallel::MPIcar::MPI_Allreduce", XS_Parallel__MPIcar_MPI_Allreduce, file);
        newXS("Parallel::MPIcar::MPI_Scatter", XS_Parallel__MPIcar_MPI_Scatter, file);
        newXS("Parallel::MPIcar::MPI_Gather", XS_Parallel__MPIcar_MPI_Gather, file);
        newXS("Parallel::MPIcar::MPI_Sendrecv", XS_Parallel__MPIcar_MPI_Sendrecv, file);
    XSRETURN_YES;
}

