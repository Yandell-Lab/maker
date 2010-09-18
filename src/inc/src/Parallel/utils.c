/* This is included into MPI.xs.  It's just here to unclutter MPI.xs a bit */

static double
constant(name, arg)
char *name;
int arg;
{
  errno = 0;
  if (strEQ(name, "MPI_2COMPLEX"))
#ifdef MPI_2COMPLEX
    return MPI_2COMPLEX;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_2DOUBLE_COMPLEX"))
#ifdef MPI_2DOUBLE_COMPLEX
    return MPI_2DOUBLE_COMPLEX;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_2DOUBLE_PRECISION"))
#ifdef MPI_2DOUBLE_PRECISION
    return MPI_2DOUBLE_PRECISION;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_2INT"))
#ifdef MPI_2INT
    return MPI_2INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_2INTEGER"))
#ifdef MPI_2INTEGER
    return MPI_2INTEGER;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_2REAL"))
#ifdef MPI_2REAL
    return MPI_2REAL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_ANY_SOURCE"))
#ifdef MPI_ANY_SOURCE
    return MPI_ANY_SOURCE;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_ANY_TAG"))
#ifdef MPI_ANY_TAG
    return MPI_ANY_TAG;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_BAND"))
#ifdef MPI_BAND
    return MPI_BAND;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_BOR"))
#ifdef MPI_BOR
    return MPI_BOR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_BSEND_OVERHEAD"))
#ifdef MPI_BSEND_OVERHEAD
    return MPI_BSEND_OVERHEAD;
#else
      goto not_there;
#endif
  if (strEQ(name, "MPI_BXOR"))
#ifdef MPI_BXOR
    return MPI_BXOR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_BYTE"))
#ifdef MPI_BYTE
    return MPI_BYTE;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_CART"))
#ifdef MPI_CART
    return MPI_CART;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_CHAR"))
#ifdef MPI_CHAR
    return MPI_CHAR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_CHARACTER"))
#ifdef MPI_CHARACTER
    return MPI_CHARACTER;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_COMM_NULL"))
#ifdef MPI_COMM_NULL
    return MPI_COMM_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_COMM_SELF"))
#ifdef MPI_COMM_SELF
    return MPI_COMM_SELF;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_COMM_WORLD"))
#ifdef MPI_COMM_WORLD
    return MPI_COMM_WORLD;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_COMPLEX"))
#ifdef MPI_COMPLEX
    return MPI_COMPLEX;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_CONGRUENT"))
#ifdef MPI_CONGRUENT
    return MPI_CONGRUENT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_DATATYPE_NULL"))
#ifdef MPI_DATATYPE_NULL
    return MPI_DATATYPE_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_DOUBLE"))
#ifdef MPI_DOUBLE
    return MPI_DOUBLE;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_DOUBLE_COMPLEX"))
#ifdef MPI_DOUBLE_COMPLEX
    return MPI_DOUBLE_COMPLEX;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_DOUBLE_INT"))
#ifdef MPI_DOUBLE_INT
    return MPI_DOUBLE_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_DOUBLE_PRECISION"))
#ifdef MPI_DOUBLE_PRECISION
    return MPI_DOUBLE_PRECISION;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_ERRHANDLER_NULL"))
#ifdef MPI_ERRHANDLER_NULL
    return MPI_ERRHANDLER_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_ERRORS_ARE_FATAL"))
#ifdef MPI_ERRORS_ARE_FATAL
    return MPI_ERRORS_ARE_FATAL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_ERRORS_RETURN"))
#ifdef MPI_ERRORS_RETURN
    return MPI_ERRORS_RETURN;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_FLOAT"))
#ifdef MPI_FLOAT
    return MPI_FLOAT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_FLOAT_INT"))
#ifdef MPI_FLOAT_INT
    return MPI_FLOAT_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_GRAPH"))
#ifdef MPI_GRAPH
    return MPI_GRAPH;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_GROUP_EMPTY"))
#ifdef MPI_GROUP_EMPTY
    return MPI_GROUP_EMPTY;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_GROUP_NULL"))
#ifdef MPI_GROUP_NULL
    return MPI_GROUP_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_HOST"))
#ifdef MPI_HOST
    return MPI_HOST;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_IDENT"))
#ifdef MPI_IDENT
    return MPI_IDENT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_INT"))
#ifdef MPI_INT
    return MPI_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_INTEGER"))
#ifdef MPI_INTEGER
    return MPI_INTEGER;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_IO"))
#ifdef MPI_IO
    return MPI_IO;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_KEYVAL_INVALID"))
#ifdef MPI_KEYVAL_INVALID
    return MPI_KEYVAL_INVALID;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LAND"))
#ifdef MPI_LAND
    return MPI_LAND;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LB"))
#ifdef MPI_LB
    return MPI_LB;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LOGICAL"))
#ifdef MPI_LOGICAL
    return MPI_LOGICAL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LONG"))
#ifdef MPI_LONG
    return MPI_LONG;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LONG_DOUBLE"))
#ifdef MPI_LONG_DOUBLE
    return MPI_LONG_DOUBLE;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LONG_DOUBLE_INT"))
#ifdef MPI_LONG_DOUBLE_INT
    return MPI_LONG_DOUBLE_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LONG_INT"))
#ifdef MPI_LONG_INT
    return MPI_LONG_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LONG_LONG_INT"))
#ifdef MPI_LONG_LONG_INT
    return MPI_LONG_LONG_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LOR"))
#ifdef MPI_LOR
    return MPI_LOR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_LXOR"))
#ifdef MPI_LXOR
    return MPI_LXOR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MAX"))
#ifdef MPI_MAX
    return MPI_MAX;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MAXLOC"))
#ifdef MPI_MAXLOC
    return MPI_MAXLOC;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MAX_ERROR_STRING"))
#ifdef MPI_MAX_ERROR_STRING
    return MPI_MAX_ERROR_STRING;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MAX_NAME_STRING"))
#ifdef MPI_MAX_NAME_STRING
    return MPI_MAX_NAME_STRING;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MAX_PROCESSOR_NAME"))
#ifdef MPI_MAX_PROCESSOR_NAME
    return MPI_MAX_PROCESSOR_NAME;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MIN"))
#ifdef MPI_MIN
    return MPI_MIN;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_MINLOC"))
#ifdef MPI_MINLOC
    return MPI_MINLOC;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_OP_NULL"))
#ifdef MPI_OP_NULL
    return MPI_OP_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_PACKED"))
#ifdef MPI_PACKED
    return MPI_PACKED;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_PROC_NULL"))
#ifdef MPI_PROC_NULL
    return MPI_PROC_NULL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_PROD"))
#ifdef MPI_PROD
    return MPI_PROD;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_REAL"))
#ifdef MPI_REAL
    return MPI_REAL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_SHORT"))
#ifdef MPI_SHORT
    return MPI_SHORT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_SHORT_INT"))
#ifdef MPI_SHORT_INT
    return MPI_SHORT_INT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_SIMILAR"))
#ifdef MPI_SIMILAR
    return MPI_SIMILAR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_STRING"))
#ifdef MPI_STRING
    return MPI_STRING;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_SUBVERSION"))
#ifdef MPI_SUBVERSION
    return MPI_SUBVERSION;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_SUM"))
#ifdef MPI_SUM
    return MPI_SUM;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_TAG_UB"))
#ifdef MPI_TAG_UB
    return MPI_TAG_UB;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UB"))
#ifdef MPI_UB
    return MPI_UB;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNDEFINED"))
#ifdef MPI_UNDEFINED
    return MPI_UNDEFINED;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNDEFINED_RANK"))
#ifdef MPI_UNDEFINED_RANK
    return MPI_UNDEFINED_RANK;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNEQUAL"))
#ifdef MPI_UNEQUAL
    return MPI_UNEQUAL;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNSIGNED"))
#ifdef MPI_UNSIGNED
    return MPI_UNSIGNED;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNSIGNED_CHAR"))
#ifdef MPI_UNSIGNED_CHAR
    return MPI_UNSIGNED_CHAR;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNSIGNED_LONG"))
#ifdef MPI_UNSIGNED_LONG
    return MPI_UNSIGNED_LONG;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_UNSIGNED_SHORT"))
#ifdef MPI_UNSIGNED_SHORT
    return MPI_UNSIGNED_SHORT;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_VERSION"))
#ifdef MPI_VERSION
    return MPI_VERSION;
#else
    goto not_there;
#endif
  if (strEQ(name, "MPI_WTIME_IS_GLOBAL"))
#ifdef MPI_WTIME_IS_GLOBAL
    return MPI_WTIME_IS_GLOBAL;
#else
    goto not_there;
#endif

  errno = EINVAL;
  return 0;

  goto not_there;  /* -Wall */
 not_there:
  errno = ENOENT;
  return 0;
}

int MPIpm_bufsize(MPI_Datatype datatype, SV* scalar, int count) {
  switch (datatype) {

  case MPI_INT:
  case MPI_INTEGER:
    return count * sizeof(int);
    break;

  case MPI_DOUBLE:
    return count * sizeof(double);
    break;

#ifndef FLOAT_HACK    
  case MPI_FLOAT:
    return count * sizeof(float);
    break;
#endif

  case MPI_CHAR:
    return (count + 1) * sizeof(char);
    break;

  case MPI_UNSIGNED_SHORT:
    croak("Unsupported datatype %d in MPIpm_bufsize()", datatype);
    break;
   
  default:
    croak("Unrecognized datatype %d in MPIpm_bufsize()", datatype);
    break;
  }
  return(-1);
}

int MPIpm_packarray(void** buf, AV* array, MPI_Datatype datatype, int count) {
    int i,len=0,avl;
    SV** sv_tmp;
    char* t;

    avl = av_len(array);
    if(avl < 0) {
	*buf = NULL;
	return 0;
    }
    if(count > 0) {
	avl = count-1;
    }

    switch(datatype) {
      case MPI_STRING:
	for(i=0 ; i<=avl ; i++) {
	    sv_tmp = av_fetch(array, i, 0);
            if(sv_tmp != NULL)
		len += SvCUR(*sv_tmp)+1+sizeof(int);
	}

	t = *buf = (void*) malloc(len);

	for(i=0;i<=avl;i++) {
	    int n;
	    sv_tmp = av_fetch(array, i, 0);
	    n = htonl(SvCUR(*sv_tmp)+1);
	    memcpy(t, &n, sizeof(int));
	    memcpy(t+sizeof(int), SvPV(*sv_tmp, PL_na), SvCUR(*sv_tmp)+1);
	    t += SvCUR(*sv_tmp)+1+sizeof(int);
	}
	return(len);
	break;
      case MPI_INT:
      case MPI_INTEGER:
        len = (avl+1) * sizeof(int);
	t = *buf = (void*) malloc(len);

	for(i=0;i<=avl;i++) {
	    int n;
	    sv_tmp = av_fetch(array, i, 0);
	    n = SvIV(*sv_tmp);
	    memcpy(t, &n, sizeof(int));
	    t += sizeof(int);
	}
	return(len);
#ifndef FLOAT_HACK    
      case MPI_FLOAT:
        len = (avl+1) * sizeof(float);
	t = *buf = (void*) malloc(len);

	for(i=0;i<=avl;i++) {
	    float n;
	    sv_tmp = av_fetch(array, i, 0);
	    n = (float)SvNV(*sv_tmp);
	    memcpy(t, &n, sizeof(float));
	    t += sizeof(float);
	}
	return(len);
#endif
      case MPI_DOUBLE:
        len = (avl+1) * sizeof(double);
	t = *buf = (void*) malloc(len);

	for(i=0;i<=avl;i++) {
	    double n;
	    sv_tmp = av_fetch(array, i, 0);
	    n = SvNV(*sv_tmp);
	    memcpy(t, &n, sizeof(double));
	    t += sizeof(double);
	}
	return(len);
      default:
	croak("Unrecognized or unsupported datatype %d in MPIpm_packscalar()",
	      datatype);
	break;
    } 
    return(-1);
}

void MPIpm_unpackarray(void* buf, AV** array, MPI_Datatype datatype, int count) 
{
    int i;
    char* t = buf;

    av_clear(*array);
    switch(datatype) {
      case MPI_STRING:
	croak("Unrecognized or unsupported datatype %d in MPIpm_packscalar()",
	      datatype);
        break;
      case MPI_INT:
      case MPI_INTEGER:
        av_extend(*array, count);
        for(i=0;i<count;i++) {
            int u;
	    memcpy(&u, t, sizeof(int));
#ifdef SEND_DEBUG
	    printf("[%d] u=%d\n", getpid(), u);
#endif
            av_push(*array, newSViv(u));
	    t += sizeof(int);
	}
	break;
#ifndef FLOAT_HACK
      case MPI_FLOAT:
	break;
#endif
      case MPI_DOUBLE:
        av_extend(*array, count);
        for(i=0;i<count;i++) {
            double u;
	    memcpy(&u, t, sizeof(double));
#ifdef SEND_DEBUG
	    printf("[%d] u=%g\n", getpid(), u);
#endif
            av_push(*array, newSVnv(u));
	    t += sizeof(double);
	}
	break;
      default:
	croak("Unrecognized or unsupported datatype %d in MPIpm_packscalar()",
	      datatype);
	break;
    } 
}

void MPIpm_packscalar(char* buf, SV* scalar, MPI_Datatype datatype) {
  switch (datatype) {

  case MPI_INT:
  case MPI_INTEGER:
    *(int*)buf = SvIV(scalar);
    break;

  case MPI_DOUBLE:
    *(double*)buf = SvNV(scalar);
    break;

#ifndef FLOAT_HACK        
  case MPI_FLOAT:
    *(float*)buf = (float)SvNV(scalar);
    break;
#endif

  case MPI_CHAR:
    memcpy(buf, SvPV(scalar, PL_na), SvCUR(scalar));
    *(char*)(buf + SvCUR(scalar)) = '\0';
    break;

  case MPI_UNSIGNED_SHORT:
    croak("Unsupported datatype %d in MPIpm_packscalar()", datatype);
    break;
   
  default:
    croak("Unrecognized datatype %d in MPIpm_packscalar()", datatype);
    break;
  }

}

void MPIpm_unpackscalar(void* buf, SV* scalar, MPI_Datatype datatype) {
  switch (datatype) {

  case MPI_INT:
  case MPI_INTEGER:
    sv_setiv(scalar, *(int*)buf);
    break;

  case MPI_DOUBLE:
    sv_setnv(scalar, *(double*)buf);
    break;

#ifndef FLOAT_HACK        
  case MPI_FLOAT:
    sv_setnv(scalar, (double)*(float*)buf);    
    break;
#endif

  case MPI_CHAR:
    sv_setpv(scalar, buf);
    break;

  case MPI_UNSIGNED_SHORT:
    croak("Unsupported datatype %d in MPIpm_unpackscalar()", datatype);
    break;
   
  default:
    croak("Unrecognized datatype %d in MPIpm_unpackscalar()", datatype);
    break;
  }

}

void MPIpm_errhandler(char* func, int mpi_errno) {
  char buf[1024];
  SV *sv_errno, *sv_errstr;

  switch (mpi_errno) {
  case MPI_SUCCESS:
    return;
  case MPI_ERR_BUFFER:
    sprintf(buf,"%s: MPI_ERR_BUFFER: Invald buffer pointer",func);
    break;
  case MPI_ERR_COUNT:
    sprintf(buf,"%s: MPI_ERR_COUNT: Invalid count argument",func);
    break;
  case MPI_ERR_TYPE:
    sprintf(buf,"%s: MPI_ERR_TYPE: Invalid datatype argument",func);
    break;
  case MPI_ERR_TAG:
    sprintf(buf,"%s: MPI_ERR_TAG: Invalid tag argument",func);
    break;
  case MPI_ERR_COMM:
    sprintf(buf,"%s: MPI_ERR_COMM: Invalid Communicator",func);
    break;
  case MPI_ERR_RANK:
    sprintf(buf,"%s: MPI_ERR_RANK: Invalid rank",func);
    break;
  case MPI_ERR_ROOT:
    sprintf(buf,"%s: MPI_ERR_ROOT: Invalid root",func);
    break;
  case MPI_ERR_GROUP:
    sprintf(buf,"%s: MPI_ERR_GROUP: Null group passed to function",func);
    break;
  case MPI_ERR_OP:
    sprintf(buf,"%s: MPI_ERR_OP: Invalid operation",func);
    break;
  case MPI_ERR_TOPOLOGY:
    sprintf(buf,"%s: MPI_ERR_TOPOLOGY: Invalid topology",func);
    break;
  case MPI_ERR_DIMS:
    sprintf(buf,"%s: MPI_ERR_DIMS: Invalid dimension argument",func);
    break;
  case MPI_ERR_ARG:
    sprintf(buf,"%s: MPI_ERR_ARG: Invalid argument",func);
    break;
  case MPI_ERR_UNKNOWN:
    sprintf(buf,"%s: MPI_ERR_UNKNOWN: Unknown error",func);
    break;
  case MPI_ERR_TRUNCATE:
    sprintf(buf,"%s: MPI_ERR_TRUNCATE: message truncated on receive",func);
    break;
  case MPI_ERR_OTHER:
    sprintf(buf,"%s: MPI_ERR_OTHER: other error",func);
    break;
  case MPI_ERR_INTERN:
    sprintf(buf,"%s: MPI_ERR_INTERN: internal error code",func);
    break;
  case MPI_ERR_IN_STATUS:
    sprintf(buf,"%s: MPI_ERR_IN_STATUS: look in status for error value",func);
    break;
  case MPI_ERR_PENDING:
    sprintf(buf,"%s: MPI_ERR_PENDING: pending request",func);
    break; 
  case MPI_ERR_REQUEST:
    sprintf(buf,"%s: MPI_ERR_REQUEST: illegal mpi_request handle",func);
    break;
  default:
    sprintf(buf,"%s: MPI Error %d occurred",func,errno);
    break;
  }

  sv_errno = perl_get_sv("$MPI:errno", TRUE);
  sv_errstr = perl_get_sv("$MPI:errstr", TRUE); 
  sv_setiv(sv_errno,mpi_errno);
  sv_setpv(sv_errstr,buf);

  if (SvTRUE(perl_get_sv("$MPI::exceptions", TRUE))) {
    die(buf);
  }
}
