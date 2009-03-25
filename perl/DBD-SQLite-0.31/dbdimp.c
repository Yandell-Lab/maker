/* $Id: dbdimp.c,v 1.44 2004/02/14 17:36:38 matt Exp $ */

#include "SQLiteXS.h"

DBISTATE_DECLARE;

#ifndef SvPV_nolen
#define SvPV_nolen(x) SvPV(x,PL_na)
#endif

#ifndef call_method
#define call_method(x,y) perl_call_method(x,y)
#endif

#ifndef call_sv
#define call_sv(x,y) perl_call_sv(x,y)
#endif

void
sqlite_init(dbistate_t *dbistate)
{
    dTHR;
    DBIS = dbistate;
}

static void
sqlite_error(SV *h, imp_xxh_t *imp_xxh, int rc, char *what)
{
    dTHR;

    SV *errstr = DBIc_ERRSTR(imp_xxh);
    sv_setiv(DBIc_ERR(imp_xxh), (IV)rc);
    sv_setpv(errstr, what);
    
    if (DBIS->debug >= 3) {
        PerlIO_printf(DBILOGFP, "sqlite error %d recorded: %s\n",
            rc, SvPV_nolen(errstr));
    }
}

int
sqlite_db_login(SV *dbh, imp_dbh_t *imp_dbh, char *dbname, char *user, char *pass)
{
    dTHR;
    int retval;
    char *errmsg = NULL;

    if (DBIS->debug >= 3) {
        PerlIO_printf(DBILOGFP, "    login '%s' (version %s, encoding %s)\n",
            dbname, sqlite_version, sqlite_encoding);
    }

    if ((imp_dbh->db = sqlite_open(dbname, 0, &errmsg)) == NULL) {
	sqlite_error(dbh, (imp_xxh_t*)imp_dbh, 1, errmsg);
        sqlite_freemem(errmsg);
        return FALSE;
    }
    DBIc_IMPSET_on(imp_dbh);

    imp_dbh->in_tran = FALSE;
    imp_dbh->no_utf8_flag = FALSE;
    imp_dbh->functions = newAV();
    imp_dbh->aggregates = newAV();
    
    imp_dbh->handle_binary_nulls = FALSE;

    sqlite_busy_timeout(imp_dbh->db, SQL_TIMEOUT);

    if (retval = sqlite_exec(imp_dbh->db, "PRAGMA empty_result_callbacks = ON",
        NULL, NULL, &errmsg)
        != SQLITE_OK)
    {
        /*  warn("failed to set pragma: %s\n", errmsg); */
	    sqlite_error(dbh, (imp_xxh_t*)imp_dbh, retval, errmsg);
        sqlite_freemem(errmsg);
        return FALSE;
    }

    if (retval = sqlite_exec(imp_dbh->db, "PRAGMA show_datatypes = ON",
        NULL, NULL, &errmsg)
        != SQLITE_OK)
    {
        /*  warn("failed to set pragma: %s\n", errmsg); */
    	sqlite_error(dbh, (imp_xxh_t*)imp_dbh, retval, errmsg);
        sqlite_freemem(errmsg);
        return FALSE;
    }

    DBIc_ACTIVE_on(imp_dbh);

    return TRUE;
}

int
sqlite_db_disconnect (SV *dbh, imp_dbh_t *imp_dbh)
{
    dTHR;
    DBIc_ACTIVE_off(imp_dbh);

    if (DBIc_is(imp_dbh, DBIcf_AutoCommit) == FALSE) {
        sqlite_db_rollback(dbh, imp_dbh);
    }

    sqlite_close(imp_dbh->db);
    imp_dbh->db = NULL;

    av_undef(imp_dbh->functions);
    imp_dbh->functions = (AV *)NULL;

    av_undef(imp_dbh->aggregates);
    imp_dbh->aggregates = (AV *)NULL;

    return TRUE;
}

void
sqlite_db_destroy (SV *dbh, imp_dbh_t *imp_dbh)
{
    dTHR;
    if (DBIc_ACTIVE(imp_dbh)) {
        sqlite_db_disconnect(dbh, imp_dbh);
    }
    DBIc_IMPSET_off(imp_dbh);
}

int
sqlite_db_rollback(SV *dbh, imp_dbh_t *imp_dbh)
{
    dTHR;
    int retval;
    char *errmsg;

    if (imp_dbh->in_tran) {
        if (retval = sqlite_exec(imp_dbh->db, "ROLLBACK TRANSACTION",
            NULL, NULL, &errmsg)
            != SQLITE_OK)
        {
	    sqlite_error(dbh, (imp_xxh_t*)imp_dbh, retval, errmsg);
            sqlite_freemem(errmsg);
            return FALSE;
        }
        imp_dbh->in_tran = FALSE;
    }

    return TRUE;
}

int
sqlite_db_commit(SV *dbh, imp_dbh_t *imp_dbh)
{
    dTHR;
    int retval;
    char *errmsg;

    if (DBIc_is(imp_dbh, DBIcf_AutoCommit)) {
        warn("commit ineffective with AutoCommit");
        return TRUE;
    }

    if (imp_dbh->in_tran) {
        if (retval = sqlite_exec(imp_dbh->db, "COMMIT TRANSACTION",
            NULL, NULL, &errmsg)
            != SQLITE_OK)
        {
	    sqlite_error(dbh, (imp_xxh_t*)imp_dbh, retval, errmsg);
            sqlite_freemem(errmsg);
            return FALSE;
        }
        imp_dbh->in_tran = FALSE;
    }
    return TRUE;
}

int
sqlite_discon_all(SV *drh, imp_drh_t *imp_drh)
{
    dTHR;
    return FALSE; /* no way to do this */
}

void
sqlite_st_parse_sql(imp_sth_t *imp_sth, char *statement)
{
    D_imp_dbh_from_sth;
    bool in_literal = FALSE;
    SV *chunk;
    int num_params = 0;

    chunk = NEWSV(0, strlen(statement));
    sv_setpv(chunk, "");

    /* warn("parsing: %s\n", statement); */

    while (*statement) {
        /*  warn("parse: %c => %s\n", *statement, SvPV_nolen(chunk)); */
        if (*statement == '\'') {
            if (in_literal) {
                /*  either end of literal, or escape */
                if (statement[1] && statement[1] == '\'') {
                    statement++;
                    sv_catpvn(chunk, "''", 2);
                }
                else {
                    sv_catpvn(chunk, "'", 1);
                    in_literal = FALSE;
                }
            }
            else {
                in_literal = TRUE;
                sv_catpvn(chunk, "'", 1);
            }
        }
        else if (*statement == '?') {
            if (in_literal) {
                sv_catpvn(chunk, "?", 1);
            }
            else {
                num_params++;
		if (!imp_dbh->no_utf8_flag) {
                    /* sv_utf8_encode(chunk); */
		}
                av_push(imp_sth->sql, chunk);
                chunk = NEWSV(0, 20);
                sv_setpvn(chunk, "", 0);
            }
        }
        else {
            sv_catpvn(chunk, statement, 1);
        }
        statement++;
    }
    av_push(imp_sth->sql, chunk);
    DBIc_NUM_PARAMS(imp_sth) = num_params;
}

int
sqlite_st_prepare (SV *sth, imp_sth_t *imp_sth,
                char *statement, SV *attribs)
{
    dTHR;
    D_imp_dbh_from_sth;

    if (!DBIc_ACTIVE(imp_dbh)) {
      die("prepare on an inactive database handle");
    }

    /* warn("prepare statement\n"); */
    imp_sth->nrow = 0;
    imp_sth->ncols = 0;
    imp_sth->params = newAV();
    imp_sth->sql = newAV();
    imp_sth->results = 0;
    imp_sth->coldata = 0;
    imp_sth->retval = SQLITE_OK;
    sqlite_st_parse_sql(imp_sth, statement);

    return TRUE;
}

char *
sqlite_quote(imp_dbh_t *imp_dbh, SV *val)
{
    STRLEN len;
    char *cval = SvPV(val, len);
    SV *ret = sv_2mortal(NEWSV(0, SvCUR(val) + 2));
    sv_setpvn(ret, "", 0);

    while (len) {
      switch (*cval) {
        case '\'':
          sv_catpvn(ret, "''", 2);
          break;
        case 0:
          if (imp_dbh->handle_binary_nulls) {
            sv_catpvn(ret, "\\0", 2);
            break;
          }
          else {
            die("attempt to quote binary null without sqlite_handle_binary_nulls on");
          }
        case '\\':
          if (imp_dbh->handle_binary_nulls) {
            sv_catpvn(ret, "\\\\", 2);
            break;
          }
        default:
          sv_catpvn(ret, cval, 1);
      }
      *cval++; len--;
    }
    return SvPV_nolen(ret);
}

char *
sqlite_decode(imp_dbh_t *imp_dbh, char *input, size_t *len)
{
  char *ret = malloc(*len);
  char *swit = ret;

  while (*input) {
    switch (*input) {
      case '\\':
        if (imp_dbh->handle_binary_nulls && input[1] && input[1] == '0') {
          *swit++ = '\0';
          *input++;
          (*len)--;
          break;
        }
        else if (imp_dbh->handle_binary_nulls && input[1] && input[1] == '\\') {
          *swit++ = '\\';
          *input++;
          (*len)--;
          break;
        }
      default:
        *swit++ = *input;
    }
    *input++;
  }
  return ret;
}

int
_sqlite_fetch_row (imp_sth_t *imp_sth)
{
    while (1)
    {
        imp_sth->retval = sqlite_step(imp_sth->vm, 
            &(imp_sth->ncols), (const char ***)&(imp_sth->results), (const char ***)&(imp_sth->coldata));
        if (imp_sth->retval == SQLITE_BUSY) {
            continue;
        }
        break;
    }
    /* warn("step got: %d\nCol1: %s\n", imp_sth->retval, imp_sth->coldata[0]); */
    return imp_sth->retval;
}

int
sqlite_st_execute (SV *sth, imp_sth_t *imp_sth)
{
    dTHR;
    D_imp_dbh_from_sth;
    SV *sql;
    I32 pos = 0;
    char *errmsg;
    int num_params = DBIc_NUM_PARAMS(imp_sth);
    I32 i;
    int retval;

    /* warn("execute\n"); */

    if (DBIc_ACTIVE(imp_sth)) {
        sqlite_st_finish(sth, imp_sth);
    }

    sql = sv_2mortal(newSVsv(AvARRAY(imp_sth->sql)[pos++]));

    for (i = 0; i < num_params; i++) {
        SV *value = av_shift(imp_sth->params);
        if (value && SvOK(value)) {
            /* warn("binding param: %s\n", SvPV_nolen(value); */
            sv_catpvn(sql, "'", 1);
            sv_catpv(sql, sqlite_quote(imp_dbh, value));
            sv_catpvn(sql, "'", 1);
            /* warn("inserting string length: %d\n", SvCUR(sql)); */
        }
        else {
            /* warn("binding NULL\n"); */
            sv_catpvn(sql, "NULL", 4);
        }
        if (value) {
            SvREFCNT_dec(value);
        }
        sv_catsv(sql, AvARRAY(imp_sth->sql)[pos++]);
    }
    /* warn("Executing: %s;\n", SvPV_nolen(sql)); */

    if ( (!DBIc_is(imp_dbh, DBIcf_AutoCommit)) && (!imp_dbh->in_tran) ) {
        if (retval = sqlite_exec(imp_dbh->db, "BEGIN TRANSACTION",
            NULL, NULL, &errmsg)
            != SQLITE_OK)
        {
            sqlite_error(sth, (imp_xxh_t*)imp_sth, retval, errmsg);
            sqlite_freemem(errmsg);
            return -2;
        }
        imp_dbh->in_tran = TRUE;
    }

    imp_sth->results = NULL;
    if (retval = sqlite_compile(imp_dbh->db, SvPV_nolen(sql), 0, &(imp_sth->vm), &errmsg) != SQLITE_OK)
    {
        sqlite_error(sth, (imp_xxh_t*)imp_sth, retval, errmsg);
        sqlite_freemem(errmsg);
        return -2;
    }
    
    if (_sqlite_fetch_row(imp_sth) == SQLITE_ERROR) {
        sqlite_finalize(imp_sth->vm, &errmsg);
        sqlite_error(sth, (imp_xxh_t*)imp_sth, imp_sth->retval, errmsg);
        sqlite_freemem(errmsg);
        return -2;
    }
    
    DBIc_NUM_FIELDS(imp_sth) = imp_sth->ncols;
    imp_sth->nrow = -1;
    
    /* warn("Execute returned %d cols\n", imp_sth->ncols); */
    if (imp_sth->ncols == 0) {
        sqlite_finalize(imp_sth->vm, 0);
        imp_sth->nrow = sqlite_changes(imp_dbh->db);
        DBIc_IMPSET_on(imp_sth);
        return imp_sth->nrow;
    }
    
    DBIc_ACTIVE_on(imp_sth);
    /* warn("exec ok - %d rows, %d cols\n", imp_sth->nrow, imp_sth->ncols); */
    DBIc_IMPSET_on(imp_sth);
    return 0;
}

int
sqlite_st_rows (SV *sth, imp_sth_t *imp_sth)
{
    return imp_sth->nrow;
}

int
sqlite_bind_ph (SV *sth, imp_sth_t *imp_sth,
                SV *param, SV *value, IV sql_type, SV *attribs,
                                int is_inout, IV maxlen)
{
    if (is_inout) {
        croak("InOut bind params not implemented");
    }
    /* warn("bind: %s => %s\n", SvPV_nolen(param), SvPV_nolen(value)); */
    if (sql_type >= SQL_NUMERIC && sql_type <= SQL_DOUBLE) {
        av_store(imp_sth->params, SvIV(param) - 1, newSVnv(SvNV(value)));
    }
    else {
        av_store(imp_sth->params, SvIV(param) - 1, SvREFCNT_inc(value));
    }
}

AV *
sqlite_st_fetch (SV *sth, imp_sth_t *imp_sth)
{
    AV *av;
    D_imp_dbh_from_sth;
    int numFields = DBIc_NUM_FIELDS(imp_sth);
    int chopBlanks = DBIc_is(imp_sth, DBIcf_ChopBlanks);
    int i;

    /* warn("current_entry == %d\nnumFields == %d\nnrow == %d",
        current_entry, numFields, imp_sth->nrow); */

    /*
    if (!DBIc_ACTIVE(imp_sth)) {
        return Nullav;
    }
    */
    
    if ((imp_sth->retval == SQLITE_DONE) || (imp_sth->retval == SQLITE_ERROR)) {
        sqlite_st_finish(sth, imp_sth);
        return Nullav;
    }
    
    if (imp_sth->nrow == -1) {
        imp_sth->nrow++;
    }
    imp_sth->nrow++;
    
    av = DBIS->get_fbav(imp_sth);
    for (i = 0; i < numFields; i++) {
        char *val = imp_sth->results[i];
        /* warn("fetching: %d == %s\n", i, val); */
        if (val != NULL) {
            size_t len = strlen(val);
            char *decoded;
            if (chopBlanks) {
                val = savepv(val);
                while((len > 0) && (val[len-1] == ' ')) {
                    len--;
                }
                val[len] = '\0';
            }
            decoded = sqlite_decode(imp_dbh, val, &len);
            sv_setpvn(AvARRAY(av)[i], decoded, len);
            Safefree(decoded);
            if (chopBlanks) Safefree(val);

            if (!imp_dbh->no_utf8_flag) {
            /* sv_utf8_encode(AvARRAY(av)[i]); */
            }
        }
        else {
            sv_setsv(AvARRAY(av)[i], &PL_sv_undef);
        }
    }
    _sqlite_fetch_row(imp_sth);
    return av;
}

int
sqlite_st_finish (SV *sth, imp_sth_t *imp_sth)
{
    /* warn("finish statement\n"); */
    if (DBIc_ACTIVE(imp_sth)) {
        char *errmsg;
        DBIc_ACTIVE_off(imp_sth);
        if (imp_sth->retval = sqlite_finalize(imp_sth->vm, &errmsg) == SQLITE_ERROR) {
            warn("finalize failed! %s\n", errmsg);
            sqlite_error(sth, (imp_xxh_t*)imp_sth, imp_sth->retval, errmsg);
            sqlite_freemem(errmsg);
            return FALSE;
        }
    }
    return TRUE;
}

void
sqlite_st_destroy (SV *sth, imp_sth_t *imp_sth)
{
    /* warn("destroy statement\n"); */
    if (DBIc_ACTIVE(imp_sth)) {
        sqlite_st_finish(sth, imp_sth);
    }
    SvREFCNT_dec((SV*)imp_sth->sql);
    SvREFCNT_dec((SV*)imp_sth->params);
    DBIc_IMPSET_off(imp_sth);
}

int
sqlite_st_blob_read (SV *sth, imp_sth_t *imp_sth,
                int field, long offset, long len, SV *destrv, long destoffset)
{
    return 0;
}

int
sqlite_db_STORE_attrib (SV *dbh, imp_dbh_t *imp_dbh, SV *keysv, SV *valuesv)
{
    dTHR;
    char *key = SvPV_nolen(keysv);
    char *errmsg;
    int retval;

    if (strEQ(key, "AutoCommit")) {
        if (SvTRUE(valuesv)) {
            /* commit tran? */
            if ( (!DBIc_is(imp_dbh, DBIcf_AutoCommit)) && (imp_dbh->in_tran) ) {
                if (retval = sqlite_exec(imp_dbh->db, "COMMIT TRANSACTION",
                    NULL, NULL, &errmsg)
                    != SQLITE_OK)
                {
		    sqlite_error(dbh, (imp_xxh_t*)imp_dbh, retval, errmsg);
                    sqlite_freemem(errmsg);
                    return TRUE;
                }
                imp_dbh->in_tran = FALSE;
            }
        }
        DBIc_set(imp_dbh, DBIcf_AutoCommit, SvTRUE(valuesv));
        return TRUE;
    }
    else if (strEQ(key, "sqlite_no_utf8_flag") || strEQ(key, "NoUTF8Flag")) {
        warn("NoUTF8Flag is deprecated due to perl unicode weirdness\n");
        if (SvTRUE(valuesv)) {
            imp_dbh->no_utf8_flag = TRUE;
        }
        else {
            imp_dbh->no_utf8_flag = FALSE;
        }
        return TRUE;
    }
    else if (strEQ(key, "sqlite_handle_binary_nulls")) {
        if (SvTRUE(valuesv)) {
          imp_dbh->handle_binary_nulls = TRUE;
        }
        else {
          imp_dbh->handle_binary_nulls = FALSE;
        }
        return TRUE;
    }
    return FALSE;
}

SV *
sqlite_db_FETCH_attrib (SV *dbh, imp_dbh_t *imp_dbh, SV *keysv)
{
    dTHR;
    char *key = SvPV_nolen(keysv);

    if (strEQ(key, "sqlite_no_utf8_flag") || strEQ(key, "NoUTF8Flag")) {
        return newSViv(imp_dbh->no_utf8_flag ? 1 : 0);
    }
    if (strEQ(key, "sqlite_version")) {
        return newSVpv(sqlite_version, strlen(sqlite_version));
    }
    if (strEQ(key, "sqlite_encoding")) {
        return newSVpv(sqlite_encoding, strlen(sqlite_encoding));
    }
    return NULL;
}

int
sqlite_st_STORE_attrib (SV *sth, imp_sth_t *imp_sth, SV *keysv, SV *valuesv)
{
    char *key = SvPV_nolen(keysv);
    return FALSE;
}

SV *
sqlite_st_FETCH_attrib (SV *sth, imp_sth_t *imp_sth, SV *keysv)
{
    char *key = SvPV_nolen(keysv);
    SV *retsv = NULL;
    int i;

    if (!imp_sth->coldata) {
        return retsv;
    }

    i = DBIc_NUM_FIELDS(imp_sth);

    if (strEQ(key, "NAME")) {
        AV *av = newAV();
        av_extend(av, i);
        retsv = sv_2mortal(newRV(sv_2mortal((SV*)av)));
        while (--i >= 0) {
            char *fieldname = imp_sth->coldata[i];
            /* warn("Name [%d]: %s\n", i, fieldname); */
            char *dot = instr(fieldname, ".");
            if (dot) /* drop table name from field name */
                fieldname = ++dot;
            av_store(av, i, newSVpv(fieldname, 0));
        }
    }
    else if (strEQ(key, "PRECISION")) {
        AV *av = newAV();
        retsv = sv_2mortal(newRV(sv_2mortal((SV*)av)));
    }
    else if (strEQ(key, "TYPE")) {
        int i_base = i;
        AV *av = newAV();
        av_extend(av, i);
        retsv = sv_2mortal(newRV(sv_2mortal((SV*)av)));
        i = i * 2;
        while (--i >= i_base) {
            char *fieldname = imp_sth->coldata[i];
            /* warn("Type [%d]: %s\n", i, fieldname); */
            char *dot = instr(fieldname, ".");
            if (dot) /* drop table name from field name */
                fieldname = ++dot;
            av_store(av, i - i_base, newSVpv(fieldname, 0));
        }
    }
    else if (strEQ(key, "NULLABLE")) {
        AV *av = newAV();
        retsv = sv_2mortal(newRV(sv_2mortal((SV*)av)));
    }
    else if (strEQ(key, "SCALE")) {
        AV *av = newAV();
        retsv = sv_2mortal(newRV(sv_2mortal((SV*)av)));
    }
    else if (strEQ(key, "NUM_OF_FIELDS")) {
        retsv = sv_2mortal(newSViv(i));
    }

    return retsv;
}

static void
sqlite_db_set_result(sqlite_func *context, SV *result, int is_error )
{
    STRLEN len;
    char *s;
    
    if ( is_error ) {
        s = SvPV(result, len);
        sqlite_set_result_error( context, s, len );
        return;
    }

    if ( !SvOK(result) ) {
        sqlite_set_result_string( context, NULL, -1 );
    } else if( SvIOK(result) ) {
        sqlite_set_result_int( context, SvIV(result));
    } else if ( !is_error && SvIOK(result) ) {
        sqlite_set_result_double( context, SvNV(result));
    } else {
        s = SvPV(result, len);
        sqlite_set_result_string( context, s, len );
    } 
}

static void
sqlite_db_func_dispatcher(sqlite_func *context, int argc, const char **argv)
{
    dSP;
    int count;
    int i;
    SV *func;
    
    func = sqlite_user_data(context);
    
    ENTER;
    SAVETMPS;
    
    PUSHMARK(SP);
    for ( i=0; i < argc; i++ ) {
        SV *arg;
        
        if ( !argv[i] ) {
            arg = &PL_sv_undef;
        } else {
            arg = sv_2mortal( newSVpv(argv[i], 0));
        }
        
        XPUSHs(arg);
    }
    PUTBACK;

    count = call_sv(func, G_SCALAR|G_EVAL);
    
    SPAGAIN;

    /* Check for an error */
    if (SvTRUE(ERRSV) ) {
        sqlite_db_set_result( context, ERRSV, 1);
        POPs;
    } else if ( count != 1 ) {
        SV *err = sv_2mortal(newSVpvf( "function should return 1 argument, got %d",
                                       count ));

        sqlite_db_set_result( context, err, 1);
        /* Clear the stack */
        for ( i=0; i < count; i++ ) {
            POPs;
        }
    } else {
        sqlite_db_set_result( context, POPs, 0 );
    }
    
    PUTBACK;
    
    FREETMPS;
    LEAVE;
}

void
sqlite_db_create_function( SV *dbh, const char *name, int argc, SV *func )
{
    D_imp_dbh(dbh);
    int rv;
    
    /* Copy the function reference */
    SV *func_sv = newSVsv(func);
    av_push( imp_dbh->functions, func_sv );

    rv = sqlite_create_function( imp_dbh->db, name, argc, 
                                 sqlite_db_func_dispatcher, func_sv );
    if ( rv != SQLITE_OK )
    {
        croak( "sqlite_create_function failed with error %s", 
               sqlite_error_string(rv) );
    }
}

typedef struct aggrInfo aggrInfo;
struct aggrInfo {
  SV *aggr_inst;
  SV *err;
  int inited;
};

static void
sqlite_db_aggr_new_dispatcher( sqlite_func *context, aggrInfo *aggr_info )
{
    dSP;
    SV *pkg = NULL;
    int count = 0;

    aggr_info->err = NULL;
    aggr_info->aggr_inst = NULL;
    
    pkg = sqlite_user_data(context);
    if ( !pkg )
        return;

    ENTER;
    SAVETMPS;
    
    PUSHMARK(SP);
    XPUSHs( sv_2mortal( newSVsv(pkg) ) );
    PUTBACK;

    count = call_method ("new", G_EVAL|G_SCALAR);
    SPAGAIN;

    aggr_info->inited = 1;

    if ( SvTRUE( ERRSV ) ) {
        aggr_info->err =  newSVpvf ("error during aggregator's new(): %s",
                                    SvPV_nolen (ERRSV));
        POPs;
    } else if ( count != 1 ) {
        int i;
        
        aggr_info->err = newSVpvf( "new() should return one value, got %d", 
                                  count );
        /* Clear the stack */
        for ( i=0; i < count; i++ ) {
            POPs;
        }
    } else {
        SV *aggr = POPs;
        if ( SvROK(aggr) ) {
            aggr_info->aggr_inst = newSVsv(aggr);
        } else{
            aggr_info->err = newSVpvf( "new() should return a blessed reference" );
        }
    }

    PUTBACK;

    FREETMPS;
    LEAVE;

    return;
}

static void
sqlite_db_aggr_step_dispatcher (sqlite_func *context,
                                int argc, const char **argv)
{
    dSP;
    int i;
    aggrInfo *aggr;

    aggr = sqlite_aggregate_context (context, sizeof (aggrInfo));
    if ( !aggr )
        return;
    
    ENTER;
    SAVETMPS;
    
    /* initialize on first step */
    if ( !aggr->inited ) {
        sqlite_db_aggr_new_dispatcher( context, aggr );
    }

    if ( aggr->err || !aggr->aggr_inst ) 
        goto cleanup;

    PUSHMARK(SP);
    XPUSHs( sv_2mortal( newSVsv( aggr->aggr_inst ) ));
    for ( i=0; i < argc; i++ ) {
        SV *arg;
        
        if ( !argv[i] ) {
            arg = &PL_sv_undef;
        } else {
            arg = sv_2mortal( newSVpv(argv[i], 0));
        }

        XPUSHs(arg);
    }
    PUTBACK;

    call_method ("step", G_SCALAR|G_EVAL|G_DISCARD);

    /* Check for an error */
    if (SvTRUE(ERRSV) ) {
      aggr->err = newSVpvf( "error during aggregator's step(): %s",
                            SvPV_nolen(ERRSV));
      POPs;
    }

 cleanup:
    FREETMPS;
    LEAVE;
}

static void
sqlite_db_aggr_finalize_dispatcher( sqlite_func *context )
{
    dSP;
    aggrInfo *aggr, myAggr;
    int count = 0;

    aggr = sqlite_aggregate_context (context, sizeof (aggrInfo));
    
    ENTER;
    SAVETMPS;
    
    if ( !aggr ) {
        /* SQLite seems to refuse to create a context structure
           from finalize() */
        aggr = &myAggr;
        aggr->aggr_inst = NULL;
        aggr->err = NULL;
        sqlite_db_aggr_new_dispatcher (context, aggr);
    } 

    if  ( ! aggr->err && aggr->aggr_inst ) {
        PUSHMARK(SP);
        XPUSHs( sv_2mortal( newSVsv( aggr->aggr_inst )) );
        PUTBACK;

        count = call_method( "finalize", G_SCALAR|G_EVAL );
        SPAGAIN;

        if ( SvTRUE(ERRSV) ) {
            aggr->err = newSVpvf ("error during aggregator's finalize(): %s",
                                  SvPV_nolen(ERRSV) ) ;
            POPs;
        } else if ( count != 1 ) {
            int i;
            aggr->err = newSVpvf( "finalize() should return 1 value, got %d",
                                  count );
            /* Clear the stack */
            for ( i=0; i<count; i++ ) {
                POPs;
            }
        } else {
            sqlite_db_set_result( context, POPs, 0 );
        }
        PUTBACK;
    }
    
    if ( aggr->err ) {
        warn( "DBD::SQLite: error in aggregator cannot be reported to SQLite: %s",               SvPV_nolen( aggr->err ) );
        
        /* sqlite_db_set_result( context, aggr->err, 1 ); */
        SvREFCNT_dec( aggr->err );
        aggr->err = NULL;
    }
    
    if ( aggr->aggr_inst ) {
         SvREFCNT_dec( aggr->aggr_inst );
         aggr->aggr_inst = NULL;
    }
    
    FREETMPS;
    LEAVE;
}

void
sqlite_db_create_aggregate( SV *dbh, const char *name, int argc, SV *aggr_pkg )
{
    D_imp_dbh(dbh);
    int rv;
    
    /* Copy the aggregate reference */
    SV *aggr_pkg_copy = newSVsv(aggr_pkg);
    av_push( imp_dbh->aggregates, aggr_pkg_copy );

    rv = sqlite_create_aggregate( imp_dbh->db, name, argc, 
                                  sqlite_db_aggr_step_dispatcher, 
                                  sqlite_db_aggr_finalize_dispatcher, 
                                  aggr_pkg_copy );
    
    if ( rv != SQLITE_OK )
    {
        croak( "sqlite_create_aggregate failed with error %s", 
               sqlite_error_string(rv) );
    }
}

/* end */
