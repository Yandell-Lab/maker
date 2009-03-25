/* $Id: SQLite.xs,v 1.5 2003/07/31 15:11:46 matt Exp $ */

#include "SQLiteXS.h"

DBISTATE_DECLARE;

MODULE = DBD::SQLite          PACKAGE = DBD::SQLite::db

PROTOTYPES: DISABLE

AV *
list_tables(dbh)
    SV *dbh
    CODE:
    {
        RETVAL = newAV();
    }
    OUTPUT:
        RETVAL

int
last_insert_rowid(dbh)
    SV *dbh
    CODE:
    {
        D_imp_dbh(dbh);
        RETVAL = sqlite_last_insert_rowid(imp_dbh->db);
    }
    OUTPUT:
        RETVAL

void
create_function(dbh, name, argc, func)
    SV *dbh
    char *name
    int argc
    SV *func
    CODE:
    {
        sqlite_db_create_function( dbh, name, argc, func );
    }

void
create_aggregate(dbh, name, argc, aggr)
    SV *dbh
    char *name
    int argc
    SV *aggr
    CODE:
    {
        sqlite_db_create_aggregate( dbh, name, argc, aggr );
    }

MODULE = DBD::SQLite          PACKAGE = DBD::SQLite::st

PROTOTYPES: DISABLE

MODULE = DBD::SQLite          PACKAGE = DBD::SQLite

INCLUDE: SQLite.xsi
