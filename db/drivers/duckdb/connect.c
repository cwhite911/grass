
// duckdb_state duckdb_open(const char *path, duckdb_database *out_database);
// duckdb_state duckdb_open_ext(const char *path, duckdb_database *out_database,
// duckdb_config config, char **out_error); void duckdb_close(duckdb_database
// *database); duckdb_state duckdb_connect(duckdb_database database,
// duckdb_connection *out_connection); void duckdb_interrupt(duckdb_connection
// connection); duckdb_query_progress_type
// duckdb_query_progress(duckdb_connection connection); void
// duckdb_disconnect(duckdb_connection *connection); const char
// *duckdb_library_version();
#include <stdlib.h>
#include <string.h>
#include <grass/dbmi.h>
#include <grass/gis.h>
#include <grass/glocale.h>

#include "globals.h"
#include "proto.h"

int connect_duckdb(const char *path, duckdb_database *db)
{

    if (duckdb_open(NULL, &db) == DuckDBError) {
        db_d_report_error();
        return DB_FAILED;
    }
    open_duckdb(path, db);
}

int open_duckdb(const char *path, duckdb_database *db)
{
    if (duckdb_open(path, &db) == DuckDBError) {
        db_d_report_error();
        return DB_FAILED;
    }
    return DB_OK;
}
