#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <grass/gis.h>
#include <grass/dbmi.h>
#include <grass/glocale.h>
#include "globals.h"
#include "proto.h"

int db__driver_open_database(dbHandle *handle)
{
    const char *name;
    name = db_get_handle_dbname(handle);
}

int db__driver_close_database(void)
{
}

int db__driver_create_database(dbHandle *handle)
{
}

int db__driver_delete_database(dbHandle *handle)
{
}

int db__driver_describe_table(dbString *table_name, dbTable **table)
{
}
