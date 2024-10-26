#include <stdio.h>
#include <grass/gis.h>
#include <grass/dbmi.h>
#include "proto.h"
#include "globals.h"

/* init error message */
void init_error(void)
{
    db_d_init_error("DuckDB");
}
