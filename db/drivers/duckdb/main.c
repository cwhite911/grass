#include <stdlib.h>

#include <grass/gis.h>
#include <grass/dbmi.h>
#include "dbdriver.h"

#include "globals.h"

// MYSQL *connection;       /* Database connection */
dbString *errMsg = NULL; /* error message */

int main(int argc, char *argv[])
{
    init_dbdriver();
    exit(db_driver(argc, argv));
}
