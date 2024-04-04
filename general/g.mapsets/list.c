#include <string.h>
#include <stdio.h>
#include <grass/gis.h>
#include <grass/glocale.h>
#include "local_proto.h"
#include <grass/parson.h>

void list_available_mapsets(const char **mapset_name, int nmapsets,
                            const char *fs)
{
    int n;

    G_message(_("Available mapsets:"));

    for (n = 0; n < nmapsets; n++) {
        fprintf(stdout, "%s", mapset_name[n]);
        if (n < nmapsets - 1) {
            if (strcmp(fs, "newline") == 0)
                fprintf(stdout, "\n");
            else if (strcmp(fs, "space") == 0)
                fprintf(stdout, " ");
            else if (strcmp(fs, "comma") == 0)
                fprintf(stdout, ",");
            else if (strcmp(fs, "tab") == 0)
                fprintf(stdout, "\t");
            else
                fprintf(stdout, "%s", fs);
        }
    }
    fprintf(stdout, "\n");
}

void list_accessible_mapsets(const char *fs)
{
    int n;
    const char *name;

    G_message(_("Accessible mapsets:"));

    for (n = 0; (name = G_get_mapset_name(n)); n++) {
        /* match each mapset to its numeric equivalent */
        fprintf(stdout, "%s", name);
        if (G_get_mapset_name(n + 1)) {
            if (strcmp(fs, "newline") == 0)
                fprintf(stdout, "\n");
            else if (strcmp(fs, "space") == 0)
                fprintf(stdout, " ");
            else if (strcmp(fs, "comma") == 0)
                fprintf(stdout, ",");
            else if (strcmp(fs, "tab") == 0)
                fprintf(stdout, "\t");
            else
                fprintf(stdout, "%s", fs);
        }
    }
    fprintf(stdout, "\n");
}

void list_accessible_mapsets_json(const char *fs)
{
    int n;
    const char *name;
    char *serialized_string = NULL;
    JSON_Value *root_value = NULL;
    JSON_Object *root_object;
    JSON_Array *root_array, *mapsets;

    // Create root json object
    root_value = json_value_init_object();
    root_object = json_value_get_object(root_value);

    // Create mapsets array
    root_array = json_value_init_array();
    json_object_set_value(root_object, "mapsets", root_array);
    mapsets = json_object_get_array(root_object, "mapsets");

    // Check that memory was allocated to root json object
    if (root_value == NULL) {
        G_fatal_error(_("Failed to initialize JSON. Out of memory?"));
    }

    // Check that memory was allocated to mapsets array
    if (mapsets == NULL) {
        G_fatal_error(_("Failed to initialize JSON array. Out of memory?"));
    }

    // Add mapsets to mapsets array
    for (n = 0; (name = G_get_mapset_name(n)); n++) {
        // Append mapset name to mapsets array
        json_array_append_string(mapsets, name);
    }

    // Add mapsets array to root object
    json_object_set_value(root_object, "mapsets", mapsets);

    // Serialize root object to string and print it to stdout
    serialized_string = json_serialize_to_string_pretty(root_value);
    puts(serialized_string);

    // Free memory
    json_free_serialized_string(serialized_string);
    json_value_free(root_value);
}

void list_avaliable_mapsets_json(const char **mapset_name, int nmapsets)
{
    int n;
    char *serialized_string = NULL;
    JSON_Value *root_value = NULL;
    JSON_Object *root_object;
    JSON_Array *root_array, *mapsets;

    root_value = json_value_init_object();
    root_object = json_value_get_object(root_value);

    // Create mapsets array
    root_array = json_value_init_array();
    json_object_set_value(root_object, "mapsets", root_array);
    mapsets = json_object_get_array(root_object, "mapsets");

    // Check that memory was allocated to root json object
    if (root_value == NULL) {
        G_fatal_error(_("Failed to initialize JSON. Out of memory?"));
    }

    if (mapsets == NULL) {
        G_fatal_error(_("Failed to initialize JSON array. Out of memory?"));
    }

    // Append mapsets to mapsets array
    for (n = 0; n < nmapsets; n++) {
        json_array_append_string(mapsets, mapset_name[n]);
    }

    // Add mapsets array to root object
    json_object_set_value(root_object, "mapsets", mapsets);

    // Serialize root object to string and print it to stdout
    serialized_string = json_serialize_to_string_pretty(root_value);
    puts(serialized_string);

    // Free memory
    json_free_serialized_string(serialized_string);
    json_value_free(root_value);
}
