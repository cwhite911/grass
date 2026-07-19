/****************************************************************************
 *
 * MODULE:       r.class.breaks
 * AUTHOR(S):    Corey White <corey.white openplains.com>
 * PURPOSE:      Compute class breaks for a raster map from a histogram of
 *               the cell values, using Jenks natural breaks, equal
 *               intervals, or quantiles
 *
 * COPYRIGHT:    (C) 2026 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/arraystats.h>
#include <grass/gis.h>
#include <grass/gjson.h>
#include <grass/glocale.h>
#include <grass/raster.h>

enum OutputFormat { PLAIN, CSV, JSON };

static int num_slots;
static DCELL min, max;
static DCELL slot_size;
static int integer_mode; /* one slot per integer value: breaks are exact */

static inline int get_slot(DCELL c)
{
    int i;

    if (integer_mode)
        i = (int)(c - min);
    else
        i = (int)floor((c - min) / slot_size);

    if (i < 0)
        i = 0;
    if (i > num_slots - 1)
        i = num_slots - 1;
    return i;
}

/* The value a slot stands for: the exact integer in integer mode, the bin
 * center otherwise. */
static inline DCELL slot_value(int i)
{
    if (integer_mode)
        return min + i;
    return min + slot_size * (i + 0.5);
}

/* Scan the map for min and max when the range file is missing. */
static void scan_range(int infile, int rows, int cols)
{
    DCELL *inbuf = Rast_allocate_d_buf();
    int row, col;
    int found = 0;

    G_message(_("Computing range"));

    for (row = 0; row < rows; row++) {
        Rast_get_d_row(infile, inbuf, row);

        for (col = 0; col < cols; col++) {
            if (Rast_is_d_null_value(&inbuf[col]))
                continue;

            if (!found || inbuf[col] < min)
                min = inbuf[col];
            if (!found || inbuf[col] > max)
                max = inbuf[col];
            found = 1;
        }

        G_percent(row, rows, 2);
    }

    G_percent(rows, rows, 2);
    G_free(inbuf);

    if (!found)
        G_fatal_error(_("All cells in input map are NULL"));
}

static grass_int64 get_histogram(int infile, int rows, int cols,
                                 grass_int64 *hist)
{
    DCELL *inbuf = Rast_allocate_d_buf();
    grass_int64 total = 0;
    int row, col;

    G_message(_("Computing histogram"));

    for (row = 0; row < rows; row++) {
        Rast_get_d_row(infile, inbuf, row);

        for (col = 0; col < cols; col++) {
            if (Rast_is_d_null_value(&inbuf[col]))
                continue;

            hist[get_slot(inbuf[col])]++;
            total++;
        }

        G_percent(row, rows, 2);
    }

    G_percent(rows, rows, 2);
    G_free(inbuf);

    return total;
}

/* Weighted analogue of AS_class_quant: break i is the value whose
 * cumulative cell count reaches rank total * (i + 1) / nclass. */
static void quant_breaks(const grass_int64 *hist, grass_int64 total, int nclass,
                         int nbreaks, double classbreaks[])
{
    double cum = 0.0;
    int slot = 0;
    int i;

    for (i = 0; i < nbreaks; i++) {
        double rank = (double)total * (i + 1) / nclass;

        while (slot < num_slots - 1 && cum + hist[slot] < rank)
            cum += hist[slot++];
        classbreaks[i] = slot_value(slot);
    }
}

/* Count cells per class; slots and breaks are both ascending, values equal
 * to a break belong to the class below it. */
static void class_frequencies(const grass_int64 *hist, int nbreaks,
                              const double classbreaks[],
                              grass_int64 *frequencies)
{
    int slot, j = 0;

    for (slot = 0; slot < num_slots; slot++) {
        if (hist[slot] == 0)
            continue;

        while (j < nbreaks && slot_value(slot) > classbreaks[j])
            j++;
        frequencies[j] += hist[slot];
    }
}

int main(int argc, char *argv[])
{
    struct GModule *module;
    struct {
        struct Option *input, *nclasses, *algorithm, *bins, *format, *nprocs;
    } opt;
    struct {
        struct Flag *b, *r;
    } flag;
    struct FPRange range;
    enum OutputFormat format;
    char *desc;
    int infile, rows, cols;
    int nclass, nbreaks, num_bins, algorithm;
    int have_range, i;
    grass_int64 *hist, *frequencies, total;
    double *values, *weights, *classbreaks;
    double finfo = 1.0;

    G_gisinit(argv[0]);

    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("statistics"));
    G_add_keyword(_("classification"));
    G_add_keyword(_("parallel"));
    module->label = _("Computes class breaks for a raster map.");
    module->description =
        _("Classifies the range of cell values into classes using Jenks "
          "natural breaks, equal intervals, or quantiles, and prints the "
          "class breaks.");

    opt.input = G_define_standard_option(G_OPT_R_INPUT);

    opt.nclasses = G_define_option();
    opt.nclasses->key = "nclasses";
    opt.nclasses->type = TYPE_INTEGER;
    opt.nclasses->required = YES;
    opt.nclasses->options = "2-";
    opt.nclasses->description = _("Number of classes to define");

    opt.algorithm = G_define_option();
    opt.algorithm->key = "algorithm";
    opt.algorithm->type = TYPE_STRING;
    opt.algorithm->required = NO;
    opt.algorithm->options = "jen,int,qua";
    opt.algorithm->answer = "jen";
    opt.algorithm->description = _("Algorithm to use for classification");
    desc = NULL;
    G_asprintf(&desc,
               "jen;%s;"
               "int;%s;"
               "qua;%s",
               _("Jenks natural breaks"), _("simple intervals"),
               _("quantiles"));
    opt.algorithm->descriptions = desc;

    opt.bins = G_define_option();
    opt.bins->key = "bins";
    opt.bins->type = TYPE_INTEGER;
    opt.bins->required = NO;
    opt.bins->options = "2-";
    opt.bins->answer = "10000";
    opt.bins->description =
        _("Number of histogram bins used to approximate the values of "
          "floating-point maps");

    opt.format = G_define_standard_option(G_OPT_F_FORMAT);
    opt.format->options = "plain,csv,json";
    opt.format->descriptions = ("plain;Human readable text output;"
                                "csv;CSV (Comma Separated Values);"
                                "json;JSON (JavaScript Object Notation);");
    opt.format->guisection = _("Print");

    opt.nprocs = G_define_standard_option(G_OPT_M_NPROCS);

    flag.b = G_define_flag();
    flag.b->key = 'b';
    flag.b->description = _("Print only class breaks (without min and max)");

    flag.r = G_define_flag();
    flag.r->key = 'r';
    flag.r->description =
        _("Print recode rules for r.recode based on the class intervals");

    G_option_exclusive(flag.b, flag.r, NULL);

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    if (G_set_omp_num_threads(opt.nprocs) < 1)
        G_fatal_error(_("<%s> is not valid number of nprocs."),
                      opt.nprocs->answer);

    if (strcmp(opt.format->answer, "json") == 0)
        format = JSON;
    else if (strcmp(opt.format->answer, "csv") == 0)
        format = CSV;
    else
        format = PLAIN;

    nclass = atoi(opt.nclasses->answer);
    nbreaks = nclass - 1;
    num_bins = atoi(opt.bins->answer);

    infile = Rast_open_old(opt.input->answer, "");

    rows = Rast_window_rows();
    cols = Rast_window_cols();

    have_range = Rast_read_fp_range(opt.input->answer, "", &range) >= 0;
    if (have_range)
        Rast_get_fp_range_min_max(&range, &min, &max);
    if (!have_range || Rast_is_d_null_value(&min) || Rast_is_d_null_value(&max))
        scan_range(infile, rows, cols);

    /* For an integer map whose range fits into the requested number of
     * bins, one bin per integer value makes the breaks exact. */
    integer_mode = Rast_get_map_type(infile) == CELL_TYPE &&
                   max - min + 1 <= (double)num_bins;
    if (integer_mode) {
        num_slots = (int)(max - min) + 1;
        slot_size = 1.0;
    }
    else {
        num_slots = num_bins;
        slot_size = (max - min) / num_slots;
        if (slot_size <= 0.0) {
            /* constant map: a single slot holds everything */
            num_slots = 1;
            slot_size = 1.0;
        }
    }

    hist = G_calloc(num_slots, sizeof(grass_int64));
    total = get_histogram(infile, rows, cols, hist);
    Rast_close(infile);

    if (total == 0)
        G_fatal_error(_("All cells in input map are NULL"));

    classbreaks = G_calloc(nclass, sizeof(double));

    algorithm = AS_option_to_algorithm(opt.algorithm);
    switch (algorithm) {
    case CLASS_JENKS: {
        values = G_malloc(num_slots * sizeof(double));
        weights = G_malloc(num_slots * sizeof(double));
        for (i = 0; i < num_slots; i++) {
            values[i] = slot_value(i);
            weights[i] = (double)hist[i];
        }
        finfo = AS_class_jenks_weighted(values, weights, num_slots, &nbreaks,
                                        classbreaks);
        G_free(values);
        G_free(weights);
        if (finfo == 0.0)
            G_fatal_error(_("Classification algorithm failed"));
        break;
    }
    case CLASS_INTERVAL: {
        double minmax[2];

        minmax[0] = min;
        minmax[1] = max;
        AS_class_interval(minmax, 2, nbreaks, classbreaks);
        break;
    }
    case CLASS_QUANT:
        quant_breaks(hist, total, nclass, nbreaks, classbreaks);
        break;
    default:
        G_fatal_error(_("Unknown algorithm '%s'"), opt.algorithm->answer);
    }

    /* jen reduces the number of breaks when the map has fewer distinct
     * values than the requested classes */
    nclass = nbreaks + 1;

    frequencies = G_calloc(nclass, sizeof(grass_int64));
    class_frequencies(hist, nbreaks, classbreaks, frequencies);
    G_free(hist);

    if (flag.r->answer) {
        /* Highest class first: r.recode lets the last matching rule win,
         * so cells equal to a break get the class below it, consistent
         * with the frequencies reported by the other output formats. */
        for (i = nclass - 1; i >= 0; i--)
            fprintf(stdout, "%f:%f:%d\n", i ? classbreaks[i - 1] : min,
                    i < nbreaks ? classbreaks[i] : max, i + 1);
        exit(EXIT_SUCCESS);
    }

    if (flag.b->answer) {
        if (format == JSON) {
            G_JSON_Value *breaks_value = G_json_value_init_array();
            G_JSON_Array *breaks_array;
            char *json_string;

            if (breaks_value == NULL)
                G_fatal_error(
                    _("Failed to initialize JSON array. Out of memory?"));
            breaks_array = G_json_array(breaks_value);

            for (i = 0; i < nbreaks; i++)
                G_json_array_append_number(breaks_array, classbreaks[i]);

            json_string = G_json_serialize_to_string_pretty(breaks_value);
            if (!json_string) {
                G_json_value_free(breaks_value);
                G_fatal_error(_("Failed to serialize JSON to pretty format."));
            }
            puts(json_string);
            G_json_free_serialized_string(json_string);
            G_json_value_free(breaks_value);
        }
        else {
            for (i = 0; i < nbreaks; i++)
                fprintf(stdout, "%s%f", i ? "," : "", classbreaks[i]);
            fprintf(stdout, "\n");
        }
        exit(EXIT_SUCCESS);
    }

    switch (format) {
    case PLAIN:
        fprintf(stdout, _("\nClassification of <%s> into %i classes\n"),
                opt.input->answer, nclass);
        fprintf(stdout, _("Using algorithm: *** %s ***\n"),
                opt.algorithm->answer);
        if (algorithm == CLASS_JENKS)
            fprintf(stdout, _("Goodness of variance fit = %f\n"), finfo);
        fprintf(stdout, "\n");
        fprintf(stdout, _("%15s%15s%15s\n\n"), "From (excl.)", "To (incl.)",
                "Cell count");

        for (i = 0; i < nclass; i++)
            fprintf(
                stdout, "%15.5f%15.5f%15lld\n", i ? classbreaks[i - 1] : min,
                i < nbreaks ? classbreaks[i] : max, (long long)frequencies[i]);

        fprintf(stdout, _("\nNote: Minimum of first class is including\n\n"));
        break;

    case CSV:
        fprintf(stdout, "from,to,cell_count\n");
        for (i = 0; i < nclass; i++)
            fprintf(stdout, "%.5f,%.5f,%lld\n", i ? classbreaks[i - 1] : min,
                    i < nbreaks ? classbreaks[i] : max,
                    (long long)frequencies[i]);
        break;

    case JSON: {
        G_JSON_Value *root_value, *breaks_value, *intervals_value,
            *interval_value;
        G_JSON_Object *root_object, *interval_object;
        G_JSON_Array *breaks_array, *intervals_array;
        char *json_string;

        root_value = G_json_value_init_object();
        if (root_value == NULL)
            G_fatal_error(
                _("Failed to initialize JSON object. Out of memory?"));
        root_object = G_json_object(root_value);

        G_json_object_set_string(root_object, "algorithm",
                                 opt.algorithm->answer);
        G_json_object_set_number(root_object, "classes", nclass);
        G_json_object_set_number(root_object, "min", min);
        G_json_object_set_number(root_object, "max", max);
        if (algorithm == CLASS_JENKS)
            G_json_object_set_number(root_object, "goodness_of_variance_fit",
                                     finfo);

        breaks_value = G_json_value_init_array();
        if (breaks_value == NULL)
            G_fatal_error(_("Failed to initialize JSON array. Out of memory?"));
        breaks_array = G_json_array(breaks_value);

        intervals_value = G_json_value_init_array();
        if (intervals_value == NULL)
            G_fatal_error(_("Failed to initialize JSON array. Out of memory?"));
        intervals_array = G_json_array(intervals_value);

        for (i = 0; i < nclass; i++) {
            interval_value = G_json_value_init_object();
            if (interval_value == NULL)
                G_fatal_error(
                    _("Failed to initialize JSON object. Out of memory?"));
            interval_object = G_json_object(interval_value);

            G_json_object_set_number(interval_object, "from",
                                     i ? classbreaks[i - 1] : min);
            G_json_object_set_number(interval_object, "to",
                                     i < nbreaks ? classbreaks[i] : max);
            G_json_object_set_number(interval_object, "cell_count",
                                     (double)frequencies[i]);

            G_json_array_append_value(intervals_array, interval_value);

            if (i < nbreaks)
                G_json_array_append_number(breaks_array, classbreaks[i]);
        }

        G_json_object_set_value(root_object, "breaks", breaks_value);
        G_json_object_set_value(root_object, "intervals", intervals_value);

        json_string = G_json_serialize_to_string_pretty(root_value);
        if (!json_string) {
            G_json_value_free(root_value);
            G_fatal_error(_("Failed to serialize JSON to pretty format."));
        }
        puts(json_string);
        G_json_free_serialized_string(json_string);
        G_json_value_free(root_value);
        break;
    }
    }

    exit(EXIT_SUCCESS);
}
