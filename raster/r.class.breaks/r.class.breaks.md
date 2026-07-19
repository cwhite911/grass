## DESCRIPTION

*r.class.breaks* computes class breaks for the cell values of a raster
map, mainly for thematic mapping. Class breaks split the range of cell
values into **nclasses** classes; the tool prints the breaks together
with the number of cells falling into each class. The breaks can be
computed with the following algorithms (**algorithm**):

- `jen`: Jenks natural breaks (default). Finds the partition of values
  minimizing the sum of within-class variances (also known as
  Fisher-Jenks optimization or optimal univariate k-means). The result
  is the exact optimum, computed with a dynamic programming algorithm
  that runs in time proportional to the number of classes times
  n log n, parallelized with OpenMP (**nprocs**). The reported goodness
  of variance fit is `1 - SSD(classified) / SSD(unclassified)`; values
  close to 1 indicate that the classes explain most of the variance.
- `int`: equal intervals of size `(max - min) / nclasses`.
- `qua`: quantiles, classes with approximately equal numbers of cells.

The tool reads the cell values into a histogram in a single pass, so
memory consumption is independent of the raster size. For integer (CELL)
maps whose value range does not exceed **bins**, the histogram holds
every occurring value and the breaks are exact. Otherwise the values are
binned into **bins** equal-width bins over the value range and each
break is a bin center, i.e. accurate to `(max - min) / bins`; raise
**bins** for more precision.

Breaks are values from the map (up to binning precision): a cell belongs
to the class whose break is the smallest break greater than or equal to
the cell value. The minimum and maximum are not part of the breaks.

With the **-b** flag only the breaks are printed, for use in scripts or
as input to other tools. With the **-r** flag rules for *r.recode* are
printed, mapping each class interval to its class number.

## NOTES

The Jenks algorithm reduces the number of classes when the map has fewer
distinct values (or occupied bins) than the requested classes.

This tool complements *v.class*, which computes class breaks for vector
attribute data using the same classification library.

## EXAMPLES

Compute five natural breaks classes for the elevation map (North
Carolina sample dataset):

```sh
g.region raster=elevation
r.class.breaks input=elevation nclasses=5
```

Print only the breaks in JSON format:

```sh
r.class.breaks input=elevation nclasses=5 format=json -b
```

Reclassify the map into five natural breaks classes with *r.recode*:

```sh
r.class.breaks input=elevation nclasses=5 -r | r.recode input=elevation \
    output=elevation_classes rules=-
```

## SEE ALSO

*[d.vect.thematic](d.vect.thematic.md), [r.quantile](r.quantile.md),
[r.recode](r.recode.md), [r.stats.quantile](r.stats.quantile.md),
[v.class](v.class.md)*

## AUTHORS

Corey White, OpenPlains Inc.
