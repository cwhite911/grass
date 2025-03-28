## DESCRIPTION

*r.fill.dir* filters and generates a depressionless elevation map and a
flow direction map from a given raster elevation map. The method adopted
to filter the elevation map and rectify it is based on the paper titled
"Extracting topographic structure from digital elevation model data for
geographic information system analysis" by S.K. Jenson and J.O. Domingue
(1988).

The procedure takes an elevation layer as input and initially fills all
the depressions with one pass across the layer. Next, the flow direction
algorithm tries to find a unique direction for each cell. If the
watershed program detects areas with pothholes, it delineates this area
from the rest of the area and once again the depressions are filled
using the neighborhood technique used by the flow direction routine. The
final output will be a depressionless elevation layer and a unique flow
direction layer.

This (D8) flow algorithm performs as follows: At each raster cell the
code determines the slope to each of the 8 surrounding cells and assigns
the flow direction to the highest slope out of the cell. If there is
more than one equal, non-zero slope then the code picks one direction
based on preferences that are hard-coded into the program. If the
highest slope is flat and in more than one direction then the code first
tries to select an alternative based on flow directions in the adjacent
cells. *r.fill.dir* iterates that process, effectively propagating flow
directions from areas where the directions are known into the area where
the flow direction cannot otherwise be resolved.

The **format** parameter is the type of format at which the user wishes
to create the flow **direction** map. The flow direction map can be
encoded in GRASS GIS aspect format, ANSWERS (Beasley et.al, 1982), or
AGNPS (Young et.al, 1985) format, so that it can be readily used as
input to other GRASS GIS modules or the aforementioned hydrological
models. The *grass* format gives the same category values as
*[r.slope.aspect](r.slope.aspect.md)* gives for aspect, i.e. angles in
degrees counter-clockwise from east in 45 degree increments. The *agnps*
format gives category values from 1-8, with 1 facing north and
increasing values in the clockwise direction. The *answers* format gives
category values from 0-360 degrees, with 0 (represented as 360) facing
east and values increasing in the counter-clockwise direction at 45
degree increments. In all cases, NULL (no data) values are used for
cells where direction cannot be determined.

In case of local problems, those unfilled areas can be stored
optionally. Each unfilled area in this maps is numbered. The **-f** flag
instructs the program to fill single-cell pits but otherwise to just
find the undrained areas and exit. With the **-f** flag set the program
writes an elevation map with just single-cell pits filled, a direction
map with unresolved problems and a map of the undrained areas that were
found but not filled. This option was included because filling DEMs was
often not the best way to solve a drainage problem. These options let
the user get a partially-fixed elevation map, identify the remaining
problems and fix the problems appropriately.

In some cases it may be necessary to run *r.fill.dir* repeatedly (using
output from one run as input to the next run) before all of problem
areas are filled.

The resulting depressionless elevation raster map can further be
processed to derive slopes and other attributes required by other
hydrological models.

As any GRASS GIS module, *r.fill.dir* respects the computational region
settings. Thus, the module can be used to generate a flow direction map
for any sub-area within the full raster map layer. Also, *r.fill.dir*
will take into account an active raster mask.

## NOTES

- The *r.fill.dir* module can be used not only to fill depression, but
  also to detect water bodies or potential water bodies based on the
  nature of the terrain and the digital elevation model used.
- Not all depressions are errors in digital elevation models. In fact,
  many are wetlands and as Jenkins and McCauley (2006) note careless use
  of depression filling may lead to unintended consequences such as loss
  of wetlands.
- Although many hydrological algorithms require depression filling,
  advanced algorithms such as those implemented in
  *[r.watershed](r.watershed.md)* and *[r.sim.water](r.sim.water.md)* do
  not require depressionless digital elevation model to work.
- The flow direction map can be visualized with
  *[d.rast.arrow](d.rast.arrow.md)*.

## EXAMPLES

Generic example: create a depressionless (sinkless) elevation map
*ansi.fill.elev* and a flow direction map *ansi.asp* for the type
"grass":

```sh
r.fill.dir input=ansi.elev output=ansi.fill.elev direction=ansi.asp
```

North Carolina sample dataset example: The LiDAR derived 1m elevation
map is sink-filled. The outcome are a depressionless elevation map, the
flow direction map and an error map.

```sh
# set computational region to elevation map
g.region raster=elev_lid792_1m -p
# generate depressionless DEM and related maps
r.fill.dir input=elev_lid792_1m output=elev_lid792_1m_filled \
           direction=elev_lid792_1m_dir areas=elev_lid792_1m_error

# generate elevation map of pixelwise differences to see obtained terrain alterations
r.mapcalc "elev_lid792_1m_diff = elev_lid792_1m_filled - elev_lid792_1m"
r.colors elev_lid792_1m_diff color=differences

# assess univariate statistics of differences
r.univar -e elev_lid792_1m_diff

# vectorize filled areas (here all fills are of positive value, see r.univar output)
r.mapcalc "elev_lid792_1m_fill_area = if(elev_lid792_1m_diff > 0.0, 1, null() )"
r.to.vect input=elev_lid792_1m_fill_area output=elev_lid792_1m_fill_area type=area

# generate shaded terrain for better visibility of results
r.relief input=elev_lid792_1m_filled output=elev_lid792_1m_filled_shade

d.mon wx0
d.shade shade=elev_lid792_1m_filled_shade color=elev_lid792_1m_filled
d.vect elev_lid792_1m_fill_area type=boundary color=red
```

![r.fill.dir example](r_fill_dir.png)  
*Figure: Sink-filled DEM (shown as shaded terrain) with areas of filling
shown as vector polygons*

## REFERENCES

- Beasley, D.B. and L.F. Huggins. 1982. ANSWERS (areal nonpoint source
  watershed environmental response simulation): User's manual. U.S.
  EPA-905/9-82-001, Chicago, IL, 54 p.
- Jenkins, D. G., and McCauley, L. A. 2006. GIS, SINKS, FILL, and
  disappearing wetlands: unintended consequences in algorithm
  development and use. In Proceedings of the 2006 ACM symposium on
  applied computing (pp. 277-282).
- Jenson, S.K., and J.O. Domingue. 1988. Extracting topographic
  structure from digital elevation model data for geographic information
  system analysis. Photogram. Engr. and Remote Sens. 54: 1593-1600.
- Young, R.A., C.A. Onstad, D.D. Bosch and W.P. Anderson. 1985.
  Agricultural nonpoint surface pollution models (AGNPS) I and II model
  documentation. St. Paul: Minn. Pollution control Agency and Washington
  D.C., USDA-Agricultural Research Service.

## SEE ALSO

*[d.rast.arrow](d.rast.arrow.md), [d.shade](d.shade.md),
[g.region](g.region.md), [r.fillnulls](r.fillnulls.md),
[r.relief](r.relief.md), [r.slope.aspect](r.slope.aspect.md)*

## AUTHORS

Fortran version: Raghavan Srinivasan, Agricultural Engineering
Department, Purdue University  
Rewrite to C with enhancements: Roger S. Miller
