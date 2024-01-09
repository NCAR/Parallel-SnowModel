You must define four things related to your SnowModel
simulation domain:

1) Identify whether you are going to provide general lat-lon
coords and the SnowModel domain will be created from that
information (==1), or you are going to provide exact coords
of the domain boundaries in the SnowModel projection you have
defined in 4) below (==2).  This first line can be either
a 1 or 2, depending on which kinds of coordinates you are
providing.

2) The grid increment you want, in meters (like 30.0).
You must provide deltax and deltay.

3) The approximate lon-lat coordinates of the lower-left and
upper-right corners of your simulation domain (the flag in
the first line of this file is set to '1'.

*NOTE: You can get these using something like Google Earth by
pointing the cursor where you want and recording the corner
positions in decimal degrees, or ArcGIS by viewing the corner
coordinates of your domain bounding box in decimal degrees.

Alternatively, you can set the flag in the first line of this
file to be '2'.  In this case the lower-left and upper-right
corner coordinates you provide will be the exact coordinates
you want the boundaries of your SnowModel simulation domain
to be, and they will be in the projection defined in 4) below.

4) The projection you want (e.g., UTM, LAEA, etc.) in PROJ.4
string format. It must be an equal-area projection, so that
you can define a grid increment in meters (you cannot run
the model on a latitude-longitude grid; actually you can,
talk to Glen if you are intersted in this).  Three common
projections are typically used for SnowModel runs: UTM
(+proj=utm), Lambert Azimuthal Equal Area (+prog=laea),
and Albers Equal Area (+prog=aea).  I usually use UTM for
small domains (10s of km), LAEA for moderate sized domains
(100s of km), and AEA for continental domains (1000s of km;
e.g., Alaska, CONUS, No Am).

An example of the required file format for this
necessary input information is provided by the file:
SM_dxdy_cornerll_proj_INPUTS.dat. Your file must look exactly
like that file, from a formatting perspective (same lines,
same order, etc.).

For example:

======================================
1

100.0  100.0

-170.000   56.000
-130.000   71.000

+proj=aea +lat_1=55.0 +lat_2=65.0 +lat_0=50.0 +lon_0=-154.0 +x_0=0.0 +y_0=0.0 +datum=NAD83 +units=m +no_defs
======================================

======================================
2

100.0  100.0

-1710000.00  -835000.00
  435000.00   490000.00

+proj=aea +lat_1=55.0 +lat_2=65.0 +lat_0=50.0 +lon_0=-154.0 +x_0=0.0 +y_0=0.0 +datum=NAD83 +units=m +no_defs
======================================

NOTES (the lines should look like this; the spaces need
to be at least one space in between the values; that's the
only requirement; you can leave any comments at the end of
the file):

======================================
  1
SPACE
  deltax  deltay
SPACE
  lower_left_lon   lower_left_lat
  upper_right_lon  upper_right_lat
SPACE
  +proj description
======================================

======================================
  2
SPACE
  deltax  deltay
SPACE
  lower_left_x_coords   lower_left_y_coords
  upper_right_x_coords  upper_right_y_coords
SPACE
  +proj description
======================================

