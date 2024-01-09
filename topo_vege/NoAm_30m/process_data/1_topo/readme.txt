This extracts dem data from Glen's 30-m (1-arcsec) North
American dem archive, and creates the topo distribution
required for your SnowModel run.

It requires access to Glen's 1-arcsec North American dem
data archive.

It is extracting 1-arcsec dem data for all of North America,
and creating a SnowModel topo dataset.

To do that it accesses the dem data over your area of interest,
projects it to your simulation projection, grids it to your
resolution, and crops it to your simulation domain.


NOTE: 1_find_subdomain_ll.f must be compiled with gfortran, pgf90,
or pgf95 (pgf77 won't work because of the floor and ceiling calls).

