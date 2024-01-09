DEM Data:

This is Glen's 1-arcsec, North American, dem dataset.
These are the US NED data from their 3DEP program. The topo
data came from here:

https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/1/TIFF/

I believe this is the latest USGS NED data.

This dataset consists of 1 degree, 1 arcsec (~30m) topo tiles
over all of North America.

There is one Alaska tile missing from the 1 arcsec archive.
Download the 2 arcsec file and regrid that to the 1 arcsec
grid.

There are 7 tiles in Canada that are missing.  To fill those
gaps, I downloaded data from the "Canadian Digital Elevation
Model, 1945-2011" dataset, and regridded and reprojected those
data to create new tiles that fill the missing tiles spaces.

Info on that dataset is available here:

https://open.canada.ca/data/en/dataset/7f245e4d-76c2-4caa-951a-45d1d2051333

https://ftp.maps.canada.ca/pub/nrcan_rncan/elevation/cdem_mnec/doc/CDEM_product_specs.pdf

https://open.canada.ca/data/en/fgpv_vpgf/7f245e4d-76c2-4caa-951a-45d1d2051333

The data are available here:

https://ftp.maps.canada.ca/pub/nrcan_rncan/elevation/cdem_mnec/

All of the processing to create this archive was done in
Glen's directory, here:

/data1/working/NoAm_30m_topo/


LANDCOVER Data:

This is the 2015 North American Land Change Monitoring System
(NALCMS) landcover dataset. The vege data has a documents
directory that says what it is and where it came from:

http://www.cec.org/north-american-land-change-monitoring-system/

http://www.cec.org/north-american-environmental-atlas/land-cover-30m-2015-landsat-and-rapideye/

http://www.cec.org/wp-content/uploads/wpallimport/files/Atlas/Files/2010nalcms30m/north_america_2015.zip

