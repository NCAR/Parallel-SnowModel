The "stations" directory provides examples of processing
met station data.  This directory also includes the MicroMet
preprocessor that can be used to fill in missing station data
(see the original MicroMet paper).

The different (re)analyses names (like: merra2, nldas2,
nora10, wrf404) provide examples of processing gridded atmospheric
(re)analyses data.

To download the (re)analysis data, edit and run
the "wget_reanalysis_met_gliston.script" in the
"/sm/met/data_download/" directory. This script has 
the capability to process era5, merra2, nldas2, and 
nora10 using Glen's ftp site. An additional script,
"wrf_404_download.py" can download WRF CONUS 404 data 
https://rda.ucar.edu/datasets/ds559.0/

