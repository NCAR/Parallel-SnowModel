If you want to build a similar procedure for topography and
vegetation datasets other than those described here, you will
need to change the information in the files:

TOPOGRAPHY:
input_files_info/topo_filepath.dat
input_files_info/topo_proj_string.dat

VEGETATION:
input_files_info/vege_filename.dat
input_files_info/vege_proj_string.dat

Also, it is likely that your new vegetation dataset will have
different assigned vegetation class numbers than the NALCMS
dataset used here. This means you will have to change the
process for converting from your input vegetation class to
SnowModel vegetation classes. That conversion is done (and
will need to be modified) in this program:

/3_merge_topo_vege/merge_topo_vege.f

There is also some checking done in the above file that
makes sure you are converting the landcover classes from
the NALCMS classes to the SnowModel classes. Essentially,
the program stops running if it finds you are not using the
NALCMS dataset. To fix this, you will have to change those
checks or you will get an error message.

It is also important to note that the programs and scripts
in /topo1/ are VERY specific to Glen's 1 degree by 1 degree
dem data-block structure.  You will have to do something very
different for your topography dataset. If you have a single
DEM file, you may be able to do processing similar to what we
do for our NALCMS vegetation dataset (located in /2_vege/).
The 300m Global topography and land cover dataset, located
here: /topo_vege/Global_300m/, also provides a data processing
example for data in a different format than Glen's 1 degree
by 1 degree dem data-block structure data.

