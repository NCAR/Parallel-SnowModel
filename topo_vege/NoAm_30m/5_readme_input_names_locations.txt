You will also need to point to where your topography
and vegetation input files are located and what they are
called. This is done by changing the file paths found within
the following.dat files:

TOPOGRAPHY:
input_files_info/topo_filepath.dat
input_files_info/topo_proj_string.dat

VEGETATION:
input_files_info/vege_filename.dat
input_files_info/vege_proj_string.dat

If you are using Glen's 1-arcsec topography and 30 m vegetation
data files, you will only need to modify the topo_filepath.dat
and vege_filename.dat, because you have likely put them in
directory locations that are different than what Glen uses.
You don't need to change the projection strings - these
already reflect the native projection of these datasets.

