After creating the file containing your input information,
the topo and vege processing programs and scripts are run by
executing the file:

process_topo_vege.script

If this file does not have execute permission, do this:

chmod u+x process_topo_vege.script

To run this script, simply type the file name and hit enter:

process_topo_vege.script

The resulting topo_vege.gdat and
SM_domain_config_info_OUTPUTS.dat files are used later to
create the MicroMet met forcing file. This is done in the
../../met/ directory and sub-directories.

To see what your topography and vegetation inputs look like,
you can open GrADS and run the plot_topo_vege.gs plotting
script.

