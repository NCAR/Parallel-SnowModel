## main.py
"""
LOAD PYTHON LIBRARIES
"""
import netCDF4 as nc
import pyproj
import numpy as np
import xarray as xr


"""
LOAD PROGRAM MODULES
"""
from load_data.input import PreprocessNetCDF as pr
from load_data.input import Introduction as intr
from load_data.input import spatialIndex
from process_data.aggregation import Output_T_P
from process_data.aggregation import Output_RELH
from process_data.aggregation import Output_wind


"""
USER INPUT ####################################################################
"""
netcdf_input_dir = '/glade/collections/rda/data/ds612.0/'
netcdf_output_dir = '/glade/scratch/rossamower/snow/snowmodel/andy_sim/tuolumne/tuolumne_100m_exdistr/met/wrf/zoomed/'
snowmodel_bounds = '/glade/scratch/rossamower/snow/snowmodel/andy_sim/tuolumne/tuolumne_100m_exdistr/topo_vege/NoAm_30m/process_data/1_topo/outputs/1_ll_sm_corners.dat'

"""
PROGRAM INPUTS ################################################################
"""

## call Introduction class to get inputs ##
pgrm_inputs = intr(netcdf_input_dir,netcdf_output_dir)

## call spatial clip ##
spat_index = spatialIndex(netcdf_input_dir,snowmodel_bounds)

## save class instances to variables
output_tstep = pgrm_inputs.timestep
output_var = pgrm_inputs.output_var
output_wy = pgrm_inputs.water_year

"""
NETCDF PREPROCESS ##############################################################
"""

if output_var == 'RELH':
    """Q2"""
    ## create instance
    Q2 = pr('Q2',output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
            spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    Q2_ds = Q2.concatenateNetCDF()

    """PSFC"""
    ## create instance
    PSFC = pr('PSFC',output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
              spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    PSFC_ds = PSFC.concatenateNetCDF()

    """T2"""
    ## create instance
    T2 = pr('T2',output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
            spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    T2_ds = T2.concatenateNetCDF()

    ## create year month list
    yr_mo_lst = Q2.month_lst
    ## write to directory
    Output_RELH(Q2_ds,PSFC_ds,T2_ds,output_var,yr_mo_lst,output_tstep,netcdf_output_dir)

elif output_var == 'WSPD' or output_var == 'WDIR':
    """U10"""
    ## create instance
    U10 = pr('U10',output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
            spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    U10_ds = U10.concatenateNetCDF()
    print('Finished with U10')

    """V10"""
    ## create instance
    V10 = pr('V10',output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
            spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    V10_ds = V10.concatenateNetCDF()
    print('Finished with V10')

    # create year month list
    yr_mo_lst = U10.month_lst
    ## write to directory
    Output_wind(U10_ds,V10_ds,output_var,yr_mo_lst,output_tstep,netcdf_output_dir)

else: ## remaining variables [T2, PREC_ACC_NC]
    """T2 OR PREC_ACC_NC"""
    ## create instance
    var_obj = pr(output_var,output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
                spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max)
    ## create dataset
    var_ds = var_obj.concatenateNetCDF()

    ## write to directory
    Output_T_P(var_ds,output_var,var_obj.month_lst,output_tstep,netcdf_output_dir)
