## main.py
"""
LOAD PYTHON LIBRARIES
"""
import netCDF4 as nc
import pyproj
import numpy as np
import xarray as xr
import sys
import os


"""
LOAD PROGRAM MODULES
"""
from load_data.user_inputs import *
from load_data.user_inputs_load import *
from load_data.input import Preprocess10yrWRF as pr10yr
from load_data.input import Preprocess30yrWRF as pr30yr
from load_data.input import Introduction as intr
from load_data.input import spatialIndex
from process_data.aggregation import Output_T_P
from process_data.aggregation import Output_RELH
from process_data.aggregation import Output_wind

"""
PROGRAM INPUTS ################################################################
"""

# call Introduction class to get inputs ##
pgrm_inputs = intr(netcdf_input_dir,netcdf_output_dir,netcdf_fmt_input,\
                       tstep_input,output_var_input,wateryr_input)
# call spatial clip ##
spat_index = spatialIndex(netcdf_input_dir,snowmodel_bounds,IsSpecifyDomain_input)


# save class instances to variables
output_tstep = pgrm_inputs.timestep
output_var = pgrm_inputs.output_var
output_wy = pgrm_inputs.water_year
wrf_fmt_flag = pgrm_inputs.netcdf_format

print('ctrl_fpath :',netcdf_ctrl_fpath)

"""
NETCDF PREPROCESS ##############################################################
"""
"""
    IF MONTHLY MET VARIABLES HAVEN'T BEEN CREATED. CREATE THOSE. OTHERWISE,
    CALCULATE MET VARIABLES FOR SNOWMODEL INPUT.
"""
if output_tstep == 1: ## 1-hrly
    hourlyStr = '1_hrly/'
elif output_tstep == 3: ## 3-hrly
    hourlyStr = '3_hrly/'
else: ## daily
    hourlyStr = 'daily/'

month_lst = [      (str(wateryr_input-1),'09'),
                   (str(wateryr_input-1),'10'),
                   (str(wateryr_input-1),'11'),
                   (str(wateryr_input-1),'12'),
                   (str(wateryr_input),'01'),
                   (str(wateryr_input),'02'),
                   (str(wateryr_input),'03'),
                   (str(wateryr_input),'04'),
                   (str(wateryr_input),'05'),
                   (str(wateryr_input),'06'),
                   (str(wateryr_input),'07'),
                   (str(wateryr_input),'08'),
                   (str(wateryr_input),'09')
                 ]

fpath = netcdf_output_dir + hourlyStr + 'year/' + str(wateryr_input) + '/month/01/'
if os.path.exists(fpath):
    print('NO')
    for i in range(0,len(month_lst)):
        fpath = netcdf_output_dir + hourlyStr + 'year/' + month_lst[i][0] + '/month/' + month_lst[i][1]
        print(fpath)
        # if output_var == 'RELH':
        print('--RELH')
        """Q2"""
        ## create instance
        Q2_ds = xr.load_dataset(fpath + '/Q2.nc')
        print('Load Q2 finished')

        """PSFC"""
        ## create instance
        PSFC_ds = xr.load_dataset(fpath + '/PSFC.nc')
        print('Load PSFC finished')

        """T2"""
        ## create instance
        T2_ds = xr.load_dataset(fpath + '/T2.nc')
        print('Load T2 finished')
        Output_RELH(Q2_ds,PSFC_ds,T2_ds,'RELH',month_lst,output_tstep,fpath)
        print('RELH output finished')
        print()

        # elif output_var == 'WSPD' or output_var == 'WDIR':
        print('--WSPD')
        print('------WDIR')
        """U10"""
        ## create instance

        U10_ds = xr.load_dataset(fpath + '/U10.nc')
        print('Finished with U10')

        """V10"""
        ## create instance

        V10_ds = xr.load_dataset(fpath + '/V10.nc')
        print('Finished with V10')

        Output_wind(U10_ds,V10_ds,'WSPD',month_lst,output_tstep,fpath,netcdf_ctrl_fpath)

        # else: ## remaining variables [T2, PREC_ACC_NC]
        #     print('------' + output_var)
        # ## create instance
        #
        #     var_ds = xr.load_dataset(fpath + '/' + output_var + '.nc')
        #
        #
        #     ## write to directory
        #     Output_T_P(var_ds,output_var,month_lst,output_tstep,fpath)




else:
    print('YES')
    pr30yr(output_var,output_wy,output_tstep,netcdf_input_dir,netcdf_output_dir,snowmodel_bounds,
                spat_index.x_ind_min,spat_index.x_ind_max,spat_index.y_ind_min,spat_index.y_ind_max,netcdf_ctrl_fpath)
    print('PRINT PROCESS INPUT VARIABLES. RE-RUN TO GET SNOWMODEL INPUT VARIABLES')

## end program
