# user_inputs.py
"""
    GENERATES USER INPUTS FOR NETCDF PREPROCESSING
"""
"""
    LOAD PYTHON LIBRARIES
"""
import json

data = {}
data['Variable'] = []
data['Variable'].append({
    'netcdf_fmt_input': '2', # netCDF format
                             # 1 = WRF-10yr
                             # 2 = WRF-30yr

# Desired SnowModel timestep
    'tstep_input': '3',
                        # 1 = hourly
                        # 3 = 3-hourly
                        # 24 = daily

# SnowModel output variable #
    'output_var_input': 'T2',
                              # PREC_ACC_NC = precipitation [mm/hr]
                              # T2 = temperature [K]
                              # WSPD = wind speed [m/s]
                              # WDIR = wind direction [degrees 0-360]
                              # RELH = relative humidity [%]
                              # GLW = long-wave downward radiation flux [W/m^2]
                              # SWDOWN = shortwave downward radiation flux [W/m^2]

# SnowModel Water Year #
    'wateryr_input': '',
                             # 1980-2021 = WRF-10yr
                             # 2000-2021 = WRF-30yr

# flag to define domain extent for preprocessing
    'IsSpecifyDomain_input': '',
                                  # 0 = CONUS
                                  # 1 = subdomains [indicated by lat / lon inputs]

# absolute path for netcdf file ctrl file
    'netcdf_ctrl_fpath': '',

# absolute file path of base directory of downloaded for netcdf files
    'netcdf_input_dir': '',

# absolute file path base directory for saving processed netcdf files
    'netcdf_output_dir': '',

# absolute file path directory of snowmodel bounds used for clipped [REQUIRED if IsSpecifyDomain_input = 1]
    'snowmodel_bounds':'',

# proj string of forcing
    'proj_string':'+proj=lcc +lat_1=30.0 +lat_2=50.0 +lat_0=39.100006 +lon_0=-97.9 +a=6370000 +b=6370000 +type=crs' # conus 404

})

with open('./json/plotting_data.json', 'w') as outfile:
    json.dump(data, outfile)
