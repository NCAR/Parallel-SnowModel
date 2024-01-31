# user_inputs_load.py
"""
    LOADS USER INPUT FOR NETCDF PREPROCESSING
"""
"""
    LOAD PYTHON LIBRARIES
"""
import json

## load user provided inputs for JSON file ##
with open('./json/plotting_data.json') as json_file:
    data = json.load(json_file)
    for p in data['Variable']:
        netcdf_fmt_input = int(p['netcdf_fmt_input'])
        tstep_input = int(p['tstep_input'])
        output_var_input = str(p['output_var_input'])
        wateryr_input = int(p['wateryr_input'])
        IsSpecifyDomain_input = int(p['IsSpecifyDomain_input'])
        netcdf_ctrl_fpath = str(p['netcdf_ctrl_fpath'])
        netcdf_input_dir = str(p['netcdf_input_dir'])
        netcdf_output_dir = str(p['netcdf_output_dir'])
        snowmodel_bounds = str(p['snowmodel_bounds'])
        proj_string = str(p['proj_string'])
