##input.py

"""
LOAD PYTHON LIBRARIES
"""
import netCDF4 as nc
import pyproj
import numpy as np
import xarray as xr
import os
import sys
import csv



class Introduction(object):
    """
       This class outputs text to the terminal that explains the program, takes in
       user inputs, and performs checks on the inputs.
    """
    def __init__(self,input_dir,output_dir):
        """
           Initializes intro object by printing user information and checking inputs
        """
        print('')
        print('--------------------------------------------------------------------------------------')
        print('This program is designed to preprocess WRF data into a format')
        print('that SnowModel expects. Below are the following assumptions about the')
        print('WRF data.')
        print('')
        print('1). WRF data comes in hourly timesteps from 00:00 to 23:00 each day')
        print('2). The name and units of relevant meterological variables include:')
        print('     - precipitation [mm/hr] --> PREC_ACC_NC')
        print('     - temperature [K] --> T2')
        print('     - u-component of wind [m/s] --> U10')
        print('     - v-component of wind [m/s] --> V10')
        print('     - wind speed [m/s] --> WSPD')
        print('     - wind direction [degrees 0-360] --> WDIR')
        print('     - specific humidity [kg/kg] --> Q2')
        print('     - SFC pressure [Pa] --> PSFC')
        print('     - relative humidity [%] --> RELH')
        print('3). The name of the input netcdf file follows the following format')
        print('     - "wrf2d_d01_PGW_{met var}_{start year}{start month}-{start year}{start month}.nc"')
        print('4). The input netCDF files contain variable data for three months')
        print('     - Ex. "wrf2d_d01_PGW_U10_201107-201109.nc"')
        print('     - contains data for July through September 2011')
        print('5). The input netCDF files can be found in the following directory')
        print('     - "/glade/collections/rda/data/ds612.0/PGW/"')
        print('     - within this directory there are year directories')
        print('     - Ex. "/glade/collections/rda/data/ds612.0/PGW/2011"')
        print('6). Acceptable water years are 2000-2013')
        print('7). Water year designation will create monthly netCDF files for each')
        print('    SnowModel input variable from July through December of the ')
        print('    following year.')
        print('     - Ex. Water year 2012 --> creates monthly netCDF files from 07/2011 - 12/2012')
        print('8). Aggregation ocurrs for daily and 3-hourly timesteps.')
        print('    Values for precipitation [PREC_ACC_NC] are summed')
        print('    Remaining SnowModel variables are averaged using mean')
        print('      - temperature [T2]')
        print('      - u-component of wind [U10]')
        print('      - v-component of wind [V10]')
        print('      - specific humidity [Q2]')
        print('      - SFC pressure [PSFC]')
        print('9). A description of how WDIR, WSPD, and RELH are calculated is provided')
        print('    in the README.md document.')
        print('--------------------------------------------------------------------------------------')
        print('')

        self.timestep = 0
        self.output_var = ''
        self.water_year = 0
        self.acceptableOutputvars = ["PREC_ACC_NC","T2","WSPD","WDIR","RELH"]
        self.acceptablewy_int = list(range(2000,2014,1))
        self.acceptable_tstep = ['1','3','24']
        self.acceptablewy_str = [str(i) for i in self.acceptablewy_int]
        self.input_directory = input_dir + 'CTRL/'
        self.output_directory = output_dir

        self.input_timestep()
        self.output_variable()
        self.water_yr()
        self.input_fpath_check()

    def input_timestep(self):
        """
           Check input timestep.
        """
        input_check = True
        ## loop through until input error doesn't exist
        while input_check == True:
            input_check = False
            self.timestep = input("Would you like to output hourly [1], 3-hourly [3], or daily [24] timesteps?\n")
            print('--------------------------------------------------------------------------------------')
            print('')
            try:
                if self.timestep in self.acceptable_tstep:
                    self.timestep = int(self.timestep)
                    return
                else:
                    input_check = True
                    raise Exception("Invalid input: Please provide either 1, 3, or 24.\n")
            except Exception as e:
                print(e)

    def output_variable(self):
        """
           Check output variable.
        """
        input_check = True
        ## loop through until input error doesn't exist
        while input_check == True:
            input_check = False
            self.output_var = input("Which output variable would you like to create? ['PREC_ACC_NC','T2','WSPD','WDIR','RELH']\n")
            print('--------------------------------------------------------------------------------------')
            print('')
            try:
                if self.output_var in self.acceptableOutputvars:
                    return
                else:
                    input_check = True
                    raise Exception("Invalid input: Please provide ['PREC_ACC_NC','T2','WSPD','WDIR','RELH'] in all CAPS\n")
            except Exception as e:
                print(e)

    def water_yr(self):
        """
           Check input timestep.
        """
        input_check = True
        ## loop through until input error doesn't exist
        while input_check == True:
            input_check = False
            self.water_year = input("What water year would you like to output?\n")
            print('--------------------------------------------------------------------------------------')
            print('')
            try:
                if (self.water_year in self.acceptablewy_str):
                    self.water_year = int(self.water_year)
                    return
                else:
                    input_check = True
                    raise Exception("Invalid input: Please provide either digits [2000 - 2013].\n")
            except Exception as e:
                print(e)

    def input_fpath_check(self):
        """
           Check input directory.
        """
        if (os.path.isdir(self.input_directory) == False):
            print('Input directory for netCDF files could not be found.\n Program will exit.')
            print('--------------------------------------------------------------------------------------')
            print('')
            sys.exit()
        if (os.path.isdir(self.output_directory) == False):
            print('Output directory for netCDF files could not be found.\n Program will exit.')
            print('--------------------------------------------------------------------------------------')
            print('')
            sys.exit()



class PreprocessNetCDF(object):
    """
       This class outputs text to the terminal that explains the program, takes in
       user inputs, and performs checks on the inputs.
    """
    def __init__(self,varname,wateryear,Hourly,wrf_dir,output_dir,snowmodel_grid_dir,
                 lon_min_idx,lon_max_idx,lat_min_idx,lat_max_idx):
        """
           Initializes intro object by printing user information and checking inputs
        """
        ## input variables ##
        self.varname = varname
        self.wateryear = wateryear
        self.Hourly = Hourly
        self.wrf_dir = wrf_dir + 'CTRL/'
        self.output_dir = output_dir
        self.sn_grid_dir = snowmodel_grid_dir
        self.lon_min_idx = lon_min_idx
        self.lon_max_idx = lon_max_idx
        self.lat_min_idx = lat_min_idx
        self.lat_max_idx = lat_max_idx

        ## constants ##
        self.wrf_starting_name = 'wrf2d_d01_CTRL_'

        ## placeholder variables ##
        self.month_lst = []
        self.lat_min = 0.0
        self.lat_max = 0.0
        self.lon_min = 0.0
        self.lon_max = 0.0

        ## convert hourly to string ##
        if self.Hourly == 1: ## 1-hrly
            self.hourlyStr = '1_hrly/'
        elif self.Hourly == 3: ## 3-hrly
            self.hourlyStr = '3_hrly/'
        else: ## daily
            self.hourlyStr = 'daily/'

        ## create index flag ##
        if ((self.lon_min_idx + self.lon_max_idx + self.lat_min_idx + self.lat_max_idx) == -39996):
            self.index_flag = False
        else:
            self.index_flag = True

        ## acceptable lists of input values
        self.acceptablevars = ["PREC_ACC_NC","T2","U10","V10","Q2","PSFC"]
        self.acceptablewy = list(range(2000,2014,1))

        ## call check vars
        self.checkVars()

        ## month list
        self.month_lst = [(str(self.wateryear-1),'07'),
                           (str(self.wateryear-1),'08'),
                           (str(self.wateryear-1),'09'),
                           (str(self.wateryear-1),'10'),
                           (str(self.wateryear-1),'11'),
                           (str(self.wateryear-1),'12'),
                           (str(self.wateryear),'01'),
                           (str(self.wateryear),'02'),
                           (str(self.wateryear),'03'),
                           (str(self.wateryear),'04'),
                           (str(self.wateryear),'05'),
                           (str(self.wateryear),'06'),
                           (str(self.wateryear),'07'),
                           (str(self.wateryear),'08'),
                           (str(self.wateryear),'09'),
                           (str(self.wateryear),'10'),
                           (str(self.wateryear),'11'),
                           (str(self.wateryear),'12')
                         ]

        self.createDirs()
        # self.spatialClip()

    def createDirs(self):
        """
            Checks to see if output directory has been made. If not, directory gets created.
        """
        for mo in self.month_lst:
            fpath = self.output_dir + self.hourlyStr + 'year/' + mo[0] + '/month/' + mo[1] + '/'
            if os.path.exists(fpath):
                pass
            else:
                os.makedirs(fpath)


    def checkVars(self):
        ## check WRF variable input
        check = True
        while check == True:
            if self.varname not in self.acceptablevars:
                self.varname = input('Invalid Input, Please indicate variable you would like to output.\n["PREC_ACC_NC","T2","U10","V10","Q2","PSFC"]')
            else:
                check = False
        print('')

        ## check water year input
        check = True
        while check == True:
            if self.wateryear not in self.acceptablewy:
                self.wateryear = int(input('Invalid Input, Please indicate water year you would like to output for the following options.\n [2000 - 2013]'))
            else:
                check = False


    def concatenateNetCDF(self):
        """
            Concatenate netCDF files based on variable and make any daily or 3-hrly aggregation
        """
        ds_lst = []
        for i in range(0,len(self.month_lst),3):
            fname = self.wrf_dir + self.month_lst[i][0] + '/' + self.wrf_starting_name + self.varname + '_'  + self.month_lst[i][0] + self.month_lst[i][1]+ '-' + self.month_lst[i][0] + self.month_lst[i+2][1] + '.nc'
            print(fname)
            ds = xr.open_dataset(fname)
            if self.index_flag == True: ## index
                ds = ds.sel(south_north=slice(self.lat_min_idx-3,self.lat_max_idx+4),
                            west_east=slice(self.lon_min_idx-3,self.lon_max_idx+4))
            ds_lst.append(ds)

        combined_ds = xr.concat(objs = ds_lst ,dim = 'Time') ## works
        print('Finished Concat')

        if self.Hourly == 3: # 3-hrly
            if self.varname == "PREC_ACC_NC":
                resample = combined_ds.resample(Time="3H").sum(dim="Time") # sum precip
                ## reindex array
                idx_arr = resample[self.varname][:,:,:].values
                resample[self.varname][1:,:,:] = idx_arr[:-1,:,:]
                # set first index = 0
                resample[self.varname][0,:,:] = 0.0
            else:
                resample = combined_ds.resample(Time="3H").mean(dim="Time")
                ## reindex array
                idx_arr = resample[self.varname][:,:,:].values
                resample[self.varname][1:,:,:] = idx_arr[:-1,:,:]
                # set first index = 0
                resample[self.varname][0,:,:] = -9999.0
            return resample

        elif self.Hourly == 24: # daily
            if self.varname == "PREC_ACC_NC":
                resample = combined_ds.resample(Time="1D").sum(dim="Time") # sum precip
            else:
                resample = combined_ds.resample(Time="1D").mean(dim="Time")
            return resample

        else: # hrly
            return combined_ds

class spatialIndex(object):

    def __init__(self,netcdf_ctr_dir,snowmodel_bounds):
        self.netcdf_ctr_dir = netcdf_ctr_dir + 'INVARIANT/RALconus4km_wrf_constants.nc'
        self.sn_grid_dir = snowmodel_bounds
        self.lat_min = 0.0
        self.lat_max = 0.0
        self.lon_max = 0.0
        self.lon_min = 0.0

        self.x_ind_min = -9999
        self.x_ind_max = -9999
        self.y_ind_min = -9999
        self.y_ind_max = -9999


        self.readBounds()
        # self.netcdfIndices()

    def readBounds(self):
        subdomain_input = input('Would you like to process WRF datasets for entire CONUS [0] or subdomains [1]\n')
        print('--------------------------------------------------------------------------------------')
        print('')
        if int(subdomain_input) == 1:
            topo_domain_bounds = input('Would you like to use snowmodel grid lat/lon from Topo_Vege preprocessing? [0 -- No; 1 -- Yes]\n')
            print('--------------------------------------------------------------------------------------')
            print('')
            if int(topo_domain_bounds) == 1: ## process lat and lon max/min
                lat_bounds = []
                lon_bounds = []
                with open(self.sn_grid_dir, newline = '') as games:
    	            game_reader = csv.reader(games, delimiter='\t')
    	            for game in game_reader:
                        str_line = str(game[0])
                        lat_bounds.append(float(str_line.split('  ')[1]))
                        lon_bounds.append(float(str_line.split('  ')[0]))

                self.lat_max = max(lat_bounds)
                self.lat_min = min(lat_bounds)
                self.lon_max = max(lon_bounds)
                self.lon_min = min(lon_bounds)

                self.netcdfIndices()
            else: ## input
                print('Functionality has not been implemented.\n Program will exit.')
                print('--------------------------------------------------------------------------------------')
                print('')
                sys.exit()
        else:
            pass

    def netcdfIndices(self):
        ds = xr.open_dataset(self.netcdf_ctr_dir)
        xlat = ds.XLAT
        xlon = ds.XLONG
        xlat_clip = xlat.where(xlat>=self.lat_min).where(xlat<=self.lat_max)
        xlon_clip = xlon.where(xlon>=self.lon_min).where(xlon<=self.lon_max)
        y_ind_list,x_ind_list = np.where(xlat_clip.where(xlon_clip>-999).values>-999)

        self.x_ind_min = np.min(x_ind_list)
        self.x_ind_max = np.max(x_ind_list)
        self.y_ind_min = np.min(y_ind_list)
        self.y_ind_max = np.max(y_ind_list)
