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
import datetime



class Introduction(object):
    """
       This class outputs text to the terminal that explains the program, takes in
       user inputs, and performs checks on the inputs.
    """
    def __init__(self,input_dir,output_dir,netcdf_format = None,timestep = None,\
                 output_var = None, water_year = None):
        """
           Initializes intro object by printing user information and checking inputs
        """
        print('')
        print('--------------------------------------------------------------------------------------')
        print('This program is designed to preprocess WRF data into a format')
        print('that SnowModel expects. Below are the following assumptions about the')
        print('WRF data.')
        print('')
        print('--------------------------------------------------------------------------------------')

        self.timestep = timestep
        self.netcdf_format = netcdf_format
        self.output_var = output_var
        self.water_year = water_year
        self.acceptableOutputvars = ["PREC_ACC_NC","T2","WSPD","WDIR","RELH"]
        self.acceptable_tstep = ['1','3','24']
        self.output_directory = output_dir
        self.input_directory = input_dir

        self.wrf_format()
        self.printing()
        self.wrf_format_vars()
        print('Netcdf format :',self.netcdf_format)
        self.input_timestep()
        print('timestep :',self.timestep)
        self.output_variable()
        print('Output var :',self.output_var)
        self.water_yr(input_dir)
        print('Water year :',self.water_year)
        print('WRF file directory :',self.input_directory)
        self.input_fpath_check()

    def wrf_format(self):
        """
            Recieve input about netCDF format
        """
        input_check = True
        while input_check == True:
            input_check = False
            if self.netcdf_format is None:
                self.netcdf_format = input("Is the WRF data you are processing formatted like the NCAR 10-yr [1] data or NCAR 30-yr [2] data?\n")
                print('--------------------------------------------------------------------------------------')
                print('')
                try:
                    if int(self.netcdf_format) in [1,2]:
                        self.netcdf_format = int(self.netcdf_format)
                    else:
                        input_check = True
                        raise Exception("Invalid input: Please provide either 1 or 2 or create new formatting capability.\n")
                except Exception as e:
                    print(e)

    def printing(self):
        if self.netcdf_format == 1: # NCAR 10-yr
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
            print('     a). 10-yr WRF data [Water years: 2000 - 2013]')
            print('         - "wrf2d_d01_PGW_{met var}_{start year}{start month}-{start year}{start month}.nc"')
            print('4). NetCDF timestep')
            print('     a). 10-yr WRF data [Water years: 2000 - 2013]')
            print('         - hourly for three months')
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
        else: # NCAR 30-yr
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
            print('     a). 30-yr WRF data [Water years: 1980 - 2021]')
            print('         - "wrf2d_d01_{year}-{month}-{day}:{hour}:{minute}:{second}.nc"')
            print('4). NetCDF timestep')
            print('     a). 30-yr WRF data [Water years: 1980 - 2021]')
            print('         - each file represents one hour')
            print('5). The input netCDF files can be found in the following directory')
            print('     - "/glade/campaign/ncar/USGS_Water/CONUS404/"')
            print('     - within this directory there are year directories')
            print('     - Ex. "/glade/campaign/ncar/USGS_Water/CONUS404/WY2011"')
            print('6). Acceptable water years are 1980-2022')
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

    def wrf_format_vars(self):
        """
            Based on wrf_format method, set the appropriate variables
        """
        if self.netcdf_format == 1: # NCAR 10-yr WRF
            self.acceptablewy_int = list(range(2000,2014,1))
            self.acceptablewy_str = [str(i) for i in self.acceptablewy_int]
        else: # NCAR 30-yr WRF
            self.acceptablewy_int = list(range(1980,2023,1))
            self.acceptablewy_str = [str(i) for i in self.acceptablewy_int]

    def input_timestep(self):
        """
           Check input timestep.
        """
        input_check = True
        ## loop through until input error doesn't exist
        while input_check == True:
            input_check = False
            if self.timestep is None:
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
            if self.output_var is None:
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

    def water_yr(self,input_dir):
        """
           Check input timestep.
        """
        input_check = True
        ## loop through until input error doesn't exist
        while input_check == True:
            input_check = False
            if self.water_year is None:
                print('NOOOO')
                self.water_year = input("What water year would you like to output?\n")
                print('--------------------------------------------------------------------------------------')
                print('')
                try:
                    if (self.water_year in self.acceptablewy_str):
                        self.water_year = int(self.water_year)
                        if self.netcdf_format == 1: # NCAR 10-yr
                            self.input_directory = input_dir + 'CTRL/'
                        else: # NCAR 30-yr
                            self.input_directory = input_dir + 'WY' + str(self.water_year)
                        return
                    else:
                        input_check = True
                        if self.netcdf_format == 1: # NCAR 10-yr
                            raise Exception("Invalid input: Please provide either digits [2000 - 2013].\n")
                        else: # NCAR 30-yr
                            raise Exception("Invalid input: Please provide either digits [1980 - 2021].\n")
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



class Preprocess10yrWRF(object):
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
        print('--------------------------------------------------------------------------------------')
        print('')

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

class Preprocess30yrWRF(object):
    """
       This class outputs text to the terminal that explains the program, takes in
       user inputs, and performs checks on the inputs.
    """
    def __init__(self,varname,wateryear,Hourly,wrf_dir,output_dir,snowmodel_grid_dir,
                 lon_min_idx,lon_max_idx,lat_min_idx,lat_max_idx,netcdf_ctrl_fpath):
        """
           Initializes intro object by printing user information and checking inputs
        """
        ## input variables ##
        self.varname = varname
        self.wateryear = wateryear
        self.Hourly = Hourly
        self.wrf_dir = wrf_dir + 'WY' + str(self.wateryear) + '/'
        self.wrf_dir_prev = wrf_dir + 'WY' + str(self.wateryear-1) + '/'
        self.output_dir = output_dir
        self.sn_grid_dir = snowmodel_grid_dir
        self.lon_min_idx = lon_min_idx
        self.lon_max_idx = lon_max_idx
        self.lat_min_idx = lat_min_idx
        self.lat_max_idx = lat_max_idx
        self.netcdf_ctrl_fpath = netcdf_ctrl_fpath

        ## constants ##
        self.wrf_starting_name = 'wrf2d_d01_'
        self.datefmt = '%Y-%m-%d_%H:%M:%S'

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
        self.acceptablevars = ["PREC_ACC_NC","T2","U10","V10","Q2","PSFC","GLW","SWDOWN","SNOW_ACC_NC","GRAUPEL_ACC_NC"]
        # self.acceptablevars = ["GLW","SWDOWN"]
        self.acceptablewy = list(range(1980,2023,1))

        ## leap years ##
        self.leapyears = [1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]

        ## call check vars
        self.checkVars()

        ## month list

        if self.wateryear in self.leapyears:
            print('LEAP YEAR')
            print()
            self.month_lst = [ (str(self.wateryear-1),'08','31'),
                               (str(self.wateryear-1),'09','30'),
                               (str(self.wateryear-1),'10','31'),
                               (str(self.wateryear-1),'11','30'),
                               (str(self.wateryear-1),'12','31'),
                               (str(self.wateryear),'01','31'),
                               (str(self.wateryear),'02','29'),
                               (str(self.wateryear),'03','31'),
                               (str(self.wateryear),'04','30'),
                               (str(self.wateryear),'05','31'),
                               (str(self.wateryear),'06','30'),
                               (str(self.wateryear),'07','31'),
                               (str(self.wateryear),'08','31'),
                               (str(self.wateryear),'09','30')
                             ]
        else:
            print('NON WATER YEAR')

            self.month_lst = [  (str(self.wateryear-1),'08','31'),
                                (str(self.wateryear-1),'09','30'),
                                (str(self.wateryear-1),'10','31'),
                                (str(self.wateryear-1),'11','30'),
                                (str(self.wateryear-1),'12','31'),
                                (str(self.wateryear),'01','31'),
                                (str(self.wateryear),'02','28'),
                                (str(self.wateryear),'03','31'),
                                (str(self.wateryear),'04','30'),
                                (str(self.wateryear),'05','31'),
                                (str(self.wateryear),'06','30'),
                                (str(self.wateryear),'07','31'),
                                (str(self.wateryear),'08','31'),
                                (str(self.wateryear),'09','30')
                            ]

        self.createDirs()
        self.excludeVars()
        self.loopingOutput()

    def excludeVars(self):
        # open spatial dataset for lat / long
        lst = []
        for file in os.listdir(self.wrf_dir):
            if file.startswith('wrf2d_d01_'):
                lst.append(self.wrf_dir + file)
        spat_ds = xr.open_dataset(lst[0])
        # create list of variables to exclude #
        self.vars_to_exclude = list(spat_ds.keys())
        # print(self.vars_to_exclude)
        for var in self.acceptablevars:
            self.vars_to_exclude.remove(var)
        self.vars_to_exclude.remove('XTIME')

    def loopingOutput(self):
        ds_coords = xr.load_dataset(self.netcdf_ctrl_fpath)
        for i in range(1,len(self.month_lst)): # run for normal wateryears
        # for i in range(3,4): # testing purposes
            print('Working on year :',self.month_lst[i][0],'-',self.month_lst[i][1])
            mf_lst = []
            if i == 1: # september [need to grab files from previous wateryear directory]
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_21:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_22:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1]+ '-' + self.month_lst[i-1][2] + '_23:00:00')
                for file in os.listdir(self.wrf_dir_prev):
                    if file.startswith(self.wrf_starting_name +self.month_lst[i][0]+'-'+self.month_lst[i][1]):
                        mf_lst.append(file.replace(self.wrf_starting_name,''))
                lst_datetime_sort = sorted([datetime.datetime.strptime(u, self.datefmt) for u in mf_lst])
                lst_str_sort = [self.wrf_dir_prev + self.wrf_starting_name + u.strftime(self.datefmt) for u in lst_datetime_sort]

            elif i == 2: # october
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_21:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_22:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1]+ '-' + self.month_lst[i-1][2] + '_23:00:00')
                for file in os.listdir(self.wrf_dir):
                    if file.startswith(self.wrf_starting_name +self.month_lst[i][0]+'-'+self.month_lst[i][1]):
                        mf_lst.append(file.replace(self.wrf_starting_name,''))
                lst_datetime_sort = sorted([datetime.datetime.strptime(u, self.datefmt) for u in mf_lst])
                lst_str_sort = []
                for date in range(0,len(lst_datetime_sort)):
                    if date < 3:
                        lst_str_sort.append(self.wrf_dir_prev + self.wrf_starting_name + lst_datetime_sort[date].strftime(self.datefmt))
                    else:
                        lst_str_sort.append(self.wrf_dir + self.wrf_starting_name + lst_datetime_sort[date].strftime(self.datefmt))


            else: # other months
                # grab last 2 timesteps from previous month #
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_21:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1] + '-' + self.month_lst[i-1][2] + '_22:00:00')
                mf_lst.append(self.month_lst[i-1][0] + '-' + self.month_lst[i-1][1]+ '-' + self.month_lst[i-1][2] + '_23:00:00')
                for file in os.listdir(self.wrf_dir):
                    if file.startswith(self.wrf_starting_name +self.month_lst[i][0]+'-'+self.month_lst[i][1]):
                        mf_lst.append(file.replace(self.wrf_starting_name,''))

                lst_datetime_sort = sorted([datetime.datetime.strptime(u, self.datefmt) for u in mf_lst])
                lst_str_sort = [self.wrf_dir + self.wrf_starting_name + u.strftime(self.datefmt) for u in lst_datetime_sort]

            combined_ds = xr.open_mfdataset(lst_str_sort,concat_dim = 'Time',combine = 'nested',drop_variables = self.vars_to_exclude)

            print('Finished loading netcdf files')


            combined_ds = combined_ds.swap_dims({'Time': 'XTIME'})
            combined_ds = combined_ds.rename_dims({'XTIME':'Time'})
            combined_ds = combined_ds.rename({'XTIME':'Time'})

            combined_ds = combined_ds.assign_coords({'XLONG':ds_coords.XLONG[0,:,:],
                                 'XLAT':ds_coords.XLAT[0,:,:]})

            # combined_ds = combined_ds.drop('XTIME')

            # Process precipitation and snow #
            combined_ds['RAIN_ACC'] = combined_ds['PREC_ACC_NC'] - (combined_ds['SNOW_ACC_NC'] + combined_ds['GRAUPEL_ACC_NC'])
            combined_ds['SNOW_ACC'] = combined_ds['SNOW_ACC_NC'] + combined_ds['GRAUPEL_ACC_NC']
            
            # convert negative values to zero #
            combined_ds['RAIN_ACC'] = combined_ds.RAIN_ACC.where(combined_ds.RAIN_ACC > 0.0, 0.0)
            combined_ds['PREC_ACC_NC'] = combined_ds.PREC_ACC_NC.where(combined_ds.PREC_ACC_NC > 0.0, 0.0)
            combined_ds['Q2'] = combined_ds.Q2.where(combined_ds.Q2 > 0.0, 0.0)


            combined_ds = combined_ds.drop(['XTIME','SNOW_ACC_NC','GRAUPEL_ACC_NC'])

            # removing unnecessary starting and ending index #
            start = combined_ds.Time[0].values
            end = combined_ds.Time[-3].values
            combined_ds = combined_ds.sel(Time = slice(start,end))

            for var in combined_ds.variables:
                if var not in ["Time","XLAT","XLONG"]:
                    print(var)
                    if var in ["GLW","SWDOWN"]:
                        # if i == 0: #october
                            # idx_arr = combined_ds[var][3::3,:,:]
                            # time_arr = combined_ds['Time'][3::3]
                        # else:
                        idx_arr = combined_ds[var][3::3,:,:]
                        time_arr = combined_ds['Time'][3::3]
                    else:
                        da = combined_ds[var]
                        ds_ = da.to_dataset(name = var)
                        resamp_time = ds_.resample(Time = "3H").mean(dim="Time")
                        if var in ['PREC_ACC_NC','RAIN_ACC','SNOW_ACC']:
                            resamp_ds_ = ds_.resample(Time = "3H",base = 21).sum(dim="Time")
                        else:
                            resamp_ds_ = ds_.resample(Time = "3H",base = 21).mean(dim="Time")
                        # if i == 0: # october
                            # idx_arr = resamp_ds_[var][1:-1,:,:].values
                        # else:
                        idx_arr = resamp_ds_[var][0:-1,:,:].values
                        time_arr = resamp_time['Time'][1:].values

                    ds__ = xr.Dataset(
                        data_vars = dict(
                            dataset = (["Time","south_north","east_west"],idx_arr),),
                                coords={"Time": ("Time",time_arr),
                                        "XLAT": (("south_north","east_west"),ds_coords.XLAT[0,:,:]),
                                        "XLONG": (("south_north","east_west"),ds_coords.XLONG[0,:,:])
                                        })
                    ds__ = ds__.rename({'dataset' : var})

                    ds__.to_netcdf(self.output_dir + self.hourlyStr + 'year/' + self.month_lst[i][0] + '/month/' + self.month_lst[i][1] + '/' + var + '.nc',
                                  encoding = {"Time":
                                  {'dtype':'float64',
                                   'units': 'hours since 1901-01-01 00:00:00',
                                   'calendar':'standard'}})




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
                # if self.wateryear in self.leapyears:
                #     print('WARNING: WATER YEAR IS A LEAP YEAR. MAKE SURE TO ADJUST DAYS OF MONTH FOR FEBRUARY!!!')
            else:
                check = False
    #
    # def getLatLong(self):
    #     ## check WRF variable input
    #     self.ds_ctrl = xr.open_dataset(self.netcdf_ctrl_fpath)
    #
    #
    #
    #
    # def concatenateNetCDF(self):
    #     """
    #         Concatenate netCDF files based on variable and make any daily or 3-hrly aggregation
    #     """
    #     print('--------------------------------------------------------------------------------------')
    #     print('Start Concat')
    #     print('')
    #
    #     ## create list of strings of dates in folder ##
    #     lst_str_unsort = []
    #     for file in os.listdir(self.wrf_dir):
    #         if file.startswith(self.wrf_starting_name):
    #             lst_str_unsort.append(file.replace(self.wrf_starting_name,''))
    #     ## convert string to datetime and sort ##
    #     lst_datetime_sort = sorted([datetime.datetime.strptime(i, self.datefmt) for i in lst_str_unsort])
    #     ## convert back to string and add netcdf_file to begginning ##
    #     lst_str_sort = [self.wrf_starting_name + i.strftime(self.datefmt) for i in lst_datetime_sort]
    #
    #     ## create list of variables to exclude when concatenating netcdf files ##
    #     ds = xr.open_dataset(self.wrf_dir+lst_str_sort[0])
    #     vars_to_exclude = list(ds.keys())
    #     vars_to_exclude.remove(self.varname)
    #     vars_to_exclude.remove('XTIME')
    #     combined_ds = xr.open_mfdataset(self.wrf_dir + self.wrf_starting_name + '*',concat_dim = 'Time',combine = 'nested',drop_variables = vars_to_exclude)
    #
    #     print('Finished Concat')
    #     print('--------------------------------------------------------------------------------------')
    #     print('')
    #
    #     print('--------------------------------------------------------------------------------------')
    #     print('Start clean dataset')
    #     print('')
    #
    #     ## swap dimensions ##
    #     combined_ds = combined_ds.swap_dims({'Time': 'XTIME'})
    #     ## rename dimensions ##
    #     combined_ds = combined_ds.rename_dims({'XTIME':'Time'})
    #     # rename coords #
    #     combined_ds = combined_ds.rename({'XTIME':'Time'})
    #     # assign lat/long coords #
    #     combined_ds = combined_ds.assign_coords({'XLONG':self.ds_ctrl.XLONG[0,:,:],
    #                                          'XLAT':self.ds_ctrl.XLAT[0,:,:]})
    #     # drop XTIME variable #
    #     combined_ds = combined_ds.drop('XTIME')
    #
    #     print('Finished clean dataset')
    #     print('--------------------------------------------------------------------------------------')
    #     print('')
    #
    #     print('--------------------------------------------------------------------------------------')
    #     print('Start resample')
    #     print('')
    #
    #
    #     if self.Hourly == 3: # 3-hrly
    #         if self.varname == "PREC_ACC_NC":
    #             resample = combined_ds.resample(Time="3H").sum(dim="Time") # sum precip
    #             ## reindex array
    #             idx_arr = resample[self.varname][:,:,:].values
    #             resample[self.varname][1:,:,:] = idx_arr[:-1,:,:]
    #             # set first index = 0
    #             resample[self.varname][0,:,:] = 0.0
    #         else:
    #             resample = combined_ds.resample(Time="3H").mean(dim="Time")
    #             ## reindex array
    #             idx_arr = resample[self.varname][:,:,:].values
    #             resample[self.varname][1:,:,:] = idx_arr[:-1,:,:]
    #             # set first index = 0
    #             resample[self.varname][0,:,:] = -9999.0
    #         return resample
    #
    #     elif self.Hourly == 24: # daily
    #         if self.varname == "PREC_ACC_NC":
    #             resample = combined_ds.resample(Time="1D").sum(dim="Time") # sum precip
    #         else:
    #             resample = combined_ds.resample(Time="1D").mean(dim="Time")
    #         return resample
    #
    #     else: # hrly
    #         return combined_ds

class spatialIndex(object):

    def __init__(self,netcdf_ctr_dir,snowmodel_bounds,IsConus = None):
        self.netcdf_ctr_dir = netcdf_ctr_dir + 'INVARIANT/RALconus4km_wrf_constants.nc'
        self.sn_grid_dir = snowmodel_bounds
        self.IsConus = IsConus # flag to process data at conus scale or specified lat / lon bounds
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
        if self.IsConus == None:
            self.IsConus = input('Would you like to process WRF datasets for entire CONUS [0] or subdomains [1]\n')
            print('--------------------------------------------------------------------------------------')
            print('')
        try:
            if self.IsConus in [0,1]:
                pass
            else:
                raise Exception("Invalid designation for domain exten [needs to be 0 or 1].\nProgram will exit.\n")
        except Exception as e:
            print(e)
            sys.exit()

        if int(self.IsConus) == 1:
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
