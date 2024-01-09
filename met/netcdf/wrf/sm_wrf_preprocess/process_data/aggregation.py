##aggregation.py

"""
LOAD PYTHON LIBRARIES
"""
import netCDF4 as nc
import pyproj
import numpy as np
import xarray as xr
import os
import sys

class Output_T_P(object):

    def __init__(self,xr_ds,varname,yr_mo_lst,Hourly,output_dir):

        self.ds = xr_ds
        self.varname = varname
        self.yr_mo_lst = yr_mo_lst
        self.Hourly = Hourly

        ## create output base directory ##
        if self.Hourly == 1:
            hourly_str = '1_hrly/'
        elif self.Hourly == 3:
            hourly_str = '3_hrly/'
        else:
            hourly_str = 'daily/'
        self.outputdir = output_dir + hourly_str + 'year/'
        print(self.outputdir)

        self.toFile()


    def toFile(self):
        # # loop through months and output
        for var in [self.varname]:
            for item in self.yr_mo_lst:
                ms = item[1]
                yr = item[0]
                print("%s-%s" % (yr,ms))
                da_month = self.ds[var].sel(Time=slice("%s-%s" % (yr,ms),"%s-%s" % (yr,ms)))
                da_month.to_dataset().to_netcdf(self.outputdir + yr + '/month/' + ms + '/' + var + '.nc',
                                                encoding = {"Time":
                                                            {'dtype':'float64',
                                                             'units': 'hours since 1901-01-01 00:00:00',
                                                             'calendar':'standard'}})

class Output_RELH(object):

    def __init__(self,xr_Q2,xr_PSFC,xr_T2,varname,yr_mo_lst,Hourly,output_dir):

        self.ds = xr_Q2
        self.ds_PSFC = xr_PSFC
        self.ds_T2 = xr_T2
        self.varname = varname
        self.yr_mo_lst = yr_mo_lst
        self.Hourly = Hourly

        ## create output base directory ##
        if self.Hourly == 1:
            hourly_str = '1_hrly/'
        elif self.Hourly == 3:
            hourly_str = '3_hrly/'
        else:
            hourly_str = 'daily/'
        self.outputdir = output_dir + hourly_str + 'year/'


        self.calcRELH()
        self.toFile()

    def calcRELH(self):
        es_T_0 = 6.11 # hPa # might need to multiply by 100!! because pressure is in Pa
        T_0 = 273.15 # K
        L = 2500000 # J/kg
        Rw = 461.52 # J/kg*K

        self.ds['PSFC'] = self.ds_PSFC['PSFC']
        self.ds['T2'] = self.ds_T2['T2']

        ## calculate partial water vapor pressure (e)
        self.ds['e'] = (self.ds['Q2'] * self.ds['PSFC']) / (0.622 + (0.378*self.ds['Q2']))

        ## calculate saturation water vapor pressure (es)
        self.ds['e_s'] = es_T_0 * np.exp((L/Rw)*((1/T_0)-(1/self.ds['T2'])))

        ## caculate relative humidity
        self.ds['RELH_raw'] = self.ds['e'] / self.ds['e_s']

        ## clean relative humidity by setting values over 100 to 100
        conditions  = [ (self.ds['RELH_raw']<= 100.0),(self.ds['RELH_raw']> 100.0) ]
        choices     = [ self.ds['RELH_raw'],100.0 ]
        dir_array = np.select(conditions, choices, default=np.nan)

        ## check for nan values in array
        if len(np.argwhere(np.isnan(dir_array))) > 0:
            print("NAN VALUE FOUND")
            print(np.argwhere(np.isnan(dir_array)))

        self.ds['RELH'] = (('Time','south_north','west_east'),dir_array)

    def toFile(self):
        # # loop through months and output
        for var in [self.varname]:
            for item in self.yr_mo_lst:
                ms = item[1]
                yr = item[0]
                print("%s-%s" % (yr,ms))
                da_month = self.ds[var].sel(Time=slice("%s-%s" % (yr,ms),"%s-%s" % (yr,ms)))
                da_month.to_dataset().to_netcdf(self.outputdir + yr + '/month/' + ms + '/' + var + '.nc',
                                                encoding = {"Time":
                                                            {'dtype':'float64',
                                                             'units': 'hours since 1901-01-01 00:00:00',
                                                             'calendar':'standard'}})

class Output_wind(object):

    def __init__(self,xr_U10,xr_V10,varname,yr_mo_lst,Hourly,output_dir):

        self.ds = xr_U10
        self.ds_V10 = xr_V10
        self.varname = varname
        self.yr_mo_lst = yr_mo_lst
        self.Hourly = Hourly

        ## create output base directory ##
        if self.Hourly == 1:
            hourly_str = '1_hrly/'
        elif self.Hourly == 3:
            hourly_str = '3_hrly/'
        else:
            hourly_str = 'daily/'
        self.outputdir = output_dir + hourly_str + 'year/'

        self.calcwind()
        self.toFile()

    def calcwind(self):
        ## add U10
        self.ds['V10'] = self.ds_V10['V10']

        ## calculate WSPD
        self.ds['WSPD']= (self.ds['V10']**2 + self.ds['U10']**2)**0.5

        ## calculate WDIR
        self.ds['V_U'] = abs(self.ds['V10']) / abs(self.ds['U10'])
        self.ds['WDEG'] = np.arctan(self.ds['V_U']) * 57.2958

## Remember direction is based on where wind is coming from. NOT WHERE IT IS COMING TO.
## A positive U and V has a wind-direction between 180 and 270.
        conditions  = [ ((self.ds['U10'] > 0.0) & (self.ds['V10'] > 0.0)), # Q1
                       ((self.ds['U10'] > 0.0) & (self.ds['V10'] < 0.0)),  # Q2
                       ((self.ds['U10'] < 0.0) & (self.ds['V10'] < 0.0)),  # Q3
                       ((self.ds['U10'] < 0.0) & (self.ds['V10'] > 0.0)), # Q4
                       ((self.ds['V10'] > 0.0) & (self.ds['U10'] == 0.0)), # 0
                       ((self.ds['V10'] == 0.0) & (self.ds['U10'] > 0.0)), # 90
                       ((self.ds['V10'] < 0.0) & (self.ds['U10'] == 0.0)),# 180
                       ((self.ds['V10'] == 0.0) & (self.ds['U10'] < 0.0)) # 270
                      ]

        choices     = [ 90.0 - self.ds['WDEG'] + 180.0, # Q1
                       90.0 + self.ds['WDEG'] + 180.0, # Q2
                       270.0 - self.ds['WDEG'] - 180.0, # Q3
                       270.0 + self.ds['WDEG'] - 180.0, # Q4
                      180.0,
                      270.0,
                      0.0,
                      90.0
                      ]
        dir_array = np.select(conditions, choices, default=np.nan)

        ## check for nan values in array
        if len(np.argwhere(np.isnan(dir_array))) > 0:
            print("NAN VALUE FOUND")
            print(np.argwhere(np.isnan(dir_array)))

        ## create 'wdir' array from dir_array
        self.ds['WDIR'] = (('Time','south_north','west_east'),dir_array)
        self.ds = self.ds.drop(['WDEG','V_U'])
        print('Finished processing wind')

    def toFile(self):
        # # loop through months and output
        for var in ['WDIR','WSPD']:
            for item in self.yr_mo_lst:
                print(var)
                ms = item[1]
                yr = item[0]
                print("%s-%s" % (yr,ms))
                da_month = self.ds[var].sel(Time=slice("%s-%s" % (yr,ms),"%s-%s" % (yr,ms)))
                da_month.to_dataset().to_netcdf(self.outputdir + yr + '/month/' + ms + '/' + var + '.nc',
                                                encoding = {"Time":
                                                            {'dtype':'float64',
                                                             'units': 'hours since 1901-01-01 00:00:00',
                                                             'calendar':'standard'}})
