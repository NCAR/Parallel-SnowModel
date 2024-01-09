#!/usr/bin/env python
""" 
Python script to download selected files from rda.ucar.edu.
After you save the file, don't forget to make it executable
i.e. - "chmod 755 <name_of_script>"
"""
import sys, os
from urllib.request import build_opener
import pandas as pd

"""
User inputs start ------------------------------------\
"""
## directory for data download ##
download_dir = '/glade/derecho/scratch/rossamower/snow/data/met/wrf_404/raw/'
## water year for data download
wateryear = 1981
"""
User input end ------------------------------------\
"""

"""
Download Invariants / CTRL file
"""

opener = build_opener()
## see if constants has been downloaded ##
wrf_constants_fname = 'wrfconstants_usgs404.nc'
base_dir = 'https://data.rda.ucar.edu/ds559.0/'
isExists = os.path.exists(download_dir + wrf_constants_fname)

print('Downloaded WRF Constants....')
if not isExists:
    filelist = [
        'https://data.rda.ucar.edu/ds559.0/INVARIANT/wrfconstants_usgs404.nc'
    ]

    for file in filelist:
        ofile = download_dir + os.path.basename(file)
        print(ofile)
        sys.stdout.write("downloading " + ofile + " ... ")
        sys.stdout.flush()
        infile = opener.open(file)
        outfile = open(ofile, "wb")
        outfile.write(infile.read())
        outfile.close()
        sys.stdout.write("done\n")
else:
    print('WRF Constants Already Downloaded')

"""
Download files for water year
"""
## create datelist ##
dates = pd.date_range(start = f'{wateryear-1}-10-01',
              end = f'{wateryear}-10-01',
              freq = '1H')
## remove last hour ##
dates = dates[0:-1]
## string replace ##
dates_str = [str(i).replace(' ','_') for i in dates]

## see if wateryear directory has been created ##
wy_dir = download_dir + f'WY{wateryear}/'
isExists = os.path.exists(wy_dir)


print(f'Downloaded Daily WRF Files for WY{wateryear}....')
if not isExists:
    os.makedirs(wy_dir)
https_prefix = f'https://data.rda.ucar.edu/ds559.0/wy{wateryear}/'
file_prefix = 'wrf2d_d01_'
filelist = [https_prefix + i[0:4] + i[5:7] + '/' + file_prefix + i  + '.nc' for i in dates_str]
print(filelist)

## check to see if files have already been downloaded ##

for file in filelist:
    ofile = wy_dir + os.path.basename(file)
    ## check to see if file is already downloaded ##
    isExists = os.path.exists(ofile)
    if not isExists:
        sys.stdout.write("downloading " + ofile + " ... ")
        sys.stdout.flush()
        infile = opener.open(file)
        outfile = open(ofile, "wb")
        outfile.write(infile.read())
        outfile.close()
        sys.stdout.write("done\n")