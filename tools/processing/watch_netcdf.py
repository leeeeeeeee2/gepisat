#!/usr/bin/python
#
# watch_netcdf.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-08-07 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script reads a netCDF file (WATCH Forcing Data) and extracts data for
# the variable of interest for each grid point. The data extraction is then
# formatted for ASCII raster file output. 
#
# NOTE: WFDEI contains 67,209 pixels of terrestrial *observations*
#       Plus an additional 27,533 pixels over Antarctica
#       For a total of 94,742 pixels
#
# ----------
# changelog:
# ----------
# 00. created [13.08.07]
# 01. added lon and lat hashes [13.08.09]
# 02. added lon, lat, and tstep iterator [13.08.09]
# 03. added raster output function [13.10.22]
# 04. added necessary modules to import [13.10.22]
# 05. finished process_raster function [13.10.23]
# --> SWdown is 1000x's units
# 06. general housekeeping [14.12.01]
#
# -----------
# references:
# ----------- 
# 1. http://gfesuite.noaa.gov/developer/netCDFPythonInterface.html
# 2. http://www.gisintersect.com/?p=186
#
###############################################################################
## IMPORT MODULES
###############################################################################
import os
import datetime
import glob
import re
#import numpy
from scipy.io import netcdf

###############################################################################
## FUNCTIONS
###############################################################################
def writeout(f, d):
    """
    Name:     writeout
    Input:    - str, file name with path (t)
              - str, data to be written to file (d)
    Output:   None.
    Features: Writes new/overwrites existing file with data string
    """
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

def process_raster(my_dir, yr, mo, d):
    """
    Name:     process_raster
    Input:    - str, output file directory (my_dir)
              - int, year (yr)
              - int, month (mo)
              - numpy.ndarray, daily WATCH SWdown data (d)
    Output:   None.
    Features: Saves daily WATCH shortwave downwelling radiation data in ASCII 
              raster format (int value x1000)
    Depends:  writeout
    """
    # Define header line for ASCII raster:
    header = (
        "NCOLS 720\n"
        "NROWS 360\n"
        "XLLCORNER -180.0\n"
        "YLLCORNER -90.0\n"
        "CELLSIZE 0.5\n"
        "NODATA_VALUE -9999\n"
        )
    #
    # Save the shape values of each for iteration purposes
    tstep, lat, lon = d.shape
    #
    # Iterate through time:
    for t in xrange(tstep):
        # Save/create time stamp:
        this_day = t + 1
        time_stamp = datetime.date(yr, mo, this_day)
        #
        # Create output file name:
        rast_out = "%s%s_%s.txt" % (my_dir, "WATCH_0.5-Raster", time_stamp)
        #
        # Save header line to file:
        writeout(rast_out, header)
        #
        # Latitude and Longitude have two decimal places.
        for i in xrange(lat):
            # Reverse column ordering:
            y = 359 - i
            #
            # Reset outline for next row:
            outline = ""
            #
            for x in xrange(lon):
                # SWdown has five decimal places. 
                # * error values are ~1e20
                # ** for raster out, multiply by 1000
                # *** set error/missing value to -9999
                swd = data[t,y,x]
                if swd < 1.0e6:
                    swd = int(1e3*swd)
                else:
                    swd = -9999
                #
                # Add value to outline:
                outline = "%s%d " % (outline, swd)
            # End line and print to file:
            outline = "%s\n" % outline.rstrip(' ')
            OUT = open(rast_out, 'a')
            OUT.write(outline)
            OUT.close()

###############################################################################
## DEFINE CONSTANTS
###############################################################################
file_directory = (
    "/Users/twdavis/Projects/data/watch/SWdown_daily_WFDEI/"
    )
var_voi = 'SWdown'

###############################################################################
## MAIN
###############################################################################
# Read netCDF files:
my_files = glob.glob(file_directory + "SWdown_daily_WFDEI_200207*")
if my_files:
    for doc in my_files:
        # Open netCDF file as a Python object:
        f = netcdf.NetCDFFile(doc, "r")
        #
        # Read the year and month from the variable file name:
        f_name = os.path.basename(f.filename)
        try:
            year_month = re.search('_(\d{6})\.', f_name).group(1)
            this_year = int(year_month[0:4])
            this_month = int(year_month[4:6])
        except AttributeError:
            print "Year and month not retrievable from file", f_name
        except ValueError:
            print "Year and/or month not read as numbers in file", f_name
        #
        # Extract data from the variable of interest:
        data = f.variables[var_voi].data
        #
        # Output data to raster format:
        process_raster(file_directory, this_year, this_month, data)
        #
        # Close netCDF file:
        f.close()
else:
    print "Found no files in directory:", file_directory
