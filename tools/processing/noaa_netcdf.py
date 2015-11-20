#!/usr/bin/python
#
# noaa_netcdf.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-22 -- created
# 2014-12-01 -- last updated
#
###############################################################################
#              DO NOT USE THIS DATA IN THE GePiSaT DATABASE!!!!
###############################################################################
#
# ------------
# description:
# ------------
# This script reads the netCDF file of 0.5 degree resolution mean monthly 
# soil moisture data based on van den Dool, Huang & Fan (2003).
#
# Downloaded from NOAA's Physical Sciences Division
# URL: http://www.esrl.noaa.gov/psd/data/gridded/data.cpcsoil.html
#
# NOTE:
#   Data is mean monthly soil moisture from Jan 1948 to Jan 2012
#   * 792 months
#   Ocean's have prescribed soil moisture set equal to 0
#
# NOTE:
#   'lat' goes from 89.75 -- -89.75 (north to south)
#   'lon' goes from 0.25 -- 359.75 (equator -> eastwards?)
#
# ----------
# changelog:
# ----------
# 00. created [14.01.22]
# 01. added raster output [14.01.23]
# 02. included 'soilw' variable attributes [14.01.23]
# 03. general housekeeping [14.12.01]
#
###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import numpy
from scipy.io import netcdf
#
###############################################################################
## FUNCTIONS
###############################################################################
def writeout(f, d):
    """
    Name:     writeout
    Input:    - string, file name with path (t)
              - string, data to be written to file (d)
    Output:   None
    Features: Writes new/overwrites existing file with data string
    """
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

###############################################################################
## CONSTANTS
###############################################################################
file_directory = (
    "/Users/twdavis/Projects/data/noaa/soil_moisture/netcdf/"
    )
my_file = glob.glob(file_directory + "*nc")[0]

###############################################################################
## MAIN
###############################################################################
# Read netcdf file:
f = netcdf.NetCDFFile(my_file, "r")

# Save variables of interest:
f_lat = f.variables['lat'].data
f_lon = f.variables['lon'].data
f_day = f.variables['time'].data
f_msm = f.variables['soilw'].data
sf_msm = f.variables['soilw'].scale_factor
offset = f.variables['soilw'].add_offset
no_value = f.variables['soilw'].missing_value

# Close netcdf file:
f.close()

# Convert lon [0 to 360] to lon [-180 to 180]
new_lon = numpy.array([])
for lon in f_lon:
    if lon > 180:
        new_lon = numpy.append(new_lon, [-1.0*(360.-lon),])
    else:
        new_lon = numpy.append(new_lon, [lon,])

# Save the shape of the data:
# NOTE: should be: sh_time=792, sh_lat=360, sh_lon=720
(sh_time, sh_lat, sh_lon) = f_msm.shape

# Initialize base time stamp (based on time variable units):
# NOTE: 'time' units are days since 1 Jan 1800
base_time = datetime.date(1800, 1, 1)

# Define header line for ASCII raster:
header = (
    "NCOLS 720\n"
    "NROWS 360\n"
    "XLLCORNER -180.0\n"
    "YLLCORNER -90.0\n"
    "CELLSIZE 0.5\n"
    "NODATA_VALUE -9999\n"
    )

# Base values for incrementing over lat & lon:
#lat_base = -89.75
lon_base = -179.75
    
# Iterate through months:
for t in xrange(sh_time):
    # Update time stamp value:
    cur_time = base_time + datetime.timedelta(days = f_day[t])
    #
    rast_out = "%s%s_%s.txt" % (file_directory, 
                                "NOAA_W_0.5-Raster", 
                                cur_time)
    #
    # Save header line to file:
    writeout(rast_out, header)
    #
    for y in xrange(sh_lat):
        #lat = lat_base + 0.5*y
        #lat_idx = numpy.where(f_lat == lat)[0][0]
        #
        # Reset outline for next row:
        outline = ""
        #
        for x in xrange(sh_lon):
            lon = lon_base + 0.5*x
            lon_idx = numpy.where(new_lon == lon)[0][0]
            msm = f_msm[t,y,lon_idx]
            if msm == no_value:
                msm = -9999.0
            else:
                msm = 1.0*msm*sf_msm + offset
            outline = "%s%d " % (outline, msm)
            #
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(rast_out, 'a')
        OUT.write(outline)
        OUT.close()
