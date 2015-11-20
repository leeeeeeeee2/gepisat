#!/usr/bin/python
#
# cru_dat.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-23 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script reads CRU TS 3.00 data file (i.e., .dat) containing the half-
# degree land-surface elevation data.
#
# Retrieved from BADC web server (username and password required).
# filename: halfdeg.elv.grid.dat
#
# NOTE:
#   'lat' goes from -89.75 -- 89.75 (south to north)
#   'lon' goes from -179.75 -- 179.75 (east to west)
#
# ----------
# changelog:
# ----------
# 00. created [14.01.23]
# 01. general housekeeping [14.12.01]
#
###############################################################################
## IMPORT MODULES
###############################################################################
import glob
import numpy

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
    "/Users/twdavis/Projects/data/cru/cru_ts_3_00/"
    )
my_file = glob.glob(file_directory + "*dat")[0]

###############################################################################
## MAIN
###############################################################################
# Open and read dat file:
# NOTE: data is read into an array with shape (360, 720)
f = numpy.loadtxt(my_file)
(sh_lat, sh_lon) = f.shape

# Define header line for ASCII raster:
header = (
    "NCOLS 720\n"
    "NROWS 360\n"
    "XLLCORNER -180.0\n"
    "YLLCORNER -90.0\n"
    "CELLSIZE 0.5\n"
    "NODATA_VALUE -9999\n"
    )

# Define/write-out raster output filename:
rast_out = "%s%s.txt" % (file_directory, "CRU-TS3.00_elv_0.5-Raster")
writeout(rast_out, header)

# Iterate through data:
for i in xrange(sh_lat):
    # Reverse latitude:
    y = sh_lat - 1 - i
    #
    # Reset outline for next row:
    outline = ""
    #
    for x in xrange(sh_lon):
        elv = f[y,x]
        # Set no-data value:
        if elv < -500:
            elv = -9999
        outline = "%s%d " % (outline, elv)
        #
    # End line and print to file:
    outline = "%s\n" % outline.rstrip(' ')
    OUT = open(rast_out, 'a')
    OUT.write(outline)
    OUT.close()