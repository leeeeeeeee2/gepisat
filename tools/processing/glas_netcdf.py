#!/usr/bin/python
#
# glas_netcdf.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-11-20 -- created
# 2014-12-01 -- last updated
#
###############################################################################
#              DO NOT USE THIS DATA IN THE GePiSaT DATABASE!!!!
###############################################################################
#
# ------------
# description:
# ------------
# This script reads the NetCDF file for the 2005 GLAS-derived canopy height 
# based on the work of Simard et al., 2011.
#
# The dataset was downloaded from ORNL DAAC wesite:
# http://webmap.ornl.gov/wcsdown/dataset.jsp?ds_id=10023
#
# The data resolution was set using the web tool to 0.5 degree (WGS84).
# The indexing of pixels starts in the top-left corner, similar to the 
# MODIS HDF files.
#
# ----------
# changelog:
# ----------
# 00. created based on cru_netcdf.py [13.11.20]
# 01. created process_poly() function [13.11.21]
# 02. general housekeeping [14.12.01]
#
# -----
# todo:
# -----
# 1. Make distinction between get_lon_lat function for top-left and bottom-left 
#    numbering (e.g., MODIS HDF versus WATCH netCDF)
#
###############################################################################
## IMPORT MODULES
###############################################################################
import glob
from scipy.io import netcdf

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

def process_raster(data, var_name, out_dir):
    """
    Name:     process_raster
    Input:    - numpy nd.array (data)
              - string, variable name for output file (var_name)
              - string, output directory w/ path (out_dir)
    Output:   None.
    Features: Saves data in ASCII raster format
    Depends:  writeout
    """
    # Define header line for ASCII raster:
    header = (
        "NCOLS 720\n"
        "NROWS 360\n"
        "XLLCORNER -180.0\n"
        "YLLCORNER -90.0\n"
        "CELLSIZE 0.5\n"
        "NODATA_VALUE 0\n"
        )
    #
    # Save the shape values of each for iteration purposes
    lat, lon = data.shape
    #
    # Create output file name:
    rast_out = "%s%s_%s_%s_%s.txt" % (
        out_dir, 
        "GLAS",
        var_name,
        "0.5-Raster", 
        "2005"
        )
    #
    # Save header line to file:
    writeout(rast_out, header)
    #
    # Iterate through current month's data:
    for y in xrange(lat):
        # Reset outline for next row:
        outline = ""
        #
        for x in xrange(lon):
            val = data[y,x]
            #
            # Convert to integer:
            # NOTE: data source is already in integer form
            val = int(val)
            #
            # Add value to outline:
            outline = "%s%d " % (outline, val)
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(rast_out, 'a')
        OUT.write(outline)
        OUT.close()

def process_poly(d, data):
    """
    Name:     process_poly
    Input:    - string, output file w/ path (d)
              - numpy nd.array (data)
    Output:   None.
    Features: Saves 0.5 deg canopy height data in polygon (shapefile) CSV format
    """
    # Open and write header line for CSV file:
    f = d + "GLAS_RH100_Poly_2005.csv"
    header = "id,lon,lat,rh100\n"
    writeout(f, header)
    #
    # Iterate through 0.5 deg lat-lon pairs:
    # --> y = latitude (0...359)
    # --> x = longitude (0...719)
    for y in xrange(360):
        for x in xrange(720):
            (lon, lat) = get_lon_lat(x, y, 0.5)
            stid = get_stationid(lon, lat)
            val = data[y,x]
            #
            # Save RH100 to file:
            OUT = open(f, 'a')
            outline = "%d,%0.3f,%0.3f,%d\n" % (stid, lon, lat, val)
            OUT.write(outline)
            OUT.close()

def get_lon_lat(x,y,r):
    """
    Name:     get_lon_lat
    Input:    - int/nd.array, longitude index (x)
              - int/nd.array, latitude index (y)
              - float, pixel resolution (r)
    Output:   float/nd.array tuple, longitude(s) and latitude(s), degrees
    Features: Returns lat-lon pair for x-y index pair (numbered from top-left 
              corner) and pixel resolution
    """
    # Offset lat, lon to pixel centroid
    lon = -180.0 + (0.5*r)
    lat = 90.0 - (0.5*r)
    #
    # Offset lat, lon based on pixel index
    lon = lon + (x*r)
    lat = lat - (y*r)
    #
    return (lon, lat)

def get_stationid(lon, lat):
    """
    Name:     get_stationid
    Input:    - float, longitude, degrees (lon)
              - float, latitude, degrees (lat)
    Output:   int, station id (st_id)
    Features: Returns the half-degree (HDG) station ID for a pixel
              numbered from 0 (bottom-left / south-west corner) to 259199 
              (top-right / north-east corner) as defined in the GePiSaT 
              database numbering scheme
    """
    st_id = 720.0*(359.0 - ((90.0 - lat)/0.5 - 0.5)) + ((lon + 180.0)/0.5 - 0.5)
    return int(st_id)

###############################################################################
## CONSTANTS
###############################################################################
file_directory = (
    "/home/user/Projects/gepisat/data/GLAS/"
    )

###############################################################################
## MAIN
###############################################################################
# Read netcdf file in the file directory:
my_file = glob.glob(file_directory + "*.nc")[0]

# Open netcdf file for reading:
f = netcdf.NetCDFFile(my_file, "r")

# Save variable of interest:
rh100 = f.variables['Band1'].data

# Close netcdf file:
f.close()

process_raster(rh100, "RH100", file_directory)
process_poly(file_directory, rh100)
