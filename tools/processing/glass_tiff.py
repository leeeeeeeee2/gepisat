#!/usr/bin/python
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-11-25 -- created
# 2014-01-28 -- last updated
#
# glas_tiff.py
# * based on:
#   georaster.py (TW Davis 2013-11-26),
#   imageproc.py (TW Davis, 2012-03-02),
#   and modis_hdf.py (TW Davis, 2013-08-21)
#
# ------------
# description:
# ------------
# This script opens a TIFF raster file for processing.
# NOTE: pixel indexing starts in the top-left (NW) corner (similar to MODIS)
#
# ----------
# changelog:
# ----------
# 00. created from imageproc.py [13.11.25]
# 01. added raster processing [13.11.26]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## LOAD MODULES ###############################################################
# /////////////////////////////////////////////////////////////////////////////
import glob
import Image
import numpy

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## FUNCTIONS ##################################################################
# /////////////////////////////////////////////////////////////////////////////
def get_1km_grid(hdg_lon, hdg_lat):
    """Converts 0.5 deg lon-lat pair to 3600 lon-lat pairs at 1 km"""
    # Initialize five one-hundreths grid:
    okm_grid = []
    #
    # Define half-degree and five one-hundredths resolutions:
    hdg_res = 0.5
    okm_res = 1.0/120.0
    #
    # Calculate the binding box at 0.5 deg:
    westing = hdg_lon - 0.5*hdg_res
    northing = hdg_lat + 0.5*hdg_res
    #
    # Initialize centroid offsetting for 1 km grid:
    okm_lon = westing + 0.5*okm_res
    okm_lat = northing - 0.5*okm_res
    #
    # Iterate over the 60x60 box:
    for y in xrange(60):
        lat = okm_lat - y*okm_res
        for x in xrange(60):
            lon = okm_lon + x*okm_res
            okm_grid.append((lon, lat))
    #
    return okm_grid

def get_1km_rh(lon, lat, data):
    """Returns array of RH values at 1 km for a single 0.5 pixel"""
    # Variable definitions:
    #   lon  :: longitude at regular 0.5 degree resolution
    #   lat  :: latitude at regular 0.5 degree resolution
    #   data :: RH data at 1 km resolution
    #
    my_grid_rh = numpy.array([])
    my_grid_pnts = get_1km_grid(lon, lat)
    my_grid_indx = grid_to_index(my_grid_pnts)
    for indx_pair in my_grid_indx:
        x,y = indx_pair
        zval = data.getpixel((x,y))
        my_grid_rh = numpy.append(my_grid_rh, [zval])
    #
    return my_grid_rh

def get_lon_lat(x,y,r):
    """Returns lat-lon pair for x-y index from top-left corner"""
    # x (longitude): 0...7199 (0...719)
    # y (latitude): 0...3599 (0...359)
    pixel_res = r   
    # Offset lat, lon to pixel centroid
    lon = -180.0 + 0.5*pixel_res
    lat = 90.0 - 0.5*pixel_res
    # Offset lat, lon based on pixel index
    lon = lon + x*pixel_res
    lat = lat - y*pixel_res
    #
    return (lon, lat)

def get_stationid(lon, lat):
    """Returns station ID for 0.5 deg pixel"""
    # Station ID is based on 0 being the bottom (south) left (west) corner
    # and 259199 being the top (north) right (east) corner as is used in 
    # the postgreSQL database naming scheme.
    st_id = (
        720.0 * (359.0 - ((90.0 - lat)/0.5 - 0.5))
        + ((lon + 180.0)/0.5 - 0.5)
        )
    return int(st_id)

def get_x_y(lon, lat, r):
    """Returns x and y indices for lon-lat pair at given resolution"""
    # Pixel resolution:
    pixel_res = r
    #
    # Solve x and y indices:
    x = (lon + 180.0)/pixel_res - 0.5
    y = (90.0 - lat)/pixel_res - 0.5
    #
    return (int(x), int(y))

def grid_to_index(grid):
    """Converts lon-lat pairs to indices for 1 km grid"""
    okm_indices = []
    my_res = 1.0/120.0
    for grid_pair in grid:
        lon, lat = grid_pair
        x, y = get_x_y(lon, lat, my_res)
        okm_indices.append((x,y))
    #
    return okm_indices

def process_hdg_poly(f, d):
    """Resample data to 0.5 deg"""
    # Variable definitions:
    #    f :: output file for saving 0.5 deg EVI
    #    d :: data file at original resolution
    #
    # Open and write header line for CSV file:
    header = "id,lon,lat,rh100_cm\n"
    writeout(f, header)
    #
    # Iterate through 0.5 deg lat-lon pairs:
    # --> y = latitude (0...359)
    # --> x = longitude (0...719)
    for y in xrange(360):
        for x in xrange(720):
            (lon, lat) = get_lon_lat(x, y, 0.5)
            stid = get_stationid(lon, lat)
            #
            # Get 1km RH values within this 0.5 cell:
            my_rh_data = get_1km_rh(lon, lat, d)
            #
            # Check that there's data in the array:
            if len(my_rh_data) > 0:
                # Calculate the average for 0.5 pixel:
                ave_rh = my_rh_data.mean()
                #
                # Convert units (meters to cm):
                ave_rh = ave_rh * 100.0
                #
            else:
                # Use the pre-defined error value
                ave_rh = -9999
                #
            # Save average data to file:
            OUT = open(f, 'a')
            outline = "%d,%0.3f,%0.3f,%d\n" % (
                stid, lon, lat, ave_rh
                )
            OUT.write(outline)
            OUT.close()

def process_hdg_raster(f, d):
    """Resample and process data to 0.5 deg raster"""
    # Variable definitions:
    #    f :: output file for saving 0.5 deg
    #    d :: data file with original data
    #
    # Open and write header line for ASCII raster:
    header = (
        "NCOLS 720\n"
        "NROWS 360\n"
        "XLLCORNER -180.0\n"
        "YLLCORNER -90.0\n"
        "CELLSIZE 0.5\n"
        "NODATA_VALUE 0\n"
        )
    writeout(f, header)
    #
    # Iterate through 0.5 deg lat-lon pairs:
    # --> y = latitude (0...359)
    # --> x = longitude (0...719)
    for y in xrange(360):
        # Reset output string:
        outline = ""
        for x in xrange(720):
            (lon, lat) = get_lon_lat(x, y, 0.5)
            #stid = get_stationid(lon, lat)
            #
            # Get RH data within this 0.5 cell:
            my_rh_data = get_1km_rh(lon, lat, d)
            #
            # Check that there's data in the array:
            if len(my_rh_data) > 0:
                # Calculate the average value for 0.5 pixel:
                ave_rh = my_rh_data.mean()
                #
                # Convert units (meters to cm):
                ave_rh = ave_rh * 100.0
                #
            else:
                # Use error value:
                ave_rh = -9999
                #
            outline = "%s%d " % (outline, ave_rh)
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(f, 'a')
        OUT.write(outline)
        OUT.close()

def writeout(f, d):
    """Writes new/overwrites existing file"""
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## CONSTANTS ##################################################################
# /////////////////////////////////////////////////////////////////////////////
#file_dir = "/Users/twdavis/Projects/data/glas/"
file_dir = "/home/user/Projects/gepisat/data/GLAS/"
output_poly = "%s%s.csv" % (file_dir, "RH100_05rs-Poly")
output_raster = "%s%s.txt" % (file_dir, "RH100_05rs-Raster")

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## MAIN #######################################################################
# /////////////////////////////////////////////////////////////////////////////
# Find tif file in the file directory:
my_file = glob.glob(file_dir + "*.tif")[0]

# Open image file:
im = Image.open(my_file).getdata()
sh_lon, sh_lat = im.size


#process_hdg_poly(output_poly, im)
process_hdg_raster(output_raster, im)

