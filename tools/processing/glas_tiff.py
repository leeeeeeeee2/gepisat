#!/usr/bin/python
#
# glas_tiff.py
# * based on:
#   georaster.py (TW Davis 2013-11-26),
#   modis_hdf.py (TW Davis, 2013-08-21), 
#   imageproc.py (TW Davis, 2012-03-02),
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-11-25 -- created
# 2014-12-01 -- last updated
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
# 02. general housekeeping [14.12.01]
#
###############################################################################
## IMPORT MODULES
###############################################################################
import glob
import Image
import numpy

###############################################################################
## FUNCTIONS
###############################################################################
def get_1km_grid(hdg_lon, hdg_lat):
    """
    Name:     get_1km_grid
    Input:    - float, half-degree grid longitude, degrees (hdg_lon)
              - float, half-degree grid latitude, degrees (hdg_lat)
    Output:   list of tuples, longitude and latitude pairs, degrees
    Features: Converts 0.5 deg lon-lat pair to 3600 lon-lat pairs at 1 km
    """
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
    """
    Name:     get_1km_rh
    Input:    - float, longitude at 0.5 degree resolution, degrees (lon)
              - float, latitude at 0.5 degree resolution, degrees (lat)
              - numpy.ndarray, canopy height data at 1 km resolution (data)
    Output:   numpy.ndarray, canopy height data at 1 km resolution w.r.t. a 
              given 0.5 pixel (my_grid_rh)
    Features: Returns array of RH values at 1 km for a single 0.5 pixel
    Depends:  - get_1km_grid
              - grid_to_index
    """
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

def get_x_y(lon, lat, r):
    """
    Name:     get_x_y
    Input:    - float, longitude, degrees (lon)
              - float, latitude, degrees (lat)
              - float, resolution (r)
    Output:   tuple, x-y index pair
    Features: Returns x and y indices for lon-lat pair at given resolution 
              based on a numbering scheme starting from the top-left corner
              (i.e., north-west corner)
    """
    x = (lon + 180.0)/r - 0.5
    y = (90.0 - lat)/r - 0.5
    return (int(x), int(y))

def grid_to_index(grid):
    """
    Name:     grid_to_index
    Input:    list of tuples, longitude-latitude pairs, degrees (grid)
    Output:   list of tuples, x-y pairs (okm_indices)
    Features: Returns x-y indices for 1 km grid lon-lat pairs
    Depends:  get_x_y
    """
    okm_indices = []
    my_res = 1.0/120.0
    for grid_pair in grid:
        lon, lat = grid_pair
        x, y = get_x_y(lon, lat, my_res)
        okm_indices.append((x,y))
    #
    return okm_indices

def process_hdg_poly(f, d):
    """
    Name:     process_hdg_poly
    Input:    - str, output file name (f)
              - object, Image data (d)
    Output:   None.
    Features: Writes polygon (shapefile) canopy height data, resampled to 0.5 
              degree resolution, to file
    Depends:  - get_lon_lat
              - get_stationid
              - get_1km_rh
              - writeout
    """
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
                ave_rh = (1e2)*ave_rh
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
    """
    Name:     process_hdg_raster
    Input:    - str, output file name (f)
              - object, Image data (d)
    Output:   None.
    Features: Writes ASCII raster canopy height data, resampled 0.5 degree 
              resolution, to file
    Depends:  - writeout
              - get_lon_lat
              - get_1km_rh
    """
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
                ave_rh = (1e2)*ave_rh
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
#file_dir = "/Users/twdavis/Projects/data/glas/"
file_dir = "/home/user/Projects/gepisat/data/GLAS/"
output_poly = "%s%s.csv" % (file_dir, "RH100_05rs-Poly")
output_raster = "%s%s.txt" % (file_dir, "RH100_05rs-Raster")

###############################################################################
## MAIN PROGRAM
###############################################################################
# Find tif file in the file directory:
my_file = glob.glob(file_dir + "*.tif")[0]

# Open image file:
im = Image.open(my_file).getdata()
sh_lon, sh_lat = im.size

#process_hdg_poly(output_poly, im)
process_hdg_raster(output_raster, im)
