#!/usr/bin/python
#
# modis_hdf.py
#
# written by Tyler W. Davis
# Imperial College London
# (Enthought Canopy Python Environment)
#
# 2013-08-21 -- created
# 2013-10-23 -- last updated
#
# ------------
# description:
# ------------
# This script reads MODIS land-surface products (in HDF file format) and 
# produces a CSV file for importing geo-referenced data into GIS.
# --> includes assigning longitude and latitude to measurements based on 
#     the gridded data's resolution
# This script can also resample the MODIS CGM (0.05 deg) grid to a lower
# resolution (e.g., 0.5 deg) and process EVI to file.
#
# ----------
# changelog:
# ----------
# 00. created in Canopy editor [13.08.21]
# 01. Updated get_lon_lat function & fixed function call in main [13.09.23]
# --> added resolution to function parameters [13.09.24]
# 02. Created get_x_y() function to convert lon-lat to indices [13.09.24]
# 03. Abstracted get_evi() function [13.09.24]
# 04. Created get_foh_grid() function [13.09.24]
# 05. Created grid_to_index() function [13.09.24]
# 06. Created get_foh_evi() function [13.09.24]
# 07. Abstracted 0.05 deg processing to process_foh() function [13.09.24]
# 08. Create process_hdg() function [13.09.24]
# --> performs averaging for EVI resampling from 0.05 to 0.5 deg resolution
# 09. Leave EVI in units of zval (i.e., 10000 * evi) [13.09.24]
# 10. Created get_stationid() function [13.09.25]
# --> added station id to process_hdg output [13.09.25]
# 11. Remove pixel filter [13.09.25]
# --> all data is needed for rasters
# --> process_hdg uses numpy.nan values for filtering averages
# 12. Separated raster and poly output formats [13.09.26]
# 13. Changed averaging method for 0.5 resampling [13.09.27]
# --> divide by 100; should scale back pixels with few points
# --> changed back to mean (performs poorly against original 0.05 data)
# 14. added file looping [13.10.23]
#
#
#  Table 1. Standard deviations of the differences between raster images of 
#           an original 0.05 raster image and 0.5 deg raster image resolutions 
#           using three different resampling techniques: Python mean (old), 
#           per 100 (old), Group stats ave. (Qgis).
#
#           Raster calc (Ra - Rb)    +/- 2 & 4 st. dev.
#           ----------------------  -------------------------
#           0.05 rast - 0.5 (old):  [-1289, -427, 433, 1295] *
#           0.05 rast - 0.5 (Qgis): [-1308, -433, 441, 1316]
#           0.05 rast - 0.5 (new):  [-1407, -423, 560, 1543]
#           -
#           0.5 (Qgis) - 0.5 (old): [-686, -244, 197, 639]   *
#           0.5 (Qgis) - 0.5 (new): [-1078, -251, 575, 1403]
#
#
# -----------
# references:
# -----------
# Pyhdf SD class: http://pysclint.sourceforge.net/pyhdf/pyhdf.SD.html
#
# -----
# todo:
# -----
# x Implement resolution resampling (0.05 to 0.5 deg) [13.09.24]
# x Reprocess FOH and HDG files with EVI in integer (100000x's) units
# x Output to ASCII Raster format [13.09.26]
#
##############################################################################
## IMPORT MODULES
##############################################################################
import os.path
import datetime
import glob
import re
import numpy
from pyhdf import SD   # installed 0.8.3 in canopy package manager

##############################################################################
## FUNCTION DEFINITIONS
##############################################################################
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

def get_x_y(lon, lat, r):
    """Returns x and y indices for lon-lat pair at 0.05 deg resolution"""
    # Pixel resolution:
    pixel_res = r
    #
    # Solve x and y indices:
    x = (lon + 180.0)/pixel_res - 0.5
    y = (90.0 - lat)/pixel_res - 0.5
    #
    return (int(x), int(y))

def writeout(f, d):
    """Writes new/overwrites existing file"""
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

def get_ts(f):
    """Reads filename to retrieve timestamp"""
    f_name = os.path.basename(f)
    try:
        f_date = re.search('A(\d{7})\.', f_name).group(1)
    except AttributeError:
        print "Search failed for date in file", f_name
    else:
        f_year = int(f_date[0:4])
        f_doy = int(f_date[4:7]) - 1
        f_timestamp = (
            datetime.date(f_year, 1, 1) + 
            datetime.timedelta(days=f_doy)
            )
        return f_timestamp

def get_evi(f):
    """Opens an HDF file and returns an SD object"""
    if os.path.isfile(f):
        # Open SD file:
        try:
            my_hdf = SD.SD(f)
        except:
            print "Unexpected error opening file", f
        else:
            # Read dataset names and select EVI:
            # * name should be "CMG 0.05 Deg Monthly EVI"
            # * also d_key[8]
            #
            #d_key = my_hdf.datasets().keys()
            #d_name = d_key[8]
            key_name = "CMG 0.05 Deg Monthly EVI"
            #d_set = f.datasets()[d_name]
            if key_name in my_hdf.datasets().keys():
                # Pull data from dataset into array:
                d_select = my_hdf.select(key_name)
                data = d_select.get()
            else:
                print "Could not open dataset", key_name
                data = None
            #
            # Close HDF file and return data:
            my_hdf.end()
            return data
    else:
        print "File does not exist!", f

def get_foh_grid(hdg_lon, hdg_lat):
    """Converts 0.5 deg lon-lat pair to 100 lon-lat pairs at 0.05 deg"""
    # Initialize five one-hundreths grid:
    foh_grid = []
    #
    # Define half-degree and five one-hundredths resolutions:
    hdg_res = 0.5
    foh_res = 0.05
    #
    # Calculate the binding box at 0.5 deg:
    westing = hdg_lon - 0.5*hdg_res
    northing = hdg_lat + 0.5*hdg_res
    #easting = hdg_lon + 0.5*hdg_res
    #southing = hdg_lat - 0.5*hdg_res
    #
    # Initialize centroid offsetting for foh_grid:
    foh_lon = westing + 0.5*foh_res
    foh_lat = northing - 0.5*foh_res
    #
    # Iterate over the 10x10 box:
    for y in xrange(10):
        lat = foh_lat - y*foh_res
        for x in xrange(10):
            lon = foh_lon + x*foh_res
            foh_grid.append((lon, lat))
    #
    return foh_grid

def grid_to_index(grid):
    """Converts lon-lat pairs to indices for 0.05 deg grid"""
    foh_indices = []
    for grid_pair in grid:
        lon, lat = grid_pair
        x, y = get_x_y(lon, lat, 0.05)
        foh_indices.append((x,y))
    #
    return foh_indices

def get_foh_evi(lon, lat, data):
    """Returns array of EVI values at 0.05 for a single 0.5 pixel"""
    # Variable definitions:
    #   lon  :: longitude at regular 0.5 degree resolution
    #   lat  :: latitude at regular 0.5 degree resolution
    #   data :: EVI data array at 0.05 degree resolution
    #
    my_grid_evi = numpy.array([])
    my_grid_pnts = get_foh_grid(lon, lat)
    my_grid_indx = grid_to_index(my_grid_pnts)
    for indx_pair in my_grid_indx:
        x,y = indx_pair
        zval = data[y][x]
        # Use NaN for invalid data (easy to remove for averaging):
        if zval < -2000:
            zval = numpy.NaN
        my_grid_evi = numpy.append(my_grid_evi, [zval])
    #
    return my_grid_evi

def process_foh_raster(f, d):
    """Processes 0.05 degree EVI data to file"""
    # Variable definitions:
    #    f :: output file to save EVI
    #    d :: data holding 0.05 EVI
    #
    # Open and write header line for ASCII raster:
    header = (
        "NCOLS 7200\n"
        "NROWS 3600\n"
        "XLLCORNER -180.00\n"
        "YLLCORNER -90.00\n"
        "CELLSIZE 0.05\n"
        "NODATA_VALUE -3000\n"
        )
    writeout(f, header)
    #
    #  Save the dimensions of the array:
    (ydim, xdim) =  d.shape
    #
    ## Extract EVI from each pixel:
    for yval in xrange(ydim):
        # Reset outline string:
        outline = ""
        for xval in xrange(xdim):
            # Row-major ordering:
            zval = d[yval][xval]
            #evi = zval / 10000.0
            #
            # Add new value to string:
            outline = "%s%d " % (outline, zval)
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(f, 'a')
        OUT.write(outline)
        OUT.close()

def process_foh_poly(f, d):
    """Processes 0.05 degree EVI data to file"""
    # Variable definitions:
    #    f :: output file to save EVI
    #    d :: data holding 0.05 EVI
    #
    # Open and write header line for CSV file:
    header = "id,lon,lat,evi\n"
    writeout(f, header)
    #
    #  Save the dimensions of the array:
    (ydim, xdim) =  d.shape
    #
    ## Extract EVI from each pixel:
    i = 0  # index (ranges from 0...25919999)
    for yval in xrange(ydim):
        for xval in xrange(xdim):
            # Row-major ordering:
            zval = d[yval][xval]
            #evi = zval / 10000.0
            #
            # Retrieve lon and lat values for this pixel:
            (lon, lat) = get_lon_lat(xval, yval, 0.05)
            #
            # Valid range for values -2000 to 10000:
            #if zval >= -2000 and zval <= 10000:
            #    # Save data to file (~170 MB/file):
            OUT = open(f, 'a')
            outline = "%d,%0.3f,%0.3f,%d\n" % (i, lon, lat, zval)
            OUT.write(outline)
            OUT.close()
            # Increment counter:
            i += 1

def process_hdg_raster(f, d):
    """Resample and process 0.05 deg EVI to 0.5 deg"""
    # Variable definitions:
    #    f :: output file for saving 0.5 deg EVI
    #    d :: data file with 0.05 deg EVI
    #
    # Open and write header line for ASCII raster:
    header = (
        "NCOLS 720\n"
        "NROWS 360\n"
        "XLLCORNER -180.0\n"
        "YLLCORNER -90.0\n"
        "CELLSIZE 0.5\n"
        "NODATA_VALUE -3000\n"
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
            # Get 0.05 EVI values within this 0.5 cell:
            my_evi_data = get_foh_evi(lon, lat, d)
            #
            # Check that there's data in the array:
            if len(my_evi_data[~numpy.isnan(my_evi_data)]) > 0:
                # Calculate the average EVI for 0.5 pixel:
                ave_evi = my_evi_data[~numpy.isnan(my_evi_data)].mean()
                #
            else:
                # Use the pre-defined error value (i.e., -3000)
                ave_evi = -3000
                #
            outline = "%s%d " % (outline, ave_evi)
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(f, 'a')
        OUT.write(outline)
        OUT.close()

def process_hdg_poly(f, d):
    """Resample and process 0.05 deg EVI to 0.5 deg"""
    # Variable definitions:
    #    f :: output file for saving 0.5 deg EVI
    #    d :: data file with 0.05 deg EVI
    #
    # Open and write header line for CSV file:
    header = "id,lon,lat,evi\n"
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
            # Get 0.05 EVI values within this 0.5 cell:
            my_evi_data = get_foh_evi(lon, lat, d)
            #
            # Check that there's data in the array:
            if len(my_evi_data[~numpy.isnan(my_evi_data)]) > 0:
                # Calculate the average EVI for 0.5 pixel:
                ave_evi = my_evi_data[~numpy.isnan(my_evi_data)].mean()
                #
            else:
                # Use the pre-defined error value (i.e., -3000)
                ave_evi = -3000
                #
            # Save averaged EVI to file:
            OUT = open(f, 'a')
            outline = "%d,%0.3f,%0.3f,%d\n" % (
                stid, lon, lat, ave_evi
                )
            OUT.write(outline)
            OUT.close()

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

##############################################################################
## MAIN
##############################################################################
# Open directory with HDF files and read file names:
#mydir = (
#    "/Users/twdavis/Qgis/resampling_modis/2002_07/"
#    )
mydir = (
    "/Users/twdavis/Projects/data/modis/vi_cgm_monthly/aqua/"
    )
myfiles = glob.glob(mydir+'MYD13C2.A2006*.hdf')
#myfile = myfiles[0]

for myfile in myfiles:
    # Get time value for this file:
    myts = get_ts(myfile)
    #
    # Get EVI data:
    data = get_evi(myfile)
    #
    if data.any():
        print "Processing month: %s" % myts
        #
        # Process original 0.05 data:
        #outfile_1 = "%s%s_%s.txt" % (mydir, "MODIS_0.05-Raster", myts)
        #process_foh_raster(outfile_1, data)
        #
        # Resample to 0.5 degree resolution:
        outfile_2 = "%s%s_%s.txt" % (mydir, "MODIS_0.5rs-Raster", myts)
        process_hdg_raster(outfile_2, data)
        #
        # Output resampled 0.5 data as poly:
        #outfile_3 = "%s%s_%s.csv" % (mydir, "MODIS_05rs-Poly", myts)
        #process_hdg_poly(outfile_3, data)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# Save data to array:
#d_array = numpy.array(d_list)
#numpy.reshape(d_array, (3600,7200), order='C')

# Plot 2D graphic (grayscale):
#plt.imshow(numpy.reshape(d_array, (3600,7200), order='C'))
#plt.gray()
#plt.show()