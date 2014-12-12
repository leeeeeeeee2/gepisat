#!/usr/bin/python
#
# modis_hdf.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-08-21 -- created
# 2014-12-11 -- last updated
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
# 15. general housekeeping [14.12.01]
# 16. started flux tower time series processing [14.12.11]
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
# 0. Finish flux station EVI time-series processing
# 1. Distinguish get_lon_lat function for top-left and bottom-left numbering
#    schemes (e.g., MODIS HDF versus WATCH netCDF).
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
## FUNCTIONS
##############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32) 
    dt3 = dt2.replace(day=1)
    return dt3

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

def get_ts(f):
    """
    Name:     get_ts
    Input:    str, filename with path (f)
    Output:   datetime.datetime
    Features: Returns the timestamp (datetime object) associated with a MODIS
              file name (e.g., MYD13C2.A2002182.005.2007199145757.hdf)
    """
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
    """
    Name:     get_evi
    Input:    str, input file name with path (f)
    Output:   numpy.ndarray, MODIS CMG 0.05 Deg Monthly EVI
    Features: Returns the numpy array of MODIS monthly EVI from a given file
    """
    if os.path.isfile(f):
        try:
            # Opens HDF4 file type:
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
            #d_set = f.datasets()[d_name]
            #
            key_name = "CMG 0.05 Deg Monthly EVI"
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
    """
    Name:     get_foh_grid
    Input:    - float, half-degree grid longitude, degrees (hdg_lon)
              - float, half-degree grid latitude, degrees (hdg_lat)
    Output:   list of tuples, lon-lat pairs at 0.05 degree
    Features: Returns the 100 lon-lat pairs at 0.05 deg corresponding to a 
              single 0.5 pixel
    """
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
    """
    Name:     grid_to_index
    Input:    list of tuples, lon-lat pairs (grid)
    Output:   list of tuples, x-y pairs
    Features: Returns the x-y indices for lon-lat pairs at 0.05 degree 
              resolution
    Depends:  get_x_y
    """
    foh_indices = []
    for grid_pair in grid:
        lon, lat = grid_pair
        x, y = get_x_y(lon, lat, 0.05)
        foh_indices.append((x,y))
    #
    return foh_indices

def get_foh_evi(lon, lat, data):
    """
    Name:     get_foh_evi
    Input:    - float, longitude at 0.5 deg resolution (lon)
              - float, latitude at 0.5 deg resolution (lat)
              - numpy.ndarray, MODIS EVI data at 0.05 deg resolution (data)
    Output:   numpy.ndarray, 0.05 deg MODIS EVI data for given 0.5 deg pixel
    Features: Returns array of EVI values at 0.05 for a single 0.5 pixel
    Depends:  - get_foh_grid
              - grid_to_index
    """
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
    """
    Name:     process_foh_raster
    Input:    - str, output file (f)
              - numpy.ndarray, 0.05 deg MODIS EVI data (d)
    Output:   None.
    Features: Writes 0.05 degree EVI data to ASCII raster file
    Depends:  writeout
    """
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
    """
    Name:     process_foh_poly
    Input:    - str, output file (f)
              - numpy.ndarray, 0.05 deg MODIS EVI data (d)
    Features: Writes 0.05 degree EVI data to polygon (shapefile CSV) file
    Depends:  - writeout
              - get_lon_lat
    """
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
    """
    Name:     process_hdg_raster
    Input:    - str, output files (f)
              - numpy.ndarray, 0.05 deg MODIS EVI data (d)
    Output:   None.
    Features: Writes 0.5 degree EVI to ASCII raster file, resampled based on
              0.05 degree MODIS EVI
    Depends:  - writeout
              - get_lon_lat
              - get_foh_evi
    """
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
    """
    Name:     process_hdg_poly
    Input:    - str, output file (f)
              - numpy.ndarray, 0.05 deg MODIS EVI data (d)
    Output:   None.
    Features: Writes 0.5 degree EVI data to vector (shapefile CSV) file based 
              on 0.05 degree MODIS EVI
    Depends:  - writeout
              - get_stationid
              - get_foh_evi
              - get_lon_lat
    """
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

def grid_centroid(my_lon, my_lat, grid_res):
    """
    Name:     grid_centroid
    Input:    - float, longitude, degrees (my_lon)
              - float, latitude, degrees (my_lat)
              - float, grid resolution, degrees (grid_res)
    Output:   tuple, longitude-latitude pair (my_centroid)
    Features: Returns the nearest grid centroid per given coordinates and 
              resolution based on the Euclidean distance to each of the four 
              surrounding grids; if any distances are equivalent, the pixel 
              north and east is selected by default
    """
    # Create lists of regular latitude and longitude:
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = int(180.0/grid_res)
    lon_dim = int(360.0/grid_res)
    lats = [lat_min + y*grid_res for y in xrange(lat_dim)]
    lons = [lon_min + x*grid_res for x in xrange(lon_dim)]
    #
    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index-1]
        try:
            bb_lon_max = lons[lon_index+1]
        except IndexError:
            bb_lon_max = lons[-1] + grid_res
        #
    # Find bounding latitude:
    centroid_lat = None
    if my_lat in lats:
        centroid_lat = my_lat
    else:
        lats.append(my_lat)
        lats.sort()
        lat_index = lats.index(my_lat)
        bb_lat_min = lats[lat_index-1]
        try:
            bb_lat_max = lats[lat_index+1]
        except IndexError:
            bb_lat_max = lats[-1] + grid_res
        #
    # Determine nearest centroid:
    # NOTE: if dist_A equals dist_B, then centroid defaults positively (NE)
    if centroid_lon and centroid_lat:
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lon and not centroid_lat:
        # Calculate the distances between lat and bounding box:
        dist_A = bb_lat_max - my_lat
        dist_B = my_lat - bb_lat_min
        if dist_A > dist_B:
            centroid_lat = bb_lat_min
        else:
            centroid_lat = bb_lat_max
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lat and not centroid_lon:
        # Calculate the distances between lon and bounding box:
        dist_A = bb_lon_max - my_lon
        dist_B = my_lon - bb_lon_min
        if dist_A > dist_B:
            centroid_lon = bb_lon_min
        else:
            centroid_lon = bb_lon_max
        my_centroid = (centroid_lon, centroid_lat)
    else:
        # Calculate distances between lat:lon and bounding box:
        # NOTE: if all distances are equal, defaults to NE grid
        dist_A = numpy.sqrt((bb_lon_max - my_lon)**2 + (bb_lat_max - my_lat)**2)
        dist_B = numpy.sqrt((bb_lon_max - my_lon)**2 + (my_lat - bb_lat_min)**2)
        dist_C = numpy.sqrt((my_lon - bb_lon_min)**2 + (bb_lat_max - my_lat)**2)
        dist_D = numpy.sqrt((my_lon - bb_lon_min)**2 + (my_lat - bb_lat_min)**2)
        min_dist = min([dist_A, dist_B, dist_C, dist_D])
        #
        # Determine centroid based on min distance:
        if dist_A == min_dist:
            my_centroid = (bb_lon_max, bb_lat_max)
        elif dist_B == min_dist:
            my_centroid = (bb_lon_max, bb_lat_min)
        elif dist_C == min_dist:
            my_centroid = (bb_lon_min, bb_lat_max)
        elif dist_D == min_dist:
            my_centroid = (bb_lon_min, bb_lat_min)
            #
    # Return nearest centroid:
    return my_centroid

##############################################################################
## MAIN PROGRAM:
##############################################################################
# Directory naming:
mac = 0
if mac:
    terra_dir = "/Users/twdavis/Projects/data/modis/vi_cgm_monthly/terra/"
    aqua_dir = "/Users/twdavis/Projects/data/modis/vi_cgm_monthly/aqua/"
    met_dir = '/Users/twdavis/Dropbox/Work/Imperial/flux/data/psql-data/flux/'
    out_dir = "/Users/twdavis/Desktop/"
else:
    terra_dir = '/usr/local/share/database/modis/evi/terra/'
    aqua_dir = '/usr/local/share/database/modis/evi/aqua/'
    met_dir = '/home/user/Dropbox/Work/Imperial/flux/data/psql-data/flux/'
    out_dir = '/home/user/Desktop/'
#
# Read Aqua and Terra HDF4 files; note multiple years of data in Aqua dir
my_files = []
for y in range(2002, 2007):
    if y == 2002:
        # Terra
        temp_files = glob.glob(terra_dir + 'MOD13C2.A2002*.hdf')
    else:
        # Aqua
        temp_files = glob.glob(aqua_dir + 'MYD13C2.A' + str(y) + '*.hdf')
    if temp_files:
            for temp_f in numpy.sort(temp_files):
                my_files.append(temp_f)

# Read flux station meta data (if processing EVI for each site)
met_files = glob.glob(met_dir + '*2002-06.csv')
met_file = met_files[0]
met_data = numpy.loadtxt(
    met_file, 
    delimiter=",", 
    skiprows=1,
    usecols=(4, 6, 7),
    dtype={'names' : ('station', 'lat', 'lon'),
           'formats' : ('S6', 'f4', 'f4')},
)

# Prepare array for EVI time series:
monthly_evi_temp = numpy.repeat(numpy.nan, 60)
i = 0
for d in met_data:
    my_tuple = ()
    my_formats = ()
    my_names = ()
    #
    my_station, my_lat, my_lon = d
    my_tuple += (my_station, my_lat, my_lon)
    my_formats += ('S6', 'f4', 'f4')
    my_names += ('station', 'lat', 'lon')
    #
    grid_lon, grid_lat = grid_centroid(my_lon, my_lat, 0.05)
    my_tuple += (grid_lat, grid_lon)
    my_formats += ('f4', 'f4')
    my_names += ('grid_lat', 'grid_lon')
    #
    grid_x, grid_y = get_x_y(grid_lon, grid_lat, 0.05)
    my_tuple += (grid_y, grid_x)
    my_formats += ('i4', 'i4')
    my_names += ('grid_y', 'grid_x')
    #
    base_date = datetime.date(2002, 1, 1)
    for evi in monthly_evi_temp:
        my_tuple += (evi,)
        my_formats += ('f4',)
        date_name = 'EVI_%d-%02d' % (base_date.year, base_date.month)
        my_names += (date_name,)
        base_date = add_one_month(base_date)
    #
    if i == 0:
        all_station_data = numpy.array(my_tuple, 
                                       dtype={'names' : my_names, 
                                              'formats' : my_formats}, 
                                       ndmin=1)
    else:
        tmp_station_data = numpy.array(my_tuple, 
                                       dtype={'names' : my_names, 
                                              'formats' : my_formats}, 
                                       ndmin=1)
        all_station_data = numpy.append(all_station_data, tmp_station_data, 
                                        axis=0)
    #
    i += 1


for my_file in my_files:
    # Get time value for this file:
    myts = get_ts(my_file)
    myts_name = 'EVI_%d-%02d' % (myts.year, myts.month)
    #
    # Get EVI data:
    data = get_evi(my_file)
    #
    if data.any():
        print "Processing month: %s" % myts
        #
        # Process original 0.05 data:
        #outfile_1 = "%s%s_%s.txt" % (mydir, "MODIS_0.05-Raster", myts)
        #process_foh_raster(outfile_1, data)
        #
        # Resample to 0.5 degree resolution:
        #outfile_2 = "%s%s_%s.txt" % (out_dir, "MODIS_0.5rs-Raster", myts)
        #process_hdg_raster(outfile_2, data)
        #
        # Output resampled 0.5 data as poly:
        #outfile_3 = "%s%s_%s.csv" % (mydir, "MODIS_05rs-Poly", myts)
        #process_hdg_poly(outfile_3, data)
        #
        # Fill all_station_data array:
        num_stations = len(all_station_data['station'])
        for j in xrange(num_stations):
            my_x = all_station_data[j]['grid_x']
            my_y = all_station_data[j]['grid_y']
            my_evi = data[my_y][my_x]
            # Valid range for values -2000 to 10000:
            if my_evi < -2000:
                my_evi = numpy.nan
            elif my_evi > 10000:
                my_evi = 10000
            all_station_data[j][myts_name] = (1e-4)*my_evi

# Print out all_station_data
out_file = 'FLUXNET_MODIS-0.05-EVI_2002-06.txt'
my_header = ''
for my_name in my_names:
    my_header += my_name
    my_header += ','
my_header = my_header.rstrip(',')
my_header += '\n'
writeout(out_dir + out_file, my_header)
out_formats = '%s,%f,%f,%0.3f,%0.3f,%d,%d,' + '%0.4f,'*59 + '%0.4f\n' 
for station_data in all_station_data:
    try:
        f = open(out_dir + out_file, 'a')
    except IOError:
        print "Cannot write to file:", out_file
    else:
        f.write(out_formats % tuple(station_data))
        f.close()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# Save data to array:
#d_array = numpy.array(d_list)
#numpy.reshape(d_array, (3600,7200), order='C')

# Plot 2D graphic (grayscale):
#plt.imshow(numpy.reshape(d_array, (3600,7200), order='C'))
#plt.gray()
#plt.show()