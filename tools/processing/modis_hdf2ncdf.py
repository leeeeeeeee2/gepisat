#!/usr/bin/python
#
# modis_hdf2ncdf.py
# * based on modis_hdf.py (2013-10)
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-08-21 -- created
# 2014-11-24 -- last updated
#
#
# ------------
# description:
# ------------
# This script reads MODIS land-surface products (in HDF file format),
# upscales the data to 0.5 degree resolution, and saves the data to 
# netCDF file format.
#
# ----------
# changelog:
# ----------
# ---- modis_hdf.py ----
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
# ---- modis_hdf2ncdf.py ----
# 15. added scipy.io.netcdf to preamble [14.05.06]
# 16. removed unnecessary functions [14.05.06]
# --> e.g. process_raster & process_poly
# 17. added new functions for netcdf meta data & lon, lat variables [14.05.06]
# 18. created upscale_evi function [14.05.07] 
# 19. added time variable & update processing to one file output [14.05.08]
# 20. updated function doc [14.09.17]
# 21. started MODIS atmospheric data processing [14.11.24]
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
#
##############################################################################
## IMPORT MODULES
##############################################################################
import os.path
import datetime
import glob
import re
import numpy
from scipy.io import netcdf
from pyhdf import SD   # installed 0.8.3 in canopy package manager

##############################################################################
## GLOBAL VARIABLES:
##############################################################################
ERROR_VAL = 1.0e+20

##############################################################################
## FUNCTION DEFINITIONS
##############################################################################
def get_lon_lat(x,y,r):
    """
    Name:     get_lon_lat
    Input:    - int, latitude index, e.g., 0--319 (x)
              - int, longitude index, e.g., 0--719 (y)
              - float, resolution, e.g., 0.5 (r)
    Output:   tuple, longitude-latitude pair (lon, lat)
    Features: Returns lon-lat pair for an x-y index pair (numbered from the 
              bottom-left corner) and pixel resolution
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
    Input:    - float, longitude (lon)
              - float, latitude (lat)
              - float, resolution (r)
    Output:   tuple, x-y indices
    Features: Returns x and y indices for a given lon-lat pair and pixel 
              resolution
    """
    x = (lon + 180.0)/r - 0.5
    y = (90.0 - lat)/r - 0.5
    #
    return (int(x), int(y))

def get_ts(f):
    """
    Name:     get_ts
    Input:    string, file name w/ path (f)
    Output:   datetime date (f_timestamp)
    Features: Returns the timestamp of a MODIS HDF file based on the filename
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
    Input:    string, filename w/ path (f)
    Output:   array (data)
    Features: Returns the data array of 0.05 deg monthly EVI for a given HDF 
              file 
    """
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
    """
    Name:     get_foh_grid
    Input:    None.
    Output:   list, list of tuples (foh_grid)
    Features: Returns a list of one hundred 0.05 degree lon-lat pairs based 
              on the class's half degree (HDG) coordinates (i.e., hdg_lon, 
              hdg_lat)
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
    Input:    list, list of lon-lat tuples (grid)
    Output:   list, list of tuples (foh_indices)
    Features: Returns a list of x-y indices based for a given list of 0.05 
              degree lon-lat pairs
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
    Input:    - float, 0.5 degree pixel longitude, degrees (lon)
              - float, 0.5 degree pixel latitude, degrees (lat)
              - numpy nd.array, 0.5 degree data (data)
    Output:   numpy nd.array (my_grid_evi)
    Features: Returns array of one hundred EVI values at 0.05 degree resolution 
              for a single 0.5 degree resolution pixel
    Depends:  - get_foh_grid
              - grid_to_index
    """
    my_grid_evi = numpy.array([])
    my_grid_pnts = get_foh_grid(lon, lat)
    my_grid_indx = grid_to_index(my_grid_pnts)
    for indx_pair in my_grid_indx:
        x,y = indx_pair
        zval = data[y][x]
        # Note: MODIS missing data is -3000
        # use NaN for invalid data (easy to remove for averaging):
        if zval < -2900:
            zval = numpy.NaN
        my_grid_evi = numpy.append(my_grid_evi, [zval])
    #
    return my_grid_evi

def upscale_evi(d):
    """
    Name:     upscale_evi
    Input:    numpy nd.array, 0.05 deg MODIS EVI (d)
    Output:   numpy nd.array (hdf_evi)
    Features: Returns an array of resampled 0.05 deg EVI to 0.5 deg resolution
              rounded to four decimal places
    Depends:  - get_lon_lat
              - get_foh_evi
    """
    # Initialize data array of floating points values:
    hdf_evi = numpy.zeros(shape=(360,720))
    #
    # Iterate through 0.5 deg lat-lon pairs:
    # --> y = latitude (0...359)
    # --> x = longitude (0...719)
    for y in xrange(360):
        for x in xrange(720):
            (lon, lat) = get_lon_lat(x, y, 0.5)
            #
            # Get 0.05 EVI values within this 0.5 cell:
            my_evi_data = get_foh_evi(lon, lat, d)
            #
            # Check that there's data in the array:
            if len(my_evi_data[~numpy.isnan(my_evi_data)]) > 0:
                # Calculate the average EVI for 0.5 pixel:
                ave_evi = my_evi_data[~numpy.isnan(my_evi_data)].mean()
                #
                # Scale result between 0 and 1:
                ave_evi = (ave_evi/10000.0)
                #
                # Adjust points outside threshold:
                if ave_evi < 0.0:
                    ave_evi = 0.0
                elif ave_evi > 1.0:
                    ave_evi = 1.0
                #
            else:
                # Assign global error value
                ave_evi = ERROR_VAL
                #
            hdf_evi[y,x] = ave_evi
    #
    # Round to four decimal places (the same as native MODIS EVI)
    numpy.around(hdf_evi, decimals=4, out=hdf_evi)
    #
    return hdf_evi

def nc_history():
    """
    Name:     nc_history
    Input:    None.
    Output:   string (my_str)
    Features: Returns a string for netCDF file history field based on the file's
              creation date
    """
    # f.history = 'created datetime()'
    my_str = "created %s" % datetime.date.today()
    return my_str

def nc_lat():
    """
    Name:     nc_lat
    Input:    None.
    Output:   numpy nd.array (my_lats)
    Features: Returns an array of latitudes from -90 to 90 at 0.5 deg resolution
    """
    my_lats = numpy.array([])
    x = 0
    for y in xrange(360):
        (lon, lat) = get_lon_lat(x,y,0.5)
        my_lats = numpy.append(my_lats, [lat,])
    return my_lats

def nc_lon():
    """
    Name:     nc_lon
    Input:    None.
    Output:   numpy nd.array (my_lons)
    Features: Returns an array of longitudes from -180 to 180 at 0.5 deg 
              resolution
    """
    my_lons = numpy.array([])
    y = 0
    for x in xrange(720):
        (lon, lat) = get_lon_lat(x,y,0.5)
        my_lons = numpy.append(my_lons, [lon,])
    return my_lons

def get_days(ts):
    """
    Name:     get_days
    Input:    datetime date (ts)
    Output:   int (delta.days)
    Features: Returns the number of days since 1 Jan 1860 for a given timestamp
    """
    # Set base time:
    base_time = datetime.date(1860, 1, 1)
    #
    # Calculate days since base time:
    delta = (ts - base_time)
    #
    # Return the number of days:
    return delta.days

##############################################################################
## MAIN
##############################################################################
# Open directory with HDF files and read file names:
#
# NOTE: change, 'mydir', 'myfiles' and the f.satellite assignments when
#       switching from Aqua to Terra MODIS data.
#  AQUA years: 2003--2011, 2013 (2012 is missing month of June)
#  TERRA years: 2001, 2002, 2012
#
mydir = (
    #"/home/user/Projects/gepisat/data/modis/aqua_evi/"
    "/Users/twdavis/Projects/data/modis/vi_cgm_monthly/terra/"
    )
myyears = (2001, 2002, 2012)
for year in myyears:
    #myfiles = glob.glob(mydir+'MYD13C2.A' + str(year) + '*.hdf') # aqua
    myfiles = glob.glob(mydir+'MOD13C2.A' + str(year) + '*.hdf') # terra
    
    # Create a new netCDF file for saving the output:
    nc_file = "%s%s%d%s.nc" % (mydir, 'ISI-MIP_', year, '-fAPAR_0.5deg')
    f = netcdf.netcdf_file(nc_file, 'w')
    #
    # Add meta data variables to file:
    f.contact1 = 'tyler.davis@imperial.ac.uk'
    f.contact2 = 'c.prentice@imperial.ac.uk'
    f.history = nc_history()
    f.institution = 'Imperial College London'
    f.note1 = (
        'EVI to fAPAR is based on Eq. 11 in Xiao et al., 2005 '
        '(Ecological Applications, vol. 15, no. 3, pp. 954-969)'
        )
    f.note2 = (
        'Upscaling of EVI from 0.05 to 0.5 degrees is based on the '
        'arithmetic mean of the one hundred 0.05 pixels with valid data'
        )
    f.note3 = (
        'Mean values for each 0.5 degree pixel are scaled between 0 and 1 '
        'and are rounded to four decimal places; '
        'any values greater than 1 are set equal to 1; '
        'any values less than 0 are set equal to 0; '
        'pixels without data are set equal to ' +
        str(ERROR_VAL)
        )
    #f.satellite = 'Aqua'
    f.satellite = 'Terra'
    f.title = 'Monthly 0.5 degree fAPAR based on MODIS EVI'
    
    # Create latitude dimension & variable:
    f.createDimension('lat', 360)
    lat = f.createVariable('lat', 'd', ('lat',) )
    lat[:] = nc_lat()
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat.axis = 'Y'
    #
    # Create longitude dimension & variable:
    f.createDimension('lon', 720)
    lon = f.createVariable('lon', 'd', ('lon',))
    lon[:] = nc_lon()
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    lon.axis = 'X'
    #
    # Create time dimension & variable
    f.createDimension('time', len(myfiles))
    t = f.createVariable('time', 'd', ('time',))
    t.standard_name = 'time'
    t.units = 'days since 1860-01-01 00:00:00'
    
    # Create fAPAR variable:
    fapar = f.createVariable('fAPAR', 'd', ('time','lat','lon'))
    fapar._FillValue = ERROR_VAL
    fapar.missing_value = ERROR_VAL
    fapar.long_name = (
        'Fractionally absorbed photosynthetically active radiation'
        )
    fapar.units = 'none'
    
    # Initialize time index:
    t_idx = 0
    
    # Iterate through all the HDF files (in order)
    for myfile in numpy.sort(myfiles):
        # Get time value for this file & update netCDF time variable:
        myts = get_ts(myfile)
        my_days = get_days(myts)
        t[t_idx] = my_days
        #
        # Open HDF file and get array of EVI data
        data = get_evi(myfile)
        #
        if data.any():
            # 
            # Save EVI to netCDF variable:
            fapar[t_idx] = upscale_evi(data)
            #
            # Update time index:
            t_idx += 1
    
    # Close and save netCDF file:
    f.close()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
if 0:
    #### Testing netCDF write ####
    nc_file = mydir + 'example.nc'
    f = netcdf.netcdf_file(nc_file, 'w')
    
    # Meta data:
    f.history = nc_history()
    f.contact1 = 'tyler.davis@imperial.ac.uk'
    f.contact2 = 'c.prentice@imperial.ac.uk'
    f.title = 'Monthly 0.5 degree fAPAR based on MODIS EVI'
    f.institution = 'Imperial College London'
    
    # Create latitude dimension & variable:
    f.createDimension('lat', 360)
    lat = f.createVariable('lat', 'd', ('lat',) )
    lat[:] = nc_lat()
    lat.units = 'degrees_north'
    
    # Create longitude dimension & variable:
    f.createDimension('lon', 720)
    lon = f.createVariable('lon', 'd', ('lon',))
    lon[:] = nc_lon()
    lon.units = 'degrees_east'
    
    # Close & save changes made to netCDF file:
    f.close()
    
    #### Testing netCDF read ####
    mydir = "/Users/twdavis/Projects/data/modis/fapar/"
    nc_file = mydir + "ISI-MIP_2001-fAPAR_0.5deg.nc"
    f = netcdf.NetCDFFile(nc_file, 'r')
    f.close()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# Save data to array:
#d_array = numpy.array(d_list)
#numpy.reshape(d_array, (3600,7200), order='C')

# Plot 2D graphic (grayscale):
#plt.imshow(numpy.reshape(d_array, (3600,7200), order='C'))
#plt.gray()
#plt.show()