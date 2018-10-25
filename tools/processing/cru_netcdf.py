#!/usr/bin/python
#
# cru_netcdf.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-11-05 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script reads NetCDF files (e.g., CRU TS 3.21) and extracts variables of
# interest (e.g., monthly average max and min air temperature, monthly average
# vapor pressure). Writes variable of interest to ASCII raster format.
#
# NOTE: CRU TS 3.21 is made of 0.5 deg. resolution pixels
#       Indexes begin at -179.75, -89.75
#       Indexes end at 179.75, 89.75
#
#       Example files:
#       1. cru_ts3.21.1901.2012.cld.dat.nc, cloudiness, %
#       2. cru_ts3.21.1901.2012.pre.dat.nc, precipitation, mm
#       3. cru_ts3.21.1901.2012.tmn.dat.nc, min. temperature, deg. C
#       4. cru_ts3.21.1901.2012.tmp.dat.nc, mean temperature, deg. C
#       5. cru_ts3.21.1901.2012.tmx.dat.nc, max. temperature, deg. C
#       6. cru_ts3.21.1901.2012.vap.dat.nc, vapor pressure, hPa
#
# NOTE 2: CRU TS3.21 data is from January 1901 to December 2012
#         Only interested in data from 1990--2012
#
# ----------
# changelog:
# ----------
# 00. created [13.11.05]
# 01. edited vpd calculation to 1990--2012 only [13.11.06]
# 02. updated process raster time_stamp [13.11.06]
# 03. added process_poly() function [13.11.21]
# 04. added get_monthly_cru & get_time_index() [14.02.18]
# 05. updated calc_vpd, process_raster, and process_poly [14.02.18]
# 06. changed the way of processing vpd (month-at-a-time) [14.02.18]
# 07. added add_years() function [14.02.18]
# 08. added get_clip & get_missing functions [14.02.18]
# 09. implemented mean monthly calcs for vpd [14.02.18]
# 10. fixed error with VPD calc [14.02.18]
# --> 237.3 not 273.
# 11. updated function documentation [14.09.17]
# 12. code housekeeping [14.12.01]
#
# -----
# todo:
# -----
# 1. Check stationid properness by implementing process_poly
# 2. Make distinction between get_lon_lat function for top-left and bottom-left
#    numbering (e.g., MODIS HDF versus WATCH netCDF)
#
###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import os

import numpy
from scipy.io import netcdf

###############################################################################
## FUNCTIONS
###############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime date
    Output:   datetime date
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3

def add_years(d, years):
    """
    Name:     add_years
    Input:    - datetime date (d)
              - int (years)
    Output:   datetime date
    Features: Adds years to the datetime preserving calendar date, if it
              exists, otherwise uses the following day (i.e., February 29
              becomes March 1)
    Ref:      G. Rees (2013) Stack Overflow
              http://stackoverflow.com/questions/15741618/add-one-year-in-
              current-date-python
    """
    d_new = None
    try:
        d_new = d.replace(year = d.year + years)
    except ValueError:
        d_new = d + (
            datetime.date(d.year + years, 1, 1) - datetime.date(d.year, 1, 1)
            )
    finally:
        return d_new

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

def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime date, base timestamp (bt)
              - datetime date, current timestamp to be found (ct)
              - numpy nd.array, days since base timestamp (aot)
    Output:   int
    Features: Finds the index in an array of CRU TS days for a given timestamp
    """
    # For CRU TS 3.21, the aot is indexed for mid-month days, e.g. 15--16th
    # therefore, to make certain that ct index preceeds the index for the
    # correct month in aot, make the day of the current month less than
    # the 15th or 16th (i.e., replace day with '1'):
    ct = ct.replace(day=1)
    #
    # Calculate the time difference between ct and bt:
    dt = (ct - bt).days
    #
    # Append dt to the aot array:
    aot = numpy.append(aot, [dt,])
    #
    # Find the first index of dt in the sorted array:
    idx = numpy.where(numpy.sort(aot)==dt)[0][0]
    return idx

def get_monthly_cru(d, ct, v):
    """
    Name:     get_monthly_cru
    Input:    - str, directory to CRU netcdf file (d)
              - datetime date, current month datetime object (ct)
              - str, variable of interest (v)
    Output:   numpy nd.array
    Depends:  get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and
              variable of interest (e.g., cld, pre, tmp)
    """
    # Search directory for netCDF file:
    my_file = glob.glob(d + "*" + v + ".dat.nc")[0]
    #
    if my_file:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")
        #
        # Save data for variables of interest:
        # NOTE: for CRU TS 3.21:
        #       variables: 'lat', 'lon', 'time', v
        #       where v is 'tmp', 'pre', 'cld', 'vap'
        # LAT:  -89.75 -- 89.75
        # LON:  -179.75 -- 179.75
        # TIME:
        #       units: days since 1900-1-1
        #       shape: (1344,)
        #       values: mid-day of each month (e.g., 15th or 16th)
        # DATA:
        #       'cld' units = %
        #       'pre' units = mm
        #       'tmp' units = deg. C
        #       'vap' units = hPa
        #       Missing value = 9.96e+36
        # Save the base time stamp:
        bt = datetime.date(1900,1,1)
        #
        # Read the time data as array:
        f_time = f.variables['time'].data
        #
        # Find the time index for the current date:
        ti = get_time_index(bt, ct, f_time)
        #
        # Get the spatial data for current time:
        f_data = f.variables[v].data[ti]
        f.close()
        return f_data

def calculate_vpd(tmin, tmax, vap):
    """
    Name:     calculate_vpd
    Input:    - numpy nd.array, mean monthly min daily air temp, deg C (tmin)
              - numpy nd.array, mean monthly max daily air temp, deg C (tmax)
              - numpy nd.array, mean monthly vapor pressure, hPa (vap)
    Output:   numpy nd.array, mean monthly vapor pressure deficit, kPa (vpd)
    Features: Returns an array of mean monthly vapor pressure deficit
    Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure
              Calculation Methods, in Evaporation and Evapotranspiration:
              Measurements and Estimations, Springer, London.
                vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
                where:
                    tc = average daily air temperature, deg C
                    ea = actual vapor pressure, kPa
    """
    # Initialize array:
    # NOTE: maintains large number for missing values
    vpd = 1e7*numpy.ones((360, 720))
    #
    # Iterate through each data point:
    lat, lon = tmin.shape
    #
    for y in xrange(lat):
        for x in xrange(lon):
            tm = tmin[y,x]
            tx = tmax[y,x]
            ea = vap[y,x]
            if tm < 1.e6 and tx < 1.e6 and ea < 1.e6:
                to = 0.5*(tm + tx)
                vpd[y,x] = (0.611*numpy.exp((17.27*to)/(to + 237.3)) - 0.10*ea)
    return vpd

def process_raster(data, var_name, out_dir, time_stamp):
    """
    Name:     process_raster
    Input:    - numpy nd.array (data)
              - string, variable name for output file (var_name)
              - string, output directory w/ path (out_dir)
              - datetime date, time stamp of data (time_stamp)
    Output:   None.
    Features: Saves monthly data in ASCII raster format (int value x100)
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
    lat, lon = data.shape
    #
    # Create output file name:
    rast_out = "%s%s_%s_%s_%s.txt" % (
        out_dir,
        "CRU-TS",
        var_name,
        "0.5-Raster",
        time_stamp
        )
    #
    # Save header line to file:
    writeout(rast_out, header)
    #
    # Iterate through current month's data:
    for i in xrange(lat):
        # Reverse column ordering:
        y = 359 - i
        #
        # Reset outline for next row:
        outline = ""
        #
        for x in xrange(lon):
            val = data[y,x]
            if val < 1.e6:
                val = int(val * 100)
            else:
                val = -9999
            # Add value to outline:
            outline = "%s%d " % (outline, val)
        # End line and print to file:
        outline = "%s\n" % outline.rstrip(' ')
        OUT = open(rast_out, 'a')
        OUT.write(outline)
        OUT.close()

def process_poly(d, var_name, data, time_stamp):
    """
    Name:     process_poly
    Input:    - string, output file w/ path (d)
              - string, variable name for output file (var_name)
              - numpy nd.array (data)
              - datetime date, time stamp of data (time_stamp)
    Output:   None.
    Features: Saves data in polygon (shapefile) CSV format
    """
    # Open and write header line for CSV file:
    # Save the shape values of each for iteration purposes
    sh_lat, sh_lon = data.shape
    #
    # Create output file name:
    poly_out = "%s%s_%s_%s_%s.csv" % (
        d,
        "CRU-TS",
        var_name,
        "0.5-Poly",
        time_stamp
        )
    #
    # Save header line to file:
    header = "id,lon,lat,vpd\n"
    writeout(poly_out, header)
    #
    # Iterate through 0.5 deg lat-lon pairs:
    # --> y = latitude (0...359)
    # --> x = longitude (0...719)
    for y in xrange(sh_lat):
        for x in xrange(sh_lon):
            val = data[y,x]
            if val < 1.e6:
                # Save station info:
                (lon, lat) = get_lon_lat(x, y, 0.5)
                stid = get_stationid(lon, lat)
                #
                # Save VPD to file:
                OUT = open(poly_out, 'a')
                outline = "%d,%0.3f,%0.3f,%0.5f\n" % (
                    stid, lon, lat, val
                    )
                OUT.write(outline)
                OUT.close()

def get_lon_lat(x,y,r):
    """
    Name:     get_lon_lat
    Input:    - int/nd.array, longitude index (x)
              - int/nd.array, latitude index (y)
              - float, pixel resolution (r)
    Output:   float/nd.array tuple, longitude(s) and latitude(s), degrees
    Features: Returns lon-lat pair for an x-y index pair (numbered from the
              bottom-left corner) and pixel resolution
    """
    # Offset lat, lon to pixel centroid
    lon = -180.0 + (0.5*r)
    lat = -90.0 + (0.5*r)
    #
    # Offset lat, lon based on pixel index
    lon = lon + (x*r)
    lat = lat + (y*r)
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
    st_id = (
        720.0*(359.0 - ((90.0 - lat)/0.5 - 0.5)) + ((lon + 180.0)/0.5 - 0.5)
    )
    return int(st_id)

def get_clip(d):
    """
    Name:     get_clip
    Input:    numpy nd.array (d)
    Output:   numpy nd.array (temp_array)
    Features: Returns an array defining land (1) and ocean (0) pixels based on
              CRU TS data (ocean pixels have error value ~1.e+36)
    """
    temp_array = d - d.max()
    temp_array = temp_array.clip(min=-1)
    temp_array = -1.0*temp_array
    return temp_array

def get_missing(d):
    """
    Name:     get_missing
    Input:    numpy nd.array, an output from the get_clip function (d)
    Output:   numpy nd.array (temp_array)
    Features: Returns an array where ocean pixels (based on the output of the
              get_clip function) are set equal to error value -9999.0"""
    temp_array = 9999.0*(d - 1.0)
    return temp_array

###############################################################################
## CONSTANTS
###############################################################################
# Data directory:
file_directory = (
    #"/Users/twdavis/Projects/data/cru/cru_ts_3_21/"
    "/usr/local/share/database/cru/"
    )
output_directory = (
    #"/Users/twdavis/Projects/data/cru/cru_ts_3_21/vpd/"
    "/home/user/Desktop/temp/"
    )

# Start and end dates to process (25 years):
start_date = datetime.date(1988,12,1)
end_date = datetime.date(2013,1,1)

###############################################################################
## MAIN PROGRAM
###############################################################################
# Initialize mean monthly VPD field (based on multiple years):
vpd_mean = numpy.zeros((360, 720))

# Process each month between start and end dates:
cur_date = start_date
while cur_date < end_date:
    # Open and ead netcdf files in the file directory:
    tmn = get_monthly_cru(file_directory, cur_date, "tmn")
    tmx = get_monthly_cru(file_directory, cur_date, "tmx")
    vap = get_monthly_cru(file_directory, cur_date, "vap")
    #
    # Get the clip and error fields:
    clip_field = get_clip(tmn)
    error_field = get_missing(clip_field)
    #
    # Calculate VPD:
    vpd = calculate_vpd(tmn, tmx, vap)
    #
    # Add new vpd to mean:
    vpd_mean = vpd_mean + vpd*clip_field
    #
    # Process data into raster:
    #process_raster(vpd, "vpd", output_directory, cur_date)
    #process_poly(output_directory, "VPD", vpd, cur_date)
    #
    # Update date field:
    #cur_date = add_one_month(cur_date)  # <- for monthly processing
    cur_date = add_years(cur_date, 1)  # <- for annual means
#
# Compute mean vpd:
vpd_mean = (0.04)*vpd_mean

# Set error values (need to be >1e6 for raster):
vpd_mean = vpd_mean + -1000.0*error_field

# Process data into raster:
process_raster(vpd_mean, "vpd", output_directory, start_date)
