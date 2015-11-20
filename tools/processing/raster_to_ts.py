#!/usr/bin/python
#
# raster_to_ts.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-09-23 -- created
# 2014-09-29 -- last updated
#
# ------------
# description:
# ------------
# This script creates a time series based on a single pixel (or list of pixels)
# read from monthly ASCII raster at 0.5 degree resolution.
#
# ----------
# changelog:
# ----------
# 01. update for daily solar radiation (Ho_day) analysis [14.09.29]
# 02. update for daily time series datetime objects [14.09.29]
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import glob
import numpy
import re

###############################################################################
## FUNCTIONS 
###############################################################################
def add_one_year(dt0):
    """
    Name:     add_one_year
    Input:    datetime date
    Output:   datetime date
    Features: Adds one year to the datetime preserving calendar date, if it 
              exists, otherwise uses the following day (i.e., February 29 
              becomes March 1)
    Ref:      G. Rees (2013) Stack Overflow
              http://stackoverflow.com/questions/15741618/add-one-year-in-
              current-date-python
    """
    try:
        return dt0.replace(year = dt0.year + 1)
    except ValueError:
        return dt0 + (
            datetime.date(dt0.year + 1, 1, 1) - datetime.date(dt0.year, 1 ,1)
        )

def get_year_days(ts):
    """
    Name:     get_year_days
    Input:    datetime date
    Output:   int
    Features: Returns the total number of days in the year
    Depends:  add_one_year
    """
    ts1 = datetime.date(ts.year, 1 , 1)
    ts2 = add_one_year(ts1)
    return (ts2 - ts1).days

def grid_centroid(my_lon, my_lat):
    """
    Name:     grid_centroid
    Input:    - float, longitude (my_lon)
              - float, latitude (my_lat)
    Output:   tuple, longitude latitude pair
    Features: Returns the nearest 0.5 deg. grid centroid per given coordinates
              based on the Euclidean distance to each of the four surrounding 
              grids; if any distances are equivalent, the pixel north and east
              is selected by default
    """
    # Create lists of regular latitude and longitude:
    grid_res = 0.5
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = 360
    lon_dim = 720
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
    # NOTE: if dist_A equals dist_B, then centroid defaults positively 
    #       i.e., north / east
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
        dist_A = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0
            )
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

def get_x_y(lon, lat, r):
    """
    Name:     get_x_y
    Input:    - float, longitude (lon)
              - float, latitude (lat)
              - float, resolution (r)
    Output:   tuple, x-y indices
    Features: Returns array indices for a given lon-lat pair and pixel 
              resolution
    """
    # Solve x and y indices:
    x = (lon + 180.0)/r - 0.5
    y = (90.0 - lat)/r - 0.5
    #
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

def save_to_file(f, d, v):
    """
    Name:     save_to_file
    Input:    - str, file name with path (f)
              - dict, time:data (d)
              - str, variable of interest (v)
    Output:   None.
    Features: Writes time series data to file
    Depends:  - get_year_days
              - writeout
    """
    # Define variable units
    units = {
        'CPA' : 'unitless',
        'CWD' : 'mm',
        'EET' : 'mm',
        'PET' : 'mm',
        'PPFD' : 'mol.m2',
        'Ho_day' : 'kJ.m2'
    }
    #
    # Create file with header:
    header = "Timestamp,Year," + v + "_" + units[v] + "\n"
    writeout(f, header)
    #
    for k in numpy.sort(d.keys()):
        ny = get_year_days(k)
        dd = k.year + (k.timetuple().tm_yday - 1.0)/ny
        out_line = "%s,%0.6f,%0.3f\n" % (k, dd, d[k])
        try:
            OUT = open(f, 'a')
            OUT.write(out_line)
        except IOError:
            print "Could not append to file", f
        else:
            OUT.close()

###############################################################################
## DEFINITIONS 
###############################################################################
mac = 0   # 1: TRUE, 0: FALSE
if mac:
    # Mac:
    analysis_dirs = []
else:
    # Linux:
    analysis_dirs = ["/home/user/Projects/gepisat/data/stash/an1/",
                     "/home/user/Projects/gepisat/data/stash/an2/",
                     "/home/user/Projects/gepisat/data/stash/an3/",
                     "/home/user/Projects/gepisat/data/stash/an4/",
                     "/home/user/Projects/gepisat/data/stash/an5/",
                     "/home/user/Projects/gepisat/data/stash/an6/"]
    output_dir = "/home/user/Projects/gepisat/data/stash/"

# Silwood Park Campus:     (51.4078,-0.6405)
# Imperial College London: (51.4974,-0.1752)
(st_lat, st_lon) = (51.4078, -0.6405)

# Get pixel information:
(px_lon, px_lat) = grid_centroid(st_lon, st_lat)
(px_x, px_y) = get_x_y(px_lon, px_lat, 0.5)

# Set analysis and variable flags:
voi = "Ho_day" # CPA, CWD, EET, PET, PPFD, Ho_day

###############################################################################
## MAIN PROGRAM 
###############################################################################
for analysis in [i+1 for i in xrange(6)]:
    # Read files for processing:
    my_files = glob.glob(analysis_dirs[analysis-1] + "*" + voi + "*.txt")
    #
    if my_files:
        # Initialize/reset time series data dictionary:
        px_data = {}
        #
        for my_file in numpy.sort(my_files):
            try:
                # Read data and search for date in filename:
                my_date = re.search('\d{4}-\d{2}-\d{2}', my_file).group(0)
                my_data = numpy.loadtxt(my_file, skiprows=6)
            except:
                print "Issues were encounted with file", my_file
            else:
                # Create a datetime object for dictionary:
                temp_yr, temp_mo, temp_dy = my_date.split('-')
                this_year = int(temp_yr)
                this_month = int(temp_mo)
                this_day = int(temp_dy)
                my_ts = datetime.date(this_year, this_month, this_day)
                # Extract value from raster and scale (div by 1000.0)
                px_val = my_data[px_y, px_x]
                #px_data[my_ts] = (px_val/1000.0)
                # For Ho_day, raster is already in units of kJ/m2
                px_data[my_ts] = (px_val)
        #
        # Save data to file:
        output_file = (
            output_dir + "Analysis-" + str(analysis) + "_" + voi + "_" + 
            str(px_lat) + "_" + str(px_lon) + ".txt"
        )
        save_to_file(output_file, px_data, voi)