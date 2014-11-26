#!/usr/bin/python
#
#
# grid_centroid.py
#
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-09-09 -- created
# 2013-09-13 -- last updated
#
# ------------
# description:
# ------------
# Given a point location, i.e., lon and lat, this script finds the four grid
# centroids surrounding the point, calculates the distances between the point 
# and centroids, and returns the closest centroid. If calculated distances are
# equal, the nearest centroid defaults north and east.
#
#  lat (-90,90)
#  |
#  .-------.-------.          N
#  |       |       |          |
#  |   C   |   A   |      W --.-- E
#  |       | x     |          |
#  .---y---z-------.          S
#  |       |       |
#  |   D   |   B   |
#  |       |       |
#  .-------.-------.--> lon (-180,180)
#
#  [A-D]: grid centroids of bounding box around given lat:lon pair
#  x: some given lon:lat pair
#
#  Examples:
#  - case x, centroid A is the nearest to the given point
#  - case y, centroid C is defaulted to as nearest
#  - case z, centroid A is defaulted to as nearest
#
# ----------
# changelog:
# ----------
# 01. regular grid variables (resolution, SW coordinate, number of pixels) have
#     been abstracted so grids of different resolutions can be easier 
#     implemented if necessary [13.09.09]
# 02. changed function name from gridwork to grid_centroid [13.09.13]
#
################################################################################
## MODULES:
################################################################################
import numpy

################################################################################
## FUNCTIONS:
################################################################################
def grid_centroid(my_lon, my_lat):
    #
    # Create lists of regular latitude and longitude:
    grid_res = 0.5
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = 360
    lon_dim = 720
    lats = [lat_min + y * grid_res for y in xrange(lat_dim)]
    lons = [lon_min + x * grid_res for x in xrange(lon_dim)]
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
    # NOTE: if dist_A equals dist_B, then centroid defaults positively (north / east)
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
        dist_A = numpy.sqrt((bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_B = numpy.sqrt((bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0)
        dist_C = numpy.sqrt((my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_D = numpy.sqrt((my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0)
        min_dist = min([dist_A, dist_B, dist_C, dist_D])
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

################################################################################
## MAIN:
################################################################################
some_lon = 0.0
some_lat = 0.0
grid_lon, grid_lat = gridwork(some_lon, some_lat)
print grid_lon, grid_lat
