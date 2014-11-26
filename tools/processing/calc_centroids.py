#!/usr/bin/python
#
# reg_points.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-08-28 -- created
# 2013-08-28 -- last updated
#
# ------------
# description:
# ------------
# This script calculates the centroids for a regular grid.
#
# ----------
# changelog:
# ----------
# 00. created [13.08.28]
# 01. updated comments [14.02.12]
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import os.path

###############################################################################
## DEFINE FUNCTIONS:
###############################################################################
def calc_centroids(gv, ev, out):
    # gv :: grid values tuple (ncol, nrow, resolution)
    # ev :: extent values tuple (northing, southing, westing, easting)
    # out :: output file (with full path)
    #
    # Create an output header:
    header = "id,lat,lon\n"
    #
    # Assign values:
    x_tot, y_tot, dx = gv
    y_max, y_min, x_min, x_max = ev
    #
    # Process centroids:
    # cat(0) = top left corner (x_min, y_max)
    for y in xrange(y_tot):
        # Latitude changes for each y (offset):
        lat = (y_max - 0.5*dx) - y*dx
        if y == 0:
            try:
                f = open(out, "w")
            except IOError:
                print "Cannot open", out, "for writing"
            else:
                f.write(header)
                f.close()
        for x in xrange(x_tot):
            # Row-major category ID:
            cat = x_tot*y + x
            # Longitude changes for each x (offset):
            lon = (x_min + 0.5*dx) + x*dx
            # Print out:
            try:
                f = open(out, "a")
            except IOError:
                # If file writing not available, print to stdout:
                print cat, lat, lon
            else:
                f.write("%d,%0.3f,%0.3f\n" % (cat,lat,lon))
                f.close()

###############################################################################
## DEFINE VARIABLES:
###############################################################################
# Define number of cells (x and y):
grid_cells_x = 720
grid_cells_y = 360

# Define resolution:
grid_resolution = 0.05

# Define spatial extents:
ext_north = 90.0
ext_south = -90.0
ext_west = -180.0
ext_east = 180.0

# Define output file and directory:
output_dir = "/home/user/Desktop/"
output_file = "wfdei_grid.csv"

###############################################################################
## MAIN:
###############################################################################
# Calculate grid cell centroids:
grid_vals = (grid_cells_x, grid_cells_y, grid_resolution)
ext_vals = (ext_north, ext_south, ext_west, ext_east)
output_location = output_dir + output_file
calc_centroids(grid_vals, ext_vals, output_location)

