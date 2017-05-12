#!/usr/bin/python
#
# utilities.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-05
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2017 Prentice Lab
#
# This file is part of the GePiSaT (Global ecosystem Production in Space and
# Time) model.
#
# GePiSaT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# GePiSaT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GePiSaT.  If not, see <http://www.gnu.org/licenses/>.
#
# ---------
# citation:
# ---------
# Davis, T.W., B.D. Stocker, X.M.P. Gilbert, T.F. Keenan, H. Wang, B.J. Evans,
# and I.C. Prentice. The Global ecosystem Production in Space and Time
# (GePiSaT) Model of the terrestrial biosphere: Part 1 - Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.

###############################################################################
# IMPORT MODULES
###############################################################################
import glob
import os
import logging

import numpy


###############################################################################
# FUNCTIONS
###############################################################################
def filter_grid_metadata(meta_path):
    """
    Filters the half-degree grid (HDG) meta data file for flux stations only
    """
    try:
        f = open(meta_path, 'r')
        header_line = f.readline()
        data = f.readlines()
    except:
        logging.error("Failed to read meta data from file %s", meta_path)
    else:
        f.close()

        # Prepare the output file:
        mf_base, mf_ext = os.path.splitext(meta_path)
        out_path = "%s_flx%s" % (mf_base, mf_ext)

        OUT = open(out_path, 'w')
        OUT.write(header_line)
        OUT.close()

        # Get list of flux station 0.5-degree grid points:
        station_list = get_station_latlon()

        # Find the indexes for lon/lat
        header_parts = header_line.rstrip("\n").split(",")
        lat_idx = header_parts.index("lat")
        lon_idx = header_parts.index("lon")

        # Read through all data and filter based on lon/lat
        for my_line in data:
            my_parts = my_line.rstrip("\n").split(",")
            pxl_lat = float(my_parts[lat_idx])
            pxl_lon = float(my_parts[lon_idx])
            if (pxl_lat, pxl_lon) in station_list:
                OUT = open(out_path, 'a')
                OUT.write(my_line)
                OUT.close()


def find_files(my_dir, my_pattern):
    """
    Name:     find_files
    Inputs:   - str, directory path (my_dir)
              - str, file name search pattern (my_pattern)
    Outputs:  list, file paths
    Features: Returns a sorted list of files found at a given directory with
              file names that match a given pattern
    """
    my_files = []
    if os.path.isdir(my_dir):
        s_str = os.path.join(my_dir, my_pattern)
        my_files = glob.glob(s_str)
    if len(my_files) > 0:
        my_files = sorted(my_files)
    return my_files


def get_station_latlon():
    """
    Returns a list of (lat, lon) pairs at 0.5-degree resolution (same as
    gridded datasets) for each flux station. Use this list to filter the
    name of gridded stations required in the GePiSaT database. Currently,
    the flux station metadata file is hard-coded: SUBJECT TO CHANGE BASED ON
    THE USER!

    TODO: make the FLUXNET metadata file a global variable that is easy for
          users to set before processing
    """
    grid_points = []

    base_dir = os.path.expanduser("~")
    data_dir = os.path.join(base_dir, "Data", "psql")
    met_file = os.path.join(data_dir, "Fluxdata_Met-Data.csv")
    if os.path.isfile(met_file):
        my_data = numpy.loadtxt(
                fname=met_file,
                dtype={'names': ('station', 'lat', 'lon'),
                       'formats': ('S6', 'f4', 'f4')},
                delimiter=',',
                usecols=(4, 6, 7),
                skiprows=1)

        for line in my_data:
            flux_site, flux_lat, flux_lon = line
            grid_lon, grid_lat = grid_centroid(flux_lon, flux_lat, 0.5)
            grid_points.append((grid_lat, grid_lon))

        grid_points = sorted(list(set(grid_points)))

    return grid_points


def grid_centroid(my_lon, my_lat, grid_res=0.5):
    """
    Name:     grid_centroid
    Input:    - float, longitude (my_lon)
              - float, latitude (my_lat)
    Output:   tuple, longitude latitude pair (my_centroid)
    Features: Returns the nearest 0.5 deg. grid centroid per given coordinates
              based on the Euclidean distance to each of the four surrounding
              grids; if any distances are equivalent, the pixel north and east
              is selected by default
    """
    # Create lists of regular latitude and longitude:
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = int(180./grid_res)
    lon_dim = int(360./grid_res)
    lats = [lat_min + y*grid_res for y in range(lat_dim)]
    lons = [lon_min + x*grid_res for x in range(lon_dim)]

    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index - 1]
        try:
            bb_lon_max = lons[lon_index + 1]
        except IndexError:
            bb_lon_max = lons[-1] + grid_res

    # Find bounding latitude:
    centroid_lat = None
    if my_lat in lats:
        centroid_lat = my_lat
    else:
        lats.append(my_lat)
        lats.sort()
        lat_index = lats.index(my_lat)
        bb_lat_min = lats[lat_index - 1]
        try:
            bb_lat_max = lats[lat_index + 1]
        except IndexError:
            bb_lat_max = lats[-1] + grid_res

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
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0)
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0)
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

    # Return nearest centroid:
    return my_centroid


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
        logging.error("Error: cannot write to file: %s", f)
    else:
        OUT.close()
