#!/usr/bin/python
#
# glasdata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-01-13
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
# (GePiSaT) Model of the terrestrial biosphere: Part 1 — Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.

###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import glob
import logging

import numpy
import PIL.Image as Image

from .utilities import writeout
from .var import VAR


###############################################################################
# CLASSES
###############################################################################
class GLASDATA:
    """
    Name:     GLASDATA
    Features: This class processes 1km TIF rasters to 0.5 deg resolution data
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    variable_file = ""
    dataset_file = ""
    data = None

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, vf, df, d):
        """
        Name:     GLASDATA.__init__
        Input:    - string, variable list output file with path (vf)
                  - string, data set output file with path (df)
                  - PIL sequence, image data (d)
        """
        # Create a class logger
        self.logger = logging.getLogger("GLASDATA")
        self.logger.info("GLASDATA class initialized")

        self.variable_file = vf
        self.dataset_file = df
        self.data = d

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def process_data(self):
        """
        Name:     GLASDATA.process_data
        Input:    None.
        Output:   None.
        Features: Writes the 0.5 degree averaged image data to variable list
                  and data set files
        Depends:  - get_lon_lat
                  - get_station
                  - get_1km_rh
        """
        # Iterate through spatial data:
        for y in range(360):
            for x in range(720):
                # Determine lon, lat, and station ID:
                lon, lat = self.get_lon_lat(x, y, 0.5)
                st_id = self.get_stationid(lon, lat)

                # Create variable class:
                station_parts = ('HDG', st_id)
                my_line = VAR(station_parts, 'RH100', 'grid')
                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                my_line.printLine(self.variable_file)
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                # Retrieve 1km data for this 0.5 deg pixel:
                my_rh_data = self.get_1km_rh(lon, lat)

                # Check that there's data in the array:
                if len(my_rh_data) > 0:
                    val = my_rh_data.mean()
                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Append to existing file:
                OUT = open(self.dataset_file, 'a')

                # Create/write output line:
                time_stamp = datetime.date(2005, 1, 1)
                outline = "%s,%s,%s,%0.2f\n" % (
                    my_line.msvIDX,
                    my_line.stationID,
                    time_stamp,
                    val
                    )
                OUT.write(outline)
                OUT.close()

    def get_1km_grid(self, hdg_lon, hdg_lat):
        """
        Name:     GLASDATA.get_1km_grid
        Input:    - float, 0.5 degree longitude (hdg_lon)
                  - float, 0.5 degree latitude (hdg_lat)
        Output:   list, list of tuples (okm_grid)
        Features: Returns a list of 3600 lon-lat pairs at 1km resolution for a
                  given 0.5 degree lon-lat pair
        """
        # Initialize five one-hundreths grid:
        okm_grid = []

        # Define half-degree and five one-hundredths resolutions:
        hdg_res = 0.5
        okm_res = 1.0/120.0

        # Calculate the binding box at 0.5 deg:
        westing = hdg_lon - 0.5*hdg_res
        northing = hdg_lat + 0.5*hdg_res

        # Initialize centroid offsetting for 1 km grid:
        okm_lon = westing + 0.5*okm_res
        okm_lat = northing - 0.5*okm_res

        # Iterate over the 60x60 box:
        for y in range(60):
            lat = okm_lat - y*okm_res
            for x in range(60):
                lon = okm_lon + x*okm_res
                okm_grid.append((lon, lat))

        return okm_grid

    def get_1km_rh(self, lon, lat):
        """
        Name:     GLASDATA.get_1km_rh
        Input:    - float, 0.5 degree longitude (lon)
                  - float, 0.5 degree latitude (lat)
        Output:   numpy nd.array
        Features: Returns array of 3600 canopy height (RH) values at 1 km
                  resolution for a single 0.5 pixel
        Depends:  - get_1km_grid
                  - grid_to_index
        """
        my_grid_rh = numpy.array([])
        my_grid_pnts = self.get_1km_grid(lon, lat)
        my_grid_indx = self.grid_to_index(my_grid_pnts)
        for indx_pair in my_grid_indx:
            x, y = indx_pair
            zval = self.data.getpixel((x, y))
            my_grid_rh = numpy.append(my_grid_rh, [zval])

        return my_grid_rh

    def get_lon_lat(self, x, y, r):
        """
        Name:     GLASDATA.get_lon_lat
        Input:    - int, longitude index (x)
                  - int, latitude index (y)
                  - float, pixel resolution (r)
        Output:   tuple, longitude and latitude pair, degrees
        Features: Returns lon-lat pair for an x-y index pair (numbered from the
                  bottom-left corner) and pixel resolution
        """
        # Offset lat, lon to pixel centroid
        lon = -180.0 + 0.5*r
        lat = 90.0 - 0.5*r

        # Offset lat, lon based on pixel index
        lon = lon + x*r
        lat = lat - y*r
        return (lon, lat)

    def get_stationid(self, lon, lat):
        """
        Name:     GLASDATA.get_stationid
        Input:    - float, 0.5 degree longitude (lon)
                  - float, 0.5 degree latitude (lat)
        Output:   int, station id (st_id)
        Features: Returns the half-degree (HDG) station ID for a pixel
                  numbered from 0 (bottom-left / south-west corner) to 259199
                  (top-right / north-east corner) as defined in the GePiSaT
                  database numbering scheme
        """
        # Station ID is based on 0 being the bottom (south) left (west) corner
        # and 259199 being the top (north) right (east) corner as is used in
        # the postgreSQL database naming scheme.
        st_id = (
            720.0 * (359.0 - ((90.0 - lat)/0.5 - 0.5)) +
            ((lon + 180.0)/0.5 - 0.5))
        return st_id

    def get_x_y(self, lon, lat, r):
        """
        Name:     GLASDATA.get_x_y
        Input:    - float, longitude (lon)
                  - float, latitude (lat)
                  - float, resolution (r)
        Output:   tuple, x-y indices
        Features: Returns x and y indices for a given lon-lat pair and pixel
                  resolution
        """
        # Solve x and y indices:
        x = (lon + 180.0)/r - 0.5
        y = (90.0 - lat)/r - 0.5

        return (int(x), int(y))

    def grid_to_index(self, grid):
        """
        Name:     GLASDATA.grid_to_index
        Input:    list, list of tuples (grid)
        Output:   list, list of tuples (okm_indices)
        Features: Returns a list of x-y indices based on a list of lon-lat
                  pairs for a 1 km grid
        Depends:  get_x_y
        """
        okm_indices = []
        my_res = 1.0/120.0
        for grid_pair in grid:
            lon, lat = grid_pair
            x, y = self.get_x_y(lon, lat, my_res)
            okm_indices.append((x, y))

        return okm_indices


###############################################################################
# FUNCTIONS
###############################################################################
def process_glas(d):
    """
    Name:     process_glas
    Input:    string, input/output file directory (d)
    Features: Processes GLAS TIF raster data into 0.5 degree canopy height and
              saves variable list and data set table output files
    Depends:  writeout
    """
    # Read TIF raster file in the file directory:
    my_file = glob.glob(d + "*.tif")[0]

    # Open file for reading and save variable of interest:
    im = Image.open(my_file).getdata()

    # Save the shape values of each for iteration purposes
    sh_lat, sh_lon = im.size

    # Prepare var output file:
    my_var_out = "GLAS_Var-List_rh100.csv"
    var_outfile = d + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)

    # Prepare data set output file:
    my_dat_out = "GLAS_Data-Set_rh100_2005.csv"
    dat_outfile = d + my_dat_out
    dat_headerline = "msvidx,stationid,datetime,data\n"
    writeout(dat_outfile, dat_headerline)

    # Initialize and run GLAS class:
    my_glas = GLASDATA(var_outfile, dat_outfile, im)
    my_glas.process_data()
