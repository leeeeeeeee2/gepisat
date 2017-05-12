#!/usr/bin/python
#
# modisdata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-12
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
import datetime
import logging
import os
import re

import numpy
from pyhdf import SD   # also provided by python-hdf4 (e.g., Linux Mint)

from .utilities import find_files
from .utilities import get_station_latlon
from .utilities import writeout
from .var import VAR


###############################################################################
# CLASSES
###############################################################################
class MODISDATA:
    """
    Name:     MODISDATA
    Features: This class processes MODIS EVI to FAPAR for the data_set database
              table
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    hdg_id = -9999      # 0.5 deg station id
    hdg_lon = -999.0    # 0.5 deg longitude
    hdg_lat = -99.00    # 0.5 deg latitude
    hdg_ave = -3000     # 0.5 deg ave pxl. value (based on 0.05 data)
    hdg_res = 0.5       # 0.5 deg pixel resolution
    foh_res = 0.05      # 0.05 deg pixel resolution
    station_id = 0      # station id
    msv_idx = ''        # msvidx
    data_time = datetime.date(1900, 1, 1)  # data time
    data_val = -9999.0  # EVI (pixel / 10000)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, x, y, t, d):
        """
        Name:     MODISDATA.__init__
        Input:    - int, latitude index, e.g., 0--359 (x)
                  - int, longitude index, e.g., 0--719 (y)
                  - datetime.date, timestamp of observation (t)
                  - numpy nd.array, MODIS data, e.g., 3600x7200 (d)
        """
        # Create a class logger
        self.logger = logging.getLogger("MODISDATA")
        self.logger.debug("MODISDATA class initialized")

        # Set datetime:
        self.data_time = t

        # Get longitude and latitude from indices:
        (self.hdg_lon, self.hdg_lat) = self.get_lon_lat(x, y, self.hdg_res)

        # Get HDG station ID:
        self.hdg_id = self.get_stationid()

        # Get station id and msvidx:
        station_parts = ('HDG', self.hdg_id)
        my_var = VAR(station_parts, 'FAPAR', 'grid')
        self.station_id = my_var.stationID
        self.msv_idx = my_var.msvIDX

        # Get 0.05 res data for this 0.5 pixel:
        my_points = self.get_foh_points(self.hdg_lon, self.hdg_lat, d)

        # Check that there's data in the array:
        if len(my_points[~numpy.isnan(my_points)]) > 0:
            # Calculate the average EVI for 0.5 pixel:
            self.hdg_ave = my_points[~numpy.isnan(my_points)].mean()
            self.data_val = self.hdg_ave / 10000.0
        else:
            # Use the pre-defined error value (i.e., -3000)
            self.hdg_ave = -3000
            self.data_val = -9999.0

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def get_lon_lat(self, x, y, r):
        """
        Name:     MODISDATA.get_lon_lat
        Input:    - int, latitude index, e.g., 0--319 (x)
                  - int, longitude index, e.g., 0--719 (y)
                  - float, resolution, e.g., 0.5 (r)
        Output:   tuple, longitude-latitude pair (lon, lat)
        Features: Returns lon-lat pair for an x-y index pair (numbered from the
                  bottom-left corner) and pixel resolution
        """
        # Offset lat, lon to pixel centroid
        lon = -180.0 + 0.5*r
        lat = 90.0 - 0.5*r

        # Calculate lat, lon based on pixel index:
        lon = lon + x*r
        lat = lat - y*r

        return (lon, lat)

    def get_stationid(self):
        """
        Name:     MODISDATA.get_stationid
        Input:    None.
        Output:   int, station id (st_id)
        Features: Returns the half-degree (HDG) station ID for a pixel
                  numbered from 0 (bottom-left / south-west corner) to 259199
                  (top-right / north-east corner) as defined in the GePiSaT
                  database numbering scheme
        """
        st_id = (
            720.0 * (359.0 - ((90.0 - self.hdg_lat)/self.hdg_res - 0.5)) +
            ((self.hdg_lon + 180.0)/self.hdg_res - 0.5))
        return int(st_id)

    def get_foh_points(self, lon, lat, d):
        """
        Name:     MODISDATA.get_foh_points
        Input:    - float, 0.5 degree pixel longitude (lon)
                  - float, 0.5 degree pixel latitude (lat)
                  - numpy nd.array, 0.05 degree data (d)
        Output:   numpy nd.array (my_grid_data)
        Features: Returns array of 100 data values at 0.05 degrees based on a
                  single 0.5 degree pixel location
        Depends:  - get_foh_grid
                  - grid_to_index
        """
        my_grid_data = numpy.array([])
        my_grid_pnts = self.get_foh_grid()
        my_grid_indx = self.grid_to_index(my_grid_pnts)
        for indx_pair in my_grid_indx:
            x, y = indx_pair
            zval = d[y][x]
            # Use NaN for invalid data (easy to remove for averaging):
            if zval < -2000:
                zval = numpy.NaN
            my_grid_data = numpy.append(my_grid_data, [zval])

        return my_grid_data

    def get_foh_grid(self):
        """
        Name:     MODISDATA.get_foh_grid
        Input:    None.
        Output:   list, list of tuples (foh_grid)
        Features: Returns a list of one hundred 0.05 degree lon-lat pairs based
                  on the class's half degree (HDG) coordinates (i.e., hdg_lon,
                  hdg_lat)
        """
        # Initialize five one-hundreths grid:
        foh_grid = []

        # Calculate the binding box at 0.5 deg:
        westing = self.hdg_lon - 0.5*self.hdg_res
        northing = self.hdg_lat + 0.5*self.hdg_res

        # Initialize centroid offsetting for foh_grid:
        foh_lon = westing + 0.5*self.foh_res
        foh_lat = northing - 0.5*self.foh_res

        # Iterate over the 10x10 box:
        for y in range(10):
            lat = foh_lat - y*self.foh_res
            for x in range(10):
                lon = foh_lon + x*self.foh_res
                foh_grid.append((lon, lat))

        return foh_grid

    def grid_to_index(self, grid):
        """
        Name:     MODISDATA.grid_to_index
        Input:    list, list of lon-lat tuples (grid)
        Output:   list, list of tuples (foh_indices)
        Features: Returns a list of x-y indices based for a given list of 0.05
                  degree lon-lat pairs
        Depends:  get_x_y
        """
        foh_indices = []
        for grid_pair in grid:
            lon, lat = grid_pair
            x, y = self.get_x_y(lon, lat, self.foh_res)
            foh_indices.append((x, y))

        return foh_indices

    def get_x_y(self, lon, lat, r):
        """
        Name:     MODISDATA.get_x_y
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

    def print_line(self, outfile):
        """
        Name:     MODISDATA.print_line
        Input:    string, output file name and path (outfile)
        Output:   None.
        Features: Writes to file the contents for data_set database table
        """
        try:
            # Append to existing file:
            OUT = open(outfile, 'a')
            #
            # Create/write output lines:
            outline = "%s,%s,%s,%0.5f\n" % (
                self.msv_idx,
                self.station_id,
                self.data_time,
                self.data_val
                )
            OUT.write(outline)
        except IOError:
            self.logger.error("could not append to file: %s", outfile)
        else:
            OUT.close()


###############################################################################
# FUNCTIONS
###############################################################################
def get_modis_ts(f):
    """
    Name:     get_modis_ts
    Input:    string, file name (f)
    Output:   datetime.date
    Features: Returns a timestamp based on the MODIS filenaming scheme
              'MYD13C2.AYEARDOY...'
    """
    try:
        f_name = os.path.basename(f)
        f_date = re.search('\.A(\d{7})\.', f_name).group(1)
    except AttributeError:
        logging.error("Search failed for date in file %s", f_name)
        f_timestamp = None
    else:
        f_year = int(f_date[0:4])
        f_doy = int(f_date[4:7]) - 1
        f_timestamp = (
            datetime.date(f_year, 1, 1) + datetime.timedelta(days=f_doy))
    finally:
        return f_timestamp


def process_modis(my_dir, voi):
    """
    Name:     process_modis
    Input:    - string, input/output file directory (d)
              - string, variable of interest (voi)
    Output:   None.
    Features: Processes 0.05 degree MODIS EVI from HDF files into 0.5 degree
              FAPAR and saves to variable list and data set table output files
    Depends:  - find_files
              - get_modis_ts
              - writeout
    Ref:      pyhdf: http://pysclint.sourceforge.net/pyhdf/
              python-hdf4: http://fhs.github.io/python-hdf4/
              * provides pyhdf package
              * requires libhdf4 and libhdf4-dev (not to be mistaken with
                libhdf4-alt used by QGIS)

    TODO: follow watchdata.py example to update this method
    """
    # Search directory for MODIS HDF files:
    my_files = find_files(my_dir, "*.hdf")
    num_files = len(my_files)

    # Prepare var output file:
    var_file = "MODIS_Var-List_fapar.csv"
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Flag for varlist:
    varlist_flag = True

    # Get list of flux station 0.5-degree grid points:
    station_list = get_station_latlon()

    if num_files > 0:
        # Define and/or create the output directory (subdir of my_dir)
        out_dir = os.path.join(my_dir, "out")
        if not os.path.exists(out_dir):
            try:
                os.makedirs(out_dir)
            except:
                logging.warning(
                    "failed to create output directory; using input directory")
                out_dir = my_dir
            else:
                logging.info("created output directory %s", out_dir)

        var_path = os.path.join(out_dir, var_file)
        writeout(var_path, var_headerline)

        # Read through each file:
        for i in range(num_files):
            doc = my_files[i]
            logging.info("Processing file (%d/%d)", i+1, num_files)
            try:
                # Try to open file for reading:
                my_hdf = SD.SD(doc)
            except:
                logging.error("Unexpected error opening file %s", doc)
                data = None
            else:
                # Check that VOI is in the dataset:
                if voi in my_hdf.datasets().keys():
                    # Pull data from dataset into array:
                    d_select = my_hdf.select(voi)
                    data = d_select.get()
                else:
                    logging.error("Could not open dataset %s", voi)
                    data = None

                # Close HDF file and return data:
                my_hdf.end()

            # Continue processing if data was retrieved:
            if data is not None and data.any():
                # Get timestamp from filename:
                my_ts = get_modis_ts(doc)

                if my_ts is not None:
                    # Prepare data set output file:
                    dat_file = "MODIS_Data-Set_%s.csv" % my_ts
                    dat_path = os.path.join(out_dir, dat_file)
                    writeout(dat_path, dat_headerline)

                    # Iterate through each 0.5 deg lat:lon pair
                    # x (longitude): 0...719
                    # y (latitude): 0...359
                    for y in range(360):
                        for x in range(720):
                            my_modis = MODISDATA(x, y, my_ts, data)
                            pxl_lat = my_modis.hdg_lat
                            pxl_lon = my_modis.hdg_lon

                            # Filter grids based on flux stations:
                            if (pxl_lat, pxl_lon) in station_list:
                                st_parts = ('HDG', my_modis.hdg_id)

                                if varlist_flag:
                                    # ~~~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~ #
                                    my_line = VAR(st_parts, 'FAPAR', 'grid')
                                    my_line.printLine(var_path)
                                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                                my_modis.print_line(dat_path)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # Close varflag after processing first doc:
            varlist_flag = 0
