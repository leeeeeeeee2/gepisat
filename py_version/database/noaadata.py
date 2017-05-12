#!/usr/bin/python
#
# noaadata.py
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

import numpy

from .utilities import find_files
from .utilities import get_station_latlon
from .utilities import writeout
from .var import VAR


###############################################################################
# FUNCTIONS
###############################################################################
def get_noaa_file(my_dir):
    """
    Name:     get_noaa_file
    Inputs:   str, directory path (my_dir)
    Outputs:  str, file path
    Features: Returns the global annual mean CO2 file from a given directory
    Depends:  find_files
    """
    # Search directory for WATCH netCDF files:
    my_files = find_files(my_dir, "co2_annmean_gl.txt")
    num_files = len(my_files)
    if num_files == 1:
        return my_files[0]
    else:
        raise IOError("Failed to find data file!")


def process_co2(my_dir):
    """
    Name:     process_co2
    Input:    string, input file directory (my_dir)
    Output:   None.
    Features: Processes NOAA ESRL global mean annual CO2 data into variable
              list and data set table output files
    Depends:  - get_noaa_file
              - get_station_latlon
              - writeout
    """
    # Prepare var output file and output file header lines:
    var_file = "NOAA_Var-List_co2.csv"
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Flag for varlist:
    varlist_flag = True

    # Get list of flux station 0.5-degree grid points:
    station_list = get_station_latlon()

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

    try:
        # Find data file and open it
        my_file = get_noaa_file(my_dir)
    except:
        logging.error("Could not open NOAA file")
    else:
        my_data = numpy.fromregex(
            my_file,
            r"\s+(\d{4})\s+(\d{3}\.{1}\d{2})\s+(\d{1}\.{1}\d{2})",
            [('year', 'int'), ('mean', 'float'), ('unc', 'float')])

        if len(my_data) > 0:
            # Process each year
            for line in my_data:
                # Extract date and value from data:
                year, val, _ = line
                st_date = datetime.date(year, 1, 1)

                logging.info("Processing %s", st_date)

                # Prepare data set output file:
                dat_file = "NOAA_Data-Set_co2_%s.csv" % year
                dat_path = os.path.join(out_dir, dat_file)
                writeout(dat_path, dat_headerline)

                # Iterate through each lat:lon pair
                #  x (longitude): 0...719
                #  y (latitude): 0...359
                sh_lat = 360
                sh_lon = 720
                latitude = [-89.75 + i*0.5 for i in range(sh_lat)]
                longitude = [-179.75 + i*0.5 for i in range(sh_lon)]
                for y in range(sh_lat):
                    pxl_lat = latitude[y]
                    for x in range(sh_lon):
                        pxl_lon = longitude[x]

                        # Filter grids based on flux stations:
                        if (pxl_lat, pxl_lon) in station_list:
                            # Calc station ID:
                            st_id = 720*y + x
                            st_parts = ('HDG', st_id)
                            st_voi = "CO2"

                            my_var = VAR(st_parts, st_voi, 'grid')
                            if varlist_flag:
                                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                                my_var.printLine(var_path)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                            # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                            OUT = open(dat_path, 'a')
                            outline = "%s,%s,%s,%0.5f\n" % (
                                my_var.msvIDX,
                                my_var.stationID,
                                st_date,
                                val)
                            OUT.write(outline)
                            OUT.close()
                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                # Close the var list flag after processing first doc:
                varlist_flag = False
