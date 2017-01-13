#!/usr/bin/python
#
# watchdata.py
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
# (GePiSaT) Model of the terrestrial biosphere: Part 1 - Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.

###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import glob
import logging
import os
import re

from scipy.io import netcdf

from .utilities import writeout
from .var import VAR


###############################################################################
# CLASSES
###############################################################################
class WATCHDATA:
    """
    Name:     WATCHDATA
    Features: This class processes WATCH observations for the data_set table
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    msv_idx = []        # specific identifier made up of station and var id
    station_id = ""     # station ID (from file name)
    data_time = datetime.date(1900, 1, 1)  # timestamp
    data_value = -9999.0      # observation data

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, st_parts, val_parts):
        """
        Name:     WATCHDATA.__init__
        Input:    - tuple, station prefix and ID (st_parts)
                  - tuple, variable name, value, and timestamp (val_parts)
        """
        # Create a class logger
        self.logger = logging.getLogger("WATCHDATA")
        self.logger.info("WATCHDATA class initialized")

        # Create a variable class & save stationid and msvidx values:
        my_var = VAR(st_parts, val_parts[0], 'grid')
        self.station_id = my_var.stationID
        self.msv_idx = my_var.msvIDX

        # Save observation and time stamp:
        self.data_value = val_parts[1]
        self.data_time = val_parts[2]

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def print_line(self, outfile):
        """
        Name:     WATCHDATA.print_line
        Input:    string, output file name and path (outfile)
        Output:   None.
        Features: Writes to file the contents for data_set database table
        """
        try:
            # Append to existing file:
            OUT = open(outfile, 'a')

            # Create/write output lines:
            outline = "%s,%s,%s,%0.5f\n" % (
                self.msv_idx,
                self.station_id,
                self.data_time,
                self.data_value
                )
            OUT.write(outline)
        except:
            self.logger.error("could not append to file: %s", outfile)
        else:
            OUT.close()


###############################################################################
# FUNCTIONS
###############################################################################
def process_watch(my_dir, voi):
    """
    Name:     process_watch
    Input:    - string, input file directory (my_dir)
              - string, variable of interest (voi)
    Output:   None.
    Features: Processes WATCH WFDEI netCDF files into variable list and data
              set table output files
    Depends:  writeout
    """
    # Search directory for WATCH netCDF files:
    s_str = os.path.join(my_dir, "*.nc")
    my_files = glob.glob(s_str)

    # Prepare var output file:
    var_file = "WFDEI_Var-List_test.csv"
    var_path = os.path.join(my_dir, var_file)
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_path, var_headerline)

    # Flag for varlist:
    varlist_flag = True

    if my_files:
        # Read through each file:
        for doc in my_files:
            if os.path.isfile(doc):
                try:
                    # Try to open file for reading:
                    f = netcdf.NetCDFFile(doc, "r")

                    # Read the year and month values from filename (YYYYMM):
                    yr_mo = re.search(
                        '_(\d{6})\.', os.path.basename(f.filename)).group(1)
                    this_year = int(yr_mo[0:4])
                    this_month = int(yr_mo[4:6])
                except IOError:
                    logging.error("Could not read file: %s", doc)
                except AttributeError:
                    logging.error("Year and month not found from file %s", doc)
                except ValueError:
                    logging.error("Year and month not numbers in file %s", doc)
                else:
                    # Prepare data set output file:
                    dat_file = "WFDEI_Data-Set_%s.csv" % yr_mo
                    dat_path = os.path.join(my_dir, dat_file)
                    dat_headerline = "msvidx,stationid,datetime,data\n"
                    writeout(dat_path, dat_headerline)

                    # Save the shape values of each variable:
                    sh_day, sh_lat, sh_lon = f.variables[voi].shape

                    # Iterate through each lat:lon pair
                    # * row-major ordering from bottom left
                    #  x (longitude): 0...719
                    #  y (latitude): 0...359
                    for y in range(sh_lat):
                        # Save latitude:
                        pxl_lat = f.variables['lat'].data[y]

                        for x in range(sh_lon):
                            # Calc station ID:
                            st_id = 720*y + x
                            st_parts = ('HDG', st_id)

                            if varlist_flag:
                                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                                my_line = VAR(st_parts, voi, 'grid')
                                my_line.printLine(var_path)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                            # Iterate through each day
                            for t in range(sh_day):
                                # Get timestamp for this day:
                                this_day = t+1
                                time_stamp = datetime.date(
                                    this_year, this_month, this_day)

                                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                                # Read SWdown for each pixel
                                # * NOTE 1: variable has five decimal places
                                # * NOTE 2: missing values are equal to ~1e20
                                # * NOTE 3: Antarctica is < -60 latitude
                                pxl_swr = f.variables[voi].data[t, y, x]
                                if pxl_swr < 1.0e6 and pxl_lat > -60:
                                    # Process pixel
                                    obs_parts = (voi, pxl_swr, time_stamp)
                                    my_data = WATCHDATA(st_parts, obs_parts)
                                    my_data.print_line(dat_path)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    # Close the var list flag after processing first doc:
                    varlist_flag = False
                    pxl_lat = None
                    pxl_swr = None
                    f.close()
