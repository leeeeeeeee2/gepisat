#!/usr/bin/python
#
# splashdata.py
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
import os
import re

import numpy

from .utilities import writeout
from .var import VAR


###############################################################################
# FUNCTIONS
###############################################################################
def process_alpha(d):
    """
    Name:     process_alpha
    Input:    string, input/output file directory (d)
    Features: Processes the Cramer-Prentice alpha raster files into variable
              list and data set table output files
    Depends:  writeout
    """
    # Read the ASCII raster files from directory:
    my_dir = d
    my_files = glob.glob(my_dir + "CP-alpha_*txt")

    # Prepare var output file:
    my_var_out = "Var-List_alpha.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_outfile, var_headerline)

    # Flag for varlist:
    varlist_flag = 1

    # Read through files:
    if my_files:
        for doc in my_files:
            try:
                # Load file, skipping header lines:
                f = numpy.loadtxt(doc, skiprows=6)

                # Read timestamp from filename:
                yr_mo = re.search(
                    '_(\d{4}-\d{2}-\d{2})\.',
                    os.path.basename(doc)
                    ).group(1)
                this_year = int(yr_mo.split('-')[0])
                this_month = int(yr_mo.split('-')[1])
                my_ts = datetime.date(this_year, this_month, 1)
            except:
                logging.error("Error reading file, %s", doc)
            else:
                # Save the data shape (360x720):
                (sh_lat, sh_lon) = f.shape

                # Prepare data set output file:
                my_dat_out = "Data-Set_alpha_%s.csv" % my_ts
                dat_outfile = my_dir + my_dat_out
                dat_headerline = "msvidx,stationid,datetime,data\n"
                writeout(dat_outfile, dat_headerline)

                # Iterate through data file:
                # NOTE: ASCII raster is organized in 360 rows and 720 cols
                # * rows ordered from 89.75 to -89.75 lat (north > south)
                # * cols ordered from -179.75 to 179.75 lon (east > west)
                for y in range(sh_lat):
                    # Reverse order latitude index for station numbering:
                    i = 359 - y
                    for x in range(sh_lon):
                        # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                        st_id = 720*i + x
                        station_parts = ('HDG', st_id)
                        my_line = VAR(station_parts, 'alpha', 'grid')
                        if varlist_flag:
                            my_line.printLine(var_outfile)
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        #
                        # Read pixel value (1000x alpha):
                        # NOTE: missing values are -9999
                        pxl = f[y, x]
                        if pxl != -9999:
                            # Scale pixel to alpha:
                            pxl = 1.0 * pxl / 1000.0

                            # ~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~ #
                            # Append to existing file:
                            OUT = open(dat_outfile, 'a')

                            # Create/write output line:
                            outline = "%s,%s,%s,%0.3f\n" % (
                                my_line.msvIDX,
                                my_line.stationID,
                                my_ts,
                                pxl
                                )
                            OUT.write(outline)
                            OUT.close()
                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                # Turn off varlist flag after processing first month:
                varlist_flag = 0
