#!/usr/bin/python
#
# splashdata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-04-28
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

from .utilities import find_files
from .utilities import writeout
from .var import VAR


###############################################################################
# CLASSES
###############################################################################
class SPLASH_DATA:
    """
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    msvIDX = None        # list of identifiers made up of station and var id
    stationID = None     # station ID (from file name)
    dateTime = datetime.datetime(1999, 1, 1, 0, 0, 0)  # timestamp
    vari = None          # list of variables (from headerline)
    data = None          # list of data (for each line of measurements)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, fileName, v, d):
        """
        Name:     SPLASH_DATA.__init__
        Input:    - string, file name with path (fileName)
                  - list, variable names, e.g., parsed header row (v)
                  - list, observation data, e.g., parsed data row (d)
        Depends:  VAR class
        """
        # Create a class logger
        self.logger = logging.getLogger("FLUXDATA_2015")
        self.logger.debug("FLUXDATA_2015 class initialized")

        # Calculate the number of fields:
        numfields = len(d)
        self.logger.debug("Read %d fields", numfields)

        # Process timestamp:
        try:
            i_time = v.index("Month")
        except:
            self.logger.error("No timestamp!")
        else:
            m_time_str = d[i_time]
            try:
                m_time = datetime.datetime.strptime(m_time_str, "%Y-%m-%d")
            except:
                self.logger.error("Timestamp format error!")
            else:
                self.dateTime = m_time

                # Process data:
                try:
                    i_alpha = v.index("alpha")
                except ValueError as err:
                    self.logger.error("Missing index! %s" % (str(err)))
                else:
                    m_alpha = float(d[i_alpha])
                    var_alpha = VAR(fileName, "alpha", "flux")

                    self.msvIDX = var_alpha.msvIDX
                    self.vari = "alpha"
                    self.data = m_alpha
                    self.stationID = var_alpha.stationID

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def print_line(self, outfile):
        """
        Name:     SPLASH_DATA.print_line
        Input:    string, output file and path (outfile)
        Output:   None.
        Features: Writes to file the contents for the data_set database table
        """
        try:
            # Append data to output file:
            OUT = open(outfile, 'a')

            # Create/write output lines:
            outline = "%s,%s,%s,%0.3f\n" % (
                self.msvIDX,
                self.stationID,
                self.dateTime,
                float(self.data)
            )
            OUT.write(outline)
        except:
            self.logger.error("could not append to file: %s", outfile)
        else:
            OUT.close()


###############################################################################
# FUNCTIONS
###############################################################################
def process_alpha(my_dir):
    """
    Name:     process_alpha
    Input:    string, input/output file directory (my_dir)
    Features: Processes the Cramer-Prentice alpha raster files into variable
              list and data set table output files
    Depends:  writeout
    """
    # Read the ASCII raster files from directory:
    my_files = find_files(my_dir, "*bioindex.csv")

    # Prepare headerlines for var and data files:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = True

    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = True

    # Check that files were found:
    n_files = len(my_files)
    if n_files > 0:
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

        # Read through each file:
        i = 0
        for doc in my_files:
            i += 1
            if os.path.isfile(doc):
                doc_name = os.path.basename(doc)

                try:
                    # Try to get filename prefix:
                    # NOTE: VAR class does this search for you ...
                    my_station = re.search('^(\S{6})_', doc_name).group(1)
                    my_dat_out = "%s_alpha_Data-Set.csv" % (my_station)
                    my_var_out = "%s_alpha_Var-List.csv" % (my_station)
                except AttributeError:
                    logging.error("Could not read file prefix from %s", doc)
                else:
                    logging.info(
                        "Processing (%d/%d): %s" % (i, n_files, my_station))

                    # Initialize data-set outfile:
                    if dat_flag:
                        dat_outfile = os.path.join(out_dir, my_dat_out)
                        writeout(dat_outfile, dat_headerline)

                    # Check to see if station has been var'ed:
                    if my_station in var_out_list:
                        var_flag = False
                    else:
                        # Initialize var outfile:
                        var_outfile = os.path.join(out_dir, my_var_out)
                        writeout(var_outfile, var_headerline)

                        # Add var to var_list and set flag to true:
                        var_out_list.append(my_station)
                        var_flag = True

                    try:
                        # Try to open file for reading:
                        f = open(doc, 'r')

                        # Read headerline and remaining content:
                        header = f.readline().rstrip("\n")
                        content = f.readlines()
                    except IOError:
                        logging.error("Could not read file: %s", doc)
                    else:
                        # Close the file after done reading:
                        f.close()

                        # Parse header line:
                        var_parts = header.split(',')

                        # ~~~~~~~~~~ VAR-LIST ~~~~~~~~~~ #
                        # For each flux variable, produce the data row:
                        if var_flag:
                            for var in var_parts:
                                my_line = VAR(doc, var, 'flux')
                                my_line.printLine(var_outfile)

                        # ~~~~~~~~~~ DATA-SET ~~~~~~~~~~ #
                        # For each timestamp, produce variable data
                        if dat_flag:
                            for row in content:
                                row = row.rstrip("\n")
                                dset = row.split(',')

                                # Create data class and print out lines:
                                my_data = SPLASH_DATA(doc, var_parts, dset)
                                my_data.print_line(dat_outfile)
