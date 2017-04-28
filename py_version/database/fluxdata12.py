#!/usr/bin/python
#
# fluxdata12.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-04-21
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
class FLUXDATA_2012:
    """
    Name:     FLUXDATA_2012
    Features: This class creates GePiSaT database data_set table entries based
              on fluxdata.org 2012 flux tower data files
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    msvIDX = []        # list of identifiers made up of station and var id
    stationID = ""     # station ID (from file name)
    dateTime = datetime.datetime(1999, 1, 1, 0, 0, 0)  # timestamp
    timeVals = []      # timestamp values
    vari = []          # list of variables (from headerline)
    data = []          # list of data (for each line of measurements)
    #
    # Variable quality control (qc) flags:
    varFlags = {'NEE_f': 'NEE_GPP_qc',
                'GPP_f': 'NEE_GPP_qc',
                'LE_f': 'LE_fqc',
                'H_f': 'H_fqc',
                'G_f': 'G_fqc',
                'Ta_f': 'Ta_fqc',
                'Ts1_f': 'Ts1_fqc',
                'Ts2_f': 'Ts2_fqc',
                'VPD_f': 'VPD_fqc',
                'Precip_f': 'Precip_fqc',
                'SWC1_f': 'SWC1_fqc',
                'SWC2_f': 'SWC2_fqc',
                'WS_f': 'WS_fqc',
                'Rg_f': 'Rg_fqc',
                'PPFD_f': 'PPFD_fqc',
                'Rn_f': 'Rn_fqc',
                'gsurf_f': 'gsurf_flag'}

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, fileName, v, d):
        """
        Name:     FLUXDATA_2012.__init__
        Input:    - string, file name with path (fileName)
                  - list, variable names, e.g., parsed header row (v)
                  - list, observation data, e.g., parsed data row (d)
        Depends:  VAR class
        """
        # Create a class logger
        self.logger = logging.getLogger("FLUXDATA_2012")
        self.logger.debug("FLUXDATA_2012 class initialized")

        # Calculate the number of fields:
        numfields = len(d)

        # Strip time values and calculate the timestamp:
        # WARNING: the midnight hour time stamp of the last day of the
        # year has the new year category and DoY = 366 or 367 which
        # will cause problems in getTS.

        # Flux data file time fields:
        # [0]: year
        # [1]: day of year
        # [2]: time (minutes)
        # [3]: datetime (fractional days)
        self.timeVals = d[0:3]
        try:
            # Extract the year with a regular expression from file name:
            year = re.search('\.(\d{4})\.', fileName).group(1)
        except AttributeError:
            # Search failed
            self.timeVals = d[0:3]
        else:
            # Add correct year to timevals:
            self.timeVals[0] = float(year)

        # Get timestamp:
        self.dateTime = self.getTS(self.timeVals)

        # Empty the class lists:
        self.msvIDX = []
        self.data = []
        self.vari = []

        # Get data (note 4:numfields skips the four time fields in CSV):
        for i in range(4, numfields):
            # If variable is core (key in varFlags dictionary):
            if (v[i] in self.varFlags.keys()):
                # Get qc value:
                qcflag = self.varFlags[v[i]]     # pull flag name from Dict
                qcindex = v.index(qcflag)   # find index corresp. to qc flag
                qcvalue = d[qcindex]        # get value of qc flag

                # Check that qc is for observation:
                if qcflag is "NEE_GPP_qc":
                    qcheck = self.check_gpp(qcvalue)
                else:
                    qcheck = self.check_vars(qcvalue)

                # Append only if core variable is an observation:
                if (qcheck == 1):
                    # Check variable against VAR class list of core vars
                    # and save only if variable is core.  Also, use this
                    # opportunity to grab the msvidx and stationid from VAR.
                    l = VAR(fileName, v[i], 'flux')
                    if l.varCore:
                        self.msvIDX.append(l.msvIDX)
                        self.vari.append(v[i])
                        self.data.append(d[i])
                        self.stationID = l.stationID

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def getTS(self, tvals):
        """
        Name:     FLUXDATA_2012.getTS
        Input:    list, fluxdata.org time field data (tvals)
                  [0] year, [1] day of year, [2] time (min), [3] datetime
        Output:   datetime.datetime
        Features: Returns the timestamp based on the three fluxdata.org time
                  fields: year, day of year, and time (minutes)
        """
        # Note: use float() because decimal point is in data file
        # Note: tvals[3] is not used (not necessary for building timestamp)
        yr = float(tvals[0])
        doy = float(tvals[1])-1
        mn = float(tvals[2])*60
        ts = (datetime.datetime(int(yr), 1, 1, 0, 0, 0) +
              datetime.timedelta(days=int(doy)) +
              datetime.timedelta(minutes=mn))
        return ts

    def check_gpp(self, qcf):
        """
        Name:     FLUXDATA_2012.check_gpp
        Input:    int, quality control flag (qcf)
        Output:   int, data type (rval)
        Features: Returns a value based on the data type defined by the
                  quality control flag for the NEE and GPP flux tower variables
                  0: missing, 1: observed, 2: gap-filled
        """
        # Initialize return value:
        rval = -1

        if (float(qcf) == -9999):
            # CHECK 1: if value is missing:
            rval = 0
        elif (float(qcf) == 1 or float(qcf) == 2):
            # CHECK 2: if value is original (obsevation):
            rval = 1
        elif (float(qcf) >= 3 and float(qcf) <= 6):
            # CHECK 3: if value has been gap-filled:
            rval = 2

        return rval

    def check_vars(self, qcf):
        """
        Name:     FLUXDATA_2012.check_vars
        Input:    int, quality control flag (qcf)
        Output:   int, data type (rval)
        Features: Returns a value based on the data type defined by the
                  quality control flag for state variables (not NEE or GPP)
                  0: missing, 1: observed, 2: gap-filled
        """
        # Initialize return value
        rval = -1
        if (float(qcf) == -9999):
            # CHECK 1: if value is missing:
            rval = 0
        elif (float(qcf) == 0):
            # CHECK 2: if value is original (observation):
            rval = 1
        elif (float(qcf) > 0 and float(qcf) <= 3):
            # CHECK 3: if value has been gap-filled:
            rval = 2

        return rval

    def print_line(self, outfile):
        """
        Name:     FLUXDATA_2012.print_line
        Input:    string, output file and path (outfile)
        Output:   None.
        Features: Writes to file the contents for the data_set database table
        """
        try:
            # Append data to output file:
            OUT = open(outfile, 'a')

            # Create/write output lines:
            for i in range(len(self.data)):
                outline = "%s,%s,%s,%0.3f\n" % (
                    self.msvIDX[i],
                    self.stationID,
                    self.dateTime,
                    float(self.data[i])
                    )
                OUT.write(outline)
        except:
            self.logger.error("could not append to file: %s", outfile)
        else:
            OUT.close()


###############################################################################
# FUNCTIONS
###############################################################################
def process_flux_2012(my_dir):
    """
    Name:     process_flux_2012
    Input:    string, input file directory (my_dir)
    Output:   None.
    Features: Processes 2012 flux tower data files into variable list and
              data set table output files
    Depends:  - find_files
              - writeout
    """
    # Search directory for fluxdata CSV files:
    my_ext = "*allvars.csv"
    my_files = find_files(my_dir, my_ext)

    # Prepare header lines and write to file:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = True

    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = True

    # Check that files were found:
    if my_files:
        # Read through each file:
        for doc in my_files:
            if os.path.isfile(doc):
                doc_name = os.path.basename(doc)
                try:
                    # Try to get filename prefix:
                    my_station = re.search('^\S{6}', doc_name).group(0)
                    my_dat = re.search('^\S{6}\.{1}\d{4}', doc_name).group(0)
                    my_dat_out = my_dat + "_Data-Set.csv"
                    my_var_out = my_station + "_Var-List.csv"
                except AttributeError:
                    logging.error("Could not read file prefix from %s", doc)
                else:
                    # Initialize data-set outfile:
                    if dat_flag:
                        dat_outfile = os.path.join(my_dir, my_dat_out)
                        writeout(dat_outfile, dat_headerline)

                    # Check to see if station has been var'ed:
                    if my_station in var_out_list:
                        var_flag = False
                    else:
                        # Initialize var outfile:
                        var_outfile = os.path.join(my_dir, my_var_out)
                        writeout(var_outfile, var_headerline)

                        # Add var to var_list and set flag to true:
                        var_out_list.append(my_station)
                        var_flag = True

                try:
                    # Try to open file for reading:
                    f = open(doc, 'r')

                    # Read the headerline and strip whitespace from the end:
                    header = f.readline()
                    header = header.rstrip("\n")

                    # Read the remaining content:
                    content = f.readlines()
                except IOError:
                    logging.error("Could not read file: %s", doc)
                else:
                    # Close the file after done reading:
                    f.close()

                    # Parse header line:
                    parts = header.split(',')

                    # ~~~~~~~~~~ VAR-LIST ~~~~~~~~~~ #
                    # For each flux variable, produce the data row:
                    if var_flag:
                        for var in parts:
                            my_line = VAR(doc, var, 'flux')
                            my_line.printLine(var_outfile)

                    # ~~~~~~~~~~ DATA-SET ~~~~~~~~~~ #
                    # For each timestamp, produce variable data
                    if dat_flag:
                        for row in content:
                            row = row.rstrip("\n")
                            dataset = row.split(',')

                            # Create data class and print out lines:
                            my_data = FLUXDATA_2012(doc, parts, dataset)
                            my_data.print_line(dat_outfile)
    else:
        logging.warning("No files found in directory: %s", my_dir)
