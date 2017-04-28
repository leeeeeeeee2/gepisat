#!/usr/bin/python
#
# fluxdata15.py
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
class FLUXDATA_2015:
    """
    Name:     FLUXDATA_2015
    Features: This class creates GePiSaT database data_set table entries based
              on fluxdata.org 2015 flux tower data files
    History:  Version 3.0
              - rough draft of processing procedure [17.01.22]
    Ref:
    http://fluxnet.fluxdata.org/data/fluxnet2015-dataset/fullset-data-product/
    Note:
        FLUXDATA FULLSET data file naming scheme:
        > HH - half-hourly data
        > HR - hourly data

        CSV Header Items of Interest:
        > TIMESTAMP_START (YYYYMMDDHHMM)
        > PPFD_IN (W m-2)
          * not available at all stations!
        > SW_IN_F (W m-2)
        > SW_IN_F_QC (0 == measured)
        > NEE_VUT_REF (umolCO2 m-2 s-1)
        > NEE_VUT_REF_QC (0 == measured)
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

    # Variable quality control (qc) flags:
    varFlags = {
        'NEE_VUT_REF': 'NEE_VUT_REF_QC',   # 0 = measured
        'SW_IN_F': 'SW_IN_F_QC',           # 0 = measured
        'PPFD_IN': 'PPFD_IN_QC'            # not defined for half-hourly data
    }

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, fileName, v, d):
        """
        Name:     FLUXDATA_2015.__init__
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

        # Empty the class lists:
        self.msvIDX = []
        self.data = []
        self.vari = []

        # Process timestamp:
        try:
            i_time = v.index("TIMESTAMP_START")
        except:
            self.logger.error("No timestamp!")
        else:
            m_time_str = d[i_time]
            try:
                m_time = datetime.datetime.strptime(m_time_str, "%Y%m%d%H%M")
            except:
                self.logger.error("Timestamp format error!")
            else:
                self.dateTime = m_time

                # Process data:
                try:
                    i_sw_in_f = v.index("SW_IN_F")
                    i_sw_in_f_qc = v.index("SW_IN_F_QC")
                    i_nee_ref = v.index("NEE_VUT_REF")
                    i_nee_ref_qc = v.index("NEE_VUT_REF_QC")
                except ValueError as err:
                    self.logger.error("Missing index! %s" % (str(err)))
                else:
                    m_sw_in_f = float(d[i_sw_in_f])
                    m_sw_in_f_qc = float(d[i_sw_in_f_qc])
                    m_nee_ref = float(d[i_nee_ref])
                    m_nee_ref_qc = float(d[i_nee_ref_qc])

                    # Check for measured quality flag:
                    if m_sw_in_f_qc == 0:
                        # Grab the msvidx and stationid from VAR.
                        var_sw = VAR(fileName, "SW_IN_F", 'flux')
                        self.msvIDX.append(var_sw.msvIDX)
                        self.vari.append("SW_IN_F")
                        self.data.append(m_sw_in_f)
                        self.stationID = var_sw.stationID

                    if m_nee_ref_qc == 0:
                        var_nee = VAR(fileName, "NEE_VUT_REF", 'flux')
                        self.msvIDX.append(var_nee.msvIDX)
                        self.vari.append("NEE_VUT_REF")
                        self.data.append(m_nee_ref)
                        self.stationID = var_nee.stationID

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def print_line(self, outfile):
        """
        Name:     FLUXDATA.print_line
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
def pearsons_r(x, y):
    """
    Name:     pearsons_r
    Inputs:   - numpy.ndarray, ordinate/observation (x)
              - numpy.ndarray, abscissa/modeled (y)
    Outputs:  float, Pearson's r
    Features: Returns Pearson's r, the correlation coefficient, between two
              data arrays
    """
    n = float(len(x))
    ssxy = (x*y).sum() - x.sum()*y.sum()/n
    ssx = numpy.power(y, 2.0).sum() - numpy.power(y.sum(), 2.0)/n
    sst = numpy.power(x, 2.0).sum() - numpy.power(x.sum(), 2.0)/n
    pr = ssxy/numpy.power(ssx*sst, 0.5)
    return pr


def process_flux_2015(my_dir):
    """
    Name:     process_flux_2015
    Input:    string, input file directory (my_dir)
    Output:   None.
    Features: Processes 2015 flux tower data files into variable list and
              data set table output files
    Depends:  - find_files
              - writeout
    """
    # Search directory for fluxdata CSV files:
    my_ext = "*_FLUXNET2015_FULLSET_HH_*"
    my_files = find_files(my_dir, my_ext)

    # Prepare header lines for output files:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = True

    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = True

    # NOTE: uncomment to filter by date
    # s_time = datetime.datetime(2002, 1, 1)    # used to check data starting
    # e_time = datetime.datetime(2007, 1, 1)    # and ending times

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
                    my_station = re.search('^FLX_(\S{6})_', doc_name).group(1)
                    my_years = re.search('_(\d{4}-\d{4})_', doc_name).group(1)
                    my_dat_out = "%s_%s_Data-Set.csv" % (my_station, my_years)
                    my_var_out = my_station + "_Var-List.csv"
                except AttributeError:
                    logging.error("Could not read file prefix from %s", doc)
                else:
                    logging.info(
                        "Processing (%d/%d): %s (%s)" % (
                            i, n_files, my_station, my_years))

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

                        # Read headerline and strip whitespace from the end:
                        header = f.readline().rstrip("\n")

                        # Read the remaining content:
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
                                my_data = FLUXDATA_2015(doc, var_parts, dset)
                                my_data.print_line(dat_outfile)
                                # NOTE: uncomment to filter by date
                                # if my_data.dateTime >= s_time:
                                #    if my_data.dateTime < e_time:
                                #        my_data.print_line(dat_outfile)
    else:
        logging.warning("No files found in directory: %s", my_dir)


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    import numpy
    import matplotlib.pyplot as plt

    # my_dir = "/usr/local/share/data/fluxnet/2015/half_hourly"
    my_dir = "/home/user/Desktop/temp/flux"
    my_ext = "*_FLUXNET2015_FULLSET_HH_*"
    my_file_paths = find_files(my_dir, my_ext)
    for my_path in sorted(my_file_paths):
        my_file = os.path.basename(my_path)
        my_station = re.search('^FLX_(\S{6})_', my_file).group(1)
        my_years = re.search('_(\d{4}-\d{4})_', my_file).group(1)
        label_str = "%s (%s)" % (my_station, my_years)

        if os.path.isfile(my_path):
            # Reader header line:
            f = open(my_path, 'r')
            header = f.readline().rstrip("\n")
            f.close()

            # Find indexes for items interested in
            h_items = header.split(",")
            try:
                i_sw_in_f = h_items.index("SW_IN_F")
                i_sw_in_f_qc = h_items.index("SW_IN_F_QC")
                i_nee_ref = h_items.index("NEE_VUT_REF")
                i_nee_ref_qc = h_items.index("NEE_VUT_REF_QC")
            except ValueError as err:
                print("Skipping %s. %s" % (my_file, str(err)))
            else:
                # Read data into structured array:
                my_data = numpy.loadtxt(
                    fname=my_path,
                    dtype={'names': ('sw_in_f', 'sw_in_f_qc',
                                     'nee_vut_ref', 'nee_vut_ref_qc'),
                           'formats': ('f4', 'f4', 'f4', 'f4')},
                    delimiter=',',
                    skiprows=1,
                    usecols=(i_sw_in_f, i_sw_in_f_qc, i_nee_ref, i_nee_ref_qc)
                )

                # Extract measurements:
                i_measured = numpy.where(
                    (my_data['sw_in_f_qc'] == 0) &
                    (my_data['sw_in_f'] != -9999) &
                    (my_data['nee_vut_ref_qc'] == 0) &
                    (my_data['nee_vut_ref'] != -9999))[0]
                m_sw_in_f = my_data['sw_in_f'][i_measured]
                m_nee_ref = my_data['nee_vut_ref'][i_measured]

                # Get correlation coefficient:
                my_r = pearsons_r(m_nee_ref, m_sw_in_f)

                # Let's see what the data look like:
                if False:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=12)
                    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
                    ax1.plot(m_sw_in_f, m_nee_ref, 'ro', label=label_str)
                    ax1.set_xlabel('SW_IN_F, W m$^{-2}$', fontsize=14)
                    ax1.set_ylabel(
                        'NEE_VUT_REF, mol m$^{-2}$ s$^{-1}$', fontsize=14)
                    ax1.legend(
                        bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                        ncol=2, mode="expand", borderaxespad=0., fontsize=14)
                    plt.show()

                print("%s: %d (r = %0.3f)" % (
                    label_str, len(i_measured), my_r))

                # Clear variables:
                h_items = None
                i_sw_in_f = None
                i_sw_in_f_qc = None
                i_nee_ref = None
                i_nee_ref_qc = None
                i_measured = None
                my_data = None
                m_sw_in_f = None
                m_nee_ref = None
