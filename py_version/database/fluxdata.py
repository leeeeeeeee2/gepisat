#!/usr/bin/python
#
# fluxdata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-01-15
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

try:
    from .utilities import writeout
    from .var import VAR
except (ImportError, ValueError):
    pass


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
        self.logger.info("FLUXDATA_2012 class initialized")

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


class FLUXDATA_2015:
    """
    Name:     FLUXDATA_2015
    Features: This class creates GePiSaT database data_set table entries based
              on fluxdata.org 2015 flux tower data files
    Ref:
    http://fluxnet.fluxdata.org/data/fluxnet2015-dataset/fullset-data-product/
    NOTE:
        FLUXDATA FULLSET data file naming scheme:
        > HH - half-hourly data
        > HR - hourly data
        CSV Header Items of Interest:
        > TIMESTAMP_START
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
    #
    # Variable quality control (qc) flags:
    varFlags = {
        'NEE_VUT_REF': 'NEE_VUT_REF_QC',   # 0 = measured
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
        self.logger.info("FLUXDATA_2015 class initialized")

        # Calculate the number of fields:
        numfields = len(d)
        self.logger.debug("Read %d fields", numfields)

        # Empty the class lists:
        self.msvIDX = []
        self.data = []
        self.vari = []

        # @TODO: convert timekeeping string (YYYYMMDDHHMM) to datetime object
        # @TODO: find list indexes of variables of interest

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
def find_files(my_dir, my_pattern):
    """
    Name:     find_files
    Inputs:   - str, directory path (my_dir)
              - str, file name search pattern (my_pattern)
    Outputs:  list, file paths
    Features: Returns a list of files found at a given directory with file
              names that match a given pattern
    """
    my_files = []
    if os.path.isdir(my_dir):
        s_str = os.path.join(my_dir, my_pattern)
        my_files = glob.glob(s_str)
    return my_files


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


def process_flux_2012(my_dir):
    """
    Name:     process_flux_2012
    Input:    string, input file directory (my_dir)
    Output:   None.
    Features: Processes 2012 flux tower data files into variable list and
              data set table output files
    Depends:  writeout
    """
    # Search directory for fluxdata CSV files:
    my_ext = "*allvars.csv"
    s_str = os.path.join(my_dir, my_ext)
    my_files = glob.glob(s_str)

    # Prepare header lines and write to file:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = 1

    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = 1

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
                        var_flag = 0
                    else:
                        # Initialize var outfile:
                        var_outfile = os.path.join(my_dir, my_var_out)
                        writeout(var_outfile, var_headerline)

                        # Add var to var_list and set flag to true:
                        var_out_list.append(my_station)
                        var_flag = 1

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


def process_flux_2015(my_dir):
    """
    Name:     process_flux_2015
    Input:    string, input file directory (my_dir)
    Output:   None.
    Features: Processes 2015 flux tower data files into variable list and
              data set table output files
    Depends:  writeout
    """
    # Search directory for fluxdata CSV files:
    my_ext = "*_FLUXNET2015_FULLSET_HH_*"
    s_str = os.path.join(my_dir, my_ext)
    my_files = glob.glob(s_str)

    # Prepare header lines for output files:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = 1

    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = 1

    # Check that files were found:
    n_files = len(my_files)
    if n_files > 0:
        # Read through each file:
        for doc in my_files:
            if os.path.isfile(doc):
                doc_name = os.path.basename(doc)
                try:
                    # Try to get filename prefix:
                    my_station = re.search('^FLX_(\S{6})_', doc_name).group(1)
                    my_dat = re.search('_(\d{4}-\d{4})_', doc_name).group(1)
                    my_dat_out = "%s_%s_Data-Set.csv" % (my_station, my_dat)
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
                        var_flag = 0
                    else:
                        # Initialize var outfile:
                        var_outfile = os.path.join(my_dir, my_var_out)
                        writeout(var_outfile, var_headerline)

                        # Add var to var_list and set flag to true:
                        var_out_list.append(my_station)
                        var_flag = 1

                    try:
                        # Try to open file for reading:
                        f = open(doc, 'r')

                        # Read headerline and strip whitespace from the end:
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
    else:
        logging.warning("No files found in directory: %s", my_dir)


###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    import numpy
    import matplotlib.pyplot as plt

    my_dir = "/usr/local/share/data/fluxnet/2015/half_hourly"
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
                i_ppfd_in = h_items.index("PPFD_IN")
            except ValueError as err:
                print("Skipping %s. %s" % (my_file, str(err)))
            else:
                # Read data into structured array:
                my_data = numpy.loadtxt(
                    fname=my_path,
                    dtype={'names': ('sw_in_f', 'sw_in_f_qc', 'ppfd_in'),
                           'formats': ('f4', 'f4', 'f4')},
                    delimiter=',',
                    skiprows=1,
                    usecols=(i_sw_in_f, i_sw_in_f_qc, i_ppfd_in)
                )

                # Extract measurements:
                i_measured = numpy.where(
                    (my_data['sw_in_f_qc'] == 0) &
                    (my_data['sw_in_f'] != -9999) &
                    (my_data['ppfd_in'] != -9999))[0]
                m_sw_in_f = my_data['sw_in_f'][i_measured]
                m_ppfd_in = my_data['ppfd_in'][i_measured]

                # Get correlation coefficient:
                my_r = pearsons_r(m_ppfd_in, m_sw_in_f)

                # Let's see what the data look like:
                if True:
                    fig = plt.figure()
                    ax1 = fig.add_subplot(111)
                    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=12)
                    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=12)
                    ax1.plot(m_ppfd_in, m_sw_in_f, 'ro', label=label_str)
                    ax1.set_ylabel('SW_IN_F, W m$^{-2}$', fontsize=14)
                    ax1.set_xlabel('PPFD_IN, W m$^{-2}$', fontsize=14)
                    ax1.legend(
                        bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                        ncol=2, mode="expand", borderaxespad=0., fontsize=14)
                    plt.show()

                # Linear regression:
                xdata = m_ppfd_in.copy()
                x = xdata[:, numpy.newaxis]
                y = m_sw_in_f.copy()
                fit_slope, fit_sse, fit_rank, fit_s = numpy.linalg.lstsq(x, y)

                # Calculate variance/st. dev of slope
                s2xy = fit_sse/(len(x) - 2.0)
                x_bar = xdata.mean()
                ssxx = (numpy.power(xdata - x_bar, 2.0)).sum()
                s2m = s2xy/ssxx
                sm = numpy.sqrt(s2m)

                print("%s: %s +/- %0.3f (r = %0.3f)" % (
                    label_str, fit_slope[0], sm, my_r))

                # Clear variables:
                h_items = None
                i_sw_in_f = None
                i_sw_in_f_qc = None
                i_ppfd_in = None
                i_measured = None
                my_data = None
                m_sw_in_f = None
                m_ppfd_in = None
                x_data = None
                x = None
                y = None

                #
                # RESULTS:
                # AR-Vir: SW_IN = 0.52 PPFD_IN
