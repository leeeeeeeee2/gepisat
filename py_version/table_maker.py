#!/usr/bin/python
#
# table_maker.py
#
# VERSION 2.02
#
# 2013-05-14 -- created
# 2015-11-18 -- last updated
#
# ---------
# citation:
# ---------
# I. C. Prentice, T. W. Davis, X. M. P. Gilbert, B. D. Stocker, B. J. Evans,
# H. Wang, and T. F. Keenan, "The Global ecosystem in Space and Time (GePiSaT)
# Model of the Terrestrial Biosphere," (in progress).
#
# ------------
# description:
# ------------
# This script reads raw observation data, e.g., flux tower and gridded
# satellite data, and outputs CSV files formatted for the PostgreSQL database
# tables.  Currently, fluxdata.org CSV files, WATCH (WFDEI) 0.5 degree
# satellite netCDF files, and MODIS 0.05 HDF files, CRU netcdf, GLAS geoTIFF,
# and CRU dat files are supported.
#
# Only two of the three tables are output for flux tower data
# (var_list and data_set).
#
# ----------
# changelog:
# ----------
# 01. Updated input file path and output file name for bitnami [13.06.17]
# 02. Fixed IOError output (filename is FILE, not f) [13.06.17]
# 03. Reading data to make data_set table [13.06.17]
# 04. Changed output for both var_list and data_set to be restricted to
#     core-vars only [13.06.18]
# 05. Changed DATA class to FLUXDATA, because it's specific to fluxdata.org
#     station data [13.06.18]
# 06. Abbreviated the list of core variables to those with qc flags [13.07.04]
# 07. Added core variable qc flag dictionary [13.07.04]
# 08. Implemented output for data-set based on only observations [13.07.04]
# 09. Added "re" module for string searching [13.07.04]
# 10. Added scipy.io for netcdf handling [13.08.09]
# 11. Changed class "LINE" to "VAR" [13.08.09]
# 12. Added "SWdown" to coreVars and variableUnits hashes [13.08.09]
# 13. Moved main to "process_flux" function call to give room for gridded data
#     processing [13.08.09]
# 14. Added process_grid function for WATCH WFDEI [13.08.19]
# 15. Updated VAR class getStation function
#     separated it based on variable type 'flux' or 'grid' [13.08.19]
# 16. Changed WFD grid station naming scheme from lat+lon combo to
#     3 char + 6 digits (e.g., "WFD000123") [13.08.27]
# --> changed "WFD" to "HDF" [13.09.23]
# 17. Reordered processing of grids; spatial first (row-major) to get var_list,
#     then temporal to get data_set [13.08.27]
# 18. Updated get_flux_station to search filename for station ID [13.08.27]
# 19. Individualized data_set output files [13.08.27]
# 20. Added varlist_flag (to run varlist output once only) [13.08.27]
# 21. Deleted lathash and lonhash (not used) [13.08.30]
# 22. Fixed station parts error (moved from var_flag condition) [13.08.30]
# 23. Updated directory names (i.e., Projects) [13.09.09]
# 24. Added pxl_lat in process_watch function [13.09.23]
# 25. Added FAPAR to VAR class [13.09.27]
# 26. Added numpy to modules list [13.09.27]
# 27. Added pyhdf.SD to modules list [13.09.27]
# 28. Created MODISDATA class [13.09.27]
# 29. Created process_modis() function [13.09.27]
# 30. Changed vars in FLUXDATA class to vari [13.09.27]
# --> vars is a built-in function name
# 31. New output handling in process_flux for dat and var files [13.10.14]
# 32. Abbreviated core vars dict but kept the numbering [13.10.16]
# 33. Added core variable check to FLUXDATA class [13.10.16]
# 34. Added VPD to core variables list in units of kPa [13.11.06]
# 35. Added add_one_month() function [13.11.06]
# 36. Created process_cru() function [13.11.06]
# 37. Created CRUDATA class [13.11.06]
# 38. Added RH100 (canopy height, m) to core variables list [13.11.20]
# 39. Created process_glas() function for canopy height [13.11.20]
# 40. Corrected errors with stationid naming for GLAS gridded data [13.11.21]
# --> Created GLASDATA class to perform the calculations
# 41. Changed GLAS canopy height source data to TIF raster [13.11.28]
# 42. Added Image to module list [13.11.28]
# 43. Renamed process_cru to process_cru_vpd [14.01.17]
# 44. New process_cru for single netCDF file [14.01.17]
# 45. Added new variables to VAR class [14.01.17]
# --> Tc (monthly mean temp, deg. C)
# --> Pre (monthly precip, mm)
# 46. Fixed calc_vpd in CRUDATA class [14.01.17]
# --> skip 89 years (not 90)
# 47. Added new variable to VAR class [14.01.22]
# --> Cld (monthly cloudiness, % [0.0-100.0])
# 48. Added new variable to VAR class [14.01.23]
# --> Elv (CRU TS3.00 elevation data)
# 49. Created process_cru_elv function [14.01.23]
# 50. Deleted RH100 from GePiSaT database [14.02.04]
# --> renumbered CO2 at varid 21
# 51. Created process_alpha() [14.02.04]
# 52. Fixed error in VPD calc [14.02.18]
# --> 237.3 not 273.3
# 53. Added get_monthly_cru & get_time_index functions [14.02.19]
# 54. Deleted CRUDATA class [14.02.19]
# 55. Added calculate_vpd function [14.02.19]
# 56. Updated process_cru_vpd function [14.02.19]
# 57. Updated process_cru function [14.02.19]
# 58. General updates to style and inline comments [14.04.16]
# 59. More housekeeping [14.09.03]
# 60. Improved initializationo of vpd in calculate_vpd() [14.09.03]
# 61. PEP8 style fixes [15.11.18]
#
# -----
# todo:
# -----
# 1. VAR class:
#    a. add new dictionary for variable names?
#       * map the "_f" to more natural names
# 2. Improve calculate_vpd
#    * create clip and error field arrays
#    * calculate vpd in a single line computation
#
###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import glob
import Image
import os.path
import re
import numpy
from pyhdf import SD
from scipy.io import netcdf


###############################################################################
# CLASSES
###############################################################################
class VAR:
    """
    Name:     VAR
    Features: This class creates output lines for GePiSaT var_list database
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    msvIDX = ""        # specific identifier made up of station and var id
    stationID = ""     # station ID (from file name)
    varID = 0          # number associated with a specific variable
    varName = ""       # name of variable
    varUnit = ""       # unit of measure for variable
    varType = ""       # flux, met, or grid
    varCore = ""       # core (1) or non-core (0)

    # Core variable, and unit dictionaries:
    # * core variable 1--17 based on fluxdata.org
    #   (those core vars not used are commented out)
    # * additional core variables added for GePiSaT (18--26)
    coreVars = {
        'NEE_f': 1,      # Net ecosystem exchange (CO2 flux)
        #'GPP_f': 2,
        #'LE_f': 3,
        #'H_f': 4,
        #'G_f': 5,
        #'Ta_f': 6,
        #'Ts1_f': 7,
        #'Ts2_f': 8,
        #'VPD_f': 9,
        #'Precip_f': 10,
        #'SWC1_f': 11,
        #'SWC2_f': 12,
        #'WS_f': 13,
        #'Rg_f': 14,
        'PPFD_f': 15,    # Photosynthetic photon flux density
        #'Rn_f': 16,
        #'gsurf_f': 17,
        'SWdown': 18,    # WATCH shortwave downwelling solar radiation
        'FAPAR': 19,     # MODIS-based fraction of absorbed PAR
        'VPD': 20,       # CRU-based vapor pressure deficit
        'CO2': 21,       # NOAA sea-surface annual atm. CO2 concen.
        'Tc': 22,        # CRU monthly mean daily air temperature
        'Pre': 23,       # CRU monthly total precipitation
        'Cld': 24,       # CRU monthly cloudiness
        'Elv': 25,       # CRU 0.5 degree pixel centroid elevations
        'alpha': 26      # Cramer-Prentice bioclimatic moisture index
        }
    #
    variableUnits = {'NEE_f': 'umolCO2 m-2 s-1',
                     'GPP_f': 'umolCO2 m-2 s-1',
                     'LE_f': 'W m-2',
                     'H_f': 'W m-2',
                     'G_f': 'W m-2',
                     'Ta_f': 'deg C',
                     'Ts1_f': 'deg C',
                     'Ts2_f': 'deg C',
                     #'VPD_f': 'hPa',
                     #'Precip_f': 'mm',
                     'SWC1_f': '%',
                     'SWC2_f': '%',
                     'WS_f': 'm s-1',
                     'Rg_f': 'W m-2',
                     'PPFD_f': 'umol m-2 s-1',
                     'Rn_f': 'W m-2',
                     'gsurf_f': 'mmol m-2 s-1',
                     'SWdown': 'W m-2',
                     'FAPAR': 'NA',
                     'VPD': 'kPa',
                     'CO2': 'ppm',
                     'Tc': 'deg C',
                     'Pre': 'mm',
                     'Cld': '%',
                     'Elv': 'm',
                     'alpha': 'NA'}

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, fileName, variable, t):
        """
        Name:     VAR.__init__
        Input:    - string/tuple, filename or file naming tuple (fileName)
                  - string, variable name (variable)
                  - string, variable type, e.g., 'flux' or 'grid' (t)
        """
        # Set variable type and name (based on input):
        self.varType = t
        self.varName = variable
        #
        # Get variable ID, core boolean, and units:
        self.varID = self.getVarID(variable)
        self.varCore = self.isCore(variable)
        self.varUnit = self.getUnits(variable)
        #
        # Get station ID and msvidx based on variable type:
        if self.varType == 'flux':
            self.stationID = self.get_flux_station(fileName)
            self.msvIDX = self.getIDX(self.varID, self.stationID)
        elif self.varType == 'grid':
            # fileName has station prefix and ID:
            self.stationID = self.get_grid_station(fileName)
            self.msvIDX = self.getIDX(self.varID, self.stationID)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def get_flux_station(self, fileName):
        """
        Name:     VAR.get_flux_station
        Input:    string, filename with path (fileName)
        Output:   string, station name (sname)
        Features: Returns the flux station name from flux data file
        """
        # Initialize station name:
        sname = ""
        try:
            # Use regular expression search on filename
            # note: group(1) returns what is inside the search ()
            sname = re.search(
                '(\w{2}-\w{2,3})\.',
                os.path.basename(fileName)
                ).group(1)
        except AttributeError:
            print "Station name not found in file:", fileName
            #
        return sname

    def get_grid_station(self, fileName):
        """
        Name:     VAR.get_grid_station
        Input:    tuple, station prefix and station ID (fileName)
        Output:   string, station name (sname)
        Features: Returns the grid station name by appending prefix and station
                  ID parts together
        """
        # Parse grid information (sent as tuple):
        s_pre, s_id = fileName
        #
        # Format ID:
        sname = "%s%06d" % (s_pre, s_id)
        return sname

    def getVarID(self, variable):
        """
        Name:     VAR.getVarID
        Input:    string, variable name (variable)
        Output:   int, variable ID (vid)
        Features: Returns the variable ID based on the coreVars dictionary
        """
        if variable in self.coreVars:
            vid = self.coreVars[variable]
        else:
            vid = 0
        return vid

    def isCore(self, variable):
        """
        Name:     VAR.isCore
        Input:    string, variable name (variable)
        Output:   boolean
        Features: Returns boolean of whether the variable name is listed in
                  the coreVars dictionary
        """
        # Set boolean for core variable:
        if variable in self.coreVars:
            core = 1
        else:
            core = 0
        return core

    def getUnits(self, variable):
        """
        Name:     VAR.getUnits
        Input:    string, variable name (variable)
        Output:   string, variable units (v)
        Features: Returns the units for a variable in the variableUnits
                  dictionary
        """
        if variable in self.variableUnits:
            v = self.variableUnits[variable]
        else:
            v = ''
        return v

    def getIDX(self, vid, sid):
        """
        Name:     VAR.getIDX
        Input:    - int, variable ID (vid)
                  - string, station name (sid)
        Output:   string, msvidx (msv)
        Features: Returns the msvidx GePiSaT database entry based on the
                  station name and variable ID
        """
        msv = "%s.%02d" % (sid, vid)
        return msv

    def printLine(self, outfile):
        """
        Name:     VAR.printLine
        Input:    string, output file name and path (outfile)
        Output:   None.
        Features: Writes to file the contents for the var_list database table
        """
        # Don't print out timestamp variables or non-core variables:
        if self.varID != 0 and self.varCore == 1:
            # Create output line:
            outline = "%s,%s,%02d,%s,%s,%s,%s\n" % (
                self.msvIDX,
                self.stationID,
                self.varID,
                self.varName,
                self.varUnit,
                self.varType,
                self.varCore
                )
            #
            # Append data to output file:
            try:
                OUT = open(outfile, 'a')
                OUT.write(outline)
            except IOError:
                print "VAR Class: could not append to file", outfile
            else:
                OUT.close()


class FLUXDATA:
    """
    Name:     FLUXDATA
    Features: This class creates GePiSaT database data_set table entries based
              on fluxdata.org flux tower data files
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
        Name:     FLUXDATA.__init__
        Input:    - string, file name with path (fileName)
                  - list, variable names, e.g., parsed header row (v)
                  - list, observation data, e.g., parsed data row (d)
        Depends:  VAR class
        """
        # Calculate the number of fields:
        numfields = len(d)
        #
        # Strip time values and calculate the timestamp:
        # WARNING: the midnight hour time stamp of the last day of the
        # year has the new year category and DoY = 366 or 367 which
        # will cause problems in getTS.
        #
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
            #
        # Get timestamp:
        self.dateTime = self.getTS(self.timeVals)
        #
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
                #
                # Check that qc is for observation:
                if qcflag is "NEE_GPP_qc":
                    qcheck = self.check_gpp(qcvalue)
                else:
                    qcheck = self.check_vars(qcvalue)
                #
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
        Name:     FLUXDATA.getTS
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
        ts = (
            datetime.datetime(int(yr), 1, 1, 0, 0, 0)
            + datetime.timedelta(days=int(doy))
            + datetime.timedelta(minutes=mn)
            )
        return ts

    def check_gpp(self, qcf):
        """
        Name:     FLUXDATA.check_gpp
        Input:    int, quality control flag (qcf)
        Output:   int, data type (rval)
        Features: Returns a value based on the data type defined by the
                  quality control flag for the NEE and GPP flux tower variables
                  0: missing, 1: observed, 2: gap-filled
        """
        # Initialize return value:
        rval = -1
        #
        # CHECK 1: if value is missing:
        if (float(qcf) == -9999):
            rval = 0
        # CHECK 2: if value is original (obsevation):
        elif (float(qcf) == 1 or float(qcf) == 2):
            rval = 1
        # CHECK 3: if value has been gap-filled:
        elif (float(qcf) >= 3 and float(qcf) <= 6):
            rval = 2
        #
        # Return:
        return rval

    def check_vars(self, qcf):
        """
        Name:     FLUXDATA.check_vars
        Input:    int, quality control flag (qcf)
        Output:   int, data type (rval)
        Features: Returns a value based on the data type defined by the
                  quality control flag for state variables (not NEE or GPP)
                  0: missing, 1: observed, 2: gap-filled
        """
        # Initialize return value
        rval = -1
        # CHECK 1: if value is missing:
        if (float(qcf) == -9999):
            rval = 0
        # CHECK 2: if value is original (observation):
        elif (float(qcf) == 0):
            rval = 1
        # CHECK 3: if value has been gap-filled:
        elif (float(qcf) > 0 and float(qcf) <= 3):
            rval = 2
        #
        # Return:
        return rval

    def print_line(self, outfile):
        """
        Name:     FLUXDATA.print_line
        Input:    string, output file and path (outfile)
        Output:   None.
        Features: Writes to file the contents for the data_set database table
        """
        # Append data to output file:
        OUT = open(outfile, 'a')
        #
        # Create/write output lines:
        for i in range(len(self.data)):
            outline = "%s,%s,%s,%0.3f\n" % (
                self.msvIDX[i],
                self.stationID,
                self.dateTime,
                float(self.data[i])
                )
            OUT.write(outline)
        #
        # Close output file:
        OUT.close()


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
        # Create a variable class & save stationid and msvidx values:
        my_var = VAR(st_parts, val_parts[0], 'grid')
        self.station_id = my_var.stationID
        self.msv_idx = my_var.msvIDX
        #
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
        # Append to existing file:
        OUT = open(outfile, 'a')
        #
        # Create/write output lines:
        outline = "%s,%s,%s,%0.5f\n" % (
            self.msv_idx,
            self.station_id,
            self.data_time,
            self.data_value
            )
        OUT.write(outline)
        OUT.close()


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
        # Set datetime:
        self.data_time = t
        #
        # Get longitude and latitude from indices:
        (self.hdg_lon, self.hdg_lat) = self.get_lon_lat(x, y, self.hdg_res)
        #
        # Get HDG station ID:
        self.hdg_id = self.get_stationid()
        #
        # Get station id and msvidx:
        station_parts = ('HDG', self.hdg_id)
        my_var = VAR(station_parts, 'FAPAR', 'grid')
        self.station_id = my_var.stationID
        self.msv_idx = my_var.msvIDX
        #
        # Get 0.05 res data for this 0.5 pixel:
        my_points = self.get_foh_points(self.hdg_lon, self.hdg_lat, d)
        #
        # Check that there's data in the array:
        if len(my_points[~numpy.isnan(my_points)]) > 0:
            # Calculate the average EVI for 0.5 pixel:
            self.hdg_ave = my_points[~numpy.isnan(my_points)].mean()
            self.data_val = self.hdg_ave / 10000.0
            #
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
        #
        # Calculate lat, lon based on pixel index:
        lon = lon + x*r
        lat = lat - y*r
        #
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
            720.0 * (359.0 - ((90.0 - self.hdg_lat)/self.hdg_res - 0.5))
            + ((self.hdg_lon + 180.0)/self.hdg_res - 0.5)
            )
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
        #
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
        #
        # Calculate the binding box at 0.5 deg:
        westing = self.hdg_lon - 0.5*self.hdg_res
        northing = self.hdg_lat + 0.5*self.hdg_res
        #
        # Initialize centroid offsetting for foh_grid:
        foh_lon = westing + 0.5*self.foh_res
        foh_lat = northing - 0.5*self.foh_res
        #
        # Iterate over the 10x10 box:
        for y in xrange(10):
            lat = foh_lat - y*self.foh_res
            for x in xrange(10):
                lon = foh_lon + x*self.foh_res
                foh_grid.append((lon, lat))
        #
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
        #
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
        #
        return (int(x), int(y))

    def print_line(self, outfile):
        """
        Name:     MODISDATA.print_line
        Input:    string, output file name and path (outfile)
        Output:   None.
        Features: Writes to file the contents for data_set database table
        """
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
        OUT.close()


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
        for y in xrange(360):
            for x in xrange(720):
                # Determine lon, lat, and station ID:
                lon, lat = self.get_lon_lat(x, y, 0.5)
                st_id = self.get_stationid(lon, lat)
                #
                # Create variable class:
                station_parts = ('HDG', st_id)
                my_line = VAR(station_parts, 'RH100', 'grid')
                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                my_line.printLine(self.variable_file)
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                #
                # Retrieve 1km data for this 0.5 deg pixel:
                my_rh_data = self.get_1km_rh(lon, lat)
                #
                # Check that there's data in the array:
                if len(my_rh_data) > 0:
                    val = my_rh_data.mean()
                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Append to existing file:
                OUT = open(self.dataset_file, 'a')
                #
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
        #
        # Define half-degree and five one-hundredths resolutions:
        hdg_res = 0.5
        okm_res = 1.0/120.0
        #
        # Calculate the binding box at 0.5 deg:
        westing = hdg_lon - 0.5*hdg_res
        northing = hdg_lat + 0.5*hdg_res
        #
        # Initialize centroid offsetting for 1 km grid:
        okm_lon = westing + 0.5*okm_res
        okm_lat = northing - 0.5*okm_res
        #
        # Iterate over the 60x60 box:
        for y in xrange(60):
            lat = okm_lat - y*okm_res
            for x in xrange(60):
                lon = okm_lon + x*okm_res
                okm_grid.append((lon, lat))
        #
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
        #
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
        #
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
            720.0 * (359.0 - ((90.0 - lat)/0.5 - 0.5))
            + ((lon + 180.0)/0.5 - 0.5)
            )
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
        #
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
        #
        return okm_indices


###############################################################################
# FUNCTIONS
###############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime date (dt0)
    Output:   datetime date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3


def get_modis_ts(f):
    """
    Name:     get_modis_ts
    Input:    string, file name (f)
    Output:   datetime.date
    Features: Returns a timestamp based on the MODIS filenaming scheme
              'MYD13C2.AYEARDOY...'
    """
    f_name = os.path.basename(f)
    try:
        f_date = re.search('A(\d{7})\.', f_name).group(1)
    except AttributeError:
        print "Search failed for date in file", f_name
    else:
        f_year = int(f_date[0:4])
        f_doy = int(f_date[4:7]) - 1
        f_timestamp = (
            datetime.date(f_year, 1, 1) +
            datetime.timedelta(days=f_doy)
            )
        return f_timestamp


def get_monthly_cru(d, ct, v):
    """
    Name:     get_monthly_cru
    Input:    - string, directory to CRU netcdf file (d)
              - datetime.date, current month datetime object (ct)
              - string, variable of interest (v)
    Output:   numpy nd.array
    Depends:  get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and
              variable of interest (e.g., cld, pre, tmp)
    """
    # Search directory for netCDF file:
    my_file = glob.glob(d + "*" + v + ".dat.nc")[0]
    #
    if my_file:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")
        #
        # Save data for variables of interest:
        # NOTE: for CRU TS 3.21:
        #       variables: 'lat', 'lon', 'time', v
        #       where v is 'tmp', 'pre', 'cld'
        # LAT:  -89.75 -- 89.75
        # LON:  -179.75 -- 179.75
        # TIME:
        #       units: days since 1900-1-1
        #       shape: (1344,)
        #       values: mid-day of each month (e.g., 15th or 16th)
        # DATA:
        #       'cld' units = %
        #       'pre' units = mm
        #       'tmp' units = deg. C
        #       Missing value = 9.96e+36
        # Save the base time stamp:
        bt = datetime.date(1900, 1, 1)
        #
        # Read the time data as array:
        f_time = f.variables['time'].data
        #
        # Find the time index for the current date:
        ti = get_time_index(bt, ct, f_time)
        #
        # Get the spatial data for current time:
        f_data = f.variables[v].data[ti]
        f.close()
        return f_data


def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime.date, base timestamp
              - datetime.date, current timestamp
              - numpy nd.array, days since base timestamp
    Output:   int
    Features: Finds the index in an array of CRU TS days for a given timestamp
    """
    # For CRU TS 3.21, the aot is indexed for mid-month days, e.g. 15--16th
    # therefore, to make certain that ct index preceeds the index for the
    # correct month in aot, make the day of the current month less than
    # the 15th or 16th (i.e., replace day with '1'):
    ct = ct.replace(day=1)
    #
    # Calculate the time difference between ct and bt:
    dt = (ct - bt).days
    #
    # Append dt to the aot array:
    aot = numpy.append(aot, [dt, ])
    #
    # Find the first index of dt in the sorted array:
    idx = numpy.where(numpy.sort(aot) == dt)[0][0]
    return idx


def calculate_vpd(tmin, tmax, vap):
    """
    Name:     calculate_vpd
    Input:    - numpy nd.array, mean monthly min daily air temp, deg C (tmin)
              - numpy nd.array, mean monthly max daily air temp, deg C (tmax)
              - numpy nd.array, mean monthly vapor pressure, hPa (vap)
    Output:   numpy nd.array, mean monthly vapor pressure deficit, kPa (vpd)
    Features: Returns an array of mean monthly vapor pressure deficit
    Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure
              Calculation Methods, in Evaporation and Evapotranspiration:
              Measurements and Estimations, Springer, London.
                vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
                where:
                    tc = average daily air temperature, deg C
                    ea = actual vapor pressure, kPa
    """
    # Initialize array:
    # NOTE: maintains missing value
    vpd = -9999.0*(numpy.zeros(shape=(360, 720)) + 1)
    #
    # Iterate through each data point:
    lat, lon = tmin.shape
    #
    for y in xrange(lat):
        for x in xrange(lon):
            tm = tmin[y, x]
            tx = tmax[y, x]
            ea = vap[y, x]
            if tm < 1.e6 and tx < 1.e6 and ea < 1.e6:
                to = 0.5*(tm + tx)
                vpd[y, x] = (
                    0.611*numpy.exp((17.27*to)/(to + 237.3)) - 0.10*ea
                )
    return vpd


def process_flux(d):
    """
    Name:     process_flux
    Input:    string, input file directory (d)
    Output:   None.
    Features: Processes flux tower data files into variable list and data set
              table output files
    Depends:  writeout
    """
    # Search directory for fluxdata CSV files:
    my_dir = d
    my_ext = "*allvars.csv"
    my_files = glob.glob(my_dir + my_ext)
    #
    # Prepare output files:
    my_var_out = "Fluxdata_Var-List_all.csv"
    var_outfile = my_dir + my_var_out
    my_dat_out = "Fluxdata_Data-Set_test.csv"
    dat_outfile = my_dir + my_dat_out
    #
    # Prepare header lines and write to file:
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_outfile, var_headerline)
    dat_headerline = "msvidx,stationid,datetime,data\n"
    writeout(dat_outfile, dat_headerline)
    #
    # Create a list for book keeping which stations have been var-processed:
    var_out_list = []
    var_flag = 1
    #
    # Statically set data-set processing flag
    # in case you want to just process var-list
    dat_flag = 1
    #
    # Check that files were found:
    if my_files:
        # Read through each file:
        for doc in my_files:
            if os.path.isfile(doc):
                try:
                    # Try to get filename prefix:
                    my_station = re.search(
                        '^\S{6}',
                        os.path.basename(doc)
                        ).group(0)
                    my_dat = re.search(
                        '^\S{6}\.{1}\d{4}',
                        os.path.basename(doc)
                        ).group(0)
                    my_dat_out = my_dat + "_Data-Set.csv"
                    my_var_out = my_station + "_Var-List.csv"
                except AttributeError:
                    print "Could not read file prefix from", doc
                else:
                    # Initialize data-set outfile:
                    if dat_flag:
                        dat_outfile = my_dir + my_dat_out
                        writeout(dat_outfile, dat_headerline)
                    #
                    # Check to see if station has been var'ed:
                    if my_station in var_out_list:
                        var_flag = 0
                    else:
                        # Initialize var outfile:
                        var_outfile = my_dir + my_var_out
                        writeout(var_outfile, var_headerline)
                        #
                        # Add var to var_list and set flag to true:
                        var_out_list.append(my_station)
                        var_flag = 1
                #
                try:
                    # Try to open file for reading:
                    f = open(doc, 'r')
                    #
                    # Read the headerline and strip whitespace from the end:
                    header = f.readline()
                    header = header.rstrip()
                    #
                    # Read the remaining content:
                    content = f.readlines()
                    #
                    # Close the file after done reading:
                    f.close()
                    #
                except IOError:
                    print "Could not read file: " + doc
                else:
                    # Parse header line:
                    parts = header.split(',')
                    #
                    # ~~~~~~~~~~ VAR-LIST ~~~~~~~~~~ #
                    # For each flux variable, produce the data row:
                    if var_flag:
                        for var in parts:
                            my_line = VAR(doc, var, 'flux')
                            my_line.printLine(var_outfile)
                            #
                    #
                    # ~~~~~~~~~~ DATA-SET ~~~~~~~~~~ #
                    # For each timestamp, produce variable data
                    if dat_flag:
                        for row in content:
                            row = row.rstrip()
                            dataset = row.split(',')
                            #
                            # Create data class and print out lines:
                            my_data = FLUXDATA(doc, parts, dataset)
                            my_data.print_line(dat_outfile)
    else:
        print "No files found in directory:", my_dir


def process_watch(d, voi):
    """
    Name:     process_watch
    Input:    - string, input file directory (d)
              - string, variable of interest (voi)
    Output:   None.
    Features: Processes WATCH WFDEI netCDF files into variable list and data
              set table output files
    Depends:  writeout
    """
    # Search directory for WATCH netCDF files:
    my_dir = d
    my_ext = "*.nc"
    my_files = glob.glob(my_dir + my_ext)
    #
    # Prepare var output file:
    my_var_out = "WFDEI_Var-List_test.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_outfile, var_headerline)
    #
    # Flag for varlist:
    varlist_flag = 1
    #
    if my_files:
        # Read through each file:
        for doc in my_files:
            if os.path.isfile(doc):
                try:
                    # Try to open file for reading:
                    f = netcdf.NetCDFFile(doc, "r")
                    #
                    # Read the year and month values from filename (YYYYMM):
                    yr_mo = re.search(
                        '_(\d{6})\.',
                        os.path.basename(f.filename)
                        ).group(1)
                    this_year = int(yr_mo[0:4])
                    this_month = int(yr_mo[4:6])
                except IOError:
                    print "Could not read file: " + doc
                except AttributeError:
                    print "Year and month not retrievable from file", doc
                except ValueError:
                    print "Year and/or month not numbers in file", doc
                else:
                    # Prepare data set output file:
                    my_dat_out = "WFDEI_Data-Set_%s.csv" % yr_mo
                    dat_outfile = my_dir + my_dat_out
                    dat_headerline = "msvidx,stationid,datetime,data\n"
                    writeout(dat_outfile, dat_headerline)
                    #
                    # Save the shape values of each variable:
                    sh_day, sh_lat, sh_lon = f.variables[voi].shape
                    #
                    # Iterate through each lat:lon pair
                    # * row-major ordering from bottom left
                    #  x (longitude): 0...719
                    #  y (latitude): 0...359
                    for y in xrange(sh_lat):
                        # Save latitude:
                        pxl_lat = f.variables['lat'].data[y]
                        #
                        for x in xrange(sh_lon):
                            # Calc station ID:
                            st_id = 720*y + x
                            station_parts = ('HDG', st_id)
                            #
                            if varlist_flag:
                                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                                my_line = VAR(station_parts, voi, 'grid')
                                my_line.printLine(var_outfile)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                            #
                            # Iterate through each day
                            for t in xrange(sh_day):
                                # Get timestamp for this day:
                                this_day = t+1
                                time_stamp = datetime.date(
                                    this_year,
                                    this_month,
                                    this_day
                                    )
                                #
                                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                                # Read SWdown for each pixel
                                # * NOTE 1: variable has five decimal places
                                # * NOTE 2: missing values are equal to ~1e20
                                # * NOTE 3: Antarctica is < -60 latitude
                                pxl_swr = f.variables[voi].data[t, y, x]
                                if pxl_swr < 1.0e6 and pxl_lat > -60:
                                    # Process pixel
                                    obs_parts = (voi, pxl_swr, time_stamp)
                                    my_data = WATCHDATA(
                                        station_parts,
                                        obs_parts
                                        )
                                    my_data.print_line(dat_outfile)
                                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                    # Close the var list flag after processing first doc:
                    varlist_flag = 0
                    f.close()


def process_modis(d, voi):
    """
    Name:     process_modis
    Input:    - string, input/output file directory (d)
              - string, variable of interest (voi)
    Output:   None.
    Features: Processes 0.05 degree MODIS EVI from HDF files into 0.5 degree
              FAPAR and saves to variable list and data set table output files
    Depends:  writeout
    """
    # Search directory for MODIS HDF files:
    my_dir = d
    my_ext = "*.hdf"
    my_files = glob.glob(my_dir + my_ext)
    #
    # Prepare var output file:
    my_var_out = "MODIS_Var-List_fapar.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_outfile, var_headerline)
    #
    # Flag for varlist:
    varlist_flag = 1
    #
    if my_files:
        # Read through each file:
        for doc in my_files:
            try:
                # Try to open file for reading:
                my_hdf = SD.SD(doc)
            except:
                print "Unexpected error opening file", doc
                data = None
            else:
                # Check that VOI is in the dataset:
                if voi in my_hdf.datasets().keys():
                    # Pull data from dataset into array:
                    d_select = my_hdf.select(voi)
                    data = d_select.get()
                else:
                    print "Could not open dataset", voi
                    data = None
                    #
                # Close HDF file and return data:
                my_hdf.end()
                #
            # Continue processing if data was retrieved:
            if data.any():
                # Get timestamp from filename:
                my_ts = get_modis_ts(doc)
                #
                # Prepare data set output file:
                my_dat_out = "MODIS_Data-Set_%s.csv" % my_ts
                dat_outfile = my_dir + my_dat_out
                dat_headerline = "msvidx,stationid,datetime,data\n"
                writeout(dat_outfile, dat_headerline)
                #
                # Iterate through each 0.5 deg lat:lon pair
                # x (longitude): 0...719
                # y (latitude): 0...359
                for y in xrange(360):
                    for x in xrange(720):
                        my_modis = MODISDATA(x, y, my_ts, data)
                        station_parts = ('HDG', my_modis.hdg_id)
                        #
                        if varlist_flag:
                            # ~~~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~~~ #
                            my_line = VAR(station_parts, 'FAPAR', 'grid')
                            my_line.printLine(var_outfile)
                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                        #
                        # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                        my_modis.print_line(dat_outfile)
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # Close varflag after processing first doc:
            varlist_flag = 0


def process_cru_elv(d):
    """
    Name:     process_cru_elv
    Input:    string, input/output file directory (d)
    Output:   None.
    Features: Processes CRU TS 3.00 dat file associated with elevation into
              variable list and data set table output files
    Depends:  writeout
    """
    my_dir = d
    #
    # Search directory for CRU dat file:
    my_file = glob.glob(my_dir + "*dat")[0]
    #
    # Open and read dat file:
    # NOTE: data is read into an array with shape (360, 720)
    #      'lat' goes from -89.75 -- 89.75 (south to north)
    #      'lon' goes from -179.75 -- 179.75 (east to west)
    f = numpy.loadtxt(my_file)
    (sh_lat, sh_lon) = f.shape
    #
    # Assign time stamp as CRU TS 3.00 date:
    time_stamp = datetime.date(2006, 6, 1)
    #
    # Prepare var output file:
    my_var_out = "CRU_Var-List_elv.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)
    #
    # Prepare data set output file:
    my_dat_out = "CRU_Data-Set_elv.csv"
    dat_outfile = my_dir + my_dat_out
    dat_headerline = "msvidx,stationid,datetime,data\n"
    writeout(dat_outfile, dat_headerline)
    #
    # Iterate through data:
    for y in xrange(sh_lat):
        for x in xrange(sh_lon):
            # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
            st_id = 720*y + x
            station_parts = ('HDG', st_id)
            my_line = VAR(station_parts, 'Elv', 'grid')
            my_line.printLine(var_outfile)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            #
            elv = f[y, x]
            # Set no-data value:
            if elv > -500:
                # Append to existing file:
                OUT = open(dat_outfile, 'a')
                # Create/write output line:
                outline = "%s,%s,%s,%0.5f\n" % (
                    my_line.msvIDX,
                    my_line.stationID,
                    time_stamp,
                    elv
                    )
                OUT.write(outline)
                OUT.close()


def process_cru_vpd(cru_dir, my_dir):
    """
    Name:     process_cru_vpd
    Input:    - string, directory for CRU input files (cru_dir)
              - string, output file directory (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF files (tmn, tmx, vap) to calculate
              monthly VPD and save variable list and data set table output
              files
    Depends:  - writeout
              - add_one_month
    """
    # Define the start and end dates you want to process (2002-2006):
    start_date = datetime.date(2002, 1, 1)
    end_date = datetime.date(2007, 1, 1)
    #
    # Set flag for varlist:
    varlist_flag = 1
    #
    # Prepare var output file:
    my_var_out = "CRU_Var-List_vpd.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)
    #
    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # Open and read netcdf files in the file directory:
        tmn = get_monthly_cru(cru_dir, cur_date, "tmn")
        tmx = get_monthly_cru(cru_dir, cur_date, "tmx")
        vap = get_monthly_cru(cru_dir, cur_date, "vap")
        #
        # Calculate VPD & save shape:
        vpd = calculate_vpd(tmn, tmx, vap)
        (sh_lat, sh_lon) = vpd.shape
        #
        # Prepare data set output file:
        my_dat_out = "CRU_Data-Set_vpd_%s.csv" % cur_date
        dat_outfile = my_dir + my_dat_out
        dat_headerline = "msvidx,stationid,datetime,data\n"
        writeout(dat_outfile, dat_headerline)
        #
        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in xrange(sh_lat):
            for x in xrange(sh_lon):
                # Calc station ID:
                st_id = 720*y + x
                station_parts = ('HDG', st_id)
                #
                my_line = VAR(station_parts, 'VPD', 'grid')
                #
                if varlist_flag:
                    # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                    my_line.printLine(var_outfile)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                #
                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Read VPD for each pixel
                # * NOTE: missing VPD are equal to -9999
                pxl_vpd = vpd[y, x]
                if pxl_vpd > -9999.0:
                    # Append to existing file:
                    OUT = open(dat_outfile, 'a')
                    #
                    # Create/write output line:
                    outline = "%s,%s,%s,%0.5f\n" % (
                        my_line.msvIDX,
                        my_line.stationID,
                        cur_date,
                        pxl_vpd
                        )
                    OUT.write(outline)
                    OUT.close()
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Turn off varlist flag after processing first month:
        varlist_flag = 0
        #
        # Increment cur_date:
        cur_date = add_one_month(cur_date)


def process_cru(v, cru_dir, my_dir):
    """
    Name:     process_cru
    Input:    - string, variable name, e.g., tmp, pre, cld (v)
              - string, directory name for CRU TS data file (cru_dir)
              - string, directory for output files (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF file by month into variable list and data
              set table output files
    Depends:  - writeout
              - add_one_month
    """
    # Define the start and end dates you want to process (2002-2006):
    start_date = datetime.date(2002, 1, 1)
    end_date = datetime.date(2007, 1, 1)
    #
    # Set flag for varlist:
    varlist_flag = 1
    #
    # Prepare var output file:
    my_var_out = "CRU_Var-List_" + v + ".csv"
    var_outfile = my_dir + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)
    #
    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # # Open and read netcdf files in the file directory:
        my_data = get_monthly_cru(cru_dir, cur_date, v)
        (sh_lat, sh_lon) = my_data.shape
        #
        # Prepare data set output file:
        my_dat_out = "CRU_Data-Set_%s_%s.csv" % (v, cur_date)
        dat_outfile = my_dir + my_dat_out
        dat_headerline = "msvidx,stationid,datetime,data\n"
        writeout(dat_outfile, dat_headerline)
        #
        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in xrange(sh_lat):
            for x in xrange(sh_lon):
                # Calc station ID:
                st_id = 720*y + x
                station_parts = ('HDG', st_id)
                #
                if v == 'tmp':
                    my_var = "Tc"
                elif v == 'pre':
                    my_var = "Pre"
                elif v == 'cld':
                    my_var = "Cld"
                my_line = VAR(station_parts, my_var, 'grid')
                #
                if varlist_flag:
                    # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                    my_line.printLine(var_outfile)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                #
                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Read each each pixel
                # * NOTE: missing values are ~ 1e7
                pxl = my_data[y, x]
                if pxl < 1.e6:
                    # Append to existing file:
                    OUT = open(dat_outfile, 'a')
                    #
                    # Create/write output line:
                    outline = "%s,%s,%s,%0.3f\n" % (
                        my_line.msvIDX,
                        my_line.stationID,
                        cur_date,
                        pxl
                        )
                    OUT.write(outline)
                    OUT.close()
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Turn off varlist flag after processing first month:
        varlist_flag = 0
        #
        # Increment cur_date:
        cur_date = add_one_month(cur_date)


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
    #
    # Open file for reading and save variable of interest:
    im = Image.open(my_file).getdata()
    #
    # Save the shape values of each for iteration purposes
    sh_lat, sh_lon = im.size
    #
    # Prepare var output file:
    my_var_out = "GLAS_Var-List_rh100.csv"
    var_outfile = d + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)
    #
    # Prepare data set output file:
    my_dat_out = "GLAS_Data-Set_rh100_2005.csv"
    dat_outfile = d + my_dat_out
    dat_headerline = "msvidx,stationid,datetime,data\n"
    writeout(dat_outfile, dat_headerline)
    #
    # Initialize and run GLAS class:
    my_glas = GLASDATA(var_outfile, dat_outfile, im)
    my_glas.process_data()


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
    #
    # Prepare var output file:
    my_var_out = "Var-List_alpha.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    writeout(var_outfile, var_headerline)
    #
    # Flag for varlist:
    varlist_flag = 1
    #
    # Read through files:
    if my_files:
        for doc in my_files:
            try:
                # Load file, skipping header lines:
                f = numpy.loadtxt(doc, skiprows=6)
                #
                # Read timestamp from filename:
                yr_mo = re.search(
                    '_(\d{4}-\d{2}-\d{2})\.',
                    os.path.basename(doc)
                    ).group(1)
                this_year = int(yr_mo.split('-')[0])
                this_month = int(yr_mo.split('-')[1])
                my_ts = datetime.date(this_year, this_month, 1)
            except:
                print "Error reading file, ", doc
            else:
                # Save the data shape (360x720):
                (sh_lat, sh_lon) = f.shape
                #
                # Prepare data set output file:
                my_dat_out = "Data-Set_alpha_%s.csv" % my_ts
                dat_outfile = my_dir + my_dat_out
                dat_headerline = "msvidx,stationid,datetime,data\n"
                writeout(dat_outfile, dat_headerline)
                #
                # Iterate through data file:
                # NOTE: ASCII raster is organized in 360 rows and 720 cols
                # * rows ordered from 89.75 to -89.75 lat (north > south)
                # * cols ordered from -179.75 to 179.75 lon (east > west)
                for y in xrange(sh_lat):
                    # Reverse order latitude index for station numbering:
                    i = 359 - y
                    for x in xrange(sh_lon):
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
                            #
                            # ~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~ #
                            # Append to existing file:
                            OUT = open(dat_outfile, 'a')
                            #
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
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

###############################################################################
# MAIN PROGRAM
###############################################################################
# Process flux data:
flux_dir = (
    "/Users/twdavis/Projects/"
    "data/flux/fluxdata.org/allTowers_indYears/"
    "2006/"
    )
#process_flux(flux_dir)

# Process WATCH data:
watch_dir = (
    "/Users/twdavis/Projects"
    "/data/watch"
    "/SWdown_daily_WFDEI/"
    )
watch_voi = 'SWdown'
#process_watch(watch_dir, watch_voi)

# Process MODIS data:
modis_dir = (
    "/Users/twdavis/Projects"
    "/data/modis/vi_cgm_monthly"
    "/aqua/"
    )
modis_voi = "CMG 0.05 Deg Monthly EVI"
#process_modis(modis_dir, modis_voi)

cru_dir = (
    "/usr/local/share/database/cru/"
    )
out_dir = (
    "/home/user/Projects/gepisat/data/vpd/"
    )
#process_cru_vpd(cru_dir, out_dir)
#process_cru_elv(cru_dir)
#process_cru("cld",cru_dir,out_dir)

glas_dir = (
    "/home/user/Projects/gepisat/data/GLAS/"
    )
#process_glas(glas_dir)

alpha_dir = (
    "/home/user/Projects/gepisat/data/alpha/raster_files/"
    )
process_alpha(alpha_dir)
