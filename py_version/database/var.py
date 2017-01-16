#!/usr/bin/python
#
# var.py
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
import logging
import re
import os


###############################################################################
# CLASSES
###############################################################################
class VAR:
    """
    Name:     VAR
    Features: This class creates output lines for GePiSaT var_list database
    History:  Version 3.0
              - added FLUXNET 2015 variables for NEE and PPFD [17.01.15]

    @TODO:  add new dictionary for variable names?
            (e.g., map the "_f" to more natural names)
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
        'NEE_f': 1,        # Net ecosystem exchange (CO2 flux); FLUXNET 2012
        'NEE_VUT_REF': 2,  # Net ecosystem exchange (CO2 flux); FLUXNET 2015
        # 'LE_f': 3,
        # 'H_f': 4,
        # 'G_f': 5,
        # 'Ta_f': 6,
        # 'Ts1_f': 7,
        # 'Ts2_f': 8,
        # 'VPD_f': 9,
        # 'Precip_f': 10,
        # 'SWC1_f': 11,
        # 'SWC2_f': 12,
        # 'WS_f': 13,
        # 'Rg_f': 14,
        'PPFD_f': 15,    # Photosynthetic photon flux density; FLUXNET 2012
        'PPFD_IN': 16,   # Photosynthetic photon flux density; FLUXNET 2015
        # 'Rn_f': 17,
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
                     'NEE_VUT_REF': 'umolCO2 m-2 s-1',
                     'LE_f': 'W m-2',
                     'H_f': 'W m-2',
                     'G_f': 'W m-2',
                     'Ta_f': 'deg C',
                     'Ts1_f': 'deg C',
                     'Ts2_f': 'deg C',
                     # 'VPD_f': 'hPa',
                     # 'Precip_f': 'mm',
                     'SWC1_f': '%',
                     'SWC2_f': '%',
                     'WS_f': 'm s-1',
                     'Rg_f': 'W m-2',
                     'PPFD_f': 'umol m-2 s-1',
                     'PPFD_IN': 'W m-2',
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
        # Create a class logger
        self.logger = logging.getLogger("VAR")
        self.logger.info("VAR class initialized")

        # Set variable type and name (based on input):
        self.varType = t
        self.varName = variable

        # Get variable ID, core boolean, and units:
        self.varID = self.getVarID(variable)
        self.varCore = self.isCore(variable)
        self.varUnit = self.getUnits(variable)

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
            self.logger.error("Station name not found in file: %s", fileName)
        else:
            self.logger.info("Found station %s", sname)
        finally:
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

            # Append data to output file:
            try:
                OUT = open(outfile, 'a')
                OUT.write(outline)
            except IOError:
                self.logger.error("could not append to file %s", outfile)
            else:
                OUT.close()
