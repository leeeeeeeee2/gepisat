#!/usr/bin/python
#
# model.py
#
# VERSION 2.2.0-dev
#
# LAST UPDATED: 2016-05-13
#
# ---------
# citation:
# ---------
# I. C. Prentice, T. W. Davis, X. M. P. Gilbert, B. D. Stocker, B. J. Evans,
# H. Wang, and T. F. Keenan, "The Global ecosystem in Space and Time (GePiSaT)
# Model of the Terrestrial Biosphere," (in progress).

###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import logging
import os

import numpy

from const import kPo
from const import kTo
from const import kL
from const import kG
from const import kR
from const import kMa
from db_setup import connectSQL
from utilities import init_summary_dict
from utilities import grid_centroid


###############################################################################
# CLASSES
###############################################################################
class DATA(object):
    """
    Name:     DATA
    Features: This class performs the data IO
    History   Version 2.2.0-dev
              - class inherits object for Python2/3 compatibility [16.05.13]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     DATA.__init__
        Input:    None.
        Features: Initialize the class
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("DATA class initialized")

        # Hidden properties:
        self._outputdir = "."
        self._stations = {}
        self._summaryfile = "summary.txt"

        # File handler (change if necessary)
        self.logger.debug("Initializing new file handler")
        self.fh = GPSQL()

        # Summary value dictionary:
        self.logger.debug("Initializing summary value dictionary")
        self.summary = init_summary_dict()

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Property Definitions
    # ////////////////////////////////////////////////////////////////////////
    @property
    def lue_file(self):
        """LUE summary file with output directory"""
        return os.path.join(self._outputdir, self._luefile)

    @lue_file.setter
    def lue_file(self, val):
        if isinstance(val, str):
            self.logger.debug("LUE file set to '%s'", val)
            self._luefile = val
        else:
            try:
                temp_val = str(val)
            except:
                self.logger.exception("Failed to read LUE file name")
                raise
            else:
                self.logger.debug("LUE file set to '%s'", temp_val)
                self._luefile = temp_val

    @property
    def output_dir(self):
        """Directory for saving export files."""
        return self._outputdir

    @output_dir.setter
    def output_dir(self, val):
        if isinstance(val, str):
            if os.path.isdir(val):
                self.logger.debug("Output directory set to %s", val)
                self._outputdir = val
            else:
                self.logger.error("Directory %s does not exist!", val)
                raise IOError("Output directory does not exist!")
        else:
            try:
                temp_val = str(val)
            except:
                self.logger.exception("Failed to read output directory.")
                raise
            else:
                if os.path.isdir(temp_val):
                    self.logger.debug("Output directory set to %s", temp_val)
                    self._outputdir = temp_val
                else:
                    self.logger.error("Directory %s does not exist!", temp_val)
                    raise IOError("Output directory does not exist!")

    @property
    def stations(self):
        """List of station names"""
        return sorted(list(self._stations.keys()))

    @stations.setter
    def stations(self, vals):
        if not isinstance(vals, list):
            self.logger.error("Stations must be given as a list!")
            raise TypeError("Stations must be given as a list!")
        else:
            for station in vals:
                self._stations[station] = {}
                lue_file = "%s_%s.txt" % (station, "LUE")
                lue_path = os.path.join(self.output_dir, lue_file)
                gpp_file = "%s_%s.txt" % (station, "GPP-daily")
                gpp_path = os.path.join(self.output_dir, gpp_file)
                self._stations[station]["lue"] = lue_path
                self._stations[station]["gpp"] = gpp_path

    @property
    def summary_file(self):
        """Summary file name with output directory"""
        return os.path.join(self._outputdir, self._summaryfile)

    @summary_file.setter
    def summary_file(self, val):
        if isinstance(val, str):
            self.logger.debug("Summary file set to '%s'", val)
            self._summaryfile = val
        else:
            try:
                temp_val = str(val)
            except:
                self.logger.exception("Failed to read summary file name")
                raise
            else:
                self.logger.debug("Summary file set to '%s'", temp_val)
                self._summaryfile = temp_val

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def create_summary_file(self, fname):
        """
        Name:     DATA.create_summary_file
        Inputs:   str, summary file name (fname)
        Outputs:  None.
        Features: Creates a summary file with the given name in the specified
                  output directory
        """
        self.summary_file = fname

        if os.path.isfile(self.summary_file):
            self.logger.warning("Overwriting '%s'", self.summary_file)
        try:
            SFILE = open(self.summary_file, 'w')
        except IOError:
            self.logger.exception("Failed to open '%s'", self.summary_file)
            raise
        else:
            # Create summary file header line:
            summary_header = ""
            for i in sorted(list(self.summary.keys())):
                field = self.summary[i]["name"]
                summary_header += str(field)
                summary_header += ","
            summary_header = summary_header.rstrip(",")
            summary_header += "\n"
            SFILE.write(summary_header)
        finally:
            SFILE.close()

    def find_date_range(self, station):
        """
        Name:     DATA.find_date_range
        Inputs:   str, station name (station)
        Outputs:  tuple, starting and ending dates
        Features: Searches for starting and ending dates for a given station's
                  available data
        """
        try:
            (sd, ed) = self.fh.get_date_range(station)
        except IOError:
            self.logger.exception(("Search failed. Check to make certain that "
                                   "the file handler has access to your input "
                                   "data."))
            raise
        except:
            self.logger.exception(("Encountered unknown error during search. "
                                   "Stopping."))
            raise
        else:
            self.logger.debug("Found date range %s to %s" % (sd, ed))
            return (sd, ed)

    def find_elv_pressure(self, station):
        """
        @TODO
        """
        elv, atm = self.fh.get_elap(station)

    def find_station_grid(self, station):
        """
        Name:     DATA.find_station_grid
        Inputs:   str, station name (station)
        Outputs:  str, grid name
        Features: Returns the grid name for a given station
        """
        try:
            grid = self.fh.flux_to_grid(station)
        except IOError:
            self.logger.exception(("Search failed. Check to make certain that "
                                   "the file handler has access to your input "
                                   "data."))
            raise
        except:
            self.logger.exception(("Encountered unknown error during search. "
                                   "Stopping."))
            raise
        else:
            self.logger.debug("Found grid %s for station %s" % (grid, station))
            return grid

    def find_stations(self):
        """
        Name:     DATA.find_stations
        Inputs:   None.
        Outputs:  None.
        Features: Searches for a list of flux tower stations
        """
        self.logger.debug("Searching for stations...")
        self.stations = self.fh.get_stations()
        self.logger.debug("Found %d stations!", len(self.stations))

    def get_gpp_file(self, station):
        """
        Name:     DATA.get_gpp_file
        Inputs:   str, station name (station)
        Outputs:  str, GPP station file
        Features: Returns the GPP file for a given station
        """
        if station in self.stations:
            return self._stations[station]["gpp"]
        else:
            self.logger.error("Could not find station '%s'", station)
            raise ValueError("Station '%s' not defined." % (station))

    def get_lue_file(self, station):
        """
        Name:     DATA.get_lue_file
        Inputs:   str, station name (station)
        Outputs:  str, LUE station file
        Features: Returns the LUE file for a given station
        """
        if station in self.stations:
            return self._stations[station]["lue"]
        else:
            self.logger.error("Could not find station '%s'", station)
            raise ValueError("Station '%s' not defined." % (station))

    def get_point_data(self, station, var, tp):
        """
        Name:     DATA.get_point_data
        Inputs:   str, station name (station)
                  str, variable name (var)
                  datetime.date, time point (tp)
        Output:   float/numpy.ndarray (my_result)
        Features: Returns data of a given variable for a given station at a
                  given time point
        """
        my_result = self.fh.get_data_point(station, var, tp)
        return my_result

    def write_summary(self):
        """
        Name:     DATA.create_summary_file
        Inputs:   None.
        Outputs:  None.
        Features: Writes summary values to file
        """
        if not os.path.isfile(self.summary_file):
            self.logger.error("Summary file has not been created!")
            raise
        else:
            try:
                SFILE = open(self.summary_file, 'a')
            except IOError:
                self.logger.exception("Failed to open '%s'", self.summary_file)
                raise
            else:
                # Create summary file header line:
                output_line = ""
                for i in sorted(list(self.summary.keys())):
                    value = self.summary[i]["val"]
                    output_line += str(value)
                    output_line += ","
                output_line = output_line.rstrip(",")
                output_line += "\n"
                SFILE.write(output_line)
            finally:
                SFILE.close()


class GPSQL(object):
    """
    Name:     GPSQL
    Features: This class performs data IO with the GePiSaT PSQL database
    History:  Version 2.2.0-dev
              - class inherits object for Python2/3 compatibility [16.05.13]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     GPSQL.__init__
        Input:    None.
        Features: Initialize the class
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("GPSQL class initialized")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def flux_to_grid(self, flux_station):
        """
        Name:     GPSQL.flux_to_grid
        Input:    str, station name (flux_station)
        Output:   str, grid station ID (grid_station)
        Features: Returns grid station ID based on the location of a given flux
                  tower
        Depends:  - connectSQL
                  - get_lon_lat
                  - grid_centroid
        """
        # Get lat and lon of flux tower:
        (fst_lon, fst_lat) = self.get_lon_lat(flux_station)

        # Determine grid centroid lon and lat:
        (grd_lon, grd_lat) = grid_centroid(fst_lon, fst_lat)

        # Get grid station name based on centroid coordinates:
        q = ("SELECT met_data.stationid "
             "FROM met_data "
             "WHERE met_data.geom = %s "
             "AND met_data.lon = %s "
             "AND met_data.lat = %s;")
        params = ("grid", grd_lon, grd_lat)

        # Connect to database:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and return results:
            try:
                cur.execute(q, params)
                grid_station = cur.fetchone()[0]
            except:
                self.logger.exception("Query failed.")
                raise IOError(
                    "Failed to fetch grid for station '%s'" % (flux_station))
            else:
                con.close()
                return grid_station
        else:
            self.logger.error("Failed to connect to GePiSaT database")
            raise IOError(
                "Failed to fetch grid for station '%s'" % (flux_station))

    def get_data_point(self, station, var, tp):
        """
        Name:     get_data_points
        Input:    - str, station name (station)
                  - str, variable name (var)
                  - datetime.date, time point (tp)
        Output:   float/numpy.ndarray (my_result)
        Features: Returns data point or array of data for a given station,
                  variable and time
        Depends:  connectSQL
        """
        msvidx = self.get_msvidx(station, var)

        # Define SQL query:
        q = ("SELECT data_set.data "
             "FROM data_set "
             "WHERE data_set.msvidx = %s "
             "AND data_set.datetime = %s;")
        params = (msvidx, tp)

        # Connect to database and start a cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and return result:
            cur.execute(q, params)
            if cur.rowcount == 1:
                my_result = cur.fetchone()[0]
            elif cur.rowcount > 1:
                my_result = numpy.array([])
                for record in cur:
                    my_result = numpy.append(my_result, record[0])
            else:
                my_result = None
                self.logger.warning("No data found!")
        else:
            self.logger.warning("Failed to connect to GePiSaT database")

        return my_result

    def get_date_range(self, station):
        """
        Name:     GPSQL.get_date_range
        Input:    str, station name (station)
        Output:   tuple, starting and ending dates
                  - datetime.date, starting date (sd)
                  - datetime.date, ending date (ed)
        Features: Returns the starting and ending dates of NEE-PPFD data pairs
                  for a given station
        Depends:  - connectSQL
                  - get_msvidx
        """
        # Get msvidx values for specified station:
        ppfdi = self.get_msvidx(station, 'PPFD_f')
        neei = self.get_msvidx(station, 'NEE_f')

        # SQL query parameters:
        params = (ppfdi, neei)

        # Define start date query:
        q1 = ("SELECT data_set.datetime "
              "FROM data_set "
              "WHERE data_set.msvidx = %s OR data_set.msvidx = %s "
              "ORDER BY data_set.datetime ASC LIMIT 1;")

        # Define end date query:
        q2 = ("SELECT data_set.datetime "
              "FROM data_set "
              "WHERE data_set.msvidx = %s "
              "OR data_set.msvidx = %s "
              "ORDER BY data_set.datetime DESC LIMIT 1;")

        # Connect to database and start a cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            try:
                # Get start date from datetime object:
                cur.execute(q1, params)
                sd = cur.fetchone()[0].date()

                # Get end date from datetime object:
                cur.execute(q2, params)
                ed = cur.fetchone()[0].date()
            except:
                self.logger.exception("Failed to fetch date range")
                raise IOError(
                    "Failed to fetch date range for station %s" % (station))
            else:
                # Make the starting date begin at day 1
                sd = sd.replace(day=1)
            finally:
                con.close()
                return (sd, ed)
        else:
            self.logger.error("Failed to connect to GePiSaT database")
            raise IOError(
                "Failed to fetch date range for station %s" % (station))

    def get_elap(self, station):
        """
        Name:     GPSQL.get_elap
        Input:    str, station name (s)
        Output:   - float, elevation, m (z)
                  - float, atmospheric pressure, Pa (patm)
        Features: Returns the elevation and atmospheric pressure for a given
                  station
        Depends:  - connectSQL
                  - flux_to_grid
                  - get_data_point
        Ref:      Allen et al. (1998)
        """
        # Define query w/ parameters:
        q = ("SELECT met_data.ele "
             "FROM met_data "
             "WHERE met_data.stationid = %s;")
        params = (station,)

        # Connect to database:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            try:
                # Execute query and return results:
                cur.execute(q, params)
                station_ele = cur.fetchone()[0]
            except:
                self.logger.exception("Query failed.")
                station_ele = -9999
            else:
                self.logger.debug(
                    "Elevation for %s is %f" % (station, station_ele))
            finally:
                con.close()
        else:
            station_ele = -9999

        # Check to see that elevation is valid:
        if float(station_ele) == -9999:
            # Find CRU Elv:
            elv_sd = datetime.date(2006, 6, 1)
            hdg_station = self.flux_to_grid(station)
            elv_data = self.get_data_point(hdg_station, 'Elv', elv_sd)
            station_ele = elv_data

        # Convert elevation to pressure, Pa:
        z = float(station_ele)
        patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))

        return (z, patm)

    def get_lon_lat(self, station):
        """
        Name:     GPSQL.get_lon_lat
        Input:    str, station name (station)
        Output:   tuple, lon-lat pair
                  - float, longitude (my_lon)
                  - float, latitude (my_lat)
        Features: Return longitude and latitude pair for a given station based
                  on the GePiSaT database meta-data table
        Depends:  connectSQL
        """
        # Define query:
        q = ("SELECT met_data.lon, met_data.lat "
             "FROM met_data "
             "WHERE met_data.stationid = %s;")
        params = (station,)

        # Connect to database and start a cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and return results:
            try:
                cur.execute(q, params)
                my_lon, my_lat = cur.fetchone()
            except:
                self.logger.exception("Querry failed.")
                raise
            else:
                self.logger.debug("Found lon and lat for %s", station)
                con.close()
                return (my_lon, my_lat)
        else:
            self.logger.error("Failed to connect to GePiSaT database!")
            raise IOError("Failed to connect to GePiSaT database!")

    def get_msvidx(self, station, variable):
        """
        Name:     GPSQL.get_msvidx
        Input:    - str, station name (station)
                  - str, variable name (variable)
        Output:   string, msvidx (result)
        Features: Returns the msvidx from the GePiSaT database based on the
                  station and variable name
        Depends:  connectSQL
        """
        # Define query:
        q = ("SELECT var_list.msvidx "
             "FROM var_list "
             "WHERE var_list.stationid = %s "
             "AND var_list.varname = %s;")
        params = (station, variable)

        # Initialize return variable:
        result = ""

        # Connect to database and start cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and fetch results:
            self.logger.debug("Querying database for msvidx ...")
            cur.execute(q, params)
            try:
                result = cur.fetchone()[0]
            except:
                self.logger.exception(
                    ("Could not find an msvidx value for station '%s' and "
                     "variable '%s'") % (station, variable))
                result = ""
            finally:
                con.close()
        else:
            self.logger.warning("Failed to connect to GePiSaT database")

        return result

    def get_stations(self):
        """
        Name:     GPSQL.get_stations
        Input:    None.
        Output:   list, station names (results)
        Features: Returns a list of flux station names from GePiSaT database
        Depends:  connectSQL
        """
        # Define query:
        q = ("SELECT stationid "
             "FROM met_data "
             "WHERE dim=0 "
             "AND geom=%s "
             "ORDER BY stationid ASC;")
        params = ("point",)

        # Initialize return list:
        results = []

        # Connect to database and start cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and fetch results:
            self.logger.debug("Querying database for stations ...")
            cur.execute(q, params)
            for record in cur:
                results.append(record[0])  # <- extract record from tuple

            con.close()
        else:
            self.logger.warning("Failed to connect to GePiSaT database")

        return results


###############################################################################
# FUNCTIONS
###############################################################################


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("data.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # DEBUG
    my_data = DATA()
    my_data.output_dir = "/home/user/Desktop/temp/out"
    my_stations = my_data.find_stations()
    for station in my_data.stations:
        grid_station = my_data.find_station_grid(station)
        print("%s,%s" % (station, grid_station))
    # my_data.create_summary_file("summary_statistics.txt")
    # my_data.write_summary()
