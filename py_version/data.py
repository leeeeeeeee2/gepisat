#!/usr/bin/python
#
# model.py
#
# VERSION 2.2.0-dev
#
# LAST UPDATED: 2016-07-06
#
# ---------
# citation:
# ---------
# I. C. Prentice, T. W. Davis, X. M. P. Gilbert, B. D. Stocker, B. J. Evans,
# H. Wang, and T. F. Keenan, "The Global ecosystem in Space and Time (GePiSaT)
# Model of the Terrestrial Biosphere," (in progress).

# @TODO: gapfilling methods

###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import logging
import os

import numpy

from db_setup import connectSQL
from solar import SOLAR_TOA
from utilities import add_one_day
from utilities import add_one_month
from utilities import elv2pres
from utilities import init_summary_dict
from utilities import grid_centroid
from utilities import simpson


###############################################################################
# CLASSES
###############################################################################
class DATA(object):
    """
    Name:     DATA
    Features: This class performs the data IO
    History   Version 2.2.0-dev
              - class inherits object for Python2/3 compatibility [16.05.13]
              - created end and start date properties [16.05.20]
              - created set working station function [16.05.20]
              - created find elv, pressure, & monthly NEE:PPFD pairs [16.05.20]
              - added properties for fh attributes [16.05.21]
              - finished gapfill monthly PPFD function [16.06.26]
              - output directory setter moved to function [16.06.26]
              - updated logging statements [16.07.06]
              - created sub to daily gpp function [16.07.06]
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
    def station(self):
        """Current working station"""
        if self.fh._station is not None:
            return self.fh._station
        else:
            self.logger.error("Trying to access station before setting!")
            raise RuntimeError("Station is not set!")

    @station.setter
    def station(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def lon(self):
        """Current working station's longitude"""
        if self.fh._station is not None and self.fh._lon is not None:
            return self.fh._lon
        else:
            self.logger.error("Station or its lon is not set!")
            raise RuntimeError("Station and/or station lon is not set!")

    @lon.setter
    def lon(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def lat(self):
        """Current working station's latitide"""
        if self.fh._station is not None and self.fh._lat is not None:
            return self.fh._lat
        else:
            self.logger.error("Station or its lat is not set!")
            raise RuntimeError("Station and/or station lat is not set!")

    @lat.setter
    def lat(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def elv(self):
        """Current working station's elevation"""
        if self.fh._station is not None and self.fh._elv is not None:
            return self.fh._elv
        else:
            self.logger.error("Station or its elevation is not set!")
            raise RuntimeError("Station and/or station elv is not set!")

    @elv.setter
    def elv(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def atmpres(self):
        """Current working station's atm. pressure"""
        if self.fh._station is not None and self.fh._atmpres is not None:
            return self.fh._atmpres
        else:
            self.logger.error("Station or its atm. pressure is not set!")
            raise RuntimeError("Station and/or station atmpres is not set!")

    @atmpres.setter
    def atmpres(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def end_date(self):
        """Current working station's end date"""
        if self.fh._enddate is not None:
            return self.fh._enddate
        else:
            print("No end date found; try setting working station.")
            self.logger.warning("End date not set; returning epoch")
            return datetime.date(2000, 1, 1)

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
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def start_date(self):
        """Current station's starting date"""
        if self.fh._startdate is not None:
            return self.fh._startdate
        else:
            print("No start date found; please set working station.")
            self.logger.warning("Start date not set; returning epoch")
            return datetime.date(2000, 1, 1)

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
        Features: Convenience function that searches for starting and ending
                  dates for a given station's available data
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
        Name:     DATA.find_elv_pressure
        Inputs:   str, station name (station)
        Outputs:  tuple of floats
                  * float, elevation, m
                  * float, atm. pressure, Pa
        Features: Convenience function for getting a station's elevation and
                  pressure.
        """
        return self.fh.get_elap(station)

    def find_monthly_nee_ppfd(self, start_date):
        """
        Name:     DATA.find_monthly_nee_ppfd
        Inputs:   datetime.date, starting date of the search month (start_date)
        Outputs:  tuple of arrays
                  * numpy.ndarray, NEE
                  * numpy.ndarray, PPFD
        Features: Returns tuple of data arrays corresponding to the current
                  flux station's NEE and PPFD pairs for a given month
        """
        if self.fh._station is not None:
            return self.fh.get_monthly_flux_pairs(self.fh._station, start_date)
        else:
            self.logger.error("No working station set.")
            raise IOError("No flux station found! Please set working station.")

    def find_station_grid(self, station):
        """
        Name:     DATA.find_station_grid
        Inputs:   str, station name (station)
        Outputs:  str, grid name
        Features: Convenience function that returns the grid name for a given
                  station
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

    def gapfill_monthly_ppfd(self, sd, to_save=False):
        """
        Name:     DATA.gapfill_monthly_ppfd
        Inputs:   - datetime.date, date (sd)
                  - bool, write to file flag (to_save)
        Outputs:  tuple
                  - numpy.ndarray, timestamps
                  - numpy.ndarray, PPFD values
        Features: Returns half-hourly gapfilled photosynthetic photon flux
                  density (PPFD) and associated timestamps for a given month
        Depends:  - add_one_month
                  - add_one_day
        """
        # Initialize monthly gapless PPFD array & associated datetimes:
        gf_ppfd = numpy.array([])
        gf_dates = numpy.array([])

        # Calculate the end date:
        end_date = add_one_month(sd)

        # Iterate through each day of the month:
        cur_date = sd
        while cur_date < end_date:
            # Gapfill daily PPFD
            (daily_time, daily_ppfd) = self.fh.gapfill_daily_ppfd(cur_date,
                                                                  to_save)

            # Append data to monthly arrays:
            gf_dates = numpy.append(gf_dates, daily_time)
            gf_ppfd = numpy.append(gf_ppfd, daily_ppfd)

            # Increment day
            cur_date = add_one_day(cur_date)

        return (gf_dates, gf_ppfd)

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

    def print_current_vals(self):
        """Convenience function for printing data properties"""
        print("Station: %s" % (self.station))
        print("  lon:   %0.3f deg" % (self.lon))
        print("  lat:   %0.3f deg" % (self.lat))
        print("  elv:   %0.3f m" % (self.elv))
        print("  P:     %0.3f kPa" % (1e-3*self.atmpres))
        print("  start: %s" % (self.start_date))
        print("  end:   %s" % (self.end_date))
        print("  grid:  %s" % (self.fh._grid))

    def set_output_directory(self, my_dir):
        """
        Name:     DATA.set_output_directory
        Inputs:   str, output directory (my_dir)
        Outputs:  None.
        Features: Assigns the output directory
        """
        if isinstance(my_dir, str):
            if os.path.isdir(my_dir):
                self.logger.debug("Output directory set to %s", my_dir)
                self.fh.set_working_directory(my_dir)
                self._outputdir = my_dir
            else:
                self.logger.error("Directory %s does not exist!", my_dir)
                raise IOError("Output directory does not exist!")
        else:
            try:
                temp_dir = my_dir.decode("utf-8")
            except:
                self.logger.exception("Failed to read output directory.")
                raise
            else:
                if os.path.isdir(temp_dir):
                    self.logger.debug("Output directory set to %s", temp_dir)
                    self.fh.set_working_directory(my_dir)
                    self._outputdir = temp_dir
                else:
                    self.logger.error("Directory %s does not exist!", temp_dir)
                    raise IOError("Output directory does not exist!")

    def set_working_station(self, station):
        """
        Name:     DATA.set_working_station
        Inputs:   str, flux station (station)
        Outputs:  None.
        Features: Prepares station parameters
        """
        if station in self.stations:
            self.logger.info("Setting up for station %s", station)
            # Set station in file handle:
            self.fh.set_working_station(station)
            # self._workingstation = station
            # self._startdate, self._enddate = self.find_date_range(station)
            # self._grid = self.find_station_grid(station)
            # self._elv, self._atmpres = self.find_elv_pressure(station)

    def sub_to_daily_gpp(self, ts, gpp, gpp_err, to_save=False):
        """
        Name:     calc_daily_gpp
        Inputs:   - numpy.ndarray, half-hourly datetimes for a month (ts)
                  - numpy.ndarray, half-hourly gap-filled GPP, umol/m2/s (gpp)
                  - numpy.ndarray, half-hourly GPP model errors (gpp_err)
                  - [optional] bool, whether to save to file (to_save)
        Output:   numpy.ndarray, multidimensional array
                  > 'Timestamp' daily timestamps
                  > 'GPP' daily GPP, mol/m2
                  > 'GPP_err' daily GPP error, mol/m2
        Features: Returns time stamped daily GPP and its associated errors
        Depends:  - add_one_day
                  - simpson
        """
        # Get starting and ending dates:
        starting_date = ts[0]
        ending_date = ts[-1]

        # Assign output file and write header line if requested:
        out_file = "%s_%s-%s_daily_GPP.txt" % (
            self.station, starting_date, ending_date)
        out_path = os.path.join(self.outputdir, out_file)
        if to_save:
            header_line = "Timestamp,GPP_mol.m2,GPP_err_mol.m2\n"
            try:
                f = open(out_path, "w")
            except:
                self.logger.exception("Cannot write to file '%s'", out_path)
            else:
                f.write(header_line)
                f.close()

        # Iterate through days:
        cur_date = starting_date
        while (cur_date < ending_date):
            # Find GPP values associated with current day:
            my_idx = numpy.where(
                (ts >= cur_date) & (ts < add_one_day(cur_date)))
            my_gpp = gpp[my_idx]
            my_gpp_err = gpp_err[my_idx]

            # Integrate to daily:
            sec_in_hh = 1800   # seconds in a half-hour (i.e., simpson's h)
            day_gpp = simpson(my_gpp, sec_in_hh)
            day_gpp_err = simpson(my_gpp_err, sec_in_hh)
            # @TODO: should this be clipped to a min of zero?

            # Convert from umol to moles:
            day_gpp *= 1e-6
            day_gpp_err *= 1e-6

            # Save results:
            if to_save:
                try:
                    f = open(out_path, 'a')
                except IOError:
                    self.logger.exception(
                        "Failed to appending to '%s'", out_path)
                else:
                    f.write("%s,%f,%f\n" % (cur_date, day_gpp, day_gpp_err))
                finally:
                    f.close()

            if cur_date == starting_date:
                gpp_daily = numpy.array(
                    (cur_date.date(), day_gpp, day_gpp_err),
                    dtype={'names': ('Timestamp', 'GPP', 'GPP_err'),
                           'formats': ('O', 'f4', 'f4')},
                    ndmin=1
                )
            else:
                temp_array = numpy.array(
                    (cur_date.date(), day_gpp, day_gpp_err),
                    dtype={'names': ('Timestamp', 'GPP', 'GPP_err'),
                           'formats': ('O', 'f4', 'f4')},
                    ndmin=1
                )
                gpp_daily = numpy.append(gpp_daily, temp_array, axis=0)

            # Increment current day:
            cur_date = add_one_day(cur_date)

        return gpp_daily

    def write_summary(self):
        """
        Name:     DATA.write_summary
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
              - added get monthly NEE:PPFD pairs [16.05.20]
              - created reset properties function [16.05.20]
              - created set working station function [16.05.20]
              - changed logger name [16.05.21]
              - created get daily PPFD function [16.06.26]
              - finished gapfill daily PPFD function [16.06.26]
              - created output directory [16.06.26]
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
        self.logger = logging.getLogger("GPSQL")
        self.logger.info("GPSQL class initialized")
        self.reset_properties()

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
        self._lon = fst_lon
        self._lat = fst_lat

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

    def gapfill_daily_ppfd(self, cur_date, to_save=False):
        """
        Name:     GPSQL.gapfill_daily_ppfd
        Inputs:   - datetime.date, date (cur_date)
                  - [OPTIONAL] bool, write results to file flag (to_save)
        Outputs:  tuple
                  - numpy.ndarray, timestamps
                  - numpy.ndarray, gapfilled PPFD values,
        Features: Returns half-hourly gapfilled photosynthetic photon flux
                  density (PPFD) and associated timestamps; write results to
                  file if requested
        """
        self.logger.debug("Gapfilling %s", cur_date)

        # Assign output file and write header line if requested:
        out_file = "%s_%s_gapfill.txt" % (self._station, cur_date)
        out_path = os.path.join(self._outputdir, out_file)
        if to_save:
            header_line = "Timestamp,ObsPPFD_umol.m2,GFPPFD_umol.m2\n"
            try:
                f = open(out_path, "w")
            except:
                self.logger.exception("Cannot write to file '%s'", out_path)
            else:
                f.write(header_line)
                f.close()

        # Initialize gapless PPFD array & associated datetimes:
        gf_ppfd = numpy.array([])
        gf_dates = numpy.array([])

        # Get daily PPFD values (in a numpy.array) from database:
        (daily_ts, daily_ppfd) = self.get_daily_ppfd(cur_date)

        # Check to see if any gaps are present in current day's observations:
        number_obs = len(daily_ppfd)
        if number_obs < 49:
            # Gaps are present:
            self.logger.debug("Gaps found.")
            gapfill_dict = {}

            # Create start and end datetime objects for this day:
            start_time = datetime.datetime(
                cur_date.year, cur_date.month, cur_date.day, 0, 0, 0)

            end_time = datetime.datetime(
                cur_date.year, cur_date.month, cur_date.day, 23, 59, 59)

            # Initialize dictionary values with half-hourly timestamp keys:
            cur_time = start_time
            while cur_time < end_time:
                my_time = "%s" % (cur_time.time())
                gapfill_dict[my_time] = 0.0

                # Add datetime object to array of monthly timestamps:
                gf_dates = numpy.append(gf_dates, [cur_time, ])

                cur_time += datetime.timedelta(minutes=30)

            # Convert date to Julian day:
            jday = cur_date.timetuple().tm_yday
            self.logger.debug("DOY = %d", jday)

            # Calculate daily ET solar radiation curve:
            et_solar = SOLAR_TOA(self._lat, self._lon)
            et_solar.calculate_daily_fluxes(jday, cur_date.year)

            # Get satellite measurement of solar rad (SWdown) [W m-2]:
            grid_srad = self.get_data_point(self._grid, 'SWdown', cur_date)

            # Convert to daily shortwave radiation [J m-2]:
            # NOTE: if None, then gridded data is not available, so set it
            #       equal to modeled
            if grid_srad is not None:
                grid_srad *= (86400.0)
            else:
                self.logger.warning("Grid %s shortwave not found!", self._grid)
                grid_srad = et_solar.ho_jm2
            self.logger.debug("Using observed radiation, %f J m^-2", grid_srad)

            # Calculate scaling factor (i.e., observed/modeled):
            if et_solar.ho_jm2 != 0:
                sfactor = grid_srad/et_solar.ho_jm2
            else:
                self.logger.warning("Zero solar irradiation at %s!", cur_date)
                sfactor = 1.0
            self.logger.debug("Scaling factor set to %f", sfactor)

            # Add scaled half-hourly PPFD to dictionary [umol m-2 s-1]:
            for i in range(48):
                val = et_solar.ppfd_hh[i]
                my_time = "%s" % et_solar.local_time[i].time()
                gapfill_dict[my_time] = sfactor*val

            # Add observations to gapfill dictionary & save for output:
            ppfd_obs = {}
            for i in range(number_obs):
                my_time = "%s" % daily_ts[i].time()
                gapfill_dict[my_time] = daily_ppfd[i]
                ppfd_obs[my_time] = daily_ppfd[i]

            # Save to PPFD time series:
            gf_ppfd = numpy.append(
                gf_ppfd,
                [gapfill_dict[x] for x in sorted(list(gapfill_dict.keys()))])

            if to_save:
                for t in sorted(list(gapfill_dict.keys())):
                    if t in ppfd_obs.keys():
                        obs = ppfd_obs[t]
                    else:
                        obs = numpy.nan
                    gfv = gapfill_dict[t]
                    dt = "%s %s" % (cur_date, t)
                    try:
                        f = open(out_path, 'a')
                    except IOError:
                        self.logger.exception(
                            "Failed to appending to '%s'", out_path)
                    else:
                        f.write("%s,%f,%f\n" % (dt, obs, gfv))
                        f.close()
        else:
            self.logger.debug("No gaps found.")

            # Append daily series
            #   NOTE: drop last entry from daily_ppfd (midnight next day)
            gf_ppfd = daily_ppfd[0:-1]
            gf_dates = daily_ts[0:-1]

            if to_save:
                for i in range(len(daily_ts) - 1):
                    dt = "%s" % daily_ts[i]
                    obs = daily_ppfd[i]
                    try:
                        f = open(out_path, 'a')
                    except IOError:
                        self.logger.exception(
                            "Failed to appending to '%s'", out_path)
                    else:
                        f.write("%s,%0.3f,%0.3f\n" % (dt, obs, obs))
                        f.close()

        return (gf_dates, gf_ppfd)

    def get_daily_ppfd(self, start_date):
        """
        Name:     GPSQL.get_daily_ppfd
        Inputs:   datetime.date, date (start_date)
        Outputs:  tuple
                  - numpy.ndarray, time stamps
                  - numpy.ndarray, PPFD values
        Features: Returns daily photosynthetic photon flux density (PPFD) and
                  associated timestamps arrays for a given date
        Depends:  connectSQL
        """
        # Define query:
        q = ("SELECT data_set.datetime, data_set.data "
             "FROM data_set "
             "WHERE data_set.msvidx = %s "
             "AND data_set.datetime BETWEEN DATE %s AND DATE %s "
             "ORDER BY data_set.datetime ASC;")
        params = (self._ppfd_msv, start_date, add_one_day(start_date))

        # Initialize return arrays:
        ppfd_vals = numpy.array([])
        time_vals = numpy.array([])

        # Connect to database and start a cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            # Execute query and store results:
            cur.execute(q, params)
            if cur.rowcount > 0:
                for record in cur:
                    time_vals = numpy.append(time_vals, record[0])
                    ppfd_vals = numpy.append(ppfd_vals, record[1])
            else:
                self.logger.warning("No data found!")

            # Close connection and return results:
            con.close()
        else:
            self.logger.warning("Failed to connect to GePiSaT database")

        return (time_vals, ppfd_vals)

    def get_data_point(self, station, var, tp):
        """
        Name:     GPSQL.get_data_points
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
            my_result = None
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
        if station == self._station:
            ppfdi = self._ppfd_msv
            neei = self._nee_msv
        else:
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
                  - elv2pres
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
            self.logger.debug("Trying grid elevation for station %s", station)
            elv_sd = datetime.date(2006, 6, 1)
            hdg_station = self.flux_to_grid(station)
            elv_data = self.get_data_point(hdg_station, 'Elv', elv_sd)
            if elv_data is None:
                station_ele = -9999
            else:
                station_ele = elv_data

        # Convert elevation to pressure, Pa:
        z = float(station_ele)
        if z == -9999:
            patm = -9999
        else:
            patm = elv2pres(z)

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

    def get_monthly_flux_pairs(self, station, start_date):
        """
        Name:     get_monthly_flux_pairs
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Output:   tuple, PPFD and NEE observations
                  - numpy.ndarray, NEE (nee_vals)
                  - numpy.ndarray, PPFD (ppfd_vals)
        Features: Returns one month of PPFD-NEE observation pairs for a given
                  station and month
        Depends:  - add_one_month
                  - connectSQL
                  - get_msvidx
        """
        if station == self._station:
            ppfd_idx = self._ppfd_msv
            nee_idx = self._nee_msv
        else:
            # Get msvidx values for specified station
            ppfd_idx = self.get_msvidx(station, 'PPFD_f')
            nee_idx = self.get_msvidx(station, 'NEE_f')

        # Increment start date one month:
        end_date = add_one_month(start_date)

        # SQL query parameters:
        # NOTE: nee msvidx comes before ppfd msvidx (i.e., 01 v. 15)
        params = (nee_idx, ppfd_idx, start_date, end_date, nee_idx, ppfd_idx)

        # Define query
        # NOTE: okay to use string concatenation because ppfd and nee are
        # function return values:
        q = (
            "SELECT * "
            "FROM crosstab('"
            "select datetime, msvidx, data from data_set "
            "where msvidx = ''%s'' "
            "or msvidx = ''%s'' "
            "and datetime between date ''%s'' and date ''%s'' "
            "order by 1, 2', "
            "'select distinct msvidx from data_set "
            "where msvidx = ''%s'' "
            "or msvidx = ''%s'' order by 1') "
            "AS ct(row_name TIMESTAMP, nee FLOAT, ppfd FLOAT);"
            ) % params

        # Initialize empty arrays:
        nee_vals = numpy.array([])
        ppfd_vals = numpy.array([])

        # Connect to database and start a cursor:
        self.logger.debug("Connecting to GePiSaT database ...")
        con = connectSQL()
        if con is not None:
            cur = con.cursor()

            try:
                # Execute query and fetch results:
                cur.execute(q)
            except:
                self.logger.exception("Search failed.")
                raise
            else:
                if cur.rowcount > 0:
                    for record in cur:
                        # NOTE: nee is ordered before ppfd because nee's
                        #       msvidx value (i.e., *.01) comes before ppfd's
                        #       msvidx value (i.e., *.15)
                        _, nee, ppfd = record

                        # Only save matched pairs:
                        if ppfd is not None and nee is not None:
                            nee_vals = numpy.append(nee_vals, nee)
                            ppfd_vals = numpy.append(ppfd_vals, ppfd)
            finally:
                con.close()
        else:
            self.logger.warning("Failed to connect to GePiSaT database")

        return (nee_vals, ppfd_vals)

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

    def reset_properties(self):
        """
        Name:     GPSQL.reset_properties
        Inputs:   None.
        Outputs:  None.
        Features: Resets the station properties
        """
        self._station = None       # flux station
        self._grid = None          # grid station
        self._lon = None           # longitude
        self._lat = None           # latitide
        self._elv = None           # elevation
        self._atmpres = None       # atmospheric pressure
        self._startdate = None     # station's data starting date
        self._enddate = None       # station's data ending date
        self._ppfd_msv = None      # database index for station PPFD
        self._nee_msv = None       # database index for station NEE

    def set_working_directory(self, my_dir):
        """
        Name:     GPSQL.set_working_directory
        Inputs:   str, existing directory (my_dir)
        Outputs:  None.
        Features: Assigns output/working directory
        """
        if os.path.isdir(my_dir):
            self._outputdir = my_dir
        else:
            raise IOError("Directory '%s' does not exist!", my_dir)

    def set_working_station(self, station):
        """
        Name:     GPSQL.set_working_station
        Inputs:   str, flux station (station)
        Outputs:  None.
        Features: Sets the appropriate class variables for a given station
        Depends:  - flux_to_grid
                  - get_date_range
                  - get_elap
        """
        self.logger.debug("Resetting working station properties ...")
        self.reset_properties()

        self.logger.debug("Setting up station %s ...", station)
        self._station = station

        # Set PPFD and NEE msvidx values:
        self._ppfd_msv = self.get_msvidx(station, "PPFD_f")
        self._nee_msv = self.get_msvidx(station, "NEE_f")
        self.logger.debug("Found PPFD msvidx %s", self._ppfd_msv)
        self.logger.debug("Found NEE msvidx %s", self._nee_msv)

        # Search for equivalent grid station (also sets lon and lat):
        try:
            self._grid = self.flux_to_grid(station)
        except IOError:
            self.logger.error(("Search failed. Check to make certain that "
                               "your credentials for the database are "
                               "correct."))
            self._grid = None
        except:
            self.logger.exception("Encountered unknown error.")
            raise
        else:
            self.logger.debug("Found grid %s", self._grid)

        # Search for starting and ending dates:
        try:
            self._startdate, self._enddate = self.get_date_range(station)
        except IOError:
            self.logger.error(("Search failed. Check to make certain that "
                               "your credentials for the database are "
                               "correct."))
            self._startdate = None
            self._enddate = None
        except:
            self.logger.exception("Encountered unknown error.")
            raise
        else:
            self.logger.debug("Found date range %s to %s" % (
                self._startdate, self._enddate))

        # Set elevation and pressure:
        self._elv, self._atmpres = self.get_elap(station)
        self.logger.debug("Found elevation, %f m, and pressure, %f Pa" % (
            self._elv, self._atmpres))


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
    output_dir = os.path.join(
        os.path.expanduser("~"), "Desktop", "temp", "out")
    my_data.set_output_directory(output_dir)
    my_stations = my_data.find_stations()
    station = my_data.stations[1]
    my_data.set_working_station(station)
    sd = my_data.start_date
    gf_time, gf_ppfd = my_data.gapfill_monthly_ppfd(sd, True)
    # nee, ppfd = my_data.find_monthly_nee_ppfd(sd)
    # for i in range(len(nee)):
    #    print("NEE: %f; PPFD: %f" % (nee[i], ppfd[i]))

    # for station in my_data.stations:
    #    my_data.set_working_station(station)
    #    my_data.print_current_vals()
    # my_data.create_summary_file("summary_statistics.txt")
    # my_data.write_summary()
