#!/usr/bin/python
#
# data.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-13
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

import numpy

from .file_handler import GPSQL
from .lue import LUE
from .resources import mkdir_p
from .utilities import add_one_day
from .utilities import add_one_month
from .utilities import init_summary_dict
from .utilities import simpson


###############################################################################
# CLASSES
###############################################################################
class DATA(object):
    """
    Name:     DATA
    Features: This class performs the data IO
    History   Version 3.0.0-dev
              - class inherits object for Python2/3 compatibility [16.05.13]
              - created end and start date properties [16.05.20]
              - created set working station function [16.05.20]
              - created find elv, pressure, & monthly NEE:PPFD pairs [16.05.20]
              - added properties for fh attributes [16.05.21]
              - finished gapfill monthly PPFD function [16.06.26]
              - output directory setter moved to function [16.06.26]
              - updated logging statements [16.07.06]
              - created sub to daily gpp function [16.07.06]
              - created find monthly meteorology function [16.07.22]
              - import LUE [16.07.22]
              - moved to gepisat package [16.07.22]
              - added warning statement for negative daily GPP [17.01.25]
              - changed datetime to date in daily GPP writeouts [17.01.27]
              - sub to daily function writes to a single file [17.01.27]
              - changed alpha to station variable [17.05.13]
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
        self.logger.debug("DATA class initialized")

        # Hidden properties:
        self._outputdir = "."
        self._stations = {}
        self._summaryfile = "summary.txt"

        # Data file handler (change if necessary)
        self.logger.debug("Initializing new data file handler")
        self.fh = GPSQL()

        # Summary value dictionary:
        self.logger.debug("Initializing summary value dictionary")
        self.summary = init_summary_dict()

        # LUE monthly data handler
        self.lue = LUE()

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

    def find_monthly_meteorology(self, start_date):
        """
        Name:     DATA.find_monthly_meteorology
        Inputs:   datetime.date, date (start_date)
        Outputs:  dict
        Features: Returns monthly gridded meteorology values
        Note:     Priestly-Taylor alpha is currently a flux-station variable
        """
        # Retrieve monthly gridded meteorology
        if self.fh._grid is not None:
            cpa = self.fh.get_monthly_alpha(self.fh._station, start_date)
            co2 = self.fh.get_annual_co2(self.fh._grid, start_date)
            fpar = self.fh.get_monthly_fapar(self.fh._grid, start_date)
            tc = self.fh.get_monthly_tair(self.fh._grid, start_date)
            vpd = self.fh.get_monthly_vpd(self.fh._grid, start_date)

            return {"alpha": cpa, "co2": co2, "fpar": fpar, "tair": tc,
                    "vpd": vpd}
        else:
            self.logger.error("No working station set!")
            raise IOError("No flux station found! Please set working station.")

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
            self.logger.debug("Setting up for station %s", station)

            # Define the light-use efficiency data output file:
            self.lue_file = "{}_LUE.txt".format(station)

            # Set station in file handle:
            self.fh.set_working_station(station)

            # Alternative setters:
            # self._workingstation = station
            # self._startdate, self._enddate = self.find_date_range(station)
            # self._grid = self.find_station_grid(station)
            # self._elv, self._atmpres = self.find_elv_pressure(station)

    def sub_to_daily_gpp(self, ts, gpp, gpp_err, to_save=False):
        """
        Name:     sub_to_daily_gpp
        Inputs:   - numpy.ndarray, half-hourly datetimes for a month (ts)
                  - numpy.ndarray, half-hourly gap-filled GPP, umol/m2/s (gpp)
                  - numpy.ndarray, half-hourly GPP model errors (gpp_err)
                  - [optional] bool, whether to save to file (to_save)
        Output:   numpy.ndarray, multidimensional array
                  > 'Timestamp' daily timestamps
                  > 'GPP' daily GPP, mol/m2
                  > 'GPP_err' daily GPP error, mol/m2
        Features: Returns time stamped daily GPP and its associated errors;
                  writes daily results to individual monthly files into
                  site-specific subdirectories (created if not existing)
        Depends:  - add_one_day
                  - mkdir_p
                  - simpson
        """
        # Get starting and ending dates
        # (gap filling starts the time series at midnight):
        starting_date = ts[0]
        ending_date = ts[-1]

        # Assign output file and write header line if requested:
        out_file = "%s_daily_GPP.txt" % (self.station)
        out_dir = os.path.join(self.output_dir, "daily_gpp")
        out_path = os.path.join(out_dir, out_file)
        if to_save:
            if not os.path.isdir(out_dir):
                try:
                    mkdir_p(out_dir)
                except:
                    self.logger.critical(
                        "Could not create output directory for daily GPP!")
                    to_save = False

            if not os.path.isfile(out_path) and to_save:
                header_line = "Timestamp,GPP_mol.m2,GPP_err_mol.m2\n"
                try:
                    f = open(out_path, "w")
                except:
                    self.logger.error("Cannot write to file '%s'", out_path)
                    to_save = False
                else:
                    f.write(header_line)
                    f.close()

        # Iterate through days:
        cur_date = starting_date
        while (cur_date < ending_date):
            # Find GPP values associated with current day:
            # @TODO: handle if numpy where returns empty array
            my_idx = numpy.where(
                (ts >= cur_date) & (ts < add_one_day(cur_date)))
            my_gpp = gpp[my_idx]
            my_gpp_err = gpp_err[my_idx]

            if (my_gpp < 0).any():
                self.logger.warning("Negative daily GPP for %s at %s" % (
                    self.station, cur_date))

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
                    f.write(
                        "%s,%f,%f\n" % (cur_date.date(), day_gpp, day_gpp_err))
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

    def update_lue(self, month, gpp, gpp_err, ppfd):
        """
        Name:     DATA.update_lue
        Inputs:   - datetime.date (month)
                  - numpy.ndarray, monthly gross-primary production (gpp)
                  - numpy.ndarray, monthly modeling errors (gpp_err)
                  - numpy.ndarray, monthly photosyn. photon flux density (ppfd)
        Features: Updates the monthly meteorology for the light-use efficiency
                  model
        """
        # STAGE 2: The New LUE Model:
        #  GPP = f(PPFD, fAPAR, VPD, CPA, Tair, Patm, CO2)
        #  - monthly variables: PPFD, fAPAR, CPA, VPD, Tair
        #  - annual variables: CO2
        #  - constants in time: Patm=f(elevation)
        try:
            monthly_met = self.find_monthly_meteorology(month)
        except IOError:
            self.logger.error("Cannot update LUE without a set station!")
            raise
        else:
            self.lue.add_station_val(self.fh._station,
                                     month,
                                     gpp,
                                     gpp_err,
                                     monthly_met['fpar'],
                                     ppfd,
                                     monthly_met['vpd'],
                                     monthly_met['alpha'],
                                     monthly_met['tair'],
                                     monthly_met['co2'],
                                     self.atmpres,
                                     self.elv)

    def write_monthly_results(self, station):
        """
        Name:     DATA.write_monthly_results
        Inputs:   str, station name (station)
        Outputs:  None.
        Features: Writes monthly light-use efficiency model values to file
        """
        self.lue.write_out_val(station, self.lue_file)

    def write_summary(self):
        """
        Name:     DATA.write_summary
        Inputs:   None.
        Outputs:  None.
        Features: Writes summary values to file

        @TODO:    Currently, DATA's summary file is not updated by PARTI_STATS
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
    my_data.print_current_vals()
    sd = my_data.start_date
    monthly_met = my_data.find_monthly_meteorology(add_one_month(sd))
    for k, v in monthly_met.items():
        print("  %s = %s" % (k, v))
    # gf_time, gf_ppfd = my_data.gapfill_monthly_ppfd(sd, to_save=False)

    # nee, ppfd = my_data.find_monthly_nee_ppfd(sd)
    # for i in range(len(nee)):
    #    print("NEE: %f; PPFD: %f" % (nee[i], ppfd[i]))

    # for station in my_data.stations:
    #    my_data.set_working_station(station)
    #    my_data.print_current_vals()
    # my_data.create_summary_file("summary_statistics.txt")
    # my_data.write_summary()
