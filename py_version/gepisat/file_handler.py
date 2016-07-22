#!/usr/bin/python
#
# file_handler.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2016-07-22
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
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

from .db_util import connectSQL
from .solar import SOLAR_TOA
from .utilities import add_one_day
from .utilities import add_one_month
from .utilities import elv2pres
from .utilities import grid_centroid


###############################################################################
# CLASSES
###############################################################################
class GPSQL(object):
    """
    Name:     GPSQL
    Features: This class performs data IO with the GePiSaT PSQL database
    History:  Version 3.0.0-dev
              - class inherits object for Python2/3 compatibility [16.05.13]
              - added get monthly NEE:PPFD pairs [16.05.20]
              - created reset properties function [16.05.20]
              - created set working station function [16.05.20]
              - changed logger name [16.05.21]
              - created get daily PPFD function [16.06.26]
              - finished gapfill daily PPFD function [16.06.26]
              - created output directory [16.06.26]
              - moved to file handler module [16.07.22]
              - created get monthly meteorology functions [16.07.22]
              - moved to gepisat package [16.07.22]
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
            self.logger.debug(
                "Using observed radiation, %f MJ m^-2", grid_srad*1e-6)

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

    def get_annual_co2(self, station, start_date):
        """
        Name:     GPSQL.get_annual_co2
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Outputs:  float, array or NoneType
        Features: Returns annual CO2
        Depends:  get_data_point
        """
        annual_sd = start_date.replace(month=1, day=1)
        return self.get_data_point(station, "CO2", annual_sd)

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

    def get_monthly_alpha(self, station, start_date):
        """
        Name:     GPSQL.get_monthly_alpha
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Output:   float or NoneType
        Features: Returns monthly gridded Priestly-Taylor alpha
        Depends:  get_data_point
        """
        return self.get_data_point(station, "alpha", start_date)

    def get_monthly_fapar(self, station, start_date):
        """
        Name:     GPSQL.get_monthly_fapar
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Output:   float or NoneType
        Features: Returns monthly gridded fractionally absorbed
                  photosyntheically active radiation
        Depends:  get_data_point
        """
        return self.get_data_point(station, "FAPAR", start_date)

    def get_monthly_flux_pairs(self, station, start_date):
        """
        Name:     GPSQL.get_monthly_flux_pairs
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

    def get_monthly_tair(self, station, start_date):
        """
        Name:     GPSQL.get_monthly_tair
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Output:   float or NoneType
        Features: Returns monthly gridded near-surface air temperature
        Depends:  get_data_point
        """
        return self.get_data_point(station, "Tc", start_date)

    def get_monthly_vpd(self, station, start_date):
        """
        Name:     GPSQL.get_monthly_vpd
        Input:    - str, station name (station)
                  - datetime.date (start_date)
        Output:   float or NoneType
        Features: Returns monthly gridded vapor pressure deficit
        Depends:  get_data_point
        """
        return self.get_data_point(station, "VPD", start_date)

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
                self.logger.error("Could not find an msvidx value "
                                  "for station '%s' "
                                  "and variable '%s'" % (station, variable))
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

        # Database indexes for station-specific:
        self._ppfd_msv = None      # photosynthetic photon flux density
        self._nee_msv = None       # net ecosystem exchange

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
                  - get_msvidx
        """
        self.logger.debug("Resetting working station properties ...")
        self.reset_properties()

        self.logger.debug("Setting up station %s ...", station)
        self._station = station

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

        # Set PPFD and NEE msvidx values:
        self._ppfd_msv = self.get_msvidx(station, "PPFD_f")
        self._nee_msv = self.get_msvidx(station, "NEE_f")

        self.logger.debug("Found PPFD msvidx %s", self._ppfd_msv)
        self.logger.debug("Found NEE msvidx %s", self._nee_msv)

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
        self.logger.debug("Found elevation, %f m", self._elv)
        self.logger.debug("Found pressure, %f Pa", self._atmpres)
