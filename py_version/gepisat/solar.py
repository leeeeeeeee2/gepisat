#!/usr/bin/python
#
# solar.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2016-07-22
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

import numpy

from .const import ke
from .const import keps
from .const import kfFEC
from .const import kGsc
from .const import komega
from .const import pir
from .utilities import dcos
from .utilities import dsin


###############################################################################
# CLASSES
###############################################################################
class SOLAR_TOA:
    """
    Name:     SOLAR_TOA
    Features: This class calculates the half-hourly extraterrestrial PPFD
              [umol m-2 s-1], and the daily solar irradiation Ho [J m-2]
              based on the SPLASH model
    History:  Version 3.0.0-dev
              - import global constants [16.04.01]
              - import utility functions [16.04.01]
              - updated methods to mirror SPLASH [16.04.01]
              - fixed xrange for Python2/3 support [16.06.26]
              - fixed bad logging statement [16.06.26]
              - moved to gepisat package [16.07.22]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lat, lon):
        """
        Name:     SOLAR_TOA.__init__
        Inputs:   - float, longitude (lon)
                  - float, latitude (lat)
        Features: Initializes point-based top-of-the-atmosphere radiation class
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("SOLAR_TOA class called")

        # Error handle and assign required public variables:
        if lat > 90.0 or lat < -90.0:
            self.logger.error(
                "Latitude outside range of validity, (-90 to 90)!")
            raise ValueError(
                "Latitude outside range of validity, (-90 to 90)!")
        else:
            self.logger.info("latitude set to %0.3f degrees", lat)
            self.lat = lat

        if lon > 180.0 or lon < -180.0:
            self.logger.error(
                "Longitude outside range of validity (-180 to 180)!")
            raise ValueError(
                "Longitude outside range of validity (-180 to 180)!")
        else:
            self.lon = lon

    def calculate_daily_fluxes(self, n, y=0):
        """
        Name:     SOLAR_TOA.calculate_daily_fluxes
        Inputs:   - int, day of year (n)
                  - [optional] int, year (y)
        Outputs:  None.
        Features: Calculates the half-hourly and daily top-of-the-atmosphere
                  solar radiation and photosynthetic photon fluxes
        """
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            kN = 365
            self.year = 2001
        elif y < 0:
            self.logger.error("year set out of range")
            raise ValueError(
                "Please use a valid Julian or Gregorian calendar year")
        else:
            kN = self.julian_day((y+1), 1, 1) - self.julian_day(y, 1, 1)
            self.year = y
        self.kN = kN

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Create local time series, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if n < 1 or n > 366:
            self.logger.error(
                "Day of year outside range of validity, (1 to 366)!")
            raise ValueError(
                "Day of year outside range of validity (1 to 366)!")
        else:
            self.day = n
            local_hh = numpy.array([0.5*i for i in range(48)])
            self.local_time = numpy.array([
                datetime.datetime(self.year, 1, 1, 0, 0, 0) +
                datetime.timedelta(days=(n-1)) +
                datetime.timedelta(hours=i) for i in local_hh
            ])

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        self.my_nu = my_nu
        self.my_lambda = my_lambda
        self.logger.debug("true anomaly, nu, set to %f degrees", my_nu)
        self.logger.debug("true lon, lambda, set to %f degrees", my_lambda)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        kee = ke**2
        my_rho = (1.0 - kee)/(1.0 + ke*dcos(my_nu))
        dr = (1.0/my_rho)**2
        self.dr = dr
        self.logger.debug(
            "relative Earth-Sun distance, rho, set to %f", my_rho)
        self.logger.debug("distance factor, dr, set to %f", dr)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        delta = numpy.arcsin(dsin(my_lambda)*dsin(keps))
        delta /= pir
        self.delta = delta
        self.logger.debug("declination, delta, set to %f", delta)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate time zone hour, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.lon < 0:
            # Swap to positive to "round down" negative numbers:
            temp_lon = -1.0*self.lon
            temp_tzh = int(temp_lon/15)
            tz_hour = -1.0*temp_tzh
        else:
            tz_hour = int(self.lon/15)
        self.tz_hour = tz_hour
        self.logger.debug("time zone hour set to %d", tz_hour)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the equation of time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Spencer (1971)
        eot = self.spencer_eot(n)
        self.eot_hour = eot
        self.logger.debug("Equation of Time set to %f", eot)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the longitude correction, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lonc = (15.0*tz_hour - self.lon)/15.0
        self.lc_hour = lonc
        self.logger.debug("longitude corrector set to %f", lonc)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate the solar time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ts_hh = local_hh + eot - lonc
        self.ts_hh = ts_hh

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate the hour angle, degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        w_hh = (360./24.)*(ts_hh - 12.0)
        self.w_hh = w_hh

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = dsin(delta)*dsin(self.lat)
        rv = dcos(delta)*dcos(self.lat)
        self.ru = ru
        self.rv = rv
        self.logger.debug("variable substitute, ru, set to %f", ru)
        self.logger.debug("variable substitute, rv, set to %f", rv)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            # Polar day (no sunset)
            self.logger.debug("polar day---no sunset")
            hs = 180.0
        elif (ru/rv) <= -1.0:
            # Polar night (no sunrise)
            self.logger.debug("polar night---no sunrise")
            hs = 0.0
        else:
            hs = -1.0*ru/rv
            hs = numpy.arccos(hs)
            hs /= pir
        self.hs = hs
        self.logger.debug("sunset angle, hs, set to %f", hs)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate the half-hourly solar radiation flux, W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        io_hh = kGsc*dr*(ru + rv*dcos(w_hh))
        io_hh = io_hh.clip(min=0)
        self.io_wm2 = io_hh

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate the half-hourly PPFD, umol/m^2/s
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_hh = kfFEC*io_hh
        self.ppfd_hh = ppfd_hh

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate the daily solar irradiation, J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 1.10.3, Duffy & Beckman (1993)
        ho = (86400.0/numpy.pi)*kGsc*dr*(ru*pir*hs + rv*dsin(hs))
        self.ho_jm2 = ho
        self.logger.info("daily ET radiation set to %f MJ/m^2", (1.0e-6)*ho)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def berger_tls(self, n):
        """
        Name:     SOLAR_TOA.berger_tls
        Input:    int, day of year
        Output:   tuple,
                  - true anomaly, degrees
                  - true longitude, degrees
        Features: Returns true anomaly and true longitude for a given day
        Depends:  - ke
                  - komega
        Ref:      Berger, A. L. (1978), Long term variations of daily
                  insolation and quaternary climatic changes, J. Atmos. Sci.,
                  35, 2362-2367.
        """
        self.logger.debug("calculating heliocentric longitudes for day %d", n)

        # Variable substitutes:
        xee = ke**2
        xec = ke**3
        xse = numpy.sqrt(1.0 - xee)

        # Mean longitude for vernal equinox:
        xlam = (ke/2.0 + xec/8.0)*(1.0 + xse)*dsin(komega)
        xlam -= xee/4.0*(0.5 + xse)*dsin(2.0*komega)
        xlam += xec/8.0*(1.0/3.0 + xse)*dsin(3.0*komega)
        xlam *= 2.0
        xlam /= pir
        self.logger.debug("mean longitude for vernal equinox set to %f", xlam)

        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        self.logger.debug("mean longitude for day of year set to %f", dlamm)

        # Mean anomaly:
        anm = (dlamm - komega)
        ranm = (anm*pir)
        self.logger.debug("mean anomaly set to %f", ranm)

        # True anomaly:
        ranv = ranm
        ranv += (2.0*ke - xec/4.0)*numpy.sin(ranm)
        ranv += 5.0/4.0*xee*numpy.sin(2.0*ranm)
        ranv += 13.0/12.0*xec*numpy.sin(3.0*ranm)
        anv = ranv/pir

        # True longitude:
        my_tls = anv + komega
        if my_tls < 0:
            my_tls += 360.0
        elif my_tls > 360:
            my_tls -= 360.0
        self.logger.debug("true longitude set to %f", my_tls)

        # True anomaly:
        my_nu = (my_tls - komega)
        if my_nu < 0:
            my_nu += 360.0
        self.logger.debug("true anomaly set to %f", my_nu)

        return(my_nu, my_tls)

    def julian_day(self, y, m, i):
        """
        Name:     SOLAR_TOA.julian_day
        Input:    - int, year (y)
                  - int, month (m)
                  - int, day of month (i)
        Output:   float, Julian Ephemeris Day
        Features: Converts Gregorian date (year, month, day) to Julian
                  Ephemeris Day
        Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical
                  Algorithms
        """
        self.logger.debug("calculating Julian day")
        if m <= 2.0:
            y -= 1.0
            m += 12.0

        a = int(y/100)
        b = 2 - a + int(a/4)

        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde

    def spencer_eot(self, n):
        """
        Name:     SOLAR_TOA.spencer_eot
        Input:    int, day of the year (n)
        Output:   float, equation of time, hours
        Features: Returns the equation of time
        Depends:  - dcos
                  - dsin
        Ref:      Spencer, J.W. (1971), Fourier series representation of the
                  position of the sun, Search, 2 (5), p. 172.
        """
        B = 2.0*numpy.pi*(n - 1.0)/self.kN
        my_eot = (7.5e-6)
        my_eot += (1.868e-3)*dcos(B)
        my_eot -= (3.2077e-2)*dsin(B)
        my_eot -= (1.4615e-2)*dcos(2.0*B)
        my_eot -= (4.0849e-2)*dsin(2.0*B)
        my_eot *= (12.0/numpy.pi)
        return(my_eot)
