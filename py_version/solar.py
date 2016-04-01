#!/usr/bin/python
#
# model.py
#
# VERSION 2.2.0-dev
#
# LAST UPDATED: 2016-04-01
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

###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import logging

import numpy

from const import ke
from const import keps
from const import kfFEC
from const import kGsc
from const import komega
from const import pir
from utilities import dcos
from utilities import dsin


###############################################################################
# CLASSES
###############################################################################
class SOLAR_TOA:
    """
    Name:     SOLAR_TOA
    Features: This class calculates the half-hourly extraterrestrial PPFD
              [umol m-2 s-1], and the daily solar irradiation Ho [J m-2]
              based on the SPLASH model
    History:  Version
              - import global constants [16.04.01]
              - import utility functions [16.04.01]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    # @TODO separate functions from init

    def __init__(self, lon, lat, n, y=0):
        """
        Name:     SOLAR_TOA.__init__
        Input:    - float, longitude (lon)
                  - float, latitude (lat)
                  - int, day of year (n)
                  - int, year (y)
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("SOLAR_TOA class called")

        # Error handle and assign required public variables:
        if y == 0:
            self.year = 2001
        else:
            self.year = y
        if lat > 90.0 or lat < -90.0:
            print("Latitude outside range of validity (-90 to 90)!")
        else:
            self.lat = lat
        if lon > 180.0 or lon < -180.0:
            print("Longitude outside range of validity (-180 to 180)!")
        else:
            self.lon = lon
        if n < 1 or n > 366:
            print("Day outside range of validity (1 to 366)!")
        else:
            self.doy = n

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Create local time series, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        local_hh = numpy.array([0.5*i for i in xrange(48)])
        self.local_time = numpy.array([
            datetime.datetime(self.year, 1, 1, 0, 0, 0) +
            datetime.timedelta(days=(n-1)) +
            datetime.timedelta(hours=i) for i in local_hh
        ])

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            self.kN = 365.
        else:
            self.kN = self.julian_day((y + 1), 1, 1) - self.julian_day(y, 1, 1)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        #self.nu_deg = my_nu
        #self.lambda_deg = my_lambda

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        my_rho = (1.0 - self.ke**2)/(1.0 + self.ke*dcos(my_nu))
        dr = (1.0/my_rho)**2
        #self.dr = dr

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        delta = numpy.arcsin(dsin(my_lambda)*dsin(keps))
        delta /= pir
        #self.delta_deg = delta

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate time zone hour, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if lon < 0:
            # Swap to positive to "round down" negative numbers:
            temp_lon = -1.0*lon
            temp_tzh = int(temp_lon/15)
            tz_hour = -1.0*temp_tzh
        else:
            tz_hour = int(lon/15)
        #self.tz_hour = tz_hour

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the equation of time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Spencer (1971)
        eot = self.spencer_eot(n)
        #self.eot_hour = eot

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the longitude correction, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lc = (15.0*tz_hour - lon)/15.0
        #self.lc_hour = lc

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate the solar time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ts_hh = local_hh + eot - lc
        #self.ts_hh = ts_hh

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate the hour angle, degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        w_hh = (360./24.)*(ts_hh - 12.0)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = dsin(delta)*dsin(lat)
        rv = dcos(delta)*dcos(lat)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            # Polar day (no sunset)
            hs = 180.0
        elif (ru/rv) <= -1.0:
            # Polar night (no sunrise)
            hs = 0.0
        else:
            hs = -1.0*ru/rv
            hs = numpy.arccos(hs)
            hs /= pir
        #self.hs_deg = hs

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate the half-hourly solar radiation flux, W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        io_hh = kGsc*dr*(ru + rv*dcos(w_hh))
        io_hh = io_hh.clip(min=0)
        #self.io_wm2 = io_hh

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
