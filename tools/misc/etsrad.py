#!/usr/bin/python
#
# etsrad.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-09-10 -- created
# 2014-10-18 -- last updated
#
# ------------
# description:
# ------------
# This script calculates extraterrestrial solar radiation.
# * daily half-hour time series
# * total daylight hours
# 
# -----------
# references:
# -----------
# 1. J.A. Duffie and W.A. Beckman, 2013. "Solar Engineering of Thermal 
#    Processes," 4th ed. John Wiley and Sons, NJ. (p. 37)
# 2. J.W. Spencer, "Fourier series representation of the position of the sun,"
#    Search, vol. 2, p. 172, 1971. Available online:
#    http://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html
# 3. P.I. Cooper, "The absorption of radiation in solar stills," Solar Energy,
#    vol. 12, no. 3, pp. 333-346, 1969. Available online:
#    http://www.sciencedirect.com/science/article/pii/0038092X69900474
# 4. W.B. Stine and M. Geyer, 2001. "Power from the Sun," available online:
#    http://www.powerfromthesun.net/book.html
# 5. H.M. Woolf, "On the computation of solar evaluation angles and the 
#    determination of sunrise and sunset times," in the National Aeronautics 
#    and Space Administration (NASA) Report TM-X-164, September 1968.
# 6. R.L. Snyder and S. Eching, 2002. "Penman-Monteith (hourly) Reference Evapo-
#    transpiration Equations for Estimating ETos and ETrs with Hourly Weather
#    Data," Regents of the University of California, revised 2008. Available 
#    online: http://biomet.ucdavis.edu/Evapotranspiration/PMhrXLS/PMhrDoc.pdf
# 7. D.W. Meek, J.L. Hatfield, T.A. Howell, S.B. Idso, and R.J. Reginato, 
#    "A generalized relationship between photosynthetically active radiation 
#    and solar radiation," Agron. J., vol. 76, pp. 939--945, 1984.
# 8. J.M. Chen, T.A. Black, D.T. Price, and R.E. Carter, "Model for calculating
#    photosynthetic photon flux densities in forest openings on slopes," 
#    Journal of Applied Meteorology, vol. 32, no. 10, pp. 1656--1665, 1993.
# 9. C. Wehrli, "Extraterrestrial solar spectrum," WRC Publication 615, 
#    Physikalisch-Meteorologisches Observatorium/World Radiation Center 
#    (PMO/WRC), Davos Dorf, Switzerland, 7 pp, 1985.
#
# ----------
# changelog:
# ----------
# 01. change lat/lon order in __init__() to lon, lat, day [13.09.13]
# 02. added srad_to_ppfd conversion factor [13.09.13]
# 03. added local_sec time series [13.09.16]
# 04. calculates daylight hours based on 1-sec time series [13.09.16]
# 05. added daylight radiation integrals & their calculation [13.09.16]
# 06. updated calc_hours() [13.09.16]
# --> daily integral of PPFD (units of umol m-2); instead of daylight average
# --> implemented Simpson's rule instead of Trapezoidal rule
# 07. added functions [14.07.16]
# --> earth_period()
# --> earth_velocity()
# --> correct_kepler()
# --> simplified_kepler()
# --> julian_day()
# --> equinox()
# --> get_lambda()
# 08. moved functions into Solar class [14.10.15]
# 09. added daylight savings factor to class input [14.10.17]
# 10. added datetime module [14.10.18]
# 11. added local_time & solar_time as SOLAR class variables [14.10.18]
# 12. added check for ds input [14.10.18] 
#
################################################################################
## IMPORT MODULES:
################################################################################
from sys import exit
import datetime
import numpy

################################################################################
## CLASSES:
################################################################################
class SOLAR:
    """
    Name:     SOLAR
    Features: This class calculates the half-hourly extraterrestrial PPFD 
              [umol m-2 s-1], and the daily extraterrestrial PPFD [umol m-2]
              based on the STASH 2.0 methodology
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    ke = 0.0167   # eccentricity for 2000 CE (Berger, 1978)
    keps = 23.44  # obliquity for 2000 CE, degrees (Berger, 1978)
    kfFEC = 2.04  # From flux to energy conversion, umol/J (Meek et al., 1984)
    kGsc = 1360.8 # Solar constant, W/m^2 (Kopp & Lean, 2011)
    komega = 283. # longitude of perihelion for 2000 CE, degrees (Berger, 1978)
    #
    # List of local time at half-hourly time step: 
    local_hh = numpy.array([0.5*i for i in xrange(48)])
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization 
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lon, lat, n, 
                 y=0, ds=0, drm='loutre', lamm = 'kepler', delm = 'loutre'):
        """
        Name:     SOLAR.__init__
        Input:    - float, longitude (lon)
                  - float, latitude (lat)
                  - int, day of year (n)
                  - int, year (optional)
                  - int, daylight savings (ds)
                  - string, distance method (drm)
                  - string, lambda method (lamm)
                  - string, delta method (delm)
        """
        # Error handle and assign required public variables:
        self.user_year = y
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.user_lat = lat
        if lon > 180.0 or lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.user_lon = lon
        if n < 1 or n > 366:
            print "Day outside range of validity (1 to 366)!"
            exit(1)
        else:
            self.user_day = n
        if ds < 0 or ds > 1:
            print "Set daylight savings time to 0 or 1"
            exit(1)
        else:
            self.ds = ds
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Create datetime series
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.local_time = numpy.array([
            datetime.datetime(y, 1 ,1 ,0, 0, 0) + 
            datetime.timedelta(days=(n-1)) + 
            datetime.timedelta(hours=i) for i in self.local_hh
        ])
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.user_year == 0:
            self.kN = 365
        else:
            self.kN = (
                self.julian_day((y + 1), 1, 1) - 
                self.julian_day(y, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if lamm == 'kepler':
            # Simplified Kepler Method:
            self.my_nu, self.my_lambda = self.map_days(n, y)
        elif lamm == 'woolf':
            # Woolf Method:
            woolf_b = (n - 1.0)*(360.0/self.kN)
            self.my_lambda = (
                279.9348 + woolf_b + 1.914827*self.dsin(woolf_b) - 
                0.079525*self.dcos(woolf_b) + 0.019938*self.dsin(2*woolf_b) - 
                0.00162*self.dcos(2*woolf_b)
            )
            if self.my_lambda < 0:
                self.my_lambda += 360.0
            elif self.my_lambda > 360:
                self.my_lambda -= 360.0
            self.my_nu = (self.my_lambda - self.komega)
            if self.my_nu < 0:
                self.my_nu += 360.0
        elif lamm == 'berger':
            # Berger'78 Method:
            xee = self.ke**2 
            xec = self.ke**3
            xse = numpy.sqrt(1.0 - xee)
            xlam = (
                (self.ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(self.komega) - 
                xee/4.0*(0.5 + xse)*self.dsin(2.0*self.komega) + 
                xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*self.komega)
                )
            xlam = numpy.degrees(2.0*xlam)
            dlamm = xlam + (n - 80.0)*(360.0/self.kN)
            anm = dlamm - self.komega
            ranm = numpy.radians(anm)
            ranv = (ranm + (2.0*self.ke - xec/4.0)*numpy.sin(ranm) + 
                5.0/4.0*xee*numpy.sin(2.0*ranm) + 
                13.0/12.0*xec*numpy.sin(3.0*ranm))
            anv = numpy.degrees(ranv)
            self.my_lambda = anv + self.komega
            if self.my_lambda < 0:
                self.my_lambda += 360.0
            elif self.my_lambda > 360:
                self.my_lambda -= 360.0
            self.my_nu = (self.my_lambda - self.komega)
            if self.my_nu < 0:
                self.my_nu += 360.0
        else:
            print "Heliocentric longitude method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if drm == 'loutre':
            # Eq. 7--9, STASH 2.0 Documentation
            my_rho = (1.0 - self.ke**2)/(1.0 + self.ke*self.dcos(self.my_nu))
            self.dr = (1.0/my_rho)**2
        elif drm == 'klein':
            # Eq. 11, STASH 2.0 Documentation
            self.dr = (
                1.0 + 2.0*self.ke*numpy.cos(2.0*numpy.pi*self.user_day/self.kN)
            )
        else:
            print "Distance factor method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if delm == 'loutre':
            # Eq. 14, STASH 2.0 Documentation (Loutre, 2002)
            self.delta = numpy.arcsin(
                self.dsin(self.my_lambda)*self.dsin(self.keps)
            )
            self.delta = self.delta*(180.0/numpy.pi)
        elif delm == 'cooper':
            # Eq. 16, STASH 2.0 Documentation (Cooper, 1969)
            self.delta = self.keps*numpy.sin(
                2.0*numpy.pi*(n + self.komega)/self.kN
            )
        elif delm == 'circle':
            # Eq. 15, STASH 2.0 Documentation
            self.delta = -1.0*self.keps*numpy.cos(
                2.0*numpy.pi*(n + 10.)/self.kN
            )
        elif delm == 'spencer':
            # Eq. 17, STASH 2.0 Documentation
            spencer_b = (n - 1.0)*(2.0*numpy.pi/self.kN)
            self.delta = (0.006918 - 
                          0.399912*numpy.cos(spencer_b) + 
                          0.070257*numpy.sin(spencer_b) - 
                          0.006758*numpy.cos(2.0*spencer_b) +
                          0.000907*numpy.sin(2.0*spencer_b) - 
                          0.0022697*numpy.cos(3.0*spencer_b) + 
                          0.00148*numpy.sin(3.0*spencer_b))
            self.delta *= (180.0/numpy.pi)
        else:
            print "Declination angle method not recognized!"
            exit(1)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 5. Calculate time zone hour, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if self.user_lon < 0:
            # Swap to positive to "round down" negative numbers:
            temp_lon = -1.0*self.user_lon
            temp_tzh = int(temp_lon/15)
            self.tz_hour = -1.0*temp_tzh
        else:
            self.tz_hour = int(self.user_lon/15)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the equation of time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Eq. 21, STASH 2.0 Documentation
        B = 2.0*numpy.pi*(n - 1.0)/self.kN
        self.eot = 24.0/(2.0*numpy.pi)*(
            (7.5e-6) + (1.868e-3)*self.dcos(B) - (3.2077e-2)*self.dsin(B) - 
            (1.4615e-2)*self.dcos(2.0*B) - (4.0849e-2)*self.dsin(2.0*B)
        )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the longitude correction, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.lc = (15.0*self.tz_hour - self.user_lon)/15.0
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8a. Calculate the solar time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.ts_hh = self.local_hh + self.eot - self.lc - self.ds
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8b. Create solar datetime series
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.solar_time = numpy.array([
            datetime.datetime(y, 1 ,1 ,0, 0, 0) + 
            datetime.timedelta(days=(n-1)) + 
            datetime.timedelta(hours=i) for i in self.ts_hh
        ])
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate the hour angle, degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.w_hh = (360./24.)*(self.ts_hh - 12.0)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(self.delta)*self.dsin(self.user_lat)
        rv = self.dcos(self.delta)*self.dcos(self.user_lat)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Note: ru/rv == tan(delta)*tan(lat)
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            self.hs = 180.0   # Polar day (no sunset)
        elif (ru/rv) <= -1.0:
            self.hs = 0.0     # Polar night (no sunrise)
        else:
            self.hs = (180.0/numpy.pi)*numpy.arccos(-1.0*ru/rv)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate the extraterrestrial solar radiation, W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.ra_hh = self.kGsc*self.dr*(ru + rv*self.dcos(self.w_hh))
        self.ra_hh = self.ra_hh.clip(min=0)
        self.et_ppfd_hh = self.ra_hh*self.kfFEC # converted to umol/m^2/s
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate the daily extraterrestrial solar radiation, J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.ra_d = (86400.0/numpy.pi)*self.kGsc*self.dr*(
            ru*(numpy.pi/180.0)*self.hs + rv*self.dsin(self.hs)
        )
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def dcos(self, x):
        """
        Name:     SOLAR.dcos
        Input:    float, angle, degrees (x)
        Output:   float, cos(x*pi/180)
        Features: Calculates the cosine of an angle given in degrees
        """
        return numpy.cos(x*numpy.pi/180.0)
    #
    def dsin(self, x):
        """
        Name:     SOLAR.dsin
        Input:    float, angle, degrees (x)
        Output:   float, sin(x*pi/180)
        Features: Calculates the sine of an angle given in degrees
        """
        return numpy.sin(x*numpy.pi/180.0)
    #
    def julian_day(self, y, m, i):
        """
        Name:     SOLAR.julian_day
        Input:    - int, year (y)
                  - int, month (m)
                  - int, day of month (i)
        Output:   float, Julian Ephemeris Day
        Features: Converts Gregorian date (year, month, day) to Julian 
                  Ephemeris Day
        Ref:      Eq. 7.1, Meeus, J. (1991), Ch.7 "Julian Day," Astronomical 
                  Algorithms
        """
        if m <= 2.0:
            y -= 1.0
            m += 12.0
        #
        a = int(y/100)
        b = 2 - a + int(a/4)
        #
        jde = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + i + b - 1524.5
        return jde
    #
    def equinox(self, year, opt=0):
        """
        Name:     SOLAR.equinox
        Input:    - int, year (year)
                  - int, option (opt)
                    0: vernal equinox     1: summer solstice
                    2: autumnal equinox   3: winter solstice
        Output:   float, day of the year
        Features: Calculates the day of the year on which seasonal dates fall
        Depends:  julian_day
        Ref:      J. Meeus (1991), Ch.26 "Equinoxes and solstices," 
                  Astronomical Algorithms
        """
        # Table 26.C (Meeus, 1991)
        periodic_terms = numpy.array([
            # A    B          C
            485, 324.96,   1934.136,
            203, 337.23,  32964.467,
            199, 342.08,     20.186,
            182,  27.85, 445267.112,
            156,  73.14,  45036.886,
            136, 171.52,  22518.443,
            77, 222.54,  65928.934,
            74, 296.72,   3034.906,
            70, 243.58,   9037.513,
            58, 119.81,  33718.147,
            52, 297.17,    150.678,
            50,  21.02,   2281.226,
            45, 247.54,  29929.562,
            44, 325.15,  31555.956,
            29,  60.93,   4443.417,
            18, 155.12,  67555.328,
            17, 288.79,   4562.452,
            16, 198.04,  62894.029,
            14, 199.76,  31436.921,
            12,  95.39,  14577.848,
            12, 287.11,  31931.756,
            12, 320.81,  34777.259,
            9, 227.73,   1222.114,
            8,  15.45,  16859.074
            ])
        #
        # Table 26.A (Meeus, 1991)
        jde_table_a = numpy.array([
            numpy.array([1721139.29189, 365242.13740,  0.06134,  
                         0.00111, -0.00071]), # March equinox
            numpy.array([1721233.25401, 365241.72562, -0.05323,  
                         0.00907,  0.00025]), # June solstice
            numpy.array([1721325.70455, 365242.49558, -0.11677, 
                         -0.00297,  0.00074]), # September equinox
            numpy.array([1721414.39987, 365242.88257, -0.00769, 
                         -0.00933, -0.00006])  # December solstice
            ])
        #
        # Table 26.B (Meeus, 1991)
        jde_table_b = numpy.array([
            numpy.array([2451623.80984, 365242.37404,  0.05169, 
                         -0.00411, -0.00057]), # March equinox
            numpy.array([2451716.56767, 365241.62603,  0.00325,  
                         0.00888, -0.00030]), # June solstice
            numpy.array([2451810.21715, 365242.01767, -0.11575,  
                         0.00337,  0.00078]), # September equinox
            numpy.array([2451900.05952, 365242.74049, -0.06223, 
                         -0.00823,  0.00032])  # December solstice
            ])
        #
        if year < 1000:
            # Use Table 26.A for years -1000 to 1000
            jde_table = jde_table_a
            y = year/1000.0
        else:
            # Use Table 26.B for years 1000 to 3000
            jde_table = jde_table_b
            y = (year - 2000.0)/1000.0
        #
        # Calculate the mean equinox based on Table 26.A or 26.B:
        jde0 = (
       	    jde_table[opt,0] +
       	    jde_table[opt,1]*y +
       	    jde_table[opt,2]*y*y +
       	    jde_table[opt,3]*y*y*y +
       	    jde_table[opt,4]*y*y*y*y
        )
        #
        # Calculate the other three terms:
        t = (jde0 - 2451545.0)/36525.0
        w = (35999.373*t) - 2.47
        dl = 1.0 + 0.0334*self.dcos(w) + 0.0007*self.dcos(2*w)
        #
        # Compute the sum of the 24 periodic terms:
        s = 0
        j = 0
        for i in xrange(24):
            s += periodic_terms[j]*self.dcos(
                (periodic_terms[j+1] + (periodic_terms[j+2]*t))
            )
            j += 3
        #
        # Calculate the JDE (Julian Ephemeris Day) of the equinox/solstice:
        jde = jde0 + ((s*0.00001)/dl)
        #
        # Calcuate the JDE for the first day of this year:
        jde_year = self.julian_day(year, 1, 1)
        #
        # Return the day of year:
        return (jde - jde_year + 1)
    #
    def earth_velocity(self, lon):
        """
        Name:     SOLAR.earth_velocity
        Input:    longitude(s) w.r.t. the perihelion (lon)
        Output:   angular velocity(ies), rad/day (w)
        Features: Calculates earth's angular velocity at a given longitude
                  relative to the perihelion
        Ref:      Kepler's Second Law of Planetary Motion
        """
        w = (
            2.0*numpy.pi/self.kN*(
                1.0 + self.ke*self.dcos(lon)
                )**2.0*(1.0 - self.ke**2)**(-1.5)
        )
        return w
    #
    def correct_kepler(self, lon, t):
        """
        Name:     SOLAR.correct_kepler
        Input:    - longitudes w.r.t. the perihelion; at least two points 
                    before and one point after the 180 degree cross-over (lon)
                  - days associates with the longitudes (t)
        Output:   numpy.ndarray, days (y)
        Features: Corrects the array of days output from simplified_kepler due 
                  to the asymptote in the arctan at 180 degrees
        """
        # Find indexes of positive and negative values:
        neg_points = numpy.where(t<0)[0]
        pos_points = numpy.where(t>=0)[0]
        #
        # Linearly extrapolate the positive position of the first negative 
        # value for:
        #   xa ... ya  <-- known value pair a
        #   xb ... yb  <-- known value pair b
        #   xi ... yi  <-- extrapolated value pair i
        # where all variables are known except for yi, such that:
        #   (yi-ya)/(yb-ya) = (xi-xa)/(xb-xa)
        # and solving for yi:
        #   yi = ya + (yb-ya)*(xi-xa)/(xb-xa)
        xa = pos_points[-2]
        xb = pos_points[-1]
        xi = neg_points[0]
        #
        ya = t[xa]
        yb = t[xb]
        yi = ya + (yb - ya)*(xi - xa)/(xb - xa)
        #
        # Calculate the difference between the pos. and neg. values at y:
        diff_y = yi - t[xi]
        #
        # Add the difference to the negative values to push them to the proper
        # positive position, i.e.:
        #
        #      +y                     +y
        #       |     o               |           o
        #       |   o                 |         o  
        #       | o                   |       o    
        #       |------------- +x  => |     o      
        #       |           o         |   o        
        #       |         o           | o          
        #       |       o             |------------ +x
        #
        #       (a) Original Data     (b) Corrected
        #
        y = numpy.copy(t)
        y[neg_points] = y[neg_points] + diff_y
        #
        return (y)
    #
    def simplified_kepler(self,lon):
        """
        Name:     SOLAR.simplified_kepler
        Input:    float, longitude w.r.t. the perihelion (lon)
        Output:   float, days traveled from the perihelion (t)
        Features: Calculates the time (days) of earth's orbit around the sun
                  for a given longitude using a simplified Kepler approach
        Depends:  correct_kepler
        Ref:      Kepler's Second Law of Planetary Motion
        """
        # Make lon into numpy array, if it is not already:
        if isinstance(lon, numpy.float) or isinstance(lon, numpy.int):
            lon = numpy.array([lon,])
        elif not isinstance(lon, numpy.ndarray):
            lon = numpy.array(lon)
        #
        # Calculate the days
        # NOTE: arctan(inf) = pi/2
        t = (
            self.kN/numpy.pi*(
                numpy.arctan(
                    self.dsin(lon)/(self.dcos(lon) + 1)
                ) - self.ke*self.dsin(lon)
            )
        )
        #
        # Correct the days, if necessary:
        if len(lon) > 1:
            t = self.correct_kepler(lon, t)
        #
        return t
    #
    def map_days(self, n, y):
        """
        Name:     SOLAR.map_days
        Input:    - int, day of the year (n)
                  - int, year (y)
        Output:   float, longitude relative to the vernal equinox (lamda_doy)
        Features: Computes earth's longitude relative to the vernal equinox 
                  for a given day of the year
        Depends:  Functions:
                  - earth_velocity
                  - equinox
                  - simplified_kepler
        """
        # Create longitude field and compute the days for orbit:
        lon_nu = numpy.array([1.0*i for i in xrange(361)])
        day_nu = self.simplified_kepler(lon_nu)
        #
        # Compute the angle of the vernal equinox w.r.t. the perihelion
        # i.e., the explementary angle of komega:
        wp = (360.0 - self.komega)
        #
        # Calculate the length of time it takes earth to travel from the
        # perihelion (t=0) to the vernal equinox:
        if wp < 180.0:
            days_to_ve = self.simplified_kepler(wp)[0]
        else:
            lon_temp = numpy.array([179.0, 180.0, 181.0, wp])
            days_to_ve = self.simplified_kepler(lon_temp)[3]
        #
        # Get day of year of vernal equinox
        if y == 0:
            day_ve = 80.0
        else:
            day_ve = self.equinox(y)
        #
        # Calculate the offset between days to and day of vernal equinox:
        offset = (day_ve - days_to_ve)
        #
        # Calculate the calendar days and set between 0 and kN:
        calendar_days = (day_nu + offset)
        calendar_days[numpy.where(calendar_days >= self.kN)] -= self.kN
        #
        # Check to see if n is listed in calendar:
        if n in calendar_days:
            icalendar = numpy.where(calendar_days==n)[0][0]
            nu_doy = lon_nu[icalendar]
        else:
            # Find the calendar day the precedes doy:
            calendar_temp = numpy.sort(numpy.append(calendar_days, [n,]))
            dbefore = calendar_temp[numpy.where(calendar_temp==n)[0][0]-1]
            #
            # Get the index of the preceding day:
            ibefore = numpy.where(calendar_days == dbefore)[0][0]
            #
            # Get the angular velocity for the longitude of the preceding day:
            vbefore = self.earth_velocity(lon_nu[ibefore])
            #
            # Calculate the delta time
            dt = (n - dbefore)
            #
            # Calculate the new longitude, degrees:
            nu_doy = lon_nu[ibefore] + (vbefore*dt)*180.0/numpy.pi
        #
        # Convert nu to lambda:
        lamda_doy = nu_doy + self.komega
        if lamda_doy >= 360:
            lamda_doy -= 360
        #
        return (nu_doy, lamda_doy)

################################################################################
## FUNCTIONS:
################################################################################


################################################################################
## MAIN:
################################################################################
# Coordinates (degrees):
#lat = 51.408
#lon = -0.64
my_lat = 32.75
my_lon = -96.78

# Julian day of the year:
my_day = 83

my_solar = SOLAR(my_lon, my_lat, my_day)
