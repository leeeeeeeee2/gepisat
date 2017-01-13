#!/usr/bin/python
#
# lue.py
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
# (GePiSaT) Model of the terrestrial biosphere: Part 1 — Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.

###############################################################################
# IMPORT MODULES
###############################################################################
import logging
import os

import numpy

from .const import kc
from .const import kco
from .const import kphio
from .const import kPo
from .const import kR
from .const import gs25
from .const import dha
from .const import kc25
from .const import ko25
from .const import dhac
from .const import dhao


###############################################################################
# CLASSES
###############################################################################
class LUE:
    """
    Name:     LUE
    Features: This class stores monthly LUE estimates and writes results to
              file.
    History   Version 3.0.0-dev
              - class separated from model [16.01.17]
              - Python 2/3 supported print statements [16.01.17]
              - moved constants to const.py [16.07.22]
              - fixed xrange for Python 3 support [16.07.22]
              - moved to gepisat package [16.07.22]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     LUE.__init__
        Features: Initializes dictionaries for the light-use efficiency model
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("LUE class initialized")

        # Dictionary of stations' monthly values & their units:
        # * this is for printing to file
        self.station_vals = {}
        self.station_units = {'Timestamp': 'NA',
                              'GPP': 'mol_m2',
                              'GPP_err': 'mol_m2',
                              'fPAR': 'NA',
                              'PPFD': 'mol_m2',
                              'VPD': 'kPa',
                              'CPA': 'NA',
                              'Tair': 'degC',
                              'CO2': 'ppm',
                              'Patm': 'Pa',
                              'elv': 'm',
                              'Iabs': 'mol_m2',
                              'ca': 'Pa',
                              'Gs': 'Pa',
                              'D': 'Pa',
                              'K': 'Pa',
                              'ns': 'NA',
                              'fa': 'NA',
                              'beta1': 'NA',
                              'beta2': 'NA',
                              'GPP_hat1': 'mol_m2',
                              'GPP_hat2': 'mol_m2'}

        # Define station value header line:
        header_vals = ['Timestamp', 'GPP', 'GPP_err', 'fPAR', 'PPFD',
                       'VPD', 'CPA', 'Tair', 'CO2', 'Patm', 'elv', 'Iabs',
                       'ca', 'Gs', 'D', 'K', 'ns', 'fa', 'beta1', 'beta2',
                       'GPP_hat1', 'GPP_hat2']
        header = ''
        for k in header_vals:
            header += k
            v = self.station_units[k]
            if v == 'NA':
                header += ','
            else:
                header += '.'
                header += v
                header += ','
        header = header.strip(',')
        header += '\n'
        self.value_header = header

        # Dictionary of stations' light-use efficiency model variables
        # * this is for modeling
        self.st_lue_vars = {}
        self.lue_var_units = {'GPP': 'mol_m2',
                              'GPP_err': 'mol_m2',
                              'Iabs': 'mol_m2',
                              'ca': 'Pa',
                              'Gs': 'Pa',
                              'D': 'Pa',
                              'K': 'Pa',
                              'ns': 'NA',
                              'fa': 'NA',
                              'beta1': 'NA',
                              'beta2': 'NA'}

        # Dictionary of stations' light-use efficiency model fit/fitness:
        # * this is for model results
        self.station_lue = {}

        # Define standard viscosity of water, Pa s
        self.n25 = self.viscosity_h2o(25.0, kPo)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def add_station_val(self, station, month, gpp, gpp_err,
                        fpar, ppfd, vpd, alpha, tair, co2, patm, elv):
        """
        Name:     LUE.add_station_val
        Input:    - string, station name (station)
                  - datetime.date, current month (month)
                  - float, monthly GPP, mol/m2 (gpp)
                  - float, associated GPP error, mol/m2 (gpp_err)
                  - float, monthly FAPAR (fpar)
                  - float, monthly PPFD, mol/m2 (ppfd)
                  - float, monthly VPD, kPa (vpd)
                  - float, monthly CPA (alpha)
                  - float, monthly air temp, degC (tmp)
                  - float, annual atm. CO2, ppm (co2)
                  - float, monthly atm. pressure, Pa (patm)
                  - float, elevation, m (elv)
        Output:   None.
        Features: Appends a set of monthly values to the value dictionary
        Depends:  - calc_gstar
                  - calc_k
                  - nxgn
                  - viscosity_h2o
        """
        # Calculate lue variables & place into tuple:
        if fpar is not None:
            iabs = fpar*ppfd                     # mol/m2, abs. PPFD
        else:
            iabs = numpy.nan

        if co2 is not None:
            ca = (1.e-6)*co2*patm                # Pa, atms. CO2
        else:
            ca = numpy.nan

        if vpd is not None:
            d = (1e3)*vpd                        # Pa, vapor pressure deficit
        else:
            d = numpy.nan

        if tair is not None:
            gs = self.calc_gstar(tair)           # Pa, photores. comp. point
            k = self.calc_k(tair, patm)          # Pa, Michaelis-Menten coef.
            ns = self.viscosity_h2o(tair, patm)  # Pa s, viscosity
            ns /= self.n25                       # unitless, water viscosity
        else:
            tair = numpy.nan
            gs = numpy.nan
            k = numpy.nan
            ns = numpy.nan

        if alpha is not None:
            fa = (alpha/1.26)**(0.25)            # unitless, func. of alpha
        else:
            fa = numpy.nan

        beta1, beta2 = self.beta_estimate(ca, d, k, gs, ns, tair, elv)
        yhat1 = self.nxtgn(iabs, ca, gs, d, k, ns, fa, beta1)
        yhat2 = self.nxtgn(iabs, ca, gs, d, k, ns, fa, beta2)

        # Place value parameters into tuples:
        val_params = (month, gpp, gpp_err, fpar, ppfd, vpd, alpha, tair, co2,
                      patm, elv, iabs, ca, gs, d, k, ns, fa, beta1, beta2,
                      yhat1, yhat2)
        var_params = (gpp, gpp_err, iabs, ca, gs, d, k, ns, fa, beta1, beta2)

        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Station Values
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize list if station key doesn't exist:
        if station not in self.station_vals.keys():
            self.station_vals[station] = numpy.array(
                val_params,
                dtype={'names': ('Timestamp', 'GPP', 'GPP_err',
                                 'fPAR', 'PPFD', 'VPD',
                                 'CPA', 'Tair', 'CO2',
                                 'Patm', 'elv', 'Iabs',
                                 'ca', 'Gs', 'D',
                                 'K', 'ns', 'fa',
                                 'beta1', 'beta2',
                                 'GPP_hat1', 'GPP_hat2'),
                       'formats': ('O', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4',
                                   'f4', 'f4')},
                ndmin=1
            )
        else:
            # Add new parameters to list:
            temp_array = numpy.array(
                val_params,
                dtype={'names': ('Timestamp', 'GPP', 'GPP_err',
                                 'fPAR', 'PPFD', 'VPD',
                                 'CPA', 'Tair', 'CO2',
                                 'Patm', 'elv', 'Iabs',
                                 'ca', 'Gs', 'D',
                                 'K', 'ns', 'fa',
                                 'beta1', 'beta2',
                                 'GPP_hat1', 'GPP_hat2'),
                       'formats': ('O', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4', 'f4',
                                   'f4', 'f4',
                                   'f4', 'f4')},
                ndmin=1
            )
            self.station_vals[station] = numpy.append(
                self.station_vals[station], temp_array, axis=0)

        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Station Variables
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Eliminate data points where VPD or air temperature are negative:
        if d > 0 and tair > 0:
            # Initialize array if station key doesn't exist:
            if station not in self.st_lue_vars.keys():
                self.st_lue_vars[station] = numpy.array(
                    var_params,
                    dtype={
                        'names': ('GPP', 'GPP_err', 'Iabs',
                                  'ca', 'Gs', 'D',
                                  'K', 'ns', 'fa',
                                  'beta1', 'beta2'),
                        'formats': ('f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4')
                        },
                    ndmin=1
                )
            else:
                temp_array = numpy.array(
                    var_params,
                    dtype={
                        'names': ('GPP', 'GPP_err', 'Iabs',
                                  'ca', 'Gs', 'D',
                                  'K', 'ns', 'fa',
                                  'beta1', 'beta2'),
                        'formats': ('f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4')
                        },
                    ndmin=1
                )
                self.st_lue_vars[station] = numpy.append(
                    self.st_lue_vars[station], temp_array, axis=0)

    def beta_estimate(self, my_ca, my_d, my_k, my_gs, my_ns, my_t, my_z):
        """
        Name:     LUE.beta_estimate
        Input:    - float, atmospheric CO2 concentration, Pa (my_ca)
                  - float, vapor pressure deficit, Pa (my_d)
                  - float, Michaelis-Menten coeff, Pa (my_k)
                  - float, photorespiratory comp point, Pa (my_gs)
                  - float, viscosity, unitless (my_ns)
                  - float, air temperature, deg C (my_t)
                  - float, elevation, m (my_z)
        Output:   tuple
                  - float, predicted beta from simple expression (beta_p1)
                  - float, predicted beta from
        Features: Returns an estimate for beta based on the Wang Han equation
        """
        if not numpy.isfinite(my_z):
            my_z = 0.

        vals = numpy.array([my_ca, my_d, my_k, my_gs, my_ns, my_t, my_z])
        if numpy.isfinite(vals).all():
            whe = numpy.exp(
                1.19 +
                0.0545*(my_t - 25.0) -       # T in deg C
                0.5*numpy.log(1e-3*my_d) -   # D in kPa
                0.0815*(1e-3*my_z)           # z in km
            )
            chi = whe/(1. + whe)

            beta_p1 = 1.6*my_ns*my_d*(chi**2)
            beta_p1 /= (1. - chi)**2
            beta_p1 /= my_k

            beta_p2 = 1.6*my_ns*my_d
            beta_p2 /= (my_k + my_gs)
            beta_p2 *= (chi*my_ca - my_gs)**2
            beta_p2 /= (my_ca**2)
            beta_p2 /= (chi - 1.)**2
        else:
            beta_p1 = numpy.nan
            beta_p2 = numpy.nan

        return (beta_p1, beta_p2)

    def calc_gstar(self, tc):
        """
        Name:     LUE.calc_gstar
        Input:    float, air temperature, degrees C (tc)
        Output:   float, gamma-star, Pa (gs)
        Features: Returns the temperature-dependent photorespiratory
                  compensation point, Gamma star (Pascals), based on constants
                  derived from Bernacchi et al. (2001) study.
        Ref:      Bernacchi et al. (2001), Improved temperature response
                  functions for models of Rubisco-limited photosynthesis,
                  Plant, Cell and Environment, 24, 253--259.
        """
        gs = gs25*numpy.exp(dha*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        return gs

    def calc_k(self, tc, patm):
        """
        Name:     LUE.calc_k
        Input:    - float, air temperature, degrees C (tc)
                  - float, atmospheric pressure, Pa (patm)
        Output:   float (mmk)
        Features: Returns the temperature & pressure dependent Michaelis-Menten
                  coefficient, K (Pascals).
        Ref:      Bernacchi et al. (2001), Improved temperature response
                  functions for models of Rubisco-limited photosynthesis,
                  Plant, Cell and Environment, 24, 253--259.
        """
        vc = kc25*numpy.exp(dhac*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        vo = ko25*numpy.exp(dhao*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        k = vc*(1 + kco*(1e-6)*patm/vo)
        return k

    def density_h2o(self, tc, p):
        """
        Name:     LUE.density_h2o
        Input:    - float, air temperature (tc), degrees C
                  - float, atmospheric pressure (p), Pa
        Output:   float, density of water, kg/m^3
        Features: Calculates density of water at a given temperature and
                  pressure using the Tumlirz Equation
        Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
                  pure water and sea water, Tech. Rept., Marine Physical
                  Laboratory, San Diego, CA.
        """
        # Calculate lambda, (bar cm^3)/g:
        my_lambda = 1788.316
        my_lambda += 21.55053*tc
        my_lambda += -0.4695911*tc*tc
        my_lambda += (3.096363e-3)*tc*tc*tc
        my_lambda += -(7.341182e-6)*tc*tc*tc*tc

        # Calculate po, bar
        po = 5918.499
        po += 58.05267*tc
        po += -1.1253317*tc*tc
        po += (6.6123869e-3)*tc*tc*tc
        po += -(1.4661625e-5)*tc*tc*tc*tc

        # Calculate vinf, cm^3/g
        vinf = 0.6980547
        vinf += -(7.435626e-4)*tc
        vinf += (3.704258e-5)*tc*tc
        vinf += -(6.315724e-7)*tc*tc*tc
        vinf += (9.829576e-9)*tc*tc*tc*tc
        vinf += -(1.197269e-10)*tc*tc*tc*tc*tc
        vinf += (1.005461e-12)*tc*tc*tc*tc*tc*tc
        vinf += -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc
        vinf += (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc
        vinf += -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

        # Convert pressure to bars (1 bar = 100000 Pa)
        pbar = (1e-5)*p

        # Calculate the specific volume (cm^3 g^-1):
        v = vinf + my_lambda/(po + pbar)

        # Convert to density (g cm^-3) -> 1000 g/kg; 1e6 cm^3/m^3 -> kg/m^3:
        rho = (1e3/v)

        return rho

    def nxtgn(self, iabs, ca, gs, d, k, ns, fa, beta):
        """
        Name:     LUE.nxtgn
        Input:    - float, 'Iabs' : mol/m2, fAPARxPPFD
                  - float, 'ca' : Pa, atmospheric CO2
                  - float, 'Gs' : Pa, photores. comp. point
                  - float, 'D' : Pa, vapor pressure deficit
                  - float, 'K' : Pa, Michaelis-Menten coeff.
                  - float, 'ns' : mPa s, viscosity of water
                  - float, 'fa' : unitless, function of alpha
                  - float, beta parameter (beta)
        Output:   float, estimate of GPP (gpp)
        Features: Returns an estimate of GPP based on the next-generation light
                  and water use efficiency model.
        """
        # Define default GPP return value:
        gpp = numpy.nan

        vals = numpy.array([iabs, ca, gs, d, k, ns, fa, beta])
        if numpy.isfinite(vals).all():
            # Define variable substitutes:
            vdcg = ca - gs
            vacg = ca + 2.*gs
            vbkg = beta*(k + gs)

            # Check for negatives:
            if vbkg > 0:
                vsr = numpy.sqrt(1.6*ns*d/(vbkg))

                # Based on the m' formulation (see Regressing_LUE.pdf)
                m = vdcg/(vacg + 3.*gs*vsr)
                mpi = m**2 - kc**(2./3.)*(m**(4./3.))

                # Check for negatives:
                if mpi > 0:
                    mp = numpy.sqrt(mpi)
                    gpp = kphio*iabs*fa*mp

        return gpp

    def viscosity_h2o(self, tc, p):
        """
        Name:     LUE.viscosity_h2o
        Input:    - float, ambient temperature (tc), degrees C
                  - float, ambient pressure (p), Pa
        Return:   float, viscosity of water (mu), Pa s
        Features: Calculates viscosity of water at a given temperature and
                  pressure.
        Depends:  density_h2o
        Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
                  Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
                  international formulation for the viscosity of H2O, J. Phys.
                  Chem. Ref. Data, Vol. 38(2), pp. 101-125.
        """
        # Define reference temperature, density, and pressure values:
        tk_ast = 647.096      # Kelvin
        rho_ast = 322.0       # kg/m^3
        mu_ast = (1e-6)       # Pa s

        # Get the density of water, kg/m^3
        rho = self.density_h2o(tc, p)

        # Calculate dimensionless parameters:
        tbar = (tc + 273.15)/tk_ast
        tbarx = tbar**(0.5)
        tbar2 = tbar**2
        tbar3 = tbar**3
        rbar = rho/rho_ast

        # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
        mu0 = 1.67752
        mu0 += 2.20462/tbar
        mu0 += 0.6366564/tbar2
        mu0 += -0.241605/tbar3
        mu0 = 1e2*tbarx/mu0

        # Create Table 3, Huber et al. (2009):
        hj0 = (0.520094, 0.0850895, -1.08374, -0.289555, 0., 0.)
        hj1 = (0.222531, 0.999115, 1.88797, 1.26613, 0., 0.120573)
        hj2 = (-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.)
        hj3 = (0.161913,  0.257399, 0., 0., 0., 0.)
        hj4 = (-0.0325372, 0., 0., 0.0698452, 0., 0.)
        hj5 = (0., 0., 0., 0., 0.00872102, 0.)
        hj6 = (0., 0., 0., -0.00435673, 0., -0.000593264)
        h = hj0 + hj1 + hj2 + hj3 + hj4 + hj5 + hj6
        h_array = numpy.reshape(numpy.array(h), (7, 6))

        # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
        mu1 = 0.
        ctbar = (1./tbar) - 1.
        for i in range(6):
            coef1 = numpy.power(ctbar, i)
            coef2 = 0.
            for j in range(7):
                coef2 += h_array[j][i]*numpy.power((rbar - 1.), j)
            mu1 += coef1*coef2
        mu1 = numpy.exp(rbar*mu1)

        # Calculate mu_bar (Eq. 2, Huber et al., 2009)
        #   assumes mu2 = 1
        mu_bar = mu0*mu1

        # Calculate mu (Eq. 1, Huber et al., 2009)
        mu = mu_bar*mu_ast    # Pa s

        return mu

    def write_out_val(self, station, out_file):
        """
        Name:     LUE.write_out_val
        Input:    - string, station name (station)
                  - string, output file (out_file)
        Output:   None.
        Features: Writes to file the monthly values associated with the light
                  use efficiency equation for a given station
        """
        self.logger.info("Writing LUE values for station %s", station)

        # Create file if it doesn't exist:
        if not os.path.isfile(out_file):
            try:
                f = open(out_file, 'w')
            except IOError:
                self.logger.error("Cannot write to file %s", out_file)
            else:
                f.write(self.value_header)
                f.close()

        # Print if station has data:
        if station in self.station_vals.keys():
            for t in self.station_vals[station]:
                try:
                    f = open(out_file, 'a')
                except IOError:
                    self.logger.error("Cannot append to file %s", out_file)
                else:
                    f.write(("%s,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,"
                             "%f,%f\n") % tuple(t))
                    f.close()
