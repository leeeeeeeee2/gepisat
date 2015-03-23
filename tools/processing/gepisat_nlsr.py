#!/usr/bin/python
#
# gepisat_nlsr.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-11-18 -- created
# 2015-03-23 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script performs non-linear least squared regression for the next-gen
# LUE model.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. added LUE class from GePiSaT model [14.11.19]
# --> added calc_statistics class function based on FLUX_PARTI
# 02. added calc_lue function from GePiSaT model [14.11.19]
# --> completely revised for next_gen_lue (m' formula)
# 03. added prediction functions for phi_o and beta [14.11.22]
# 04. updated plotting (PDF file) [14.11.24]
# 05. added function of alpha to next_gen_lue & calc_lue [14.11.24]
# 06. general housekeeping [14.12.11]
# 07. added get_lue() for quick processing with plot_mo_lue.py [14.12.12]
# 08. updated r-squared calc using scipy.stats.linregress [15.01.29]
# 09. reworking for scipy v0.9.0
#     --> assuming constant phio (Long et al., 1993)
# 10. updated LUE class [15.02.01]
#     --> added calc_gstar, calc_k, viscosity_h2o
#     --> added dictionary for LUE variables & units
# 11. imported scipy.optimize.leastsq [15.02.07]
# 12. created func and func_der functions [15.02.07]
# 13. moved next_gen_lue, lue_resid and lue_jacob into LUE class [15.02.12]
# 14. changed density equation to Tumlirz (Fisher et al., 1975) [15.02.12]
# 15. error estimate to fitted beta parameter [15.02.13]
# 16. added checks for negatives in next_gen_lue and lue_der funcs [15.02.15]
# 17. created get_elv function [15.02.16]
# 18. created wang_han_eq function [15.02.16]
# 19. added elevation to LUE class variable inputs [15.02.16]
# 20. updated beta_estimate function & added it to st_lue_vars [15.02.16]
# 21. plotting of GPP based on seasonal beta estimates [15.02.19]
# 22. created model_fitness function w/ R2 and Willmott IA [15.02.19]
# 23. updated write_out_val [15.02.20]
# 24. added LUE variables & GPP estimates to station_vals  [15.02.19]
# 25. added calc_mgs_beta to LUE class [15.02.27]
# 26. overhaul on model_fitness function [15.02.27]
#     --> definitions/distinctions between R-squared, r, and IOA
# 27. new plots for seasonal, MGS beta estimates + basic v nxtgn [15.02.27]
# 28. added ground-state simple-formula global beta [15.03.17]
# 29. changed calc_mgs_beta to calc_mgs_params [15.03.17]
# 30. deleted calc_hessian function [15.03.17]
# 31. started on crop plots (C4 + irrigation) [15.03.17]
# 32. removed fa from mean growing season calculation [15.03.18]
# 33. work on wetland over-estimation [15.03.18]
# 34. new gamma star equation for temp and press [15.03.23]
# 
# ~~~~~
# todo:
# ~~~~~
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import glob
import numpy
import os.path
from scipy.optimize import leastsq
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

###############################################################################
## CLASSES:
###############################################################################
class LUE:
    """
    Name:     LUE
    Features: This class stores monthly LUE estimates and writes results to
              file.
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    kc = 0.41          # Jmax cost coefficient
    kphio = 0.093      # intrinsic quantum efficiency (Long et al., 1993)
    kbeta = 244.033    # ground-state simple-formula global beta
    kPo = 101325.      # standard atmosphere, Pa (Allen, 1973)
    kTo = 25.          # base temperature, deg C (Prentice, unpublished)
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     LUE.__init__
        Features: Initializes dictionaries for the light-use efficiency model
        """
        # Dictionary of stations' monthly values & their units:
        # * this is for printing to file
        self.station_vals = {}
        self.station_units = {'Timestamp' : 'NA', 
                              'GPP' : 'mol_m2', 
                              'GPP_err' : 'mol_m2', 
                              'fPAR' : 'NA', 
                              'PPFD' : 'mol_m2', 
                              'VPD' : 'kPa', 
                              'CPA' : 'NA', 
                              'Tair' : 'degC', 
                              'CO2' : 'ppm', 
                              'Patm' : 'Pa',
                              'elv' : 'm',
                              'Iabs' : 'mol_m2', 
                              'ca' : 'Pa',
                              'Gs' : 'Pa', 
                              'D' : 'Pa', 
                              'K' : 'Pa', 
                              'ns' : 'NA', 
                              'fa' : 'NA',
                              'GPP_hat' : 'mol_m2'}
        #
        # Define station value header line:
        header_vals = ['Timestamp', 'GPP', 'GPP_err', 'fPAR', 'PPFD',  
                       'VPD', 'CPA', 'Tair', 'CO2', 'Patm', 'elv', 'Iabs', 
                       'ca', 'Gs', 'D', 'K', 'ns', 'fa', 'GPP_hat']
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
        #
        # Dictionary of stations' light-use efficiency model variables
        # * this is for modeling
        self.st_lue_vars = {}
        self.lue_var_units = {'GPP' : 'mol_m2',
                              'GPP_err' : 'mol_m2',
                              'Iabs' : 'mol_m2', 
                              'ca' : 'Pa',
                              'Gs' : 'Pa', 
                              'D' : 'Pa', 
                              'K' : 'Pa', 
                              'ns' : 'NA', 
                              'fa' : 'NA'}
        #
        # Dictionary of stations' light-use efficiency model fit/fitness:
        # * this is for model results
        self.station_lue = {}
        #
        # Define standard viscosity of water, Pa s
        self.n25 = self.viscosity_h2o(self.kTo, self.kPo)
    #
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
                  - viscosity_h2o
                  - kPo
        """
        # Check for missing alpha (due to new STASH code):
        if alpha is None:
            alpha = -9999.0
        #
        # Calculate lue variables & place into tuple:
        iabs = fpar*ppfd                    # mol/m2, abs. PPFD
        ca = (1.e-6)*co2*patm               # Pa, atms. CO2
        gs = self.calc_gstar(tair, patm)    # Pa, photores. comp. point
        d = (1e3)*vpd                       # Pa, vapor pressure deficit
        k = self.calc_k(tair, patm)         # Pa, Michaelis-Menten coef.
        ns = self.viscosity_h2o(tair, patm) # Pa s, viscosity
        ns /= self.n25                      # unitless, water viscosity
        fa = (alpha/1.26)**(0.25)           # unitless, func. of alpha
        yhat = self.nxtgn(iabs, ca, gs, d, k, ns, fa)
        #
        # Place value parameters into tuples:
        val_params = (month, gpp, gpp_err, fpar, ppfd, vpd, alpha, tair, co2, 
                      patm, elv, iabs, ca, gs, d, k, ns, fa, yhat)
        var_params = (gpp, gpp_err, iabs, ca, gs, d, k, ns, fa)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Station Values
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize list if station key doesn't exist:
        if station not in self.station_vals.keys():
            self.station_vals[station] = numpy.array(
                val_params,
                dtype={'names' : ('Timestamp', 'GPP', 'GPP_err', 
                                  'fPAR', 'PPFD', 'VPD', 
                                  'CPA', 'Tair', 'CO2', 
                                  'Patm', 'elv', 'Iabs', 
                                  'ca', 'Gs', 'D', 
                                  'K', 'ns', 'fa', 
                                  'GPP_hat'),
                       'formats' : ('O', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4')},
                ndmin=1
                )
        else:
            # Add new parameters to list:
            temp_array = numpy.array(
                val_params,
                dtype={'names' : ('Timestamp', 'GPP', 'GPP_err', 
                                  'fPAR', 'PPFD', 'VPD', 
                                  'CPA', 'Tair', 'CO2', 
                                  'Patm', 'elv', 'Iabs', 
                                  'ca', 'Gs', 'D', 
                                  'K', 'ns', 'fa', 
                                  'GPP_hat'),
                       'formats' : ('O', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4', 'f4', 'f4',
                                    'f4')},
                ndmin=1
                )
            self.station_vals[station] = numpy.append(
                self.station_vals[station], 
                temp_array, 
                axis=0
                )
        #
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
                        'names' : ('GPP', 'GPP_err', 'Iabs', 
                                   'ca', 'Gs', 'D', 
                                   'K', 'ns', 'fa'),
                        'formats' : ('f4', 'f4', 'f4', 
                                     'f4', 'f4', 'f4', 
                                     'f4', 'f4', 'f4')
                        },
                    ndmin=1
                    )
            else:
                temp_array = numpy.array(
                    var_params,
                    dtype={
                        'names' : ('GPP', 'GPP_err', 'Iabs', 
                                   'ca', 'Gs', 'D', 
                                   'K', 'ns', 'fa'),
                        'formats' : ('f4', 'f4', 'f4', 
                                     'f4', 'f4', 'f4', 
                                     'f4', 'f4', 'f4')
                        },
                    ndmin=1
                    )
                self.st_lue_vars[station] = numpy.append(
                    self.st_lue_vars[station], 
                    temp_array, 
                    axis=0
                    )
    #
    def calc_gstar(self, tc, patm):
        """
        Name:     LUE.calc_gstar
        Input:    - float, air temperature, degrees C (tc)
                  - float, atmospheric pressure, Pa (patm)
        Output:   float, gamma-star, Pa (gs)
        Features: Returns the temperature and pressure dependent 
                  photorespiratory compensation point, Gamma star (Pascals), 
                  based on constants derived from Bernacchi et al. (2001).
        Ref:      Bernacchi et al. (2001), Improved temperature response 
                  functions for models of Rubisco-limited photosynthesis, 
                  Plant, Cell and Environment, 24, 253--259.
        """
        # Define constants
        gsc = 7.472      # empirical constant
        kco = 2.09476e5  # ppm, US Standard Atmosphere
        dha = 37830      # J/mol
        kR = 8.3145      # J/mol/K
        tk = tc + 273.15
        #
        gs = (5e-7)*kco*patm*numpy.exp(gsc - dha/(kR*tk))
        return gs
    #
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
        # Define constants
        kc25 = 39.97     # Pa, assuming 25 deg C & 98.716 kPa
        ko25 = (2.748e4) # Pa, assuming 25 deg C & 98.716 kPa
        dhac = 79430     # J/mol
        dhao = 36380     # J/mol
        kR = 8.3145      # J/mol/K
        kco = 2.09476e5  # ppm, US Standard Atmosphere
        #
        vc = kc25*numpy.exp(dhac*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        vo = ko25*numpy.exp(dhao*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        k = vc*(1 + kco*(1e-6)*patm/vo)
        return k
    #
    def calc_mgs_params(self, station):
        """
        Name:     LUE.calc_mgs_params
        Input:    str, station name (station)
        Ouput:    numpy.ndarray
        Features: Returns the mean-growing season parameters (i.e., ca, d, k,
                  gs, ns, and tair)
        """
        # NOTE: vars in st_lue_vars are above freezing (i.e., growing season)
        if station in self.st_lue_vars.keys():
            # Extract vectors of data:
            v_tair = self.station_vals[station]['Tair']    # deg C
            v_d = self.station_vals[station]['D']          # Pa
            v_ca = self.station_vals[station]['ca']        # Pa
            v_patm = self.station_vals[station]['Patm']    # Pa
            #
            # Extract constants:
            my_patm = v_patm[0]
            #
            # Find indexes where non-freezing:
            mgs_idx = numpy.where(v_tair > 0)
            #
            # Calculate mean-growing season values:
            mgs_tair = v_tair[mgs_idx].mean()
            mgs_d = v_d[mgs_idx].mean()
            mgs_ca = v_ca[mgs_idx].mean()
            #
            mgs_gs = self.calc_gstar(mgs_tair, my_patm)
            mgs_k = self.calc_k(mgs_tair, my_patm)
            mgs_ns = self.viscosity_h2o(mgs_tair, my_patm)
            mgs_ns /= self.n25
            #
            mgs_params = (station,
                          mgs_ca, mgs_d, mgs_k, 
                          mgs_gs, mgs_ns, mgs_tair)
        else:
            mgs_params = (station,
                          numpy.nan, numpy.nan, numpy.nan, 
                          numpy.nan, numpy.nan, numpy.nan)
        #
        mgs_return = numpy.array(
            mgs_params,
            dtype={'names' : ('station',
                              'mgs_ca', 'mgs_d', 'mgs_k', 
                              'mgs_gs', 'mgs_ns', 'mgs_tair'),
                   'formats' : ('S6',
                                'f4', 'f4', 'f4', 
                                'f4', 'f4', 'f4')},
            ndmin=1
            )
        return mgs_return
    #
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
        #
        # Calculate po, bar
        po = 5918.499
        po += 58.05267*tc
        po += -1.1253317*tc*tc
        po += (6.6123869e-3)*tc*tc*tc
        po += -(1.4661625e-5)*tc*tc*tc*tc
        #
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
        #
        # Convert pressure to bars (1 bar = 100000 Pa)
        pbar = (1e-5)*p
        #
        # Calculate the specific volume (cm^3 g^-1):
        v = vinf + my_lambda/(po + pbar)
        #
        # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
        rho = (1e3/v)
        #
        return rho
    #
    def next_gen_lue(self, x, beta):
        """
        Name:     LUE.next_gen_lue
        Input:    - numpy.ndarray, (k,M) shape of predictors (x)
                    > 'Iabs' : mol/m2, fAPARxPPFD
                    > 'ca' : Pa, atmospheric CO2
                    > 'Gs' : Pa, photores. comp. point
                    > 'D' : Pa, vapor pressure deficit
                    > 'K' : Pa, Michaelis-Menten coeff.
                    > 'ns' : mPa s, viscosity of water
                    > 'fa' : unitless, function of alpha
                  - float, beta parameter (beta)
        Output:   numpy.ndarray (gpp)
        Features: Returns array of GPP based on the next-generation light and
                  water use efficiency model.
        Depends:  - kc
                  - kphio
                  
        @TODO: eliminate bad points instead of fixing them
        """
        # Define variable substitutes:
        vdcg = x['ca'] - x['Gs']
        vacg = x['ca'] + 2.*x['Gs']
        vbkg = beta*(x['K'] + x['Gs'])
        #
        # Check for negatives in sqrt:
        vbkg_temp = [1e-6 + 0.*i if i<0 else i for i in vbkg]
        vbkg = numpy.array(vbkg_temp)
        #
        vsr = numpy.sqrt(1.6*x['ns']*x['D']/(vbkg))
        #
        # Based on the m' formulation (see Regressing_LUE.pdf)
        m = vdcg/(vacg + 3.*x['Gs']*vsr)
        temp_part = numpy.power(m, 2.) - self.kc**(2./3.)*numpy.power(m, 4./3.)
        # 
        # Check for negatives in temp part:
        temp_array = [1e-6 + 0.*i if i<0 else i for i in temp_part]
        temp_part = numpy.array(temp_array)
        #
        gpp = self.kphio*x['Iabs']*x['fa']*numpy.sqrt(temp_part)
        #
        return gpp
    #
    def nxtgn(self, iabs, ca, gs, d, k, ns, fa):
        """
        Name:     LUE.nxtgn
        Input:    - float, 'Iabs' : mol/m2, fAPARxPPFD
                  - float, 'ca' : Pa, atmospheric CO2
                  - float, 'Gs' : Pa, photores. comp. point
                  - float, 'D' : Pa, vapor pressure deficit
                  - float, 'K' : Pa, Michaelis-Menten coeff.
                  - float, 'ns' : mPa s, viscosity of water
                  - float, 'fa' : unitless, function of alpha
        Output:   float, estimate of GPP (gpp)
        Features: Returns an estimate of GPP based on the next-generation light 
                  and water use efficiency model.
        Depends:  - kc
                  - kphio
                  - kbeta
        """
        # Define default GPP return value:
        gpp = numpy.nan
        #
        # Define variable substitutes:
        vdcg = ca - gs
        vacg = ca + 2.*gs
        vbkg = self.kbeta*(k + gs)
        #
        # Check for negatives:
        if vbkg > 0:
            vsr = numpy.sqrt(1.6*ns*d/(vbkg))
            #
            # Based on the m' formulation (see Regressing_LUE.pdf)
            m = vdcg/(vacg + 3.*gs*vsr)
            mpi = m**2 - self.kc**(2./3.)*(m**(4./3.))
            # 
            # Check for negatives:
            if mpi > 0:
                mp = numpy.sqrt(mpi)
                gpp = self.kphio*iabs*fa*mp
        #
        return gpp
    #
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
        #
        # Get the density of water, kg/m^3
        rho = self.density_h2o(tc, p)
        #
        # Calculate dimensionless parameters:
        tbar = (tc + 273.15)/tk_ast
        tbarx = tbar**(0.5)
        tbar2 = tbar**2
        tbar3 = tbar**3
        rbar = rho/rho_ast
        #
        # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
        mu0 = 1.67752 
        mu0 += 2.20462/tbar 
        mu0 += 0.6366564/tbar2 
        mu0 += -0.241605/tbar3
        mu0 = 1e2*tbarx/mu0
        #
        # Create Table 3, Huber et al. (2009):
        hj0 = (0.520094, 0.0850895, -1.08374, -0.289555, 0., 0.)
        hj1 = (0.222531, 0.999115, 1.88797, 1.26613, 0., 0.120573)
        hj2 = (-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.)
        hj3 = (0.161913,  0.257399, 0., 0., 0., 0.)
        hj4 = (-0.0325372, 0., 0., 0.0698452, 0., 0.)
        hj5 = (0., 0., 0., 0., 0.00872102, 0.)
        hj6 = (0., 0., 0., -0.00435673, 0., -0.000593264)
        h = hj0 + hj1 + hj2 + hj3 + hj4 + hj5 + hj6
        h_array = numpy.reshape(numpy.array(h), (7,6))
        #
        # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
        mu1 = 0.
        ctbar = (1./tbar) - 1.
        for i in xrange(6):
            coef1 = numpy.power(ctbar, i)
            coef2 = 0.
            for j in xrange(7):
                coef2 += h_array[j][i]*numpy.power((rbar - 1.), j)
            mu1 += coef1*coef2
        mu1 = numpy.exp(rbar*mu1)
        #
        # Calculate mu_bar (Eq. 2, Huber et al., 2009)
        #   assumes mu2 = 1
        mu_bar = mu0*mu1
        #
        # Calculate mu (Eq. 1, Huber et al., 2009)
        mu = mu_bar*mu_ast    # Pa s
        #
        return mu
    #
    def write_out_val(self, station, out_file):
        """
        Name:     LUE.write_out_val
        Input:    - string, station name (station)
                  - string, output file (out_file)
        Output:   None.
        Features: Writes to file the monthly values associated with the light 
                  use efficiency equation for a given station
        """
        # Create file if it doesn't exist:
        if not os.path.isfile(out_file):
            try:
                f = open(out_file, 'w')
            except IOError:
                print "Cannot write to file:", out_file
            else:
                f.write(self.value_header)
                f.close()
        #
        # Print if station has data:
        if station in self.station_vals.keys():
            for t in self.station_vals[station]:
                try:
                    f = open(out_file, 'a')
                except IOError:
                    print "Cannot append to file:", out_file
                else:
                    f.write(("%s,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f,%f,"
                             "%f,%f,%f,%f\n") % tuple(t))
                    f.close()
    #
    def write_out_lue(self, out_file):
        """
        Name:     LUE.write_out_lue
        Input:    string, output file (out_file)
        Output:   None.
        Features: Writes to file the calculated light use efficiency for all 
                  stations
                  
        @TODO: update for basic and next-gen models?
        """
        # Create file if it doesn't exist:
        if not os.path.isfile(out_file):
            temp_string = '*_max,*_min,*_ave,*_std,*_skw,*_krt'
            lue_head = ("Station,r_sq,phio_est,phio_opt,beta_est,beta_opt," +
                        "phio_err,beta_err,phio_t,beta_t,phio_p,beta_p," +
                        temp_string.replace('*', 'GPP') +
                        "," +
                        temp_string.replace('*', 'Iabs') +
                        "," +
                        temp_string.replace('*', 'ca') +
                        "," +
                        temp_string.replace('*', 'Gs') +
                        "," +
                        temp_string.replace('*', 'D') +
                        "," +
                        temp_string.replace('*', 'K') +
                        "," +
                        temp_string.replace('*', 'eta') +
                        "\n")
            try:
                f = open(out_file, 'w')
            except IOError:
                print "Cannot write to file:", out_file
            else:
                f.write(lue_head)
                f.close()
        #
        # Print each station if it has data:
        for station in numpy.sort(self.station_lue.keys()):
            t = (station,) + self.station_lue[station]
            try:
                f = open(out_file, 'a')
            except IOError:
                print "Cannot append to file:", out_file
            else:
                f.write("%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f,"
                        "%f,%f,%f,%f,%f,%f\n" % t)
                f.close()
    #
    # ~~~~~~~~~~~~~~~~~ NOT USED FUNCTIONS ~~~~~~~~~~~~~~~~~~
    #
    def lue_resid(self, p, x, y):
        """
        Name:     LUE.lue_reside
        Inputs:   - float, model estimate (p)
                  - numpy.ndarray, model predictors (x)
                  - numpy.ndarray, model observations (y)
        Features: Returns the weighted residuals of the light-use efficiency 
                  model for a given model estimate
        Depends:  next_gen_lue
        
        NOT USED
        """
        sigma = x['GPP_err']
        y_hat = self.next_gen_lue(x, p)
        resid = (y_hat - y)/sigma
        return resid
    #
    def lue_jacob(self, p, x, y):
        """
        Name:     LUE.lue_jacob
        Inputs:   - float, model estimate (p)
                  - numpy.ndarray, model predictors (x)
                  - numpy.ndarray, model observations (y)
        Outputs:  numpy.ndarray (jacob)
        Features: Returns the Jacobian (first partial derivative) for the 
                  light-use efficiency equation (w.r.t. beta)
        Depends:  - kc
                  - kphio
        
        NOT USED
        """
        # Variable substitutes:
        kc23 = self.kc**(2./3.)
        xh = numpy.sqrt(1.6)
        b2 = p**2
        ced = x['ns']*x['D']
        v1 = xh*x['Gs']*x['ns']*x['D']
        v2 = x['ca'] - x['Gs']
        v22 = v2**2
        v243 = v2**(4./3.)
        v3 = x['K'] + x['Gs']
        bv3 = p*v3
        v4 = ced/bv3
        xv4 = numpy.sqrt(v4)
        v5 = 3.*xh*x['Gs']*xv4 + 2.*x['Gs'] + x['ca']
        v52 = v5**2
        v53 = v5**3
        v543 = v5**(4./3.)
        v573 = v5**(7./3.)
        s1 = self.kphio*x['fa']*x['Iabs']
        f1 = v3*xv4
        f2 = v22/v52
        f3 = v243/v543
        f4 = kc23*v1*v243
        f5 = v1*v22
        f6 = f2 - kc23*f3
        #
        # Check for negatives in f6:
        temp_list = [1e-6 + 0.*i if i < 0 else i for i in f6]
        f6 = numpy.array(temp_list)
        #
        # First partial derivative
        n1 = 3.*f5
        n2 = 2.*f4
        d1 = b2*f1*v53
        d2 = b2*f1*v573
        d3 = 2.*numpy.sqrt(f6)
        jacob = s1*(n1/d1 - n2/d2)
        jacob /= d3
        #
        return(jacob)
    #
    def calc_statistics(self, my_array):
        """
        Name:     LUE.calc_statistics
        Input:    numpy.ndarray (my_array)
        Output:   tuple, statistical quantities
                  - float, max value (max_val)
                  - float, min value (min_val)
                  - float, mean value (ave_val)
                  - float, standard deviation (std_val)
                  - float, skewness (skew_val)
                  - float, kurtosis (kurt_val)
        Features: Returns the basic/advanced statistics for an array of values
        
        NOT USED
        """
        # Make sure my_array is a numpy array:
        if not isinstance(my_array, numpy.ndarray):
            my_array = numpy.array(my_array)
        #
        # Initialize return array:
        my_stats = numpy.array(tuple([-9999., -9999., -9999., 
                                      -9999., -9999., -9999.]),
                               dtype={'names' : ('max', 'min', 'ave', 
                                                 'std', 'skw', 'krt'),
                                      'formats' : ('f4', 'f4', 'f4', 
                                                   'f4', 'f4', 'f4')},
                               ndmin=1)
        #
        # Make sure my_array is not empty or crashes on skew/kurt:
        if my_array.any() and len(my_array) > 1:
            # Max, min, mean, st dev, skew, kurtosis (offset from normal)
            max_val = my_array.max()
            min_val = my_array.min()
            ave_val = my_array.mean()
            std_val = my_array.std()
            #
            # Address divide by zero issues:
            if std_val == 0:
                std_val = 1e-4
            #
            skew_val = (
                sum((my_array - ave_val)**3)/
                ((len(my_array) - 1)*std_val**3)
                )
            kurt_val = (
                sum((my_array - ave_val)**4)/
                ((len(my_array) - 1)*std_val**4) - 3
                )
        else:
            # Maintain initial quantity values:
            max_val = -9999.
            min_val = -9999.
            ave_val = -9999.
            std_val = -9999.
            skew_val = -9999.
            kurt_val = -9999.
            #
        my_stats[0] = (max_val, min_val, ave_val, std_val, skew_val, kurt_val)
        return my_stats
    #
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
                  - float, predicted beta from precise expression (beta_p2)
        Features: Returns an estimate for beta based on the Wang Han equation
        
        NOT USED
        """
        if not numpy.isfinite(my_z):
            my_z = 0.
        #
        whe = numpy.exp(
            1.19 
            + 0.0545*(my_t - 25.)        # T in deg C
            - 0.5*numpy.log(1e-3*my_d)   # D in kPa
            - 0.0815*(1e-3*my_z)         # z in km
        )
        chi = whe/(1. + whe)
        #
        beta_p1 = 1.6*my_ns*my_d*(chi**2)
        beta_p1 /= (1. - chi)**2
        beta_p1 /= my_k
        #
        beta_p2 = 1.6*my_ns*my_d
        beta_p2 /= (my_k + my_gs)
        beta_p2 *= (chi*my_ca - my_gs)**2
        beta_p2 /= (my_ca**2)
        beta_p2 /= (chi - 1.)**2
        #
        return (beta_p1, beta_p2)
    #
    

###############################################################################
## FUNCTIONS:
###############################################################################
def get_elv(meta_dir, st_name):
    """
    Name:     get_elv
    Input:    - str, metadata file w/ path (meta_dir)
              - str, station name (st_name)
    Output:   float, elevation, m (my_elv)
    Features: Returns the station elevation from GePiSaT database metadata file
    """
    my_files = glob.glob(meta_dir + '*Met-Data*')
    if my_files:
        my_file = my_files[0]
        my_data = numpy.loadtxt(fname=my_file,
                                dtype={'names' : ('stationid', 'ele'),
                                       'formats' : ('S6', 'f4')},
                                delimiter=',',
                                skiprows=1,
                                usecols=(4,8),
                                converters={4 : numpy.str,
                                            8 : numpy.float})
        my_idx = numpy.where(my_data['stationid'] == st_name)[0]
        if my_idx:
            my_elv = my_data['ele'][my_idx[0]]
            if my_elv == -9999:
                my_elv = numpy.nan
        else:
            my_elv = numpy.nan
    else:
        my_elv = numpy.nan
    #
    return my_elv

def get_veg_type(meta_dir, st_name):
    """
    Name:     get_veg_type
    Inputs:   - str, meta data file w/ path (meta_dir)
              - str, station name (st_name)
    Output:   str, vegation type, short name 
    Features: Returns the vegation type for a given station
    """
    my_files = glob.glob(meta_dir + '*Met-Data*')
    if my_files:
        my_file = my_files[0]
        my_data = numpy.loadtxt(fname=my_file,
                                dtype={'names' : ('stationid', 'classid'),
                                       'formats' : ('S6', 'S3')},
                                delimiter=',',
                                skiprows=1,
                                usecols=(4,9),
                                converters={4 : numpy.str,
                                            9 : numpy.str})
        my_idx = numpy.where(my_data['stationid'] == st_name)
        try:
            my_veg = my_data['classid'][my_idx][0]
        except:
            my_veg = ''
    else:
        my_veg = ''
    #
    return my_veg

def get_lue(lue_class, station):
    """
    Name:     get_lue
    Input:    - LUE class (lue_class)
              - string, station name (station)
    Output:   tuple
              - numpy.ndarray, Iabs (x_data)
              - numpy.ndarray, GPP (y_data)
              - float, fitted phi_o (st_phio)
              - float, fitted beta (st_beta)
              - float, R-squared (st_rsqr)
    Features: Fits the next generation LUE model to the monthly flux tower 
              data and returns the model fit
    Depends:  - next_gen_lue
              - calc_gstar
              - calc_k
              - viscosity_h2o
    """
    # Initialize fitness parameters:
    st_rsqr = -9999.              # model coef. of determination
    st_beta = -9999.              # beta parameter
    st_phio = lue_class.kphio     # intrinsic quantum efficiency parameter
    #
    if station in lue_class.station_vals.keys():
        # Get predictors and observations:
        x_data = numpy.copy(lue_class.st_lue_vars[station])
        y_data = x_data['GPP']
        #
        # Estimate beta:
        est_beta = lue_class.beta_estimate(station)
        #
        # Curve fit:
        plsq, cov, infodict, mesg, ier = leastsq(lue_class.lue_resid, 
                                                 est_beta,
                                                 args=(x_data, y_data), 
                                                 Dfun=lue_class.lue_jacob, 
                                                 full_output=True)
        #
        beta_fit = plsq[0]
        y_fit = lue_class.next_gen_lue(x_data, beta_fit)
        #
        # Calculate fitness:
        m = len(y_data)
        ssxx = numpy.power(y_fit, 2).sum() - m*(y_fit.mean()**2)
        ssyy = numpy.power(y_data, 2).sum() - m*(y_data.mean()**2)
        ssxy = (y_fit*y_data).sum() - m*y_fit.mean()*y_data.mean()
        rsqr = ssxy**2/(ssxx*ssyy)
        #
        st_beta = beta_fit
        st_rsqr = rsqr
    #
    return (x_data, y_data, st_phio, st_beta, st_rsqr)

def make_plot(my_fit, my_obs, my_text):
    """
    Name:     make_plot
    Input:    - numpy.ndarray, modelled GPP (my_fit)
              - numpy.ndarray, observed GPP (my_obs)
              - str, plot text (my_text)
    Output:   None
    Features: Creates a plot of predicted versus observed GPP
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plot_max = numpy.concatenate((my_fit, my_obs)).max()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(my_obs, my_fit, 'ro', label='Fitted') 
    ax1.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax1.set_ylabel('Modeled GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.set_xlabel('Observed GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                ncol=3, mode="expand", borderaxespad=0., fontsize=14)
    ax1.text(0.05, 0.95, my_text, transform=ax1.transAxes, fontsize=14, 
            verticalalignment='top', bbox=props)
    plt.show()

def wang_han_eq(my_d, my_k, my_n, my_t, my_z, my_gs, my_ca):
    """
    Name:     wang_han_eq
    Inputs:   - float, vapor pressure deficit, Pa (my_d)
              - float, Michaelis-Menten coef, Pa (my_k)
              - float, viscosity, unitless (my_n)
              - float, air temperature, deg C (my_t)
              - float, elevation, m (my_z)
              - float, photorespiratory comp. point, Pa (my_gs)
              - float, atmospheric CO2 concentration, Pa (my_ca)
    Output:   tuple
              - float, predicted chi from Wang Han w/o z (chi2)
              - float, predicted chi from Wang Han (chi1)
              - float, predicted beta fro. simple expression (beta_p1)
              - float, predicted beta from precise equation (beta_p2)
    Features: Returns the estimate for beta based on the Wang Han equation
    """
    # Wang Han's equation (new method):
    whe = numpy.exp(
        1.19 
        + 0.0545*(my_t - 25.)         # T in deg C
        - 0.5*numpy.log(1e-3*my_d)    # D in kPa
        - 0.0815*(1e-3*my_z)          # z in km
    )
    chi1 = whe/(1. + whe)
    #
    whe2 = numpy.exp(
        1.19 
        + 0.0545*(my_t - 25.)         # T in deg C
        - 0.5*numpy.log(1e-3*my_d)    # D in kPa
    )
    chi2 = whe2/(1. + whe2)
    #
    # VPD-based chi estimate (old method):
    #chi2 = 0.5 - (0.5 - 0.9)*(2.5e3 - my_d)/(2.5e3)
    #
    # Simple beta estimate:
    beta_p1 = 1.6*my_n*my_d*(chi1**2)
    beta_p1 /= (1. - chi1)**2
    beta_p1 /= my_k
    #
    # Precise beta estimate:
    beta_p2 = 1.6*my_n*my_d
    beta_p2 /= (my_k + my_gs)
    beta_p2 *= (chi1*my_ca - my_gs)**2
    beta_p2 /= (my_ca**2)
    beta_p2 /= (chi1 - 1.)**2
    #
    return (chi2, chi1, beta_p1, beta_p2)

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

def make_bplot(my_data, my_names, x_lab):
    """
    Name:     make_bplot
    Input:    - list (my_data)
              - str, xtick label names (my_names)
              - str, x axis label (x_lab)
    Features: Boxplot of data
    """
    fig, ax1 = plt.subplots(figsize=(6,12))
    fig.canvas.set_window_title('A Boxplot Example')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    #
    bp = plt.boxplot(my_data, notch=0, sym='+', vert=0, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    #
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    ax1.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                   alpha=0.5)
    ax1.set_axisbelow(True)
    ax1.set_xlabel(x_lab, fontsize=16)
    ytickNames = plt.setp(ax1, yticklabels=my_names)
    plt.setp(ytickNames, rotation=0, fontsize=6) # was fontsize=8
    plt.show()

def model_fitness(fit, obs):
    """
    Name:     model_fitness
    Input:    - numpy.ndarray, fitted data (fit)
              - numpy.ndarray, observed data (obs)
    Output:   - float, coefficient of determination, (rsqr)
              - float, correlation coefficient (r)
              - float, index of agreement (ioa)
    Features: Returns model fitness parameters (rsqr, r, ioa)
    """
    # 1. Calculate the R-square
    # VERSION 1: ratio of regression sum of squares (SSR) to total sum of 
    #            squares (SST), where SSR = SST - SSE (error sum of squares):
    sst = obs - obs.mean()
    sst = numpy.power(sst, 2.0)
    sst = sst.sum()
    sse = obs - fit
    sse = numpy.power(sse, 2.0)
    sse = sse.sum()
    rsqr = (sst - sse)/sst
    #
    # VERSION 2: correlation coefficient, r
    n = float(len(obs))
    ssxy = (obs*fit).sum() - obs.sum()*fit.sum()/n
    ssx = numpy.power(fit, 2.0).sum() - numpy.power(fit.sum(), 2.0)/n
    sst = numpy.power(obs, 2.0).sum() - numpy.power(obs.sum(), 2.0)/n
    r = ssxy/numpy.power(ssx*sst, 0.5)
    #
    # 2. Calculate the index of agreement
    mae = numpy.abs(fit - obs)
    smae = mae.sum()
    mao = numpy.abs(obs - obs.mean())
    smao = mao.sum()
    tsmo = 2.0*smao
    #
    if smae <= tsmo:
        ioa = 1.0 - smae/tsmo
    else:
        ioa = tsmo/smae - 1.0
    #
    return (rsqr, r, ioa)

def make_one_plot(my_obs, my_fit, my_txt):
    """
    Name:     make_one_plot
    Input:    - numpy.ndarray, observed GPP (my_obs)
              - numpy.ndarray, modelled GPP (my_fit)
              - str, plot text for model (my_txt)
    Output:   None
    Features: Creates a plot of predicted versus observed GPP
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plot_max = numpy.concatenate((my_fit, my_obs)).max()
    #
    fig = plt.figure(figsize=(8,8), dpi=180)
    #
    ax1 = fig.add_subplot(111)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(my_obs, my_fit, 'ro', label='Basic LUE formula') 
    ax1.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax1.set_ylabel('Modeled GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.set_xlabel('Observed GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)
    ax1.text(0.05, 0.95, my_txt, transform=ax1.transAxes, fontsize=14, 
             verticalalignment='top', bbox=props)
    #
    plt.show()

def make_two_plots(my_obs, my_fit1, my_fit2, my_txt1, my_txt2, v=1):
    """
    Name:     make_two_plots
    Input:    - numpy.ndarray, observed GPP (my_obs)
              - numpy.ndarray, modelled 1 GPP (my_fit1)
              - numpy.ndarray, modelled 2 GPP (my_fit2)
              - str, plot text for model 1 (my_txt1)
              - str, plot text for model 2 (my_txt2)
              - int, version number
    Output:   None
    Features: Creates two plots of predicted versus observed GPP
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plot_max = numpy.concatenate((my_fit1, my_fit2, my_obs)).max()
    #
    fig = plt.figure(figsize=(14,8), dpi=180)
    #
    ax1 = fig.add_subplot(121)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    if v == 1:
        ax1.plot(my_obs, my_fit1, 'ro', label='Simple $\\beta$ formula') 
    elif v == 2:
        ax1.plot(my_obs, my_fit1, 'ro', label='Basic LUE model') 
    elif v == 3:
        ax1.plot(my_obs, my_fit1, 'ro', label='Seasonal')
    ax1.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax1.set_ylabel('Modeled GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.set_xlabel('Observed GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)
    ax1.text(0.05, 0.95, my_txt1, transform=ax1.transAxes, fontsize=14, 
             verticalalignment='top', bbox=props)
    #
    ax2 = fig.add_subplot(122)
    plt.setp(ax2.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=14)
    if v == 1:
        ax2.plot(my_obs, my_fit2, 'ro', label='Precise $\\beta$ formula') 
    elif v == 2:
        ax2.plot(my_obs, my_fit2, 'ro', label='Next-Gen LUE model') 
    elif v == 3:
        ax2.plot(my_obs, my_fit2, 'ro', label='Mean Growing Season')
    ax2.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax2.set_ylabel('Modeled GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax2.set_xlabel('Observed GPP, mol CO$_2$ m$^{-2}$', fontsize=16)
    ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)
    ax2.text(0.05, 0.95, my_txt2, transform=ax2.transAxes, fontsize=14, 
             verticalalignment='top', bbox=props)
    #
    plt.show()

def elv2pres(z):
    """
    Name:     elv2pres
    Input:    float, elevation above sea level (z), m
    Output:   float, atmospheric pressure, Pa
    Features: Calculates atm. pressure for a given elevation
    Ref:      Allen et al. (1998)
    """
    # Def constants:
    kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
    kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kTo = 298.15   # base temperature, K (Prentice, unpublished)
    kG = 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    #
    p = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
    return p

def calculate_vpd(tmp, vap):
    """
    Name:     calculate_vpd
    Input:    - float, mean monthly daily air temp, deg C (tmp)
              - float, mean monthly vapor pressure, hPa (vap)
    Output:   float, mean monthly vapor pressure deficit, kPa (vpd)
    Features: Returns mean monthly vapor pressure deficit
    Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure 
              Calculation Methods, in Evaporation and Evapotranspiration: 
              Measurements and Estimations, Springer, London.
                vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
                where:
                    tc = average daily air temperature, deg C
                    ea = actual vapor pressure, kPa
    """
    vpd = (0.611*numpy.exp((17.27*tmp)/(tmp + 237.3)) - 0.10*vap)
    return vpd

def plot_two_obs(my_obs1, my_fit1, my_obs2, my_fit2, my_txt1, my_txt2):
    """
    Name:     plot_two_obs
    Input:    - numpy.ndarray, observed 1 GPP (my_obs1)
              - numpy.ndarray, observed 2 GPP (my_obs2)
              - numpy.ndarray, modelled 1 GPP (my_fit1)
              - numpy.ndarray, modelled 2 GPP (my_fit2)
              - str, plot text for model 1 (my_txt1)
              - str, plot text for model 2 (my_txt2)
    Output:   None
    Features: Creates two plots of predicted versus observed GPP
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plot_max = numpy.concatenate((my_fit1, my_fit2, my_obs1, my_obs2)).max()
    #
    fig = plt.figure(figsize=(14,8), dpi=180)
    #
    ax1 = fig.add_subplot(121)
    plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(my_obs1, my_fit1, 'ro') 
    ax1.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax1.set_ylabel('Modeled GPP, mol C m$^{-2}$ mo$^{-1}$', fontsize=16)
    ax1.set_xlabel('Observed GPP, mol C m$^{-2}$ mo$^{-1}$', fontsize=16)
    ax1.text(0.05, 0.95, my_txt1, transform=ax1.transAxes, fontsize=14, 
             verticalalignment='top', bbox=props)
    #
    ax2 = fig.add_subplot(122)
    plt.setp(ax2.get_xticklabels(), rotation=0, fontsize=14)
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=14)
    ax2.plot(my_obs2, my_fit2, 'ro') 
    ax2.plot([0., plot_max], [0., plot_max], '--k', label='1:1 Line')
    ax2.set_ylabel('Modeled GPP, mol C m$^{-2}$ mo$^{-1}$', fontsize=16)
    ax2.set_xlabel('Observed GPP, mol C m$^{-2}$ mo$^{-1}$', fontsize=16)
    ax2.text(0.05, 0.95, my_txt2, transform=ax2.transAxes, fontsize=14, 
             verticalalignment='top', bbox=props)
    #
    plt.show()

###############################################################################
## MAIN PROGRAM:
###############################################################################
mac = False
if mac:
    my_dir = '/Users/twdavis/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    met_dir = '/Users/twdavis/Dropbox/Work/Imperial/flux/data/psql-data/flux/'
    out_dir = '/Users/twdavis/Desktop/'
else:
    my_dir = '/home/user/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    met_dir = '/home/user/Dropbox/Work/Imperial/flux/data/psql-data/flux/'
    out_dir = '/home/user/Desktop/'

my_files = glob.glob(my_dir + 'GePiSaT_All_Data-26*')
my_data = numpy.loadtxt(
    fname=my_files[0],
    dtype={'names': ('station', 'timestamp', 'gpp', 'gpp_err', 'fapar',
                     'ppfd', 'vpd', 'alpha', 'tair', 'co2', 'patm'),
           'formats' : ('S6', 'O', 'f4', 'f4', 'f4',
                        'f4', 'f4', 'f4', 'f4', 'f4', 'f4')},
    delimiter=',',
    skiprows=1,
    converters={
        0 : numpy.str,
        1 : lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'),
        2 : numpy.float,                       # mol/m2
        3 : numpy.float,                       # mol/m2
        4 : numpy.float,                       # unitless
        5 : numpy.float,                       # mol/m2
        6 : numpy.float,                       # kPa
        7 : numpy.float,                       # unitless
        8 : numpy.float,                       # deg C
        9 : numpy.float,                       # ppm
        10 : numpy.float                       # Pa
    }
)


# Create a LUE class as in GePiSaT model
my_lue = LUE()
for line in my_data:
    (st_name, st_time, st_gpp, st_gpp_err, st_fpar, 
     st_ppfd, st_vpd, st_cpa, st_tair, st_co2, st_patm) = line
    st_elv = get_elv(met_dir, st_name)
    #
    my_lue.add_station_val(st_name, 
                           st_time.date(), 
                           st_gpp,
                           st_gpp_err, 
                           st_fpar, 
                           st_ppfd,
                           st_vpd,
                           st_cpa,
                           st_tair,
                           st_co2,
                           st_patm,
                           st_elv)

# List of each station:
all_stations = numpy.sort(numpy.array(list(set(my_data['station']))))

# Test calc_mgs_params function:
# @TODO: get mean growing season values 
station = all_stations[0]
my_lue.calc_mgs_params(station)

# Test write_out_val function:
station = 'AT-Neu'
out_file = out_dir + station + "_lue_vals.txt"
my_lue.write_out_val(station, out_file)


##############################################################################
#                                                                            #
#                                WORKSPACE                                   #
#                                                                            #
##############################################################################
# Plot crops

# Dictionary of stations, years, veg and irrigation:
my_crop_dict = {
    'BE-Lon' : {2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : False, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'CH-Oe2' : {2005 : {'C4' : False, 'Irr' : False}},
    'DE-Geb' : {2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : False, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'DE-Kli' : {2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : False, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'DK-Fou' : {2005 : {'C4' : False, 'Irr' : False}},
    'DK-Ris' : {2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : False, 'Irr' : False}},
    'ES-ES2' : {2004 : {'C4' : False, 'Irr' : True},
                2005 : {'C4' : False, 'Irr' : True},
                2006 : {'C4' : False, 'Irr' : True}},
    'FR-Gri' : {2005 : {'C4' : True, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'IE-Ca1' : {2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : False, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'IT-BCi' : {2004 : {'C4' : True, 'Irr' : True},
                2005 : {'C4' : True, 'Irr' : True},
                2006 : {'C4' : True, 'Irr' : True}},
    'IT-Cas' : {2006 : {'C4' : True, 'Irr' : True}},
    'NL-Lan' : {2005 : {'C4' : True, 'Irr' : False},
                2006 : {'C4' : True, 'Irr' : False}},
    'NL-Lut' : {2006 : {'C4' : False, 'Irr' : False}},
    'NL-Mol' : {2005 : {'C4' : False, 'Irr' : False},
                2006 : {}},
    'UK-ESa' : {2003 : {},
                2004 : {},
                2005 : {}},
    'UK-Her' : {2006 : {'C4' : False, 'Irr' : False}},
    'US-ARM' : {2003 : {'C4' : False},
                2004 : {'C4' : False},
                2005 : {'C4' : True},
                2006 : {'C4' : False}},
    'US-Bo1' : {2002 : {'C4' : False, 'Irr' : False},
                2003 : {'C4' : True, 'Irr' : False},
                2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : True, 'Irr' : False},
                2006 : {'C4' : False, 'Irr' : False}},
    'US-Ne1' : {2002 : {'C4' : True, 'Irr' : True},
                2003 : {'C4' : True, 'Irr' : True},
                2004 : {'C4' : True, 'Irr' : True},
                2005 : {'C4' : True, 'Irr' : True}},
    'US-Ne2' : {2002 : {'C4' : False, 'Irr' : True},
                2003 : {'C4' : True, 'Irr' : True},
                2004 : {'C4' : False, 'Irr' : True},
                2005 : {'C4' : True, 'Irr' : True}},
    'US-Ne3' : {2002 : {'C4' : False, 'Irr' : False},
                2003 : {'C4' : True, 'Irr' : False},
                2004 : {'C4' : False, 'Irr' : False},
                2005 : {'C4' : True, 'Irr' : False}}
}

crop_stations = numpy.sort(my_crop_dict.keys())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# WETLAND REGRESSIONS
wet_stations = [
    'AU-Fog', 'CA-Mer', 'CZ-wet', 'FI-Kaa', 
    'PL-wet', 'SE-Deg', 'SE-Faj', 'UK-AMo', 
]

wetlands = {}
for station in wet_stations:
    xdata = numpy.copy(my_lue.st_lue_vars[station])
    yobs = xdata['GPP']
    x = yobs[:, numpy.newaxis]
    #
    # Mean Growing Season GPP
    iabs_data = xdata['Iabs']
    fa_data = xdata['fa']
    mgs_data = my_lue.calc_mgs_params(station)
    mgs_ymod = [my_lue.nxtgn(
        iabs_data[i], 
        mgs_data['mgs_ca'], 
        mgs_data['mgs_gs'], 
        mgs_data['mgs_d'],
        mgs_data['mgs_k'],
        mgs_data['mgs_ns'],
        fa_data[i]
        )[0] for i in xrange(len(iabs_data))]
    y = numpy.array(mgs_ymod)
    #
    fit_slope, fit_sse, fit_rank, fit_s = numpy.linalg.lstsq(x, y)
    wetlands[station] = (fit_slope[0], mgs_data['mgs_tair'][0])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Plot MGS versus seasonally-varying parameter estimates of GPP
# Presentation stations
# station = 'NL-Loo'  nlloo_obs = numpy.copy(yobs)
#                     nlloo_mod = numpy.copy(mgs_ymod)
#                     nlloo_txt = my_txt2
#
# station = 'BE-Vie'  bevie_obs = numpy.copy(yobs)
#                     bevie_mod = numpy.copy(mgs_ymod)
#                     bevie_txt = my_txt2
#
# station = 'DE-Hai'  dehai_obs = numpy.copy(yobs)
#                     dehai_mod = numpy.copy(mgs_ymod)
#                     dehai_txt = my_txt2
#
# station = 'CH-Oe1'  choe1_obs = numpy.copy(yobs)
#                     choe1_mod = numpy.copy(mgs_ymod)
#                     choe1_txt = my_txt2
#
# station = 'US-Bo1'  usbo1_obs = numpy.copy(yobs)
#                     usbo1_mod = numpy.copy(mgs_ymod)
#                     usbo1_txt = my_txt2
#
# station = 'SE-Deg'  sedeg_obs = numpy.copy(yobs)
#                     sedeg_mod = numpy.copy(mgs_ymod)
#                     sedeg_txt = my_txt2
#
# 93, 17, 10, 97, 0.20, 0.25
# plot_two_obs(nlloo_obs, nlloo_mod, bevie_obs, bevie_mod, nlloo_txt, bevie_txt)
# plot_two_obs(dehai_obs, dehai_mod, choe1_obs, choe1_mod, dehai_txt, choe1_txt)
# plot_two_obs(usbo1_obs, usbo1_mod, sedeg_obs, sedeg_mod, usbo1_txt, sedeg_txt)

fig_file = out_dir + 'GePiSaT_nxgn_GPP_seas-v-mgs.pdf'
pp = PdfPages(fig_file)

for station in all_stations:
    xdata = numpy.copy(my_lue.st_lue_vars[station])
    yobs = xdata['GPP']
    #
    # Seasonal GPP
    sea_ymod = my_lue.next_gen_lue(xdata, my_lue.kbeta)
    my_rsq1, my_r1, my_ioa1 = model_fitness(sea_ymod, yobs)
    #
    # Mean Growing Season GPP
    iabs_data = xdata['Iabs']
    fa_data = xdata['fa']
    mgs_data = my_lue.calc_mgs_params(station)
    mgs_ymod = [my_lue.nxtgn(
        iabs_data[i], 
        mgs_data['mgs_ca'], 
        mgs_data['mgs_gs'], 
        mgs_data['mgs_d'],
        mgs_data['mgs_k'],
        mgs_data['mgs_ns'],
        fa_data[i]
        )[0] for i in xrange(len(iabs_data))]
    mgs_ymod = numpy.array(mgs_ymod)
    my_rsq2, my_r2, my_ioa2 = model_fitness(mgs_ymod, yobs)
    #
    # Plot
    my_veg_type = get_veg_type(met_dir, station)
    my_txt1 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n"
               "$R^2=%0.3f$\n$r=%0.3f$\n"
               "$IA=%0.3f$") % (station, my_veg_type, my_rsq1,
                                my_r1, my_ioa1)
    my_txt2 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n"
               "$R^2=%0.3f$\n$r=%0.3f$\n"
               "$IA=%0.3f$") % (station, my_veg_type, my_rsq2,
                                my_r2, my_ioa2)
    #
    if len(yobs) > 2:
        make_two_plots(yobs, sea_ymod, mgs_ymod, my_txt1, my_txt2, 3)
        pp.savefig()
        plt.close()

d = pp.infodict()
d['Title'] = 'GePiSaT next-gen GPP plots of seasonal versus mean-growning season parameters'
d['Author'] = 'Tyler W. Davis'
pp.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Beni's data:
# fAPAR = 1.0 for all months
# alpha = 1.26 for all months
# co2 = 376 ppm for all months
# elv = 450 m
my_dir = '/Users/twdavis/Dropbox/Work/Imperial/collaborations/stocker/p-model_intercomparison/'
my_files = glob.glob(my_dir + '*2002*.txt')
my_mo_ppfd = numpy.loadtxt(
    fname = my_files[1],
    dtype={'names': ('month', 'ppfd'),
           'formats' : ('i4', 'f4')},
    delimiter=',',
    skiprows=1,
    converters={0 : numpy.int, 
                1 : numpy.float}
)
my_mo_tc = numpy.loadtxt(
    fname = my_files[3],
    dtype={'names': ('month', 'tc'),
           'formats' : ('i4', 'f4')},
    delimiter=',',
    skiprows=1,
    converters={0 : numpy.int, 
                1 : numpy.float}
)
my_mo_vap = numpy.loadtxt(
    fname = my_files[4],
    dtype={'names': ('month', 'vap'),
           'formats' : ('i4', 'f4')},
    delimiter=',',
    skiprows=1,
    converters={0 : numpy.int, 
                1 : numpy.float}
)

my_lue = LUE()
st_name = 'CH-Oe1'
st_gpp = 0.0
st_gpp_err = 0.0
st_fpar = 1.0
st_cpa = 1.26
st_co2 = 376.
st_elv = 450.
st_patm = elv2pres(st_elv)
for i in xrange(12):
    st_time = datetime.date(2002, (i + 1), 1)
    st_tair = my_mo_tc['tc'][i]
    st_vpd = calculate_vpd(st_tair, my_mo_vap['vap'][i])
    st_ppfd = my_mo_ppfd['ppfd'][i]
    #
    my_lue.add_station_val(st_name, 
                           st_time, 
                           st_gpp,
                           st_gpp_err, 
                           st_fpar, 
                           st_ppfd,
                           st_vpd,
                           st_cpa,
                           st_tair,
                           st_co2,
                           st_patm,
                           st_elv)

my_x = my_lue.st_lue_vars[st_name]
my_y = my_lue.next_gen_lue(my_x, 244.033)

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
ax1.plot(numpy.array([i for i in xrange(12)]), my_y, 'k-o') 
ax1.set_ylabel('GPP, mol C m$^{-2}$ mo$^{-1}$', fontsize=16)
ax1.set_xlabel('Month', fontsize=16)
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
ax1.plot(numpy.array([i for i in xrange(12)]), my_x['D'], 'k-o') 
ax1.set_ylabel('', fontsize=16)
ax1.set_xlabel('Month', fontsize=16)
plt.show()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Estimate beta from Wang Han equation:
out_file = out_dir + 'Wang_Han_beta_estimates.txt'
header = ('Timestamp,station,Tair_degC,D_kPa,K_Pa,'
          'ns,elv_km,Gs_Pa,ca_Pa,chi2,chi1,beta1,beta2\n')
writeout(out_file, header)
beta1_dict = {}
beta2_dict = {}
chi_dict = {}
monthly_dict = {}
for station in all_stations:
    z_data = numpy.copy(my_lue.station_vals[station])
    m = z_data.shape[0]
    for i in xrange(m):
        elv = z_data['elv'][i]                  # m
        tair = z_data['Tair'][i]                # deg C
        patm = z_data['Patm'][i]                # Pa
        d = 1e3*z_data['VPD'][i]                # Pa
        k = my_lue.calc_k(tair, patm)           # Pa
        ns = my_lue.viscosity_h2o(tair, patm)   # Pa s
        ns /= my_lue.n25                        # unitless
        gs = my_lue.calc_gstar(tair)            # Pa
        ca = z_data['CO2'][i]                   # ppm
        ca *= (1.e-6)*patm                      # Pa
        #
        cest2, cest1, best1, best2 = wang_han_eq(d, k, ns, tair, elv, gs, ca)
        #
        if numpy.isfinite(best1) and numpy.isfinite(best2):
            my_time = z_data['Timestamp'][i]
            my_month = my_time.month
            if my_month in monthly_dict.keys():
                monthly_dict[my_month] += (best2,)
            else:
                monthly_dict[my_month] = (best2,)
            #
            if station in chi_dict.keys():
                chi_dict[station] += (cest1,)
            else:
                chi_dict[station] = (cest1,)
            #
            if station in beta1_dict.keys():
                beta1_dict[station] += (best1,)
            else:
                beta1_dict[station] = (best1,)
            #
            if station in beta2_dict.keys():
                beta2_dict[station] += (best2,)
            else:
                beta2_dict[station] = (best2,)
            #
            out_line = '%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % (
                z_data['Timestamp'][i].date(), station, tair, d, k, ns, 
                elv, gs, ca, cest2, cest1, best1, best2)
            #
            try:
                OUT = open(out_file, 'a')
                OUT.write(out_line)
            except IOError:
                print "Error: cannot write to file: ", out_file
            else:
                OUT.close()

# Read estimates:
# Note: 0.5 < chi < 1.0
est_data = numpy.loadtxt(
    fname=out_file,
    dtype={
        'names': (
            'timestamp', 'station', 'tair', 'D', 'K',
            'ns', 'elv', 'Gs', 'ca', 'chi2', 'chi1', 'beta1', 'beta2'
        ),
        'formats' : (
            'O', 'S6', 'f4', 'f4', 'f4', 'f4', 'f4',
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4'
        )
    },
    delimiter=',',
    skiprows=1,
    converters={
        0 : lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'),
        1 : numpy.str,
        2 : numpy.float,                       # deg C
        3 : numpy.float,                       # kPa
        4 : numpy.float,                       # Pa
        5 : numpy.float,                       # unitless
        6 : numpy.float,                       # km
        7 : numpy.float,                       # Pa
        8 : numpy.float,                       # Pa
        9 : numpy.float,                       # unitless
        10 : numpy.float,                      # unitless
        11 : numpy.float,                      # unitless
        12 : numpy.float                       # unitless
    }
)



fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
ax1.plot(est_data['chi1'], est_data['chi2'], 'ko') 
ax1.set_ylabel('$\\chi_o$ without elevation', fontsize=16)
ax1.set_xlabel('$\\chi_o$ full equation', fontsize=16)
plt.show()


monthly_names = numpy.sort(monthly_dict.keys())
monthly_chi = []
monthly_best1 = []
monthly_best2 = []
for name in monthly_names:
    monthly_best2.append(monthly_dict[name])

month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
make_bplot(monthly_chi, month_names, '$\\chi_o$')
make_bplot(monthly_best1, month_names, 'Simple $\\beta$')
make_bplot(monthly_best2, month_names, 'Precise $\\beta$')

beta1_names = numpy.sort(beta1_dict.keys())
beta1_data = []
for name in beta1_names:
    beta1_data.append(beta1_dict[name])

beta2_names = numpy.sort(beta2_dict.keys())
beta2_data = []
for name in beta2_names:
    beta2_data.append(beta2_dict[name])

chi_names = numpy.sort(chi_dict.keys())
chi_data = []
for name in chi_names:
    chi_data.append(chi_dict[name])

make_bplot(beta1_data, beta1_names, 'Simple $\\beta$')
make_bplot(beta2_data, beta2_names, 'Precise $\\beta$')
make_bplot(chi_data, chi_names, '$\\chi_o$')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

##
## PLOT GPP BASED ON MEAN-GROWING SEASON BETA ESTIMATES FOR EACH SITE
##    depends: (1) model_fitness, (2) get_veg_type, (3) make_two_plots
##
fig_file = out_dir + 'GePiSaT_LUE_beta_mgs.pdf'
pp = PdfPages(fig_file)

for station in all_stations:
    x_data = numpy.copy(my_lue.st_lue_vars[station])
    y_data = x_data['GPP']
    #
    b_est1, b_est2 = my_lue.calc_mgs_beta(station)
    y_fit1 = my_lue.next_gen_lue(x_data, b_est1)
    y_fit2 = my_lue.next_gen_lue(x_data, b_est2)
    #
    my_rsq1a, my_rsq1b, my_ioa1 = model_fitness(y_fit1, y_data)
    my_rsq2a, my_rsq2b, my_ioa2 = model_fitness(y_fit2, y_data)
    #
    # @TODO: get vegetation type for station & add to plot
    my_veg_type = get_veg_type(met_dir, station)
    my_txt1 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n$\\beta=%0.2f$\n"
               "$R^2=%0.3f$\n$r=%0.3f$\n"
               "$IA=%0.3f$") % (station, my_veg_type, b_est1, my_rsq1a,
                                my_rsq1b, my_ioa1)
    my_txt2 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n$\\beta=%0.2f$\n"
               "$R^2=%0.3f$\n$r=%0.3f$\n"
               "$IA=%0.3f$") % (station, my_veg_type, b_est2, my_rsq2a,
                                my_rsq2b, my_ioa2)
    #
    if len(y_data) > 2:
        make_two_plots(y_data, y_fit1, y_fit2, my_txt1, my_txt2)
        pp.savefig()
        plt.close()

d = pp.infodict()
d['Title'] = 'GePiSaT GPP Plots with mean-growing season estimates of beta'
d['Author'] = 'Tyler W. Davis'
pp.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  
##
## PLOT GPP BASED ON SEASONAL BETA ESTIMATES FOR EACH SITE
##   depends: (1) get_veg_type, (2) model_fitness, (3) make_two_plots
##
fig_file = out_dir + 'GePiSaT_LUE_beta_seasonal.pdf'
pp = PdfPages(fig_file)

for station in all_stations:
    x_data = numpy.copy(my_lue.st_lue_vars[station])
    y_data = x_data['GPP']
    #
    my_veg_type = get_veg_type(met_dir, station)
    #
    b_est1 = x_data['beta1']
    b_ave1 = b_est1.mean()
    b_std1 = b_est1.std()
    y_fit1 = my_lue.next_gen_lue(x_data, b_est1)
    my_rsq1a, my_rsq1b, my_ioa1 = model_fitness(y_fit1, y_data)
    my_txt1 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n"
              "$\\beta=%0.2f\\pm %0.2f$\n"
              "$R^2=%0.3f$\n$r=%0.3f$\n"
              "$IA=%0.3f$") % (station, my_veg_type, b_ave1, b_std1, my_rsq1a, 
                               my_rsq1b, my_ioa1)
    #
    b_est2 = x_data['beta2']
    b_ave2 = b_est2.mean()
    b_std2 = b_est2.std()
    y_fit2 = my_lue.next_gen_lue(x_data, b_est2)
    my_rsq2a, my_rsq2b, my_ioa2 = model_fitness(y_fit2, y_data)
    my_txt2 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n"
              "$\\beta=%0.2f\\pm %0.2f$\n"
              "$R^2=%0.3f$\n$r=%0.3f$\n"
              "$IA=%0.3f$") % (station, my_veg_type, b_ave2, b_std2, my_rsq2a, 
                               my_rsq2b, my_ioa2)
    #
    if len(y_data) > 2:
        make_two_plots(y_data, y_fit1, y_fit2, my_txt1, my_txt2)
        pp.savefig()
        plt.close()

d = pp.infodict()
d['Title'] = 'GePiSaT GPP Plots with seasonal estimates of beta'
d['Author'] = 'Tyler W. Davis'
pp.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

##
## PLOT GPP BASED ON BASIC V. NEXT-GEN LUE FORMULA
##   depends: (1) get_veg_type, (2) model_fitness, (3) make_two_plots
##
fig_file = out_dir + 'GePiSaT_LUE_basic_v_nextgen.pdf'
pp = PdfPages(fig_file)

for station in all_stations:
    lue_data = numpy.copy(my_lue.st_lue_vars[station])
    x_data = lue_data['Iabs']
    x_data = numpy.vstack([x_data, numpy.zeros(len(x_data))]).T
    y_data = lue_data['GPP']
    #
    my_veg_type = get_veg_type(met_dir, station)
    #
    # Basic LUE:
    try:
        lue_est = numpy.linalg.lstsq(x_data, y_data)[0][0]
        #slope, intrcp, r, p, sterr = scipy.stats.linregress(x_data, y_data)
    except:
        lue_est = numpy.nan
    #
    y_fit1 = lue_est*lue_data['Iabs']
    my_rsq1a, my_rsq1b, my_ioa1 = model_fitness(y_fit1, y_data)
    my_txt1 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n"
              "$\\phi=%0.4f$\n"
              "$R^2=%0.3f$\n$r=%0.3f$\n"
              "$IA=%0.3f$") % (station, my_veg_type, lue_est, my_rsq1a, 
                               my_rsq1b, my_ioa1)
    
    #
    # Next-Gen LUE:
    b_est1, b_est2 = my_lue.calc_mgs_beta(station)
    y_fit2 = my_lue.next_gen_lue(lue_data, b_est2)
    my_rsq2a, my_rsq2b, my_ioa2 = model_fitness(y_fit2, y_data)
    my_txt2 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n$\\beta=%0.2f$\n"
               "$R^2=%0.3f$\n$r=%0.3f$\n"
               "$IA=%0.3f$") % (station, my_veg_type, b_est2, 
                                my_rsq2a, my_rsq2b, my_ioa2)
    #
    if len(y_data) > 2:
        make_two_plots(y_data, y_fit1, y_fit2, my_txt1, my_txt2, v=2)
        pp.savefig()
        plt.close()

d = pp.infodict()
d['Title'] = 'GePiSaT Basic v. Next-Gen LUE Plots'
d['Author'] = 'Tyler W. Davis'
pp.close()
