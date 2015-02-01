#!/usr/bin/python
#
# gepisat_nlsr.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-11-18 -- created
# 2015-02-01 -- last updated
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
# 
# ~~~~~
# todo:
# ~~~~~
# 1. add LUE class and other functions to GePiSaT's model.py
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import glob
import numpy
import os.path
from scipy.optimize import curve_fit
import scipy.stats
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
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, mod='NXGEN'):
        """
        Name:     LUE.__init__
        Input:    str, model type parameter (BASIC or NXGEN)
        Features: Initializes dictionaries for the light-use efficiency model
        """
        # Dictionary of stations' monthly meteorological values & their units:
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
                              'Patm' : 'Pa'}
        #
        # Define station value header line:
        header_vals = ['Timestamp', 'GPP', 'GPP_err', 'fPAR', 'PPFD',  
                       'VPD', 'CPA', 'Tair', 'CO2', 'Patm']
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
        # also set basic boolean
        self.st_lue_vars = {}
        if mod == 'BASIC':
            self.basic = True
            self.lue_var_units = {'Iabs' : 'mol_m2'}
        elif mod == 'NXGEN':
            self.basic = False
            self.lue_var_units = {'GPP' : 'mol_m2',
                                  'GPP_err' : 'mol_m2',
                                  'Iabs' : 'mol_m2', 
                                  'ca' : 'Pa',
                                  'Gs' : 'Pa', 
                                  'D' : 'Pa', 
                                  'K' : 'Pa', 
                                  'ns' : 'NA', 
                                  'fa' : 'NA'}
        else:
            print 'Light-use efficiency model type not recognized!'
            exit(1)
        #
        # Dictionary of stations' light-use efficiency model fit/fitness:
        self.station_lue = {}
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def add_station_val(self, station, month, gpp, gpp_err, 
                        fpar, ppfd, vpd, alpha, tair, co2, patm):
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
        Output:   None.
        Features: Appends a set of monthly values to the value dictionary
        Depends:  - calc_gstar
                  - calc_k
                  - viscosity_h2o
        """
        # Check for missing alpha (due to new STASH code):
        if alpha is None:
            alpha = -9999.0
        #
        # Place value parameters into a tuple:
        val_params = (month, gpp, gpp_err, fpar, ppfd, vpd, alpha, tair, co2, 
                      patm)
        #
        # Calculate lue variables & place into tuple:
        iabs = fpar*ppfd                    # mol/m2, abs. PPFD
        ca = (1.e-6)*co2*patm               # Pa, atms. CO2
        gs = self.calc_gstar(tair)          # Pa, photores. comp. point
        d = (1e3)*vpd                       # Pa, vapor pressure deficit
        k = self.calc_k(tair, patm)         # Pa, Michaelis-Menten coef.
        no = 1.00                           # mPa s, viscosity at 20 degC
        ns = self.viscosity_h2o(tair)/no    # unitless, water viscosity
        fa = (alpha/1.26)**(0.25)           # unitless, func. of alpha
        var_params = (gpp, gpp_err, iabs, ca, gs, d, k, ns, fa)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Meteorological Values
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # Initialize list if station key doesn't exist:
        if station not in self.station_vals.keys():
            self.station_vals[station] = numpy.array(
                val_params,
                dtype={'names' : ('Timestamp', 'GPP', 'GPP_err', 'fPAR', 'PPFD', 
                                  'VPD', 'CPA', 'Tair', 'CO2', 'Patm'),
                       'formats' : ('O', 'f4', 'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 'f4', 'f4')},
                        ndmin=1
                )
        else:
            # Add new parameters to list:
            temp_array = numpy.array(
                val_params,
                dtype={'names' : ('Timestamp', 'GPP', 'GPP_err', 'fPAR', 'PPFD', 
                                  'VPD', 'CPA', 'Tair', 'CO2', 'Patm'),
                       'formats' : ('O', 'f4', 'f4', 'f4', 'f4', 
                                    'f4', 'f4', 'f4', 'f4', 'f4')},
                        ndmin=1
                )
            self.station_vals[station] = numpy.append(
                self.station_vals[station], 
                temp_array, 
                axis=0
                )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~
        # LUE Model Variables
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
        # Define constants
        gs25 = 4.220  # Pa, assuming 25 deg C & 98.716 kPa)
        dha = 37830   # J/mol
        kR = 8.3145   # J/mol/K
        #
        gs = gs25*numpy.exp(dha*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
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
        kc = kc25*numpy.exp(dhac*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        ko = ko25*numpy.exp(dhao*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
        k = kc*(1 + kco*(1e-6)*patm/ko)
        return k
    #
    def viscosity_h2o(self, tc):
        """
        Name:     LUE.viscosity_h2o
        Input:    float, air temperature, degrees C (tc)
        Output:   float, viscosity, mPa s (n)
        Features: Returns the temperature-dependent viscosity of water (mPa s) 
                  based on the Vogel Equation
        """
        n = 2.4263e-2*numpy.exp(5.78919e2/((tc + 273.15) - 1.37546e2))
        return n
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
                    f.write(
                        ("%s,%0.4f,%0.4f,%0.4f,%0.4f,"
                        "%0.4f,%0.3f,%0.2f,%0.2f,%0.2f\n") % t)
                    f.close()
    #
    def write_out_lue(self, out_file):
        """
        Name:     LUE.write_out_lue
        Input:    string, output file (out_file)
        Output:   None.
        Features: Writes to file the calculated light use efficiency for all 
                  stations
                  
        @TODO: update for basic and next-gen models
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
        """
        # Make sure my_array is a numpy array:
        if not isinstance(my_array, numpy.ndarray):
            my_array = numpy.array(my_array)
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
        return (max_val, min_val, ave_val, std_val, skew_val, kurt_val)

###############################################################################
## FUNCTIONS:
###############################################################################
def next_gen_lue(x, beta):
    """
    Name:     next_gen_lue
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
    """
    # Define constants:
    kc = 0.41      # Jmax cost parameter
    phi_o = 0.093  # Long et al. (1993)
    #
    if beta >= 0:
        # Correct divide by zero errors:
        if beta == 0:
            beta += 1e-6
        #
        # Define variable substitutes:
        vdg = x['ca'] - x['Gs']
        vag = x['ca'] + 2.0*x['Gs']
        vsr = numpy.sqrt(
            1.6*x['ns']*x['D']/(beta*(x['K'] + x['Gs']))
        )
        #
        # Based on the m' formulation (see Regressing_LUE.pdf)
        temp_part = (vdg/(vag + 3.*x['Gs']*vsr))**2 - kc**(2./3.)*(
            vdg/(vag + 3.*x['Gs']*vsr))**(4./3.)
        if (temp_part <= 0).any():
            gpp = x['Iabs']
            gpp /= x['Iabs']
            gpp *= 1e20
        else:
            gpp = phi_o*x['Iabs']*x['fa']*numpy.sqrt(temp_part)
    else:
        gpp = x['Iabs']
        gpp /= x['Iabs']
        gpp *= 1e20
    return gpp

def predict_params(d, eta, gs, k):
    """
    Name:     predict_params
    Input:    -numpy.ndarray, atmos. CO2 conc. (ca)          DON'T NEED
              -numpy.ndarray, vap. press. deficit (d)
              -numpy.ndarray, water visc. (eta)
              -numpy.ndarray, GPP (gpp)                      DON'T NEED
              -numpy.ndarray, photo. resp. comp. point (gs)
              -numpy.ndarray, abs. PAR (iabs)                DON'T NEED
              -numpy.ndarray, Michaelis-Menten coef. (k)
    Output:   float
    Features: Returns statistically-based estimates of beta
    """
    # Predict phio (statistical relationship)
    #phio_p = ((1.761e-1)*gs['min'] + (1.702e-2)*gpp['std'] -
    #          (6.894e-3)*k['min'] - (3.531e-4)*iabs['ave'])
    #
    # Predict beta (based on Colin's method)
    my_d = 0.5*(d['max'] + d['min'])
    my_k = k['ave']
    my_g = 0.5*(gs['max'] + gs['min'])
    my_n = 0.5*(eta['max'] + eta['min'])
    #
    chi = 0.5 - (0.5 - 0.9)*(2.5e3 - my_d)/(2.5e3)
    beta_p = 1.6*my_d/(my_k + my_g)
    beta_p *= ((chi - my_g)/(1.0 - chi))**2
    beta_p *= my_n
    beta_p *= 0.0017154
    #
    return beta_p[0]

def calc_lue(lue_class, station):
    """
    Name:     calc_lue
    Input:    - LUE class (lue_class)
              - string, station name (station)
    Output:   None.
    Features: Fits the next generation LUE model to the monthly flux tower 
              data and saves the fit to LUE class
    Depends:  next_gen_lue
    """
    # Initialize fitness parameters:
    st_rsqr = -9999.     # model coef. of determination
    st_phio = 0.093      # intrinsic quantum efficiency parameter
    st_beta = -9999.     # beta parameter
    st_phio_err = 0.     # \ standard errors 
    st_beta_err = -9999. # /  of the estimates
    st_phio_t = 0.       # \ t-values 
    st_beta_t = -9999.   # /  of the estimates
    st_phio_p = 0.       # \ p-values
    st_beta_p = -9999.   # /  of the estimates
    #
    temp_stats = numpy.array(tuple([-9999., -9999., -9999., 
                                    -9999., -9999., -9999.]),
                             dtype={'names' : ('max', 'min', 'ave', 
                                               'std', 'skw', 'krt'),
                                    'formats' : ('f4', 'f4', 'f4', 
                                                 'f4', 'f4', 'f4')},
                             ndmin=1)
    #
    # Initialize parameter statistics:
    # max_val, min_val, ave_val, std_val, skew_val, kurt_val
    gpp_stats = numpy.copy(temp_stats)
    iabs_stats = numpy.copy(temp_stats)
    ca_stats = numpy.copy(temp_stats)
    gs_stats = numpy.copy(temp_stats)
    k_stats = numpy.copy(temp_stats)
    ns_stats = numpy.copy(temp_stats)
    vpd_stats = numpy.copy(temp_stats)
    fa_stats = numpy.copy(temp_stats)
    #
    if (station in my_lue.station_vals.keys() and 
        station in my_lue.st_lue_vars.keys()):
        #
        x_data = my_lue.st_lue_vars[station]
        y_data = my_lue.station_vals[station]['GPP']
        #
        # Calculate data statistics:
        gpp_stats[0] = my_lue.calc_statistics(y_data)
        iabs_stats[0] = my_lue.calc_statistics(x_data['Iabs'])
        ca_stats[0] = my_lue.calc_statistics(x_data['ca'])
        gs_stats[0] = my_lue.calc_statistics(x_data['Gs'])
        vpd_stats[0] = my_lue.calc_statistics(x_data['D'])
        k_stats[0] = my_lue.calc_statistics(x_data['K'])
        ns_stats[0] = my_lue.calc_statistics(x_data['ns'])
        fa_stats[0] = my_lue.calc_statistics(x_data['fa'])
        #
        # Estimate parameter starting values:
        est_beta = predict_params(vpd_stats, ns_stats, gs_stats, k_stats)
        #
        # Curve fit:
        try:
            fit_opt, fit_cov = curve_fit(next_gen_lue, 
                                         x_data, 
                                         y_data, 
                                         p0=est_beta)
            #print 'fit_opt', fit_opt
            #print 'fit_cov', fit_cov
        except:
            st_beta = -9999.
        else:
            st_beta = fit_opt[0]
            #
            fit_var = fit_cov[0][0]
            if numpy.isfinite(fit_var) and not (fit_var < 0):
                # Get parameter standard errors:
                st_beta_err = numpy.sqrt(fit_var)
                #
                # Calculate t-values:
                st_beta_t = st_beta/st_beta_err
                # 
                # Calculate p-values:
                st_beta_p = scipy.stats.t.pdf(-abs(st_beta_t), num_rows)
                #
                # Calculate r-squared:
                dum_x = next_gen_lue(x_data, st_beta)
                dum_slope, dum_intrcp, dum_r, dum_p, dum_sterr = (
                    scipy.stats.linregress(dum_x, y_data))
                st_rsqr = dum_r**2
    #
    # Save fit to LUE class
    params = (tuple([st_rsqr, st_phio, st_phio, est_beta, st_beta, 
                     st_phio_err, st_beta_err, st_phio_t, st_beta_t, 
                     st_phio_p, st_beta_p]) + 
                     tuple(gpp_stats[0]) +
                     tuple(iabs_stats[0]) +
                     tuple(ca_stats[0]) +
                     tuple(gs_stats[0]) +
                     tuple(vpd_stats[0]) +
                     tuple(k_stats[0]) +
                     tuple(eta_stats[0]))
    lue_class.station_lue[station] = params

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
    st_rsqr = -9999.     # model coef. of determination
    st_phio = -9999.     # intrinsic quantum efficiency parameter
    st_beta = -9999.     # beta parameter
    st_phio_err = -9999. # \ standard errors 
    st_beta_err = -9999. # /  of the estimates
    st_phio_t = -9999.   # \ t-values 
    st_beta_t = -9999.   # /  of the estimates
    st_phio_p = -9999.   # \ p-values
    st_beta_p = -9999.   # /  of the estimates
    #
    # Initialize data statistics:
    temp_stats = numpy.array(tuple([-9999., -9999., -9999., -9999., -9999., 
                                    -9999.]),
                             dtype={'names' : ('max', 'min', 'ave', 'std', 
                                               'skw', 'krt'),
                                    'formats' : ('f4', 'f4', 'f4', 'f4', 'f4', 
                                                 'f4')},
                             ndmin=1)
    #
    gpp_stats = numpy.copy(temp_stats)
    iabs_stats = numpy.copy(temp_stats)
    ca_stats = numpy.copy(temp_stats)
    gs_stats = numpy.copy(temp_stats)
    k_stats = numpy.copy(temp_stats)
    eta_stats = numpy.copy(temp_stats)
    vpd_stats = numpy.copy(temp_stats)
    #
    if station in lue_class.station_vals.keys():
        num_rows = len(lue_class.station_vals[station])
        for i in xrange(num_rows):
            (st_time, st_gpp, st_gpp_err, st_fpar, st_ppfd, st_vpd, st_cpa, 
            st_tair, st_co2, st_patm) = lue_class.station_vals[station][i]
            #
            # Calculate other necessary parameters for regression:
            st_ca = (1.e-6)*st_co2*st_patm       # Pa, atms. CO2
            st_gs = calc_gstar(st_tair)          # Pa, photores. comp. point
            st_k = calc_k(st_tair, st_patm)      # Pa, Michaelis-Menten coef.
            st_eta = viscosity_h2o(st_tair)      # mPa s, water viscosity
            st_iabs = st_fpar*st_ppfd            # mol/m2, abs. PPFD
            st_fa = (st_cpa/1.26)**(0.25)        # unitless, func. of alpha
            #
            # Filter variables out of range:
            if st_vpd < 0:
                st_vpd = numpy.nan
            #
            if i == 0:
                x_data = numpy.array(
                    tuple([st_iabs, st_ca, st_gs, st_vpd, st_k, st_eta, st_fa]),
                    dtype={'names' : ('Iabs', 'ca', 'Gs', 'D', 'K', 'eta', 'fa'),
                        'formats' : ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')},
                    ndmin=1
                )
                y_data = numpy.array([st_gpp,])
            else:
                x_temp = numpy.array(
                    tuple([st_iabs, st_ca, st_gs, st_vpd, st_k, st_eta, st_fa]),
                    dtype={'names' : ('Iabs', 'ca', 'Gs', 'D', 'K', 'eta', 'fa'),
                        'formats' : ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')},
                    ndmin=1
                )
                x_data = numpy.append(x_data, x_temp, axis=0)
                y_data = numpy.append(y_data, [st_gpp,])
        #
        # Remove nans from data sets:
        st_idx = numpy.where(~numpy.isnan(x_data['D']))[0]
        x_data = x_data[st_idx,]
        y_data = y_data[st_idx]
        num_rows = len(st_idx)
        #
        # Calculate predictor statistics:
        gpp_stats[0] = lue_class.calc_statistics(y_data)
        iabs_stats[0] = lue_class.calc_statistics(x_data['Iabs'])
        ca_stats[0] = lue_class.calc_statistics(x_data['ca'])
        gs_stats[0] = lue_class.calc_statistics(x_data['Gs'])
        vpd_stats[0] = lue_class.calc_statistics(x_data['D'])
        k_stats[0] = lue_class.calc_statistics(x_data['K'])
        eta_stats[0] = lue_class.calc_statistics(x_data['eta'])
        #
        est_phio, est_beta = predict_params(ca_stats, vpd_stats, eta_stats, 
                                            gpp_stats, gs_stats, iabs_stats, 
                                            k_stats)
        # Curve fit:
        try:
            fit_opt, fit_cov = curve_fit(next_gen_lue, 
                                        x_data, 
                                        y_data, 
                                        p0=[est_phio, est_beta])
        except:
            st_phio = -9999.
            st_beta = -9999.
        else:
            st_phio, st_beta = fit_opt
            #
            try:
                fit_var = numpy.diag(fit_cov)
            except ValueError:
                fit_var = [0.0, 0.0]
            else:
                if numpy.isfinite(fit_var).all() and not (fit_var < 0).any():
                    # Get parameter standard errors:
                    (st_phio_err, st_beta_err) = numpy.sqrt(fit_var)
                    #
                    # Calculate t-values:
                    st_phio_t = st_phio/st_phio_err
                    st_beta_t = st_beta/st_beta_err
                    # 
                    # Calculate p-values:
                    st_phio_p, st_beta_p = scipy.stats.t.pdf(
                        -abs(numpy.array([st_phio_t, st_beta_t])),
                        num_rows
                    )
                    #
                    # Calculate r-squared:
                    dum_x = next_gen_lue(x_data, st_phio, st_beta)
                    dum_slope, dum_intrcp, dum_r, dum_p, dum_sterr = (
                        scipy.stats.linregress(dum_x, y_data))
                    st_rsqr = dum_r**2
    #
    return (x_data, y_data, st_phio, st_beta, st_rsqr)

###############################################################################
## MAIN PROGRAM:
###############################################################################
mac = 0
if mac:
    my_dir = '/Users/twdavis/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    out_dir = '/Users/twdavis/Desktop/'
else:
    my_dir = '/home/user/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    out_dir = '/home/user/Desktop/'

my_files = glob.glob(my_dir + 'GePiSaT_All_Data-26*')
my_data = numpy.loadtxt(
    fname=my_files[0],
    dtype={
        'names': (
            'station', 'timestamp', 'gpp', 'gpp_err', 'fapar',
            'ppfd', 'vpd', 'alpha', 'tair', 'co2', 'patm'
        ),
        'formats' : (
            'S6', 'O', 'f4', 'f4', 'f4',
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4'
        )
    },
    delimiter=',',
    skiprows=1,
    converters={
        0 : numpy.str,
        1 : lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'),
        2 : numpy.float,                       # mol/m2
        3 : numpy.float,                       # mol/m2
        4 : numpy.float,                       # unitless
        5 : numpy.float,                       # mol/m2
        #6 : lambda x: numpy.float(x)*1e3,      # Pa
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
    (st_name, st_time, st_gpp, st_gpp_err, st_fpar, st_ppfd, st_vpd, st_cpa, 
     st_tair, st_co2, st_patm) = line
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
                           st_patm)

# List of each station:
all_stations = numpy.array(list(set(my_data['station'])))

# Test calc_lue function:
for station in all_stations:
    calc_lue(my_lue, station)

my_lue.write_out_lue(out_dir + "GePiSaT_nxgn.txt")




# For plot_mo_lue.py basic v advanced plot
#x_data, y_data, st_phio, st_beta, st_rsqr = get_lue(my_lue, 'ES-ES1')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEGIN my_ls_fit FUNCTION:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station = 'CZ-wet'
x_data = my_lue.st_lue_vars[station]
y_data = my_lue.st_lue_vars[station]['GPP']
w_data = numpy.power(my_lue.st_lue_vars[station]['GPP_err'], -1)
w_data /= w_data.sum()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEGIN calc_lue FUNCTION:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if 0:
    # Initialize multi-page PDF:
    fig_file = out_dir + 'GePiSaT_nxtgn_figs.pdf'
    pp = PdfPages(fig_file)
    
    temp_stats = numpy.array(tuple([-9999., -9999., -9999., 
                                    -9999., -9999., -9999.]),
                            dtype={'names' : ('max', 'min', 'ave', 
                                            'std', 'skw', 'krt'),
                                    'formats' : ('f4', 'f4', 'f4', 
                                                'f4', 'f4', 'f4')},
                            ndmin=1)
    #
    #station = 'ES-ES1'
    for station in numpy.sort(all_stations):
        # Initialize return values:
        st_rsqr = -9999.     # model coef. of determination
        st_phio = 0.093      # intrinsic quantum efficiency parameter
        st_beta = -9999.     # beta parameter
        st_phio_err = 0.     # \ standard errors 
        st_beta_err = -9999. # /  of the estimates
        st_phio_t = 0.       # \ t-values 
        st_beta_t = -9999.   # /  of the estimates
        st_phio_p = 0.       # \ p-values
        st_beta_p = -9999.   # /  of the estimates
        #
        # Initialize parameter statistics:
        # max_val, min_val, ave_val, std_val, skew_val, kurt_val
        gpp_stats = numpy.copy(temp_stats)
        iabs_stats = numpy.copy(temp_stats)
        ca_stats = numpy.copy(temp_stats)
        gs_stats = numpy.copy(temp_stats)
        k_stats = numpy.copy(temp_stats)
        eta_stats = numpy.copy(temp_stats)
        vpd_stats = numpy.copy(temp_stats)
        #cpa_stats = numpy.copy(temp_stats) #?
        #
        if station in my_lue.station_vals.keys():
            num_rows = len(my_lue.station_vals[station])
            for i in xrange(num_rows):
                (st_time, st_gpp, st_gpp_err, st_fpar, st_ppfd, st_vpd, st_cpa, 
                st_tair, st_co2, st_patm) = my_lue.station_vals[station][i]
                #
                # Calculate other necessary parameters for regression:
                st_ca = (1.e-6)*st_co2*st_patm       # Pa, atms. CO2
                st_gs = calc_gstar(st_tair)          # Pa, photores. comp. point
                st_k = calc_k(st_tair, st_patm)      # Pa, Michaelis-Menten coef.
                st_eta = viscosity_h2o(st_tair)      # mPa s, water viscosity
                st_iabs = st_fpar*st_ppfd            # mol/m2, abs. PPFD
                st_fa = (st_cpa/1.26)**(0.25)        # unitless, func. of alpha
                #
                # Filter variables out of range:
                if st_vpd < 0:
                    st_vpd = numpy.nan
                #
                if i == 0:
                    x_data = numpy.array(
                        tuple([st_iabs, st_ca, st_gs, st_vpd, st_k, st_eta, st_fa]),
                        dtype={'names' : ('Iabs', 'ca', 'Gs', 'D', 'K', 'eta', 'fa'),
                            'formats' : ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')},
                        ndmin=1
                    )
                    y_data = numpy.array([st_gpp,])
                else:
                    x_temp = numpy.array(
                        tuple([st_iabs, st_ca, st_gs, st_vpd, st_k, st_eta, st_fa]),
                        dtype={'names' : ('Iabs', 'ca', 'Gs', 'D', 'K', 'eta', 'fa'),
                            'formats' : ('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')},
                        ndmin=1
                    )
                    x_data = numpy.append(x_data, x_temp, axis=0)
                    y_data = numpy.append(y_data, [st_gpp,])
            #
            # Remove nans from data sets:
            st_idx = numpy.where(~numpy.isnan(x_data['D']))[0]
            x_data = x_data[st_idx,]
            y_data = y_data[st_idx]
            num_rows = len(st_idx)
            #
            # Calculate data statistics:
            gpp_stats[0] = my_lue.calc_statistics(y_data)
            iabs_stats[0] = my_lue.calc_statistics(x_data['Iabs'])
            ca_stats[0] = my_lue.calc_statistics(x_data['ca'])
            gs_stats[0] = my_lue.calc_statistics(x_data['Gs'])
            vpd_stats[0] = my_lue.calc_statistics(x_data['D'])
            k_stats[0] = my_lue.calc_statistics(x_data['K'])
            eta_stats[0] = my_lue.calc_statistics(x_data['eta'])
            #
            # Estimate parameter starting values:
            est_beta = predict_params(vpd_stats, eta_stats, gs_stats, k_stats)
            #
            # Curve fit:
            try:
                fit_opt, fit_cov = curve_fit(next_gen_lue, 
                                            x_data, 
                                            y_data, 
                                            p0=est_beta)
            except:
                st_beta = -9999.
            else:
                st_beta = fit_opt[0]
                #
                fit_var = fit_cov[0][0]
                if numpy.isfinite(fit_var) and not (fit_var < 0):
                    # Get parameter standard errors:
                    st_beta_err = numpy.sqrt(fit_var)
                    #
                    # Calculate t-values:
                    st_beta_t = st_beta/st_beta_err
                    # 
                    # Calculate p-values:
                    st_beta_p = scipy.stats.t.pdf(-abs(st_beta_t), num_rows)
                    #
                    # Calculate r-squared:
                    dum_x = next_gen_lue(x_data, st_beta)
                    dum_slope, dum_intrcp, dum_r, dum_p, dum_sterr = (
                        scipy.stats.linregress(dum_x, y_data))
                    st_rsqr = dum_r**2
                #
                # Plot fit:
                text_str = ("$\\mathrm{%s}$, $\\mathrm{w/}$ $f_\\alpha$\n"
                            "$\\phi_o=%0.4f\\pm%0.4f$\n"
                            "$\\beta=%0.2f\pm%0.2f$\n$R^2=%0.3f$") % (station, 
                            st_phio, st_phio_err, st_beta, st_beta_err, st_rsqr)
                fig = plt.figure();
                ax1 = fig.add_subplot(111)
                ax1.plot(next_gen_lue(x_data, st_phio, st_beta), y_data, 'ro', 
                        numpy.sort(y_data), numpy.sort(y_data), '--k', label='LUE')
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                ax1.text(0.05, 0.95, text_str, transform=ax1.transAxes, fontsize=12, 
                        verticalalignment='top', bbox=props)
                ax1.set_ylabel('Observed GPP, mol CO$_2$ m$^{-2}$')
                ax1.set_xlabel('Modeled GPP, mol CO$_2$ m$^{-2}$')
                pp.savefig()
                plt.close()
                #
    # Close PDF file
    d = pp.infodict()
    d['Title'] = 'GePiSaT Next-Gen LUE Plots with alpha'
    d['Author'] = 'Tyler W. Davis'
    pp.close()
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # END CALC_LUE
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # PLOT
    fig_file = out_dir + 'GePiSaT_nxtgn_figs.pdf'
    pp = PdfPages(fig_file)
    
    text_str = ("$\\mathrm{%s}$\n$\\phi_o=%0.4f\\pm%0.4f$\n"
                "$\\beta=%0.2f\pm%0.2f$\n$R^2=%0.3f$") % (
                station, st_phio, st_phio_err, st_beta, st_beta_err, st_rsqr
    )
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(next_gen_lue(x_data, st_phio, st_beta), y_data, 'ro', 
            numpy.sort(y_data), numpy.sort(y_data), '--k', label='LUE')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, text_str, transform=ax1.transAxes, fontsize=12, 
            verticalalignment='top', bbox=props)
    ax1.set_ylabel('Observed GPP, mol CO$_2$ m$^{-2}$')
    ax1.set_xlabel('Modeled GPP, mol CO$_2$ m$^{-2}$')
    #ax1.set_title(station)
    plt.show()
    
    pp.savefig()
    
    pp.close()
    
