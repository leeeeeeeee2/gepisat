#!/usr/bin/python
#
# flux_parti.py
#
# VERSION 2.2.0-dev
#
# LAST UPDATES: 2016-04-01
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

import numpy
from scipy.optimize import curve_fit
import scipy.stats

from utilities import calc_statistics
from utilities import goodness_of_fit
from utilities import pearsons_r
from utilities import peirce_dev


###############################################################################
# CLASSES
###############################################################################
class STAGE1:
    """
    Name:     STAGE1
    Features: This class performs flux partitioning of monthly NEE & PPFD
              observations
    History   Version 2.2.0-dev
              - class separated from model [16.01.17]
              - Python 2/3 supported print statements [16.01.17]
              - moved utility functions to utilities script [16.04.01]
                * calc statistics
                * goodness of fit
                * pearsons r
                * peirce dev
              - reduced function calls by adding modes [16.04.01]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    # Flux tower name & month being processed:
    name = ""
    month = datetime.date(1999, 1, 1)

    # PPFD and NEE arrays:
    # "h" and "l" indicate hyperbolic and linear models
    # "ro" indicates data has removed outliers
    ppfd_obs = numpy.array([])
    ppfd_obs_h_ro = numpy.array([])
    ppfd_obs_l_ro = numpy.array([])
    nee_obs = numpy.array([])
    nee_obs_h_ro = numpy.array([])
    nee_obs_l_ro = numpy.array([])
    nee_model_h = numpy.array([])
    nee_model_h_ro = numpy.array([])
    nee_model_l = numpy.array([])
    nee_model_l_ro = numpy.array([])

    # PPFD Statistics (working, observation backup, & model outlier backups):
    max_ppfd = 0.0
    min_ppfd = 0.0
    ave_ppfd = 0.0
    std_ppfd = 0.0
    skw_ppfd = 0.0
    krt_ppfd = 0.0

    max_ppfd_obs = 0.0
    min_ppfd_obs = 0.0
    ave_ppfd_obs = 0.0
    std_ppfd_obs = 0.0
    skw_ppfd_obs = 0.0
    krt_ppfd_obs = 0.0

    max_ppfd_ro_h = 0.0
    min_ppfd_ro_h = 0.0
    ave_ppfd_ro_h = 0.0
    std_ppfd_ro_h = 0.0
    skw_ppfd_ro_h = 0.0
    krt_ppfd_ro_h = 0.0

    max_ppfd_ro_l = 0.0
    min_ppfd_ro_l = 0.0
    ave_ppfd_ro_l = 0.0
    std_ppfd_ro_l = 0.0
    skw_ppfd_ro_l = 0.0
    krt_ppfd_ro_l = 0.0

    # NEE Statistics (working, observation backup, & model outlier backups):
    max_nee = 0.0
    min_nee = 0.0
    ave_nee = 0.0
    std_nee = 0.0
    skw_nee = 0.0
    krt_nee = 0.0

    max_nee_obs = 0.0
    min_nee_obs = 0.0
    ave_nee_obs = 0.0
    std_nee_obs = 0.0
    skw_nee_obs = 0.0
    krt_nee_obs = 0.0

    max_nee_ro_h = 0.0
    min_nee_ro_h = 0.0
    ave_nee_ro_h = 0.0
    std_nee_ro_h = 0.0
    skw_nee_ro_h = 0.0
    krt_nee_ro_h = 0.0

    max_nee_ro_l = 0.0
    min_nee_ro_l = 0.0
    ave_nee_ro_l = 0.0
    std_nee_ro_l = 0.0
    skw_nee_ro_l = 0.0
    krt_nee_ro_l = 0.0

    # Correlation coefficients (NEE v PPFD):
    r_obs = 0.0
    r_ro_h = 0.0
    r_ro_l = 0.0

    # Linear model coefficients and their standard errors [alpha, R]:
    # includes student's t statistic and p-value for model parameters
    # * alpha :: [unitless]
    # * R     :: [umol m-2 s-1]
    lm_estimates = [1.0, 1.0]
    lm_estimates_obs = [-9999.0, -9999.0]
    lm_estimates_ro = [-9999.0, -9999.0]
    lm_optimized = [1.0, 1.0]
    lm_optim_err = [1.0, 1.0]
    lm_optim_t = [0.0, 0.0]
    lm_optim_p = [0.0, 0.0]
    lm_optimized_ro = [1.0, 1.0]
    lm_optim_err_ro = [1.0, 1.0]
    lm_optim_t_ro = [0.0, 0.0]
    lm_optim_p_ro = [0.0, 0.0]

    # Hyperbolic model coefficients and their standard errors [Foo, alpha, R]:
    # includes student's t statistic and p-value for model parameters
    # * Foo   :: [umol m-2 s-1]
    # * alpha :: [unitless]
    # * R     :: [umol m-2 s-1]
    hm_estimates = [1.0, 1.0, 1.0]
    hm_estimates_obs = [-9999.0, -9999.0, -9999.0]
    hm_estimates_ro = [-9999.0, -9999.0, -9999.0]
    hm_optimized = [1.0, 1.0, 1.0]
    hm_optim_err = [1.0, 1.0, 1.0]
    hm_optim_t = [0.0, 0.0, 0.0]
    hm_optim_p = [0.0, 0.0, 0.0]
    hm_optimized_ro = [1.0, 1.0, 1.0]
    hm_optim_err_ro = [1.0, 1.0, 1.0]
    hm_optim_t_ro = [0.0, 0.0, 0.0]
    hm_optim_p_ro = [0.0, 0.0, 0.0]

    # Linear model fitting statistics:
    lm_rsqr = 0.0
    lm_rsqr_ro = 0.0
    lm_rmse = 0.0
    lm_rmse_ro = 0.0
    lm_mse = 0.0
    lm_mse_ro = 0.0
    lm_outliers = 0.0
    lm_resid_var = 0.0
    lm_resid_var_ro = 0.0

    # Hyperbolic model fitting statistics:
    hm_rsqr = 0.0
    hm_rsqr_ro = 0.0
    hm_rmse = 0.0
    hm_rmse_ro = 0.0
    hm_mse = 0.0
    hm_mse_ro = 0.0
    hm_outliers = 0
    hm_resid_var = 0.0
    hm_resid_var_ro = 0.0

    # Model selection:
    mod_select = 0

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, ppfd, nee, tower, mo):
        """
        Name:     STAGE1.__init__
        Input:    - list, monthly PPFD observations (ppfd)
                  - list, monthly NEE observations (nee)
                  - string, flux tower name (tower)
                  - datetime.date, current month (mo)
        Features: Initialize the class with observation data, calculate basic
                  statistics, estimate model parameters, and calculate
                  Pearson's correlation coefficient
        Depends:  - calc_statistics
                  - pearsons_r
                  - save_estimates
                  - save_stats
                  - update_guess
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.info("STAGE1 class called")

        # Save tower name & month being processed:
        self.name = tower
        self.month = mo

        # Save PPFD and NEE lists:
        self.ppfd_obs = ppfd
        self.nee_obs = nee

        # Calculate statistics for observations:
        (
            self.max_ppfd,
            self.min_ppfd,
            self.ave_ppfd,
            self.std_ppfd,
            self.skw_ppfd,
            self.krt_ppfd) = calc_statistics(ppfd)
        (
            self.max_nee,
            self.min_nee,
            self.ave_nee,
            self.std_nee,
            self.skw_nee,
            self.krt_nee) = calc_statistics(nee)
        self.save_stats(obs=1, h=-1)

        # Update model parameters:
        self.update_guess()
        self.save_estimates(obs=1, h=-1)

        # Calculate Pearson's correlation coefficient:
        self.r_obs = pearsons_r(nee, ppfd)

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def calc_gpp(self, ppfd):
        """
        Name:     STAGE1.calc_gpp
        Input:    numpy.ndarray, monthly PPFD (ppfd)
        Output:   tuple, modeled GPP and associated error
                  - numpy.ndarray, modeled GPP (gpp)
                  - numpy.ndarray, associated error (gpp_err)
        Features: Returns GPP (umol m-2 s-1) and associated error based on the
                  best model
        Ref:      Chapters 11.2 & 11.3, GePiSaT Documentation
        """
        self.logger.debug(
            "calculating GPP for model selection %d", self.mod_select)

        # Get model selection:
        if (self.mod_select == 1 or self.mod_select == 2):
            model = "H"
        elif (self.mod_select == 3 or self.mod_select == 4):
            model = "L"

        # Get outlier selection:
        if (self.mod_select == 1 or self.mod_select == 3):
            outlier = 0
        elif (self.mod_select == 2 or self.mod_select == 4):
            outlier = 1

        if self.mod_select == 0:
            # Default to hyperbolic model with outliers removed estimates
            (foo, alpha, r) = self.hm_estimates_ro
            (foo_err, alpha_err, r_err) = (0., 0., 0.)

            if foo == 0 or alpha == 0:
                # Use linear estimates instead:
                (alpha, r) = self.lm_estimates_ro
                (alpha_err, r_err) = (0., 0.)

                # Calculate GPP and GPP err:
                gpp = (alpha*ppfd)
                gpp_err = (ppfd*alpha_err)
            else:
                # Variable substitutes:
                apf = alpha*ppfd + foo
                afp = alpha*foo*ppfd
                afpp = alpha*foo*(ppfd**2)

                # Calculate GPP and its associated error:
                gpp = afp/apf
                gpp_err = numpy.sqrt(
                    (alpha_err**2)*((foo*ppfd*apf - afpp)/(apf**2))**2 +
                    (foo_err**2)*((alpha*ppfd*apf - afp)/(apf**2))**2
                )

            # Clip out negative values of GPP:
            gpp = gpp.clip(min=0)
        elif model.upper() == "H":
            # Retrieve parameters for hyperbolic model:
            if outlier == 0:
                (foo, alpha, r) = self.hm_optimized
                (foo_err, alpha_err, r_err) = self.hm_optim_err
            elif outlier == 1:
                (foo, alpha, r) = self.hm_optimized_ro
                (foo_err, alpha_err, r_err) = self.hm_optim_err_ro

            # Variable substitutes:
            apf = alpha*ppfd + foo
            afp = alpha*foo*ppfd
            afpp = alpha*foo*(ppfd**2)

            # Calculate GPP and its associated error:
            gpp = afp/apf
            gpp_err = numpy.sqrt(
                (alpha_err**2)*((foo*ppfd*apf - afpp)/(apf**2))**2 +
                (foo_err**2)*((alpha*ppfd*apf - afp)/(apf**2))**2
            )
        elif model.upper() == "L":
            # Retrieve parameters for linear model:
            if outlier == 0:
                (alpha, r) = self.lm_optimized
                (alpha_err, r_err) = self.lm_optim_err
            elif outlier == 1:
                (alpha, r) = self.lm_optimized_ro
                (alpha_err, r_err) = self.lm_optim_err_ro

            # Calculate GPP:
            gpp = (alpha*ppfd)
            gpp_err = (ppfd*alpha_err)

        # Return GPP
        return (gpp, gpp_err)

    def calc_model(self, mode, ro):
        """
        Name:     STAGE1.calc_model
        Input:    - str, model mode (mode)
                    'H' for hyperbolic
                    'L' for linear
                  - bool, remove outliers (ro)
        Output:   None.
        Features: Calculates NEE and the fitness statistics using the
                  optimization parameters for the hyperbolic or linear model
                  based on observation data with or without outliers removed
        Depends:  - model_h
                  - model_l
                  - goodness_of_fit
        """
        self.logger.debug("Model mode %s", mode)
        self.logger.debug("Outlier removal: %s", ro)

        if mode == "H":
            if ro:
                self.logger.debug(
                    "calculating hyperbolic model without outliers")
                self.nee_model_h_ro = self.model_h(self.ppfd_obs_h_ro,
                                                   self.hm_optimized_ro[0],
                                                   self.hm_optimized_ro[1],
                                                   self.hm_optimized_ro[2])

                if -9999 in self.hm_optimized_ro:
                    self.hm_mse_ro = -9999.0
                    self.hm_rmse_ro = -9999.0
                    self.hm_rsqr_ro = -9999.0
                else:
                    (self.hm_mse_ro,
                     self.hm_rmse_ro,
                     self.hm_rsqr_ro) = goodness_of_fit(self.nee_model_h_ro,
                                                        self.nee_obs_h_ro, 3)
            else:
                self.logger.debug("calculating hyperbolic model")
                self.nee_model_h = self.model_h(self.ppfd_obs,
                                                self.hm_optimized[0],
                                                self.hm_optimized[1],
                                                self.hm_optimized[2])

                if -9999 in self.hm_optimized:
                    self.hm_mse = -9999.0
                    self.hm_rmse = -9999.0
                    self.hm_rsqr = -9999.0
                else:
                    (self.hm_mse,
                     self.hm_rmse,
                     self.hm_rsqr) = goodness_of_fit(self.nee_model_h,
                                                     self.nee_obs, 3)
        elif mode == "L":
            if ro:
                self.logger.debug("calculating linear model without outliers")
                self.nee_model_l_ro = self.model_l(self.ppfd_obs_l_ro,
                                                   self.lm_optimized_ro[0],
                                                   self.lm_optimized_ro[1])

                if -9999 in self.lm_optimized_ro:
                    self.lm_mse_ro = -9999.0
                    self.lm_rmse_ro = -9999.0
                    self.lm_rsqr_ro = -9999.0
                else:
                    (self.lm_mse_ro,
                     self.lm_rmse_ro,
                     self.lm_rsqr_ro) = goodness_of_fit(self.nee_model_l_ro,
                                                        self.nee_obs_l_ro, 2)
            else:
                self.logger.debug("calculating linear model")
                self.nee_model_l = self.model_l(self.ppfd_obs,
                                                self.lm_optimized[0],
                                                self.lm_optimized[1])

                if -9999 in self.lm_optimized:
                    self.lm_mse = -9999.0
                    self.lm_rmse = -9999.0
                    self.lm_rsqr = -9999.0
                else:
                    (self.lm_mse,
                     self.lm_rmse,
                     self.lm_rsqr) = goodness_of_fit(self.nee_model_l,
                                                     self.nee_obs, 2)

    def model_h(self, x, Foo, alpha, R):
        """
        Name:     STAGE1.model_h
        Input:    - numpy.ndarray, monthly PPFD (x)
                  - float, hyperbolic parameter (Foo)
                  - float, hyperbolic parameter (alpha)
                  - float, hyperbolic parameter (R)
        Output:   numpy.ndarray, modeled NEE
        Features: Returns array of NEE based on the hyperbolic model in the
                  form: Y = (ax+b)/(bx+c)
        Ref:      Eq. 15, GePiSaT Documentation
        """
        a = 1.0*(alpha*R) - 1.0*(alpha*Foo)
        b = 1.0*(Foo*R)
        c = 1.0*alpha
        d = 1.0*Foo

        numerator = 1.0*(a*x) + b
        denominator = 1.0*(c*x) + d

        # Compensate for zero division by adding a small number to each value:
        if 0 in denominator:
            denominator = [elem + 1.0e-6 for elem in denominator]
            denominator = numpy.array(denominator)

        return (numerator/denominator)

    def model_l(self, x, alpha, R):
        """
        Name:     STAGE1.model_l
        Input:    - numpy.ndarray, monthly PPFD (x)
                  - float, linear parameter (alpha)
                  - float, linear parameter (R)
        Output:   - numpy.ndarray, modeled NEE
        Features: Returns array of NEE based on the linear model in the form:
                  Y = ax+b
        Ref:      Eq. 14, GePiSaT Documentation
        """
        a = -1.0*alpha
        b = R

        return (a*x + b)

    def model_selection(self):
        """
        Name:     STAGE1.model_selection
        Input:    None.
        Output:   None.
        Features: Saves the model (i.e., none, hyperbolic with
                  observations, hyperbolic with outliers removed, linear with
                  observations, linear with outliers removed) that ``best''
                  represents the data based on thesholds of fitness, validity
                  of parameter values, and parameter significance
        """
        # Define threshold values:
        r2_thresh = 0.5       # minimum R-squared fitness
        p_thresh = 0.05       # minimum parameter significance (95%)

        # Initialize selection dictionary:
        select_dict = {'obs_h': 1,
                       'ro_h': 2,
                       'obs_l': 3,
                       'ro_l': 4}

        # Initialize model failure dictionary:
        model_dict = {'obs_h': 0,
                      'ro_h': 0,
                      'obs_l': 0,
                      'ro_l': 0}

        # Initialize R-squared dictionary:
        r2_dict = {'obs_h': self.hm_rsqr,
                   'obs_l': self.lm_rsqr,
                   'ro_h': self.hm_rsqr_ro,
                   'ro_l': self.lm_rsqr_ro}

        # Initialize P-value dictionary:
        p_dict = {'obs_h': numpy.array(self.hm_optim_p),
                  'ro_h': numpy.array(self.hm_optim_p_ro),
                  'obs_l': numpy.array(self.lm_optim_p),
                  'ro_l': numpy.array(self.lm_optim_p_ro)}

        # Initialize parameter value dictionary:
        value_dict = {
            'obs_h': {'foo': self.hm_optimized[0],
                      'alpha': self.hm_optimized[1],
                      'r': self.hm_optimized[2]
                      },
            'ro_h': {'foo': self.hm_optimized_ro[0],
                     'alpha': self.hm_optimized_ro[1],
                     'r': self.hm_optimized_ro[2]
                     },
            'obs_l': {'alpha': self.lm_optimized[0],
                      'r': self.lm_optimized[1]
                      },
            'ro_l': {'alpha': self.lm_optimized_ro[0],
                     'r': self.lm_optimized_ro[1]
                     }
        }

        # Initialize model selection:
        best_model = 0

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIRST CRITERIA--- R-SQUARED THRESHOLD
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (numpy.array(r2_dict.values()) < r2_thresh).any():
            # Eliminate those with poor R2:
            m_failed = [k for k, v in r2_dict.items() if v < r2_thresh]
            for m in m_failed:
                model_dict[m] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # SECOND CRITERIA--- PARAMETER RANGE OF VALIDITY
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # -- rectangular hypberbola:
        if (value_dict['obs_h']['foo'] > 100 or
                value_dict['obs_h']['foo'] < 0 or
                value_dict['obs_h']['alpha'] > 1 or
                value_dict['obs_h']['alpha'] < 0):
            model_dict['obs_h'] = 1
        if (value_dict['ro_h']['foo'] > 100 or
                value_dict['ro_h']['foo'] < 0 or
                value_dict['ro_h']['alpha'] > 1 or
                value_dict['ro_h']['alpha'] < 0):
            model_dict['ro_h'] = 1

        # -- linear:
        if (value_dict['obs_l']['alpha'] > 1 or
                value_dict['obs_l']['alpha'] < 0):
            model_dict['obs_l'] = 1
        if (value_dict['ro_l']['alpha'] > 1 or
                value_dict['ro_l']['alpha'] < 0):
            model_dict['ro_l'] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # THIRD CRITERIA--- PARAMETER SIGNIFICANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (p_dict['obs_h'] > p_thresh).any():
            model_dict['obs_h'] = 1
        if (p_dict['ro_h'] > p_thresh).any():
            model_dict['ro_h'] = 1
        if (p_dict['obs_l'] > p_thresh).any():
            model_dict['obs_l'] = 1
        if (p_dict['ro_l'] > p_thresh).any():
            model_dict['ro_l'] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ANALYSIS:
        # Let's see which models passed the selection criteria and how many
        # of them there are
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        mods = [k for k, v in model_dict.items() if v == 0]
        howmany = len(mods)

        if howmany > 0:
            # At least one model passed the initial criteria
            if howmany == 1:
                # There is one best choice:
                best_model = select_dict[mods[0]]
            elif howmany == 2:
                # Choose the better of the two:
                if 1.01*r2_dict[mods[1]] >= 0.99*r2_dict[mods[0]]:
                    best_model = select_dict[mods[1]]
                else:
                    best_model = select_dict[mods[0]]
            else:
                # Gather R-squared and worst p-values for remaining models:
                mods_rs = [r2_dict[m] for m in mods]
                mods_ps = [p_dict[m].max() for m in mods]

                # Get the max R-squared and min p-value:
                max_rs = max(mods_rs)
                min_ps = min(mods_ps)

                # How many models meet this criteria:
                mods_temp = [
                    m for m in mods if
                    mods_rs[mods.index(m)] == max_rs and
                    mods_ps[mods.index(m)] == min_ps
                    ]
                n_temp = len(mods_temp)

                if n_temp == 0:
                    # There is no "best" model, go through default heirarchy:
                    if 'ro_h' in mods:
                        best_model = select_dict['ro_h']
                    elif 'ro_l' in mods:
                        best_model = select_dict['ro_l']
                    elif 'obs_l' in mods:
                        best_model = select_dict['obs_l']
                    else:
                        best_model = select_dict['obs_h']
                elif n_temp == 1:
                    # Select the best model:
                    best_model = select_dict[mods_temp[0]]
                elif n_temp == 2:
                    # Choose the better of the two:
                    if (1.01*r2_dict[mods_temp[1]] >=
                            0.99*r2_dict[mods_temp[0]]):
                        best_model = select_dict[mods_temp[1]]
                    else:
                        best_model = select_dict[mods_temp[0]]
                else:
                    # There are 3 or 4 models
                    # that all having the same max R2 and min p-value.
                    # Default choices: ro_h, ro_l, obs_l, obs_h
                    if 'ro_h' in mods_temp:
                        best_model = select_dict['ro_h']
                    elif 'ro_l' in mods_temp:
                        best_model = select_dict['ro_l']
                    elif 'obs_l' in mods_temp:
                        best_model = select_dict['obs_l']
                    else:
                        best_model = select_dict['obs_h']

        # Save model selection:
        self.mod_select = best_model

    def partition(self, to_write, rm_out):
        """
        Name:     STAGE1.partition
        Input:    - numpy.ndarray, month of NEE observations (nee)
                  - numpy.ndarray, month of PPFD observations (ppfd)
                  - int, write to file boolean (to_write)
                  - int, outlier removal boolean (rm_out)
                  - str, station name (tower)
                  - datetime.date (month)
        Output:   FLUX_PARTI class object (my_class)
        Features: Returns a class with flux partitioning results based on
                  NEE:PPFD observations; processes outliers (optional);
                  saves results to file (optional)
        """
        ## ~~~~~~~~~~ ANALYZE OBSERVATIONS ~~~~~~~~~~ ##
        # ---------------------
        # Model H optimization:
        # ---------------------
        try:
            mh_opt, mh_cov = curve_fit(self.model_h,
                                       self.ppfd_obs,
                                       self.nee_obs,
                                       p0=self.hm_estimates)
        except:
            self.hm_optimized = [-9999.0, -9999.0, -9999.0]
            self.hm_optim_err = [-9999.0, -9999.0, -9999.0]
        else:
            (self.hm_optimized[0],
             self.hm_optimized[1],
             self.hm_optimized[2]) = mh_opt

            # Extract variance values from matrix:
            try:
                mh_var = numpy.diag(mh_cov)
            except ValueError:
                mh_var = [0, 0, 0]
            else:
                if numpy.isfinite(mh_var).all() and not (mh_var < 0).any():
                    (self.hm_optim_err[0],
                     self.hm_optim_err[1],
                     self.hm_optim_err[2]) = numpy.sqrt(mh_var)

                    # Calculate the t-value:
                    (self.hm_optim_t[0],
                     self.hm_optim_t[1],
                     self.hm_optim_t[2]) = (numpy.array(self.hm_optimized) /
                                            numpy.array(self.hm_optim_err))

                    # Calculate the p-value:
                    hm_df = len(self.nee_obs) - 3.0  # degrees of freedom
                    (self.hm_optim_p[0],
                     self.hm_optim_p[1],
                     self.hm_optim_p[2]) = scipy.stats.t.pdf(
                        -abs(numpy.array(self.hm_optim_t)), hm_df)
                else:
                    self.hm_optim_err[0] = 0.0
                    self.hm_optim_err[1] = 0.0
                    self.hm_optim_err[2] = 0.0
        finally:
            self.calc_model(mode="H", ro=False)

        # ---------------------
        # Model L optimization:
        # ---------------------
        try:
            ml_opt, ml_cov = curve_fit(self.model_l,
                                       self.ppfd_obs,
                                       self.nee_obs,
                                       p0=self.lm_estimates)
        except:
            self.lm_optimized = [-9999.0, -9999.0]
            self.lm_optim_err = [-9999.0, -9999.0]
        else:
            (self.lm_optimized[0],
             self.lm_optimized[1]) = ml_opt

            # Extract variance values from matrix:
            try:
                ml_var = numpy.diag(ml_cov)
            except ValueError:
                ml_var = [0, 0]
            else:
                if numpy.isfinite(ml_var).all() and not (ml_var < 0).any():
                    (self.lm_optim_err[0],
                     self.lm_optim_err[1]) = numpy.sqrt(ml_var)

                    # Calculate the t-value:
                    (self.lm_optim_t[0],
                     self.lm_optim_t[1]) = (numpy.array(self.lm_optimized) /
                                            numpy.array(self.lm_optim_err))

                    # Calculate the p-value:
                    lm_df = len(self.nee_obs) - 2  # degrees of freedom
                    (self.lm_optim_p[0],
                     self.lm_optim_p[1]) = scipy.stats.t.pdf(
                        -abs(numpy.array(self.lm_optim_t)), lm_df)
                else:
                    self.lm_optim_err[0] = 0.0
                    self.lm_optim_err[1] = 0.0
        finally:
            self.calc_model(mode="L", ro=False)

        ## ~~~~~~~~~~ WRITE OBSERVATION RESULTS ~~~~~~~~~~ ##
        if to_write:
            output_file = "out/%s_%s.txt" % (self.name, self.month)

            header0 = "MH_guess,%f,%f,%f\n" % (self.hm_estimates[1],
                                               self.hm_estimates[2],
                                               self.hm_estimates[0])

            header1 = "ML_guess,%f,%f\n" % (self.lm_estimates[0],
                                            self.lm_estimates[1])

            header2 = "MH_opt,%f,%f,%f,%f,%f\n" % (self.hm_optimized[1],
                                                   self.hm_optimized[2],
                                                   self.hm_optimized[0],
                                                   self.hm_rmse,
                                                   self.hm_rsqr)

            header3 = "ML_opt,%f,%f,,%f,%f\n" % (self.lm_optimized[0],
                                                 self.lm_optimized[1],
                                                 self.lm_rmse,
                                                 self.lm_rsqr)

            OUTFILE = open(output_file, 'w')
            OUTFILE.write(header0)
            OUTFILE.write(header1)
            OUTFILE.write(header2)
            OUTFILE.write(header3)
            OUTFILE.write("ppfd_obs,nee_obs,nee_mh,nee_ml\n")
            for po, no, nh, nl in map(None,
                                      self.ppfd_obs,
                                      self.nee_obs,
                                      self.nee_model_h,
                                      self.nee_model_l):
                outline = "%s,%s,%s,%s\n" % (po, no, nh, nl)
                OUTFILE.write(outline)
            OUTFILE.close()
        ## ~~~~~~~~~~ REMOVE OUTLIERS AND RE-ANALYZE ~~~~~~~~~~ ##
        if rm_out:
            # ---------------------
            # Model H optimization:
            # ---------------------
            self.remove_outliers(mode="H")
            try:
                mh_opt_ro, mh_cov_ro = curve_fit(self.model_h,
                                                 self.ppfd_obs_h_ro,
                                                 self.nee_obs_h_ro,
                                                 p0=self.hm_estimates)
            except:
                self.hm_optimized_ro = [-9999.0, -9999.0, -9999.0]
                self.hm_optim_err_ro = [-9999.0, -9999.0, -9999.0]
            else:
                (self.hm_optimized_ro[0],
                 self.hm_optimized_ro[1],
                 self.hm_optimized_ro[2]) = mh_opt_ro

                # Extract variance values from matrix:
                try:
                    mh_var_ro = numpy.diag(mh_cov_ro)
                except ValueError:
                    mh_var_ro = [0, 0, 0]
                else:
                    if (numpy.isfinite(mh_var_ro).all() and not
                            (mh_var_ro < 0).any()):
                        (self.hm_optim_err_ro[0],
                         self.hm_optim_err_ro[1],
                         self.hm_optim_err_ro[2]) = numpy.sqrt(mh_var_ro)

                        # Calculate the t-value:
                        (self.hm_optim_t_ro[0],
                         self.hm_optim_t_ro[1],
                         self.hm_optim_t_ro[2]) = (
                            numpy.array(self.hm_optimized_ro) /
                            numpy.array(self.hm_optim_err_ro))

                        # Calculate the p-value:
                        hm_df_ro = len(self.nee_obs) - self.hm_outliers - 3
                        (self.hm_optim_p_ro[0],
                         self.hm_optim_p_ro[1],
                         self.hm_optim_p_ro[2]) = scipy.stats.t.pdf(
                            -abs(numpy.array(self.hm_optim_t_ro)), hm_df_ro)
                    else:
                        self.hm_optim_err_ro[0] = 0.0
                        self.hm_optim_err_ro[1] = 0.0
                        self.hm_optim_err_ro[2] = 0.0
            finally:
                self.calc_model(mode="H", ro=True)

            # ---------------------
            # Model L optimization:
            # ---------------------
            self.remove_outliers(mode="L")
            try:
                ml_opt_ro, ml_cov_ro = curve_fit(self.model_l,
                                                 self.ppfd_obs_l_ro,
                                                 self.nee_obs_l_ro,
                                                 p0=self.lm_estimates)
            except:
                self.lm_optimized_ro = [-9999.0, -9999.0]
                self.lm_optim_err_ro = [-9999.0, -9999.0]
            else:
                (self.lm_optimized_ro[0],
                 self.lm_optimized_ro[1]) = ml_opt_ro

                # Extract variance values from matrix:
                try:
                    ml_var_ro = numpy.diag(ml_cov_ro)
                except ValueError:
                    ml_var_ro = [0, 0]
                else:
                    if (numpy.isfinite(ml_var_ro).all() and not
                            (ml_var_ro < 0).any()):
                        (self.lm_optim_err_ro[0],
                         self.lm_optim_err_ro[1]) = numpy.sqrt(ml_var_ro)

                        # Calculate the t-value:
                        (self.lm_optim_t_ro[0],
                         self.lm_optim_t_ro[1]) = (
                            numpy.array(self.lm_optimized_ro) /
                            numpy.array(self.lm_optim_err_ro))

                        # Calculate the p-value:
                        lm_df_ro = len(self.nee_obs) - self.lm_outliers - 2
                        (self.lm_optim_p_ro[0],
                         self.lm_optim_p_ro[1]) = scipy.stats.t.pdf(
                            -abs(numpy.array(self.lm_optim_t_ro)), lm_df_ro)
                    else:
                        self.lm_optim_err_ro[0] = 0.0
                        self.lm_optim_err_ro[1] = 0.0
            finally:
                self.calc_model(mode="L", ro=True)

            ## ~~~~~~~~~~ WRITE OUTLIER-FREE RESULTS ~~~~~~~~~~ ##
            if to_write:
                output_file = "out/%s_%s_ro.txt" % (self.name, self.month)
                header0 = "MH_guess,%0.5f,%0.2f,%0.2f\n" % (
                    self.hm_estimates[1],
                    self.hm_estimates[2],
                    self.hm_estimates[0])
                header1 = "ML_guess,%0.5f,%0.2f\n" % (
                    self.lm_estimates[0],
                    self.lm_estimates[1])
                header2 = "MH_opt,%0.5f,%0.2f,%0.2f,%0.2f,%0.3f\n" % (
                    self.hm_optimized_ro[1],
                    self.hm_optimized_ro[2],
                    self.hm_optimized_ro[0],
                    self.hm_rmse_ro,
                    self.hm_rsqr_ro)
                header3 = "ML_opt,%0.5f,%0.2f,,%0.2f,%0.3f\n" % (
                    self.lm_optimized_ro[0],
                    self.lm_optimized_ro[1],
                    self.lm_rmse_ro,
                    self.lm_rsqr_ro)
                OUTFILE = open(output_file, 'w')
                OUTFILE.write(header0)
                OUTFILE.write(header1)
                OUTFILE.write(header2)
                OUTFILE.write(header3)
                OUTFILE.write("ppfd_obs_h,nee_obs_h,nee_mod_h,"
                              "ppfd_obs_l,nee_obs_l,nee_mod_l\n")
                for (poh, noh, nhh, pol, nol, nll) in map(None,
                                                          self.ppfd_obs_h_ro,
                                                          self.nee_obs_h_ro,
                                                          self.nee_model_h_ro,
                                                          self.ppfd_obs_l_ro,
                                                          self.nee_obs_l_ro,
                                                          self.nee_model_l_ro):
                    outline = "%s,%s,%s,%s,%s,%s\n" % (
                        poh, noh, nhh, pol, nol, nll)
                    OUTFILE.write(outline)
                OUTFILE.close()

        # Run model selection criteria:
        self.model_selection()

    def remove_outliers(self, mode):
        """
        Name:     STAGE1.remove_outliers
        Input:    str, model mode (mode)
                    "H" hyperbolic
                    "L" linear
        Output:   int, success flag (rval)
        Features: Returns flag indicating the successful removal of outliers
                  from model observations; saves outlier-free data, new
                  statistics, and updated model estimates
        Depends:  - peirce_dev
                  - calc_statistics
                  - save_stats
                  - update_guess
                  - save_estimates
                  - pearsons_r
        """
        if mode == "H":
            self.logger.debug("removing outliers from hyperbolic model")
            peirce_m = 3
            mse = self.hm_mse
            nee_model = self.nee_model_h
        elif mode == "L":
            self.logger.debug("removing outliers from linear model")
            peirce_m = 2
            mse = self.lm_mse
            nee_model = self.nee_model_l

        # Set Peirce values:
        peirce_cap_n = len(self.nee_obs)
        peirce_lc_n = 1

        # Calculate tolerance
        peirce_x2 = peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
        peirce_delta2 = mse*peirce_x2

        # Calculate the square errors:
        sq_errors = (self.nee_obs - nee_model)**2.0

        # Find if/where exceedance occurs:
        outliers_index = numpy.where(sq_errors > peirce_delta2)[0]
        outliers_found = len(outliers_index)

        # Run check again if no outliers are found in first attempt:
        if (outliers_found == 0):
            peirce_lc_n = 2
            peirce_x2 = peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
            peirce_delta2 = mse*peirce_x2
            outliers_index = numpy.where(sq_errors > peirce_delta2)[0]
            outliers_found = len(outliers_index)

            # Reset n
            peirce_lc_n = 1

        # Increment n until it is greater than number of outliers found:
        while (peirce_lc_n <= outliers_found):
            peirce_lc_n += 1

            # Check that n < N:
            if peirce_lc_n >= peirce_cap_n:
                peirce_lc_n = outliers_found + 1.0
            else:
                peirce_x2 = peirce_dev(
                    peirce_cap_n, peirce_lc_n, peirce_m)
                peirce_delta2 = mse*peirce_x2
                outliers_index = numpy.where(sq_errors > peirce_delta2)[0]
                outliers_found = len(outliers_index)

        ppfd_obs_ro = numpy.delete(self.ppfd_obs, outliers_index)
        nee_obs_ro = numpy.delete(self.nee_obs, outliers_index)
        if ppfd_obs_ro.any() and nee_obs_ro.any():
            # Recalculate statistics and update guess values:
            (self.max_ppfd,
             self.min_ppfd,
             self.ave_ppfd,
             self.std_ppfd,
             self.skw_ppfd,
             self.krt_ppfd) = calc_statistics(ppfd_obs_ro)

            (self.max_nee,
             self.min_nee,
             self.ave_nee,
             self.std_nee,
             self.skw_nee,
             self.krt_nee) = calc_statistics(nee_obs_ro)

            pear = pearsons_r(ppfd_obs_ro, nee_obs_ro)
            rval = 1
        else:
            pear = -9999
            rval = 0

        if mode == "H":
            self.hm_outliers = outliers_found

            # Remove outliers found:
            self.ppfd_obs_h_ro = ppfd_obs_ro
            self.nee_obs_h_ro = nee_obs_ro
            self.r_ro_h = pear

            self.save_stats(obs=0, h=1)
            self.update_guess()
            self.save_estimates(obs=0, h=1)
        elif mode == "L":
            self.lm_outliers = outliers_found

            # Remove outliers found:
            self.ppfd_obs_l_ro = ppfd_obs_ro
            self.nee_obs_l_ro = nee_obs_ro
            self.r_ro_l = pear

            self.save_stats(obs=0, h=0)
            self.update_guess()
            self.save_estimates(obs=0, h=0)

        return rval

    def save_estimates(self, obs, h):
        """
        Name:     STAGE1.save_estimates
        Input:    - int, observation flag (obs)
                  - int, hyperbolic model flag (h)
        Output:   None.
        Features: Saves the current model estimates to either observation
                  (obs==1) or outliers removed (obs==0) for either the
                  hyperbolic (h==1) or linear (h==0) model
        """
        if obs == 1:
            # Backup Model Estimates for Observation Data:
            self.logger.debug("saving observation estimates")
            self.lm_estimates_obs = self.lm_estimates
            self.hm_estimates_obs = self.hm_estimates
        elif obs == 0:
            # Backup Model Estimates for Data with Removed Outliers (RO):
            if h == 1:
                # Outliers Based on Hyperbolic Model
                self.logger.debug("saving hyperbolic model estimates")
                self.hm_estimates_ro = self.hm_estimates
            elif h == 0:
                # Outliers Based on Linear Model
                self.logger.debug("saving linear model estimates")
                self.lm_estimates_ro = self.lm_estimates

    def save_stats(self, obs, h):
        """
        Name:     STAGE1.save_stats
        Inputs:   - int, observation flag (obs)
                  - int, hyperbolic model flag (h)
        Output:   None.
        Features: Saves the current statistical parameters to either
                  observation (obs==1) or outliers removed (obs==0) for either
                  the hyperbolic (h==1) or linear (h==0) model
        """
        if obs == 1:
            # Backup Observation Data:
            self.logger.debug("saving stats for observations")
            self.max_ppfd_obs = self.max_ppfd
            self.min_ppfd_obs = self.min_ppfd
            self.ave_ppfd_obs = self.ave_ppfd
            self.std_ppfd_obs = self.std_ppfd
            self.skw_ppfd_obs = self.skw_ppfd
            self.krt_ppfd_obs = self.krt_ppfd
            self.max_nee_obs = self.max_nee
            self.min_nee_obs = self.min_nee
            self.ave_nee_obs = self.ave_nee
            self.std_nee_obs = self.std_nee
            self.skw_nee_obs = self.skw_nee
            self.krt_nee_obs = self.krt_nee
        elif obs == 0:
            # Backup Data with Removed Outliers (RO):
            if h:
                # Outliers Based on Hyperbolic Model
                self.logger.debug("saving stats for hyperbolic regression")
                self.max_ppfd_ro_h = self.max_ppfd
                self.min_ppfd_ro_h = self.min_ppfd
                self.ave_ppfd_ro_h = self.ave_ppfd
                self.std_ppfd_ro_h = self.std_ppfd
                self.skw_ppfd_ro_h = self.skw_ppfd
                self.krt_ppfd_ro_h = self.krt_ppfd
                self.max_nee_ro_h = self.max_nee
                self.min_nee_ro_h = self.min_nee
                self.ave_nee_ro_h = self.ave_nee
                self.std_nee_ro_h = self.std_nee
                self.skw_nee_ro_h = self.skw_nee
                self.krt_nee_ro_h = self.krt_nee
            else:
                # Outliers Based on Linear Model
                self.logger.debug("saving stats for linear regression")
                self.max_ppfd_ro_l = self.max_ppfd
                self.min_ppfd_ro_l = self.min_ppfd
                self.ave_ppfd_ro_l = self.ave_ppfd
                self.std_ppfd_ro_l = self.std_ppfd
                self.skw_ppfd_ro_l = self.skw_ppfd
                self.krt_ppfd_ro_l = self.krt_ppfd
                self.max_nee_ro_l = self.max_nee
                self.min_nee_ro_l = self.min_nee
                self.ave_nee_ro_l = self.ave_nee
                self.std_nee_ro_l = self.std_nee
                self.skw_nee_ro_l = self.skw_nee
                self.krt_nee_ro_l = self.krt_nee

    def summary_statistics(self):
        """
        Name:     STAGE1.summary_statistics
        Input:    None.
        Output:   string, summary of statistics (sum_stat)
        Features: Returns output string that summarizes all variables in this
                  class
        """
        # Fields:
        # 01. name .......... :: station name
        # 02. month ......... :: year, month, day starting the analysis period
        # 03. n_obs ......... :: number of observations
        # 04. n_h ........... :: outliers identified in model h
        # 05. n_l ........... :: outliers identified in model l
        # 06. foo_est_obs_h   :: Foo estimate for observations for model H
        # 07. foo_opt_obs_h   :: Foo optimized for observations for model H
        # 08. foo_err_obs_h   :: Foo std error for observations for model H
        # 09. foo_t_obs_h     :: Foo t-statistic for observations for model H
        # 10. foo_p_obs_h     :: Foo p-value for observations for model H
        # 11. foo_est_ro_h    :: Foo estimate for obs w/o outliers model H
        # 12. foo_opt_ro_h .. :: Foo optimized for observations w/o outliers
        #                        for model H
        # 13. foo_err_ro_h .. :: Foo standard error for observations w/o
        #                        outliers for model H
        # 14. foo_t_ro_h      :: Foo t-value for obs w/o outliers for model H
        # 15. foo_p_ro_h      :: Foo p-value for obs w/o outliers for model H
        # 16. alpha_est_obs_h :: alpha estimate for observations for model H
        # 17. alpha_opt_obs_h :: alpha optimized for observations for model H
        # 18. alpha_err_obs_h :: alpha standard error for observations for
        #                        model H
        # 19. alpha_t_obs_h   :: alpha t-value for observations for model L
        # 20. alpha_p_obs_h   :: alpha p-value for observations for model L
        # 21. alpha_est_ro_h  :: alpha estimate for obs w/o outliers model H
        # 22. alpha_opt_ro_h  :: alpha optimized for observations w/o outliers
        #                        for model H
        # 23. alpha_err_ro_h  :: alpha standard error for observations w/o
        #                        outliers for model H
        # 24. alpha_t_ro_h    :: alpha t-value for observations for model L
        # 25. alpha_p_ro_h    :: alpha p-value for observations for model L
        # 26. alpha_est_obs_l :: alpha estimate for observations for model L
        # 27. alpha_opt_obs_l :: alpha optimized for observations for model L
        # 28. alpha_err_obs_l :: alpha standard error for observations for
        #                        model L
        # 29. alpha_t_obs_l   :: alpha t-value for observations for model L
        # 30. alpha_p_obs_l   :: alpha p-value for observations for model L
        # 31. alpha_est_ro_l  :: alpha estimate for obs w/o outliers model L
        # 32. alpha_opt_ro_l  :: alpha optimized for observations w/o outliers
        #                        for model L
        # 33. alpha_err_ro_l  :: alpha standard error for observations w/o
        #                        outliers for model L
        # 34. alpha_t_ro_l    :: alpha t-value for obs w/o outliers model L
        # 35. alpha_p_ro_l    :: alpha p-value for obs w/o outliers model L
        # 36. r_est_obs_h ... :: R estimate for observations for model H
        # 37. r_opt_obs_h ... :: R optimized for observations for model H
        # 38. r_err_obs_h ... :: R standard error for observations for model H
        # 39. r_t_obs_h ..... :: R t-value for observations for model H
        # 40. r_p_obs_h ..... :: R p-value for observations for model H
        # 41. r_est_ro_h .... :: R estimate for obs w/o outliers model H
        # 42. r_opt_ro_h .... :: R optimized for observations w/o outliers for
        #                        model H
        # 43. r_err_ro_h .... :: R standard error for observations w/o outliers
        #                        for model H
        # 44. r_t_ro_h ...... :: R t-value for obs w/o outliers for model H
        # 45. r_p_ro_h ...... :: R p-value for obs w/o outliers for model H
        # 46. r_est_obs_l ... :: R estimate for observations for model L
        # 47. r_opt_obs_l ... :: R optimized for observations for model L
        # 48. r_err_obs_l ... :: R standard error for observations for model L
        # 49. r_t_obs_l ..... :: R t-value for observations for model L
        # 50. r_p_obs_l ..... :: R p-value for observations for model L
        # 51. r_est_ro_l .... :: R estimate for obs w/o outliers model L
        # 52. r_opt_ro_l .... :: R optimized for observations w/o outliers for
        #                        model L
        # 53. r_err_ro_l .... :: R standard error for observations w/o outliers
        #                        for model L
        # 54. r_t_ro_l ...... :: R t-value for obs w/o outliers for model L
        # 55. r_p_ro_l ...... :: R p-value for obs w/o outliers for model L
        # 56. r2_obs_h ...... :: r-squared for observations for model H
        # 57. r2_ro_h ....... :: r-squared for observations w/o outliers for
        #                        model H
        # 58. rmse_obs_h .... :: RMSE for observations for model H
        # 59. rmse_ro_h ..... :: RMSE for observations w/o outliers for model H
        # 60. r2_obs_l ...... :: r-squared for observations for model L
        # 61. r2_ro_l ....... :: r-squared for observations w/o outliers for
        #                        model L
        # 62. rmse_obs_l .... :: RMSE for observations for model L
        # 63. rmse_ro_l ..... :: RMSE for observations w/o outliers for model L
        # 64. min_ppfd_obs .. :: minimum of PPFD observations
        # 65. max_ppfd_obs .. :: maximum of PPFD observations
        # 66. ave_ppfd_obs .. :: average of PPFD observations
        # 67. std_ppfd_obs .. :: standard deviation of PPFD observations
        # 68. skw_ppfd_obs .. :: skew of PPFD observations
        # 69. krt_ppfd_obs .. :: kurtosis of PPFD observations
        # 70. min_ppfd_ro_h   :: minimum of PPFD obs w/o outliers (model H)
        # 71. max_ppfd_ro_h   :: maximum of PPFD obs w/o outliers (model H)
        # 72. ave_ppfd_ro_h   :: average of PPFD obs w/o outliers (model H)
        # 73. std_ppfd_ro_h   :: standard deviation of PPFD obs
        # 74. skw_ppfd_ro_h   :: skew of PPFD obs w/o outliers (model H)
        # 75. krt_ppfd_ro_h   :: kurtosis of PPFD obs w/o outliers (model H)
        # 76. min_ppfd_ro_l   :: minimum of PPFD obs w/o outliers (model L)
        # 77. max_ppfd_ro_l   :: maximum of PPFD obs w/o outliers (model L)
        # 78. ave_ppfd_ro_l   :: average of PPFD obs w/o outliers (model L)
        # 79. std_ppfd_ro_l   :: standard deviation of PPFD obs
        # 80. skw_ppfd_ro_l   :: skew of PPFD obs w/o outliers (model L)
        # 81. krt_ppfd_ro_l   :: kurtosis of PPFD obs w/o outliers (model L)
        # 82. min_nee_obs ... :: minimum of NEE observations
        # 83. max_nee_obs ... :: maximum of NEE observations
        # 84. ave_nee_obs ... :: average of NEE observations
        # 85. std_nee_obs ... :: standard deviation of NEE observations
        # 86. skw_nee_obs ... :: skew of NEE observations
        # 87. krt_nee_obs ... :: kurtosis of NEE observations
        # 88. min_nee_ro_h .. :: minimum of NEE obs w/o outliers (model H)
        # 89. max_nee_ro_h .. :: maximum of NEE obs w/o outliers (model H)
        # 90. ave_nee_ro_h .. :: average of NEE obs w/o outliers (model H)
        # 91. std_nee_ro_h .. :: standard deviation of NEE obs
        # 92. skw_nee_ro_h .. :: skew of NEE obs w/o outliers (model H)
        # 93. krt_nee_ro_h .. :: kurtosis of NEE obs w/o outliers (model H)
        # 94. min_nee_ro_l .. :: minimum of NEE obs w/o outliers (model L)
        # 95. max_nee_ro_l .. :: maximum of NEE obs w/o outliers (model L)
        # 96. ave_nee_ro_l .. :: average of NEE obs w/o outliers (model L)
        # 97. std_nee_ro_l .. :: standard deviation of NEE obs
        # 98. skw_nee_ro_l .. :: skew of NEE obs w/o outliers (model L)
        # 99. krt_nee_ro_l .. :: kurtosis of NEE obs w/o outliers (model L)
        #100. pearsonr_obs .. :: pearson's r NEE & PPFD obs
        #101. pearsonr_ro_h . :: pearson's r NEE & PPFD w/o outliers (model H)
        #102. pearsonr_ro_l . :: pearson's r NEE & PPFD w/o outliers (model L)
        #103. mod_select .... :: the best model (0--4)
        #
        sum_stat = (
            "%s,%s,%d,%d,%d,"                             # 01--05
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 06--10
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 11--15
            "%f,%f,%f,%f,%f,"                             # 16--20
            "%f,%f,%f,%f,%f,"                             # 21--25
            "%f,%f,%f,%f,%f,"                             # 26--30
            "%f,%f,%f,%f,%f,"                             # 31--35
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 36--40
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 41--45
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 46--50
            "%0.5f,%0.5f,%0.5f,%f,%f,"                    # 51--55
            "%0.3f,%0.3f,%0.2f,%0.2f,"                    # 56--59
            "%0.3f,%0.3f,%0.2f,%0.2f,"                    # 60--63
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 64--69
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 70--75
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 76--81
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 82--87
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 88--93
            "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,"        # 94--99
            "%0.3f,%0.3f,%0.3f,%d\n"                      # 100--103
            ) % (
            self.name, self.month,                                   # 01--02
            len(self.nee_obs), self.hm_outliers, self.lm_outliers,   # 03--05
            self.hm_estimates_obs[0], self.hm_optimized[0],          # 06--07
            self.hm_optim_err[0],     self.hm_optim_t[0],            # 08--09
            self.hm_optim_p[0],       self.hm_estimates_ro[0],       # 10--11
            self.hm_optimized_ro[0],  self.hm_optim_err_ro[0],       # 12--13
            self.hm_optim_t_ro[0],    self.hm_optim_p_ro[0],         # 14--15
            self.hm_estimates_obs[1], self.hm_optimized[1],          # 16--17
            self.hm_optim_err[1],     self.hm_optim_t[1],            # 18--19
            self.hm_optim_p[1],       self.hm_estimates_ro[1],       # 20--21
            self.hm_optimized_ro[1],  self.hm_optim_err_ro[1],       # 22--23
            self.hm_optim_t_ro[1],    self.hm_optim_p_ro[1],         # 24--25
            self.lm_estimates_obs[0], self.lm_optimized[0],          # 26--27
            self.lm_optim_err[0],     self.lm_optim_t[0],            # 28--29
            self.lm_optim_p[0],       self.lm_estimates_ro[0],       # 30--31
            self.lm_optimized_ro[0],  self.lm_optim_err_ro[0],       # 32--33
            self.lm_optim_t_ro[0],    self.lm_optim_p_ro[0],         # 34--35
            self.hm_estimates_obs[2], self.hm_optimized[2],          # 36--37
            self.hm_optim_err[2],     self.hm_optim_t[2],            # 38--39
            self.hm_optim_p[2],       self.hm_estimates_ro[2],       # 40--41
            self.hm_optimized_ro[2],  self.hm_optim_err_ro[2],       # 42--43
            self.hm_optim_t_ro[2],    self.hm_optim_p_ro[2],         # 44--45
            self.lm_estimates_obs[1], self.lm_optimized[1],          # 46--47
            self.lm_optim_err[1],     self.lm_optim_t[1],            # 48--49
            self.lm_optim_p[1],       self.lm_estimates_ro[1],       # 50--51
            self.lm_optimized_ro[1],  self.lm_optim_err_ro[1],       # 52--53
            self.lm_optim_t_ro[1],    self.lm_optim_p_ro[1],         # 54--55
            self.hm_rsqr,             self.hm_rsqr_ro,               # 56--57
            self.hm_rmse,             self.hm_rmse_ro,               # 58--59
            self.lm_rsqr,             self.lm_rsqr_ro,               # 60--61
            self.lm_rmse,             self.lm_rmse_ro,               # 62--63
            self.min_ppfd_obs,        self.max_ppfd_obs,             # 64--65
            self.ave_ppfd_obs,        self.std_ppfd_obs,             # 66--67
            self.skw_ppfd_obs,        self.krt_ppfd_obs,             # 68--69
            self.min_ppfd_ro_h,       self.max_ppfd_ro_h,            # 70--71
            self.ave_ppfd_ro_h,       self.std_ppfd_ro_h,            # 72--73
            self.skw_ppfd_ro_h,       self.krt_ppfd_ro_h,            # 74--75
            self.min_ppfd_ro_l,       self.max_ppfd_ro_l,            # 76--77
            self.ave_ppfd_ro_l,       self.std_ppfd_ro_l,            # 78--79
            self.skw_ppfd_ro_l,       self.krt_ppfd_ro_l,            # 80--81
            self.min_nee_obs,         self.max_nee_obs,              # 82--83
            self.ave_nee_obs,         self.std_nee_obs,              # 84--85
            self.skw_nee_obs,         self.krt_nee_obs,              # 86--87
            self.min_nee_ro_h,        self.max_nee_ro_h,             # 88--89
            self.ave_nee_ro_h,        self.std_nee_ro_h,             # 90--91
            self.skw_nee_ro_h,        self.krt_nee_ro_h,             # 92--93
            self.min_nee_ro_l,        self.max_nee_ro_l,             # 94--95
            self.ave_nee_ro_l,        self.std_nee_ro_l,             # 96--97
            self.skw_nee_ro_l,        self.krt_nee_ro_l,             # 98--99
            self.r_obs,               self.r_ro_h,                  # 100--101
            self.r_ro_l,              self.mod_select              # 102--103
            )

        return sum_stat

    def update_guess(self):
        """
        Name:     STAGE1.update_guess
        Input:    None.
        Output:   None.
        Features: Updates the model parameter estimation values based on the
                  current data statistics
        """
        self.logger.debug("updating initial model guess values")

        # Linear model [alpha, R]:
        try:
            self.lm_estimates[0] = (
                0.672*(self.max_nee - self.min_nee) /
                (self.max_ppfd - self.min_ppfd)
                )
        except ZeroDivisionError:
            self.lm_estimates[0] = 0.672*(self.max_nee - self.min_nee)/1.0e-3
        else:
            if abs(self.lm_estimates[0]) <= 5.0e-4:
                self.lm_estimates[0] = 1.0e-3

        self.lm_estimates[1] = (
            0.899*self.std_nee
            + 0.827*self.ave_nee
            + 0.00628*self.ave_ppfd
            - 0.008*self.std_ppfd
        )

        # Hyperbolic model [Foo, alpha, R]:
        self.hm_estimates[0] = 3.83*self.std_nee
        try:
            self.hm_estimates[1] = (
                1.96*(self.max_nee - self.min_nee) /
                (self.max_ppfd - self.min_ppfd)
                )
        except ZeroDivisionError:
            self.hm_estimates[1] = 1.96*(self.max_nee - self.min_nee)/1.0e-3
        else:
            if abs(self.hm_estimates[1]) <= 5.0e-4:
                self.hm_estimates[1] = 1.0e-3
        self.hm_estimates[2] = 0.69*self.std_nee

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("stage1.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)
