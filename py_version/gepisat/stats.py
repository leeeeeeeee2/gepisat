#!/usr/bin/python
#
# stats.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-01-24
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
import scipy
from scipy.optimize import curve_fit

from .utilities import calc_statistics
from .utilities import get_hm_estimates
from .utilities import get_lm_estimates
from .utilities import get_outliers
from .utilities import goodness_of_fit
from .utilities import hyp_model
from .utilities import init_summary_dict
from .utilities import lin_model
from .utilities import pearsons_r


###############################################################################
# CLASSES
###############################################################################
class PARTI_STATS(object):
    """
    Name:     PARTI_STATS
    Features: This class holds all the statistics values for FLUX_PARTI
              (see init_summary_dict in utilities for complementary list)
    History   Version 3.0.0-dev
              - created [16.05.20]
              - moved to its own module [16.07.22]
              - created summary string property [16.07.22]
              - moved to gepisat package [16.07.22]
              - updated write fit params function [17.01.22]
              - fixed header items in write fit params function [17.01.24]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     PARTI_STATS.__init__
        Inputs:   None.
        Features: Initializes the class
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.debug("PARTI_STATS class called")

        # Initialize empty stats dictionary:
        self.summary = init_summary_dict()

        # Initialize empty data arrays:
        self.nee = numpy.array([])
        self.ppfd = numpy.array([])

        # Initialize outlier index tuples:
        self.hyp_outliers = (numpy.array([]), numpy.array([]))
        self.lin_outliers = (numpy.array([]), numpy.array([]))

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Property Definitions
    # ////////////////////////////////////////////////////////////////////////
    @property
    def name(self):
        """Station name"""
        k = self.get_summary_idx("name")
        return self.summary[k]['val']

    @name.setter
    def name(self, val):
        k = self.get_summary_idx("name")
        self.summary[k]['val'] = val

    @property
    def date(self):
        """Date for station data (i.e., month)"""
        k = self.get_summary_idx("month")
        return self.summary[k]['val']

    @date.setter
    def date(self, val):
        k = self.get_summary_idx("month")
        if isinstance(val, datetime.date):
            self.summary[k]['val'] = val
        else:
            raise TypeError("Date attribute must be datetime date type")

    @property
    def num_obs(self):
        """Number of observations"""
        k = self.get_summary_idx("n_obs")
        return self.summary[k]['val']

    @num_obs.setter
    def num_obs(self, val):
        k = self.get_summary_idx("n_obs")
        self.summary[k]['val'] = int(val)

    @property
    def nh_outliers(self):
        """Number of outliers from the hyperbolic model"""
        k = self.get_summary_idx("n_h")
        return self.summary[k]['val']

    @nh_outliers.setter
    def nh_outliers(self, val):
        k = self.get_summary_idx("n_h")
        self.summary[k]['val'] = int(val)

    @property
    def nl_outliers(self):
        """Number of outliers from the linear model"""
        k = self.get_summary_idx("n_l")
        return self.summary[k]['val']

    @nl_outliers.setter
    def nl_outliers(self, val):
        k = self.get_summary_idx("n_l")
        self.summary[k]['val'] = int(val)

    @property
    def max_ppfd_obs(self):
        """Maximum PPFD observation"""
        k = self.get_summary_idx("max_ppfd_obs")
        return self.summary[k]['val']

    @max_ppfd_obs.setter
    def max_ppfd_obs(self, val):
        k = self.get_summary_idx("max_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def min_ppfd_obs(self):
        """Minimum PPFD observation"""
        k = self.get_summary_idx("min_ppfd_obs")
        return self.summary[k]['val']

    @min_ppfd_obs.setter
    def min_ppfd_obs(self, val):
        k = self.get_summary_idx("min_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def ave_ppfd_obs(self):
        """Averge PPFD observation"""
        k = self.get_summary_idx("ave_ppfd_obs")
        return self.summary[k]['val']

    @ave_ppfd_obs.setter
    def ave_ppfd_obs(self, val):
        k = self.get_summary_idx("ave_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def std_ppfd_obs(self):
        """Standard deviation of PPFD observation"""
        k = self.get_summary_idx("std_ppfd_obs")
        return self.summary[k]['val']

    @std_ppfd_obs.setter
    def std_ppfd_obs(self, val):
        k = self.get_summary_idx("std_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def skw_ppfd_obs(self):
        """Skew of PPFD observation"""
        k = self.get_summary_idx("skw_ppfd_obs")
        return self.summary[k]['val']

    @skw_ppfd_obs.setter
    def skw_ppfd_obs(self, val):
        k = self.get_summary_idx("skw_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def krt_ppfd_obs(self):
        """Kurtosis of PPFD observation"""
        k = self.get_summary_idx("krt_ppfd_obs")
        return self.summary[k]['val']

    @krt_ppfd_obs.setter
    def krt_ppfd_obs(self, val):
        k = self.get_summary_idx("krt_ppfd_obs")
        self.summary[k]['val'] = float(val)

    @property
    def min_ppfd_ro_h(self):
        """Minimum of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("min_ppfd_ro_h")
        return self.summary[k]['val']

    @min_ppfd_ro_h.setter
    def min_ppfd_ro_h(self, val):
        k = self.get_summary_idx("min_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def max_ppfd_ro_h(self):
        """Maximum of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("max_ppfd_ro_h")
        return self.summary[k]['val']

    @max_ppfd_ro_h.setter
    def max_ppfd_ro_h(self, val):
        k = self.get_summary_idx("max_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def ave_ppfd_ro_h(self):
        """Average of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("ave_ppfd_ro_h")
        return self.summary[k]['val']

    @ave_ppfd_ro_h.setter
    def ave_ppfd_ro_h(self, val):
        k = self.get_summary_idx("ave_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def std_ppfd_ro_h(self):
        """Std deviations of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("std_ppfd_ro_h")
        return self.summary[k]['val']

    @std_ppfd_ro_h.setter
    def std_ppfd_ro_h(self, val):
        k = self.get_summary_idx("std_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def skw_ppfd_ro_h(self):
        """Skew of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("skw_ppfd_ro_h")
        return self.summary[k]['val']

    @skw_ppfd_ro_h.setter
    def skw_ppfd_ro_h(self, val):
        k = self.get_summary_idx("skw_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def krt_ppfd_ro_h(self):
        """Kurtosis of PPFD with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("krt_ppfd_ro_h")
        return self.summary[k]['val']

    @krt_ppfd_ro_h.setter
    def krt_ppfd_ro_h(self, val):
        k = self.get_summary_idx("krt_ppfd_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def min_ppfd_ro_l(self):
        """Minimum of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("min_ppfd_ro_l")
        return self.summary[k]['val']

    @min_ppfd_ro_l.setter
    def min_ppfd_ro_l(self, val):
        k = self.get_summary_idx("min_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def max_ppfd_ro_l(self):
        """Maximum of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("max_ppfd_ro_l")
        return self.summary[k]['val']

    @max_ppfd_ro_l.setter
    def max_ppfd_ro_l(self, val):
        k = self.get_summary_idx("max_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def ave_ppfd_ro_l(self):
        """Average of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("ave_ppfd_ro_l")
        return self.summary[k]['val']

    @ave_ppfd_ro_l.setter
    def ave_ppfd_ro_l(self, val):
        k = self.get_summary_idx("ave_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def std_ppfd_ro_l(self):
        """Standard deviation of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("std_ppfd_ro_l")
        return self.summary[k]['val']

    @std_ppfd_ro_l.setter
    def std_ppfd_ro_l(self, val):
        k = self.get_summary_idx("std_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def skw_ppfd_ro_l(self):
        """Skew of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("skw_ppfd_ro_l")
        return self.summary[k]['val']

    @skw_ppfd_ro_l.setter
    def skw_ppfd_ro_l(self, val):
        k = self.get_summary_idx("skw_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def krt_ppfd_ro_l(self):
        """Kurtosis of PPFD with outliers removed by linear model"""
        k = self.get_summary_idx("krt_ppfd_ro_l")
        return self.summary[k]['val']

    @krt_ppfd_ro_l.setter
    def krt_ppfd_ro_l(self, val):
        k = self.get_summary_idx("krt_ppfd_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def min_nee_obs(self):
        """Minimum of NEE observations"""
        k = self.get_summary_idx("min_nee_obs")
        return self.summary[k]['val']

    @min_nee_obs.setter
    def min_nee_obs(self, val):
        k = self.get_summary_idx("min_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def max_nee_obs(self):
        """Maximum of NEE observations"""
        k = self.get_summary_idx("max_nee_obs")
        return self.summary[k]['val']

    @max_nee_obs.setter
    def max_nee_obs(self, val):
        k = self.get_summary_idx("max_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def ave_nee_obs(self):
        """Average of NEE observations"""
        k = self.get_summary_idx("ave_nee_obs")
        return self.summary[k]['val']

    @ave_nee_obs.setter
    def ave_nee_obs(self, val):
        k = self.get_summary_idx("ave_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def std_nee_obs(self):
        """Standard deviation of NEE observations"""
        k = self.get_summary_idx("std_nee_obs")
        return self.summary[k]['val']

    @std_nee_obs.setter
    def std_nee_obs(self, val):
        k = self.get_summary_idx("std_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def skw_nee_obs(self):
        """Skew of NEE observations"""
        k = self.get_summary_idx("skw_nee_obs")
        return self.summary[k]['val']

    @skw_nee_obs.setter
    def skw_nee_obs(self, val):
        k = self.get_summary_idx("skw_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def krt_nee_obs(self):
        """Kurtosis of NEE observations"""
        k = self.get_summary_idx("krt_nee_obs")
        return self.summary[k]['val']

    @krt_nee_obs.setter
    def krt_nee_obs(self, val):
        k = self.get_summary_idx("krt_nee_obs")
        self.summary[k]['val'] = float(val)

    @property
    def min_nee_ro_h(self):
        """Minimum of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("min_nee_ro_h")
        return self.summary[k]['val']

    @min_nee_ro_h.setter
    def min_nee_ro_h(self, val):
        k = self.get_summary_idx("min_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def max_nee_ro_h(self):
        """Maximum of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("max_nee_ro_h")
        return self.summary[k]['val']

    @max_nee_ro_h.setter
    def max_nee_ro_h(self, val):
        k = self.get_summary_idx("max_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def ave_nee_ro_h(self):
        """Average of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("ave_nee_ro_h")
        return self.summary[k]['val']

    @ave_nee_ro_h.setter
    def ave_nee_ro_h(self, val):
        k = self.get_summary_idx("ave_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def std_nee_ro_h(self):
        """Std. deviation of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("std_nee_ro_h")
        return self.summary[k]['val']

    @std_nee_ro_h.setter
    def std_nee_ro_h(self, val):
        k = self.get_summary_idx("std_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def skw_nee_ro_h(self):
        """Skew of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("skw_nee_ro_h")
        return self.summary[k]['val']

    @skw_nee_ro_h.setter
    def skw_nee_ro_h(self, val):
        k = self.get_summary_idx("skw_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def krt_nee_ro_h(self):
        """Kurtosis of NEE with outliers removed by hyperbolic model"""
        k = self.get_summary_idx("krt_nee_ro_h")
        return self.summary[k]['val']

    @krt_nee_ro_h.setter
    def krt_nee_ro_h(self, val):
        k = self.get_summary_idx("krt_nee_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def min_nee_ro_l(self):
        """Minimum of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("min_nee_ro_l")
        return self.summary[k]['val']

    @min_nee_ro_l.setter
    def min_nee_ro_l(self, val):
        k = self.get_summary_idx("min_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def max_nee_ro_l(self):
        """Maximum of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("max_nee_ro_l")
        return self.summary[k]['val']

    @max_nee_ro_l.setter
    def max_nee_ro_l(self, val):
        k = self.get_summary_idx("max_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def ave_nee_ro_l(self):
        """Average of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("ave_nee_ro_l")
        return self.summary[k]['val']

    @ave_nee_ro_l.setter
    def ave_nee_ro_l(self, val):
        k = self.get_summary_idx("ave_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def std_nee_ro_l(self):
        """Standard deviation of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("std_nee_ro_l")
        return self.summary[k]['val']

    @std_nee_ro_l.setter
    def std_nee_ro_l(self, val):
        k = self.get_summary_idx("std_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def skw_nee_ro_l(self):
        """Skew of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("skw_nee_ro_l")
        return self.summary[k]['val']

    @skw_nee_ro_l.setter
    def skw_nee_ro_l(self, val):
        k = self.get_summary_idx("skw_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def krt_nee_ro_l(self):
        """Kurtosis of NEE with outliers removed by linear model"""
        k = self.get_summary_idx("krt_nee_ro_l")
        return self.summary[k]['val']

    @krt_nee_ro_l.setter
    def krt_nee_ro_l(self, val):
        k = self.get_summary_idx("krt_nee_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def foo_est_obs_h(self):
        """Estimated value for hyperbolic foo from observations"""
        k = self.get_summary_idx("foo_est_obs_h")
        return self.summary[k]['val']

    @foo_est_obs_h.setter
    def foo_est_obs_h(self, val):
        k = self.get_summary_idx("foo_est_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_opt_obs_h(self):
        """Optimized value of hyperbolic foo from observations"""
        k = self.get_summary_idx("foo_opt_obs_h")
        return self.summary[k]['val']

    @foo_opt_obs_h.setter
    def foo_opt_obs_h(self, val):
        k = self.get_summary_idx("foo_opt_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_err_obs_h(self):
        """Optimized value error of hyperbolic foo from observations"""
        k = self.get_summary_idx("foo_err_obs_h")
        return self.summary[k]['val']

    @foo_err_obs_h.setter
    def foo_err_obs_h(self, val):
        k = self.get_summary_idx("foo_err_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_t_obs_h(self):
        """T-value of optimized hyperbolic foo from observations"""
        k = self.get_summary_idx("foo_t_obs_h")
        return self.summary[k]['val']

    @foo_t_obs_h.setter
    def foo_t_obs_h(self, val):
        k = self.get_summary_idx("foo_t_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_p_obs_h(self):
        """p-value of optimized hyperbolic foo from observations"""
        k = self.get_summary_idx("foo_p_obs_h")
        return self.summary[k]['val']

    @foo_p_obs_h.setter
    def foo_p_obs_h(self, val):
        k = self.get_summary_idx("foo_p_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_est_ro_h(self):
        """Estimated value of hyperbolic foo without outliers"""
        k = self.get_summary_idx("foo_est_ro_h")
        return self.summary[k]['val']

    @foo_est_ro_h.setter
    def foo_est_ro_h(self, val):
        k = self.get_summary_idx("foo_est_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_opt_ro_h(self):
        """Optimized value of hyperbolic foo without outliers"""
        k = self.get_summary_idx("foo_opt_ro_h")
        return self.summary[k]['val']

    @foo_opt_ro_h.setter
    def foo_opt_ro_h(self, val):
        k = self.get_summary_idx("foo_opt_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_err_ro_h(self):
        """Optimized value error of hyperbolic foo without outliers"""
        k = self.get_summary_idx("foo_err_ro_h")
        return self.summary[k]['val']

    @foo_err_ro_h.setter
    def foo_err_ro_h(self, val):
        k = self.get_summary_idx("foo_err_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_t_ro_h(self):
        """T-value of optimized hyperbolic foo without outliers"""
        k = self.get_summary_idx("foo_t_ro_h")
        return self.summary[k]['val']

    @foo_t_ro_h.setter
    def foo_t_ro_h(self, val):
        k = self.get_summary_idx("foo_t_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def foo_p_ro_h(self):
        """p-value of optimized hyperbolic foo without outliers"""
        k = self.get_summary_idx("foo_p_ro_h")
        return self.summary[k]['val']

    @foo_p_ro_h.setter
    def foo_p_ro_h(self, val):
        k = self.get_summary_idx("foo_p_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_est_obs_h(self):
        """Estimated value of hyperbolic alpha from observations"""
        k = self.get_summary_idx("alpha_est_obs_h")
        return self.summary[k]['val']

    @alpha_est_obs_h.setter
    def alpha_est_obs_h(self, val):
        k = self.get_summary_idx("alpha_est_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_opt_obs_h(self):
        """Optimized value of hyperbolic alpha from observations"""
        k = self.get_summary_idx("alpha_opt_obs_h")
        return self.summary[k]['val']

    @alpha_opt_obs_h.setter
    def alpha_opt_obs_h(self, val):
        k = self.get_summary_idx("alpha_opt_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_err_obs_h(self):
        """Optimized value error of hyperbolic alpha from observations"""
        k = self.get_summary_idx("alpha_err_obs_h")
        return self.summary[k]['val']

    @alpha_err_obs_h.setter
    def alpha_err_obs_h(self, val):
        k = self.get_summary_idx("alpha_err_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_t_obs_h(self):
        """T-value of optimized hyperbolic alpha from observations"""
        k = self.get_summary_idx("alpha_t_obs_h")
        return self.summary[k]['val']

    @alpha_t_obs_h.setter
    def alpha_t_obs_h(self, val):
        k = self.get_summary_idx("alpha_t_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_p_obs_h(self):
        """p-value of optimized hyperbolic alpha from observations"""
        k = self.get_summary_idx("alpha_p_obs_h")
        return self.summary[k]['val']

    @alpha_p_obs_h.setter
    def alpha_p_obs_h(self, val):
        k = self.get_summary_idx("alpha_p_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_est_ro_h(self):
        """Estimated value of hyperbolic alpha without outliers"""
        k = self.get_summary_idx("alpha_est_ro_h")
        return self.summary[k]['val']

    @alpha_est_ro_h.setter
    def alpha_est_ro_h(self, val):
        k = self.get_summary_idx("alpha_est_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_opt_ro_h(self):
        """Optimized value of hyperbolic alpha without outliers"""
        k = self.get_summary_idx("alpha_opt_ro_h")
        return self.summary[k]['val']

    @alpha_opt_ro_h.setter
    def alpha_opt_ro_h(self, val):
        k = self.get_summary_idx("alpha_opt_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_err_ro_h(self):
        """Optimized value error of hyperbolic alpha without outliers"""
        k = self.get_summary_idx("alpha_err_ro_h")
        return self.summary[k]['val']

    @alpha_err_ro_h.setter
    def alpha_err_ro_h(self, val):
        k = self.get_summary_idx("alpha_err_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_t_ro_h(self):
        """T-value of hyperbolic alpha without outliers"""
        k = self.get_summary_idx("alpha_t_ro_h")
        return self.summary[k]['val']

    @alpha_t_ro_h.setter
    def alpha_t_ro_h(self, val):
        k = self.get_summary_idx("alpha_t_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_p_ro_h(self):
        """p-value of hyperbolic alpha without outliers"""
        k = self.get_summary_idx("alpha_p_ro_h")
        return self.summary[k]['val']

    @alpha_p_ro_h.setter
    def alpha_p_ro_h(self, val):
        k = self.get_summary_idx("alpha_p_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_est_obs_l(self):
        """Estimated value of linear alpha from observations"""
        k = self.get_summary_idx("alpha_est_obs_l")
        return self.summary[k]['val']

    @alpha_est_obs_l.setter
    def alpha_est_obs_l(self, val):
        k = self.get_summary_idx("alpha_est_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_opt_obs_l(self):
        """Optimized value of linear alpha from observations"""
        k = self.get_summary_idx("alpha_opt_obs_l")
        return self.summary[k]['val']

    @alpha_opt_obs_l.setter
    def alpha_opt_obs_l(self, val):
        k = self.get_summary_idx("alpha_opt_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_err_obs_l(self):
        """Optimized value error of linear alpha from observations"""
        k = self.get_summary_idx("alpha_err_obs_l")
        return self.summary[k]['val']

    @alpha_err_obs_l.setter
    def alpha_err_obs_l(self, val):
        k = self.get_summary_idx("alpha_err_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_t_obs_l(self):
        """T-value of linear alpha from observations"""
        k = self.get_summary_idx("alpha_t_obs_l")
        return self.summary[k]['val']

    @alpha_t_obs_l.setter
    def alpha_t_obs_l(self, val):
        k = self.get_summary_idx("alpha_t_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_p_obs_l(self):
        """p-value of linear alpha from observations"""
        k = self.get_summary_idx("alpha_p_obs_l")
        return self.summary[k]['val']

    @alpha_p_obs_l.setter
    def alpha_p_obs_l(self, val):
        k = self.get_summary_idx("alpha_p_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_est_ro_l(self):
        """Estimated value of linear alpha without outliers"""
        k = self.get_summary_idx("alpha_est_ro_l")
        return self.summary[k]['val']

    @alpha_est_ro_l.setter
    def alpha_est_ro_l(self, val):
        k = self.get_summary_idx("alpha_est_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_opt_ro_l(self):
        """Optimized value of linear alpha without outliers"""
        k = self.get_summary_idx("alpha_opt_ro_l")
        return self.summary[k]['val']

    @alpha_opt_ro_l.setter
    def alpha_opt_ro_l(self, val):
        k = self.get_summary_idx("alpha_opt_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_err_ro_l(self):
        """Optimized value error of linear alpha without outliers"""
        k = self.get_summary_idx("alpha_err_ro_l")
        return self.summary[k]['val']

    @alpha_err_ro_l.setter
    def alpha_err_ro_l(self, val):
        k = self.get_summary_idx("alpha_err_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_t_ro_l(self):
        """T-value of optimized linear alpha without outliers"""
        k = self.get_summary_idx("alpha_t_ro_l")
        return self.summary[k]['val']

    @alpha_t_ro_l.setter
    def alpha_t_ro_l(self, val):
        k = self.get_summary_idx("alpha_t_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def alpha_p_ro_l(self):
        """p-value of optimized linear alpha without outliers"""
        k = self.get_summary_idx("alpha_p_ro_l")
        return self.summary[k]['val']

    @alpha_p_ro_l.setter
    def alpha_p_ro_l(self, val):
        k = self.get_summary_idx("alpha_p_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_est_obs_h(self):
        """Estimated value of hyperbolic r from observations"""
        k = self.get_summary_idx("r_est_obs_h")
        return self.summary[k]['val']

    @r_est_obs_h.setter
    def r_est_obs_h(self, val):
        k = self.get_summary_idx("r_est_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_opt_obs_h(self):
        """Optimized value of hyperbolic r from observations"""
        k = self.get_summary_idx("r_opt_obs_h")
        return self.summary[k]['val']

    @r_opt_obs_h.setter
    def r_opt_obs_h(self, val):
        k = self.get_summary_idx("r_opt_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_err_obs_h(self):
        """Optimized value error of hyperbolic r from observations"""
        k = self.get_summary_idx("r_err_obs_h")
        return self.summary[k]['val']

    @r_err_obs_h.setter
    def r_err_obs_h(self, val):
        k = self.get_summary_idx("r_err_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_t_obs_h(self):
        """T-value of optimized hyperbolic r from observations"""
        k = self.get_summary_idx("r_t_obs_h")
        return self.summary[k]['val']

    @r_t_obs_h.setter
    def r_t_obs_h(self, val):
        k = self.get_summary_idx("r_t_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_p_obs_h(self):
        """p-value of optimized hyperbolic r from observations"""
        k = self.get_summary_idx("r_p_obs_h")
        return self.summary[k]['val']

    @r_p_obs_h.setter
    def r_p_obs_h(self, val):
        k = self.get_summary_idx("r_p_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_est_ro_h(self):
        """Estimated value of hyperbolic r without outliers"""
        k = self.get_summary_idx("r_est_ro_h")
        return self.summary[k]['val']

    @r_est_ro_h.setter
    def r_est_ro_h(self, val):
        k = self.get_summary_idx("r_est_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_opt_ro_h(self):
        """Optimized value of hyperbolic r without outliers"""
        k = self.get_summary_idx("r_opt_ro_h")
        return self.summary[k]['val']

    @r_opt_ro_h.setter
    def r_opt_ro_h(self, val):
        k = self.get_summary_idx("r_opt_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_err_ro_h(self):
        """Optimized value error of hyperbolic r without outliers"""
        k = self.get_summary_idx("r_err_ro_h")
        return self.summary[k]['val']

    @r_err_ro_h.setter
    def r_err_ro_h(self, val):
        k = self.get_summary_idx("r_err_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_t_ro_h(self):
        """T-value of optimized hyperbolic r without outliers"""
        k = self.get_summary_idx("r_t_ro_h")
        return self.summary[k]['val']

    @r_t_ro_h.setter
    def r_t_ro_h(self, val):
        k = self.get_summary_idx("r_t_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_p_ro_h(self):
        """p-value of optimized hyperbolic r without outliers"""
        k = self.get_summary_idx("r_p_ro_h")
        return self.summary[k]['val']

    @r_p_ro_h.setter
    def r_p_ro_h(self, val):
        k = self.get_summary_idx("r_p_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r_est_obs_l(self):
        """Estimated value of linear r from observations"""
        k = self.get_summary_idx("r_est_obs_l")
        return self.summary[k]['val']

    @r_est_obs_l.setter
    def r_est_obs_l(self, val):
        k = self.get_summary_idx("r_est_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_opt_obs_l(self):
        """Optimized value of linear r from observations"""
        k = self.get_summary_idx("r_opt_obs_l")
        return self.summary[k]['val']

    @r_opt_obs_l.setter
    def r_opt_obs_l(self, val):
        k = self.get_summary_idx("r_opt_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_err_obs_l(self):
        """Optimized value error of linear r from observations"""
        k = self.get_summary_idx("r_err_obs_l")
        return self.summary[k]['val']

    @r_err_obs_l.setter
    def r_err_obs_l(self, val):
        k = self.get_summary_idx("r_err_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_t_obs_l(self):
        """T-value of optimized linear r from observations"""
        k = self.get_summary_idx("r_t_obs_l")
        return self.summary[k]['val']

    @r_t_obs_l.setter
    def r_t_obs_l(self, val):
        k = self.get_summary_idx("r_t_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_p_obs_l(self):
        """p-value of optimized linear r from observations"""
        k = self.get_summary_idx("r_p_obs_l")
        return self.summary[k]['val']

    @r_p_obs_l.setter
    def r_p_obs_l(self, val):
        k = self.get_summary_idx("r_p_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_est_ro_l(self):
        """Estimated value of linear r without outliers"""
        k = self.get_summary_idx("r_est_ro_l")
        return self.summary[k]['val']

    @r_est_ro_l.setter
    def r_est_ro_l(self, val):
        k = self.get_summary_idx("r_est_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_opt_ro_l(self):
        """Optimized value of linear r without outliers"""
        k = self.get_summary_idx("r_opt_ro_l")
        return self.summary[k]['val']

    @r_opt_ro_l.setter
    def r_opt_ro_l(self, val):
        k = self.get_summary_idx("r_opt_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_err_ro_l(self):
        """Optimized value error of linear r without outliers"""
        k = self.get_summary_idx("r_err_ro_l")
        return self.summary[k]['val']

    @r_err_ro_l.setter
    def r_err_ro_l(self, val):
        k = self.get_summary_idx("r_err_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_t_ro_l(self):
        """T-value of optimized linear r without outliers"""
        k = self.get_summary_idx("r_t_ro_l")
        return self.summary[k]['val']

    @r_t_ro_l.setter
    def r_t_ro_l(self, val):
        k = self.get_summary_idx("r_t_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r_p_ro_l(self):
        """p-value of optimized linear r without outliers"""
        k = self.get_summary_idx("r_p_ro_l")
        return self.summary[k]['val']

    @r_p_ro_l.setter
    def r_p_ro_l(self, val):
        k = self.get_summary_idx("r_p_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def pearson_r_obs(self):
        """Pearson's correlation coefficient for observations"""
        k = self.get_summary_idx("pearson_r_obs")
        return self.summary[k]['val']

    @pearson_r_obs.setter
    def pearson_r_obs(self, val):
        k = self.get_summary_idx("pearson_r_obs")
        self.summary[k]['val'] = float(val)

    @property
    def pearson_r_ro_h(self):
        """Pearson's correlation coefficient without hyperbolic outliers"""
        k = self.get_summary_idx("pearson_r_ro_h")
        return self.summary[k]['val']

    @pearson_r_ro_h.setter
    def pearson_r_ro_h(self, val):
        k = self.get_summary_idx("pearson_r_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def pearson_r_ro_l(self):
        """Pearson's correlation coefficient without linear outliers"""
        k = self.get_summary_idx("pearson_r_ro_l")
        return self.summary[k]['val']

    @pearson_r_ro_l.setter
    def pearson_r_ro_l(self, val):
        k = self.get_summary_idx("pearson_r_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def r2_obs_h(self):
        """Hyperbolic model coefficient of determination for observations"""
        k = self.get_summary_idx("r2_obs_h")
        return self.summary[k]['val']

    @r2_obs_h.setter
    def r2_obs_h(self, val):
        k = self.get_summary_idx("r2_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def r2_ro_h(self):
        """Hyperbolic model coefficient of determination without outliers"""
        k = self.get_summary_idx("r2_ro_h")
        return self.summary[k]['val']

    @r2_ro_h.setter
    def r2_ro_h(self, val):
        k = self.get_summary_idx("r2_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def rmse_obs_h(self):
        """Hyperbolic model root-mean-squared-error for observations"""
        k = self.get_summary_idx("rmse_obs_h")
        return self.summary[k]['val']

    @rmse_obs_h.setter
    def rmse_obs_h(self, val):
        k = self.get_summary_idx("rmse_obs_h")
        self.summary[k]['val'] = float(val)

    @property
    def rmse_ro_h(self):
        """Hyperbolic model root-mean-squared-error without outliers"""
        k = self.get_summary_idx("rmse_ro_h")
        return self.summary[k]['val']

    @rmse_ro_h.setter
    def rmse_ro_h(self, val):
        k = self.get_summary_idx("rmse_ro_h")
        self.summary[k]['val'] = float(val)

    @property
    def r2_obs_l(self):
        """Linear model coefficient of determination for observations"""
        k = self.get_summary_idx("r2_obs_l")
        return self.summary[k]['val']

    @r2_obs_l.setter
    def r2_obs_l(self, val):
        k = self.get_summary_idx("r2_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def r2_ro_l(self):
        """Linear model coefficient of determination without outliers"""
        k = self.get_summary_idx("r2_ro_l")
        return self.summary[k]['val']

    @r2_ro_l.setter
    def r2_ro_l(self, val):
        k = self.get_summary_idx("r2_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def rmse_obs_l(self):
        """Linear model root-mean-squared-error for observations"""
        k = self.get_summary_idx("rmse_obs_l")
        return self.summary[k]['val']

    @rmse_obs_l.setter
    def rmse_obs_l(self, val):
        k = self.get_summary_idx("rmse_obs_l")
        self.summary[k]['val'] = float(val)

    @property
    def rmse_ro_l(self):
        """Linear model root-mean-squared-error without outliers"""
        k = self.get_summary_idx("rmse_ro_l")
        return self.summary[k]['val']

    @rmse_ro_l.setter
    def rmse_ro_l(self, val):
        k = self.get_summary_idx("rmse_ro_l")
        self.summary[k]['val'] = float(val)

    @property
    def model_select(self):
        """Best model index based on fitting"""
        k = self.get_summary_idx("model_select")
        return self.summary[k]['val']

    @model_select.setter
    def model_select(self, val):
        k = self.get_summary_idx("model_select")
        self.summary[k]['val'] = val

    @property
    def summary_string(self):
        """String of full summary statistics"""
        sstring = ""
        for k in sorted(list(self.summary.keys())):
            val = self.summary[k]['val']
            sstring += "{},".format(val)
        sstring = sstring.rstrip(",")
        sstring += "\n"
        return sstring

    @summary_string.setter
    def summary_string(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def find_outliers(self):
        """
        Name:     PARTI_STATS.find_outliers
        Inputs:   None.
        Outputs:  None.
        Features: Convenience function for removing outliers and updating
                  model initial guess values
        Depends:  - remove_linear_outliers
                  - remove_hyperbolic_outliers
                  - update_guess
        """
        self.remove_linear_outliers()
        self.remove_hyperbolic_outliers()
        self.update_guess()

    def fit_hm_obs(self):
        """
        Name:     PARTI_STATS.fit_hm_obs
        Inputs:   None.
        Outputs:  None.
        Features: Curve fits observations to hyperbolic model
        """
        # Define initial parameter values:
        hm_guess = (self.foo_est_obs_h, self.alpha_est_obs_h, self.r_est_obs_h)
        try:
            self.logger.debug("fitting hyperbolic model to observations")
            opt, cov = curve_fit(hyp_model, self.ppfd, self.nee, p0=hm_guess)
        except:
            self.logger.exception("curve fitting failed!")
        else:
            self.logger.debug("optimized hyperbolic values: %s", opt)
            (f_opt, a_opt, r_opt) = opt
            self.foo_opt_obs_h = f_opt
            self.alpha_opt_obs_h = a_opt
            self.r_opt_obs_h = r_opt

            # Calculate the goodness of fit:
            fit = hyp_model(self.ppfd, f_opt, a_opt, r_opt)
            _, rmse, rsqr = goodness_of_fit(fit, self.nee, 3)
            self.rmse_obs_h = rmse
            self.r2_obs_h = rsqr

            # Extract variance values from matrix:
            try:
                self.logger.debug("extracting variance from covariance matrix")
                var = numpy.diag(cov)
            except ValueError:
                self.logger.exception("extraction failed!")
            else:
                if numpy.isfinite(var).all() and not (var < 0).any():
                    (f_err, a_err, r_err) = numpy.sqrt(var)
                    self.foo_err_obs_h = f_err
                    self.alpha_err_obs_h = a_err
                    self.r_err_obs_h = r_err

                    # Calculate the t-value:
                    f_t = f_opt/f_err
                    a_t = a_opt/a_err
                    r_t = r_opt/r_err
                    self.foo_t_obs_h = f_t
                    self.alpha_t_obs_h = a_t
                    self.r_t_obs_h = r_t

                    # Calculate the p-value:
                    dof = self.num_obs - 3.0  # degrees of freedom
                    self.foo_p_obs_h = scipy.stats.t.pdf(-abs(f_t), dof)
                    self.alpha_p_obs_h = scipy.stats.t.pdf(-abs(a_t), dof)
                    self.r_p_obs_h = scipy.stats.t.pdf(-abs(r_t), dof)

    def fit_hm_ro(self):
        """
        Name:     PARTI_STATS.fit_hm_ro
        Inputs:   None.
        Outputs:  None.
        Features: Curve fits observations without outliers to hyperbolic model
        """
        # Remove outliers from observations:
        ppfd_ro = numpy.delete(self.ppfd, self.hyp_outliers)
        nee_ro = numpy.delete(self.nee, self.hyp_outliers)

        # Define initial parameter values:
        hm_guess = (self.foo_est_ro_h, self.alpha_est_ro_h, self.r_est_ro_h)
        try:
            self.logger.debug("fitting hyperbolic model without outliers")
            opt, cov = curve_fit(hyp_model, ppfd_ro, nee_ro, p0=hm_guess)
        except:
            self.logger.exception("curve fitting failed!")
        else:
            self.logger.debug("optimized hyperbolic values: %s", opt)
            (f_opt, a_opt, r_opt) = opt
            self.foo_opt_ro_h = f_opt
            self.alpha_opt_ro_h = a_opt
            self.r_opt_ro_h = r_opt

            # Calculate the goodness of fit:
            fit = hyp_model(ppfd_ro, f_opt, a_opt, r_opt)
            _, rmse, rsqr = goodness_of_fit(fit, nee_ro, 3)
            self.rmse_ro_h = rmse
            self.r2_ro_h = rsqr

            # Extract variance values from matrix:
            try:
                self.logger.debug("extracting variance from covariance matrix")
                var = numpy.diag(cov)
            except ValueError:
                self.logger.exception("extraction failed!")
            else:
                if numpy.isfinite(var).all() and not (var < 0).any():
                    (f_err, a_err, r_err) = numpy.sqrt(var)
                    self.foo_err_ro_h = f_err
                    self.alpha_err_ro_h = a_err
                    self.r_err_ro_h = r_err

                    # Calculate the t-value:
                    f_t = f_opt/f_err
                    a_t = a_opt/a_err
                    r_t = r_opt/r_err
                    self.foo_t_ro_h = f_t
                    self.alpha_t_ro_h = a_t
                    self.r_t_ro_h = r_t

                    # Calculate the p-value:
                    dof = self.num_obs - self.nh_outliers - 3.0
                    self.foo_p_ro_h = scipy.stats.t.pdf(-abs(f_t), dof)
                    self.alpha_p_ro_h = scipy.stats.t.pdf(-abs(a_t), dof)
                    self.r_p_ro_h = scipy.stats.t.pdf(-abs(r_t), dof)

    def fit_lm_obs(self):
        """
        Name:     PARTI_STATS.fit_lm_obs
        Inputs:   None.
        Outputs:  None.
        Features: Curve fits observations to linear model
        """
        # Define the initial guess values:
        lm_guess = (self.alpha_est_obs_l, self.r_est_obs_l)
        try:
            self.logger.debug("fitting linear model to observations")
            opt, cov = curve_fit(lin_model, self.ppfd, self.nee, p0=lm_guess)
        except:
            self.logger.exception("curve fitting failed!")
        else:
            self.logger.debug("optimized linear values %s", opt)
            a_opt, r_opt = opt
            self.alpha_opt_obs_l = a_opt
            self.r_opt_obs_l = r_opt

            # Calculate the goodness of fit:
            fit = lin_model(self.ppfd, a_opt, r_opt)
            _, rmse, rsqr = goodness_of_fit(fit, self.nee, 2)
            self.rmse_obs_l = rmse
            self.r2_obs_l = rsqr

            # Extract variance values from matrix:
            try:
                self.logger.debug("extracting variance from covariance matrix")
                var = numpy.diag(cov)
            except ValueError:
                self.logger.exception("extraction failed!")
            else:
                if numpy.isfinite(var).all() and not (var < 0).any():
                    a_err, r_err = numpy.sqrt(var)
                    self.alpha_err_obs_l = a_err
                    self.r_err_obs_l = r_err

                    # Calculate the t-value:
                    a_t = a_opt/a_err
                    r_t = r_opt/r_err
                    self.alpha_t_obs_l = a_t
                    self.r_t_obs_l = r_t

                    # Calculate the p-value:
                    dof = self.num_obs - 2  # degrees of freedom
                    self.alpha_p_obs_l = scipy.stats.t.pdf(-abs(a_t), dof)
                    self.r_p_obs_l = scipy.stats.t.pdf(-abs(r_t), dof)

    def fit_lm_ro(self):
        """
        Name:     PARTI_STATS.fit_lm_ro
        Inputs:   None.
        Outputs:  None.
        Features: Curve fits observations without outliers to linear model
        """
        # Remove outliers from observations:
        ppfd_ro = numpy.delete(self.ppfd, self.lin_outliers)
        nee_ro = numpy.delete(self.nee, self.lin_outliers)

        # Define the initial guess values:
        lm_guess = (self.alpha_est_ro_l, self.r_est_ro_l)
        try:
            self.logger.debug("fitting linear model without outliers")
            opt, cov = curve_fit(lin_model, ppfd_ro, nee_ro, p0=lm_guess)
        except:
            self.logger.exception("curve fitting failed!")
        else:
            self.logger.debug("optimized linear values %s", opt)
            a_opt, r_opt = opt
            self.alpha_opt_ro_l = a_opt
            self.r_opt_ro_l = r_opt

            # Calculate the goodness of fit:
            fit = lin_model(ppfd_ro, a_opt, r_opt)
            _, rmse, rsqr = goodness_of_fit(fit, nee_ro, 2)
            self.rmse_ro_l = rmse
            self.r2_ro_l = rsqr

            # Extract variance values from matrix:
            try:
                self.logger.debug("extracting variance from covariance matrix")
                var = numpy.diag(cov)
            except ValueError:
                self.logger.exception("extraction failed!")
            else:
                if numpy.isfinite(var).all() and not (var < 0).any():
                    a_err, r_err = numpy.sqrt(var)
                    self.alpha_err_ro_l = a_err
                    self.r_err_ro_l = r_err

                    # Calculate the t-value:
                    a_t = a_opt/a_err
                    r_t = r_opt/r_err
                    self.alpha_t_ro_l = a_t
                    self.r_t_ro_l = r_t

                    # Calculate the p-value:
                    dof = self.num_obs - self.nl_outliers - 2
                    self.alpha_p_ro_l = scipy.stats.t.pdf(-abs(a_t), dof)
                    self.r_p_ro_l = scipy.stats.t.pdf(-abs(r_t), dof)

    def get_initial_guess(self):
        """
        Name:     PARTI_STATS.get_initial_guess
        Inputs:   None.
        Outputs:  None.
        Features: Calculates the initial guess values for the linear and
                  hyperbolic models based on observation data statistics
        Depends:  - get_lm_estimates
                  - get_hm_estimates
        """
        lm_alpha, lm_r = get_lm_estimates(
            self.max_nee_obs,
            self.min_nee_obs,
            self.ave_nee_obs,
            self.std_nee_obs,
            self.max_ppfd_obs,
            self.min_ppfd_obs,
            self.ave_ppfd_obs,
            self.std_ppfd_obs)
        self.alpha_est_obs_l = lm_alpha
        self.r_est_obs_l = lm_r

        hm_foo, hm_alpha, hm_r = get_hm_estimates(
            self.max_nee_obs,
            self.min_nee_obs,
            self.std_nee_obs,
            self.max_ppfd_obs,
            self.min_ppfd_obs)
        self.foo_est_obs_h = hm_foo
        self.alpha_est_obs_h = hm_alpha
        self.r_est_obs_h = hm_r

    def get_summary_idx(self, name):
        """
        Name:     PARTI_STATS.get_summary_idx
        Inputs:   str, field name (name)
        Outputs:  int, dictionary key (k)
        Features: Returns the dictionary key corresponding to the given name
                  field
        """
        found = False
        for k, v in self.summary.items():
            if v['name'] == name:
                found = True
                return k
        if not found:
            raise ValueError(
                "Name '%s' not found in summary dictionary!" % (name))

    def print_summary(self):
        """
        Features: Convenience function for printing the full summary dict
        """
        for k in sorted(list(self.summary.keys())):
            name = self.summary[k]['name']
            val = self.summary[k]['val']
            print("%s: %s" % (name, val))

    def remove_hyperbolic_outliers(self):
        """
        Name:     PARTI_STATS.remove_hyperbolic_outliers
        Inputs:   None.
        Outputs:  None.
        Features: Finds outlier indexes and updates statistics
        Depends:  - hyp_model
                  - get_outliers
                  - calc_statistics
                  - pearsons_r
        """
        # Calculate model fitted values:
        fit = hyp_model(self.ppfd, self.foo_opt_obs_h, self.alpha_opt_obs_h,
                        self.r_opt_obs_h)

        # Get outlier indexes:
        self.hyp_outliers = get_outliers(fit, self.nee, 3)
        self.nh_outliers = len(self.hyp_outliers[0])

        # Remove outliers from observations:
        ppfd_ro = numpy.delete(self.ppfd, self.hyp_outliers)
        nee_ro = numpy.delete(self.nee, self.hyp_outliers)

        # Re-calculate statistics:
        max_p, min_p, ave_p, std_p, skw_p, krt_p = calc_statistics(ppfd_ro)
        self.max_ppfd_ro_h = max_p
        self.min_ppfd_ro_h = min_p
        self.ave_ppfd_ro_h = ave_p
        self.std_ppfd_ro_h = std_p
        self.skw_ppfd_ro_h = skw_p
        self.krt_ppfd_ro_h = krt_p

        max_n, min_n, ave_n, std_n, skw_n, krt_n = calc_statistics(nee_ro)
        self.max_nee_ro_h = max_n
        self.min_nee_ro_h = min_n
        self.ave_nee_ro_h = ave_n
        self.std_nee_ro_h = std_n
        self.skw_nee_ro_h = skw_n
        self.krt_nee_ro_h = krt_n

        self.logger.debug(("calculating correlation coefficient for "
                           "observations without hyperbolic outliers"))
        self.pearson_r_ro_h = pearsons_r(ppfd_ro, nee_ro)

    def remove_linear_outliers(self):
        """
        Name:     PARTI_STATS.remove_linear_outliers
        Inputs:   None.
        Outputs:  None.
        Features: Finds outlier indexes and updates statistics
        Depends:  - lin_model
                  - get_outliers
                  - calc_statistics
                  - pearsons_r
        """
        # Calculate model fitted values:
        fit = lin_model(self.ppfd, self.alpha_opt_obs_l, self.r_opt_obs_l)

        # Get outlier indexes:
        self.lin_outliers = get_outliers(fit, self.nee, 2)
        self.nl_outliers = len(self.lin_outliers[0])

        # Remove outliers from observations:
        ppfd_ro = numpy.delete(self.ppfd, self.lin_outliers)
        nee_ro = numpy.delete(self.nee, self.lin_outliers)

        # Re-calculate statistics:
        max_p, min_p, ave_p, std_p, skw_p, krt_p = calc_statistics(ppfd_ro)
        self.max_ppfd_ro_l = max_p
        self.min_ppfd_ro_l = min_p
        self.ave_ppfd_ro_l = ave_p
        self.std_ppfd_ro_l = std_p
        self.skw_ppfd_ro_l = skw_p
        self.krt_ppfd_ro_l = krt_p

        max_n, min_n, ave_n, std_n, skw_n, krt_n = calc_statistics(nee_ro)
        self.max_nee_ro_l = max_n
        self.min_nee_ro_l = min_n
        self.ave_nee_ro_l = ave_n
        self.std_nee_ro_l = std_n
        self.skw_nee_ro_l = skw_n
        self.krt_nee_ro_l = krt_n

        self.logger.debug(("calculating correlation coefficient for "
                           "observations without linear outliers"))
        self.pearson_r_ro_l = pearsons_r(ppfd_ro, nee_ro)

    def select_best_model(self):
        """
        Name:     PARTI_STATS.select_best_model
        Input:    None.
        Output:   None.
        Features: Saves the model (i.e., none, hyperbolic with
                  observations, hyperbolic with outliers removed, linear with
                  observations, linear with outliers removed) that ``best''
                  represents the data based on thesholds of fitness, validity
                  of parameter values, and parameter significance
        """
        self.logger.debug("Selecting the best model")

        # Define threshold values:
        r2_thresh = 0.5       # minimum R-squared fitness
        p_thresh = 0.05       # minimum parameter significance (95%)

        # Initialize selection dictionary:
        select_dict = {'obs_h': 1, 'ro_h': 2, 'obs_l': 3, 'ro_l': 4}

        # Initialize model failure dictionary:
        model_dict = {'obs_h': 0, 'ro_h': 0, 'obs_l': 0, 'ro_l': 0}

        # Initialize R-squared dictionary:
        r2_dict = {'obs_h': self.r2_obs_h, 'ro_h': self.r2_ro_h,
                   'obs_l': self.r2_obs_l, 'ro_l': self.r2_ro_l}

        # Initialize P-value dictionary:
        p_dict = {
            'obs_h': numpy.array(
                [self.foo_p_obs_h, self.alpha_p_obs_h, self.r_p_obs_h]),
            'ro_h': numpy.array(
                [self.foo_p_ro_h, self.alpha_p_ro_h, self.r_p_ro_h]),
            'obs_l': numpy.array([self.alpha_p_obs_l, self.r_p_obs_l]),
            'ro_l': numpy.array([self.alpha_p_ro_l, self.r_p_ro_l])}

        # Initialize parameter value dictionary:
        value_dict = {
            'obs_h': {'foo': self.foo_opt_obs_h,
                      'alpha': self.alpha_opt_obs_h,
                      'r': self.r_opt_obs_h},
            'ro_h': {'foo': self.foo_opt_ro_h,
                     'alpha': self.alpha_opt_ro_h,
                     'r': self.r_opt_ro_h},
            'obs_l': {'alpha': self.alpha_opt_obs_l,
                      'r': self.r_opt_obs_l},
            'ro_l': {'alpha': self.alpha_opt_ro_l,
                     'r': self.r_opt_ro_l}}

        # Initialize model selection:
        best_model = 0

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIRST CRITERIA--- R-SQUARED THRESHOLD
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (numpy.array(list(r2_dict.values())) < r2_thresh).any():
            # Eliminate those with poor R2:
            m_failed = [k for k, v in r2_dict.items() if v < r2_thresh]
            for m in m_failed:
                self.logger.debug("%s failed r-squared test", m)
                model_dict[m] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # SECOND CRITERIA--- PARAMETER RANGE OF VALIDITY
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # -- rectangular hypberbola:
        if (value_dict['obs_h']['foo'] > 100 or
                value_dict['obs_h']['foo'] < 0 or
                value_dict['obs_h']['alpha'] > 1 or
                value_dict['obs_h']['alpha'] < 0):
            self.logger.debug("obs_h failed validity test")
            model_dict['obs_h'] = 1
        if (value_dict['ro_h']['foo'] > 100 or
                value_dict['ro_h']['foo'] < 0 or
                value_dict['ro_h']['alpha'] > 1 or
                value_dict['ro_h']['alpha'] < 0):
            self.logger.debug("ro_h failed validity test")
            model_dict['ro_h'] = 1

        # -- linear:
        if (value_dict['obs_l']['alpha'] > 1 or
                value_dict['obs_l']['alpha'] < 0):
            self.logger.debug("obs_l failed validity test")
            model_dict['obs_l'] = 1
        if (value_dict['ro_l']['alpha'] > 1 or
                value_dict['ro_l']['alpha'] < 0):
            self.logger.debug("ro_l failed validity test")
            model_dict['ro_l'] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # THIRD CRITERIA--- PARAMETER SIGNIFICANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (p_dict['obs_h'] > p_thresh).any():
            self.logger.debug("obs_h failed significance test")
            model_dict['obs_h'] = 1
        if (p_dict['ro_h'] > p_thresh).any():
            self.logger.debug("ro_h failed significance test")
            model_dict['ro_h'] = 1
        if (p_dict['obs_l'] > p_thresh).any():
            self.logger.debug("obs_l failed significance test")
            model_dict['obs_l'] = 1
        if (p_dict['ro_l'] > p_thresh).any():
            self.logger.debug("ro_l failed significance test")
            model_dict['ro_l'] = 1

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ANALYSIS:
        # Let's see which models passed the selection criteria and how many
        # of them there are
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        mods = [k for k, v in model_dict.items() if v == 0]
        howmany = len(mods)
        self.logger.debug("Good models found: %d", howmany)

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
        else:
            best_model = 0

        # Save model selection:
        self.logger.debug("Best model is %s", best_model)
        self.model_select = best_model

    def update_guess(self):
        """
        Name:     PARTI_STATS.update_guess
        Inputs:   None.
        Outputs:  None.
        Features: Calculates the initial guess values for the linear and
                  hyperbolic models based on observations without outliers
                  statistics
        Depends:  - get_lm_estimates
                  - get_hm_estimates
        """
        lm_alpha, lm_r = get_lm_estimates(
            self.max_nee_ro_l,
            self.min_nee_ro_l,
            self.ave_nee_ro_l,
            self.std_nee_ro_l,
            self.max_ppfd_ro_l,
            self.min_ppfd_ro_l,
            self.ave_ppfd_ro_l,
            self.std_ppfd_ro_l)
        self.alpha_est_ro_l = lm_alpha
        self.r_est_ro_l = lm_r

        hm_foo, hm_alpha, hm_r = get_hm_estimates(
            self.max_nee_ro_h,
            self.min_nee_ro_h,
            self.std_nee_ro_h,
            self.max_ppfd_ro_h,
            self.min_ppfd_ro_h)
        self.foo_est_ro_h = hm_foo
        self.alpha_est_ro_h = hm_alpha
        self.r_est_ro_h = hm_r

    def set_observations(self, nee, ppfd):
        """
        Name:     PARTI_STATS.set_observations
        Inputs:   - numpy.ndarray, NEE observations (nee)
                  - numpy.ndarray, PPFD observations (ppfd)
        Outputs:  None.
        Features: Saves observation data and statistics
        Depends:  calc_statistics
        """
        # Set the number of observations:
        if len(nee) == len(ppfd):
            self.nee = nee
            self.ppfd = ppfd
            self.num_obs = len(nee)
            self.logger.debug("Saving %d observations", len(nee))
        else:
            raise ValueError("NEE and PPFD array lengths must match!")

        # Calculate the observation stats:
        self.logger.debug("Calculating observation statistics")
        (max_p, min_p, ave_p, std_p, skw_p, krt_p) = calc_statistics(ppfd)
        self.max_ppfd_obs = max_p
        self.min_ppfd_obs = min_p
        self.ave_ppfd_obs = ave_p
        self.std_ppfd_obs = std_p
        self.skw_ppfd_obs = skw_p
        self.krt_ppfd_obs = krt_p

        (max_n, min_n, ave_n, std_n, skw_n, krt_n) = calc_statistics(nee)
        self.max_nee_obs = max_n
        self.min_nee_obs = min_n
        self.ave_nee_obs = ave_n
        self.std_nee_obs = std_n
        self.skw_nee_obs = skw_n
        self.krt_nee_obs = krt_n

        self.logger.debug("Calculating observation correlation coefficient")
        self.pearson_r_obs = pearsons_r(nee, ppfd)

    def write_fit_params(self, output_dir):
        """
        Name:     PARTI_STATS.write_fit_params
        Inputs:   str, output directory path (output_dir)
        Outputs:  None.
        Features: Writes the partitioning parameters and data to file
        Note:     Output is read by plot_partitioning.R (in tools/plotting/)
                  Headers are organized as: alpha, R, foo, RMSE, R2 in order
                  to comply with the optim_params R function
        """
        # ~~~~~~~~~~ OBSERVATIONS ~~~~~~~~~~
        obs_file = "%s_%s_obs.txt" % (self.name, self.date)
        obs_path = os.path.join(output_dir, obs_file)

        head0 = "MH_guess,%f,%f,%f\n" % (
            self.alpha_est_obs_h, self.r_est_obs_h, self.foo_est_obs_h)

        head1 = "ML_guess,%f,%f\n" % (
            self.alpha_est_obs_l, self.r_est_obs_l)

        head2 = "MH_opt,%f,%f,%f,%f,%f\n" % (
            self.alpha_opt_obs_h, self.r_opt_obs_h, self.foo_opt_obs_h,
            self.rmse_obs_h, self.r2_obs_h)

        head3 = "ML_opt,%f,%f,,%f,%f\n" % (
            self.alpha_opt_obs_l, self.r_opt_obs_l, self.rmse_obs_l,
            self.r2_obs_l)

        try:
            self.logger.debug("opening %s for writing", obs_path)
            OUTFILE = open(obs_path, 'w')

            self.logger.debug("calculating model fitted values")
            fit_h = hyp_model(
                self.ppfd, self.foo_opt_obs_h, self.alpha_opt_obs_h,
                self.r_opt_obs_h)
            fit_l = lin_model(
                self.ppfd, self.alpha_opt_obs_l, self.r_opt_obs_l)
        except:
            self.logger.exception("failed to write %s", obs_path)
        else:
            OUTFILE.write(head0)
            OUTFILE.write(head1)
            OUTFILE.write(head2)
            OUTFILE.write(head3)
            OUTFILE.write("ppfd_obs,nee_obs,nee_mh,nee_ml\n")

            for po, no, nh, nl in map(None, self.ppfd, self.nee, fit_h, fit_l):
                outline = "%s,%s,%s,%s\n" % (po, no, nh, nl)
                OUTFILE.write(outline)

            OUTFILE.close()

        # ~~~~~~~~~~ OUTLIERS REMOVED ~~~~~~~~~~
        ro_file = "%s_%s_ro.txt" % (self.name, self.date)
        ro_path = os.path.join(output_dir, ro_file)

        head4 = "MH_guess,%f,%f,%f\n" % (
            self.alpha_est_obs_h, self.r_est_obs_h, self.foo_est_obs_h)
        head5 = "ML_guess,%f,%f\n" % (
            self.alpha_est_obs_l, self.r_est_obs_l)
        head6 = "MH_opt,%f,%f,%f,%f,%f\n" % (
            self.alpha_opt_ro_h, self.r_opt_ro_h, self.foo_opt_ro_h,
            self.rmse_ro_h, self.r2_ro_h)
        head7 = "ML_opt,%f,%f,,%f,%f\n" % (
            self.alpha_opt_ro_l, self.r_opt_ro_l,
            self.rmse_ro_l, self.r2_ro_l)

        try:
            OUTFILE = open(ro_path, 'w')

            # Remove outliers:
            ppfd_rol = numpy.delete(self.ppfd, self.lin_outliers)
            nee_rol = numpy.delete(self.nee, self.lin_outliers)
            ppfd_roh = numpy.delete(self.ppfd, self.hyp_outliers)
            nee_roh = numpy.delete(self.nee, self.hyp_outliers)

            # Calculate model fitted parameters:
            fit_h = hyp_model(ppfd_roh, self.foo_opt_ro_h, self.alpha_opt_ro_h,
                self.r_opt_ro_h)
            fit_l = lin_model(ppfd_rol, self.alpha_opt_ro_l, self.r_opt_ro_l)
        except:
            self.logger.exception("failed to open file %s", ro_path)
        else:
            OUTFILE.write(head4)
            OUTFILE.write(head5)
            OUTFILE.write(head6)
            OUTFILE.write(head7)
            OUTFILE.write("ppfd_obs_h,nee_obs_h,nee_mod_h,"
                          "ppfd_obs_l,nee_obs_l,nee_mod_l\n")
            for (poh, noh, nhh, pol, nol, nll) in map(None, ppfd_roh, nee_roh,
                                                      fit_h, ppfd_rol, nee_rol,
                                                      fit_l):
                outline = "%s,%s,%s,%s,%s,%s\n" % (
                    poh, noh, nhh, pol, nol, nll)
                OUTFILE.write(outline)
            OUTFILE.close()
