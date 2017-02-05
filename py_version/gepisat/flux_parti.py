#!/usr/bin/python
#
# flux_parti.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-02-05
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
import logging
import os

import numpy

from .data import DATA
from .stats import PARTI_STATS


###############################################################################
# CLASSES
###############################################################################
class FLUX_PARTI(object):
    """
    Name:     FLUX_PARTI
    Features: This class performs flux partitioning of monthly NEE & PPFD
              observations
    History   Version 3.0.0-dev
              - class separated from model [16.01.17]
              - Python 2/3 supported print statements [16.01.17]
              - moved utility functions to utilities script [16.04.01]
                * calc statistics
                * goodness of fit
                * pearsons r
                * peirce dev
              - reduced function calls by adding modes [16.04.01]
              - added object inheritance for Python2/3 support [16.05.20]
              - moved the vast majority of function and attributes to
                utilities script and PARTI_STATS class [16.05.20]
              - changed partition variable to outliers [16.07.06]
              - updated calc gpp function [16.07.06]
              - created best model class property [16.07.06]
              - moved to gepisat package [16.07.22]
              - created a write partition function [17.01.22]
              - added array length check in partition [17.01.25]
              - return NaNs when no observation pairs present [17.02.05]
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, ppfd, nee, tower, mo):
        """
        Name:     FLUX_PARTI.__init__
        Input:    - list, monthly PPFD observations (ppfd)
                  - list, monthly NEE observations (nee)
                  - string, flux tower name (tower)
                  - datetime.date, current month (mo)
        Features: Initializes the observation data statistics and initial
                  model fitting parameters
        """
        # Create a class logger
        self.logger = logging.getLogger(__name__)
        self.logger.debug(
            "Initializing FLUX_PARTI class for %s @ %s" % (tower, mo))

        # Initialize stats class:
        self.stats = PARTI_STATS()
        self.stats.name = tower
        self.stats.date = mo

        # Calculate observation statistics:
        self.stats.set_observations(nee, ppfd)

        # Update initial model parameters:
        self.stats.get_initial_guess()

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Property Definitions
    # ////////////////////////////////////////////////////////////////////////
    @property
    def best_model(self):
        """Best model based on fitness statistics"""
        return self.stats.model_select

    @best_model.setter
    def best_model(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    @property
    def summary_stats(self):
        """The full summary statistics"""
        return self.stats.summary_string

    @summary_stats.setter
    def summary_stats(self, val):
        raise NotImplementedError("You can not set this parameter this way.")

    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def calc_gpp(self, ppfd):
        """
        Name:     FLUX_PARTI.calc_gpp
        Input:    numpy.ndarray, monthly PPFD (ppfd)
        Output:   tuple, modeled GPP and associated error
                  - numpy.ndarray, modeled GPP (gpp)
                  - numpy.ndarray, associated error (gpp_err)
        Features: Returns GPP (umol m-2 s-1) and associated error based on the
                  best fitted model
        """
        self.logger.debug("calculating GPP for model %d", self.best_model)

        # Get model selection:
        if (self.best_model == 1 or self.best_model == 2):
            model = "H"
        elif (self.best_model == 3 or self.best_model == 4):
            model = "L"

        # Get outlier selection:
        if (self.best_model == 1 or self.best_model == 3):
            outlier = 0
        elif (self.best_model == 2 or self.best_model == 4):
            outlier = 1

        if self.stats.num_obs == 0:
            # No NEE:PPFD; therefore, no model fits, return NaNs
            gpp = numpy.ones(len(ppfd)) * numpy.nan
            gpp_err = numpy.zeros(len(ppfd))
        elif self.best_model == 0:
            # Default to hyperbolic model with outliers removed estimates
            foo = self.stats.foo_est_obs_h
            foo_err = 0.0
            alpha = self.stats.alpha_est_obs_h
            alpha_err = 0.0

            if foo == 0 or alpha == 0:
                self.logger.debug(
                    "No best model, defaulting to linear observed model")

                # Use linear estimates instead:
                alpha = self.stats.alpha_est_obs_l
                alpha_err = 0.

                # Calculate GPP and GPP err:
                gpp = (alpha*ppfd)
                gpp_err = (ppfd*alpha_err)
            else:
                self.logger.debug(
                    "No best model, defaulting to hyperbolic observed model")

                # Variable substitutes:
                apf = alpha*ppfd + foo
                afp = alpha*foo*ppfd
                afpp = alpha*foo*(ppfd**2)

                # Calculate GPP and its associated error:
                gpp = afp/apf
                gpp_err = numpy.sqrt(
                    (alpha_err**2)*((foo*ppfd*apf - afpp)/(apf**2))**2 +
                    (foo_err**2)*((alpha*ppfd*apf - afp)/(apf**2))**2)

            # Clip out negative values of GPP:
            gpp = gpp.clip(min=0)
        elif model.upper() == "H":
            # Retrieve parameters for hyperbolic model:
            if outlier == 0:
                self.logger.debug("Using optimized hyperbolic model")
                foo = self.stats.foo_opt_obs_h
                foo_err = self.stats.foo_err_obs_h
                alpha = self.stats.alpha_opt_obs_h
                alpha_err = self.stats.alpha_err_obs_h
                # (foo, alpha, r) = self.hm_optimized
                # (foo_err, alpha_err, r_err) = self.hm_optim_err
            elif outlier == 1:
                self.logger.debug(
                    "Using optimized hyperbolic model with outliers removed")
                foo = self.stats.foo_opt_ro_h
                foo_err = self.stats.foo_err_ro_h
                alpha = self.stats.alpha_opt_ro_h
                alpha_err = self.stats.alpha_err_ro_h
                # (foo, alpha, r) = self.hm_optimized_ro
                # (foo_err, alpha_err, r_err) = self.hm_optim_err_ro

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
                self.logger.debug("Using optimized linear model")
                alpha = self.stats.alpha_opt_obs_l
                alpha_err = self.stats.alpha_err_obs_l
                # (alpha, r) = self.lm_optimized
                # (alpha_err, r_err) = self.lm_optim_err
            elif outlier == 1:
                self.logger.debug(
                    "Using optimized linear model with outliers removed")
                alpha = self.stats.alpha_opt_ro_l
                alpha_err = self.stats.alpha_err_ro_l
                # (alpha, r) = self.lm_optimized_ro
                # (alpha_err, r_err) = self.lm_optim_err_ro

            # Calculate GPP:
            gpp = (alpha*ppfd)
            gpp_err = (ppfd*alpha_err)

        # Return GPP
        return (gpp, gpp_err)

    def partition(self, outliers=False):
        """
        Name:     FLUX_PARTI.partition
        Input:    [optional] bool, whether to remove outliers (rm_out)
        Output:   None.
        Features: Performs flux partitioning based on NEE:PPFD observations;
                  processes outliers (optional);
                  saves results to file (optional)
        """
        if self.stats.num_obs > 3:
            self.logger.debug("fitting hyperbolic and linear models")
            self.stats.fit_hm_obs()
            self.stats.fit_lm_obs()

            if outliers:
                self.logger.debug("finding outliers ...")
                self.stats.find_outliers()
                self.logger.debug(
                    "fitting hyperbolic and linear models without outliers")
                self.stats.fit_hm_ro()
                self.stats.fit_lm_ro()

            self.stats.select_best_model()

    def write_partition(self, out_dir):
        """
        Name:     FLUX_PARTI.write_partition
        Inputs:   str, output directory (out_dir)
        Outputs:  None.
        Features: Writes partitioning parameters to file

        @TODO: to omit write outs if partitioning was not performed?
        """
        if os.path.isdir(out_dir):
            self.stats.write_fit_params(out_dir)
        else:
            self.logger.error("Output directory does not exist!")


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("flux_parti.log")
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # DEBUG
    my_data = DATA()
    my_data.output_dir = os.path.join(
        os.path.expanduser("~"), "Desktop", "temp", "out")
    my_stations = my_data.find_stations()
    station = my_data.stations[1]
    my_data.set_working_station(station)
    sd = my_data.start_date
    nee, ppfd = my_data.find_monthly_nee_ppfd(sd)
    my_parter = FLUX_PARTI(ppfd, nee, station, sd)
    my_parter.partition(True)
    my_parter.stats.print_summary()
