#!/usr/bin/python
#
# model.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-13
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
#
# ------------
# description:
# ------------
# This script presents the Python application code for the Global ecosystem
# in Space and Time (GePiSaT) model of the terrestrial biosphere.
#
# GePiSaT takes a simplistic approach to modelling terrestrial gross primary
# production (GPP) by making the best use of in situ observations while
# defensibly representing the principal ecophysiological processes that govern
# GPP, including:
# (1) the eddy-covariance method of partitioning net CO2 and solar radiative
#     fluxes into monthly quantities of GPP and respiration; and
# (2) the optimality principle of vegetation minimizing the summed costs
#     associated with maintaining carbon fixation and water transport
#     capabilities.
#
# This script is used to connect to the gepisat database, query and process
# data, and perform model analyses.
#
# Currently this includes the following:
#
# **STAGE 1**
# -- Acquires monthly PPFD and NEE observation pairs from postgreSQL database
# -- Performs GPP and Re partitioning of NEE and PPFD based on models defined
#    by Ruimy et al. 1995
# -- Removes outliers in PPFD:NEE observations based on model fitting using
#    Peirce's criterion (Peirce 1852; Gould 1855)
# -- Outputs original and outlier-free monthly datasets w/ fitting params
#    > [STATION]_[YYYY-MM-DD].txt, original
#    > [STATION]_[YYYY-MM-DD]_ro.txt, outlier-free
# -- Outputs partitioning statistics to file
#    > summary_statistics.txt
# -- Selects the best model parameters for modeling
#    > selection criteria are based on: model fitness, parameter significance
#
# **STAGE 2**
# -- Returns WFDEI shortwave radiation grid measurement based on flux tower
#    location and day of the year
# -- Calculates extraterrestrial solar radiation (in units of PPFD)
# -- Gap-fills daily PPFD observations with scaled extraterrestrial solar
# -- Outputs observations and gap-filled time series
#    > [STATION]-GF_[YYYY-MM-DD].txt
# -- Calculates half-hourly GPP based on partitioning
# -- Integrates monthly gap-filled PPFD and GPP
# -- Returns fAPAR based on flux tower location and month
# -- Curve-fits light-use efficiency based on PEM: GPP = e fPAR PPFD
# -- Outputs station-specific LUE parameters
#    > [STATION]_LUE.txt
# -- @TODO: Outputs summary of LUE for each month for each station
#    > LUE_All_Stations.txt
#
# **STAGE 3**
# -- @TODO: Global prediction
#
# -----
# todo:
# -----
# - FLUX_PARTI class
#   + Why not use the optimized obs parameters as the guess for ro?
#       * if optim fails---you'll have crumby initial guesses
#   + add variance of residuals to summary statistics (4 fields)
#       * goodness_of_fit() function
# - partition function:
#   + check data validity after outliers are removed before reprocessing
# - model_select()
#   + Consider implementing either Shapiro-Wilks or Anderson test of
#     normality on model residuals
#   --> scipy.stats.shapiro(my_resids)
#       * returns test statistic and p-value
#       * Note: p-value indicates significantly different from normal
#   --> scipy.stats.anderson(my_resids, 'norm')
#       * returns test statistic, array of critical values, and sig values
#       * if test statistic is larger than the critical value at the
#         significance value you are interested in, then the null can be
#         rejected, i.e., not normal
#   + Or use Willmott's revised index of agreement
# - Consider checking for system closure (each station, each month?)
#   + short equation: Rn + G + LE + H = 0
#
# NOTE: IT-SR2 and ZM-Mon not processed.
#
###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import logging
import os

from gepisat.data import DATA
from gepisat.flux_parti import FLUX_PARTI
from gepisat.utilities import add_one_month
from gepisat.utilities import simpson


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == "__main__":
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.FileHandler("run.log", mode='w')
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # Initialize data class for handling all data IO operations, which
    # includes connecting to the GePiSaT database:
    my_data = DATA()

    # Define output directory:
    output_dir = os.path.join(
        os.path.expanduser("~"), "Desktop", "temp", "out_flux")
    part_out_dir = os.path.join(output_dir, "partition")
    my_data.set_output_directory(output_dir)

    # Get list of all flux station names:
    my_data.find_stations()         # use this function to search the DB
    # my_data.stations = ['CZ-wet', ]   # or manually define stations here

    # Initialize summary statistics file:
    my_data.create_summary_file("summary_statistics.txt")

    # Iterate through stations:
    for station in my_data.stations:
        try:
            # Initialize station data:
            my_data.set_working_station(station)
        except:
            root_logger.error(
                "Failed to set working station; skipping %s", station)
        else:
            print("Processing %s" % (station))
            root_logger.info("Processing %s", station)
            # Process each month in time:
            sd_of_interest = datetime.date(2003, 1, 1)
            ed_of_interest = datetime.date(2015, 1, 1)
            if my_data.start_date > sd_of_interest:
                sd = my_data.start_date
            else:
                sd = sd_of_interest
            if my_data.end_date < ed_of_interest:
                ed = my_data.end_date
            else:
                ed = ed_of_interest
            while sd < ed:
                try:
                    # Get PPFD and NEE array pairs [umol m-2 s-1]:
                    mo_nee, mo_ppfd = my_data.find_monthly_nee_ppfd(sd)
                except:
                    root_logger.error(
                        "Failed to find monthly NEE/PPFD; skipping %s %s" % (
                            station, sd))
                else:
                    # Initialize a FLUX_PARTI class:
                    my_parter = FLUX_PARTI(mo_ppfd, mo_nee, station, sd)

                    # Perform GPP partitioning:
                    my_parter.partition(outliers=True)

                    # OPTIONAL: write the partition results for plotting
                    if True:
                        my_parter.write_partition(part_out_dir)

                    # Perform half-hourly PPFD gapfilling (umol m-2 s-1):
                    gf_time, gf_ppfd = my_data.gapfill_monthly_ppfd(
                        sd, to_save=False)

                    # Calculate half-hourly GPP (umol m-2 s-1)
                    gf_gpp, gf_gpp_err = my_parter.calc_gpp(gf_ppfd)

                    # OPTIONAL: calculate daily GPP (mol m-2):
                    if True:
                        gpp_daily = my_data.sub_to_daily_gpp(
                            gf_time, gf_gpp, gf_gpp_err, to_save=True)

                    # Continue processing if partitioning was successful:
                    if my_parter.best_model > 0:
                        # Integrate PPFD & GPP [umol m-2]; dt = 30 min
                        ppfd_mo = simpson(gf_ppfd.clip(min=0), 1800)
                        gpp_mo = simpson(gf_gpp.clip(min=0), 1800)
                        gpp_mo_err = simpson(gf_gpp_err, 1800)

                        # Convert units from [umol m-2] to [mol m-2]:
                        ppfd_mo *= (1e-6)
                        gpp_mo *= (1e-6)
                        gpp_mo_err *= (1e-6)

                        # STAGE 2: Update LUE parameters:
                        my_data.update_lue(sd, gpp_mo, gpp_mo_err, ppfd_mo)

                    # Save class summary statistics:
                    SFILE = open(my_data.summary_file, 'a')
                    SFILE.write(my_parter.summary_stats)
                    SFILE.close()
                finally:
                    # Increment date:
                    sd = add_one_month(sd)

            # STAGE 2: Write monthly LUE parameters to file:
            my_data.write_monthly_results(station)
