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
#    **STAGE 1**
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
#    **STAGE 2**
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
# -- Outputs summary of LUE for each month for each station
#    > LUE_All_Stations.txt
#
#    **STAGE 3**
# -- TBA
#
# ----------
# changelog:
# ----------
# VERSION 0.0.0
#  - file created based on tableQuery.py (v.0.15) [13.07.05]
#  - import module datetime [13.07.08]
# VERSION 0.0.1
#  - added station and date to ppfd and nee query [13.07.09]
#  - filtered ppfd and nee pairs (only when both are observed) [13.07.09]
#  - added add_one_month function [13.07.09]
#  - added recthyp function [13.07.09]
#  - import numpy and scipy modules [13.07.09]
# VERSION 0.0.2
#  - changed requirements on ppfd and nee pairs [13.07.10]
#    --> ppfd has to be positive
#  - added goodness_of_fit function [13.07.10]
# VERSION 0.0.3
#  - added get_dates function [13.07.11]
#  - added update_guess function [13.07.11]
#  - added check on monthly_pair; need at least two observations for
#    processing regressions [13.07.11]
#    --> updated to >2 observations (required for model_h) [13.07.16]
#    --> updated to >3 observations otherwise mse divides by zero [13.07.19]
# VERSION 0.0.4
#  - added station iteration [13.07.12]
# VERSION 0.0.5
#  - updated naming convensions for PEP8 [13.07.16]
#  - updated get_dates function [13.07.16]
#    --> made startind date begin at day 1
#  - added PPFD & NEE list sepration in monthly_ppfd_nee function [13.07.16]
# VERSION 0.1.0
#  - changed the way model_h handles zero division errors [13.07.17]
#    --> add +1e-6 to denominator
#  - added standard error lists for optimization parameters and added
#    calculation for them [13.07.17]
#  - import sys (noticed is was used in connect_sql function) [13.07.17]
#  - changed datasets from lists to numpy.ndarrays [13.07.17]
#  - implemented class MyClass [13.07.17]
#  - import scipy.special [13.07.17]
#  - added peirceDev function [13.07.17]
#  - implemented  outlier removal [13.07.17]
#  - now outputs two files (w & w/o outliers) [13.07.17]
# VERSION 0.2.0
#  - significant formatting updates [13.07.18]
#  - added summary_statistics function and output [13.07.18]
#  - important update to crosstab query in monthly_ppfd_nee function;
#    added category_sql [13.07.18]
#  - ppfd and nee records flipped in monthly_ppfd_nee tuples [13.07.18]
#    --> for some reason ppfd and nee are received in the main code in the
#        opposite order as sent
# VERSION 0.2.1
#  - updated update_guess:
#    --> based on first round of summary statistics [13.07.20]
#    --> based on second round of summary statistics [13.07.22]
#    --> also updated decimal places in the write to file sections [13.07.22]
#  - added two more decimal places (%0.5f) to alpha in output [13.07.20]
#    --> added additional decimal place to rmse values [13.07.22]
#  - changed r2 -> r2.adj in goodness_of_fit function [13.07.22]
# VERSION 0.3.0
#  - changed class name from "MyClass" to "FLUX_PARTI" [13.09.11]
#  - started get_daily_flux and add_one_day functions [13.09.11]
#  - added con.close() commands at the end of sql functions [13.09.11]
# VERSION 0.4.0
#  - created summary_file_init() function [13.09.12]
#  - created partition() function w/ write_out and rm_out booleans [13.09.12]
#  - updated update_guess ZeroDivision exception handling; changed
#    denominator only [13.09.12]
#  - updated connect_sql() to look for credential file [13.09.12]
#  - added os.path to module list [13.09.12]
# VERSION 0.5.0
#  - created get_lon_lat() function [13.09.13]
#  - added grid_centroid() function [13.09.13]
#  - created flux_to_grid() function [13.09.13]
#  - added SOLAR class [13.09.13]
#  - created gapfill_ppfd() function [13.09.13]
#  - created get_data_point() function [13.09.13]
# VERSION 0.5.1
#  - updated solar class [13.09.16]
#    --> includes integrals of daily radiation
#  - updated gapfill_ppfd() [13.09.16]
#    --> calculate julian day from datetime.date object
# VERSION 0.5.2
#  - updated solar class [13.09.17]
#    --> daily PPFD integral now in units of umol m-2
#    --> implemented Simpson's rule instead of Trapezoidal rule
#  - updated gapfill_ppfd() [13.09.17]
#    --> convert grid (WATCH) PPFD to daily units (umol m-2)
#  - updated FLUX_PARTI class [13.09.17]
#    --> added calc_gpp function w/ three user parameters
# VERSION 0.5.3
#  - fixed typo in simpson function [13.09.18]
#  - convert units of umol m-2 to mol m-2 in integrations [13.09.18]
# VERSION 0.6.0
#  - added fapar retrieval [13.10.01]
# VERSION 0.7.0
#  - added lue output (int gpp, int ppfd, fAPAR) [13.10.02]
#    --> moved to LUE class [13.10.03]
# VERSION 0.8.0
#  - added new class LUE [13.10.03]
# VERSION 1.0.0
#  - re-named script to "model.py" [13.10.11]
# VERSION 1.0.1-run01
#  - updated for new bitnami db [13.10.18]
# VERSION 1.0.2-run02
#  - updated goodness_of_fit for minimum value check (=4) [13.10.18]
# VERSION 1.0.3-run03
#  - updated remove_mh|ml_outliers functions [13.10.18]
#    --> checks that n does not exceed N
# VERSION 1.0.4-run04
#  - updated calc_statistics [13.10.18]
#    --> check that my_array is not empty
# VERSION 1.0.5-run05
#  - updated partition function [13.10.18]
#    --> added TypeError exception in curve_fit outlier-free data (model H&L)
# VERSION 1.0.6-run06
#  - updated peirce_dev [13.10.20]
#    --> if N <=1, return x2 = 0 (otherwise divide by zero error)
#    --> nan check on Lamda calc
#    --> return 0 for negative x2 vals
#  - updated get_stations() [13.10.20]
#    --> added where clause for dim=0 (don't select grid stations)
# VERSION 1.0.7-run07
#  - added check for NaN in curve_fit covariance array [13.10.20]
# VERSION 1.0.8-run08
#  - generalized exceptions for curve_fit (all cases) [13.10.20]
# VERSION 1.0.9-run09
#  - updated remove_ml|mh_outliers & peirce_dev functions [13.11.04]
#    --> unnecessary but kept
# VERSION 1.0.10-run11
#  - updated the write out in partition function [13.11.04]
#    --> changed zip to map and padded with "None"
# VERSION 1.1.0-run12
#  - added obs and ro versions of NEE and PPFD statistics [13.11.12]
#  - created save_stats() function for saving obs and ro stats [13.11.12]
#  - created save_estimates() function for obs and ro params [13.11.12]
#  - updated summary statistics fields (from 50 to 79) [13.11.12]
# VERSION 1.2.0-run13
#  - added VPD to LUE class & calc_lue() [13.11.18]
# VERSION 1.3.0-run14
#  - added RH100 (canopy height, m) to LUE class & calc_lue() [13.11.20]
# VERSION 1.3.0-run15
#  - updated canopy height data (based on original TIFF) [13.11.27]
# VERSION 1.4.0
#  - fixed error in calc_gpp in FLUX_PARTI class [14.01.10]
#    --> hm_optimized(_ro) parameters given in wrong order
#  - uncommented peirce_x2 and peirce_d2 in FLUX_PARTI [14.01.10]
#    --> values are self referenced in remove_mh|ml_outliers functions
# VERSION 1.4.0-run16 [14.01.11]
# VERSION 1.5.0
#  - GPP calculation error [14.01.12]
#    --> source: http://chemwiki.ucdavis.edu/Analytical_Chemistry/            \
#        Quantifying_Nature/Significant_Digits/Propagation_of_Error
#    --> added to calc_gpp function in FLUX_PARTI
#  - added gpp_err as second return value [14.01.12]
#    --> gpp_err is then aggregated using simpson
#        ? square-root of the sum of the squared error was too small
#    --> updated LUE add_station_val and write_out to accommodate gpp_err
# VERSION 1.5.0-run17 [14.01.13]
# VERSION 1.6.0
#  - removed canopy height (RH100) from LUE [14.02.04]
#  - added alpha, Tc, CO2, and Patm to LUE [14.02.04]
#  - implemented new LUE [14.02.04]
# VERSION 1.6.1-run18
#  - updated get_stations (added geom check) [14.02.05]
# VERSION 1.6.2-run19
#  - updated calc_k, calc_kc, calc_ko, calc_gstar functions [14.02.19]
#    --> use constant partial pressures for kc, ko and gstar
#  - updated lue_model [14.02.19]
#    --> no more V1-V4 shananagens
# VERSION 1.6.3
#  - added [lm/hm]_resid_var[_ro] to variable list in FLUX_PARTI [14.03.10]
# VERSION 1.7.0-run20
#  - added Pearson's r function [14.03.11]
#  - included Pearson's r for obs of NEE & PPFD to summary stats [14.03.11]
#  - updated calc methods in calc_statistics function [14.03.11]
#    --> take advantage of numpy arrays for .max(), .min(), .mean(), and .std()
#  - updated optim err filters in partition() [14.03.11]
#    --> isfinite(cov).all() and not (cov<0).any()
#  - added model selection criteria [14.03.11]
#  - updated calculate_gpp() function for model selection [14.03.11]
#  - removed decimal restriction in SQBA param in LUE class [14.03.11]
# VERSION 1.7.1-run21
#  - updated Pearson's R formulation [14.03.12]
#  - updated model selection criteria [14.03.12]
#    --> parameter boundaries
#    --> fitness (R2 >= 0.2)
#    --> / outliers versus observations
#    --> \ simplest model (linear versus hyperbolic)
#    --> R2 comparison
# VERSION 1.7.2-run22
#  - removed decimal limitation on alpha in sum_stat [14.03.13]
#  - removed decimal limitation on all params in partition() [14.03.13]
# VERSION 1.7.3
#  - fixed check on numpy.diag for parameter st error calculation [14.03.18]
#  - changed summary stats decimal places for params to 5 [14.03.18]
# VERSION 1.7.4
#  - moved flux_to_grid and get_msvidx outside sd<ed loop [14.03.20]
# VERSION 1.8.0-run23
#  - added t-value and p-value to summary stats [14.03.23]
#    --> there are now 103 fields in summary stats
#  - added t-value and p-value calcs to partition() [14.03.23]
#    --> calculate t-value:
#        * t-value = parameter / std. error
#    --> calculate p-value:
#        * scipy.stats.t.pdf(-abs(t-value), df)
#          where df = degrees of freedom (len(x)-len(args))
#  - placed numpy.diag in try-catch block in partition [14.03.23]
#  - set output files to be placed in subdir 'out' [14.03.23]
# VERSION 1.8.1
#  - set r2 threshold to 0 for model selection [14.03.27]
# VERSION 1.8.2
#  - amended model selection linear model alpha validity range [14.04.01]
#    --> 0 < alpha < 1
#  - new model_selection criteria [14.04.03]
#    --> check three criteria: R2 thresh, p-value thresh, param validity
#    --> if more than one model meets all criteria, check for difference
#        between linear and hyperbolic (more than 1% difference) like in
#        Ruimy et al. 1995
#    --> outlier versus observation
#    --> problem finding model ...
# VERSION 1.8.3-run24
#  - added tower & month variables to partition() [14.04.04]
#  - added tower & month to FLUX_PARTI class [14.04.04]
#    --> initialized in partition() function
#    --> removed name & month from summary_statistics function call;
#        no longer needed
#  - updated model_selection criteria [14.04.04]
#    --> better account for default models when selection fails to find "best"
#  - added divide by zero check in calculate_stats in FLUX_PARTI [14.04.04]
#  - removed numpy_ones (unnecessary) from skew & kurt calcs [14.04.04]
# VERSION 1.8.4-run25
#  - moved gap_fill outside conditionals to run all months [14.06.24]
#    --> undone after model run
#  - added else to create blank FLUX_PARTI class for poor months [14.06.24]
#  - moved summary stats outside conditional to run all months [14.06.24]
# VERSION 1.9.0
#  - updated class & function doc [14.09.25]
#  - updated SOLAR class based on STASH 2.0 methods [14.09.26]
#  - updated calc_lue and lue_model functions for Colin's method [14.09.26]
#  - added density_h2o and viscosity_h2o functions [14.09.26]
# VERSION 1.9.1
#  - updated komega to float --- may influence cooper delta [14.09.29]
# VERSION 1.9.2
#  - updated Woolf's method in Solar---check lambda [14.10.07]
#  - added Berger's method to Solar class [14.10.07]
#  - corrected LUE model---Colin's method [14.10.07]
# VERSION 1.9.3
#  - added Spencer's method to Solar class [14.10.15]
# VERSION 1.9.4-run26
#  - modified runtime such that sum. stats writes for all stations [14.10.22]
#  - modified my_array check in FLUX_PARTI calc_statistics [14.10.22]
#  - added check for None type in add_station_val in LUE class [14.10.22]
#  - added check for alpha = -9999 in calc_lue [14.10.22]
# VERSION 2.0.0
#  - created version numbers [15.01.29]
#  - updated SOLAR class [15.01.29]
#    --> removed approx. methods (e.g., for dr, lambda, nu, and delta)
#    --> renamed class variables
#    --> suppressed a number of class variables (that aren't used)
#    --> ``stripped'' version of the SOLAR class in etsrad.py
#  - new LUE class (from gepisat_nlsr.py) [15.01.29]
#  - added next_gen_lue function (from gepisat_nlsr.py) [15.01.29]
#  - new calc_k function (from gepisat_nlsr.py) [15.01.29]
#  - updated calc_gstar [15.01.29]
#  - new calc_lue function (from gepisat_nlsr.py) [15.01.29]
#  - new viscosity_h2o function [15.01.29]
#    --> based on Vogel equation
#    --> removed density_h2o function
#  - updated get_pressure with STASH constants [15.01.29]
#  - moved patm calculation out of time-loop [15.01.29]
#  - changed units for VPD from kPa to Pa [15.01.29]
#  - updated gapfill_ppfd() function [15.01.30]
#    --> now outputs array of associated datetime objects
#  - updated calc_gpp() in FLUX_PARTI class for model select 0 [15.01.30]
#    --> try hyperbolic model estimates with outliers removed first or default
#        to linear model estimates with outliers removed
#  - updated Pearson's r calculation in FLUXPARTI class [15.01.30]
#    --> use scipy.stats.linregress()
#  - updated the order of functions [15.01.30]
#    --> base and dependant function sections
#    --> otherwise alphabetized
# VERSION 2.1.0
#  - updated LUE class [15.02.20]
#    --> removed basic_lue() function
#    --> moved calc_gstar function to LUE class
#    --> moved calc_k function to LUE class
#    --> moved next_gen_lue function to LUE class
#    --> moved predict_params to LUE class & renamed beta_estimate
#    --> moved viscosity_h2o to LUE class & updated with Huber method
#    --> removed calc_lue function (GPP estimates now made in LUE class)
#  - updated get_pressure function to also return elevation [15.02.20]
#  - fixed time issue in gapfill_ppfd [15.02.20]
#  - created calc_daily_gpp function [15.02.20]
# VERSION 2.1.1
#  - created gapfill_ppfd_day [15.02.23]
#  - renamed gapfill_ppfd to gapfill_ppfd_month [15.02.23]
#    --> moved daily processing & file writing to gapfill_ppfd_day
#  - updated SOLAR class with local_time class variable [15.02.23]
#  - created daily_gpp function [15.02.23]
#  - PEP8 style fixes [15.11.18]
# VERSION 2.2.0-dev
#  - separated individual class files [16.01.17]
#  - Python 2/3 consistency checks [16.01.17]
#  - imports connectSQL from db_setup [16.01.17]
#
# -----
# todo:
# -----
# - Create a constants file
# - Separate utility functions to their own file & import them
# - Add logging
# - Consider creating a data class to handle the interfact between the main
#   model and the necessary input variables
#   + this will allow alternative means of getting input data (e.g., read from
#     file) to be implemented
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
# X gapfill_ppfd
#    X what to do with sfactor when ppfd_integral = 0?
#    --> currently sfactor set to 1.0 (arbitrarily)
#    --> DOESN'T MATTER, POLAR NIGHT, ppfd_hh WILL BE ARRAY OF ZEROS
# X Implement specific model runs, e.g.:
#    X calculate daily GPP for a given day and station
#    X create gapfill_ppfd_day & gapfill_ppfd_month
#
###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import os.path

# import psycopg2
import numpy
import scipy.special
import scipy.stats
from scipy.optimize import curve_fit

from db_setup import connectSQL
from flux_parti import FLUX_PARTI
from lue import LUE


###############################################################################
# FUNCTIONS
###############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Base Functions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def add_one_day(dt0):
    """
    Name:     add_one_day
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt1)
    Features: Adds one day to datetime
    """
    dt1 = dt0 + datetime.timedelta(days=1)
    return dt1


def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3


def grid_centroid(my_lon, my_lat):
    """
    Name:     grid_centroid
    Input:    - float, longitude (my_lon)
              - float, latitude (my_lat)
    Output:   tuple, longitude latitude pair (my_centroid)
    Features: Returns the nearest 0.5 deg. grid centroid per given coordinates
              based on the Euclidean distance to each of the four surrounding
              grids; if any distances are equivalent, the pixel north and east
              is selected by default
    """
    # Create lists of regular latitude and longitude:
    grid_res = 0.5
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = 360
    lon_dim = 720
    lats = [lat_min + y * grid_res for y in xrange(lat_dim)]
    lons = [lon_min + x * grid_res for x in xrange(lon_dim)]

    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index-1]
        try:
            bb_lon_max = lons[lon_index+1]
        except IndexError:
            bb_lon_max = lons[-1] + grid_res

    # Find bounding latitude:
    centroid_lat = None
    if my_lat in lats:
        centroid_lat = my_lat
    else:
        lats.append(my_lat)
        lats.sort()
        lat_index = lats.index(my_lat)
        bb_lat_min = lats[lat_index-1]
        try:
            bb_lat_max = lats[lat_index+1]
        except IndexError:
            bb_lat_max = lats[-1] + grid_res

    # Determine nearest centroid:
    # NOTE: if dist_A equals dist_B, then centroid defaults positively
    #       i.e., north / east
    if centroid_lon and centroid_lat:
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lon and not centroid_lat:
        # Calculate the distances between lat and bounding box:
        dist_A = bb_lat_max - my_lat
        dist_B = my_lat - bb_lat_min
        if dist_A > dist_B:
            centroid_lat = bb_lat_min
        else:
            centroid_lat = bb_lat_max
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lat and not centroid_lon:
        # Calculate the distances between lon and bounding box:
        dist_A = bb_lon_max - my_lon
        dist_B = my_lon - bb_lon_min
        if dist_A > dist_B:
            centroid_lon = bb_lon_min
        else:
            centroid_lon = bb_lon_max
        my_centroid = (centroid_lon, centroid_lat)
    else:
        # Calculate distances between lat:lon and bounding box:
        # NOTE: if all distances are equal, defaults to NE grid
        dist_A = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        min_dist = min([dist_A, dist_B, dist_C, dist_D])

        # Determine centroid based on min distance:
        if dist_A == min_dist:
            my_centroid = (bb_lon_max, bb_lat_max)
        elif dist_B == min_dist:
            my_centroid = (bb_lon_max, bb_lat_min)
        elif dist_C == min_dist:
            my_centroid = (bb_lon_min, bb_lat_max)
        elif dist_D == min_dist:
            my_centroid = (bb_lon_min, bb_lat_min)

    # Return nearest centroid:
    return my_centroid


def simpson(my_array, h):
    """
    Name:     simpson
    Input:    - numpy.ndarray (my_array)
              - int, step-length (h)
    Features: Returns the numerical integration over an array by applying
              Simpson's (composite) rule for numerical integration, i.e.
              int(fx) = (h/3)*[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) ... f(xn)]
    Note:     n should be even, but if it is not, the first and last values of
              radiation curve are typically zero, so no harm done
    """
    n = len(my_array)
    s = my_array[0] + my_array[-1]

    for i in xrange(1, n, 2):
        s += 4.0*my_array[i]
    for j in xrange(2, n-1, 2):
        s += 2.0*my_array[j]
    s = (s*h)/3.0
    return s


def summary_file_init(summary_file):
    """
    Name:     summary_file_init
    Input:    str, output file (summary_file)
    Output:   None.
    Features: Creates a new summary file with header line
    """
    summary_header = (
        "%s,%s,%s,%s,%s,"        # 01--05
        "%s,%s,%s,%s,%s,"        # 06--10
        "%s,%s,%s,%s,%s,"        # 11--15
        "%s,%s,%s,%s,%s,"        # 16--20
        "%s,%s,%s,%s,%s,"        # 21--25
        "%s,%s,%s,%s,%s,"        # 26--30
        "%s,%s,%s,%s,%s,"        # 31--35
        "%s,%s,%s,%s,%s,"        # 36--40
        "%s,%s,%s,%s,%s,"        # 41--45
        "%s,%s,%s,%s,%s,"        # 46--50
        "%s,%s,%s,%s,%s,"        # 51--55
        "%s,%s,%s,%s,"           # 56--59
        "%s,%s,%s,%s,"           # 60--63
        "%s,%s,%s,%s,%s,%s,"     # 64--69
        "%s,%s,%s,%s,%s,%s,"     # 70--75
        "%s,%s,%s,%s,%s,%s,"     # 76--81
        "%s,%s,%s,%s,%s,%s,"     # 82--87
        "%s,%s,%s,%s,%s,%s,"     # 88--93
        "%s,%s,%s,%s,%s,%s,"     # 94--99
        "%s,%s,%s,%s\n"          # 100--103
        ) % (
        "name", "month",   "n_obs", "n_h",    "n_l",               # 01--05
        "foo_est_obs_h",   "foo_opt_obs_h",   "foo_err_obs_h",     # 06--08
        "foo_t_obs_h",     "foo_p_obs_h",                          # 09--10
        "foo_est_ro_h",    "foo_opt_ro_h",    "foo_err_ro_h",      # 11--13
        "foo_t_ro_h",      "foo_p_ro_h",                           # 14--15
        "alpha_est_obs_h", "alpha_opt_obs_h", "alpha_err_obs_h",   # 16--18
        "alpha_t_obs_h",   "alpha_p_obs_h",                        # 19--20
        "alpha_est_ro_h",  "alpha_opt_ro_h",  "alpha_err_ro_h",    # 21--23
        "alpha_t_ro_h",    "alpha_p_ro_h",                         # 24--25
        "alpha_est_obs_l", "alpha_opt_obs_l", "alpha_err_obs_l",   # 26--28
        "alpha_t_obs_l",   "alpha_p_obs_l",                        # 29--30
        "alpha_est_ro_l",  "alpha_opt_ro_l",  "alpha_err_ro_l",    # 31--33
        "alpha_t_ro_l",    "alpha_p_ro_l",                         # 34--35
        "r_est_obs_h",     "r_opt_obs_h",     "r_err_obs_h",       # 36--38
        "r_t_obs_h",       "r_p_obs_h",                            # 39--40
        "r_est_ro_h",      "r_opt_ro_h",      "r_err_ro_h",        # 41--43
        "r_t_ro_h",        "r_p_ro_h",                             # 44--45
        "r_est_obs_l",     "r_opt_obs_l",     "r_err_obs_l",       # 46--48
        "r_t_obs_l",       "r_p_obs_l",                            # 49--50
        "r_est_ro_l",      "r_opt_ro_l",      "r_err_ro_l",        # 51--53
        "r_t_ro_l",        "r_p_ro_l",                             # 54--55
        "r2_obs_h",        "r2_ro_h",                              # 56--57
        "rmse_obs_h",      "rmse_ro_h",                            # 58--59
        "r2_obs_l",        "r2_ro_l",                              # 60--61
        "rmse_obs_l",      "rmse_ro_l",                            # 62--63
        "min_ppfd_obs",    "max_ppfd_obs",    "ave_ppfd_obs",      # 64--66
        "std_ppfd_obs",    "skw_ppfd_obs",    "krt_ppfd_obs",      # 67--69
        "min_ppfd_ro_h",   "max_ppfd_ro_h",   "ave_ppfd_ro_h",     # 70--72
        "std_ppfd_ro_h",   "skw_ppfd_ro_h",   "krt_ppfd_ro_h",     # 73--75
        "min_ppfd_ro_l",   "max_ppfd_ro_l",   "ave_ppfd_ro_l",     # 76--78
        "std_ppfd_ro_l",   "skw_ppfd_ro_l",   "krt_ppfd_ro_l",     # 79--81
        "min_nee_obs",     "max_nee_obs",     "ave_nee_obs",       # 82--84
        "std_nee_obs",     "skw_nee_obs",     "krt_nee_obs",       # 85--87
        "min_nee_ro_h",    "max_nee_ro_h",    "ave_nee_ro_h",      # 88--90
        "std_nee_ro_h",    "skw_nee_ro_h",    "krt_nee_ro_h",      # 91--93
        "min_nee_ro_l",    "max_nee_ro_l",    "ave_nee_ro_l",      # 94--96
        "std_nee_ro_l",    "skw_nee_ro_l",    "krt_nee_ro_l",      # 97--99
        "pearson_r_obs",   "pearson_r_ro_h",  "pearson_r_ro_l",  # 100--102
        "model_select"                                           # 103
        )
    # Open file for writing:
    try:
        SFILE = open(summary_file, 'w')
    except IOError:
        print("Cannot open file '%s' for writting" % (summary_file))
        raise
    else:
        SFILE.write(summary_header)
        SFILE.close()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Dependant Functions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calc_daily_gpp(ts, gpp, gpp_err):
    """
    Name:     calc_daily_gpp
    Inputs:   - numpy.ndarray, half-hourly time stamps for a given month (ts)
              - numpy.ndarray, half-hourly gap-filled GPP, umol/m2/s (gpp)
              - numpy.ndarray, half-hourly GPP model errors (gpp_err)
    Output:   numpy.ndarray, multidimensional array
              > 'Timestamp' daily timestamps
              > 'GPP' daily GPP, mol/m2
              > 'GPP_err' daily GPP error, mol/m2
    Features: Returns time stamped daily GPP and its associated errors
    Depends:  - add_one_day
              - simpson

    @TODO: write to file option?
    """
    # Get starting and ending dates:
    starting_date = ts[0]
    ending_date = ts[-1]

    # Iterate through days:
    cur_date = starting_date
    while (cur_date < ending_date):
        # Find GPP values associated with current day:
        my_idx = numpy.where(
            (ts >= cur_date) & (ts < add_one_day(cur_date)))[0]
        my_gpp = [gpp[i] for i in my_idx]
        my_gpp_err = [gpp_err[i] for i in my_idx]

        # Integrate to daily:
        # @TODO: should this be clipped to a min of zero?
        day_gpp = simpson(numpy.array(my_gpp), 1800)
        day_gpp_err = simpson(numpy.array(my_gpp_err), 1800)

        # Convert from umol to moles:
        day_gpp *= 1e-6
        day_gpp_err *= 1e-6

        # Save results:
        if cur_date == starting_date:
            gpp_daily = numpy.array(
                (cur_date.date(), day_gpp, day_gpp_err),
                dtype={'names': ('Timestamp', 'GPP', 'GPP_err'),
                       'formats': ('O', 'f4', 'f4')},
                ndmin=1
            )
        else:
            temp_array = numpy.array(
                (cur_date.date(), day_gpp, day_gpp_err),
                dtype={'names': ('Timestamp', 'GPP', 'GPP_err'),
                       'formats': ('O', 'f4', 'f4')},
                ndmin=1
            )
            gpp_daily = numpy.append(gpp_daily, temp_array, axis=0)

        # Increment current day:
        cur_date = add_one_day(cur_date)

    return gpp_daily


def daily_gpp(station, cur_date):
    """
    Name:     daily_gpp
    Inputs:   - str, flux station name (station)
              - datetime.date, current date (cur_date)
    Output:   float, daily GPP, mol/m^2
              > 'nan' for months without observations (i.e., no partitioning)
    Features: Calculates the GPP and associated error for a given day and
              station
    Depends:  - monthly_ppfd_nee
              - partition
              - gapfill_ppfd_day
              - calc_daily_gpp
    """
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Get the monthly partitioning parameters
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Set datetime to start of the month:
    month_date = cur_date.replace(day=1)

    # Get monthly NEE and PPFD observations from database:
    mo_nee, mo_ppfd = monthly_ppfd_nee(station, month_date)

    # Calculate partitioning parameters:
    if mo_nee and mo_ppfd:
        mo_parti = partition(mo_nee, mo_ppfd, to_write=False, rm_out=True,
                             tower=station, month=month_date)

        # ~~~~~~~~~~~~~~~~~~~
        # Calculate daily GPP
        # ~~~~~~~~~~~~~~~~~~~
        # Get gapfilled PPFD [umol m^-2 m^-1]:
        (hh_ts, hh_ppfd) = gapfill_ppfd_day(station=station,
                                            cur_date=cur_date,
                                            to_write=False,
                                            out_file='/dev/null')

        # Convert PPFD to GPP [umol m^-2 s^-1]:
        hh_gpp, hh_gpp_err = mo_parti.calc_gpp(hh_ppfd)

        # Calculate daily GPP (mol/m2):
        day_ts, day_gpp, day_gpp_err = calc_daily_gpp(hh_ts, hh_gpp,
                                                      hh_gpp_err)[0]
    else:
        day_gpp = numpy.nan
        day_gpp_err = numpy.nan

    return (day_gpp, day_gpp_err)


def flux_to_grid(flux_station):
    """
    Name:     flux_to_grid
    Input:    str, station name (flux_station)
    Output:   int, grid station ID (grid_station)
    Features: Returns grid station ID based on the location of a given flux
              tower
    Depends:  - connectSQL
              - get_lon_lat
              - grid_centroid
    """
    # Get lat and lon of flux tower:
    (fst_lon, fst_lat) = get_lon_lat(flux_station)

    # Determine grid centroid lon and lat:
    (grd_lon, grd_lat) = grid_centroid(fst_lon, fst_lat)

    # Get grid station name based on centroid coordinates:
    params = ("grid", grd_lon, grd_lat)
    q = (
        "SELECT met_data.stationid "
        "FROM met_data "
        "WHERE met_data.geom = %s "
        "AND met_data.lon = %s "
        "AND met_data.lat = %s;"
        )

    # Connect to database:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and return results:
        cur.execute(q, params)
        grid_station = cur.fetchone()[0]
        con.close()

        return grid_station


def gapfill_ppfd_day(station, cur_date, to_write, out_file):
    """
    Name:     gapfill_ppfd_day
    Input:    - str, station name (station)
              - datetime.date (cur_date)
              - int, boolean for writing to file (to_write)
              - str, output file name (out_file)
    Output:   - numpy.ndarray of datetime objects (monthly_timestamp_hh)
              - numpy.ndarray of gapfilled PPFD (monthly_ppfd_gapless)
    Features: Returns array of half-hourly gapless PPFD (umol m-2 s-1) for a
              given station and day and associated timestamps; prints
              observations and gapfilled PPFD to file
    Depends:  - flux_to_grid
              - get_daily_ppfd
              - get_data_point
              - get_lon_lat
              - get_msvidx
              - SOLAR
    """
    # Check if output file exists:
    if to_write:
        if not os.path.isfile(out_file):
            header = "Timestamp,PPFDobs,PPFDgf\n"
            try:
                f = open(out_file, "w")
            except IOError:
                print("Cannot write to file '%s'" % (out_file))
            else:
                f.write(header)
                f.close()

    # Initialize monthly gapless PPFD array & associated timestamps:
    daily_ppfd_gapless = numpy.array([])
    daily_timestamp_hh = numpy.array([])

    # Get daily PPFD values (in a numpy.array) from database:
    (daily_ts, daily_ppfd) = get_daily_ppfd(station, cur_date)

    # Check to see if any gaps are present in current day's observations:
    number_obs = len(daily_ppfd)
    if number_obs < 49:
        # Gaps are present:
        gapfill_dict = {}

        # Create start and end datetime objects for this day:
        start_time = datetime.datetime(cur_date.year,
                                       cur_date.month,
                                       cur_date.day,
                                       0, 0, 0)

        end_time = datetime.datetime(cur_date.year,
                                     cur_date.month,
                                     cur_date.day,
                                     23, 59, 59)

        # Initialize dictionary values with half-hourly timestamp keys:
        cur_time = start_time
        while cur_time < end_time:
            my_time = "%s" % cur_time.time()
            gapfill_dict[my_time] = 0.0

            # Add datetime object to array of monthly timestamps:
            daily_timestamp_hh = numpy.append(daily_timestamp_hh, [cur_time, ])

            cur_time = cur_time + datetime.timedelta(minutes=30)

        # Convert date to Julian day:
        jday = cur_date.timetuple().tm_yday

        # Calculate daily ET solar radiation curve:
        (flux_lon, flux_lat) = get_lon_lat(station)
        et_solar = SOLAR(flux_lon, flux_lat, jday, cur_date.year)

        # Get satellite measurement of solar rad (SWdown) [W m-2]:
        grid_station = flux_to_grid(station)
        grid_msvidx = get_msvidx(grid_station, 'SWdown')
        grid_srad = get_data_point(grid_msvidx, cur_date)

        # Convert to daily shortwave radiation [J m-2]:
        # NOTE: if None, then gridded data not available, set equal to modeled
        if grid_srad is not None:
            grid_srad *= (86400.0)
        else:
            grid_srad = et_solar.ho_jm2

        # Calculate scaling factor (i.e., observed/modeled):
        if et_solar.ho_jm2 != 0:
            sfactor = grid_srad/et_solar.ho_jm2
        else:
            sfactor = 1.0

        # Add scaled half-hourly PPFD to dictionary [umol m-2 s-1]:
        for i in xrange(48):
            val = et_solar.ppfd_hh[i]
            my_time = "%s" % et_solar.local_time[i].time()
            gapfill_dict[my_time] = sfactor*val

        # Add observations to gapfill dictionary & save for output:
        ppfd_obs = {}
        for i in xrange(number_obs):
            my_time = "%s" % daily_ts[i].time()
            gapfill_dict[my_time] = daily_ppfd[i]
            ppfd_obs[my_time] = daily_ppfd[i]

        # Save to PPFD time series:
        daily_ppfd_gapless = numpy.append(
            daily_ppfd_gapless,
            [gapfill_dict[x] for x in sorted(gapfill_dict.keys())]
        )

        # Write to file:
        if to_write:
            for t in sorted(gapfill_dict.keys()):
                if t in ppfd_obs.keys():
                    obs = ppfd_obs[t]
                else:
                    obs = -9999.0
                gfv = gapfill_dict[t]
                dt = "%s %s" % (cur_date, t)
                try:
                    f = open(out_file, 'a')
                except IOError:
                    print("Cannot append to file '%s'" % (out_file))
                else:
                    f.write("%s,%f,%f\n" % (dt, obs, gfv))
                    f.close()
    else:
        # No gaps; append daily series
        #   NOTE: drop last entry from daily_ppfd (midnight next day)
        daily_ppfd_gapless = daily_ppfd[0:-1]
        daily_timestamp_hh = daily_ts[0:-1]

        # Write to file:
        if to_write:
            for i in xrange(len(daily_ts) - 1):
                dt = "%s" % daily_ts[i]
                obs = daily_ppfd[i]
                try:
                    f = open(out_file, 'a')
                except IOError:
                    print("Cannot append to file '%s'" % (out_file))
                else:
                    f.write("%s,%0.3f,%0.3f\n" % (dt, obs, obs))
                    f.close()

    return (daily_timestamp_hh, daily_ppfd_gapless)


def gapfill_ppfd_month(station, start_date, to_write):
    """
    Name:     gapfill_ppfd_month
    Input:    - str, station name (station)
              - datetime.date (start_date)
              - int, write to file boolean (to_write)
    Output:   - numpy.ndarray of datetime objects (monthly_timestamp_hh)
              - numpy.ndarray of gapfilled PPFD (monthly_ppfd_gapless)
    Features: Returns array of half-hourly gapless PPFD (umol m-2 s-1) for a
              given station and month and associated timestamps; prints
              gap-filled PPFD to file
    Depends:  - add_one_day
              - add_one_month
              - gapfill_ppfd_day
    """
    # Initialize monthly gapless PPFD array & associated timestamps:
    monthly_ppfd_gapless = numpy.array([])
    monthly_timestamp_hh = numpy.array([])

    # Calculate the end date:
    end_date = add_one_month(start_date)

    # Create output file (if to_write):
    out_file = "out/%s-GF_%s.txt" % (station, start_date)

    # Iterate through each day of the month:
    cur_date = start_date
    while cur_date < end_date:
        # Gapfill daily PPFD
        (gf_daily_time, gf_daily_ppfd) = gapfill_ppfd_day(station=station,
                                                          cur_date=cur_date,
                                                          to_write=to_write,
                                                          out_file=out_file)
        # Append data to monthly arrays:
        monthly_timestamp_hh = numpy.append(monthly_timestamp_hh,
                                            gf_daily_time)
        monthly_ppfd_gapless = numpy.append(monthly_ppfd_gapless,
                                            gf_daily_ppfd)

        # Increment day
        cur_date = add_one_day(cur_date)

    return (monthly_timestamp_hh, monthly_ppfd_gapless)


def get_daily_ppfd(station, start_date):
    """
    Name:     get_daily_ppfd
    Input:    - str, station name (station)
              - datetime.date, date of interest (start_date)
    Output:   tuple, arrays of time stamps and PPFD data
              - numpy.ndarray, timestamps (time_vals)
              - numpy.ndarray, associated PPFD (ppfd_vals)
    Features: Returns the half-hourly timestamps and PPFD observations for a
              given day
    Depends:  - add_one_day
              - connectSQL
              - get_msvidx
    """
    # Get msvidx value for PPFD:
    ppfd_idx = get_msvidx(station, 'PPFD_f')

    # SQL query parameters:
    params = (ppfd_idx, start_date, add_one_day(start_date))

    # Define query:
    q = (
        "SELECT data_set.datetime, data_set.data "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "AND data_set.datetime BETWEEN DATE %s AND DATE %s "
        "ORDER BY data_set.datetime ASC;"
        )

    # Connect to database and start a cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and store results:
        cur.execute(q, params)
        ppfd_vals = numpy.array([])
        time_vals = numpy.array([])
        if cur.rowcount > 0:
            for record in cur:
                time_vals = numpy.append(time_vals, record[0])
                ppfd_vals = numpy.append(ppfd_vals, record[1])

        # Close connection and return results:
        con.close()

        return (time_vals, ppfd_vals)


def get_data_point(msvidx, time_point):
    """
    Name:     get_data_points
    Input:    - str, msvidx (msvidx)
              - datetime.date (time_point)
    Output:   float/numpy.ndarray (my_result)
    Features: Returns data point or array of data for a given msvidx (i.e.,
              station and variable) and time
    Depends:  connectSQL
    """
    # SQL query params:
    params = (msvidx, time_point)

    # Define SQL query:
    q = (
        "SELECT data_set.data "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "AND data_set.datetime = %s;"
        )

    # Connect to database and start a cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and return result:
        cur.execute(q, params)
        if cur.rowcount == 1:
            my_result = cur.fetchone()[0]
        elif cur.rowcount > 1:
            my_result = numpy.array([])
            for record in cur:
                my_result = numpy.append(my_result, record[0])
        else:
            my_result = None
            print("No data found in function get_data_point")

        return my_result


def get_dates(station):
    """
    Name:     get_dates
    Input:    str, station name (station)
    Output:   tuple, starting and ending dates
              - datetime.date, starting date (sd)
              - datetime.date, ending date (ed)
    Features: Returns the starting and ending dates of NEE-PPFD data pairs for
              a given station
    Depends:  - connectSQL
              - get_msvidx
    """
    # Get msvidx values for specified station:
    ppfdi = get_msvidx(station, 'PPFD_f')
    neei = get_msvidx(station, 'NEE_f')

    # SQL query parameters:
    params = (ppfdi, neei)

    # Define start date query:
    q1 = (
        "SELECT data_set.datetime "
        "FROM data_set "
        "WHERE data_set.msvidx = %s OR data_set.msvidx = %s "
        "ORDER BY data_set.datetime ASC LIMIT 1;"
        )

    # Define end date query:
    q2 = (
        "SELECT data_set.datetime "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "OR data_set.msvidx = %s "
        "ORDER BY data_set.datetime DESC LIMIT 1;"
        )

    # Connect to database and start a cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Get start date from datetime object:
        cur.execute(q1, params)
        sd = cur.fetchone()[0].date()

        # Get end date from datetime object:
        cur.execute(q2, params)
        ed = cur.fetchone()[0].date()

        # Make the starting date begin at day 1
        sd = sd.replace(day=1)

        # Return results:
        con.close()

        return (sd, ed)


def get_lon_lat(station):
    """
    Name:     get_lon_lat
    Input:    str, station name (station)
    Output:   tuple, lon-lat pair
              - float, longitude (my_lon)
              - float, latitude (my_lat)
    Features: Return longitude and latitude pair for a given station based on
              the GePiSaT database meta-data table
    Depends:  connectSQL
    """
    # Query paramters:
    params = (station,)

    # SQL query:
    q = (
        "SELECT met_data.lon, met_data.lat "
        "FROM met_data "
        "WHERE met_data.stationid = %s;"
        )

    # Connect to database and start a cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and return results:
        cur.execute(q, params)
        my_lon, my_lat = cur.fetchone()
        con.close()
        return (my_lon, my_lat)


def get_msvidx(station, variable):
    """
    Name:     get_msvidx
    Input:    - str, station name (station)
              - str, variable name (variable)
    Output:   string, msvidx (result)
    Features: Returns the msvidx from the GePiSaT database based on the station
              and variable name
    Depends:  connectSQL
    """
    # Define query:
    q = (
        "SELECT var_list.msvidx "
        "FROM var_list "
        "WHERE var_list.stationid = %s "
        "AND var_list.varname = %s;"
        )

    # SQL query parameters:
    params = (station, variable)

    # Connect to database and star cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and fetch results:
        cur.execute(q, params)
        try:
            result = cur.fetchone()[0]
        except:
            print(("Could not return an msvidx value for station '%s' and "
                   "variable '%s'") % (station, variable))
            result = ""
        finally:
            con.close()
            return result


def get_pressure(s):
    """
    Name:     get_pressure
    Input:    str, station name (s)
    Output:   - float, elevation, m (z)
              - float, atmospheric pressure, Pa (patm)
    Features: Returns the atmospheric pressure based on the elevation of a
              given station
    Depends:  - connectSQL
              - flux_to_grid
              - get_data_point
              - get_msvidx
    Ref:      Allen et al. (1998)
    """
    # Define constants:
    kPo = 101325    # standard atmosphere, Pa (Allen, 1973)
    kTo = 288.15    # base temperature, K (Berberan-Santos et al., 1997)
    kL = 0.0065     # temperature lapse rate, K/m (Allen, 1973)
    kG = 9.80665    # gravitational acceleration, m/s^2 (Allen, 1973)
    kR = 8.3143     # universal gas constant, J/mol/K (Allen, 1973)
    kMa = 0.028963  # molecular weight of dry air, kg/mol (Tsilingiris, 2008)

    # Define query w/ parameters:
    params = (s,)
    q = (
        "SELECT met_data.ele "
        "FROM met_data "
        "WHERE met_data.stationid = %s;"
        )

    # Connect to database:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and return results:
        cur.execute(q, params)
        station_ele = cur.fetchone()[0]
        con.close()
    else:
        station_ele = -9999

    # Check to see that elevation is valid:
    if float(station_ele) == -9999:
        # Find CRU Elv:
        elv_sd = datetime.date(2006, 6, 1)
        hdg_station = flux_to_grid(s)
        elv_msvidx = get_msvidx(hdg_station, 'Elv')
        elv_data = get_data_point(elv_msvidx, elv_sd)
        station_ele = elv_data

    # Convert elevation to pressure, Pa:
    z = float(station_ele)
    patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))

    return (z, patm)


def get_stations():
    """
    Name:     get_stations
    Input:    None.
    Output:   list, station names (results)
    Features: Returns a list of flux station names from GePiSaT database
    Depends:  connectSQL
    """
    # Define query:
    q = (
        "SELECT stationid "
        "FROM met_data "
        "WHERE dim=0 "
        "AND geom=%s "
        "ORDER BY stationid ASC;"
        )

    params = ("point",)

    # Connect to database and start cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and fetch results:
        cur.execute(q, params)
        results = []
        for record in cur:
            results.append(record[0])  # <- extract record from tuple

        con.close()

        return results


def monthly_ppfd_nee(station, start_date):
    """
    Name:     monthly_ppfd_nee
    Input:    - str, station name (station)
              - datetime.date (start_date)
    Output:   tuple, PPFD and NEE observations
              - numpy.ndarray, PPFD (ppfd_vals)
              - numpy.ndarray, NEE (nee_vals)
    Features: Returns one month of PPFD-NEE observation pairs for a given
              station and month
    Depends:  - add_one_month
              - connectSQL
              - get_msvidx
    """
    # Get msvidx values for specified station
    ppfd_idx = get_msvidx(station, 'PPFD_f')
    nee_idx = get_msvidx(station, 'NEE_f')

    # Increment start date one month:
    end_date = add_one_month(start_date)

    # SQL query parameters:
    params = (ppfd_idx, nee_idx, start_date, end_date, ppfd_idx, nee_idx)

    # Define query
    # NOTE: okay to use string concatenation % because ppfd and nee are
    # function return values:
    q = (
        "SELECT * "
        "FROM crosstab('"
        "select data_set.datetime, data_set.msvidx, data_set.data "
        "from data_set "
        "where data_set.msvidx = ''%s'' "
        "or data_set.msvidx = ''%s'' "
        "and data_set.datetime between date ''%s'' and date ''%s'' "
        "order by 1,2', "
        "'select distinct data_set.msvidx from data_set "
        "where data_set.msvidx = ''%s'' "
        "or data_set.msvidx = ''%s'' order by 1') "
        "AS ct(row_name TIMESTAMP, ppfd FLOAT, nee FLOAT);"
        ) % params

    # Connect to database and start a cursor:
    con = connectSQL()
    if con is not None:
        cur = con.cursor()

        # Execute query and fetch results:
        cur.execute(q)
        ppfd_vals = numpy.array([])
        nee_vals = numpy.array([])
        if cur.rowcount > 0:
            for record in cur:
                ppfd = record[1]
                nee = record[2]

                # Only save matched pairs:
                if (ppfd and nee):
                    nee_vals = numpy.append(nee_vals, nee)
                    ppfd_vals = numpy.append(ppfd_vals, ppfd)

        con.close()

        return (ppfd_vals, nee_vals)


def partition(nee, ppfd, to_write, rm_out, tower, month):
    """
    Name:     partition
    Input:    - numpy.ndarray, month of NEE observations (nee)
              - numpy.ndarray, month of PPFD observations (ppfd)
              - int, write to file boolean (to_write)
              - int, outlier removal boolean (rm_out)
              - str, station name (tower)
              - datetime.date (month)
    Output:   FLUX_PARTI class object (my_class)
    Features: Returns a class with flux partitioning results based on NEE-PPFD
              observations; processes outliers (optional); saves results to
              file (optional)
    Depends:  FLUX_PARTI
    """
    # Create a FLUX_PARTI class with PPFD and NEE arrays:
    my_class = FLUX_PARTI(ppfd, nee, tower, month)

    ## ~~~~~~~~~~ ANALYZE OBSERVATIONS ~~~~~~~~~~ ##
    # ---------------------
    # Model H optimization:
    # ---------------------
    try:
        mh_opt, mh_cov = curve_fit(
            my_class.model_h,
            my_class.ppfd_obs,
            my_class.nee_obs,
            p0=my_class.hm_estimates
            )
    except:
        my_class.hm_optimized = [-9999.0, -9999.0, -9999.0]
        my_class.hm_optim_err = [-9999.0, -9999.0, -9999.0]
    else:
        (
            my_class.hm_optimized[0],
            my_class.hm_optimized[1],
            my_class.hm_optimized[2]) = mh_opt

        # Extract variance values from matrix:
        try:
            mh_var = numpy.diag(mh_cov)
        except ValueError:
            mh_var = [0, 0, 0]
        else:
            if numpy.isfinite(mh_var).all() and not (mh_var < 0).any():
                (
                    my_class.hm_optim_err[0],
                    my_class.hm_optim_err[1],
                    my_class.hm_optim_err[2]
                    ) = numpy.sqrt(mh_var)

                # Calculate the t-value:
                (
                    my_class.hm_optim_t[0],
                    my_class.hm_optim_t[1],
                    my_class.hm_optim_t[2]
                    ) = (numpy.array(my_class.hm_optimized) /
                         numpy.array(my_class.hm_optim_err))

                # Calculate the p-value:
                hm_df = len(my_class.nee_obs) - 3.0  # degrees of freedom
                (my_class.hm_optim_p[0],
                 my_class.hm_optim_p[1],
                 my_class.hm_optim_p[2]) = scipy.stats.t.pdf(
                    -abs(numpy.array(my_class.hm_optim_t)), hm_df)
            else:
                my_class.hm_optim_err[0] = 0.0
                my_class.hm_optim_err[1] = 0.0
                my_class.hm_optim_err[2] = 0.0
    finally:
        my_class.calc_model_h()

    # ---------------------
    # Model L optimization:
    # ---------------------
    try:
        ml_opt, ml_cov = curve_fit(
            my_class.model_l,
            my_class.ppfd_obs,
            my_class.nee_obs,
            p0=my_class.lm_estimates
            )
    except:
        my_class.lm_optimized = [-9999.0, -9999.0]
        my_class.lm_optim_err = [-9999.0, -9999.0]
    else:
        (
            my_class.lm_optimized[0],
            my_class.lm_optimized[1]) = ml_opt

        # Extract variance values from matrix:
        try:
            ml_var = numpy.diag(ml_cov)
        except ValueError:
            ml_var = [0, 0]
        else:
            if numpy.isfinite(ml_var).all() and not (ml_var < 0).any():
                (
                    my_class.lm_optim_err[0],
                    my_class.lm_optim_err[1]
                    ) = numpy.sqrt(ml_var)

                # Calculate the t-value:
                (
                    my_class.lm_optim_t[0],
                    my_class.lm_optim_t[1]
                    ) = (numpy.array(my_class.lm_optimized) /
                         numpy.array(my_class.lm_optim_err))

                # Calculate the p-value:
                lm_df = len(my_class.nee_obs) - 2  # degrees of freedom
                (my_class.lm_optim_p[0],
                 my_class.lm_optim_p[1]) = scipy.stats.t.pdf(
                    -abs(numpy.array(my_class.lm_optim_t)), lm_df)
            else:
                my_class.lm_optim_err[0] = 0.0
                my_class.lm_optim_err[1] = 0.0
    finally:
        my_class.calc_model_l()

    ## ~~~~~~~~~~ WRITE OBSERVATION RESULTS ~~~~~~~~~~ ##
    if to_write:
        output_file = "out/%s_%s.txt" % (tower, month)
        header0 = "MH_guess,%f,%f,%f\n" % (
            my_class.hm_estimates[1],
            my_class.hm_estimates[2],
            my_class.hm_estimates[0]
            )
        header1 = "ML_guess,%f,%f\n" % (
            my_class.lm_estimates[0],
            my_class.lm_estimates[1]
            )
        header2 = "MH_opt,%f,%f,%f,%f,%f\n" % (
            my_class.hm_optimized[1],
            my_class.hm_optimized[2],
            my_class.hm_optimized[0],
            my_class.hm_rmse,
            my_class.hm_rsqr
            )
        header3 = "ML_opt,%f,%f,,%f,%f\n" % (
            my_class.lm_optimized[0],
            my_class.lm_optimized[1],
            my_class.lm_rmse,
            my_class.lm_rsqr
            )
        OUTFILE = open(output_file, 'w')
        OUTFILE.write(header0)
        OUTFILE.write(header1)
        OUTFILE.write(header2)
        OUTFILE.write(header3)
        OUTFILE.write("ppfd_obs,nee_obs,nee_mh,nee_ml\n")
        for po, no, nh, nl in map(None, my_class.ppfd_obs,
                                  my_class.nee_obs,
                                  my_class.nee_model_h,
                                  my_class.nee_model_l):
            outline = "%s,%s,%s,%s\n" % (po, no, nh, nl)
            OUTFILE.write(outline)
        OUTFILE.close()
    ## ~~~~~~~~~~ REMOVE OUTLIERS AND RE-ANALYZE ~~~~~~~~~~ ##
    if rm_out:
        # ---------------------
        # Model H optimization:
        # ---------------------
        my_class.remove_mh_outliers()
        try:
            mh_opt_ro, mh_cov_ro = curve_fit(
                my_class.model_h,
                my_class.ppfd_obs_h_ro,
                my_class.nee_obs_h_ro,
                p0=my_class.hm_estimates
                )
        except:
            my_class.hm_optimized_ro = [-9999.0, -9999.0, -9999.0]
            my_class.hm_optim_err_ro = [-9999.0, -9999.0, -9999.0]
        else:
            (
                my_class.hm_optimized_ro[0],
                my_class.hm_optimized_ro[1],
                my_class.hm_optimized_ro[2]) = mh_opt_ro

            # Extract variance values from matrix:
            try:
                mh_var_ro = numpy.diag(mh_cov_ro)
            except ValueError:
                mh_var_ro = [0, 0, 0]
            else:
                if (numpy.isfinite(mh_var_ro).all() and not
                        (mh_var_ro < 0).any()):
                    (
                        my_class.hm_optim_err_ro[0],
                        my_class.hm_optim_err_ro[1],
                        my_class.hm_optim_err_ro[2]
                        ) = numpy.sqrt(mh_var_ro)

                    # Calculate the t-value:
                    (
                        my_class.hm_optim_t_ro[0],
                        my_class.hm_optim_t_ro[1],
                        my_class.hm_optim_t_ro[2]
                        ) = (numpy.array(my_class.hm_optimized_ro) /
                             numpy.array(my_class.hm_optim_err_ro))

                    # Calculate the p-value:
                    hm_df_ro = len(my_class.nee_obs) - my_class.hm_outliers - 3
                    (my_class.hm_optim_p_ro[0],
                     my_class.hm_optim_p_ro[1],
                     my_class.hm_optim_p_ro[2]) = scipy.stats.t.pdf(
                        -abs(numpy.array(my_class.hm_optim_t_ro)), hm_df_ro)
                else:
                    my_class.hm_optim_err_ro[0] = 0.0
                    my_class.hm_optim_err_ro[1] = 0.0
                    my_class.hm_optim_err_ro[2] = 0.0
        finally:
            my_class.calc_model_h_ro()

        # ---------------------
        # Model L optimization:
        # ---------------------
        my_class.remove_ml_outliers()
        try:
            ml_opt_ro, ml_cov_ro = curve_fit(
                my_class.model_l,
                my_class.ppfd_obs_l_ro,
                my_class.nee_obs_l_ro,
                p0=my_class.lm_estimates
                )
        except:
            my_class.lm_optimized_ro = [-9999.0, -9999.0]
            my_class.lm_optim_err_ro = [-9999.0, -9999.0]
        else:
            (
                my_class.lm_optimized_ro[0],
                my_class.lm_optimized_ro[1]) = ml_opt_ro

            # Extract variance values from matrix:
            try:
                ml_var_ro = numpy.diag(ml_cov_ro)
            except ValueError:
                ml_var_ro = [0, 0]
            else:
                if (numpy.isfinite(ml_var_ro).all() and not
                        (ml_var_ro < 0).any()):
                    (
                        my_class.lm_optim_err_ro[0],
                        my_class.lm_optim_err_ro[1]
                        ) = numpy.sqrt(ml_var_ro)

                    # Calculate the t-value:
                    (
                        my_class.lm_optim_t_ro[0],
                        my_class.lm_optim_t_ro[1]
                        ) = (numpy.array(my_class.lm_optimized_ro) /
                             numpy.array(my_class.lm_optim_err_ro))

                    # Calculate the p-value:
                    lm_df_ro = len(my_class.nee_obs) - my_class.lm_outliers - 2
                    (my_class.lm_optim_p_ro[0],
                     my_class.lm_optim_p_ro[1]) = scipy.stats.t.pdf(
                        -abs(numpy.array(my_class.lm_optim_t_ro)), lm_df_ro)
                else:
                    my_class.lm_optim_err_ro[0] = 0.0
                    my_class.lm_optim_err_ro[1] = 0.0
        finally:
            my_class.calc_model_l_ro()

        ## ~~~~~~~~~~ WRITE OUTLIER-FREE RESULTS ~~~~~~~~~~ ##
        if to_write:
            output_file = "out/%s_%s_ro.txt" % (tower, month)
            header0 = "MH_guess,%0.5f,%0.2f,%0.2f\n" % (
                my_class.hm_estimates[1],
                my_class.hm_estimates[2],
                my_class.hm_estimates[0])
            header1 = "ML_guess,%0.5f,%0.2f\n" % (
                my_class.lm_estimates[0],
                my_class.lm_estimates[1])
            header2 = "MH_opt,%0.5f,%0.2f,%0.2f,%0.2f,%0.3f\n" % (
                my_class.hm_optimized_ro[1],
                my_class.hm_optimized_ro[2],
                my_class.hm_optimized_ro[0],
                my_class.hm_rmse_ro,
                my_class.hm_rsqr_ro)
            header3 = "ML_opt,%0.5f,%0.2f,,%0.2f,%0.3f\n" % (
                my_class.lm_optimized_ro[0],
                my_class.lm_optimized_ro[1],
                my_class.lm_rmse_ro,
                my_class.lm_rsqr_ro)
            OUTFILE = open(output_file, 'w')
            OUTFILE.write(header0)
            OUTFILE.write(header1)
            OUTFILE.write(header2)
            OUTFILE.write(header3)
            OUTFILE.write(
                "ppfd_obs_h,nee_obs_h,nee_mod_h,"
                "ppfd_obs_l,nee_obs_l,nee_mod_l\n"
                )
            for (
                    poh, noh, nhh, pol, nol, nll
                    ) in map(None, my_class.ppfd_obs_h_ro,
                             my_class.nee_obs_h_ro,
                             my_class.nee_model_h_ro,
                             my_class.ppfd_obs_l_ro,
                             my_class.nee_obs_l_ro,
                             my_class.nee_model_l_ro):
                outline = "%s,%s,%s,%s,%s,%s\n" % (
                    poh, noh, nhh,
                    pol, nol, nll
                    )
                OUTFILE.write(outline)
            OUTFILE.close()

    # Run model selection criteria:
    my_class.model_selection()

    # Return the FLUX_PARTI class object:
    return my_class

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == "__main__":
    # Define output directory:
    output_dir = "out"

    # Get list of all flux station names:
    #stations = get_stations()
    stations = ['CZ-wet', ]

    # Initialize summary statistics file:
    summary_file = os.path.join(output_dir, "summary_statistics.txt")
    summary_file_init(summary_file)

    # Create/initialize LUE class instance:
    my_lue = LUE()
    lue_out_file = os.path.join(output_dir, "LUE_All-Stations.txt")

    # Iterate through stations:
    for station in stations:
        # Initialize station's LUE & daily GPP output file:
        lue_file = "%s_%s.txt" % (station, "LUE")
        lue_file = os.path.join(output_dir, lue_file)
        gpp_file = "%s_%s.txt" % (station, "GPP-daily")
        gpp_file = os.path.join(output_dir, gpp_file)

        # Get first/last dates for station data:
        sd, ed = get_dates(station)

        # Get flux station's corresponding 0.5-degree grid station:
        hdg_station = flux_to_grid(station)

        # Get station's variable msvidx values:
        co2_msvidx = "US-MLO.21"
        cpa_msvidx = get_msvidx(hdg_station, 'alpha')
        fpar_msvidx = get_msvidx(hdg_station, 'FAPAR')
        tair_msvidx = get_msvidx(hdg_station, 'Tc')
        vpd_msvidx = get_msvidx(hdg_station, 'VPD')

        # Get station's elevation & atmospheric pressure:
        elv, patm = get_pressure(station)

        # Process each month in time:
        while sd < ed:
            # Get PPFD and NEE array pairs [umol m-2 s-1]:
            #   NOTE: accepts in opposite order as sent
            monthly_nee, monthly_ppfd = monthly_ppfd_nee(station, sd)

            # Process if enough data was found:
            if (len(monthly_ppfd) > 3 and len(monthly_nee) > 3):
                # Perform GPP partitioning (returns FLUX_PARTI class object):
                monthly_parti = partition(monthly_nee,
                                          monthly_ppfd,
                                          to_write=False,
                                          rm_out=True,
                                          tower=station,
                                          month=sd)

                # Perform half-hourly PPFD gapfilling (umol m-2 s-1):
                (gf_time, gf_ppfd) = gapfill_ppfd_month(station,
                                                        sd,
                                                        to_write=False)

                # Calculate half-hourly GPP (umol m-2 s-1)
                gf_gpp, gf_gpp_err = monthly_parti.calc_gpp(gf_ppfd)

                # Calculate daily GPP (mol m-2):
                # gpp_daily = calc_daily_gpp(gf_time, gf_gpp, gf_gpp_err)

                # Continue processing if partitioning was successful:
                if monthly_parti.mod_select > 0:
                    # The new LUE model:
                    # GPP = f(PPFD, fAPAR, VPD, CPA, Tair, Patm, CO2)
                    #  - monthly variables: PPFD, fAPAR, CPA, VPD, Tair
                    #  - annual variables: CO2
                    #  - constants in time: Patm=f(elevation)

                    # Integrate PPFD & GPP [umol m-2]; dt=30 min (1800 s)
                    ppfd_month = simpson(gf_ppfd.clip(min=0), 1800)
                    gpp_month = simpson(gf_gpp.clip(min=0), 1800)
                    gpp_month_err = simpson(gf_gpp_err, 1800)

                    # Convert units from [umol m-2] to [mol m-2]:
                    ppfd_month *= (1e-6)
                    gpp_month *= (1e-6)
                    gpp_month_err *= (1e-6)

                    # Retrieve annual CO2:
                    annual_sd = sd.replace(month=1)
                    co2_annual = get_data_point(co2_msvidx, annual_sd)  # ppm

                    # Retrieve monthly gridded data:
                    cpa_month = get_data_point(cpa_msvidx, sd)      # unitless
                    fpar_month = get_data_point(fpar_msvidx, sd)    # unitless
                    tair_month = get_data_point(tair_msvidx, sd)    # deg C
                    vpd_month = get_data_point(vpd_msvidx, sd)      # kPa

                    # Add LUE parameters to LUE class:
                    my_lue.add_station_val(station,
                                           sd,
                                           gpp_month,               # mol/m2
                                           gpp_month_err,           # mol/m2
                                           fpar_month,              # unitless
                                           ppfd_month,              # mol/m2
                                           vpd_month,               # kPa
                                           cpa_month,               # unitless
                                           tair_month,              # degC
                                           co2_annual,              # ppm
                                           patm,                    # Pa
                                           elv)                     # m
            else:
                # Create an 'empty' class:
                monthly_parti = FLUX_PARTI(
                    monthly_ppfd, monthly_nee, station, sd)

            # Save class summary statistics:
            SFILE = open(summary_file, 'a')
            SFILE.write(monthly_parti.summary_statistics())
            SFILE.close()

            # Increment date:
            sd = add_one_month(sd)

        # Write monthly LUE parameters to file:
        my_lue.write_out_val(station, lue_file)
