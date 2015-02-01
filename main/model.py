#!/usr/bin/python
#
# model.py 
# __version__ 2.0.0
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-07-05 -- created
# 2015-01-31 -- last updated
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
#    --> source: http://chemwiki.ucdavis.edu/Analytical_Chemistry/               \
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
#
# -----
# todo:
# -----
# xx. gapfill_ppfd
#     a. what to do with sfactor when ppfd_integral = 0?
#     --> currently sfactor set to 1.0 (arbitrarily)
#     --> DOESN'T MATTER, POLAR NIGHT, ppfd_hh WILL BE ARRAY OF ZEROS 
#
# 01. FLUX_PARTI class
#     a. Why not use the optimized obs parameters as the guess for ro?
#        * if optim fails---you'll have crumby initial guesses
#     b. add variance of residuals to summary statistics (4 fields)
#        * goodness_of_fit() function
#
# 02. partition()
#     a. check data validity after outliers are removed before reprocessing
#
# 03. model_select()
#     a. Consider implementing either Shapiro-Wilks or Anderson test of 
#        normality on model residuals
#     --> scipy.stats.shapiro(my_resids)
#         * returns test statistic and p-value
#         * Note: p-value indicates signifantly different from normal
#     --> scipy.stats.anderson(my_resids, 'norm')
#         * returns test statistic, array of critical values, and sig values
#         * if test statistic is larger than the critical value at the 
#           significance value you are interested in, then the null can be
#           rejected, i.e., not normal
#
# 04. Check for system closure (each station, each month?)
#     --> short equation: Rn + G + LE + H = 0
#
################################################################################
## IMPORT MODULES 
################################################################################
from sys import exit
import datetime
import os.path
import psycopg2
import numpy
import scipy.special
import scipy.stats
from scipy.optimize import curve_fit

################################################################################
## CLASSES 
################################################################################
class FLUX_PARTI:
    """
    Name:     FLUX_PARTI
    Features: This class performs flux partitioning of monthly NEE & PPFD 
              observations
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    # Flux tower name & month being processed:
    name = ""
    month = datetime.date(1999,1,1)
    #
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
    #
    # PPFD Statistics (working, observation backup, & model outlier backups):
    max_ppfd = 0.0
    min_ppfd = 0.0
    ave_ppfd = 0.0
    std_ppfd = 0.0
    skw_ppfd = 0.0
    krt_ppfd = 0.0
    #
    max_ppfd_obs = 0.0
    min_ppfd_obs = 0.0
    ave_ppfd_obs = 0.0
    std_ppfd_obs = 0.0
    skw_ppfd_obs = 0.0
    krt_ppfd_obs = 0.0
    #
    max_ppfd_ro_h = 0.0
    min_ppfd_ro_h = 0.0
    ave_ppfd_ro_h = 0.0
    std_ppfd_ro_h = 0.0
    skw_ppfd_ro_h = 0.0
    krt_ppfd_ro_h = 0.0
    #
    max_ppfd_ro_l = 0.0
    min_ppfd_ro_l = 0.0
    ave_ppfd_ro_l = 0.0
    std_ppfd_ro_l = 0.0
    skw_ppfd_ro_l = 0.0
    krt_ppfd_ro_l = 0.0
    #
    # NEE Statistics (working, observation backup, & model outlier backups):
    max_nee = 0.0
    min_nee = 0.0
    ave_nee = 0.0
    std_nee = 0.0
    skw_nee = 0.0
    krt_nee = 0.0
    #
    max_nee_obs = 0.0
    min_nee_obs = 0.0
    ave_nee_obs = 0.0
    std_nee_obs = 0.0
    skw_nee_obs = 0.0
    krt_nee_obs = 0.0
    #
    max_nee_ro_h = 0.0
    min_nee_ro_h = 0.0
    ave_nee_ro_h = 0.0
    std_nee_ro_h = 0.0
    skw_nee_ro_h = 0.0
    krt_nee_ro_h = 0.0
    #
    max_nee_ro_l = 0.0
    min_nee_ro_l = 0.0
    ave_nee_ro_l = 0.0
    std_nee_ro_l = 0.0
    skw_nee_ro_l = 0.0
    krt_nee_ro_l = 0.0
    #
    # Correlation coefficients (NEE v PPFD):
    r_obs = 0.0
    r_ro_h = 0.0
    r_ro_l = 0.0
    #
    # Peirce's criterion values:
    #peirce_cap_n = 0.0
    #peirce_lc_n = 1.0
    #peirce_m = 1.0
    peirce_x2 = 0.0
    peirce_delta2 = 0.0
    #
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
    #
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
    #
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
    #
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
    #
    # Model selection:
    mod_select = 0
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, p, n, t, m):
        """
        Name:     FLUX_PARTI.__init__
        Input:    - list, monthly PPFD observations (p)
                  - list, monthly NEE observations (n)
                  - string, flux tower name (t)
                  - datetime.date, current month (m)
        Features: Initialize the class with observation data, calculate basic 
                  statistics, estimate model parameters, and calculate Pearson's
                  correlation coefficient
        """
        # Save tower name & month being processed:
        self.name = t
        self.month = m
        #
        # Save PPFD and NEE lists:
        self.ppfd_obs = p
        self.nee_obs = n
        #
        # Calculate statistics for observations:
        (
            self.max_ppfd, 
            self.min_ppfd, 
            self.ave_ppfd, 
            self.std_ppfd, 
            self.skw_ppfd, 
            self.krt_ppfd) = self.calc_statistics(p)
        (
            self.max_nee, 
            self.min_nee, 
            self.ave_nee, 
            self.std_nee, 
            self.skw_nee, 
            self.krt_nee) = self.calc_statistics(n)
        self.save_stats(obs=1, h=-1)
        #
        # Update model parameters:
        self.update_guess()
        self.save_estimates(obs=1, h=-1)
        #
        # Calculate Pearson's correlation coefficient:
        self.r_obs = self.pearsons_r(n,p)
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def calc_statistics(self, my_array):
        """
        Name:     FLUX_PARTI.calc_statistics
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
            max_val = 0.0
            min_val = 0.0
            ave_val = 0.0
            std_val = 0.0
            skew_val = 0.0
            kurt_val = 0.0
            #
        return (max_val, min_val, ave_val, std_val, skew_val, kurt_val)
    #
    def pearsons_r(self, x, y):
        """
        Name:     FLUX_PARTI.pearsons_r
        Input:    - numpy.ndarray (x)
                  - numpy.ndarray (y)
        Output:   float, Pearson's r (pearsonsr)
        Features: Returns Pearson's correlation between two arrays
        """
        # Make certain data are numpy arrays:
        if not isinstance(x, numpy.ndarray):
            x = numpy.array(x)
        if not isinstance(y, numpy.ndarray):
            y = numpy.array(y)
        #
        # Initialize Pearson's r:
        pearsonsr = -9999.0
        #
        # Make certain both arrays have equal lengths:
        if (len(x) == len(y)):
            try:
                slope, intrcp, pearsonsr, p, sterr = (
                    scipy.stats.linregress(x, y)
                )
                #sxy = ((x - x.mean())*(y - y.mean())).sum()
                #sxx = ((x - x.mean())**2.0).sum()
                #syy = ((y - y.mean())**2.0).sum()
                #pearsonsr = sxy/numpy.sqrt(sxx)/numpy.sqrt(syy)
            except:
                print "Error calculating Pearson's r"
                return -9999.0
            else:
                return pearsonsr
    #
    def save_stats(self, obs, h):
        """
        Name:     FLUX_PARTI.save_stats
        Inputs:   - int, observation flag (obs)
                  - int, hyperbolic model flag (h)
        Output:   None.
        Features: Saves the current statistical parameters to either observation
                  (obs==1) or outliers removed (obs==0) for either the 
                  hyperbolic (h==1) or linear (h==0) model
        """
        if obs == 1:
            # Backup Observation Data:
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
    #
    def update_guess(self):
        """
        Name:     FLUX_PARTI.update_guess
        Input:    None.
        Output:   None.
        Features: Updates the model parameter estimation values based on the 
                  current data statistics
        """
        # Linear model [alpha, R]:
        try:
            self.lm_estimates[0] = (
                0.672*(self.max_nee - self.min_nee)/
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
        #
        # Hyperbolic model [Foo, alpha, R]:
        self.hm_estimates[0] = 3.83*self.std_nee
        try:
            self.hm_estimates[1] = (
                1.96*(self.max_nee - self.min_nee)/
                (self.max_ppfd - self.min_ppfd)
                )
        except ZeroDivisionError:
            self.hm_estimates[1] = 1.96*(self.max_nee - self.min_nee)/1.0e-3
        else:
            if abs(self.hm_estimates[1]) <= 5.0e-4:
                self.hm_estimates[1] = 1.0e-3
        self.hm_estimates[2] = 0.69*self.std_nee
    #
    def save_estimates(self, obs, h):
        """
        Name:     FLUX_PARTI.save_estimates
        Input:    - int, observation flag (obs)
                  - int, hyperbolic model flag (h)
        Output:   None.
        Features: Saves the current model estimates to either observation
                  (obs==1) or outliers removed (obs==0) for either the 
                  hyperbolic (h==1) or linear (h==0) model
        """
        if obs == 1:
            # Backup Model Estimates for Observation Data:
            self.lm_estimates_obs = self.lm_estimates
            self.hm_estimates_obs = self.hm_estimates
        elif obs == 0:
            # Backup Model Estimates for Data with Removed Outliers (RO):
            if h == 1:
                # Outliers Based on Hyperbolic Model
                self.hm_estimates_ro = self.hm_estimates
            elif h == 0:
                # Outliers Based on Linear Model
                self.lm_estimates_ro = self.lm_estimates
    #
    def model_h(self, x, Foo, alpha, R):
        """
        Name:     FLUX_PARTI.model_h
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
        #
        numerator = 1.0*(a*x) + b
        denominator = 1.0*(c*x) + d
        #
        # Compensate for zero division by adding a small number to each value:
        if 0 in denominator:
            denominator = [elem + 1.0e-6 for elem in denominator]
            denominator = numpy.array(denominator)
            #
        return (numerator/denominator)
    #
    def model_l(self, x, alpha, R):
        """
        Name:     FLUX_PARTI.model_l
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
    #
    def calc_model_h(self):
        """
        Name:     FLUX_PARTI.calc_model
        Input:    None.
        Output:   None.
        Features: Calculates NEE and the fitness statistics using the 
                  optimization parameters for the hyperbolic model based on the 
                  observation data
        """
        self.nee_model_h = self.model_h(self.ppfd_obs, 
                                        self.hm_optimized[0], 
                                        self.hm_optimized[1], 
                                        self.hm_optimized[2])
        #
        if -9999 in self.hm_optimized:
            self.hm_mse = -9999.0
            self.hm_rmse = -9999.0
            self.hm_rsqr = -9999.0
        else:
            (self.hm_mse, 
            self.hm_rmse, 
            self.hm_rsqr) = self.goodness_of_fit(self.nee_model_h, 
                                                 self.nee_obs, 
                                                 3)
    #
    def calc_model_h_ro(self):
        """
        Name:     FLUX_PARTI.calc_model_h_ro
        Input:    None.
        Output:   None.
        Features: Calculates NEE and the fitness statistics using the 
                  optimization parameters for the hyperbolic model based on the 
                  observation data with outliers removed
        """
        self.nee_model_h_ro = self.model_h(self.ppfd_obs_h_ro, 
                                           self.hm_optimized_ro[0], 
                                           self.hm_optimized_ro[1], 
                                           self.hm_optimized_ro[2])
        #
        if -9999 in self.hm_optimized_ro:
            self.hm_mse_ro = -9999.0
            self.hm_rmse_ro = -9999.0
            self.hm_rsqr_ro = -9999.0
        else:
            (self.hm_mse_ro, 
            self.hm_rmse_ro, 
            self.hm_rsqr_ro) = self.goodness_of_fit(self.nee_model_h_ro, 
                                                    self.nee_obs_h_ro, 
                                                    3)
    #
    def calc_model_l(self):
        """
        Name:     FLUX_PARTI.calc_model_l
        Input:    None.
        Output:   None.
        Features: Calculates NEE and the fitness statistics using the 
                  optimization parameters for the linear model based on the 
                  observation data
        """
        self.nee_model_l = self.model_l(self.ppfd_obs, 
                                        self.lm_optimized[0], 
                                        self.lm_optimized[1])
        #
        if -9999 in self.lm_optimized:
            self.lm_mse = -9999.0
            self.lm_rmse = -9999.0
            self.lm_rsqr = -9999.0
        else:
            (self.lm_mse, 
            self.lm_rmse, 
            self.lm_rsqr) = self.goodness_of_fit(self.nee_model_l, 
                                                 self.nee_obs, 
                                                 2)
    #
    def calc_model_l_ro(self):
        """
        Name:     FLUX_PARTI.calc_model_l_ro
        Input:    None.
        Output:   None.
        Features: Calculates NEE and the fitness statistics using the 
                  optimization parameters for the linear model based on the 
                  observation data with outliers removed
        """
        self.nee_model_l_ro = self.model_l(self.ppfd_obs_l_ro, 
                                           self.lm_optimized_ro[0], 
                                           self.lm_optimized_ro[1])
        #
        if -9999 in self.lm_optimized_ro:
            self.lm_mse_ro = -9999.0
            self.lm_rmse_ro = -9999.0
            self.lm_rsqr_ro = -9999.0
        else:
            (self.lm_mse_ro, 
            self.lm_rmse_ro, 
            self.lm_rsqr_ro) = self.goodness_of_fit(self.nee_model_l_ro, 
                                                    self.nee_obs_l_ro, 
                                                    2)
    #
    def goodness_of_fit(self, modvals, obsvals, nparams):
        """
        Name:     FLUX_PARTI.goodness_of_fit
        Input:    - numpy.ndarray, modeled values (modvals)
                  - numpy.ndarray, observed values (obsvals)
                  - int, number of model parameters (nparams)
        Output:   tuple, goodness of fit statistics
                  - float, mean squared error (mse)
                  - float, root-mean squared error (rmse)
                  - float, adjusted coefficient of determination (r2_adj)
        Features: Returns the mean squared error, RMSE, and R-squared for  
                  given modeled values
        """
        # Initialize return values:
        mse = 0.0
        rmse = 0.0
        r2_adj = 0.0
        #
        # Check that both have the same number of values and that the length
        # is greater than 4, no divide by zero issues:
        if (len(obsvals) > 4 and len(modvals) == len(obsvals)):
            # Sum of the squared error (SSE):
            sse = sum((obsvals - modvals)**2.0)
            #
            # Mean squared error:
            mse = float(sse)/(len(obsvals) - nparams)
            #
            # Total sum of the squares (SST):
            sst = sum((obsvals - float(sum(obsvals))/len(obsvals))**2.0)
            #
            # R-squared:
            #r2 = 1.0 - float(sse)/sst
            r2_adj = 1.0 - (
                float(sse)/sst*(len(obsvals) - 1.0)/
                float(len(obsvals) - nparams - 1.0)
                )
            #
            # RMSE:
            rmse = numpy.sqrt(float(sse)/len(obsvals))
            #
        return (mse, rmse, r2_adj)
    #
    def remove_mh_outliers(self):
        """
        Name:     FLUX_PARTI.remove_mh_outliers
        Input:    None.
        Output:   int, success flag (rval)
        Features: Returns flag indicating the successful removal of outliers 
                  from the hyperbolic model observations; saves outlier-free
                  data, new statistics, and updated model estimates
        Ref:      Chapter 11.5, GePiSaT Documentation
        """
        # Set Peirce values:
        peirce_cap_n = len(self.nee_obs)
        peirce_lc_n = 1
        peirce_m = 3
        #
        # Calculate tolerance
        self.peirce_x2 = self.peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
        self.peirce_delta2 = self.hm_mse*self.peirce_x2
        #
        # Calculate the square errors:
        sq_errors = (self.nee_obs - self.nee_model_h)**2.0
        #
        # Find if/where exceedance occurs:
        outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
        outliers_found = len(outliers_index)
        #
        # Run check again if no outliers are found in first attempt:
        if (outliers_found == 0):
            peirce_lc_n = 2
            self.peirce_x2 = self.peirce_dev(
                peirce_cap_n, peirce_lc_n, peirce_m
                )
            self.peirce_delta2 = self.hm_mse*self.peirce_x2
            outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
            outliers_found = len(outliers_index)
            # Reset n
            peirce_lc_n = 1
            #
        # Increment n until it is greater than number of outliers found:
        while (peirce_lc_n <= outliers_found):
            peirce_lc_n += 1
            # Check that n < N:
            if peirce_lc_n >= peirce_cap_n:
                peirce_lc_n = outliers_found + 1.0
            else:
                self.peirce_x2 = self.peirce_dev(
                    peirce_cap_n, peirce_lc_n, peirce_m)
                self.peirce_delta2 = self.hm_mse*self.peirce_x2
                outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
                outliers_found = len(outliers_index)
            #
        self.hm_outliers = outliers_found
        #
        # Remove outliers found:
        self.ppfd_obs_h_ro = numpy.delete(self.ppfd_obs, outliers_index)
        self.nee_obs_h_ro = numpy.delete(self.nee_obs, outliers_index)
        #
        # Check that both arrays still have observations:
        if self.ppfd_obs_h_ro.any() and self.nee_obs_h_ro.any():
            # Recalculate statistics and update guess values:
            (
                self.max_ppfd, 
                self.min_ppfd, 
                self.ave_ppfd, 
                self.std_ppfd, 
                self.skw_ppfd, 
                self.krt_ppfd) = self.calc_statistics(self.ppfd_obs_h_ro)
            (
                self.max_nee, 
                self.min_nee, 
                self.ave_nee, 
                self.std_nee, 
                self.skw_nee, 
                self.krt_nee) = self.calc_statistics(self.nee_obs_h_ro)
            self.save_stats(obs=0, h=1)
            self.update_guess()
            self.save_estimates(obs=0, h=1)
            self.r_ro_h = self.pearsons_r(
                self.ppfd_obs_h_ro, 
                self.nee_obs_h_ro
                )
            rval = 1
        else:
            rval = 0
        return rval
    #
    def remove_ml_outliers(self):
        """
        Name:     FLUX_PARTI.remove_ml_outliers
        Input:    None.
        Output:   int, success flag (rval)
        Features: Returns flag indicating the successful removal of outliers 
                  from the linear model observations; saves outlier-free
                  data, new statistics, and updated model estimates
        Ref:      Chapter 11.5, GePiSaT Documentation
        """
        # Set Peirce values:
        peirce_cap_n = len(self.nee_obs)
        peirce_lc_n = 1
        peirce_m = 2
        #
        # Calculate tolerance
        self.peirce_x2 = self.peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
        self.peirce_delta2 = self.hm_mse*self.peirce_x2
        #
        # Calculate the square errors:
        sq_errors = (self.nee_obs - self.nee_model_l)**2.0
        #
        # Find if/where exceedance occurs:
        outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
        outliers_found = len(outliers_index)
        #
        # Run check again if no outliers are found in first attempt:
        if (outliers_found == 0):
            peirce_lc_n = 2
            self.peirce_x2 = self.peirce_dev(
                peirce_cap_n, peirce_lc_n, peirce_m
                )
            self.peirce_delta2 = self.hm_mse*self.peirce_x2
            outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
            outliers_found = len(outliers_index)
            # Reset n
            peirce_lc_n = 1
            #
        # Increment n until it is greater than number of outliers found:
        while (peirce_lc_n <= outliers_found):
            peirce_lc_n += 1
            # Check that n < N:
            if peirce_lc_n >= peirce_cap_n:
                peirce_lc_n = outliers_found + 1.0
            else:
	        self.peirce_x2 = self.peirce_dev(
	           peirce_cap_n, peirce_lc_n, peirce_m
	           )
                self.peirce_delta2 = self.hm_mse*self.peirce_x2
                outliers_index = numpy.where(sq_errors > self.peirce_delta2)[0]
                outliers_found = len(outliers_index)
            #
        self.lm_outliers = outliers_found
        #
        # Remove outliers found:
        self.ppfd_obs_l_ro = numpy.delete(self.ppfd_obs, outliers_index)
        self.nee_obs_l_ro = numpy.delete(self.nee_obs, outliers_index)
        #
        # Make sure that arrays still have observations:
        if self.ppfd_obs_l_ro.any() and self.nee_obs_l_ro.any():
            # Recalculate statistics and update guess values:
            (
                self.max_ppfd, 
                self.min_ppfd, 
                self.ave_ppfd, 
                self.std_ppfd, 
                self.skw_ppfd, 
                self.krt_ppfd) = self.calc_statistics(self.ppfd_obs_h_ro)
            (
                self.max_nee, 
                self.min_nee, 
                self.ave_nee, 
                self.std_nee, 
                self.skw_nee, 
                self.krt_nee) = self.calc_statistics(self.nee_obs_h_ro)
            self.save_stats(obs=0, h=0)
            self.update_guess()
            self.save_estimates(obs=0, h=0)
            self.r_ro_l = self.pearsons_r(
                self.ppfd_obs_l_ro,
                self.nee_obs_l_ro
                )
            rval = 1
        else:
            rval = 0
        return rval
    #
    def peirce_dev(self, peirce_cap_n, peirce_lc_n, peirce_m):
        """
        Name:     FLUX_PARTI.peirce_dev
        Input:    - int, total number of observations (peirce_cap_n)
                  - int, number of outliers to be removed (peirce_lc_n)
                  - int, number of model unknowns (peirce_m)
        Output:   float, squared error threshold (x2)
        Features: Returns the squared threshold error deviation for outlier 
                  identification using Peirce's criterion based on Gould's
                  methodology
        Ref:      Chapter 11.5, GePiSaT Documentation
        """
        # Assign floats to input variables:
        N = float(peirce_cap_n)
        n = float(peirce_lc_n)
        m = float(peirce_m)
        #
        # Check the total number of observations:
        if N > 1:
            # Calculate Q (Nth root of Gould's equation B):
            # Note: 1/N exponent is factored to each individual term to prevent
            # OverflowError with large N (e.g., >142)
            Q = (n**(n/N)*(N-n)**((N-n)/N))/N
            #
            # Initialize R values (as floats):
            Rnew = 1.0  # <- Tried values between 1 and 10 and all seem stable
            Rold = 0.0  # <- Necessary to prompt while loop
            #
            while ( abs(Rnew - Rold) > (N*2.0e-16) ):
                # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
                ldiv = Rnew**n
                if ldiv == 0:
                    ldiv = 1.0e-6
                Lamda = ((Q**N)/(ldiv))**(1.0/(N - n))
                #
                # Calculate x-squared (straight-forward Gould's equation C):
                x2 = 1.0 + (N - m - n)/n*(1.0 - Lamda**2.0)
                #
                # Return 0 for negative x2 values:
                if x2 < 0:
                    x2 = 0
                    Rold = Rnew
                else:
                    # Use x-squared to update R (Gould's equation D):
                    Rold = Rnew
                    Rnew = (
                        numpy.exp((x2 - 1)/2.0)*
                        scipy.special.erfc(numpy.sqrt(x2)/numpy.sqrt(2.0))
                        )
                    #
        else:
            x2 = 0.0
        return x2
    #
    def summary_statistics(self):
        """
        Name:     FLUX_PARTI.summary_statistics
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
        # 08. foo_err_obs_h   :: Foo standard error for observations for model H
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
                self.name, self.month,                                 # 01--02
                len(self.nee_obs), self.hm_outliers, self.lm_outliers, # 03--05
                self.hm_estimates_obs[0], self.hm_optimized[0],        # 06--07
                self.hm_optim_err[0],     self.hm_optim_t[0],          # 08--09
                self.hm_optim_p[0],       self.hm_estimates_ro[0],     # 10--11
                self.hm_optimized_ro[0],  self.hm_optim_err_ro[0],     # 12--13
                self.hm_optim_t_ro[0],    self.hm_optim_p_ro[0],       # 14--15
                self.hm_estimates_obs[1], self.hm_optimized[1],        # 16--17
                self.hm_optim_err[1],     self.hm_optim_t[1],          # 18--19
                self.hm_optim_p[1],       self.hm_estimates_ro[1],     # 20--21
                self.hm_optimized_ro[1],  self.hm_optim_err_ro[1],     # 22--23
                self.hm_optim_t_ro[1],    self.hm_optim_p_ro[1],       # 24--25
                self.lm_estimates_obs[0], self.lm_optimized[0],        # 26--27
                self.lm_optim_err[0],     self.lm_optim_t[0],          # 28--29
                self.lm_optim_p[0],       self.lm_estimates_ro[0],     # 30--31
                self.lm_optimized_ro[0],  self.lm_optim_err_ro[0],     # 32--33
                self.lm_optim_t_ro[0],    self.lm_optim_p_ro[0],       # 34--35
                self.hm_estimates_obs[2], self.hm_optimized[2],        # 36--37
                self.hm_optim_err[2],     self.hm_optim_t[2],          # 38--39
                self.hm_optim_p[2],       self.hm_estimates_ro[2],     # 40--41
                self.hm_optimized_ro[2],  self.hm_optim_err_ro[2],     # 42--43
                self.hm_optim_t_ro[2],    self.hm_optim_p_ro[2],       # 44--45
                self.lm_estimates_obs[1], self.lm_optimized[1],        # 46--47
                self.lm_optim_err[1],     self.lm_optim_t[1],          # 48--49
                self.lm_optim_p[1],       self.lm_estimates_ro[1],     # 50--51
                self.lm_optimized_ro[1],  self.lm_optim_err_ro[1],     # 52--53
                self.lm_optim_t_ro[1],    self.lm_optim_p_ro[1],       # 54--55
                self.hm_rsqr,             self.hm_rsqr_ro,             # 56--57
                self.hm_rmse,             self.hm_rmse_ro,             # 58--59
                self.lm_rsqr,             self.lm_rsqr_ro,             # 60--61
                self.lm_rmse,             self.lm_rmse_ro,             # 62--63
                self.min_ppfd_obs,        self.max_ppfd_obs,           # 64--65
                self.ave_ppfd_obs,        self.std_ppfd_obs,           # 66--67
                self.skw_ppfd_obs,        self.krt_ppfd_obs,           # 68--69
                self.min_ppfd_ro_h,       self.max_ppfd_ro_h,          # 70--71
                self.ave_ppfd_ro_h,       self.std_ppfd_ro_h,          # 72--73
                self.skw_ppfd_ro_h,       self.krt_ppfd_ro_h,          # 74--75
                self.min_ppfd_ro_l,       self.max_ppfd_ro_l,          # 76--77
                self.ave_ppfd_ro_l,       self.std_ppfd_ro_l,          # 78--79
                self.skw_ppfd_ro_l,       self.krt_ppfd_ro_l,          # 80--81
                self.min_nee_obs,         self.max_nee_obs,            # 82--83
                self.ave_nee_obs,         self.std_nee_obs,            # 84--85
                self.skw_nee_obs,         self.krt_nee_obs,            # 86--87
                self.min_nee_ro_h,        self.max_nee_ro_h,           # 88--89
                self.ave_nee_ro_h,        self.std_nee_ro_h,           # 90--91
                self.skw_nee_ro_h,        self.krt_nee_ro_h,           # 92--93
                self.min_nee_ro_l,        self.max_nee_ro_l,           # 94--95
                self.ave_nee_ro_l,        self.std_nee_ro_l,           # 96--97
                self.skw_nee_ro_l,        self.krt_nee_ro_l,           # 98--99
                self.r_obs,               self.r_ro_h,               # 100--101
                self.r_ro_l,              self.mod_select            # 102--103
                )
        #
        return sum_stat
    #
    def model_selection(self):
        """
        Name:     FLUX_PARTI.model_selection
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
        #
        # Initialize selection dictionary:
        select_dict = {
            'obs_h' : 1,
            'ro_h'  : 2,
            'obs_l' : 3,
            'ro_l'  : 4
            }
        #
        # Initialize model failure dictionary:
        model_dict = {
            'obs_h' : 0,
            'ro_h'  : 0,
            'obs_l' : 0,
            'ro_l'  : 0
            }
        #
        # Initialize R-squared dictionary:
        r2_dict = {
            'obs_h' : self.hm_rsqr,
            'obs_l' : self.lm_rsqr,
            'ro_h'  : self.hm_rsqr_ro,
            'ro_l'  : self.lm_rsqr_ro
            }
        #
        # Initialize P-value dictionary:
        p_dict = {
            'obs_h' : numpy.array(self.hm_optim_p),
            'ro_h'  : numpy.array(self.hm_optim_p_ro),
            'obs_l' : numpy.array(self.lm_optim_p),
            'ro_l'  : numpy.array(self.lm_optim_p_ro)
        }
        #
        # Initialize parameter value dictionary:
        value_dict = {
            'obs_h' : {'foo'   : self.hm_optimized[0],
                       'alpha' : self.hm_optimized[1],
                       'r'     : self.hm_optimized[2]
                       },
            'ro_h'  : {'foo'   : self.hm_optimized_ro[0],
                       'alpha' : self.hm_optimized_ro[1],
                       'r'     : self.hm_optimized_ro[2]
                       },
            'obs_l' : {'alpha' : self.lm_optimized[0],
                       'r'     : self.lm_optimized[1]
                       },
            'ro_l'  : {'alpha' : self.lm_optimized_ro[0],
                       'r'     : self.lm_optimized_ro[1]
                       }
        }
        # Initialize model selection:
        best_model = 0
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # FIRST CRITERIA--- R-SQUARED THRESHOLD
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (numpy.array(r2_dict.values()) < r2_thresh).any():
            # Eliminate those with poor R2:
            m_failed = [k for k,v in r2_dict.items() if v < r2_thresh]
            for m in m_failed:
                model_dict[m] = 1
        #
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
        #
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
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ANALYSIS:
        # Let's see which models passed the selection criteria and how many
        # of them there are
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        mods = [k for k,v in model_dict.items() if v == 0]
        howmany = len(mods)
        #
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
                #
                # Get the max R-squared and min p-value:
                max_rs = max(mods_rs)
                min_ps = min(mods_ps)
                #
                # How many models meet this criteria:
                mods_temp = [
                    m for m in mods if 
                    mods_rs[mods.index(m)] == max_rs and
                    mods_ps[mods.index(m)] == min_ps 
                    ]
                n_temp = len(mods_temp)
                #
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
                    if 1.01*r2_dict[mods_temp[1]] >= 0.99*r2_dict[mods_temp[0]]:
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
        #
        # Save model selection:
        self.mod_select = best_model
    #
    def calc_gpp(self, ppfd):
        """
        Name:     FLUX_PARTI.calc_gpp
        Input:    numpy.ndarray, monthly PPFD (ppfd)
        Output:   tuple, modeled GPP and associated error
                  - numpy.ndarray, modeled GPP (gpp)
                  - numpy.ndarray, associated error (gpp_err)
        Features: Returns GPP (umol m-2 s-1) and associated error based on the 
                  best model
        Ref:      Chapters 11.2 & 11.3, GePiSaT Documentation
        """
        # Get model selection:
        if (self.mod_select == 1 or self.mod_select == 2):
            model = "H"
        elif (self.mod_select == 3 or self.mod_select == 4):
            model = "L"
        #
        # Get outlier selection:
        if (self.mod_select == 1 or self.mod_select == 3):
            outlier = 0
        elif (self.mod_select == 2 or self.mod_select == 4):
            outlier = 1
        #
        if self.mod_select == 0:
            # Default to hyperbolic model with outliers removed estimates
            (foo, alpha, r) = self.hm_estimates_ro
            (foo_err, alpha_err, r_err) = (0., 0., 0.)
            #
            if foo == 0 or alpha == 0:
                # Use linear estimates instead:
                (alpha, r) = self.lm_estimates_ro
                (alpha_err, r_err) = (0., 0.)
                #
                # Calculate GPP and GPP err:
                gpp = (alpha*ppfd)
                gpp_err = (ppfd*alpha_err)
                #
            else:
                # Variable substitutes:
                apf = alpha*ppfd + foo
                afp = alpha*foo*ppfd
                afpp = alpha*foo*(ppfd**2)
                #
                # Calculate GPP and its associated error:
                gpp = afp/apf
                gpp_err = numpy.sqrt(
                    (alpha_err**2)*((foo*ppfd*apf - afpp)/(apf**2))**2 + 
                    (foo_err**2)*((alpha*ppfd*apf - afp)/(apf**2))**2
                )
            #
            # Clip out negative values of GPP:
            gpp = gpp.clip(min=0)
            #
        elif model.upper() == "H":
            # Retrieve parameters for hyperbolic model:
            if outlier == 0:
                (foo, alpha, r) = self.hm_optimized
                (foo_err, alpha_err, r_err) = self.hm_optim_err
            elif outlier == 1:
                (foo, alpha, r) = self.hm_optimized_ro
                (foo_err, alpha_err, r_err) = self.hm_optim_err_ro
                #
            # Variable substitutes:
            apf = alpha*ppfd + foo
            afp = alpha*foo*ppfd
            afpp = alpha*foo*(ppfd**2)
            #
            # Calculate GPP and its associated error:
            gpp = afp/apf
            gpp_err = numpy.sqrt(
                (alpha_err**2)*((foo*ppfd*apf - afpp)/(apf**2))**2 + 
                (foo_err**2)*((alpha*ppfd*apf - afp)/(apf**2))**2
            )
            #
        elif model.upper() == "L":
            # Retrieve parameters for linear model:
            if outlier == 0:
                (alpha, r) = self.lm_optimized
                (alpha_err, r_err) = self.lm_optim_err
            elif outlier == 1:
                (alpha, r) = self.lm_optimized_ro
                (alpha_err, r_err) = self.lm_optim_err_ro
                #
            # Calculate GPP:
            gpp = (alpha*ppfd)
            gpp_err = (ppfd*alpha_err)
            #
        # Return GPP
        return (gpp, gpp_err)
#
class SOLAR:
    """
    Name:     SOLAR
    Features: This class calculates the half-hourly extraterrestrial PPFD 
              [umol m-2 s-1], and the daily solar irradiation Ho [J m-2]
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
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization 
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self, lon, lat, n, y=0):
        """
        Name:     SOLAR.__init__
        Input:    - float, longitude (lon)
                  - float, latitude (lat)
                  - int, day of year (n)
                  - int, year (y)
        """
        # Error handle and assign required public variables:
        self.year = y
        if lat > 90.0 or lat < -90.0:
            print "Latitude outside range of validity (-90 to 90)!"
            exit(1)
        else:
            self.lat = lat
        if lon > 180.0 or lon < -180.0:
            print "Longitude outside range of validity (-180 to 180)!"
            exit(1)
        else:
            self.lon = lon
        if n < 1 or n > 366:
            print "Day outside range of validity (1 to 366)!"
            exit(1)
        else:
            self.doy = n
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 0. Create local time series, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        local_hh = numpy.array([0.5*i for i in xrange(48)])
        #self.local_hh = local_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 1. Calculate number of days in year (kN), days
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if y == 0:
            self.kN = 365.
        else:
            self.kN = (
                self.julian_day((y + 1), 1, 1) - 
                self.julian_day(y, 1, 1)
            )
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 2. Calculate heliocentric longitudes (nu and lambda), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger (1978)
        my_nu, my_lambda = self.berger_tls(n)
        #self.nu_deg = my_nu
        #self.lambda_deg = my_lambda
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 3. Calculate distance factor (dr), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Berger et al. (1993)
        my_rho = (1.0 - self.ke**2)/(1.0 + self.ke*self.dcos(my_nu))
        dr = (1.0/my_rho)**2
        #self.dr = dr
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 4. Calculate declination angle (delta), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Woolf (1968)
        delta = numpy.arcsin(self.dsin(my_lambda)*self.dsin(self.keps))
        delta = delta*(180.0/numpy.pi)
        #self.delta_deg = delta
        #
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
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 6. Calculate the equation of time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Spencer (1971)
        eot = self.spencer_eot(n)
        #self.eot_hour = eot
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 7. Calculate the longitude correction, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        lc = (15.0*tz_hour - lon)/15.0
        #self.lc_hour = lc
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 8. Calculate the solar time, hours
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ts_hh = local_hh + eot - lc
        #self.ts_hh = ts_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 9. Calculate the hour angle, degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        w_hh = (360./24.)*(ts_hh - 12.0)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 10. Calculate variable substitutes (u and v), unitless
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ru = self.dsin(delta)*self.dsin(lat)
        rv = self.dcos(delta)*self.dcos(lat)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 11. Calculate the sunset hour angle (hs), degrees
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Note: ru/rv == tan(delta)*tan(lat)
        # Eq. 3.22, Stine & Geyer (2001)
        if (ru/rv) >= 1.0:
            hs = 180.0   # Polar day (no sunset)
        elif (ru/rv) <= -1.0:
            hs = 0.0     # Polar night (no sunrise)
        else:
            hs = (180.0/numpy.pi)*numpy.arccos(-1.0*ru/rv)
        #self.hs_deg = hs
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 12. Calculate the solar radiation flux, W/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        io_hh = self.kGsc*dr*(ru + rv*self.dcos(w_hh))
        io_hh = io_hh.clip(min=0)
        #self.io_wm2 = io_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 13. Calculate the half-hourly PPFD, umol/m^2/s
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ppfd_hh = self.kfFEC*io_hh
        self.ppfd_hh = ppfd_hh
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # 14. Calculate the daily solar irradiation, J/m^2
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ho = (86400.0/numpy.pi)*self.kGsc*dr*(
            ru*(numpy.pi/180.0)*hs + rv*self.dsin(hs)
        )
        self.ho_jm2 = ho
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def berger_tls(self, n):
        """
        Name:     SOLAR.berger_tls
        Input:    int, day of year
        Output:   tuple, 
                  - true anomaly, degrees
                  - true longitude, degrees
        Features: Returns true anomaly and true longitude for a given day
        Depends:  - ke
                  - komega
                  - kN
                  - dsin()
        Ref:      Berger, A. L. (1978), Long term variations of daily insolation
                  and quaternary climatic changes, J. Atmos. Sci., 35, 2362-
                  2367.
        """
        # Variable substitutes:
        xee = self.ke**2 
        xec = self.ke**3
        xse = numpy.sqrt(1.0 - xee)
        #
        # Mean longitude for vernal equinox:
        xlam = (
            (self.ke/2.0 + xec/8.0)*(1.0 + xse)*self.dsin(self.komega) - 
            xee/4.0*(0.5 + xse)*self.dsin(2.0*self.komega) + 
            xec/8.0*(1.0/3.0 + xse)*self.dsin(3.0*self.komega)
            )
        xlam = numpy.degrees(2.0*xlam)
        #
        # Mean longitude for day of year:
        dlamm = xlam + (n - 80.0)*(360.0/self.kN)
        #
        # Mean anomaly:
        anm = dlamm - self.komega
        ranm = numpy.radians(anm)
        #
        # True anomaly:
        ranv = (ranm + (2.0*self.ke - xec/4.0)*numpy.sin(ranm) + 
            5.0/4.0*xee*numpy.sin(2.0*ranm) + 
            13.0/12.0*xec*numpy.sin(3.0*ranm))
        anv = numpy.degrees(ranv)
        #
        # True longitude:
        my_tls = anv + self.komega
        if my_tls < 0:
            my_tls += 360.0
        elif my_tls > 360:
            my_tls -= 360.0
        #
        # True anomaly:
        my_nu = (my_tls - self.komega)
        if my_nu < 0:
            my_nu += 360.0
        #
        return(my_nu, my_tls)
    #
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
    def spencer_eot(self, n):
        """
        Name:     SOLAR.spencer_eot
        Input:    int, day of the year (n)
        Output:   float, equation of time, hours
        Features: Returns the equation of time
        Ref:      Spencer, J.W. (1971), Fourier series representation of the 
                  position of the sun, Search, 2 (5), p. 172. 
        """
        B = 2.0*numpy.pi*(n - 1.0)/self.kN
        my_eot = 12.0/(numpy.pi)*(
            (7.5e-6) + (1.868e-3)*self.dcos(B) - (3.2077e-2)*self.dsin(B) - 
            (1.4615e-2)*self.dcos(2.0*B) - (4.0849e-2)*self.dsin(2.0*B)
        )
        return(my_eot)
#
class LUE:
    """
    Name:     LUE
    Features: This class stores monthly LUE estimates and writes results to
              file.
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Variable Definitions
    # ////////////////////////////////////////////////////////////////////////
    # Dictionary of station's monthly GPP, PPFD, fAPAR vals:
    station_vals = {}
    #
    # Dictionary of station's LUE:
    station_lue = {}
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     LUE.__init__
        Input:    None.
        Features: Initializes empty dictionaries for LUE
        """
        self.station_vals = {}
        self.station_lue = {}
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def add_station_val(self, station, month, gpp, gpp_err, 
                        fpar, ppfd, vpd, alpha, tmp, co2, patm):
        """
        Name:     LUE.add_station_val
        Input:    - string, station name (station)
                  - datetime.date, current month (month)
                  - float, monthly GPP (gpp)
                  - float, associated GPP error (gpp_err)
                  - float, monthly FAPAR (fpar)
                  - float, monthly PPFD (ppfd)
                  - float, monthly VPD (vpd)
                  - float, monthly CPA (alpha)
                  - float, monthly air temp (tmp)
                  - float, annual atm. CO2 (co2)
                  - float, monthly atm. pressure (patm)
        Output:   None.
        Features: Appends a set of monthly values to the value dictionary
        """
        # Check for missing alpha (due to new STASH code):
        if alpha is None:
            alpha = -9999.0
        #
        # Place parameters into a tuple:
        params = (month, gpp, gpp_err, fpar, ppfd, vpd, alpha, tmp, co2, patm)
        #
        # Initialize list if station key doesn't exist:
        if station not in self.station_vals.keys():
            self.station_vals[station] = []
        #
        # Add new parameters to list:
        self.station_vals[station].append(params)
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
            lue_head = (
                "Timestamp,GPP.mol_m2,GPP_err,fAPAR,PPFD.mol_m2,VPD.Pa,"
                "ALPHA,Tc.deg_C,CO2.ppm,Patm.Pa\n"
                )
            try:
                f = open(out_file, 'w')
            except IOError:
                print "Cannot write to file:", out_file
            else:
                f.write(lue_head)
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

################################################################################
## FUNCTIONS 
################################################################################
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

def basic_lue(x, a):
    """
    Name:     basic_lue
    Input:    - numpy.ndarray (x)
              - float, LUE parameter (a)
    Output:   numpy.ndarray, modeled GPP
    Features: Returns an array after applying the basic light-use efficiency 
              model, GPP = LUE * (fPAR * PPFD)
    """
    return x*a

def calc_gstar(tc):
    """
    Name:     calc_gstar
    Input:    float, air temperature, degrees C (tc)
    Output:   float, (gs)
    Features: Returns the temperature-dependent photorespiratory compensation
              point (Pascals) based on constants derived from Bernacchi et al.
              (2001) study.
    Ref:      Bernacchi et al. (2001)
    """
    # Define constants
    gs25 = 4.220  # Pa, assuming 25 deg C & 98.716 kPa)
    dha = 37830   # J/mol
    kR = 8.3145   # J/mol/K
    #
    gs = gs25*numpy.exp(dha*(tc - 25.0)/(298.15*kR*(tc + 273.15)))
    return gs

def calc_k(tc, patm):
    """
    Name:     calc_k
    Input:    - float, air temperature, degrees C (tc)
              - float, atmospheric pressure, Pa (patm)
    Output:   float, Michaelis-Menton K, Pa (mmk)
    Features: Returns the temperature and pressure dependent Michaelis-Menten
              coefficient (Pascals) based on constants derived from Bernacchi 
              et al. (2001) study.
    Ref:      Bernacchi et al. (2001)
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

def connect_sql():
    """
    Name:     connectSQL
    Input:    None.
    Output:   psycopg2 connection (con)
    Features: Connects to postgreSQL database and returns connection handle
    """
    # Open credentials file:
    cred_file = "./user.txt"
    if os.path.isfile(cred_file):
        f = open(cred_file, "r")
        cred_data = f.readline()
        if cred_data:
            cred_data = cred_data.rstrip()
            my_user, my_pass = cred_data.split(',')
    else:
        my_user, my_pass = ('postgres', 'password1')
    #
    # Initialize connection variable:
    con = None
    #
    # Test database connection:
    try:
        con = psycopg2.connect(
            database='gepisat', 
            #database='flux_test',
            user=my_user, 
            host='localhost', 
            password=my_pass
            )
    #
    except psycopg2.DatabaseError, e:
        print 'Error %s' % e
        exit(1)
    #
    return con

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
    #
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
        #
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
        #
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
        #
        # Determine centroid based on min distance:
        if dist_A == min_dist:
            my_centroid = (bb_lon_max, bb_lat_max)
        elif dist_B == min_dist:
            my_centroid = (bb_lon_max, bb_lat_min)
        elif dist_C == min_dist:
            my_centroid = (bb_lon_min, bb_lat_max)
        elif dist_D == min_dist:
            my_centroid = (bb_lon_min, bb_lat_min)
            #
    # Return nearest centroid:
    return my_centroid

def next_gen_lue(x, phi_o, beta):
    """
    Name:     next_gen_lue
    Input:    - numpy.ndarray, (k,M) shape of predictors (x)
                > 'Iabs' : mol/m2, fAPARxPPFD
                > 'ca' : Pa, atmospheric CO2
                > 'Gs' : Pa, photores. comp. point
                > 'D' : Pa, vapor pressure deficit
                > 'K' : Pa, Michaelis-Menten coeff.
                > 'eta' : mPa s, viscosity of water
                > 'fa' : unitless, function of alpha
              - float, intrinsic quantum efficiency (phi_o)
              - float, beta parameter (beta)
    Output:   numpy.ndarray (gpp)
    Features: Returns array of GPP based on the next-generation light and
              water use efficiency model.
    """
    # Define constants:
    kc = 0.41    # Jmax cost parameter
    #
    # Check beta parameter (divide by zero):
    if beta == 0:
        beta += 1e-6
    #
    # Define variable substitutes:
    vdg = x['ca'] - x['Gs']
    vag = x['ca'] + 2.0*x['Gs']
    vsr = numpy.sqrt(1.6*x['eta']*x['D']/(beta*(x['K'] + x['Gs'])))
    #
    # Based on the m' formulation (see Regressing_LUE.pdf)
    gpp = phi_o*x['Iabs']*x['fa']*numpy.sqrt(
        (vdg/(vag + 3.*x['Gs']*vsr))**2 - 
        kc**(2./3.)*(vdg/(vag + 3.*x['Gs']*vsr))**(4./3.)
    )
    return gpp

def predict_params(ca, d, eta, gpp, gs, iabs, k):
    """
    Name:     predict_params
    Input:    -numpy.ndarray, atmos. CO2 conc. (ca)
              -numpy.ndarray, vap. press. deficit (d)
              -numpy.ndarray, water visc. (eta)
              -numpy.ndarray, GPP (gpp)
              -numpy.ndarray, photo. resp. comp. point (gs)
              -numpy.ndarray, abs. PAR (iabs)
              -numpy.ndarray, Michaelis-Menten coef. (k)
    Output:   tuple
    Features: Returns statistically-based estimate of phio and theoretical
              estimate of beta
    """
    # Predict phio (statistical relationship)
    phio_p = ((1.761e-1)*gs['min'] + (1.702e-2)*gpp['std'] -
              (6.894e-3)*k['min'] - (3.531e-4)*iabs['ave'])
    #
    # Predict beta (based on Colin's method)
    my_d = 0.5*(d['max'] + d['min'])
    my_k = k['ave']
    my_g = 0.5*(gs['max'] + gs['min'])
    my_n = 0.5*(eta['max'] + eta['min'])
    #
    # Estimate chi based on VPD (linear interpolation):
    chi = 0.5 - (0.5 - 0.9)*(2.5e3 - my_d)/(2.5e3)
    beta_p = 1.6*my_d/(my_k + my_g)
    beta_p *= ((chi - my_g)/(1.0 - chi))**2
    beta_p *= my_n
    beta_p *= 0.0017154
    #
    return (phio_p[0], beta_p[0])

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
    #
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
        )  % (
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
        print "Cannot open file", summary_file, "for writing"
    else:
        SFILE.write(summary_header)
        SFILE.close()

def viscosity_h2o(tc):
    """
    Name:     viscosity_h2o
    Input:    float, air temperature, degrees C (tc)
    Output:   float, viscosity of water, mPa s (n)
    Features: Returns the temperature-dependent viscosity of water (mPa s) 
              based on the Vogel Equation
    """
    n = 2.4263e-2*numpy.exp(5.78919e2/((tc + 273.15) - 1.37546e2))
    return n

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         Dependant Functions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calc_lue(lue_class, station):
    """
    Name:     calc_lue
    Input:    - LUE class (lue_class)
              - str, station name (station)
    Output:   None.
    Features: Fits the next generation LUE model to the monthly flux tower 
              data and saves the fit to LUE class
    Depends:  - next_gen_lue
              - calc_gstar
              - calc_k
              - viscosity_h2o
              - predict_params
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
    temp_stats = numpy.array(tuple([-9999., -9999., -9999., 
                                    -9999., -9999., -9999.]),
                             dtype={'names' : ('max', 'min', 'ave', 
                                               'std', 'skw', 'krt'),
                                    'formats' : ('f4', 'f4', 'f4', 
                                                 'f4', 'f4', 'f4')},
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
        print y_data
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
            #print "fit_opt:", fit_opt
            #print "fit_cov:", fit_cov
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
    # Save fit to LUE class
    params = (tuple([st_rsqr, est_phio, st_phio, est_beta, st_beta, 
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

def flux_to_grid(flux_station):
    """
    Name:     flux_to_grid
    Input:    str, station name (flux_station)
    Output:   int, grid station ID (grid_station)
    Features: Returns grid station ID based on the location of a given flux 
              tower
    Depends:  - connect_sql
              - get_lon_lat
              - grid_centroid
    """
    # Get lat and lon of flux tower:
    (fst_lon, fst_lat) = get_lon_lat(flux_station)
    #
    # Determine grid centroid lon and lat:
    (grd_lon, grd_lat) = grid_centroid(fst_lon, fst_lat)
    #
    # Get grid station name based on centroid coordinates:
    params = ("grid", grd_lon, grd_lat)
    q = (
        "SELECT met_data.stationid "
        "FROM met_data "
        "WHERE met_data.geom = %s "
        "AND met_data.lon = %s "
        "AND met_data.lat = %s;"
        )
    #
    # Connect to database:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return results:
    cur.execute(q, params)
    grid_station = cur.fetchone()[0]
    con.close()
    return grid_station

def gapfill_ppfd(station, start_date, to_write):
    """
    Name:     gapfill_ppfd
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
              - flux_to_grid
              - get_daily_ppfd
              - get_data_point
              - get_lon_lat
              - get_msvidx
              - SOLAR
    """
    # Initialize monthly gapless PPFD array & associated timestamps:
    monthly_ppfd_gapless = numpy.array([])
    monthly_timestamp_hh = numpy.array([])
    #
    # Get flux station lon/lat coordinates:
    (flux_lon, flux_lat) = get_lon_lat(station)
    #
    # Get grid station ID, lon/lat coordinates, and msvidx for 'SWdown':
    grid_station = flux_to_grid(station)
    (grid_lon, grid_lat) = get_lon_lat(grid_station)
    grid_msvidx = get_msvidx(grid_station, 'SWdown')
    #
    # Calculate the end date:
    end_date = add_one_month(start_date)
    #
    # Create output file (if to_write):
    out_file = "out/%s-GF_%s.txt" % (station, start_date)
    if to_write:
        header = "Timestamp, PPFDobs, PPFDgf\n"
        try:
            f = open(out_file, "w")
        except IOError:
            print "Cannot write to file:", out_file
        else:
            f.write(header)
            f.close()
    #
    # Iterate through each day of the month:
    cur_date = start_date
    while cur_date < end_date:
        # Reset dictionary for gap-filled time series:
        gapfill_dict = {}
        #
        # Create start and end datetime objects for this day:
        start_time = datetime.datetime(cur_date.year, 
                                       cur_date.month, 
                                       cur_date.day, 
                                       0, 0, 0)
        end_time = datetime.datetime(cur_date.year, 
                                     cur_date.month, 
                                     cur_date.day, 
                                     23, 59, 59)
        #
        # Initialize dictionary values with half-hourly timestamp keys:
        while start_time < end_time:
            my_time = "%s" % start_time.time()
            gapfill_dict[my_time] = 0.0
            start_time = start_time + datetime.timedelta(minutes=30)
            #
            # Add datetime object to array of monthly timestamps:
            monthly_timestamp_hh = numpy.append(monthly_timestamp_hh, 
                                                [start_time,])
            #
        # Get daily PPFD values (in a numpy.array) from database:
        #   NOTE: last index is start of next day
        (daily_ts, daily_ppfd) = get_daily_ppfd(station, cur_date)
        #
        # Check to see if any gaps are present in current day's observations:
        number_obs = len(daily_ppfd)
        if number_obs < 49:
            # Gaps are present:
            #
            # Convert date to Julian day:
            jday = cur_date.timetuple().tm_yday
            #
            # Calculate daily ET solar radiation curve:
            et_solar = SOLAR(flux_lon, flux_lat, jday, cur_date.year)
            #
            # Get satellite measurement of solar rad (SWdown) [W m-2]:
            grid_srad = get_data_point(grid_msvidx, cur_date)
            #
            # Convert to daily shortwave radiation [J m-2]:
            grid_srad_d = (86400.0)*grid_srad
            #
            # Calculate scaling factor (i.e., observed/modeled):
            if et_solar.ho_jm2 != 0:
                sfactor = grid_srad_d/et_solar.ho_jm2
            else:
                sfactor = 1.0
            #
            # Add scaled half-hourly PPFD to dictionary [umol m-2 s-1]:
            start_time = datetime.datetime(2000, 1, 1, 0, 0, 0)
            for val in et_solar.ppfd_hh:
                my_time = "%s" % start_time.time()
                gapfill_dict[my_time] = (val*sfactor)
                start_time = start_time + datetime.timedelta(minutes=30)
            #
            # Add observations to dictionary:
            ppfd_obs = {}
            for i in xrange(len(daily_ts)):
                my_time = "%s" % daily_ts[i].time()
                gapfill_dict[my_time] = daily_ppfd[i]
                ppfd_obs[my_time] = daily_ppfd[i]
            #
            # Save to monthly PPFD time series:
            monthly_ppfd_gapless = numpy.append(
                monthly_ppfd_gapless, 
                [gapfill_dict[x] for x in sorted(gapfill_dict.keys())]
            )
            #
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
                        print "Cannot append to file:", out_file
                    else:
                        f.write("%s,%0.3f,%0.3f\n" % (dt, obs, gfv))
                        f.close()
        else:
            # No gaps; append daily series to monthly
            #   NOTE: drop last entry from daily_ppfd (midnight next day)
            monthly_ppfd_gapless = numpy.append(
                monthly_ppfd_gapless, 
                daily_ppfd[0:-1]
                )
            #
            # Write to file:
            if to_write:
                for i in xrange(len(daily_ts) - 1):
                    dt = "%s" % daily_ts[i]
                    obs = daily_ppfd[i]
                    try:
                        f = open(out_file, 'a')
                    except IOError:
                        print "Cannot append to file:", out_file
                    else:
                        f.write("%s,%0.3f,%0.3f\n" % (dt, obs, obs))
                        f.close()
        #
        # Increment day
        cur_date = add_one_day(cur_date)
        #
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
              - connect_sql
              - get_msvidx
    """
    # Get msvidx value for PPFD:
    ppfd_idx = get_msvidx(station, 'PPFD_f')
    #
    # SQL query parameters:
    params = (ppfd_idx, start_date, add_one_day(start_date))
    #
    # Define query:
    q = (
        "SELECT data_set.datetime, data_set.data "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "AND data_set.datetime BETWEEN DATE %s AND DATE %s "
        "ORDER BY data_set.datetime ASC;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and store results:
    cur.execute(q, params)
    ppfd_vals = numpy.array([])
    time_vals = numpy.array([])
    if cur.rowcount > 0:
        for record in cur:
            time_vals = numpy.append(time_vals, record[0])
            ppfd_vals = numpy.append(ppfd_vals, record[1])
            #
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
    Depends:  connect_sql
    """
    # SQL query params:
    params = (msvidx, time_point)
    #
    # Define SQL query:
    q = (
        "SELECT data_set.data "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "AND data_set.datetime = %s;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
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
        print "No data found in function get_data_point"
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
    Depends:  - connect_sql
              - get_msvidx
    """
    # Get msvidx values for specified station:
    ppfdi = get_msvidx(station, 'PPFD_f')
    neei = get_msvidx(station, 'NEE_f')
    #
    # SQL query parameters:
    params = (ppfdi, neei)
    #
    # Define start date query:
    q1 = (
        "SELECT data_set.datetime "
        "FROM data_set "
        "WHERE data_set.msvidx = %s OR data_set.msvidx = %s "
        "ORDER BY data_set.datetime ASC LIMIT 1;"
        )
    #
    # Define end date query:
    q2 = (
        "SELECT data_set.datetime "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "OR data_set.msvidx = %s "
        "ORDER BY data_set.datetime DESC LIMIT 1;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Get start date from datetime object:
    cur.execute(q1, params)
    sd = cur.fetchone()[0].date()
    #
    # Get end date from datetime object:
    cur.execute(q2, params)
    ed = cur.fetchone()[0].date()
    #
    # Make the starting date begin at day 1
    sd = sd.replace(day=1)
    #
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
    Depends:  connect_sql
    """
    # Query paramters:
    params = (station,)
    #
    # SQL query:
    q = (
        "SELECT met_data.lon, met_data.lat "
        "FROM met_data "
        "WHERE met_data.stationid = %s;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
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
    Depends:  connect_sql
    """
    # Define query:
    q = (
        "SELECT var_list.msvidx "
        "FROM var_list "
        "WHERE var_list.stationid = %s "
        "AND var_list.varname = %s;"
        )
    #
    # SQL query parameters:
    params = (station, variable)
    #
    # Connect to database and star cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and fetch results:
    cur.execute(q, params)
    try:
        result = cur.fetchone()[0]
    except:
        print "Could not return an msvidx value for station", 
        print station, "and variable", variable
        result = ""
    finally:
        con.close()
        return result

def get_pressure(s):
    """
    Name:     get_pressure
    Input:    str, station name (s)
    Output:   float, atmospheric pressure, Pa (patm)
    Features: Returns the atmospheric pressure based on the elevation of a 
              given station
    Depends:  - connect_sql
              - flux_to_grid
              - get_data_point
              - get_msvidx
    Ref:      Allen et al. (1998)
    """
    # Define constants:
    kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo = 298.15   # base temperature, K (Prentice, unpublished)
    kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG = 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    #
    # Define query w/ parameters:
    params = (s,)
    q = (
        "SELECT met_data.ele "
        "FROM met_data "
        "WHERE met_data.stationid = %s;"
        )
    #
    # Connect to database:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return results:
    cur.execute(q, params)
    station_ele = cur.fetchone()[0]
    con.close()
    #
    # Check to see that elevation is valid:
    if float(station_ele) == -9999:
        # Find CRU Elv:
        elv_sd = datetime.date(2006, 6, 1)
        hdg_station = flux_to_grid(s)
        elv_msvidx = get_msvidx(hdg_station, 'Elv')
        elv_data = get_data_point(elv_msvidx, elv_sd)
        station_ele = elv_data
    #
    # Convert elevation to pressure, Pa:
    z = float(station_ele)
    patm = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
    #
    return patm

def get_stations():
    """
    Name:     get_stations
    Input:    None.
    Output:   list, station names (results)
    Features: Returns a list of flux station names from GePiSaT database
    Depends:  connect_sql
    """
    # Define query:
    q = (
        "SELECT stationid "
        "FROM met_data "
        "WHERE dim=0 "
        "AND geom=%s "
        "ORDER BY stationid ASC;"
        )
    #
    params = ("point",)
    #
    # Connect to database and start cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and fetch results:
    cur.execute(q,params)
    results = []
    for record in cur:
        results.append(record[0])  # <- extract record from tuple
        #
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
              - connect_sql
              - get_msvidx
    """
    # Get msvidx values for specified station
    ppfd_idx = get_msvidx(station, 'PPFD_f')
    nee_idx = get_msvidx(station, 'NEE_f')
    #
    # Increment start date one month:
    end_date = add_one_month(start_date)
    #
    # SQL query parameters:
    params = (ppfd_idx, nee_idx, start_date, end_date, ppfd_idx, nee_idx)
    #
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
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
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
                #
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
    #
    ## ~~~~~~~~~~ ANALYZE OBSERVATIONS ~~~~~~~~~~ ##
    # ---------------------
    # Model H optimization:
    # ---------------------
    try:
        mh_opt, mh_cov = curve_fit(
            my_class.model_h, 
            my_class.ppfd_obs, 
            my_class.nee_obs, 
            p0 = my_class.hm_estimates
            )
    except:
        my_class.hm_optimized = [-9999.0, -9999.0, -9999.0]
        my_class.hm_optim_err = [-9999.0, -9999.0, -9999.0]
    else:
        (
            my_class.hm_optimized[0], 
            my_class.hm_optimized[1], 
            my_class.hm_optimized[2]) = mh_opt
        #
        # Extract variance values from matrix:
        try:
            mh_var = numpy.diag(mh_cov)
        except ValueError:
            mh_var = [0,0,0]
        else:
            if numpy.isfinite(mh_var).all() and not (mh_var < 0).any():
                (
                    my_class.hm_optim_err[0], 
                    my_class.hm_optim_err[1], 
                    my_class.hm_optim_err[2]
                    ) = numpy.sqrt(mh_var)
                #
                # Calculate the t-value:
                (
                    my_class.hm_optim_t[0],
                    my_class.hm_optim_t[1],
                    my_class.hm_optim_t[2]
                    ) = (
                        numpy.array(my_class.hm_optimized)/
                        numpy.array(my_class.hm_optim_err)
                        )
                #
                # Calculate the p-value:
                hm_df = len(my_class.nee_obs) - 3.0  # degrees of freedom
                (
                    my_class.hm_optim_p[0],
                    my_class.hm_optim_p[1],
                    my_class.hm_optim_p[2]
                    ) = scipy.stats.t.pdf(
                        -abs(numpy.array(my_class.hm_optim_t)),
                        hm_df)
            else:
                my_class.hm_optim_err[0] = 0.0
                my_class.hm_optim_err[1] = 0.0
                my_class.hm_optim_err[2] = 0.0
    finally:
        my_class.calc_model_h()
        #
    # ---------------------
    # Model L optimization:
    # ---------------------
    try:
        ml_opt, ml_cov = curve_fit(
            my_class.model_l, 
            my_class.ppfd_obs, 
            my_class.nee_obs, 
            p0 = my_class.lm_estimates
            )
    except:
        my_class.lm_optimized = [-9999.0, -9999.0]
        my_class.lm_optim_err = [-9999.0, -9999.0]
    else:
        (
            my_class.lm_optimized[0], 
            my_class.lm_optimized[1]) = ml_opt
        #
        # Extract variance values from matrix:
        try:
            ml_var = numpy.diag(ml_cov)
        except ValueError:
            ml_var = [0,0]
        else:
            if numpy.isfinite(ml_var).all() and not (ml_var < 0).any():
                (
                    my_class.lm_optim_err[0], 
                    my_class.lm_optim_err[1]
                    ) = numpy.sqrt(ml_var)
                #
                # Calculate the t-value:
                (
                    my_class.lm_optim_t[0],
                    my_class.lm_optim_t[1]
                    ) = (
                        numpy.array(my_class.lm_optimized)/
                        numpy.array(my_class.lm_optim_err)
                        )
                #
                # Calculate the p-value:
                lm_df = len(my_class.nee_obs) - 2 # degrees of freedom
                (
                    my_class.lm_optim_p[0],
                    my_class.lm_optim_p[1]
                    ) = scipy.stats.t.pdf(
                        -abs(numpy.array(my_class.lm_optim_t)),
                        lm_df)
            else:
                my_class.lm_optim_err[0] = 0.0
                my_class.lm_optim_err[1] = 0.0
    finally:
        my_class.calc_model_l()
        #
    ## ~~~~~~~~~~ WRITE OBSERVATION RESULTS ~~~~~~~~~~ ##
    if to_write == 1:
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
    if rm_out == 1:
        # ---------------------
        # Model H optimization: 
        # ---------------------
        my_class.remove_mh_outliers()
        try:
            mh_opt_ro, mh_cov_ro = curve_fit(
                my_class.model_h, 
                my_class.ppfd_obs_h_ro, 
                my_class.nee_obs_h_ro, 
                p0 = my_class.hm_estimates
                )
        except:
            my_class.hm_optimized_ro = [-9999.0, -9999.0, -9999.0]
            my_class.hm_optim_err_ro = [-9999.0, -9999.0, -9999.0]
        else:
            (
                my_class.hm_optimized_ro[0], 
                my_class.hm_optimized_ro[1], 
                my_class.hm_optimized_ro[2]) = mh_opt_ro
            #
            # Extract variance values from matrix:
            try:
                mh_var_ro = numpy.diag(mh_cov_ro)
            except ValueError:
                mh_var_ro = [0,0,0]
            else:
                if numpy.isfinite(mh_var_ro).all() and not (mh_var_ro < 0).any():
                    (
                        my_class.hm_optim_err_ro[0], 
                        my_class.hm_optim_err_ro[1], 
                        my_class.hm_optim_err_ro[2]
                        ) = numpy.sqrt(mh_var_ro)
                    #
                    # Calculate the t-value:
                    (
                        my_class.hm_optim_t_ro[0],
                        my_class.hm_optim_t_ro[1],
                        my_class.hm_optim_t_ro[2]
                        ) = (
                            numpy.array(my_class.hm_optimized_ro)/
                            numpy.array(my_class.hm_optim_err_ro)
                            )
                    #
                    # Calculate the p-value:
                    hm_df_ro = len(my_class.nee_obs) - my_class.hm_outliers - 3
                    (
                        my_class.hm_optim_p_ro[0],
                        my_class.hm_optim_p_ro[1],
                        my_class.hm_optim_p_ro[2]
                        ) = scipy.stats.t.pdf(
                            -abs(numpy.array(my_class.hm_optim_t_ro)),
                            hm_df_ro)
                else:
                    my_class.hm_optim_err_ro[0] = 0.0
                    my_class.hm_optim_err_ro[1] = 0.0
                    my_class.hm_optim_err_ro[2] = 0.0
        finally:
            my_class.calc_model_h_ro()
            #
        # ---------------------
        # Model L optimization:
        # --------------------- 
        my_class.remove_ml_outliers()
        try:
            ml_opt_ro, ml_cov_ro = curve_fit(
                my_class.model_l, 
                my_class.ppfd_obs_l_ro, 
                my_class.nee_obs_l_ro, 
                p0 = my_class.lm_estimates
                )
        except:
            my_class.lm_optimized_ro = [-9999.0, -9999.0]
            my_class.lm_optim_err_ro = [-9999.0, -9999.0]
        else:
            (
                my_class.lm_optimized_ro[0], 
                my_class.lm_optimized_ro[1]) = ml_opt_ro
            #
            # Extract variance values from matrix:
            try:
                ml_var_ro = numpy.diag(ml_cov_ro)
            except ValueError:
                ml_var_ro = [0,0]
            else:
                if numpy.isfinite(ml_var_ro).all() and not (ml_var_ro < 0).any():
                    (
                        my_class.lm_optim_err_ro[0], 
                        my_class.lm_optim_err_ro[1]
                        ) = numpy.sqrt(ml_var_ro)
                    #
                    # Calculate the t-value:
                    (
                        my_class.lm_optim_t_ro[0],
                        my_class.lm_optim_t_ro[1]
                        ) = (
                            numpy.array(my_class.lm_optimized_ro)/
                            numpy.array(my_class.lm_optim_err_ro)
                            )
                    #
                    # Calculate the p-value:
                    lm_df_ro = len(my_class.nee_obs) - my_class.lm_outliers - 2
                    (
                        my_class.lm_optim_p_ro[0],
                        my_class.lm_optim_p_ro[1]
                        ) = scipy.stats.t.pdf(
                            -abs(numpy.array(my_class.lm_optim_t_ro)),
                            lm_df_ro)
                else:
                    my_class.lm_optim_err_ro[0] = 0.0
                    my_class.lm_optim_err_ro[1] = 0.0
        finally:
            my_class.calc_model_l_ro()
            #
        ## ~~~~~~~~~~ WRITE OUTLIER-FREE RESULTS ~~~~~~~~~~ ##
        if to_write == 1:
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
            #
    # Run model selection criteria:
    my_class.model_selection()
    #
    # Return the FLUX_PARTI class object:
    return my_class

################################################################################
## MAIN PROGRAM
################################################################################
# Get list of all flux station names:
#stations = get_stations()
stations = ['CZ-wet',]

# Initialize summary statistics file:
summary_file = "out/summary_statistics.txt"
summary_file_init(summary_file)

# Create/initialize LUE class instance:
my_lue = LUE()
lue_out_file = "out/LUE_All-Stations.txt"

# Iterate through stations:
for station in stations:
    # Initialize station's LUE & daily GPP output file:
    lue_file = "out/%s_%s.txt" % (station, "LUE")
    gpp_file = "out/%s_%s.txt" % (station, "GPP-daily")
    #
    # Get first/last dates for station data:
    sd, ed = get_dates(station)
    #
    # Get flux station's corresponding 0.5-degree grid station:
    hdg_station = flux_to_grid(station)
    #
    # Get station's variable msvidx values:
    co2_msvidx = "US-MLO.21"
    cpa_msvidx = get_msvidx(hdg_station, 'alpha')
    fpar_msvidx = get_msvidx(hdg_station, 'FAPAR')
    tair_msvidx = get_msvidx(hdg_station, 'Tc')
    vpd_msvidx = get_msvidx(hdg_station, 'VPD')
    #
    # Calculate station's atmospheric pressure based on elevation:
    patm = get_pressure(station)
    #
    # Process each month in time:
    while sd < ed:
        # Get PPFD and NEE array pairs [umol m-2 s-1]:
        #   NOTE: accepts in opposite order as sent
        monthly_nee, monthly_ppfd = monthly_ppfd_nee(station, sd)
        #
        # Process if enough data was found:
        if (len(monthly_ppfd) > 3 and len(monthly_nee) > 3):
            # Perform GPP partitioning (returns FLUX_PARTI class object):
            monthly_parti = partition(monthly_nee, 
                                      monthly_ppfd, 
                                      to_write=0, 
                                      rm_out=1,
                                      tower=station,
                                      month=sd)
            #
            # Perform half-hourly PPFD gapfilling (umol m-2 s-1):
            (gf_time, gf_ppfd) = gapfill_ppfd(station, sd, to_write=1)
            #
            # Calculate half-hourly GPP (umol m-2 s-1)
            gf_gpp, gf_gpp_err = monthly_parti.calc_gpp(gf_ppfd)
            #
            # @TODO: daily GPP calculation & writeout
            #
            # Continue processing if partitioning was successful:
            if monthly_parti.mod_select > 0:
                #
                # The new LUE model:
                # GPP = f(PPFD, fAPAR, VPD, CPA, Tair, Patm, CO2)
                #       monthly variables: PPFD, fAPAR, CPA, VPD, Tair
                #       annual variables: CO2
                #       constants in time: Patm=f(elevation)
                #
                # Integrate PPFD & GPP [umol m-2]; dt=30 min (1800 s)
                ppfd_month = simpson(gf_ppfd.clip(min=0), 1800)
                gpp_month = simpson(gf_gpp.clip(min=0), 1800)
                gpp_month_err = simpson(gf_gpp_err, 1800)
                #
                # Convert units from [umol m-2] to [mol m-2]:
                ppfd_month = (1e-6)*ppfd_month
                gpp_month = (1e-6)*gpp_month
                gpp_month_err = (1e-6)*gpp_month_err
                #
                # Retrieve annual CO2:
                annual_sd = sd.replace(month=1)
                co2_annual = get_data_point(co2_msvidx, annual_sd) # ppm
                #
                # Retrieve monthly gridded data:
                cpa_month = get_data_point(cpa_msvidx, sd)         # unitless
                fpar_month = get_data_point(fpar_msvidx, sd)       # unitless
                tair_month = get_data_point(tair_msvidx, sd)       # deg C
                vpd_month = get_data_point(vpd_msvidx, sd)         # kPa
                #
                # Add LUE parameters to LUE class:
                my_lue.add_station_val(station, 
                                       sd, 
                                       gpp_month,                  # mol/m2
                                       gpp_month_err,              # mol/m2
                                       fpar_month,                 # unitless
                                       ppfd_month,                 # mol/m2
                                       vpd_month,                  # kPa
                                       cpa_month,                  # unitless
                                       tair_month,                 # degC
                                       co2_annual,                 # ppm
                                       patm)                       # Pa
            #
        else:
            # Create an 'empty' class:
            monthly_parti = FLUX_PARTI(monthly_ppfd, monthly_nee, station, sd)
        #
        # Save class summary statistics:
        SFILE = open(summary_file, 'a')
        SFILE.write(monthly_parti.summary_statistics())
        SFILE.close()
        #
        # Increment date:
        sd = add_one_month(sd)
    #
    # Write monthly LUE parameters to file:
    my_lue.write_out_val(station, lue_file)
    #
    # Calculate LUE for station:
    calc_lue(my_lue, station)

# Write station LUE to file
my_lue.write_out_lue(lue_out_file)
