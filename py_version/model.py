#!/usr/bin/python
#
# model.py
#
# VERSION 2.2.0-dev
#
# LAST UPDATED: 2016-07-06
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
#    --> due to the order of the msvidx (NEE comes before PPFD) [16.05.20]
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
#  - moved partition function to STAGE1 class [16.04.01]
#  - separated function from init for SOLAR class [16.04.01]
#  - up to calculating GPP and daily integrations
#
# -----
# todo:
# -----
# x Create a constants file
# x Separate utility functions to their own file & import them
# x Add logging
# x Consider creating a data class to handle the interfact between the main
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
import os.path

from data import DATA
from flux_parti import FLUX_PARTI
from lue import LUE
from utilities import add_one_month
from utilities import simpson


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == "__main__":
    # Initialize data class for handling all data IO operations, which
    # includes connecting to the GePiSaT database:
    my_data = DATA()

    # Define output directory:
    output_dir = os.path.join(
        os.path.expanduser("~"), "Desktop", "temp", "out")
    my_data.set_output_directory(output_dir)

    # Get list of all flux station names:
    # my_data.find_stations()         # use this function to search the DB
    my_data.stations = ['CZ-wet', ]   # or manually define stations here

    # Initialize summary statistics file:
    my_data.create_summary_file("summary_statistics.txt")

    # Create/initialize LUE class instance:
    my_lue = LUE()

    # Iterate through stations:
    for station in my_data.stations:
        # Initialize station data:
        my_data.set_working_station(station)

        # Process each month in time:
        sd = my_data.start_date
        ed = my_data.end_date
        while sd < ed:
            # Get PPFD and NEE array pairs [umol m-2 s-1]:
            monthly_nee, monthly_ppfd = my_data.find_monthly_nee_ppfd(sd)

            # Initialize a FLUX_PARTI class:
            my_parter = FLUX_PARTI(monthly_ppfd, monthly_nee, station, sd)

            # Process if enough data was found:
            if (len(monthly_ppfd) > 3 and len(monthly_nee) > 3):

                # Perform GPP partitioning:
                my_parter.partition(outliers=True)

                # Perform half-hourly PPFD gapfilling (umol m-2 s-1):
                (gf_time, gf_ppfd) = my_data.gapfill_monthly_ppfd(
                    sd, to_save=True)

                # Calculate half-hourly GPP (umol m-2 s-1)
                gf_gpp, gf_gpp_err = my_parter.calc_gpp(gf_ppfd)

                # Calculate daily GPP (mol m-2):
                if False:
                    gpp_daily = my_data.sub_to_daily_gpp(
                        gf_time, gf_gpp, gf_gpp_err, to_save=True)

                # Continue processing if partitioning was successful:
                if my_parter.best_model > 0:
                    # The new LUE model:
                    # GPP = f(PPFD, fAPAR, VPD, CPA, Tair, Patm, CO2)
                    #  - monthly variables: PPFD, fAPAR, CPA, VPD, Tair
                    #  - annual variables: CO2
                    #  - constants in time: Patm=f(elevation)

                    # !!!!! STOPPED HERE - @TODO - !!!!!!!

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
                    co2_msvidx = "US-MLO.21"
                    cpa_msvidx = get_msvidx(hdg_station, 'alpha')
                    fpar_msvidx = get_msvidx(hdg_station, 'FAPAR')
                    tair_msvidx = get_msvidx(hdg_station, 'Tc')
                    vpd_msvidx = get_msvidx(hdg_station, 'VPD')
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

            # Save class summary statistics:
            SFILE = open(my_data.summary_file, 'a')
            SFILE.write(my_parter.summary_statistics())
            SFILE.close()

            # Increment date:
            sd = add_one_month(sd)

        # Write monthly LUE parameters to file:
        my_lue.write_out_val(station, lue_file)
