#!/usr/bin/python
#
# utilities.py
#
# LAST UPDATED: 2016-06-10
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
import scipy.stats

from const import kPo
from const import kTo
from const import kL
from const import kG
from const import kR
from const import kMa
from const import pir


###############################################################################
# FUNCTIONS
###############################################################################
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


def calc_statistics(my_array):
    """
    Name:     calc_statistics
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

    # Make sure my_array is not empty or crashes on skew/kurt:
    if my_array.any() and len(my_array) > 1:
        # Max, min, mean, st dev, skew, kurtosis (offset from normal)
        max_val = my_array.max()
        min_val = my_array.min()
        ave_val = my_array.mean()
        std_val = my_array.std()

        # Address divide by zero issues:
        if std_val == 0:
            std_val = 1e-4

        skew_val = (
            sum((my_array - ave_val)**3)/((len(my_array) - 1)*std_val**3))
        kurt_val = (
            sum((my_array - ave_val)**4)/((len(my_array) - 1)*std_val**4) - 3)
    else:
        logging.warning("Insufficient data for stats calculation")
        # Maintain initial quantity values:
        max_val = 0.0
        min_val = 0.0
        ave_val = 0.0
        std_val = 0.0
        skew_val = 0.0
        kurt_val = 0.0

    return (max_val, min_val, ave_val, std_val, skew_val, kurt_val)


def dcos(x):
    """
    Name:     dcos
    Input:    float/nd.array, angle, degrees (x)
    Output:   float/nd.array, cos(x*pi/180)
    Features: Calculates the cosine of an array of angles in degrees
    """
    if isinstance(x, float):
        logging.debug("calculating cosine of %f degrees", x)
    elif isinstance(x, numpy.ndarray):
        logging.debug("calculating cosine of numpy array of length %d", x.size)
    return numpy.cos(x*pir)


def dsin(x):
    """
    Name:     dsin
    Input:    float/nd.array, angle, degrees (x)
    Output:   float/nd.array, sin(x*pi/180)
    Features: Calculates the sine of an array of angles in degrees
    """
    if isinstance(x, float):
        logging.debug("calculating sine of %f degrees", x)
    elif isinstance(x, numpy.ndarray):
        logging.debug("calculating cosine of numpy array of length %d", x.size)
    return numpy.sin(x*pir)


def elv2pres(z):
    """
    Name:     elv2pres
    Input:    float, elevation above sea level (z), m
    Output:   float, atmospheric pressure, Pa
    Features: Calculates atm. pressure for a given elevation
    Depends:  Global constants:
              - kPo, base press.   - kTo, base temp.
              - kL, lapse rate     - kMa, mole. wt. of air
              - kG, std. gravity   - kR, gas const.
    Ref:      Allen et al. (1998)
    """
    logging.debug("estimating atmospheric pressure at %f m", z)
    p = kPo*(1.0 - kL*z/kTo)**(kG*kMa/(kR*kL))
    return p


def get_hm_estimates(max_n, min_n, std_n, max_p, min_p):
    """
    Name:     get_hm_estimates
    Inputs:   - float, maximum NEE value (max_n)
              - float, minimum NEE value (min_n)
              - float, std. dev. of NEE values (std_n)
              - float, maximum PPFD value (max_p)
              - float, minimum PPFD value (min_p)
    Outputs:  tuple of estimates
    Features: Returns estimate values for hyperbolic model fitting based
              on given data statistics
    """
    hm_foo = 3.83*std_n
    hm_r = 0.69*std_n
    try:
        hm_alpha = 1.96*(max_n - min_n)/(max_p - min_p)
    except ZeroDivisionError:
        hm_alpha = 1.96*(max_n - min_n)/1.0e-3
    else:
        if abs(hm_alpha) <= 5.0e-4:
            hm_alpha = 1.0e-3

    return (hm_foo, hm_alpha, hm_r)


def get_lm_estimates(max_n, min_n, ave_n, std_n, max_p, min_p, ave_p, std_p):
    """
    Name:     get_lm_estimates
    Inputs:   - float, maximum NEE value (max_n)
              - float, minimum NEE value (min_n)
              - float, average NEE value (ave_n)
              - float, std. deviation of NEE values (std_n)
              - float, maximum PPFD value (max_p)
              - float, minimum PPFD value (min_p)
              - float, average PPFD value (ave_p)
              - float, std. deviation of PPFD values (std_p)
    Outputs:  tuple of estimates
    Features: Returns estimate values for linear model fitting based on
              given data statistics
    """
    lm_r = 0.899*std_n + 0.827*ave_n + 0.00628*ave_p - 0.008*std_p
    try:
        lm_alpha = 0.672*(max_n - min_n)/(max_p - min_p)
    except ZeroDivisionError:
        lm_alpha = 0.672*(max_n - min_n)/1.0e-3
    else:
        if abs(lm_alpha) <= 5.0e-4:
            lm_alpha = 1.0e-3

    return (lm_alpha, lm_r)


def get_outliers(fit, obs, np):
    """
    Name:     get_outliers
    Input:    - numpy.ndarray, model fitted values (fit)
              - numpy.ndarray, observations (obs)
              - int, number of parameters (np)
    Output:   tuple of indexes
    Features: Returns indexes of outliers found using Peirce's criterion
    Depends:  - goodness_of_fit
              - peirce_dev
    """
    # Calculate mean-square error:
    mse, rmse, rsqr = goodness_of_fit(fit, obs, np)

    # Calculate the square errors:
    sq_errors = (obs - fit)**2.0

    # Set Peirce values:
    peirce_cap_n = len(obs)
    peirce_lc_n = 1
    peirce_m = np

    # Calculate tolerance
    peirce_x2 = peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
    peirce_delta2 = mse*peirce_x2

    # Find if/where exceedance occurs:
    outliers_index = numpy.where(sq_errors > peirce_delta2)
    outliers_found = len(outliers_index[0])

    # Run check again if no outliers are found in first attempt:
    if (outliers_found == 0):
        peirce_lc_n = 2
        peirce_x2 = peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
        peirce_delta2 = mse*peirce_x2
        outliers_index = numpy.where(sq_errors > peirce_delta2)
        outliers_found = len(outliers_index[0])

        # Reset n
        peirce_lc_n = 1

    # Increment n until it is greater than number of outliers found:
    while (peirce_lc_n <= outliers_found):
        peirce_lc_n += 1

        # Check that n < N:
        if peirce_lc_n >= peirce_cap_n:
            peirce_lc_n = outliers_found + 1.0
        else:
            peirce_x2 = peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m)
            peirce_delta2 = mse*peirce_x2
            outliers_index = numpy.where(sq_errors > peirce_delta2)
            outliers_found = len(outliers_index[0])

    logging.debug("Found %d outliers", outliers_found)
    return outliers_index


def goodness_of_fit(fit, obs, nparams):
    """
    Name:     goodness_of_fit
    Input:    - numpy.ndarray, modeled values (fit)
              - numpy.ndarray, observed values (obs)
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

    # Check that both have the same number of values and that the length
    # is greater than 4, no divide by zero issues:
    if len(obs) > 4 and len(fit) == len(obs):
        n = len(obs)

        # Sum of the squared error (SSE):
        sse = sum((obs - fit)**2.0)

        # Mean squared error:
        mse = float(sse)/(n - nparams)

        # Total sum of the squares (SST):
        sst = sum((obs - float(sum(obs))/n)**2.0)

        # R-squared:
        # r2 = 1.0 - float(sse)/sst
        r2_adj = 1.0 - (float(sse)/sst*(n - 1.0) / float(n - nparams - 1.0))

        # RMSE:
        rmse = numpy.sqrt(float(sse)/n)

    return (mse, rmse, r2_adj)


def grid_centroid(my_lon, my_lat, grid_res=0.5):
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
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = int(180./grid_res)
    lon_dim = int(360./grid_res)
    lats = [lat_min + y*grid_res for y in range(lat_dim)]
    lons = [lon_min + x*grid_res for x in range(lon_dim)]

    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index - 1]
        try:
            bb_lon_max = lons[lon_index + 1]
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
        bb_lat_min = lats[lat_index - 1]
        try:
            bb_lat_max = lats[lat_index + 1]
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
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0)
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0)
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0)
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


def init_summary_dict():
    """
    Name:     init_summary_dict
    Inputs:   None.
    Outputs:  dict, summary fields
    Features: Returns an empty summary dictionary
    """
    summary_dict = {
        1: {"name": "name", "val": ""},
        2: {"name": "month", "val": datetime.date(1999, 1, 1)},
        3: {"name": "n_obs", "val": 0},
        4: {"name": "n_h", "val": 0},
        5: {"name": "n_l", "val": 0},
        6: {"name": "foo_est_obs_h", "val": -9999.0},
        7: {"name": "foo_opt_obs_h", "val": -9999.0},
        8: {"name": "foo_err_obs_h", "val": 0.0},
        9: {"name": "foo_t_obs_h", "val": 0.0},
        10: {"name": "foo_p_obs_h", "val": 0.0},
        11: {"name": "foo_est_ro_h", "val": -9999.0},
        12: {"name": "foo_opt_ro_h", "val": -9999.0},
        13: {"name": "foo_err_ro_h", "val": 0.0},
        14: {"name": "foo_t_ro_h", "val": 0.0},
        15: {"name": "foo_p_ro_h", "val": 0.0},
        16: {"name": "alpha_est_obs_h", "val": -9999.0},
        17: {"name": "alpha_opt_obs_h", "val": -9999.0},
        18: {"name": "alpha_err_obs_h", "val": 0.0},
        19: {"name": "alpha_t_obs_h", "val": 0.0},
        20: {"name": "alpha_p_obs_h", "val": 0.0},
        21: {"name": "alpha_est_ro_h", "val": -9999.0},
        22: {"name": "alpha_opt_ro_h", "val": -9999.0},
        23: {"name": "alpha_err_ro_h", "val": 0.0},
        24: {"name": "alpha_t_ro_h", "val": 0.0},
        25: {"name": "alpha_p_ro_h", "val": 0.0},
        26: {"name": "alpha_est_obs_l", "val": -9999.0},
        27: {"name": "alpha_opt_obs_l", "val": -9999.0},
        28: {"name": "alpha_err_obs_l", "val": 0.0},
        29: {"name": "alpha_t_obs_l", "val": 0.0},
        30: {"name": "alpha_p_obs_l", "val": 0.0},
        31: {"name": "alpha_est_ro_l", "val": -9999.0},
        32: {"name": "alpha_opt_ro_l", "val": -9999.0},
        33: {"name": "alpha_err_ro_l", "val": 0.0},
        34: {"name": "alpha_t_ro_l", "val": 0.0},
        35: {"name": "alpha_p_ro_l", "val": 0.0},
        36: {"name": "r_est_obs_h", "val": -9999.0},
        37: {"name": "r_opt_obs_h", "val": -9999.0},
        38: {"name": "r_err_obs_h", "val": 0.0},
        39: {"name": "r_t_obs_h", "val": 0.0},
        40: {"name": "r_p_obs_h", "val": 0.0},
        41: {"name": "r_est_ro_h", "val": -9999.0},
        42: {"name": "r_opt_ro_h", "val": -9999.0},
        43: {"name": "r_err_ro_h", "val": 0.0},
        44: {"name": "r_t_ro_h", "val": 0.0},
        45: {"name": "r_p_ro_h", "val": 0.0},
        46: {"name": "r_est_obs_l", "val": -9999.0},
        47: {"name": "r_opt_obs_l", "val": -9999.0},
        48: {"name": "r_err_obs_l", "val": 0.0},
        49: {"name": "r_t_obs_l", "val": 0.0},
        50: {"name": "r_p_obs_l", "val": 0.0},
        51: {"name": "r_est_ro_l", "val": -9999.0},
        52: {"name": "r_opt_ro_l", "val": -9999.0},
        53: {"name": "r_err_ro_l", "val": 0.0},
        54: {"name": "r_t_ro_l", "val": 0.0},
        55: {"name": "r_p_ro_l", "val": 0.0},
        56: {"name": "r2_obs_h", "val": -9999.0},
        57: {"name": "r2_ro_h", "val": -9999.0},
        58: {"name": "rmse_obs_h", "val": -9999.0},
        59: {"name": "rmse_ro_h", "val": -9999.0},
        60: {"name": "r2_obs_l", "val": -9999.0},
        61: {"name": "r2_ro_l", "val": -9999.0},
        62: {"name": "rmse_obs_l", "val": -9999.0},
        63: {"name": "rmse_ro_l", "val": -9999.0},
        64: {"name": "min_ppfd_obs", "val": 0.0},
        65: {"name": "max_ppfd_obs", "val": 0.0},
        66: {"name": "ave_ppfd_obs", "val": 0.0},
        67: {"name": "std_ppfd_obs", "val": 0.0},
        68: {"name": "skw_ppfd_obs", "val": 0.0},
        69: {"name": "krt_ppfd_obs", "val": 0.0},
        70: {"name": "min_ppfd_ro_h", "val": 0.0},
        71: {"name": "max_ppfd_ro_h", "val": 0.0},
        72: {"name": "ave_ppfd_ro_h", "val": 0.0},
        73: {"name": "std_ppfd_ro_h", "val": 0.0},
        74: {"name": "skw_ppfd_ro_h", "val": 0.0},
        75: {"name": "krt_ppfd_ro_h", "val": 0.0},
        76: {"name": "min_ppfd_ro_l", "val": 0.0},
        77: {"name": "max_ppfd_ro_l", "val": 0.0},
        78: {"name": "ave_ppfd_ro_l", "val": 0.0},
        79: {"name": "std_ppfd_ro_l", "val": 0.0},
        80: {"name": "skw_ppfd_ro_l", "val": 0.0},
        81: {"name": "krt_ppfd_ro_l", "val": 0.0},
        82: {"name": "min_nee_obs", "val": 0.0},
        83: {"name": "max_nee_obs", "val": 0.0},
        84: {"name": "ave_nee_obs", "val": 0.0},
        85: {"name": "std_nee_obs", "val": 0.0},
        86: {"name": "skw_nee_obs", "val": 0.0},
        87: {"name": "krt_nee_obs", "val": 0.0},
        88: {"name": "min_nee_ro_h", "val": 0.0},
        89: {"name": "max_nee_ro_h", "val": 0.0},
        90: {"name": "ave_nee_ro_h", "val": 0.0},
        91: {"name": "std_nee_ro_h", "val": 0.0},
        92: {"name": "skw_nee_ro_h", "val": 0.0},
        93: {"name": "krt_nee_ro_h", "val": 0.0},
        94: {"name": "min_nee_ro_l", "val": 0.0},
        95: {"name": "max_nee_ro_l", "val": 0.0},
        96: {"name": "ave_nee_ro_l", "val": 0.0},
        97: {"name": "std_nee_ro_l", "val": 0.0},
        98: {"name": "skw_nee_ro_l", "val": 0.0},
        99: {"name": "krt_nee_ro_l", "val": 0.0},
        100: {"name": "pearson_r_obs", "val": 0.0},
        101: {"name": "pearson_r_ro_h", "val": 0.0},
        102: {"name": "pearson_r_ro_l", "val": 0.0},
        103: {"name": "model_select", "val": None}
    }
    return summary_dict


def hyp_model(x, Foo, alpha, R):
    """
    Name:     hyp_model
    Inputs:   - numpy.ndarray, monthly PPFD (x)
              - float, hyperbolic parameter (Foo)
              - float, hyperbolic parameter (alpha)
              - float, hyperbolic parameter (R)
    Outputs:  numpy.ndarray, modeled NEE
    Features: Returns array of NEE based on the hyperbolic flux partitioning
              model in the form: f(x) = (ax + b)/(bx + c)
    Ref:      Eq. 15, GePiSaT Documentation
    """
    a = 1.0*(alpha*R) - 1.0*(alpha*Foo)
    b = 1.0*(Foo*R)
    c = 1.0*alpha
    d = 1.0*Foo

    numerator = 1.0*(a*x) + b
    denominator = 1.0*(c*x) + d

    # Compensate for zero division by adding a small number to each value:
    denominator[numpy.where(denominator == 0)] += 1.0e-6

    return (numerator/denominator)


def lin_model(x, alpha, R):
    """
    Name:     lin_model
    Input:    - numpy.ndarray, monthly PPFD (x)
              - float, linear parameter (alpha)
              - float, linear parameter (R)
    Output:   - numpy.ndarray, modeled NEE
    Features: Returns array of NEE based on the linear flux partitioning model
              in the form: f(x) = ax + b
    Ref:      Eq. 14, GePiSaT Documentation
    """
    a = -1.0*alpha
    b = R
    return (a*x + b)


def pearsons_r(x, y):
    """
    Name:     pearsons_r
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

    # Initialize Pearson's r:
    pearsonsr = -9999.0

    # Make certain both arrays have enough points for regression and are of
    # equal lengths:
    if len(x) > 2 and len(x) == len(y):
        try:
            slope, intrcp, pearsonsr, p, sterr = scipy.stats.linregress(x, y)
        except:
            logging.exception("Error calculating Pearson's r")
        else:
            logging.debug("Pearson's r = %f", pearsonsr)

    return pearsonsr


def peirce_dev(peirce_cap_n, peirce_lc_n, peirce_m):
    """
    Name:     peirce_dev
    Input:    - int, total number of observations (peirce_cap_n)
              - int, number of outliers to be removed (peirce_lc_n)
              - int, number of model unknowns (peirce_m)
    Output:   float, squared error threshold (x2)
    Features: Returns the squared threshold error deviation for outlier
              identification using Peirce's criterion based on Gould's
              methodology
    """
    # Assign floats to input variables:
    N = float(peirce_cap_n)
    n = float(peirce_lc_n)
    m = float(peirce_m)

    # Check the total number of observations:
    if N > 1:
        # Calculate Q (Nth root of Gould's equation B):
        # Note: 1/N exponent is factored to each individual term to prevent
        # OverflowError with large N (e.g., >142)
        Q = (n**(n/N)*(N-n)**((N-n)/N))/N

        # Initialize R values (as floats):
        Rnew = 1.0  # <- Tried values between 1 and 10 and all seem stable
        Rold = 0.0  # <- Necessary to prompt while loop

        while (abs(Rnew - Rold) > (N*2.0e-16)):
            # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
            ldiv = Rnew**n
            if ldiv == 0:
                ldiv = 1.0e-6
            Lamda = ((Q**N)/(ldiv))**(1.0/(N - n))

            # Calculate x-squared (straight-forward Gould's equation C):
            x2 = 1.0 + (N - m - n)/n*(1.0 - Lamda**2.0)

            # Return 0 for negative x2 values:
            if x2 < 0:
                x2 = 0
                Rold = Rnew
            else:
                # Use x-squared to update R (Gould's equation D):
                Rold = Rnew
                Rnew = (
                    numpy.exp((x2 - 1)/2.0) *
                    scipy.special.erfc(numpy.sqrt(x2)/numpy.sqrt(2.0))
                    )
    else:
        x2 = 0.0
    return x2
