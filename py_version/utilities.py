#!/usr/bin/python
#
# utilities.py
#
# LAST UPDATED: 2016-05-13
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

from const import pir


###############################################################################
# FUNCTIONS
###############################################################################
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
            sum((my_array - ave_val)**3) /
            ((len(my_array) - 1)*std_val**3)
            )
        kurt_val = (
            sum((my_array - ave_val)**4) /
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


def goodness_of_fit(modvals, obsvals, nparams):
    """
    Name:     goodness_of_fit
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

    # Check that both have the same number of values and that the length
    # is greater than 4, no divide by zero issues:
    if len(obsvals) > 4 and len(modvals) == len(obsvals):
        # Sum of the squared error (SSE):
        sse = sum((obsvals - modvals)**2.0)

        # Mean squared error:
        mse = float(sse)/(len(obsvals) - nparams)

        # Total sum of the squares (SST):
        sst = sum((obsvals - float(sum(obsvals))/len(obsvals))**2.0)

        # R-squared:
        # r2 = 1.0 - float(sse)/sst
        r2_adj = 1.0 - (
            float(sse)/sst*(len(obsvals) - 1.0) /
            float(len(obsvals) - nparams - 1.0)
            )

        # RMSE:
        rmse = numpy.sqrt(float(sse)/len(obsvals))

    return (mse, rmse, r2_adj)


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

    @TODO:    write getters and setters in flux_parti for this dictionary
    """
    summary_dict = {
        1: {"name": "name", "val": ""},
        2: {"name": "month", "val": datetime.date(1999, 1, 1)},
        3: {"name": "n_obs", "val": 0},
        4: {"name": "n_h", "val": 0},
        5: {"name": "n_l", "val": 0},
        6: {"name": "foo_est_obs_h", "val": -9999.0},
        7: {"name": "foo_opt_obs_h", "val": 1.0},
        8: {"name": "foo_err_obs_h", "val": 0.0},
        9: {"name": "foo_t_obs_h", "val": 0.0},
        10: {"name": "foo_p_obs_h", "val": 0.0},
        11: {"name": "foo_est_ro_h", "val": -9999.0},
        12: {"name": "foo_opt_ro_h", "val": 1.0},
        13: {"name": "foo_err_ro_h", "val": 0.0},
        14: {"name": "foo_t_ro_h", "val": 0.0},
        15: {"name": "foo_p_ro_h", "val": 0.0},
        16: {"name": "alpha_est_obs_h", "val": -9999.0},
        17: {"name": "alpha_opt_obs_h", "val": 1.0},
        18: {"name": "alpha_err_obs_h", "val": 0.0},
        19: {"name": "alpha_t_obs_h", "val": 0.0},
        20: {"name": "alpha_p_obs_h", "val": 0.0},
        21: {"name": "alpha_est_ro_h", "val": -9999.0},
        22: {"name": "alpha_opt_ro_h", "val": 1.0},
        23: {"name": "alpha_err_ro_h", "val": 0.0},
        24: {"name": "alpha_t_ro_h", "val": 0.0},
        25: {"name": "alpha_p_ro_h", "val": 0.0},
        26: {"name": "alpha_est_obs_l", "val": -9999.0},
        27: {"name": "alpha_opt_obs_l", "val": 1.0},
        28: {"name": "alpha_err_obs_l", "val": 0.0},
        29: {"name": "alpha_t_obs_l", "val": 0.0},
        30: {"name": "alpha_p_obs_l", "val": 0.0},
        31: {"name": "alpha_est_ro_l", "val": -9999.0},
        32: {"name": "alpha_opt_ro_l", "val": 1.0},
        33: {"name": "alpha_err_ro_l", "val": 0.0},
        34: {"name": "alpha_t_ro_l", "val": 0.0},
        35: {"name": "alpha_p_ro_l", "val": 0.0},
        36: {"name": "r_est_obs_h", "val": -9999.0},
        37: {"name": "r_opt_obs_h", "val": 1.0},
        38: {"name": "r_err_obs_h", "val": 0.0},
        39: {"name": "r_t_obs_h", "val": 0.0},
        40: {"name": "r_p_obs_h", "val": 0.0},
        41: {"name": "r_est_ro_h", "val": -9999.0},
        42: {"name": "r_opt_ro_h", "val": 1.0},
        43: {"name": "r_err_ro_h", "val": 0.0},
        44: {"name": "r_t_ro_h", "val": 0.0},
        45: {"name": "r_p_ro_h", "val": 0.0},
        46: {"name": "r_est_obs_l", "val": -9999.0},
        47: {"name": "r_opt_obs_l", "val": 1.0},
        48: {"name": "r_err_obs_l", "val": 0.0},
        49: {"name": "r_t_obs_l", "val": 0.0},
        50: {"name": "r_p_obs_l", "val": 0.0},
        51: {"name": "r_est_ro_l", "val": -9999.0},
        52: {"name": "r_opt_ro_l", "val": 1.0},
        53: {"name": "r_err_ro_l", "val": 0.0},
        54: {"name": "r_t_ro_l", "val": 0.0},
        55: {"name": "r_p_ro_l", "val": 0.0},
        56: {"name": "r2_obs_h", "val": 0.0},
        57: {"name": "r2_ro_h", "val": 0.0},
        58: {"name": "rmse_obs_h", "val": 0.0},
        59: {"name": "rmse_ro_h", "val": 0.0},
        60: {"name": "r2_obs_l", "val": 0.0},
        61: {"name": "r2_ro_l", "val": 0.0},
        62: {"name": "rmse_obs_l", "val": 0.0},
        63: {"name": "rmse_ro_l", "val": 0.0},
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
        103: {"name": "model_select", "val": 0}
    }
    return summary_dict


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

    # Make certain both arrays have equal lengths:
    if len(x) == len(y):
        try:
            slope, intrcp, pearsonsr, p, sterr = scipy.stats.linregress(x, y)
        except:
            logging.exception("Error calculating Pearson's r")
            return -9999.0
        else:
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
