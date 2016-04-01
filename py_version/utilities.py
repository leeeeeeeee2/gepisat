#!/usr/bin/python
#
# utilities.py
#
# LAST UPDATED: 2016-04-01
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
import logging

import numpy
import scipy.stats

from utilities import pir


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
