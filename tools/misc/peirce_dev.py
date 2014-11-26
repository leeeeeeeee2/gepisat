#!/usr/bin/python
#
# peirce_dev.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-07-15 -- created
# 2013-10-20 -- last updated
#
# -----------
# definition:
# -----------
# This script calculates the critical deviation of errors from a sample of N 
# observations with k outliers and m unknowns.
#
# This script is based on the algorithms given by Gould 1855 and is in-part
# based on the pseudocode by K. Thomsen
# [http://mathforum.org/kb/thread.jspa?forumID=13&threadID=1841790&messageID=6449606]
#
# ----------
# changelog:
# ----------
# 01. Had to change the way Q was calculated (factor N exponent inside)
#     because it was throwing an OverflowError otherwise (N**N for N > 144)
# 02. Published on wikipedia.org [13.07.16]
# 03. Check for N > 1 [13.10.20]
# 04. Check for nan in Lamda calculation [13.10.20]
# 05. Return 0 for negative x2 values [13.10.20]
#
# -----------
# references:
# -----------
# 01. B. Peirce, "Criterion for the rejection of doubtful observations," The 
#     Astronomical Journal, issue 45, vol. 2, no. 21, pp. 161--163, 1852.
#
# 02. B.A. Gould, "On Peirce's criterion for the rejection of doubtful 
#     observations, with tables for facilitating its application," The 
#     Astronomical Journal, issue 83, vol. 4, no. 11, pp. 81--87, 1855.
#
################################################################################
## LIBRARIES ###################################################################
################################################################################
import numpy
import scipy.special

################################################################################
## FUNCTIONS ###################################################################
################################################################################
def peirce_dev(N, n, m):
    """Calculate threshold error deviation for outliers.""" 
    # Methods based on Gould 1855
    # N :: total number of observations
    # n :: number of outliers to be removed
    # m :: number of unknowns in the relationship (e.g., regression parameters)
    #
    # Assign floats to input variables:
    N = float(N)
    n = float(n)
    m = float(m)
    #
    # Check number of observations:
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
        while ( abs(Rnew-Rold) > (N*2.0e-16) ):
            # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
            ldiv = Rnew**n
            if ldiv == 0:
                ldiv = 1.0e-6
            Lamda = ((Q**N)/(ldiv))**(1.0/(N-n))
            #
            # Calculate x-squared (straight-forward Gould's equation C):
            x2 = 1.0 + (N-m-n)/n * (1.0-Lamda**2.0)
            #
            # If x2 goes negative, return as 0.0:
            if x2 < 0:
                x2 = 0.0
                Rold = Rnew
            else:
                # Use x-squared to update R (Gould's equation D):
                Rold = Rnew
                Rnew = numpy.exp((x2-1)/2.0)*scipy.special.erfc(
                    numpy.sqrt(x2)/numpy.sqrt(2.0)
                    )
                #
    else:
        x2 = 0.0
    return x2
