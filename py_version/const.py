#!/usr/bin/python
#
# const.py
#
# LAST UPDATED: 2016-05-20
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# I. C. Prentice, T. W. Davis, X. M. P. Gilbert, B. D. Stocker, B. J. Evans,
# H. Wang, and T. F. Keenan, "The Global ecosystem in Space and Time (GePiSaT)
# Model of the Terrestrial Biosphere," (in progress).

###############################################################################
# IMPORT MODULES:
###############################################################################
import numpy

###############################################################################
# GLOBAL CONSTANTS:
###############################################################################
ke = 0.0167     # eccentricity for 2000 CE (Berger, 1978)
keps = 23.44    # obliquity for 2000 CE, degrees (Berger, 1978)
kfFEC = 2.04    # from flux to energy conversion, umol/J (Meek et al., 1984)
kG = 9.80665    # gravitational acceleration, m/s^2 (Allen, 1973)
kGsc = 1360.8   # solar constant, W/m^2 (Kopp & Lean, 2011)
kL = 0.0065     # temperature lapse rate, K/m (Allen, 1973)
kMa = 0.028963  # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
komega = 283.0  # longitude of perihelion for 2000 CE, degrees (Berger, 1978)
kPo = 101325    # standard atmosphere, Pa (Allen, 1973)
kR = 8.31447    # universal gas constant, J/mol/K (Moldover et al., 1988)
kTo = 288.15    # base temperature, K (Berberan-Santos et al., 1997)
pir = (numpy.pi/180.0)
