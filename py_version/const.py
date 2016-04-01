#!/usr/bin/python
#
# const.py
#
# LAST UPDATED: 2016-04-01
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
ke = 0.0167      # eccentricity for 2000 CE (Berger, 1978)
keps = 23.44     # obliquity for 2000 CE, degrees (Berger, 1978)
kfFEC = 2.04     # from flux to energy conversion, umol/J (Meek et al., 1984)
kGsc = 1360.8    # solar constant, W/m^2 (Kopp & Lean, 2011)
komega = 283.0   # longitude of perihelion for 2000 CE, degrees (Berger, 1978)
pir = (numpy.pi/180.0)
