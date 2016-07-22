#!/usr/bin/python
#
# const.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2016-07-22
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
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

kc = 0.41       # Jmax cost coefficient
kphio = 0.093   # quantum efficiency (Long et. al., 1993)
gs25 = 4.220    # gamma-star at 25C, Pa (assuming 25 deg C & 98.716 kPa)
dha = 37830     # gamma-star activation energy, J/mol
kc25 = 39.97    # Michaelis-Menten Kc, Pa (assuming 25 deg C & 98.716 kPa)
ko25 = 2.748e4  # Michaelis-Menten, Pa (assuming 25 deg C & 98.716 kPa)
dhac = 79430    # Michaelis-Menten CO2 activation energy, J/mol
dhao = 36380    # Michaelis-Menten O2 activation, J/mol
kco = 2.09476e5   # atmos. CO2 concentration, ppm (US Standard Atmosphere)
