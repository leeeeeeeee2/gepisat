#!/usr/bin/python
#
# utilities.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-04-21
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2017 Prentice Lab
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
# ---------
# citation:
# ---------
# Davis, T.W., B.D. Stocker, X.M.P. Gilbert, T.F. Keenan, H. Wang, B.J. Evans,
# and I.C. Prentice. The Global ecosystem Production in Space and Time
# (GePiSaT) Model of the terrestrial biosphere: Part 1 - Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.

###############################################################################
# IMPORT MODULES
###############################################################################
import glob
import os
import logging


###############################################################################
# FUNCTIONS
###############################################################################
def find_files(my_dir, my_pattern):
    """
    Name:     find_files
    Inputs:   - str, directory path (my_dir)
              - str, file name search pattern (my_pattern)
    Outputs:  list, file paths
    Features: Returns a sorted list of files found at a given directory with
              file names that match a given pattern
    """
    my_files = []
    if os.path.isdir(my_dir):
        s_str = os.path.join(my_dir, my_pattern)
        my_files = glob.glob(s_str)
    if len(my_files) > 0:
        my_files = sorted(my_files)
    return my_files


def writeout(f, d):
    """
    Name:     writeout
    Input:    - string, file name with path (t)
              - string, data to be written to file (d)
    Output:   None
    Features: Writes new/overwrites existing file with data string
    """
    try:
        OUT = open(f, 'w')
        OUT.write(d)
    except IOError:
        logging.error("Error: cannot write to file: %s", f)
    else:
        OUT.close()
