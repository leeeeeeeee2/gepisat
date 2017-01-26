#!/usr/bin/python
#
# resources.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-01-25
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
import errno
import logging
import os


##############################################################################
# FUNCTIONS
##############################################################################
def mkdir_p(path):
    """
    Name:     mkdir_p
    Inputs:   str, directory path (path)
    Outputs:  None.
    Features: Makes directories, including intermediate directories as
              required (i.e., directory tree)
    Ref:      tzot (2009) "mkdir -p functionality in python,"
              StackOverflow, Online (Accessed 2017-01-25):
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    try:
        logging.debug('building directory tree')
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            logging.debug('directory exists')
            pass
        else:
            logging.error('failed to create directory!')
            raise
