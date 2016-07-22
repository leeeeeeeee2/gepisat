#!/usr/bin/python
#
# db_util.py
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
# ---------
# citation:
# ---------
# I. C. Prentice, T. W. Davis, X. M. P. Gilbert, B. D. Stocker, B. J. Evans,
# H. Wang, and T. F. Keenan, "The Global ecosystem in Space and Time (GePiSaT)
# Model of the Terrestrial Biosphere," (in progress).
#
###############################################################################
# IMPORT MODULES:
###############################################################################
import logging
import os.path

import psycopg2


###############################################################################
# DEFINE FUNCTIONS:
###############################################################################
def connectSQL():
    """
    Name:     connectSQL
    Input:    None.
    Output:   psycopg2 connection (con)
    Features: Connects to postgreSQL database and returns connection handle
              or None type if connection fails
    """
    # Open credentials file:
    cred_file = "user.txt"
    my_user = "postgres"
    my_pass = "bitnami"
    if os.path.isfile(cred_file):
        try:
            f = open(cred_file, "r")
        except:
            logging.exception("Failed to read crendential file")
        else:
            logging.debug("Reading credential file")
            cred_data = f.readline()
            if cred_data:
                cred_data = cred_data.rstrip()
                my_user, my_pass = cred_data.split(',')

    # Initialize connection variable:
    con = None

    # Test database connection:
    try:
        con = psycopg2.connect(
            database='gepisat',
            #database='test',
            user=my_user,
            host='localhost',
            password=my_pass
        )
    except psycopg2.DatabaseError:
        logging.exception("Failed to connect to the database")
    except:
        logging.exception(
            "Encountered unknown error while connecting to database")
    else:
        logging.debug("Database connection created")
    finally:
        return con
