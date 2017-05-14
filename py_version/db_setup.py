#!/usr/bin/python
#
# db_setup.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-13
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
#
# ------------
# description:
# ------------
# This script connects to the flux project postgreSQL database and
# provides the following functions:
# 1. Create/drop database tables (met_data, var_list, data_set)
# 2. Read data from file and populates database tables
#
# NOTE: database name is hardcoded in connectSQL()
#
# ----------
# changelog:
# ----------
# - Added check for table existance [13.05.28]
# - Changed varchar(n) datatype to text [13.06.17]
#   http://www.depesz.com/2010/03/02/charx-vs-varcharx-vs-varchar-vs-text/
# - Added popdataset() function now that we have fluxtower data
# - Abstracted user name and password to file; defaults to user postgres
#   and pass bitnami (update as necessary) [13.08.27]
# - Housekeeping [13.08.27]
# - Removed alter table owner from table create queries [13.08.27]
#   * Think about adding this command in relation to "my_user"
# - Renamed script from "postgreSQL.py" to "db_setup.py" [13.08.27]
# - Changed import module os to os.path [13.08.27]
# - Added dependency check on createdataset() [13.09.12]
# - Updated getdata function [13.09.12]
#   * takes a filename (with path)
#   * returns t regarless of errors (None-type if error occurs)
# - Updated popmeta, popvar, and popdata functions [13.09.12]
#   * check for t
#   * else print errors
# - Gridded data station ids have 9 characters; update database definition of
#   stationid and msvidx [13.09.12]
#   * stationid varchar(12)
#   * msvidx varchar(15)
# - Added glob to modules list [13.09.12]
# - Creates reset_db, clean_db, and clean_table functions [13.10.01]
# - Created db_size() function [13.10.02]
# - Added get_var|data_files() functions [14.01.10]
# - Added BITNAMI directory structure to get data functions [14.02.05]
# - General housekeeping [14.09.02]
# - PEP8 style fixes [15.11.18]
# - added Python 3 support [16.01.15]
# - added logging [16.01.15]
#
# -------------
# database size
# -------------
# empty tables: 7.28 MB
# + HDG_Met-Data: 58.86 MB
# + WFDEI_Var-List: 97.15 MB
# + WFDEI_Data-List (2002-2006): 17660.64 MB
#
# -----
# todo:
# -----
# - getdata() function needs better file handling
#   * how to handle multiple input files?
# - Add database name (e.g. 'gepisat') to the user.txt file
#
###############################################################################
# IMPORT MODULES:
###############################################################################
import logging
import os.path
import sys

from database.utilities import find_files
from gepisat.db_util import connectSQL


###############################################################################
# DEFINE FUNCTIONS:
###############################################################################
def getversion():
    """
    Name:     get_version
    Input:    None.
    Output:   None.
    Features: Connects to postgreSQL database and prints the postgreSQL version
    Depends:  connectSQL
    """
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()

        # Print version and close connection:
        logging.debug("Checking database version")
        cur.execute("SELECT version()")
        ver = cur.fetchone()
        logging.info("%s", ver)
        con.close()


def getdata(myfile):
    """
    Name:     getdata
    Input:    str, input file with path (myfile)
    Output:   tuple (t)
    Features: Returns tuple of input file contents or None type if data file
              is unreadable
    """
    # Check that file exists:
    if os.path.isfile(myfile):
        # Open input file for reading:
        try:
            f = open(myfile, 'r')
        except IOError:
            logging.error("Failed to read data file %s", myfile)
            t = None
        except:
            logging.exception("Failed at reading data file %s", myfile)
            t = None
        else:
            # Read content and save to tuple:
            logging.debug("Reading data file %s", myfile)
            content = f.readlines()[1:]  # skip header line
            t = tuple(content)
        return t


def droptable(tname):
    """
    Name:     droptable
    Input:    string, table name (tname)
    Output:   None.
    Features: Drops the given table from the postgreSQL database if it exists
              and checks and drops any tables that depend on it
    Depends:  connectSQL
    """
    # Define query:
    q = "DROP TABLE IF EXISTS " + tname

    # Check dependencies:
    if (tname == "met_data"):
        # Check to see if var_list exists:
        isvarlist = existsSQL("var_list")
        if (isvarlist):
            logging.warning("Dropping dependent table 'var_list'")
            droptable("var_list")
    elif (tname == "var_list"):
        # Check to see if data_set exists:
        isdataset = existsSQL("data_set")
        if (isdataset):
            logging.warning("Dropping dependent table 'data_set'")
            droptable("data_set")

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()
        logging.debug("Executing SQL query")
        cur.execute(q)
        con.commit()
        con.close()


def createmetdata():
    """
    Name:     createmetdata
    Input:    None.
    Output:   None.
    Features: Creates the met_data postgreSQL table, dropping an existing table
              if it exists
    Depends:  - connectSQL
              - droptable
    """
    # Define create table query
    q = (
        "CREATE TABLE met_data ("
        "mapid character(3),"
        "map text,"
        "countryid character(2),"
        "country text,"
        "stationid varchar(12) NOT NULL,"
        "station text,"
        "lat real NOT NULL,"
        "lon real NOT NULL,"
        "ele real,"
        "classid character(3),"
        "class text,"
        "climateid varchar(4),"
        "climate text,"
        "data_years text,"
        "years_data integer NOT NULL DEFAULT 0,"
        "network text,"
        "url text,"
        "geom text,"
        "coor text,"
        "dim integer NOT NULL,"
        "res text,"
        "CONSTRAINT set_pk1 PRIMARY KEY (stationid)"
        ") WITH ("
        "OIDS=FALSE"
        "); "
        )

    # Drop met_data if it exists:
    droptable("met_data")

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()
        logging.debug("Executing SQL query")
        cur.execute(q)
        con.commit()
        con.close()


def createvarlist():
    """
    Name:     createvarlist
    Input:    None.
    Output:   None.
    Features: Creates the var_list postgreSQL table, dropping an existing table
              if it exists and checking that necessary dependency table exist
    Depends:  - connectSQL
              - droptable
              - existsSQL
    """
    # Define create table query:
    q = (
        "CREATE TABLE var_list ("
        "msvidx varchar(15) NOT NULL,"
        "stationid varchar(12) REFERENCES met_data(stationid),"
        "varid integer NOT NULL,"
        "varname text,"
        "varunit text,"
        "vartype character(4) NOT NULL,"
        "varcore integer,"
        "UNIQUE(stationid, varid),"
        "CONSTRAINT set_pk2 PRIMARY KEY (msvidx)"
        ") WITH ("
        "OIDS = FALSE"
        "); "
        )

    # Drop var_list if it exists:
    droptable("var_list")

    # Check to see if met_data exists and proceed if dependency exists:
    ismetdata = existsSQL("met_data")
    if ismetdata:
        logging.debug("Creating database connection")
        con = connectSQL()
        if con is not None:
            logging.debug("Creating connection cursor")
            cur = con.cursor()
            logging.debug("Executing SQL query")
            cur.execute(q)
            con.commit()
            con.close()
    else:
        logging.warning("Cannot proceed, missing dependency 'met_data' table")


def createdataset():
    """
    Name:     createdataset
    Input:    None.
    Output:   None.
    Features: Creates the data_set table, dropping an existing table if it
              exists and checking that necessary dependency tables exist
    Depends:  - droptable
              - existsSQL
              - connectSQL
    """
    # Define create table query:
    q = (
        "CREATE TABLE data_set ("
        "msvidx varchar(15) REFERENCES var_list(msvidx),"
        "stationid varchar(12) REFERENCES met_data(stationid),"
        "datetime timestamp,"
        "data float,"
        "UNIQUE(msvidx, datetime)"
        ") WITH ("
        "OIDS = FALSE"
        "); "
        )

    # Drop data_set if it exists:
    droptable("data_set")

    # Check to see if var_list table exists and proceed if dependency exists:
    isvarlist = existsSQL("var_list")
    if isvarlist:
        logging.debug("Creating database connection")
        con = connectSQL()
        if con is not None:
            logging.debug("Creating connection cursor")
            cur = con.cursor()
            logging.debug("Executing SQL query")
            cur.execute(q)
            con.commit()
            con.close()


def popmetdata(filename):
    """
    Name:     popmetdata
    Input:    string, input file name (filename)
    Output:   None.
    Features: Populates the met_data table
    Depends:  - getdata
              - existsSQL
              - connectSQL
    """
    # Get table data
    t = getdata(filename)
    if t:
        # Make sure the number of columns is right before processing
        ncols = len(t[0].rstrip().split(','))
        if ncols != 21:
            logging.error(
                "Found %d cols, met_data table must have 21 cols!", ncols)
            sys.exit(1)

        # Check that table exists:
        ismetdata = existsSQL("met_data")
        if (ismetdata):
            logging.debug("Creating database connection")
            con = connectSQL()
            if con is not None:
                logging.debug("Creating connection cursor")
                cur = con.cursor()

                # Add data to table (met_data has 21 columns):
                logging.debug("Reading table values from file %s...", filename)
                for r in t:
                    l = (
                        "('%s','%s','%s','%s','%s','%s','%s','%s',"
                        "'%s','%s','%s','%s','%s','%s','%s','%s',"
                        "'%s','%s','%s','%s','%s')"
                        ) % tuple(r.rstrip().split(','))
                    cur.execute("INSERT INTO met_data VALUES " + l)
                logging.debug("...complete")

                # Commit changes to database:
                con.commit()
                con.close()
        else:
            logging.warning("Table 'met_data' does not exist!")
    else:
        logging.warning("No data found in file %s", filename)


def popvarlist(filename):
    """
    Name:     popvarlist
    Input:    string, input file name (filename)
    Output:   None.
    Features: Populates the var_list table
    Depends:  - getdata
              - existsSQL
              - connectSQL
    """
    # Get table data:
    t = getdata(filename)
    if t:
        # Make sure number of cols is right before processing:
        ncols = len(t[0].rstrip().split(','))
        if ncols != 7:
            logging.error(
                "Found %d cols, var_list table must have 7 cols!", ncols)
            sys.exit(1)

        # Check that table exists:
        isvarlist = existsSQL("var_list")
        if (isvarlist):
            logging.debug("Creating database connection")
            con = connectSQL()
            if con is not None:
                logging.debug("Creating connection cursor")
                cur = con.cursor()

                # Add data to table (var_list has 7 columns):
                logging.debug("Reading table values from file %s...", filename)
                for r in t:
                    l = (
                        "('%s','%s','%s','%s','%s','%s','%s')"
                        ) % tuple(r.rstrip().split(','))
                    cur.execute("INSERT INTO var_list VALUES " + l)
                logging.debug("...complete")

                # Commit changes to database and close:
                con.commit()
                con.close()
        else:
            logging.warning("Table 'var_list' does not exist!")
    else:
        logging.warning("No data found in file %s", filename)


def popdataset(filename):
    """
    Name:     popdataset
    Input:    string, input file name (filename)
    Output:   None.
    Features: Populates the data_set table
    Depends:  - getdata
              - existsSQL
              - connectSQL
    """
    # Get table data:
    t = getdata(filename)
    if t:
        # Make certain the number of columns is correct:
        ncols = len(t[0].rstrip().split(','))
        if ncols != 4:
            logging.error(
                "Found %d cols, data_set table must have 4 cols!", ncols)
            sys.exit(1)

        # Check that the table exists:
        isdataset = existsSQL("data_set")
        if (isdataset):
            logging.debug("Creating database connection")
            con = connectSQL()
            if con is not None:
                logging.debug("Creating connection cursor")
                cur = con.cursor()

                # Add data to table:
                logging.debug("Reading table values from file %s...", filename)
                for r in t:
                    l = "('%s','%s','%s','%s')" % tuple(r.rstrip().split(','))
                    cur.execute("INSERT INTO data_set VALUES " + l)
                logging.debug("...complete")

                # Commit changes to database and close:
                con.commit()
                con.close()
        else:
            logging.warning("Table data_set does not exist!")
    else:
        logging.warning("No data found in file %s", filename)


def existsSQL(tname):
    """
    Name:     existsSQL
    Input:    string, table name (tname)
    Output:   boolean
    Features: Returns boolean if the table exists in the postgreSQL database
    Depends:  connectSQL
    """
    # Define a query:
    q = (
        "SELECT * FROM information_schema.tables "
        "WHERE table_name='%s'"
        ) % tname

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()
        logging.debug("Executing SQL query")
        cur.execute(q)
        myresult = bool(cur.rowcount)
        con.close()

        return myresult


def reset_db():
    """
    Name:     reset_db
    Input:    None.
    Output:   None.
    Features: Deletes existing and creates the met_data, var_list, and data_set
              postgreSQL tables
    Depends:  - existsSQL
              - droptable
              - createmetadata
              - createvarlist
              - createdataset
    """
    # Check which tables exist:
    ismd = existsSQL('met_data')
    isvl = existsSQL('var_list')
    isds = existsSQL('data_set')

    # Drop those that do:
    if isds:
        logging.warning("Dropping 'data_set'...")
        droptable('data_set')
    if isvl:
        logging.warning("Dropping 'var_list'...")
        droptable('var_list')
    if ismd:
        logging.warning("Dropping 'met_data'...")
        droptable('met_data')

    # Recreate all three tables:
    logging.info("Creating 'met_data'...")
    createmetdata()

    logging.info("Creating 'var_list'...")
    createvarlist()

    logging.info("Creating 'data_set'...")
    createdataset()


def clean_db():
    """
    Name:     clean_db
    Input:    None.
    Output:   None.
    Features: Deletes all rows (data) from met_data, var_list, and data_set
              tables if they exist (does not drop tables)
    Depends:  - existsSQL
              - connectSQL
    """
    # Check which tables exist:
    ismd = existsSQL('met_data')
    isvl = existsSQL('var_list')
    isds = existsSQL('data_set')

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()

        # Clean tables:
        if isds:
            logging.debug("Cleaning 'data_set' table")
            q = "DELETE FROM %s;" % ("data_set")
            cur.execute(q)
        if isvl:
            logging.debug("Cleaning 'var_list' table")
            q = "DELETE FROM %s;" % ("var_list")
            cur.execute(q)
        if ismd:
            logging.debug("Cleaning 'met_data' table")
            q = "DELETE FROM %s;" % ("met_data")
            cur.execute(q)

        # Commit changes and close:
        con.commit()
        con.close()


def clean_table(tname):
    """
    Name:     clean_table
    Input:    string, table name (tname)
    Output:   None.
    Features: Deletes rows (data) from a given table (does not drop table)
    Depends:  - existsSQL
              - connectSQL
    """
    # Define query:
    q = (
        "DELETE FROM %s;"
        ) % tname

    # Check which tables exist:
    istable = existsSQL(tname)

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()

        # Clean table:
        if istable:
            logging.debug("Cleaning '%s' table", tname)
            cur.execute(q)
        else:
            logging.warning("Table '%s' does not exist", tname)

        # Commit changes and close:
        con.commit()
        con.close()


def db_size(db_name='gepisat'):
    """
    Name:     db_size
    Input:    [optional] str, database name (db_name)
    Output:   int, database size (in bytes)
    Features: Returns the disk size of a postgreSQL database (bytes)
    Depends:  connectSQL
    Note:     Hard-coded database name
    """
    # Define db name as query parameter:
    params = (db_name,)

    # Define query:
    q = (
        "SELECT pg_database_size(%s) "
        "As fulldbsize;"
        )

    logging.debug("Creating database connection")
    con = connectSQL()
    if con is not None:
        logging.debug("Creating connection cursor")
        cur = con.cursor()

        # Execute query, fetch results:
        logging.debug("Executing SQL query")
        cur.execute(q, params)
        my_result = cur.fetchone()
        con.close()

        try:
            my_size = int(my_result[0])
        except ValueError:
            logging.exception("Return value %s not an int", my_result[0])
        except TypeError:
            try:
                my_size = int(my_result)
            except:
                logging.exception("Failed to convert result %s", my_result)
        except:
            logging.exception("Failed to convert result %s", my_result)

        return my_size


def get_var_files(base_path):
    """
    Convenience function for getting GePiSaT database files.

    Name:     get_var_files
    Input:    str, base path (base_path)
    Output:   list, file names (files_list)
    Features: Returns a list of var_list file names (with path);
              NOTE: you need to update these paths!
    Depends:  find_files
    """
    path_list = [
        base_path,
        os.path.join(base_path, "cru", "flux_stations"),
        os.path.join(base_path, "flux"),
        os.path.join(base_path, "modis", "flux_stations"),
        os.path.join(base_path, "noaa", "flux_stations"),
        os.path.join(base_path, "watch", "flux_stations"),
    ]

    files_list = []
    for path in path_list:
        logging.debug("Checking for files in %s", os.path.dirname(path))
        tmp_list = find_files(path, "*Var-List*csv")
        files_list += tmp_list

    return files_list


def get_data_files(base_path):
    """
    Convenience function for getting GePiSaT database files.

    Name:     get_data_files
    Input:    str, base path (base_path)
    Output:   list, file names (files_list)
    Features: Returns a list of data_set file names (with paths)
    Depends:  find_files
    Note:     Hard-coded base_path
    """
    path_list = [
        base_path,
        os.path.join(base_path, "cru", "flux_stations"),
        os.path.join(base_path, "flux"),
        os.path.join(base_path, "modis", "flux_stations"),
        os.path.join(base_path, "noaa", "flux_stations"),
        os.path.join(base_path, "watch", "flux_stations"),
    ]

    files_list = []
    for path in path_list:
        logging.debug("Checking for files in %s", os.path.dirname(path))
        tmp_list = find_files(path, "*Data-Set*csv")
        files_list += tmp_list

    return files_list

###############################################################################
# MAIN PROGRAM:
###############################################################################
if __name__ == "__main__":
    # Create a root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_handler = logging.StreamHandler()
    rec_format = "%(asctime)s:%(levelname)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt='%H%M%S')
    root_handler.setFormatter(formatter)
    root_logger.addHandler(root_handler)

    # Test database connection:
    # getversion()

    # Reset or clean database:
    reset_db()
    # clean_db()
    # clean_table("data_set")

    # Debug: check if tables exist:
    root_logger.debug("met: %s", existsSQL('met_data'))
    root_logger.debug("var: %s", existsSQL('var_list'))
    root_logger.debug("dat: %s", existsSQL('data_set'))

    # Check database size:
    my_dbsize = db_size()
    root_logger.info("Database starting size: %0.3f MB" % (1e-6*my_dbsize))

    # Define directories:
    home_dir = os.path.expanduser("~")
    data_dir = os.path.join(home_dir, "Data", "psql")

    # for testing file searches
    run_md = True
    run_var = True
    run_ds = True
    to_test = False

    # Create and populate met_data table:
    if run_md:
        root_logger.info("Processing 'met_data' table...")
        md_files = find_files(data_dir, "*Met-Data*csv")
        for x in md_files:
            if to_test:
                root_logger.info("found %s", os.path.basename(x))
            else:
                popmetdata(x)
                root_logger.info("added %s (%0.3f MB)",
                                 os.path.basename(x), 1e-6*db_size())
        root_logger.info("... 'met_data' complete")

    if run_var:
        # Create and populate var_list table :
        # NOTE: some cols depend on met_data, so create/pop it first
        root_logger.info("Processing 'var_list' table...")
        vl_files = get_var_files(data_dir)
        for y in vl_files:
            if to_test:
                root_logger.info("found %s", os.path.basename(y))
            else:
                popvarlist(y)
                root_logger.info("added %s (%0.3f MB)",
                                 os.path.basename(y), 1e-6*db_size())
        root_logger.info("... 'var_list' complete")

    if run_ds:
        # Create and populate data_set table
        root_logger.info("Processing 'data_set' table...")
        ds_files = get_data_files(data_dir)
        for z in sorted(ds_files):
            if to_test:
                root_logger.info("found %s", os.path.basename(z))
            else:
                popdataset(z)
                root_logger.info("added %s (%0.3f MB)",
                                 os.path.basename(z), 1e-6*db_size())
        root_logger.info("... 'data_set' complete")
