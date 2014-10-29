#!/usr/bin/python
#
# db_setup.py (BITNAMI version)
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-05-15 -- created (as postgreSQL.py)
# 2014-10-22 -- last updated
#
# ------------
# description:
# ------------
# This script connects to the flux project postgreSQL database and
# provides the following functions:
# 1. Create/drop database tables (met_data, var_list, data_set)
# 2. Read data from file and populates database tables
#
# NOTE 1: getdata() function requires hardcoding data directory
#
# NOTE 2: database name is hardcoded in connectSQL() and db_size()
#
# ----------
# changelog:
# ----------
# 01. Added check for table existance [13.05.28]
# 02. Changed varchar(n) datatype to text [13.06.17]
#     http://www.depesz.com/2010/03/02/charx-vs-varcharx-vs-varchar-vs-text/ 
# 03. Added popdataset() function now that we have fluxtower data
# 04. Abstracted user name and password to file; defaults to user postgres
#     and pass bitnami (update as necessary) [13.08.27]
# 05. Housekeeping [13.08.27]
# 06. Removed alter table owner from table create queries [13.08.27]
#     * Think about adding this command in relation to "my_user"
# 07. Renamed script from "postgreSQL.py" to "db_setup.py" [13.08.27]
# 08. Changed import module os to os.path [13.08.27]
# 09. Added dependency check on createdataset() [13.09.12]
# 10. Updated getdata function [13.09.12]
#     * takes a filename (with path)
#     * returns t regarless of errors (None-type if error occurs)
# 11. Updated popmeta, popvar, and popdata functions [13.09.12]
#     * check for t
#     * else print errors
# 12. Gridded data station ids have 9 characters; update database definition of 
#     stationid and msvidx [13.09.12]
#     * stationid varchar(12)
#     * msvidx varchar(15)
# 13. Added glob to modules list [13.09.12]
# 14. Creates reset_db, clean_db, and clean_table functions [13.10.01]
# 15. Created db_size() function [13.10.02]
# 16. Added get_var|data_files() functions [14.01.10]
# 17. Added BITNAMI directory structure to get data functions [14.02.05]
# 18. General housekeeping [14.09.02]
#
# -----
# todo:
# -----
# 01. getdata() function needs better file handling
#     * how to handle multiple input files?
# 02. Add database name (e.g. 'gepisat') to the user.txt file
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import glob
import os.path
import sys
import psycopg2

###############################################################################
## DEFINE FUNCTIONS ###########################################################
###############################################################################
def connectSQL():
    """
    Name:     connectSQL
    Input:    None.
    Output:   psycopg2 connection (con)
    Features: Connects to postgreSQL database and returns connection handle
    """
    # Open credentials file:
    cred_file = "./user.txt"
    if os.path.isfile(cred_file):
        f = open(cred_file, "r")
        cred_data = f.readline()
        if cred_data:
            cred_data = cred_data.rstrip()
            my_user, my_pass = cred_data.split(',')
    else:
        my_user, my_pass = ('postgres', 'bitnami')
    #
    # Initialize connection variable:
    con = None
    #
    # Test database connection:
    try:
        con = psycopg2.connect(
            database='gepisat', 
            #database='test',
            user=my_user, 
            host='localhost', 
            password=my_pass
            )
    #
    except psycopg2.DatabaseError, e:
        print 'Error %s' % e
        sys.exit(1)
    #
    return con

def getversion():
    """
    Name:     get_version
    Input:    None.
    Output:   None.
    Features: Connects to postgreSQL database and prints the postgreSQL version
    Depends:  connectSQL
    """
    # Connect to the database and create cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Execute query:
    cur.execute("SELECT version()")
    #
    # Print version:
    ver = cur.fetchone()
    print ver
    #
    # Close connection
    con.close()

def getdata(myfile):
    """
    Name:     getdata
    Input:    string, input file with path (myfile)
    Output:   tuple (t)
    Features: Returns tuple of input file contents
    """
    # Check that file exists:
    if os.path.isfile(myfile):
        # Open input file for reading:
        try:
            f = open(myfile, 'r')
        except IOError as e:
            print "Error reading file!\n" + e
            t = None
        except:
            print "Unexpected error: ", sys.exc_info()[0]
            t = None
        else:
            # Read content and save to tuple:
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
    # Define querry:
    q = "DROP TABLE IF EXISTS " + tname
    # 
    # Check dependencies:
    if (tname == "met_data"):
        # Check to see if var_list exists:
        isvarlist = existsSQL("var_list")
        #
        # Warn and drop dependent table:
        if (isvarlist):
            print "Warning, dropping dependent table var_list"				
            droptable("var_list")
    elif (tname == "var_list"):
        # Check to see if data_set exists:
        isdataset = existsSQL("data_set")
        #
        # Warn and drop dependent table:
        if (isdataset):
            print "Warning, dropping dependent table data_set"
            droptable("data_set")
    #
    # Connect to database and create cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Execute querry:
    cur.execute(q)
    #
    # Commit changes to database and close:
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
    # Define create table querry
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
    #
    # Drop met_data if it exists:
    droptable("met_data")
    #
    # Connect to database:
    con = connectSQL()
    cur = con.cursor()
    #
    # Create met_data table:
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
    # Define create table querry:
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
    #
    # Drop var_list if it exists:
    droptable("var_list")
    #
    # Check to see if met_data exists
    ismetdata = existsSQL("met_data")
    #
    # Proceed if dependency exists:
    if ismetdata:
        # Connect to database:
        con = connectSQL()
        cur = con.cursor()
        # Create var_list table:
        cur.execute(q)
        con.commit()
        con.close()
    else:
        print "Cannot proceed, missing dependency 'met_data' table"

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
    # Define create table querry:
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
    #
    # Drop data_set if it exists:
    droptable("data_set")
    #
    # Check to see if var_list table exists:
    isvarlist = existsSQL("var_list")
    #
    # Proceed if dependency exists:
    if isvarlist:
        # Connect to database:
        con = connectSQL()
        cur = con.cursor()
        # Create data_set table:
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
        if (ncols != 21):
            print "Input for met_data table must have 21 columns!"
            sys.exit(1)
            #
        # Check that table exists:
        ismetdata = existsSQL("met_data")
        if (ismetdata):
            # Connect to database and create cursor:
            con = connectSQL()
            cur = con.cursor()
            #
            # Add data to table (met_data has 21 columns):
            for r in t:
                l = (
                    "('%s','%s','%s','%s','%s','%s','%s','%s',"
                    "'%s','%s','%s','%s','%s','%s','%s','%s',"
                    "'%s','%s','%s','%s','%s')"
                    ) % tuple(r.rstrip().split(','))
                cur.execute("INSERT INTO met_data VALUES " + l)
                #
            # Commit changes to database:
            con.commit()
            con.close()
        else:
            print "Table met_data does not exist!"
    else:
        print "No data read from file", filename

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
        if (ncols != 7):
            print "Input for var_list table must have seven columns!"
            sys.exit(1)
            #
        # Check that table exists:
        isvarlist = existsSQL("var_list")
        if (isvarlist):
            # Connect to the database and create a cursor:
            con = connectSQL()
            cur = con.cursor()
            #
            # Add data to table (var_list has 7 columns):
            for r in t:
                l = (
                    "('%s','%s','%s','%s','%s','%s','%s')"
                    ) % tuple(r.rstrip().split(','))
                cur.execute("INSERT INTO var_list VALUES " + l)
                #
            # Commit changes to database and close:
            con.commit()
            con.close()
        else:
            print "Table var_list does not exist!"
    else:
        print "No data read from file", filename

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
        if (ncols != 4):
            print "Input for data_set table must have four columns!"
            sys.exit(1)
            #
        # Check that the table exists:
        isdataset = existsSQL("data_set")
        if (isdataset):
            # Connect to database and create a cursor:
            con = connectSQL()
            cur = con.cursor()
            #
            # Add data to table:
            for r in t:
                l = "('%s','%s','%s','%s')" % tuple(r.rstrip().split(','))
                cur.execute("INSERT INTO data_set VALUES " + l)
                #
            # Commit changes to database and close:
            con.commit()
            con.close()
        else:
            print "Table data_set does not exist!"
    else:
        print "No data read from file", filename

def existsSQL(tname):
    """
    Name:     existsSQL
    Input:    string, table name (tname)
    Output:   boolean
    Features: Returns boolean if the table exists in the postgreSQL database
    Depends:  connectSQL
    """
    # Define a querry:
    q = (
        "SELECT * FROM information_schema.tables "
        "WHERE table_name='%s'"
        ) % tname
    #
    # Connect to database and start a cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Execute querry:
    cur.execute(q)
    #
    # Fetch results:
    myresult = bool(cur.rowcount)
    #
    # Close connection
    con.close()
    #
    # Return boolean:
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
    #
    # Drop those that do:
    if isds:
        print "Dropping 'data_set'..."
        droptable('data_set')
    if isvl:
        print "Dropping 'var_list'..."
        droptable('var_list')
    if ismd:
        print "Dropping 'met_data'..."
        droptable('met_data')
    #
    # Recreate all three tables:
    print "Creating 'met_data'..."
    createmetdata()
    print "Creating 'var_list'..."
    createvarlist()
    print "Creating 'data_set'..."
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
    # Initialize query:
    q = ""
    #
    # Check which tables exist:
    ismd = existsSQL('met_data')
    isvl = existsSQL('var_list')
    isds = existsSQL('data_set')
    #
    # Connect to database and start a cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Clean tables:
    if isds:
        q = "DELETE FROM %s;" % ("data_set") 
        cur.execute(q)
    if isvl:
        q = "DELETE FROM %s;" % ("var_list")
        cur.execute(q)
    if ismd:
        q = "DELETE FROM %s;" % ("met_data")
        cur.execute(q)
    #
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
    #
    # Check which tables exist:
    istable = existsSQL(tname)
    #
    # Connect to database and start a cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Clean table:
    if istable:
        cur.execute(q)
    #
    # Commit changes and close:
    con.commit()
    con.close()

def db_size():
    """
    Name:     db_size
    Input:    None.
    Output:   string, database size
    Features: Returns the disk size of a postgreSQL database (bytes)
    Depends:  connectSQL
    Note:     Hard-coded database name
    """
    # Define db name as query parameter:
    params = ('gepisat',)
    #params = ('test',)
    #
    # Define query:
    q = (
        "SELECT pg_database_size(%s) "
        "As fulldbsize;"
        )
    #
    # Connect to db and start a cursor:
    con = connectSQL()
    cur = con.cursor()
    #
    # Execute query, fetch results:
    cur.execute(q, params)
    my_result = cur.fetchone()
    con.close()
    #
    return my_result[0]

def get_var_files():
    """
    Name:     get_var_files
    Input:    None.
    Output:   list, file names (files_list)
    Features: Returns a list of var_list file names (with path)
    Note:     Hard-coded base_path
    """
    # Merl:
    base_path = "/usr/local/share/database/"
    # BITNAMI:
    #base_path = "/database/files/var_list/"
    #
    path_list = [
        #base_path + "cru/vpd/",
        #base_path + "cru/tc/",
        #base_path + "cru/pre/"
        base_path + "cru/cld/",
        #base_path + "cru/elv/",
        #base_path + "noaa/co2/",
        #base_path + "fluxtowers/station_data/",
        #base_path + "glas/canopy_height/",
        #base_path + "modis/evi/",
        #base_path + "watch/swdown/",
        #base_path + "alpha/"
    ]
    #
    files_list = []
    #
    for path in path_list:
        tmp_list = glob.glob(path + "*Var-List*")
        files_list = files_list + tmp_list
        #
    return files_list

def get_data_files():
    """
    Name:     get_data_files
    Input:    None.
    Output:   list, file names (files_list)
    Features: Returns a list of data_set file names (with paths)
    Note:     Hard-coded base_path
    """
    # Merl:
    base_path = "/usr/local/share/database/"
    # BITNAMI:
    #base_path = "/database/files/data_set/"
    #
    path_list = [
        #base_path + "cru/vpd/"
        #base_path + "cru/tc/",
        #base_path + "cru/pre/"
        #base_path + "cru/cld/",
        #base_path + "cru/elv/",
        #base_path + "noaa/co2/",
        #base_path + "fluxtowers/station_data/",
        #base_path + "glas/canopy_height/",
        #base_path + "modis/evi/aqua/",
        #base_path + "modis/evi/terra/",
        #base_path + "watch/swdown/",
        base_path + "alpha/",
    ]
    #
    files_list = []
    #
    for path in path_list:
        tmp_list = glob.glob(path + "*Data-Set*")
        files_list = files_list + tmp_list
        #
    return files_list

###############################################################################
## MAIN PROGRAM ###############################################################
###############################################################################
# Test database connection:
#getversion()

# Reset or clean database:
#reset_db()
#clean_db()
#clean_table("data_set")

# Check database size:
#my_dbsize = db_size()
#print "Start,%s" % my_dbsize

# Create and populate met_data table:
#met_dir = "/database/files/met_data/point/"
#print "Processing 'met_data'"
#md_files = glob.glob(met_dir + "*Met-Data*")
#for x in sorted(md_files):
#    #print os.path.basename(x)
#    popmetdata(x)
#    print "%s,%s" % (os.path.basename(x), db_size())

# Create and populate var_list table :
# * NOTE: some cols depend on met_data, so create/pop it first
#print "Processing 'var_list'"
#vl_files = get_var_files()
#for y in sorted(vl_files):
#    #print " ...", os.path.basename(y) 
#    popvarlist(y)
#    print "%s,%s" % (os.path.basename(y), db_size())
#
# Create and populate data_set table
#print "Processing 'data_set'"
ds_files = get_data_files()
for z in sorted(ds_files):
    #print "%s" % os.path.basename(z)
    popdataset(z)
    print "%s,%s" % (os.path.basename(z), db_size())
#
#print "met: %s" % existsSQL('met_data')
#print "var: %s" % existsSQL('var_list')
#print "dat: %s" % existsSQL('data_set')
