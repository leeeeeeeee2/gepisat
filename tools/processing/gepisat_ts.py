#!/usr/bin/python
#
# timeseries.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-03-20 -- created
# 2014-10-22 -- last updated
#
# ------------
# description:
# ------------
# This script connects to the GePiSaT database, reads observations, and
# saves them to file as a time series for analysis.
#
# Note that most of the functions presented here are based on those written 
# in model.py
#
# ----------
# changelog:
# ----------
# 00. created [14.03.20]
# 01. updated to retrieve NERC proposal sites [14.06.26]
# 02. updated function doc [14.10.07]
# 03. updated outputs to include coordinates, precip and cld [14.10.07]
# 04. outputs all climatology [14.10.22]
# 05. changed output director to same as model.py [14.10.22]
# 06. changed default return value in get_data_point [14.10.22]
# ->  type 'None' to -9999
#
################################################################################
## IMPORT MODULES 
################################################################################
import datetime
import numpy
import os.path
import sys
import psycopg2

################################################################################
## FUNCTIONS 
################################################################################
def connect_sql():
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
            #database='flux_test',
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

def get_stations():
    """
    Name:     get_stations
    Input:    None.
    Output:   list, station names (results)
    Features: Returns a list of flux station names from GePiSaT database
    Depends:  connect_sql
    """
    # Define query:
    q = (
        "SELECT stationid "
        "FROM met_data "
        "WHERE dim=0 "
        "AND geom=%s "
        "ORDER BY stationid ASC;"
        )
    #
    params = ("point",)
    #
    # Connect to database and start cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and fetch results:
    cur.execute(q,params)
    results = []
    for record in cur:
        results.append(record[0])  # <- extract record from tuple
        #
    con.close()
    return results

def flux_to_grid(flux_station):
    """
    Name:     flux_to_grid
    Input:    string, station name (flux_station)
    Output:   int, grid station ID (grid_station)
    Features: Returns grid station ID based on the location of a given flux 
              tower
    Depends:  - get_lon_lat
              - grid_centroid
              - connect_sql
    """
    # Get lat and lon of flux tower:
    (fst_lon, fst_lat) = get_lon_lat(flux_station)
    #
    # Determine grid centroid lon and lat:
    (grd_lon, grd_lat) = grid_centroid(fst_lon, fst_lat)
    #
    # Get grid station name based on centroid coordinates:
    params = ("grid", grd_lon, grd_lat)
    q = (
        "SELECT met_data.stationid "
        "FROM met_data "
        "WHERE met_data.geom = %s "
        "AND met_data.lon = %s "
        "AND met_data.lat = %s;"
        )
    #
    # Connect to database:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return results:
    cur.execute(q, params)
    grid_station = cur.fetchone()[0]
    con.close()
    return grid_station

def site_to_grid(fst_lon, fst_lat):
    """
    Name:     site_to_grid
    Input:    - float, longitude (fst_lon)
              - float, latitude (fst_lat)
    Output:   int, grid station ID (grid_station)
    Features: Returns grid station ID based on a given location
    Depends:  - grid_centroid
              - connect_sql
    """
    # Determine grid centroid lon and lat:
    (grd_lon, grd_lat) = grid_centroid(fst_lon, fst_lat)
    #
    # Get grid station name based on centroid coordinates:
    params = ("grid", grd_lon, grd_lat)
    q = (
        "SELECT met_data.stationid "
        "FROM met_data "
        "WHERE met_data.geom = %s "
        "AND met_data.lon = %s "
        "AND met_data.lat = %s;"
        )
    #
    # Connect to database:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return results:
    cur.execute(q, params)
    grid_station = cur.fetchone()[0]
    con.close()
    return grid_station

def get_lon_lat(station):
    """
    Name:     get_lon_lat
    Input:    string, station name (station)
    Output:   tuple, lon-lat pair
              - float, longitude (my_lon)
              - float, latitude (my_lat)
    Features: Return longitude and latitude pair for a given station based on 
              the GePiSaT database meta-data table
    Depends:  connect_sql
    """
    # Query paramters:
    params = (station,)
    #
    # SQL query:
    q = (
        "SELECT met_data.lon, met_data.lat "
        "FROM met_data "
        "WHERE met_data.stationid = %s;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return results:
    cur.execute(q, params)
    my_lon, my_lat = cur.fetchone()
    con.close()
    return (my_lon, my_lat)

def grid_centroid(my_lon, my_lat):
    """
    Name:     grid_centroid
    Input:    - float, longitude (my_lon)
              - float, latitude (my_lat)
    Output:   tuple, longitude latitude pair (my_centroid)
    Features: Returns the nearest 0.5 deg. grid centroid per given coordinates
              based on the Euclidean distance to each of the four surrounding 
              grids; if any distances are equivalent, the pixel north and east
              is selected by default
    """
    # Create lists of regular latitude and longitude:
    grid_res = 0.5
    lat_min = -90 + 0.5*grid_res
    lon_min = -180 + 0.5*grid_res
    lat_dim = 360
    lon_dim = 720
    lats = [lat_min + y * grid_res for y in xrange(lat_dim)]
    lons = [lon_min + x * grid_res for x in xrange(lon_dim)]
    #
    # Find bounding longitude:
    centroid_lon = None
    if my_lon in lons:
        centroid_lon = my_lon
    else:
        lons.append(my_lon)
        lons.sort()
        lon_index = lons.index(my_lon)
        bb_lon_min = lons[lon_index-1]
        try:
            bb_lon_max = lons[lon_index+1]
        except IndexError:
            bb_lon_max = lons[-1] + grid_res
        #
    # Find bounding latitude:
    centroid_lat = None
    if my_lat in lats:
        centroid_lat = my_lat
    else:
        lats.append(my_lat)
        lats.sort()
        lat_index = lats.index(my_lat)
        bb_lat_min = lats[lat_index-1]
        try:
            bb_lat_max = lats[lat_index+1]
        except IndexError:
            bb_lat_max = lats[-1] + grid_res
        #
    # Determine nearest centroid:
    # NOTE: if dist_A equals dist_B, then centroid defaults positively 
    #       i.e., north / east
    if centroid_lon and centroid_lat:
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lon and not centroid_lat:
        # Calculate the distances between lat and bounding box:
        dist_A = bb_lat_max - my_lat
        dist_B = my_lat - bb_lat_min
        if dist_A > dist_B:
            centroid_lat = bb_lat_min
        else:
            centroid_lat = bb_lat_max
        my_centroid = (centroid_lon, centroid_lat)
    elif centroid_lat and not centroid_lon:
        # Calculate the distances between lon and bounding box:
        dist_A = bb_lon_max - my_lon
        dist_B = my_lon - bb_lon_min
        if dist_A > dist_B:
            centroid_lon = bb_lon_min
        else:
            centroid_lon = bb_lon_max
        my_centroid = (centroid_lon, centroid_lat)
    else:
        # Calculate distances between lat:lon and bounding box:
        # NOTE: if all distances are equal, defaults to NE grid
        dist_A = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_B = numpy.sqrt(
            (bb_lon_max - my_lon)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        dist_C = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (bb_lat_max - my_lat)**2.0
            )
        dist_D = numpy.sqrt(
            (my_lon - bb_lon_min)**2.0 + (my_lat - bb_lat_min)**2.0
            )
        min_dist = min([dist_A, dist_B, dist_C, dist_D])
        #
        # Determine centroid based on min distance:
        if dist_A == min_dist:
            my_centroid = (bb_lon_max, bb_lat_max)
        elif dist_B == min_dist:
            my_centroid = (bb_lon_max, bb_lat_min)
        elif dist_C == min_dist:
            my_centroid = (bb_lon_min, bb_lat_max)
        elif dist_D == min_dist:
            my_centroid = (bb_lon_min, bb_lat_min)
            #
    # Return nearest centroid:
    return my_centroid

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
        print "Error: cannot write to file: ", f
    else:
        OUT.close()

def get_msvidx(station, variable):
    """
    Name:     get_msvidx
    Input:    - string, station name (station)
              - string, variable name (variable)
    Output:   string, msvidx (result)
    Features: Returns the msvidx from the GePiSaT database based on the station
              and variable name
    Depends:  connect_sql
    """
    # Define query:
    q = (
        "SELECT var_list.msvidx "
        "FROM var_list "
        "WHERE var_list.stationid = %s "
        "AND var_list.varname = %s;"
        )
    #
    # SQL query parameters:
    params = (station, variable)
    #
    # Connect to database and star cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and fetch results:
    cur.execute(q, params)
    try:
        result = cur.fetchone()[0]
    except:
        print "Could not return an msvidx value for station", 
        print station, "and variable", variable
        result = ""
    finally:
        con.close()
        return result

def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime.date (dt0)
    Output:   datetime.date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32) 
    dt3 = dt2.replace(day=1)
    return dt3

def get_data_point(msvidx, time_point):
    """
    Name:     get_data_points
    Input:    - string, msvidx (msvidx)
              - datetime.date (time_point)
    Output:   float/numpy.ndarray (my_result)
    Features: Returns data point or array of data for a given msvidx (i.e., 
              station and variable) and time
    Depends:  connect_sql
    """
    # SQL query params:
    params = (msvidx, time_point)
    #
    # Define SQL query:
    q = (
        "SELECT data_set.data "
        "FROM data_set "
        "WHERE data_set.msvidx = %s "
        "AND data_set.datetime = %s;"
        )
    #
    # Connect to database and start a cursor:
    con = connect_sql()
    cur = con.cursor()
    #
    # Execute query and return result:
    cur.execute(q, params)
    if cur.rowcount == 1:
        my_result = cur.fetchone()[0]
    elif cur.rowcount > 1:
        my_result = numpy.array([])
        for record in cur:
            my_result = numpy.append(my_result, record[0])
    else:
        my_result = -9999
        print "No data found in function get_data_point"
        print "... msvidx =", msvidx
        print "... time = %s" % time_point
    return my_result

################################################################################
## DEFINITIONS
################################################################################
# Proposed leaf temperature sites:
stations = {
    'FI-Hyy':(61.8474, 24.2948),
    'UK-Wytham':(51.777526, -1.338931),
    'CN-Ail': (24.54, 101.02),
    'ES-CanBalasc': (41.41667, 2.0667),      # Sanchez-Conta et al.
    'CN-Yuanjiang': (23.692, 101.856),       # from Colin
    'CN-XTBG': (21.95, 101.02)
}

################################################################################
## MAIN 
################################################################################
header = (
    "Station,Month,St_Lon,St_Lat,Grid_Lon,Grid_Lat,"
    "Tair_C,VPD_kPa,SF,Pre_mm,fAPAR,alpha\n"
)
out_file = "out/timeseries.txt"
writeout(out_file, header)

# Get list of all flux station names:
stations = get_stations()
#
# For use with station dict:
#for station in stations.keys():
#
for station in stations:
    # For use with station dict:
    #(st_lat, st_lon) = stations[station]
    #
    # Get lon-lat of flux tower:
    (st_lon, st_lat) = get_lon_lat(station)
    (grid_lon, grid_lat) = grid_centroid(st_lon, st_lat)
    #
    # Get the half-degree station corresponding to flux tower's location:
    hdg_station = flux_to_grid(station)
    #
    # For use with station dict:
    #hdg_station = site_to_grid(st_lon, st_lat)
    #
    # Get variable indexes:
    fpar_msvidx = get_msvidx(hdg_station, 'FAPAR')
    vpd_msvidx = get_msvidx(hdg_station, 'VPD')
    alpha_msvidx = get_msvidx(hdg_station, 'alpha')
    tair_msvidx = get_msvidx(hdg_station, 'Tc')
    #
    ##### ADDED 2014-06-25 #####
    pre_msvidx = get_msvidx(hdg_station, 'Pre')
    cld_msvidx = get_msvidx(hdg_station, 'Cld')
    #
    # Starting and ending dates:
    sd = datetime.date(2002, 1, 1)
    ed = datetime.date(2007, 1, 1)
    #
    while sd < ed:
        fpar_month = get_data_point(fpar_msvidx, sd)
        vpd_month = get_data_point(vpd_msvidx, sd)
        alpha_month = get_data_point(alpha_msvidx, sd)
        tair_month = get_data_point(tair_msvidx, sd)
        #
        ##### ADDED 2014-10-07 #####
        pre_month = get_data_point(pre_msvidx, sd)
        cld_month = get_data_point(cld_msvidx, sd)
        sf_month = 1.0 - 0.01*cld_month
        #
        try:
            OUT = open(out_file, 'a')
            OUT.write(
                "%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (
                    station,     # tower name
                    sd,          # date
                    st_lon,      # station lon
                    st_lat,      # station lat
                    grid_lon,    # grid lon.
                    grid_lat,    # grid lat.
                    tair_month,  # air temp.
                    vpd_month,   # VPD
                    sf_month,    # sun frac
                    pre_month,   # precip.
                    fpar_month,  # fAPAR
                    alpha_month  # CPA
                    )
                )
        except IOError:
            print "Could not append to file", out_file
        else:
            OUT.close()
        #
        sd = add_one_month(sd)
    #