#!/usr/bin/python
#
# crudata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-05-05
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
import datetime
import logging
import os

import numpy
from scipy.io import netcdf

from .utilities import find_files
from .utilities import get_station_latlon
from .utilities import writeout
from .var import VAR


###############################################################################
# FUNCTIONS
###############################################################################
def add_one_month(dt0):
    """
    Name:     add_one_month
    Input:    datetime date (dt0)
    Output:   datetime date (dt3)
    Features: Adds one month to datetime
    Ref:      A. Balogh (2010), ActiveState Code
              http://code.activestate.com/recipes/577274-subtract-or-add-a-
              month-to-a-datetimedate-or-datet/
    """
    dt1 = dt0.replace(day=1)
    dt2 = dt1 + datetime.timedelta(days=32)
    dt3 = dt2.replace(day=1)
    return dt3


def calculate_vpd(vap, tm, tx=None):
    """
    Name:     calculate_vpd
    Input:    - numpy nd.array, mean monthly vapor pressure, hPa (vap)
              - numpy nd.array, mean monthly min daily air temp
                or mean monthly mean daily air temp, deg C (tm)
              - [optional] numpy nd.array, mean monthly max daily air temp,
                deg C (tx)
                * USE ONLY IF USING MIN DAILY AIR TEMP
    Output:   numpy nd.array, mean monthly vapor pressure deficit, kPa (vpd)
    Features: Returns an array of mean monthly vapor pressure deficit based
              on CRU TS mean monthly vapor pressure and air temperatures;
              allows for either the max/min daily air temperatures or mean
              daily air temperature. For mean daily air temperature,
              leave 'tx' as NoneType.
    Ref:      Eq. 5.1, Abtew and Meleese (2013), Ch. 5 Vapor Pressure
              Calculation Methods, in Evaporation and Evapotranspiration:
              Measurements and Estimations, Springer, London.
                vpd = 0.611*exp[ (17.27 tc)/(tc + 237.3) ] - ea
                where:
                    tc = average daily air temperature, deg C
                    ea = actual vapor pressure, kPa

    TODO:  _ create clip and error field arrays
           _ calculate vpd in a single line computation
    """
    # Initialize array:
    # NOTE: maintains missing value
    vpd = -9999.0*numpy.ones(shape=(360, 720))

    # Iterate through each data point:
    lat, lon = tm.shape

    for y in range(lat):
        for x in range(lon):
            ea = vap[y, x]
            if tx is not None:
                tmin = tm[y, x]
                tmax = tx[y, x]
                tc = 0.5*(tmin + tmax)
            else:
                tc = tm[y, x]

            if tc < 1.e6 and ea < 1.e6:
                vpd[y, x] = (
                    0.611*numpy.exp((17.27*tc)/(tc + 237.3)) - 0.10*ea)
    return vpd


def get_cru_file(my_dir, my_date, my_voi):
    """
    Name:     get_cru_file
    Inputs:   - str, directory to CRU netcdf file (my_dir)
              - datetime.date, current month datetime object (my_date)
              - str, variable of interest (my_voi)
    Outputs:  str,
    Features: Searches for decadal, then non-specific, CRU-TS file for a given
              time stamp and variable of interest
    """
    # Year for processing:
    my_year = my_date.year

    # Try decadal files first, then generic file search:
    if my_year > 1990 and my_year < 2001:
        s_val = "*1991.2000.%s.dat.nc" % (my_voi)
    elif my_year > 2000 and my_year < 2011:
        s_val = "*2001.2010.%s.dat.nc" % (my_voi)
    elif my_year > 2010 and my_year < 2016:
        s_val = "*2011.2015.%s.dat.nc" % (my_voi)
    else:
        s_val = "*%s.dat.nc" % (my_voi)

    # Search directory for netCDF file:
    my_files = find_files(my_dir, s_val)
    num_files = len(my_files)

    if num_files == 0:
        raise IndexError(
            "Failed to find CRU files for year %s and variable %s" % (
                my_year, my_voi))
    elif num_files == 1:
        my_file = my_files[0]
    else:
        raise IndexError(
            "Found multiple CRU files for year %s and variable %s" % (
                my_year, my_voi))

    return my_file


def get_cru_lat(my_dir, my_date, my_voi):
    """
    Name:     get_cru_lat
    Inputs:   - str, directory to CRU TS netcdf file (my_dir)
              - datetime.date, current month datetime object (my_date)
              - str, variable of interest, for CRU file only (my_voi)
    Outputs:  numpy.ndarray, latitudes (degrees)
    Features: Returns an array of latitudes from a CRU TS netCDF file
    """
    try:
        my_file = get_cru_file(my_dir, my_date, my_voi)
    except:
        raise
    else:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")

        # Save data for variables of interest:
        # NOTE: for CRU TS 4.00:
        #       variables: 'lat', 'lon', 'time', voi, 'stn'
        #       where
        #           lat (latitude): -89.75 -- 89.75
        #           > units: degrees_north
        #           lon (longitude): -179.75 -- 179.75
        #           > units: degrees_east

        # Create a copy of the latitude data:
        f_data = f.variables['lat'].data.copy()
        f.close()

        return f_data


def get_cru_lon(my_dir, my_date, my_voi):
    """
    Name:     get_cru_lon
    Inputs:   - str, directory to CRU TS netcdf file (my_dir)
              - datetime.date, current month datetime object (my_date)
              - str, variable of interest, for CRU file only (my_voi)
    Outputs:  numpy.ndarray, longitudes (degrees)
    Features: Returns an array of longitudes from a CRU TS netCDF file
    """
    try:
        my_file = get_cru_file(my_dir, my_date, my_voi)
    except:
        raise
    else:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")

        # Save data for variables of interest:
        # NOTE: for CRU TS 4.00:
        #       variables: 'lat', 'lon', 'time', voi, 'stn'
        #       where
        #           lat (latitude): -89.75 -- 89.75
        #           > units: degrees_north
        #           lon (longitude): -179.75 -- 179.75
        #           > units: degrees_east

        # Create a copy of the latitude data:
        f_data = f.variables['lon'].data.copy()
        f.close()

        return f_data


def get_monthly_cru(my_dir, my_date, my_voi):
    """
    Name:     get_monthly_cru
    Input:    - string, directory to CRU netcdf file (my_dir)
              - datetime.date, current month datetime object (my_date)
              - string, variable of interest (my_voi)
    Output:   numpy nd.array
    Depends:  - find_files
              - get_cru_file
              - get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and
              variable of interest (e.g., cld, pre, tmp)
    """
    try:
        my_file = get_cru_file(my_dir, my_date, my_voi)
    except:
        raise
    else:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")

        # Save data for variables of interest:
        # NOTE: for CRU TS 4.00:
        #       variables: 'lat', 'lon', 'time', voi, 'stn'
        #       where
        #           lat (latitude): -89.75 -- 89.75
        #           lon (longitude): -179.75 -- 179.75
        #           time:
        #           > units: days since 1900-1-1
        #           > shape: (1344,)
        #           > values: mid-day of each month (e.g., 15th or 16th)
        #           voi: 'tmp', 'pre', 'cld', etc.
        #           > 'cld' units = %
        #           > 'pre' units = mm
        #           > 'tmp' units = deg. C
        #           > Missing value = 9.96e+36

        # Save the base time stamp:
        bt = datetime.date(1900, 1, 1)

        # Read the time data as array:
        f_time = f.variables['time'].data

        # Find the time index for the current date:
        ti = get_time_index(bt, my_date, f_time)

        # Get the spatial data for current time:
        f_data = f.variables[my_voi].data[ti].copy()

        f_time = None
        f.close()

        return f_data


def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime.date, base timestamp
              - datetime.date, current timestamp
              - numpy nd.array, days since base timestamp
    Output:   int
    Features: Finds the index in an array of CRU TS days for a given timestamp
    """
    # For CRU TS 3.21, the aot is indexed for mid-month days, e.g. 15--16th
    # therefore, to make certain that ct index preceeds the index for the
    # correct month in aot, make the day of the current month less than
    # the 15th or 16th (i.e., replace day with '1'):
    ct = ct.replace(day=1)

    # Calculate the time difference between ct and bt:
    dt = (ct - bt).days

    # Append dt to the aot array:
    aot = numpy.append(aot, [dt, ])

    # Find the first index of dt in the sorted array:
    idx = numpy.where(numpy.sort(aot) == dt)[0][0]

    return idx


def process_cru(my_voi, my_dir):
    """
    Name:     process_cru
    Input:    - str, variable name, e.g., tmp, pre, cld (my_voi)
              - str, directory name for CRU TS data file (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF file by month into variable list and data
              set table output files
    Depends:  - add_one_month
              - get_station_latlon
              - writeout
    """
    # Define the start and end dates you want to process:
    # NOTE: this is a hard-coded date range
    start_date = datetime.date(1991, 1, 1)
    end_date = datetime.date(2015, 1, 1)

    # Set flag for varlist:
    varlist_flag = True

    # Get list of flux station 0.5-degree grid points:
    station_list = get_station_latlon()

    # Prepare var output file and output file header lines:
    var_file = "CRU_Var-List_" + my_voi + ".csv"
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Define and/or create the output directory (subdir of my_dir)
    out_dir = os.path.join(my_dir, "out")
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except:
            logging.warning(
                "failed to create output directory; using input directory")
            out_dir = my_dir
        else:
            logging.info("created output directory %s", out_dir)

    var_path = os.path.join(out_dir, var_file)
    writeout(var_path, var_headerline)

    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # Open and read netcdf files in the file directory:
        try:
            # Open netCDF file for reading:
            cru_file = get_cru_file(my_dir, cur_date, my_voi)
            f = netcdf.NetCDFFile(cru_file, "r")
        except IndexError as e:
            logging.error("Failed to find CRU TS file. %s", str(e))
        except:
            logging.error("Failed to read CRU TS file.")
        else:
            # Save the data shape:
            sh_time, sh_lat, sh_lon = f.variables[my_voi].shape

            # Find the time index for the current date:
            b_time = datetime.date(1900, 1, 1)
            f_time = f.variables['time'].data
            ti = get_time_index(b_time, cur_date, f_time)

            # Prepare data set output file:
            dat_file = "CRU_Data-Set_%s_%s.csv" % (my_voi, cur_date)
            dat_path = os.path.join(out_dir, dat_file)
            writeout(dat_path, dat_headerline)

            # Iterate through each lat:lon pair
            # * row-major ordering from bottom left
            #  x (longitude): 0...719
            #  y (latitude): 0...359
            for y in range(sh_lat):
                pxl_lat = f.variables['lat'].data[y]
                for x in range(sh_lon):
                    pxl_lon = f.variables['lon'].data[x]

                    # Filter grids based on flux stations:
                    if (pxl_lat, pxl_lon) in station_list:
                        # Calc station ID:
                        st_id = 720*y + x
                        station_parts = ('HDG', st_id)

                        if my_voi == 'tmp':
                            my_var = "Tc"
                        elif my_voi == 'pre':
                            my_var = "Pre"
                        elif my_voi == 'cld':
                            my_var = "Cld"
                        my_line = VAR(station_parts, my_var, 'grid')

                        if varlist_flag:
                            # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                            my_line.printLine(var_path)
                            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                        # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                        # Read each each pixel
                        # NOTE: missing values are ~ 1e7
                        pxl_voi = f.variables[my_voi].data[ti, y, x]
                        if pxl_voi < 1.e6:
                            # Append to existing file:
                            OUT = open(dat_path, 'a')
                            outline = "%s,%s,%s,%0.3f\n" % (
                                my_line.msvIDX,
                                my_line.stationID,
                                cur_date,
                                pxl_voi)
                            OUT.write(outline)
                            OUT.close()
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
            # Turn off varlist flag after processing first month:
            varlist_flag = False

            # Increment cur_date:
            cur_date = add_one_month(cur_date)

            f_time = None
            pxl_lat = None
            pxl_lon = None
            pxl_voi = None
            f.close()


def process_cru_elv(my_dir):
    """
    Name:     process_cru_elv
    Input:    string, input/output file directory (my_dir)
    Output:   None.
    Features: Processes CRU TS 3.22 data file associated with elevation into
              variable list and data set table output files
    Depends:  - VAR class
              - find_files
              - get_station_latlon
              - writeout
    """
    var_file = "CRU_Var-List_elv.csv"
    dat_file = "CRU_Data-Set_elv.csv"
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Define and/or create the output directory (subdir of my_dir)
    out_dir = os.path.join(my_dir, "out")
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except:
            logging.warning(
                "failed to create output directory; using input directory")
            out_dir = my_dir
        else:
            logging.info("created output directory %s", out_dir)

    var_path = os.path.join(out_dir, var_file)
    dat_path = os.path.join(out_dir, dat_file)

    # Get list of flux station 0.5-degree grid points:
    station_list = get_station_latlon()

    try:
        # Search directory for CRU dat file:
        my_files = find_files(my_dir, "*.dat")
        my_file = my_files[0]
    except IndexError:
        logging.error("No CRU TS elevation data file found!")
    else:
        # Open and read dat file:
        # NOTE: data is read into an array with shape (360, 720)
        #      'lat' goes from -89.75 -- 89.75 (south to north)
        #      'lon' goes from -179.75 -- 179.75 (east to west)
        f = numpy.loadtxt(my_file)
        (sh_lat, sh_lon) = f.shape
        latitude = [-89.75 + i*0.5 for i in range(sh_lat)]
        longitude = [-179.75 + i*0.5 for i in range(sh_lon)]

        # Assign time stamp as CRU TS 3.00 date:
        time_stamp = datetime.date(2006, 6, 1)

        # Prepare var and data output files:
        writeout(var_path, var_headerline)
        writeout(dat_path, dat_headerline)

        # Iterate through data:
        for y in range(sh_lat):
            pxl_lat = latitude[y]
            for x in range(sh_lon):
                pxl_lon = longitude[x]

                # Filter grids based on flux stations:
                if (pxl_lat, pxl_lon) in station_list:
                    # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                    st_id = 720*y + x
                    station_parts = ('HDG', st_id)
                    my_line = VAR(station_parts, 'Elv', 'grid')
                    my_line.printLine(var_path)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                    elv = f[y, x]
                    # Set no-data value:
                    if elv > -500:
                        # Append to existing file:
                        OUT = open(dat_path, 'a')
                        # Create/write output line:
                        outline = "%s,%s,%s,%0.5f\n" % (
                            my_line.msvIDX,
                            my_line.stationID,
                            time_stamp,
                            elv)
                        OUT.write(outline)
                        OUT.close()


def process_cru_vpd(my_dir):
    """
    Name:     process_cru_vpd
    Input:    str, directory for CRU input files (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF files (tmp and vap) to calculate
              monthly VPD and save variable list and data set table output
              files
    Depends:  - add_one_month
              - calculate_vpd
              - writeout
    """
    # Define the start and end dates you want to process:
    start_date = datetime.date(1991, 1, 1)
    end_date = datetime.date(2015, 1, 1)

    # Set flag for varlist:
    varlist_flag = True

    # Get lon/lat arrays (throws exception if file not found):
    latitude = get_cru_lat(my_dir, start_date, 'tmp')
    longitude = get_cru_lon(my_dir, start_date, 'tmp')

    # Get list of flux station 0.5-degree grid points:
    station_list = get_station_latlon()

    # Define var output file and headerlines:
    var_file = "CRU_Var-List_vpd.csv"
    var_headerline = "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
    dat_headerline = "msvidx,stationid,datetime,data\n"

    # Define and/or create the output directory (subdir of my_dir)
    out_dir = os.path.join(my_dir, "out")
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except:
            logging.warning(
                "failed to create output directory; using input directory")
            out_dir = my_dir
        else:
            logging.info("created output directory %s", out_dir)

    # Prepare var output file:
    var_path = os.path.join(out_dir, var_file)
    writeout(var_path, var_headerline)

    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # Open and read netcdf files in the file directory:
        tmp = get_monthly_cru(my_dir, cur_date, "tmp")
        vap = get_monthly_cru(my_dir, cur_date, "vap")

        # Calculate VPD & save shape:
        vpd = calculate_vpd(vap, tmp)
        (sh_lat, sh_lon) = vpd.shape

        # Prepare data set output file:
        dat_file = "CRU_Data-Set_vpd_%s.csv" % cur_date
        dat_path = os.path.join(out_dir, dat_file)
        writeout(dat_path, dat_headerline)

        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in range(sh_lat):
            pxl_lat = latitude[y]
            for x in range(sh_lon):
                pxl_lon = longitude[x]

                # Filter grids based on flux stations:
                if (pxl_lat, pxl_lon) in station_list:
                    # Calc station ID:
                    st_id = 720*y + x
                    station_parts = ('HDG', st_id)

                    my_line = VAR(station_parts, 'VPD', 'grid')

                    if varlist_flag:
                        # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                        my_line.printLine(var_path)
                        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                    # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                    # Read VPD for each pixel
                    # NOTE: missing VPD are equal to -9999
                    pxl_vpd = vpd[y, x]
                    if pxl_vpd > -9999.0:
                        # Append to existing file:
                        OUT = open(dat_path, 'a')

                        # Create/write output line:
                        outline = "%s,%s,%s,%0.5f\n" % (
                            my_line.msvIDX,
                            my_line.stationID,
                            cur_date,
                            pxl_vpd
                            )
                        OUT.write(outline)
                        OUT.close()
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Turn off varlist flag after processing first month:
        varlist_flag = False

        # Increment cur_date:
        cur_date = add_one_month(cur_date)
