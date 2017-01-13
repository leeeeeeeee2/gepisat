#!/usr/bin/python
#
# crudata.py
#
# VERSION 3.0.0-dev
# LAST UPDATED: 2017-01-13
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
import glob
import logging
import os

import numpy
from scipy.io import netcdf

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
    Features: Returns an array of mean monthly vapor pressure deficit
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
    vpd = -9999.0*(numpy.zeros(shape=(360, 720)) + 1)

    # Iterate through each data point:
    lat, lon = tm.shape

    for y in range(lat):
        for x in range(lon):
            ea = vap[y, x]
            if tm is not None:
                tmin = tm[y, x]
                tmax = tx[y, x]
                tc = 0.5*(tmin + tmax)
            else:
                tc = tm

            if tc < 1.e6 and ea < 1.e6:
                vpd[y, x] = (
                    0.611*numpy.exp((17.27*tc)/(tc + 237.3)) - 0.10*ea)
    return vpd


def get_monthly_cru(d, ct, v):
    """
    Name:     get_monthly_cru
    Input:    - string, directory to CRU netcdf file (d)
              - datetime.date, current month datetime object (ct)
              - string, variable of interest (v)
    Output:   numpy nd.array
    Depends:  get_time_index
    Features: Returns 360x720 monthly CRU TS dataset for a given month and
              variable of interest (e.g., cld, pre, tmp)
    """
    # Search directory for netCDF file:
    s_val = "*%s.dat.nc" % (v)
    s_str = os.path.join(d, s_val)
    my_files = glob.glob(s_str)

    if len(my_files) == 1:
        my_file = my_files[0]
    else:
        my_file = None

    if my_file:
        # Open netCDF file for reading:
        f = netcdf.NetCDFFile(my_file, "r")

        # Save data for variables of interest:
        # NOTE: for CRU TS 3.21:
        #       variables: 'lat', 'lon', 'time', v
        #       where v is 'tmp', 'pre', 'cld'
        # LAT:  -89.75 -- 89.75
        # LON:  -179.75 -- 179.75
        # TIME:
        #       units: days since 1900-1-1
        #       shape: (1344,)
        #       values: mid-day of each month (e.g., 15th or 16th)
        # DATA:
        #       'cld' units = %
        #       'pre' units = mm
        #       'tmp' units = deg. C
        #       Missing value = 9.96e+36
        # Save the base time stamp:
        bt = datetime.date(1900, 1, 1)

        # Read the time data as array:
        f_time = f.variables['time'].data

        # Find the time index for the current date:
        ti = get_time_index(bt, ct, f_time)

        # Get the spatial data for current time:
        f_data = f.variables[v].data[ti].copy()

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


def process_cru(v, cru_dir, my_dir):
    """
    Name:     process_cru
    Input:    - string, variable name, e.g., tmp, pre, cld (v)
              - string, directory name for CRU TS data file (cru_dir)
              - string, directory for output files (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF file by month into variable list and data
              set table output files
    Depends:  - add_one_month
              - get_monthly_cru
              - writeout
    """
    # Define the start and end dates you want to process (2002-2006):
    start_date = datetime.date(2002, 1, 1)
    end_date = datetime.date(2007, 1, 1)

    # Set flag for varlist:
    varlist_flag = 1

    # Prepare var output file:
    my_var_out = "CRU_Var-List_" + v + ".csv"
    var_outfile = my_dir + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)

    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # # Open and read netcdf files in the file directory:
        my_data = get_monthly_cru(cru_dir, cur_date, v)
        (sh_lat, sh_lon) = my_data.shape

        # Prepare data set output file:
        my_dat_out = "CRU_Data-Set_%s_%s.csv" % (v, cur_date)
        dat_outfile = my_dir + my_dat_out
        dat_headerline = "msvidx,stationid,datetime,data\n"
        writeout(dat_outfile, dat_headerline)

        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in range(sh_lat):
            for x in range(sh_lon):
                # Calc station ID:
                st_id = 720*y + x
                station_parts = ('HDG', st_id)

                if v == 'tmp':
                    my_var = "Tc"
                elif v == 'pre':
                    my_var = "Pre"
                elif v == 'cld':
                    my_var = "Cld"
                my_line = VAR(station_parts, my_var, 'grid')

                if varlist_flag:
                    # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                    my_line.printLine(var_outfile)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Read each each pixel
                # * NOTE: missing values are ~ 1e7
                pxl = my_data[y, x]
                if pxl < 1.e6:
                    # Append to existing file:
                    OUT = open(dat_outfile, 'a')

                    # Create/write output line:
                    outline = "%s,%s,%s,%0.3f\n" % (
                        my_line.msvIDX,
                        my_line.stationID,
                        cur_date,
                        pxl
                        )
                    OUT.write(outline)
                    OUT.close()
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Turn off varlist flag after processing first month:
        varlist_flag = 0

        # Increment cur_date:
        cur_date = add_one_month(cur_date)


def process_cru_elv(d):
    """
    Name:     process_cru_elv
    Input:    string, input/output file directory (d)
    Output:   None.
    Features: Processes CRU TS 3.22 data file associated with elevation into
              variable list and data set table output files
    Depends:  - VAR class
              - writeout
    """
    my_dir = d

    try:
        # Search directory for CRU dat file:
        s_str = os.path.join(my_dir, "*.dat")
        my_file = glob.glob(s_str)[0]
    except IndexError:
        logging.error("No CRU TS elevation data file found!")
    else:
        # Open and read dat file:
        # NOTE: data is read into an array with shape (360, 720)
        #      'lat' goes from -89.75 -- 89.75 (south to north)
        #      'lon' goes from -179.75 -- 179.75 (east to west)
        f = numpy.loadtxt(my_file)
        (sh_lat, sh_lon) = f.shape

        # Assign time stamp as CRU TS 3.00 date:
        time_stamp = datetime.date(2006, 6, 1)

        # Prepare var output file:
        my_var_out = "CRU_Var-List_elv.csv"
        var_outfile = os.path.join(my_dir, my_var_out)
        var_headerline = (
            "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
            )
        writeout(var_outfile, var_headerline)

        # Prepare data set output file:
        my_dat_out = "CRU_Data-Set_elv.csv"
        dat_outfile = os.path.join(my_dir, my_dat_out)
        dat_headerline = "msvidx,stationid,datetime,data\n"
        writeout(dat_outfile, dat_headerline)

        # Iterate through data:
        for y in range(sh_lat):
            for x in range(sh_lon):
                # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                st_id = 720*y + x
                station_parts = ('HDG', st_id)
                my_line = VAR(station_parts, 'Elv', 'grid')
                my_line.printLine(var_outfile)
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                elv = f[y, x]
                # Set no-data value:
                if elv > -500:
                    # Append to existing file:
                    OUT = open(dat_outfile, 'a')
                    # Create/write output line:
                    outline = "%s,%s,%s,%0.5f\n" % (
                        my_line.msvIDX,
                        my_line.stationID,
                        time_stamp,
                        elv
                        )
                    OUT.write(outline)
                    OUT.close()


def process_cru_vpd(cru_dir, my_dir):
    """
    Name:     process_cru_vpd
    Input:    - string, directory for CRU input files (cru_dir)
              - string, output file directory (my_dir)
    Output:   None.
    Features: Processes CRU TS netCDF files (tmn, tmx, vap) to calculate
              monthly VPD and save variable list and data set table output
              files
    Depends:  - add_one_month
              - calculate_vpd
              - writeout
    """
    # Define the start and end dates you want to process (2002-2006):
    start_date = datetime.date(2002, 1, 1)
    end_date = datetime.date(2007, 1, 1)

    # Set flag for varlist:
    varlist_flag = 1

    # Prepare var output file:
    my_var_out = "CRU_Var-List_vpd.csv"
    var_outfile = my_dir + my_var_out
    var_headerline = (
        "msvidx,stationid,varid,varname,varunit,vartype,varcore\n"
        )
    writeout(var_outfile, var_headerline)

    # Process each month between start and end dates:
    cur_date = start_date
    while cur_date < end_date:
        # Open and read netcdf files in the file directory:
        tmn = get_monthly_cru(cru_dir, cur_date, "tmn")
        tmx = get_monthly_cru(cru_dir, cur_date, "tmx")
        vap = get_monthly_cru(cru_dir, cur_date, "vap")

        # Calculate VPD & save shape:
        vpd = calculate_vpd(tmn, tmx, vap)
        (sh_lat, sh_lon) = vpd.shape

        # Prepare data set output file:
        my_dat_out = "CRU_Data-Set_vpd_%s.csv" % cur_date
        dat_outfile = my_dir + my_dat_out
        dat_headerline = "msvidx,stationid,datetime,data\n"
        writeout(dat_outfile, dat_headerline)

        # Iterate through each lat:lon pair
        # * row-major ordering from bottom left
        #  x (longitude): 0...719
        #  y (latitude): 0...359
        for y in range(sh_lat):
            for x in range(sh_lon):
                # Calc station ID:
                st_id = 720*y + x
                station_parts = ('HDG', st_id)

                my_line = VAR(station_parts, 'VPD', 'grid')

                if varlist_flag:
                    # ~~~~~~~~~~~~~~ VAR-LIST ~~~~~~~~~~~~~~ #
                    my_line.printLine(var_outfile)
                    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

                # ~~~~~~~~~~~~~~~~ DATA-SET ~~~~~~~~~~~~~~~~ #
                # Read VPD for each pixel
                # * NOTE: missing VPD are equal to -9999
                pxl_vpd = vpd[y, x]
                if pxl_vpd > -9999.0:
                    # Append to existing file:
                    OUT = open(dat_outfile, 'a')

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
        varlist_flag = 0

        # Increment cur_date:
        cur_date = add_one_month(cur_date)
