#!/usr/bin/python
#
# plot_daily_gpp.py
#
# LAST UPDATED: 2017-01-29
#
# ------------
# description:
# ------------
#
###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import os
import re

import numpy
import matplotlib.pyplot as plt


###############################################################################
# FUNCTIONS
###############################################################################
def get_elv(my_path, st_name):
    """
    Name:     get_elv
    Input:    - str, path metadata file (my_path)
              - str, station name (st_name)
    Output:   float, elevation, m (my_elv)
    Features: Returns the station elevation from GePiSaT database metadata file
    """
    try:
        my_data = numpy.loadtxt(fname=my_path,
                                dtype={'names': ('stationid', 'ele'),
                                       'formats': ('S6', 'f4')},
                                delimiter=',',
                                skiprows=1,
                                usecols=(4, 8),
                                converters={4: numpy.str,
                                            8: numpy.float})
    except:
        return numpy.nan
    else:
        my_idx = numpy.where(my_data['stationid'] == st_name)[0]
        if my_idx:
            my_elv = my_data['ele'][my_idx[0]]
            if my_elv == -9999:
                my_elv = numpy.nan
        else:
            my_elv = numpy.nan

        return my_elv


def get_veg_type(my_path, st_name):
    """
    Name:     get_veg_type
    Inputs:   - str, path meta data file (my_path)
              - str, station name (st_name)
    Output:   str, vegation type, short name
    Features: Returns the vegation type for a given station
    """
    try:
        my_data = numpy.loadtxt(fname=my_path,
                                dtype={'names': ('stationid', 'classid'),
                                       'formats': ('S6', 'S3')},
                                delimiter=',',
                                skiprows=1,
                                usecols=(4, 9),
                                converters={4: numpy.str,
                                            9: numpy.str})
    except:
        return ''
    else:
        my_idx = numpy.where(my_data['stationid'] == st_name)
        try:
            my_veg = my_data['classid'][my_idx][0]
        except:
            my_veg = ''

        return my_veg


def get_station_name(my_path):
    """
    """
    my_station = ""
    my_file = os.path.basename(my_path)
    try:
        my_station = re.search('^\S{6}', my_file).group(0)
    except:
        return ""
    else:
        return my_station


def load_daily_gpp(my_path):
    """
    Name:     load_daily_gpp
    Inputs:   str, file path (my_path)
    Outputs:  numpy.ndarray
              > datetime.date, 'date'
              > float, 'gpp'
              > float, 'gpp_err'
    Features: Returns daily GPP file (output from GePiSaT model) as a
              structured array
    """
    my_data = numpy.array([])
    if os.path.isfile(my_path):
        my_data = numpy.loadtxt(
            fname=my_path,
            dtype={'names': ('date', 'gpp', 'gpp_err'),
                   'formats': ('O', 'f4', 'f4')},
            delimiter=',',
            skiprows=1,
            converters={0: lambda x: datetime.datetime.strptime(
                x, '%Y-%m-%d').date(),
                        1: numpy.float,
                        2: numpy.float}
        )
    return my_data


def plot_daily_gpp(daily_data):
    """
    """
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # plot_max = numpy.concatenate((my_fit, my_obs)).max()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    # plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(
        daily_data['date'],
        daily_data['gpp'],
        '-k'
        )
    # ax1.fill_between(
    #    daily_data['date'],
    #    daily_data['gpp'] + daily_data['gpp_err'],
    #    daily_data['gpp'] - daily_data['gpp_err'],
    #    facecolor='r'
    #    )
    ax1.set_ylabel('Daily GPP, mol CO$_2$ m$^{-2}$', fontsize=12)
    ax1.set_xlabel('Date', fontsize=12)
    # ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=3, mode="expand", borderaxespad=0., fontsize=14)
    # ax1.text(0.05, 0.95, my_text, transform=ax1.transAxes, fontsize=14,
    #        verticalalignment='top', bbox=props)
    plt.show()


def plot_monthly_gpp(daily_data):
    """
    """
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # plot_max = numpy.concatenate((my_fit, my_obs)).max()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    # plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.errorbar(daily_data['date'],
                 daily_data['gpp'],
                 yerr=daily_data['gpp_err'],
                 ecolor='k',
                 fmt='o'
                 )
    ax1.set_ylabel('GPP, mol CO$_2$ m$^{-2}$', fontsize=14)
    ax1.set_xlabel('Date', fontsize=14)
    # ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=3, mode="expand", borderaxespad=0., fontsize=14)
    # ax1.text(0.05, 0.95, my_text, transform=ax1.transAxes, fontsize=14,
    #        verticalalignment='top', bbox=props)
    plt.show()

###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == "__main__":
    my_file = "AT-Neu_daily_GPP.txt"
    flx_met_path = "Fluxdata_Met-Data.csv"
    my_data = load_daily_gpp(my_file)
    my_station = get_station_name(my_file)
    my_veg_type = get_veg_type(flx_met_path, my_station)
    my_elv = get_elv(flx_met_path, my_station)
    my_txt1 = ("$\\mathrm{%s}$ ($\\mathrm{%s}$)\n$z=%0.1f$ m") % (
        my_station, my_veg_type, my_elv)
    print "%s" % (my_txt1)
    # plot_daily_gpp(my_data)
