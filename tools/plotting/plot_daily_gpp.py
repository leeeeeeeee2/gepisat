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
def get_meta(my_path, st_name):
    """
    Name:     get_meta
    Input:    - str, path metadata file (my_path)
              - str, station name (st_name)
    Output:   float, elevation, m (my_elv)
    Features: Returns the station meta data acquired from the GePiSaT
              database metadata file
    """
    try:
        my_data = numpy.loadtxt(fname=my_path,
                                dtype={'names': ('stationid', 'ele', 'classid',
                                                 'climateid'),
                                       'formats': ('S6', 'f4', 'S3', 'S3')},
                                delimiter=',',
                                skiprows=1,
                                usecols=(4, 8, 9, 11),
                                converters={4: numpy.str,
                                            8: numpy.float,
                                            9: numpy.str,
                                            11: numpy.str})
    except:
        return numpy.array([])
    else:
        my_idx = numpy.where(my_data['stationid'] == st_name)[0]
        if my_idx:
            return my_data[my_idx]
        else:
            return numpy.array([])


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


def plot_daily_gpp(daily_data, my_text):
    """
    """
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
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
    ax1.text(0.05, 0.95, my_text, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', bbox=props)
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
    my_file = "FI-Sod_daily_GPP.txt"
    flx_met_path = "Fluxdata_Met-Data.csv"
    my_data = load_daily_gpp(my_file)
    my_station = get_station_name(my_file)
    my_meta = get_meta(flx_met_path, my_station)
    if my_meta.size == 1:
        my_elv = my_meta['ele']
        my_veg = my_meta['classid']
        my_clim = my_meta['climateid']
        if my_elv == -9999:
            my_txt = ("$\\mathrm{%s}$ ($\\mathrm{%s}/\\mathrm{%s}$)") % (
                my_station, my_veg, my_clim)
        else:
            my_txt = (
                "$\\mathrm{%s}$ ($\\mathrm{%s}/\\mathrm{%s}$)\n"
                "$z=%0.1f$ m") % (my_station, my_veg, my_clim, my_elv)
        print("%s" % (my_txt))
    else:
        my_txt = ""
    plot_daily_gpp(my_data, my_txt)
