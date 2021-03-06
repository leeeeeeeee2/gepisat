#!/usr/bin/python
#
# plot_daily_gpp.py
#
# LAST UPDATED: 2017-02-05
#
# ------------
# description:
# ------------
# To combine PDF files, try: `pdfjoin --rotateoversize false *.pdf`
#
# @TODO plot monthly
#
###############################################################################
# IMPORT MODULES
###############################################################################
import datetime
import glob
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
        if len(my_idx) == 1:
            return my_data[my_idx]
        else:
            return numpy.array([])


def get_station_name(my_path):
    """
    Name:     get_station_name
    Inputs:   str, path and file name (my_path)
    Outputs:  str, station name
    Features: Return the station name extracted from the output file name
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
            converters={
                0: lambda x: datetime.datetime.strptime(x, '%Y-%m-%d').date(),
                1: numpy.float,
                2: numpy.float}
        )
    return my_data


def plot_daily_gpp(daily_data, left_text, right_text, to_save, out_file):
    """
    Name:     plot_daily_gpp
    Inputs:   - numpy.ndarray, structured array with daily GPP (daily_data)
              - str, text to show in the top-left margin (left_text)
              - str, text to show in the top-right margin (right_text)
              - bool, whether to save the figure to file (to_save)
              - str, path to output figure (out_file)
    Outputs:  None.
    Features: Produces a plot of daily GPP with text given in the top-left and
              top-right margin space
    @TODO:    _ set the x-axis range from 2002-01-01 to 2007-01-01
    """
    fig = plt.figure(figsize=(8, 6), dpi=180)
    ax1 = fig.add_subplot(111)
    # plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    # plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.plot(daily_data['date'], daily_data['gpp'], '-')
    ax1.set_ylabel('GPP, mol CO$_2$ m$^{-2}$ day$^{-1}$', fontsize=12)
    ax1.set_xlabel('Date', fontsize=12)
    ax1.set_xlim([datetime.date(2001, 10, 1), datetime.date(2007, 4, 1)])
    ax1.text(0, 1.1, left_text, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', horizontalalignment='left')
    ax1.text(1, 1.1, right_text, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', horizontalalignment='right')

    if to_save:
        fig.savefig(out_file, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == "__main__":
    # Define input files and directories:
    data_dir = os.path.join(os.path.expanduser("~"), "Data", "psql", "flux")
    flx_met_path = os.path.join(data_dir, "Fluxdata_Met-Data.csv")

    #
    # PLOT DAILY GPP
    #
    out_dir = os.path.join(
        os.path.expanduser("~"), "Desktop", "temp", "out", "daily_gpp")
    s_pattern = os.path.join(out_dir, "*daily_GPP.txt")
    my_files = glob.glob(s_pattern)
    for my_file in sorted(my_files):
        my_data = load_daily_gpp(my_file)
        my_station = get_station_name(my_file)
        my_meta = get_meta(flx_met_path, my_station)

        if my_meta.size == 1:
            my_elv = my_meta['ele'][0]
            my_veg = my_meta['classid'][0]
            my_clim = my_meta['climateid'][0]
            if my_elv == -9999:
                my_txt1 = ("$\\mathrm{%s}$\n") % (my_station)
                my_txt2 = (
                    "$\\mathrm{%s:}$ $\\mathrm{%s}$\n"
                    "$\\mathrm{%s:}$ $\\mathrm{%s}$") % (
                        "Veg", my_veg, "Clim", my_clim)
            else:
                my_txt1 = (
                    "$\\mathrm{%s}$\n$z=%0.1f$ m") % (my_station, my_elv)
                my_txt2 = (
                    "$\\mathrm{%s:}$ $\\mathrm{%s}$\n"
                    "$\\mathrm{%s:}$ $\\mathrm{%s}$") % (
                        "Veg", my_veg, "Clim", my_clim)
        else:
            my_txt1 = ("$\\mathrm{%s}$\n") % (my_station)
            my_txt2 = ""

        fig_file = "%s_daily_gpp.pdf" % (my_station)
        fig_path = os.path.join(out_dir, "figs", fig_file)
        plot_daily_gpp(my_data, my_txt1, my_txt2, True, fig_path)
