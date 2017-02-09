#!/usr/bin/python
#
# plot_monthly_gpp.py
#
# LAST UPDATED: 2017-02-08
#
# ------------
# description:
# ------------
# To combine PDF files, try: `pdfjoin --rotateoversize false *.pdf`
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


def gapfill_months(data):
    """
    Name:     gapfill_months
    Inputs:   numpy.ndarray, structured array (data)
    Outputs:  numpy.ndarray, gapfilled structured array
    Features: Returns structured array with months missing GPP gap-filled with
              NaNs
    """
    if data.size > 1:
        # Initialize return data structured array:
        rdata = numpy.array(
            (data['date'][0], data['gpp'][0], data['gpp_err'][0]),
            dtype={'names': ('date', 'gpp', 'gpp_err'),
                   'formats': ('O', 'f4', 'f4')},
            ndmin=1)

        first_mo = data['date'][0]
        last_mo = data['date'][-1]
        cur_date = add_one_month(first_mo)
        while cur_date <= last_mo:
            idx = numpy.where(data['date'] == cur_date)[0]
            if len(idx) == 1:
                temp_array = numpy.array(
                    (data['date'][idx][0],
                     data['gpp'][idx][0],
                     data['gpp_err'][idx][0]),
                    dtype={'names': ('date', 'gpp', 'gpp_err'),
                           'formats': ('O', 'f4', 'f4')},
                    ndmin=1)
            else:
                temp_array = numpy.array(
                    (cur_date, numpy.nan, 0.0),
                    dtype={'names': ('date', 'gpp', 'gpp_err'),
                           'formats': ('O', 'f4', 'f4')},
                    ndmin=1)
            rdata = numpy.append(rdata, temp_array, axis=0)
            cur_date = add_one_month(cur_date)
    else:
        rdata = data
    return rdata


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


def load_gepisat_data(my_path):
    """
    Name:     load_gepisat_data
    Inputs:   str, file path to concatenated station data (my_path)
    Outputs:  numpy.ndarray, structured array
    Features: Returns a structured array of time-stamped station GPP and
              associated errors from concatenated station LUE files
    """
    my_data = numpy.array([])
    if os.path.isfile(my_path):
        my_data = numpy.loadtxt(
            fname=my_path,
            dtype={'names': ('station', 'date', 'gpp', 'gpp_err'),
                   'formats': ('S6', 'O', 'f4', 'f4')},
            delimiter=',',
            skiprows=1,
            usecols=(0, 1, 2, 3),
            converters={
                0: numpy.str,
                1: lambda x: datetime.datetime.strptime(x, '%Y-%m-%d').date(),
                2: numpy.float,
                3: numpy.float}
        )
    return my_data


def load_monthly_gpp(my_path):
    """
    Name:     load_monthly_gpp
    Inputs:   str, file path (my_path)
    Outputs:  numpy.ndarray
              > datetime.date, 'date'
              > float, 'gpp'
              > float, 'gpp_err'
    Features: Returns monthly GPP file (output from GePiSaT model) as a
              structured array; gapfills missing months with NaNs
    Depends:  - add_one_month
              - gapfill_months
    @TODO: make gapfilling an option passed as function argument
    """
    my_data = numpy.array([])
    if os.path.isfile(my_path):
        data = numpy.loadtxt(
            fname=my_path,
            dtype={'names': ('date', 'gpp', 'gpp_err'),
                   'formats': ('O', 'f4', 'f4')},
            delimiter=',',
            skiprows=1,
            usecols=(0, 1, 2),
            converters={
                0: lambda x: datetime.datetime.strptime(x, '%Y-%m-%d').date(),
                1: numpy.float,
                2: numpy.float}
        )
        my_data = gapfill_months(data)

    return my_data


def plot_monthly_gpp(mo_data, left_text, right_text, to_save, out_file):
    """
    Name:     plot_monthly_gpp
    Inputs:   - numpy.ndarray, structured array with monthly GPP (mo_data)
              - str, text to show in the top-left margin (left_text)
              - str, text to show in the top-right margin (right_text)
              - bool, whether to save the figure to file (to_save)
              - str, path to output figure (out_file)
    Outputs:  None.
    Features: Produces a plot of monthly GPP with text given in the top-left
              and top-right margin space
    Note:     use scatter to make certain there are dots for each errorbar!
    """
    fig = plt.figure(figsize=(8, 6), dpi=180)
    ax1 = fig.add_subplot(111)
    # plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
    # plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
    ax1.scatter(mo_data['date'], mo_data['gpp'], edgecolors='none')
    ax1.errorbar(mo_data['date'], mo_data['gpp'], yerr=mo_data['gpp_err'],
                 fmt='--o', ecolor='r')
    ax1.set_ylabel('GPP, mol CO$_2$ m$^{-2}$ month$^{-1}$', fontsize=12)
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


def plot_gepisat_gpp(xdata, ydata, left_text, right_text, to_save, out_file):
    """
    Name:     plot_gepisat_gpp
    Inputs:   - numpy.ndarray, structured array with 2015 GPP (xdata)
              - numpy.ndarray, structured array with 2012 GPP (ydata)
              - str, text to show in the top-left margin (left_text)
              - str, text to show in the top-right margin (right_text)
              - bool, whether to save the figure to file (to_save)
              - str, path to output figure (out_file)
    Outputs:  None.
    Features: Produces a plot of monthly GPP comparing 2012 with 2015 data
              with text given in the top-left and top-right margin space
    """
    fig = plt.figure(figsize=(8, 6), dpi=180)
    ax1 = fig.add_subplot(111)

    ax1.plot(xdata['date'], xdata['gpp'], ':^', lw=1, label="2015")
    ax1.plot(ydata['date'], ydata['gpp'], 'r:v', lw=1, label="2012")

    ax1.set_ylabel('GPP, mol CO$_2$ m$^{-2}$ mo$^{-1}$', fontsize=12)
    ax1.set_xlabel('Date', fontsize=12)
    ax1.set_xlim([datetime.date(2001, 10, 1), datetime.date(2007, 4, 1)])
    ax1.text(0, 1.1, left_text, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', horizontalalignment='left')
    ax1.text(1, 1.1, right_text, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', horizontalalignment='right')
    ax1.legend(bbox_to_anchor=(0.3, 1.02, 1., .102), loc=3,
               ncol=2, mode=None, borderaxespad=0., fontsize=12)
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
    # PLOT MONTHLY GPP
    #
    to_plot = False
    if to_plot:
        out_dir = os.path.join(
            os.path.expanduser("~"), "Desktop", "temp", "out", "lue")
        s_pattern = os.path.join(out_dir, "*_LUE.txt")
        my_files = glob.glob(s_pattern)
        for my_file in sorted(my_files):
            my_data = load_monthly_gpp(my_file)

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

            fig_file = "%s_monthly_gpp.pdf" % (my_station)
            fig_path = os.path.join(out_dir, "figs", fig_file)
            plot_monthly_gpp(my_data, my_txt1, my_txt2, True, fig_path)

    #
    # COMPARE 2012 & 2015 RESULTS
    #
    to_compare = True
    if to_compare:
        gepisat_dir = os.path.join(
            os.path.expanduser("~"), "Desktop", "temp", "out")

        gep12_file = "GePiSaT_All_Data-26.txt"
        gep15_file = "GePiSaT_All_Data-30.txt"

        gep12_path = os.path.join(gepisat_dir, gep12_file)
        gep15_path = os.path.join(gepisat_dir, gep15_file)

        gep12_data = load_gepisat_data(gep12_path)
        gep15_data = load_gepisat_data(gep15_path)

        for station in sorted(list(set(gep12_data['station']))):
            if station in gep15_data['station']:
                gep12_st_idx = numpy.where(gep12_data['station'] == station)
                gep15_st_idx = numpy.where(gep15_data['station'] == station)

                gep12_st_data = gep12_data[gep12_st_idx]
                gep15_st_data = gep15_data[gep15_st_idx]

                gep12_st_gf_data = gapfill_months(gep12_st_data)
                gep15_st_gf_data = gapfill_months(gep15_st_data)

                my_meta = get_meta(flx_met_path, station)
                if my_meta.size == 1:
                    my_elv = my_meta['ele'][0]
                    my_veg = my_meta['classid'][0]
                    my_clim = my_meta['climateid'][0]
                    if my_elv == -9999:
                        my_txt1 = ("$\\mathrm{%s}$\n") % (station)
                        my_txt2 = (
                            "$\\mathrm{%s:}$ $\\mathrm{%s}$\n"
                            "$\\mathrm{%s:}$ $\\mathrm{%s}$") % (
                                "Veg", my_veg, "Clim", my_clim)
                    else:
                        my_txt1 = (
                            "$\\mathrm{%s}$\n$z=%0.1f$ m") % (
                                station, my_elv)
                        my_txt2 = (
                            "$\\mathrm{%s:}$ $\\mathrm{%s}$\n"
                            "$\\mathrm{%s:}$ $\\mathrm{%s}$") % (
                                "Veg", my_veg, "Clim", my_clim)
                else:
                    my_txt1 = ("$\\mathrm{%s}$\n") % (station)
                    my_txt2 = ""

                fig_file = "%s_monthly_gpp_comp.pdf" % (station)
                fig_path = os.path.join(gepisat_dir, "compare_figs", fig_file)
                plot_gepisat_gpp(gep15_st_gf_data,
                                 gep12_st_gf_data,
                                 my_txt1, my_txt2,
                                 True, fig_path)
            else:
                print("no pairing for station %s" % station)
