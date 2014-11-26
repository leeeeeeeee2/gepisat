# Enthought Canopy v 1.4
#
# daily_gpp.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-06-11 -- created
# 2014-09-17 -- last updated
#
# ------------
# description:
# ------------
# This script reads the output from GePiSaT for gapfilling PPFD,
# calculates the daily PPFD totals, finds the associated flux partitioning
# parameters, and calculates daily GPP.
#
# ----------
# changelog:
# ----------
# 01. updated for Linux (not Mac) and for v.25 output [14.06.25]
# 02. separated model_select=0 and added check for zeros [14.06.25]
# 03. added sort to ppfd file iteration [14.06.25]
# 04. added check for negative GPP (set equal to zero) [14.06.25]
# 05. updated function doc [14.09.17]
#
################################################################################
#### IMPORT MODULES ############################################################
################################################################################
import datetime
import glob
import numpy
import os.path
import re

################################################################################
#### GLOBAL VARIABLES ##########################################################
################################################################################
TIME_FORMAT = '%Y-%m-%d %H:%M:%S'

################################################################################
#### FUNCTIONS #################################################################
################################################################################
def add_one_day(dt0):
    """
    Name:     add_one_day
    Input:    datetime date
    Output:   datetime date
    Features: Adds one day to datetime
    """
    dt1 = dt0 + datetime.timedelta(days=1)
    return dt1

def get_flux_station(fileName):
        """
        Name:     get_flux_station
        Input:    string, file name w/ path (fileName)
        Output:   string, station name (sname)
        Features: Returns flux station name based on file name
        """
        # Initialize station name:
        sname = ""
        try:
            # Use regular expression search on filename
            # note: group(1) returns what is inside the search ()
            sname = re.search(
                '(\w{2}-\w{2,3})-', 
                os.path.basename(fileName)
                ).group(1)
        except AttributeError:
            print "Station name not found in file:", fileName
            #
        return sname

def calc_daily(sday, d):
    """
    Name:     calc_daily
    Input:    - datetime date, starting time stamp (sday)
              - dict, time stamped data (d)
    Output:   tuple, daily GPP and associated error (total_vals, total_errs)
    Features: Returns the daily totals (mol/m^2) of half-hourly calculated GPP  
              and associated errors (umol/m^2/s) for a given day
    Depends:  - add_one_day
              - simpson
    """
    # Get daily indexes:
    my_keys = [i for i in d.keys() if i >= sday and i <= add_one_day(sday)]
    #
    # Save the associated data for each of the keys:
    my_vals = numpy.array([])
    my_errs = numpy.array([])
    for k in numpy.sort(my_keys):
        my_vals = numpy.append(my_vals, [d[k][0],])
        my_errs = numpy.append(my_errs, [d[k][1],])
    #
    # Calculate the sum of the daily values (umol m-2):
    # NOTE: half-hourly data, dt = 30 mins = 1800 s
    total_vals = simpson(my_vals, 1800)
    total_errs = simpson(my_errs, 1800)
    #
    # Convert to moles:
    total_vals = (total_vals/1000000.)
    total_errs = (total_errs/1000000.)
    #
    return (total_vals, total_errs)

def calc_gpp(stat_array, my_dict):
    """
    Name:     calc_gpp
    Input:    - numpy nd.array, statistics (stat_array)
              - dict, half-hourly PPFD observations (my_dict)
    Output:   dict, half-hourly GPP and associated error (gpp_dict)
    Features: Returns a dictionary of half-hourly GPP and associated error based
              on flux partitioning of half-hourly PPFD; partitioning parameters
              are from a statistics file (i.e., output from GePiSaT model)---
              in the absence of model parameters, estimates for the hyperbolic 
              model or linear model with outliers removed are used
    """
    # Check the model selection for this month and station:
    # NOTE: there are five options:
    # 0. no model (use estimates)
    # 1. hyperbolic with observations
    # 2. hyperbolic with outliers removed
    # 3. linear with observations
    # 4. linear with outliers removed
    my_model = stat_array['model_select'][0]
    #
    # Initialize variables:
    gpp_dict = {}
    foo = 0
    foo_err = 0
    alpha = 0
    alpha_err = 0
    gpp = 0
    gpp_err = 0
    #
    if my_model == 0:
        # Try hyperbolic first:
        foo = stat_array['foo_est_ro_h'][0]
        foo_err = 0
        alpha = stat_array['alpha_est_ro_h'][0]
        alpha_err = 0
        #
        if foo == 0 or alpha == 0:
            # Use linear instead:
            alpha = stat_array['alpha_est_ro_l'][0]
            alpha_err = 0
            #
            # Calculate GPP and GPP err:
            for ts,ppfd in my_dict.iteritems():
                gpp = (alpha*ppfd)
                gpp_err = (ppfd*alpha_err)
                #
                # Check for negative GPP:
                if gpp < 0:
                    gpp = 0
                gpp_dict[ts] = (gpp, gpp_err)
        else:
            # Calculate GPP and GPP err:
            for ts,ppfd in my_dict.iteritems():
                gpp = (1.0*alpha*foo*ppfd)/(alpha*ppfd + foo)
                gpp_err = numpy.sqrt(
                    (
                        (
                            (foo*ppfd*(alpha*ppfd+foo) - alpha*foo*ppfd**2)/
                            (alpha*ppfd+foo)**2
                        )**2 * (alpha_err**2)
                    ) + (
                        (
                            (alpha*ppfd*(alpha*ppfd+foo) - alpha*foo*ppfd)/
                            (alpha*ppfd+foo)**2
                        )**2 * (foo_err**2)
                    )
                )
                #
                # Check for negative GPP:
                if gpp < 0:
                    gpp = 0
                gpp_dict[ts] = (gpp, gpp_err)
        #
    elif my_model == 1 or my_model == 2:
        # Hyperbolic
        if my_model == 1:
            foo = stat_array['foo_opt_obs_h'][0]
            foo_err = stat_array['foo_err_obs_h'][0]
            alpha = stat_array['alpha_opt_obs_h'][0]
            alpha_err = stat_array['alpha_err_obs_h'][0]
        elif my_model == 2:
            foo = stat_array['foo_opt_ro_h'][0]
            foo_err = stat_array['foo_err_ro_h'][0]
            alpha = stat_array['alpha_opt_ro_h'][0]
            alpha_err = stat_array['alpha_err_ro_h'][0]
        #
        # Calculate GPP & GPP_err
        for ts,ppfd in my_dict.iteritems():
            gpp = (1.0*alpha*foo*ppfd)/(alpha*ppfd + foo)
            gpp_err = numpy.sqrt(
                (
                    (
                        (foo*ppfd*(alpha*ppfd+foo) - alpha*foo*ppfd**2)/
                        (alpha*ppfd+foo)**2
                    )**2 * (alpha_err**2)
                ) + (
                    (
                        (alpha*ppfd*(alpha*ppfd+foo) - alpha*foo*ppfd)/
                        (alpha*ppfd+foo)**2
                    )**2 * (foo_err**2)
                )
            )
            #
            # Check for negative GPP:
            if gpp < 0:
                gpp = 0
            gpp_dict[ts] = (gpp, gpp_err)
        #
    elif my_model == 3 or my_model == 4:
        # Linear
        if my_model == 3:
            alpha = stat_array['alpha_opt_obs_l'][0]
            alpha_err = stat_array['alpha_err_obs_l'][0]
        elif my_model == 4:
            alpha = stat_array['alpha_opt_ro_l'][0]
            alpha_err = stat_array['alpha_err_ro_l'][0]
        #
        # Calculate GPP & GPP error:
        for ts, ppfd in my_dict.iteritems():
            gpp = (alpha*ppfd)
            gpp_err = (ppfd*alpha_err)
            #
            # Check for negative GPP:
            if gpp < 0:
                gpp = 0
            gpp_dict[ts] = (gpp, gpp_err)
    #
    return(gpp_dict)

def simpson(my_array, h):
    """
    Name:     simpson
    Input:    - numpy nd.array (my_array)
              - int, width of time step, seconds (h)
    Output:   float, daily integral of values (s)
    Features: Returns the numerical integral of values using Simpson's rule
    """
    n = len(my_array)
    s = my_array[0] + my_array[-1]
    #
    for i in xrange(1, n, 2):
        s += 4.0 * my_array[i]
    for j in xrange(2, n-1, 2):
        s += 2.0 * my_array[j]
    s = s * h / 3.0
    return s

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

################################################################################
#### MAIN ######################################################################
################################################################################
# Define directory for gap filled PPFD files
mac = 0
if mac:
    ppfd_dir = (
        '/Users/twdavis/Projects/gepisat/results/2002-06/attempt22b/gapfill/'
    )
    stats_dir = (
        '/Users/twdavis/Dropbox'
        '/Work/Imperial/flux/results/2002-06/summary_stats/'
    )
    out_dir = (
    	'/Users/twdavis/Dropbox'
    	'/Work/Imperial/flux/results/2002-06/gpp/daily_gpp/'
    )
else:
    ppfd_dir = (
        '/home/user/Projects/gepisat/results/2002-06/attempt25/gapfill/'
    )
    stats_dir = (
        '/home/user/Dropbox'
        '/Work/Imperial/flux/results/2002-06/summary_stats/'
    )
    out_dir = (
        '/home/user/Projects/gepisat/results/2002-06/attempt25/gpp/'
    )

# Open and read stats:
stats_file = glob.glob(stats_dir + "*-v25.txt")[0]
if stats_file:
    stats = numpy.loadtxt(
        fname=stats_file,
        dtype={
            'names': (
                'name','month',
                'foo_opt_obs_h','foo_err_obs_h',
                'foo_est_ro_h','foo_opt_ro_h', 'foo_err_ro_h',
                'alpha_opt_obs_h', 'alpha_err_obs_h',
                'alpha_est_ro_h', 'alpha_opt_ro_h', 'alpha_err_ro_h',
                'alpha_opt_obs_l', 'alpha_err_obs_l',
                'alpha_opt_ro_l', 'alpha_err_ro_l',
                'model_select'
            ),
            'formats' : (
                'S6', 'S10',
                'f4', 'f4',
                'f4', 'f4', 'f4',
                'f4', 'f4',
                'f4', 'f4', 'f4',
                'f4', 'f4',
                'f4', 'f4',
                'i4'
            )
        },
        delimiter=',',
        skiprows=1,
        usecols = (
            0, 1,
            6, 7,
            10, 11, 12,
            16, 17,
            20, 21, 22,
            26, 27,
            31, 32,
            102
        )
    )

# Open and read each gapfill file:
ppfd_files = glob.glob(ppfd_dir + "*-GF_*txt")
if ppfd_files:
    for my_file in numpy.sort(ppfd_files):
        my_data = numpy.loadtxt(
            fname=my_file,
            dtype={
                'names': ('timestamp','ppfd_gf'),
                'formats': ('S19', 'f4')
            },
            delimiter=',',
            skiprows=1,
            usecols = (0,2)
        )
        #
        # Create a dictionary for gap filled PPFD data:
        my_ppfd = {}
        for t in my_data:
            dt = datetime.datetime.strptime(t['timestamp'], TIME_FORMAT)
            my_ppfd[dt] = t['ppfd_gf']
        #
        # Save starting and ending times for this month:
        start_time = datetime.datetime.strptime(
            my_data[0]['timestamp'], 
            TIME_FORMAT
        )
        end_time = datetime.datetime.strptime(
            my_data[-1]['timestamp'], 
            TIME_FORMAT
        )
        #
        # Get station name (from filename):
        my_station = get_flux_station(my_file)
        my_month = '%s' % start_time.date()
        #
        # Create an output file for daily GPP:
        out_file = out_dir + my_station + '_daily_gpp.txt'
        if not os.path.isfile(out_file):
            headerline = 'Timestamp,GPP_mol_m2,GPP_err\n'
            writeout(out_file, headerline)
        #
        # Get station's stats:
        my_stats = stats[
            numpy.where(
                (stats['name'] == my_station) & (stats['month'] == my_month)
            )
        ]
        #
        # Calculate GPP for this station and this month:
        my_gpp = calc_gpp(my_stats, my_ppfd)
        #
        # Calculate daily totals (umol m-2 d-1):
        while start_time < end_time:
            (gpp_val, gpp_err) = calc_daily(start_time, my_gpp)
            #
            # Write to file:
            try:
                OUT = open(out_file, 'a')
                OUT.write(
                    '%s,%0.6f,%0.6f\n' % (
                        start_time.date(), 
                        gpp_val, 
                        gpp_err
                    )
                )
            except IOError:
                print "Error: cannot write to file: ", out_file
            else:
                OUT.close()
            #
            start_time = add_one_day(start_time)
        #
