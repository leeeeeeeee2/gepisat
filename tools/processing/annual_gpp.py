#!/usr/bin/python
#
# annual_gpp.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-01-10 -- created
# 2015-01-10 -- last updated
#
# -----------
# description
# -----------
# This script processes the GPP output from GePiSaT into monthly and annual
# totals. Interpolation of monthly GPP is performed for cases where one or two
# consecutive months are missing in an attempt to complete the time series; 
# only complete time series of monthly GPP are integrated to annual totals. 

###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import numpy

###############################################################################
## FUNCTIONS
###############################################################################
def print_annual_gpp(f, d):
    """
    Name:     print_annual_gpp
    Input:    - str, output file name
              - numpy.ndarray, structured monthly GPP array
    Output:   None.
    Features: Saves annual gpp data to file
    Depends:  - writeout
    """
    # Create output file:
    my_header = 'Station,Year,GPP_mol_m2\n'
    writeout(f, my_header)
    #
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    #
    # Write each line to file:
    for j in xrange(d.shape[0]):
        month_of_gpp = numpy.array(list(d[month_names][j]))
        if ~numpy.isnan(month_of_gpp).any():
            annual_gpp = month_of_gpp.sum()
            line = tuple([d['station'][j], d['year'][j], annual_gpp])
            try:
                OUT = open(f, 'a')
                OUT.write('%s,%i,%f\n' % line)
            except IOError:
                print "Error: cannot write to file: ", f
            else:
                OUT.close()

def print_monthly_gpp(f, d):
    """
    Name:     print_monthly_gpp
    Input:    - str, output file name
              - numpy.ndarray, structured monthly GPP array
    Output:   None.
    Features: Saves monthly gpp data to file
    Depends:  - writeout
    """
    # Create output file:
    my_header = 'Station,Year,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec\n'
    writeout(f, my_header)
    #
    # Write each line to file:
    for line in d:
        line = tuple(line)
        try:
            OUT = open(f, 'a')
            OUT.write('%s,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n' % line)
        except IOError:
            print "Error: cannot write to file: ", f
        else:
            OUT.close()

def writeout(f, d):
    """
    Name:     writeout
    Input:    - str, file name with path (t)
              - str, data to be written to file (d)
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

###############################################################################
## MAIN PROGRAM 
###############################################################################
mac = 0
if mac:
    lue_dir = '/Users/twdavis/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    out_dir = '/Users/twdavis/Desktop/'
else:
    lue_dir = '/home/user/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
    out_dir = '/home/user/Desktop/'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load LUE data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lue_files = glob.glob(lue_dir + '*26.txt')
lue_file = lue_files[0]

lue_data = numpy.loadtxt(
    lue_file, 
    delimiter=",", 
    skiprows=1,
    usecols=(0,1,2,4,5),
    dtype={'names' : ('station', 'time', 'gpp', 'fapar', 'ppfd'),
           'formats' : ('S6', 'O', 'f4', 'f4', 'f4')},
    converters={0: numpy.str,
                1: lambda x: datetime.datetime.strptime(x, '%Y-%m-%d'),
                2: numpy.float,
                4: numpy.float,
                5: numpy.float}
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build data matrix for holding monthly GPP per station, per year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_stations = [i for i in numpy.sort(list(set(lue_data['station'])))]
all_years = [2002, 2003, 2004, 2005, 2006]
# * note monthly names match datetime.datetime.strftime('%b')
my_names = tuple(['station', 'year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
my_formats = tuple(['S6', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
                    'f4', 'f4', 'f4', 'f4'])
i = 0
for station in all_stations:
    for year in all_years:
        my_line = tuple([station, year]) + tuple(numpy.repeat(numpy.nan, 12))
        if i == 0:
            gpp_monthly_all = numpy.array(my_line, 
                                          dtype={'names' : my_names,
                                                 'formats' : my_formats},
                                          ndmin=1)
            i += 1
        else:
            gpp_monthly_temp = numpy.array(my_line, 
                                           dtype={'names' : my_names,
                                                  'formats' : my_formats},
                                           ndmin=1)
            gpp_monthly_all = numpy.append(gpp_monthly_all, 
                                           gpp_monthly_temp, 
                                           axis=0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Save monthly GPP for each station for each year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for lue in lue_data:
    # Station row matches station name and year
    station_row = numpy.where(
        (gpp_monthly_all['station'] == lue['station']) & 
        (gpp_monthly_all['year'] == lue['time'].year)
        )[0][0]
    # Station column matches month
    station_col = lue['time'].strftime('%b')
    gpp_monthly_all[station_row][station_col] = lue['gpp']

all_out = out_dir + 'GePiSaT_GPP-mo_All.txt'
print_monthly_gpp(all_out, gpp_monthly_all)

all_annual_out = out_dir + 'GePiSaT_GPP-an_All.txt'
print_annual_gpp(all_annual_out, gpp_monthly_all)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter empty months
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
gpp_monthly_filtered = numpy.array(gpp_monthly_all[0],
                                   dtype={'names' : my_names,
                                          'formats' : my_formats},
                                   ndmin=1)
for j in xrange(1, gpp_monthly_all.shape[0]):
    month_of_gpp = numpy.array(list(gpp_monthly_all[month_names][j]))
    if ~numpy.isnan(month_of_gpp).all():
        temp = numpy.array(gpp_monthly_all[j],
                           dtype={'names' : my_names,
                                  'formats' : my_formats},
                           ndmin=1)
        gpp_monthly_filtered = numpy.append(gpp_monthly_filtered, temp, 0)

filtered_out = out_dir + 'GePiSaT_GPP-mo_Filtered.txt'
print_monthly_gpp(filtered_out, gpp_monthly_filtered)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gap-fill single missing months
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gpp_monthly_gapfill = numpy.copy(gpp_monthly_filtered)
for k in xrange(0, gpp_monthly_gapfill.shape[0]):
    for m in xrange(len(month_names)):
        my_month = month_names[m]
        if numpy.isnan(gpp_monthly_gapfill[k][my_month]):
            # Check preceeding and next monthly GPP
            if m == 11:
                m_next = month_names[0]
                m_before = month_names[10]
            elif m == 0:
                m_next = month_names[1]
                m_before = month_names[11]
            else:
                m_next = month_names[m + 1]
                m_before = month_names[m - 1]
            #
            # Check that GPP of m_next and m_before are not NANs
            gpp_next = gpp_monthly_gapfill[k][m_next]
            gpp_before = gpp_monthly_gapfill[k][m_before]
            if ~numpy.isnan(gpp_next) and ~numpy.isnan(gpp_before):
                gpp_gapfill = 0.5*(gpp_next + gpp_before)
                gpp_monthly_gapfill[k][my_month] = gpp_gapfill

gapfill_out = out_dir + 'GePiSaT_GPP-mo_1-Gapfill.txt'
print_monthly_gpp(gapfill_out, gpp_monthly_gapfill)

gapfill_annual_out = out_dir + 'GePiSaT_GPP-an_1-Gapfill.txt'
print_annual_gpp(gapfill_annual_out, gpp_monthly_gapfill)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gap-fill two consecutive missing months
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k in xrange(0, gpp_monthly_gapfill.shape[0]):
    for m in xrange(len(month_names)):
        # Check get current, next and previous months:
        m_now = month_names[m]
        if m == 11:
            m_pre_2 = month_names[9]
            m_pre_1 = month_names[10]
            m_nxt_1 = month_names[0]
            m_nxt_2 = month_names[1]
        elif m == 0:
            m_pre_2 = month_names[10]
            m_pre_1 = month_names[11]
            m_nxt_1 = month_names[1]
            m_nxt_2 = month_names[2]
        else:
            m_nxt_1 = month_names[m + 1]
            m_pre_1 = month_names[m - 1]
            #
            if m == 1:
                m_pre_2 = month_names[11]
                m_nxt_2 = month_names[3]
            elif m == 10:
                m_pre_2 = month_names[8]
                m_nxt_2 = month_names[0]
            else:
                m_pre_2 = month_names[m - 2]
                m_nxt_2 = month_names[m + 2]
        #
        # Get current, next and previous GPP:
        gpp_pre_2 = gpp_monthly_gapfill[k][m_pre_2]
        gpp_pre_1 = gpp_monthly_gapfill[k][m_pre_1]
        gpp_now = gpp_monthly_gapfill[k][m_now]
        gpp_nxt_1 = gpp_monthly_gapfill[k][m_nxt_1]
        gpp_nxt_2 = gpp_monthly_gapfill[k][m_nxt_2]
        #
        # Check if NaN                         #  MO1 GPP1  # MO1 GPP2
        isnan_pre2 = numpy.isnan(gpp_pre_2)    #1  a   A    #
        isnan_pre1 = numpy.isnan(gpp_pre_1)    #2  x1  y1   #1  a   A
        isnan_now = numpy.isnan(gpp_now)       #3  x2  y2   #2  x1  y1
        isnan_nxt1 = numpy.isnan(gpp_nxt_1)    #4  b   B    #3  x2  y2
        isnan_nxt2 = numpy.isnan(gpp_nxt_2)    #            #4  b   B
        #
        # 1. Current and previous NaNs
        if isnan_now and isnan_pre1 and not isnan_nxt1 and not isnan_pre2:
            A = gpp_pre_2
            B = gpp_nxt_1
            y1 = A + (B-A)*(1.0)/(3.0)
            y2 = A + (B-A)*(2.0)/(3.0)
            #
            gpp_monthly_gapfill[k][m_pre_1] = y1
            gpp_monthly_gapfill[k][m] = y2
            #
        # 2. Current and next NaNs
        if isnan_now and isnan_nxt1 and not isnan_nxt2 and not isnan_pre1:
            A = gpp_pre_1
            B = gpp_nxt_2
            y1 = A + (B-A)*(1.0)/(3.0)
            y2 = A + (B-A)*(2.0)/(3.0)
            #
            gpp_monthly_gapfill[k][m_nxt_1] = y2
            gpp_monthly_gapfill[k][m] = y1
            #
gapfill2_out = out_dir + 'GePiSaT_GPP-mo_2-Gapfill.txt'
print_monthly_gpp(gapfill2_out, gpp_monthly_gapfill)

gapfill2_annual_out = out_dir + 'GePiSaT_GPP-an_2-Gapfill.txt'
print_annual_gpp(gapfill2_annual_out, gpp_monthly_gapfill)