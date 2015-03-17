#
# brazil_flux.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-03-16 -- created
# 2015-03-17 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script reads the LBA Brazilflux data (downloaded from ORNL website) to
# check data pairs (i.e., PPFD + NEE) within the nine sites.
#
###############################################################################
## IMPORT MODULES:
###############################################################################
import datetime
import glob
import numpy
import matplotlib.pyplot as plt

###############################################################################
## FUNCTIONS:
###############################################################################
def find_files(d):
    """
    Name:     find_files
    Input:    str, parent directory
    Output:   dict
    Features: Returns all the files from subdirs for Brazilflux
    """
    my_subdirs = ['MAN_K34', 'RON_FNS', 'SP_PDG', 
                  'STM_K77', 'TOC_BAN', 'PA_CAX', 
                  'RON_RJA', 'STM_K67', 'STM_K83', 
                  'MAN_K34_QAQC', 'RON_FNS_QAQC', 'SP_PDG_QAQC', 
                  'STM_K77_QAQC', 'TOC_BAN_QAQC', 'PA_CAX_QAQC', 
                  'RON_RJA_QAQC', 'STM_K67_QAQC', 'STM_K83_QAQC']
    my_search = '*Avg_hour*csv'
    my_return = {}
    #
    for sd in my_subdirs:
        my_pattern = d + sd + '/' + my_search
        my_files = glob.glob(my_pattern)
        if my_files:
            my_return[sd] = my_files[0]
    #
    return my_return

###############################################################################
## MAIN:
###############################################################################
mac = True
if mac:
    my_dir = '/Users/twdavis/Projects/data/flux/brazil_flux/data/'
else:
    my_dir = '/home/user/Desktop/CD32_BRAZIL_FLUX_NETWORK_1174/data/'
my_files = find_files(my_dir)
my_pairs = {}
for name in numpy.sort(my_files.keys()):
    my_file = my_files[name]
    my_data = numpy.loadtxt(
        fname=my_file,
        dtype={'names': ('y', 'n', 'hh', 'par', 'Fc', 'fcraw', 'NEE'),
               'formats' : ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4')},
        delimiter=',',
        skiprows=6,
        usecols=(0, 1, 2, 12, 22, 23, 26)
    )
    #
    p0 = my_data.shape[0]
    p1 = numpy.where((my_data['par'] != -9999.) & (my_data['Fc'] != -9999.))
    p2 = numpy.where((my_data['par'] != -9999.) & (my_data['fcraw'] != -9999.))
    p3 = numpy.where((my_data['par'] != -9999.) & (my_data['NEE'] != -9999.))
    #
    my_pairs[name] = (p0, p1[0].shape[0], p2[0].shape[0], p3[0].shape[0])

print "Station,PAR:Fc,PAR:Fcraw,PAR:NEE,TOTAL"
for pair in numpy.sort(my_pairs.keys()):
    print "%s,%d,%d,%d,%d" %(pair, 
                             my_pairs[pair][1], 
                             my_pairs[pair][2],
                             my_pairs[pair][3],
                             my_pairs[pair][0])

my_ts = [datetime.datetime(my_data['y'][i], 1, 1) for i in xrange(p0)]
for i in xrange(p0):
    my_ts[i] += datetime.timedelta(days=int(my_data['n'][i] - 1))
    my_ts[i] += datetime.timedelta(hours=int(my_data['hh'][i]))
my_ts = numpy.array(my_ts)

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.setp(ax1.get_xticklabels(), rotation=0, fontsize=14)
plt.setp(ax1.get_yticklabels(), rotation=0, fontsize=14)
ax1.plot(my_ts[p1], my_data['par'][p1], 'k-') 
ax1.set_ylabel('PAR', fontsize=16)
ax1.set_xlabel('Time', fontsize=16)

ax2 = fig.add_subplot(122)
plt.setp(ax2.get_xticklabels(), rotation=0, fontsize=14)
plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=14)
ax2.plot(my_data['par'][p2], my_data['Fc'][p2], 'ro')
ax2.set_ylabel('Fc', fontsize=16)
ax2.set_xlabel('PAR', fontsize=16)

plt.show()

my_data = numpy.loadtxt(
    fname=my_files[0],
    dtype={
        'names': (
            'year', 'day', 'hour', 'min', 'ta', 'taed', 'wd', 'wded', 
            'pressed', 'press', 'rg', 'rr', 'par', 'rpar', 'Rn', 'FG', 'wsed',
            'ws', 'H', 'Hraw', 'LE', 'Leraw', 'Fc', 'fcraw', 'CO2', 'sCO2',
            'NEE', 'NEEf', 'mrs', 'ust', 'rh', 'prec', 'H2O', 'FH2O', 'U', 
            'Ued	V', 'Ved', 'ee', 'ees', 'dpt', 'tsavg', 'eqtemp', 'abshu', 
            'slopee', 'radtop', 'rgs', 'rgsout', 'rgl', 'rglout', 'stdW', 
            'ang', 'Tau', 'zl', 'tprof1', 'tprof2', 'tprof3', 'tprof4', 
            'tprof5', 'tprof6', 'tprof7', 'tprof8', 'tprof9', 'tprof10', 
            'avgprofT', 'msoil1', 'msoil2', 'msoil3', 'msoil4', 'msoil5', 
            'msoil6', 'msoil7', 'msoil8', 'msoil9', 'msoil10', 'totaltet', 
            'pCO2_1', 'pCO2_2', 'pCO2_3', 'pCO2_4', 'pCO2_5', 'pCO2_6', 
            'pCO2_7', 'pCO2_8', 'pCO2_9', 'pCO2_10', 'avgsto', 'H2O1', 
            'H2O2', 'H2O3', 'H2O4', 'H2O5', 'H2O6', 'H2O7', 'H2O8', 'H2O9', 
            'H2O10', 'avgprofW', 'Wind1', 'Wind2', 'Wind3', 'Wind4', 'Wind5', 
            'WAvg', 'tsoil1', 'tsoil2', 'tsoil3', 'tsoil4', 'tsoil5'
        ),
        'formats' : (
            'i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4',
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 
            'f4', 'f4', 'f4', 'f4', 'f4', 'f4'
        )
    },
    delimiter=',',
    skiprows=3
)