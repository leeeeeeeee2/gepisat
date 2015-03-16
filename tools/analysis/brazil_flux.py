import datetime
import glob
import numpy
import matplotlib.pyplot as plt

#
# CHECK BRAZIL SITES
#
my_dir = '/home/user/Desktop/CD32_BRAZIL_FLUX_NETWORK_1174/data/MAN_K34_QAQC/'
my_files = glob.glob(my_dir + '*Avg_hour*csv')
my_data = numpy.loadtxt(
    fname=my_files[0],
    dtype={'names': ('y', 'n', 'hh', 'par', 'Fc', 'fcraw', 'NEE'),
           'formats' : ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4')},
    delimiter=',',
    skiprows=6,
    usecols=(0, 1, 2, 12, 22, 23, 26)
)

p0 = my_data.shape[0]
p1 = numpy.where((my_data['par'] != -9999.) & (my_data['Fc'] != -9999.))
p2 = numpy.where((my_data['par'] != -9999.) & (my_data['fcraw'] != -9999.))

print "PAR/Fc    :", len(p1[0]), "/", p0
print "PAR/fcraw :", len(p2[0]), "/", p0

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