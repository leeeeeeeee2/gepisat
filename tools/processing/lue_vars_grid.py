#!/usr/bin/python
#
# lue_vars_grid.py
#
# written by Tyler W. Davis
# Imperial College London
#
# 2015-02-24 -- created
# 2015-02-25 -- last updated
#
# ------------
# description:
# ------------
# This script creates netCDF files of gridded datasets required for the 
# calculation of GePiSaT's predicted GPP, i.e.:
#
#                     GPP = phi_o * m' * fa * Iabs
#
#               where:  phi_o is the intrinsic quantum efficiency
#                       m' is the CO2 substrate limitation term
#                       fa is the index of plant-available moisture
#                       Iabs is absorbed light [mol m^-2]
#
# The required gridded datasets are:
#   Iabs ... absorbed light, mol m^-2
#   VPD .... vapor pressure deficit, kPa
#   Gs ..... photorespiratory compensation point, Pa
#   K ...... Michaelis-Menten coefficient, Pa
#   ns ..... relative viscosity of water
#   fa ..... index of plant-available moisture
#   beta ... ratio of unit costs carboxylation:transpiration
#
# The required gridded observations (inputs) are:
#   elv .... elevation AMSV, meters
#   evi .... monthly enhanced vegatation index
#   swdown . daily downwelling shortwave radiation, W m^-2
#   tmp .... monthly mean daily air temperature, degrees C
#   vap .... monthly mean atmospheric vapor pressure, hPa
#
# The sources of inputs are:
#   elv .... CRU TS3.00 [netCDF] (0.5 deg x 0.5 deg)
#   evi .... MODIS MYD13C2 & MOD13C2 [HDF] (0.05 deg x 0.05 deg)
#   swdown . WATCH WFDEI [netCDF] (0.5 deg x 0.5 deg)
#   tmp .... CRU TS3.22 [netCDF] (0.5 deg x 0.5 deg)
#   vap .... CRU TS3.22 [netCDF] (0.5 deg x 0.5 deg)
#
# ----------
# changelog:
# ----------
#
###############################################################################
## IMPORT MODULES
###############################################################################
import datetime
import glob
import numpy
from scipy.io import netcdf
from pyhdf import SD

###############################################################################
## FUNCTIONS
###############################################################################
def calc_gstar(tair_file, my_month):
    """
    Name:     calc_gstar
    Input:    - str, CRU air temperature file (tair_file)
              - datetime.date, month of interest (my_month)
    Output:   numpy.ndarray, gamma-star, Pa
    Features: Returns the 360x720 gridded temperature-dependent 
              photorespiratory compensation point
              NODATA: -9999.0
    Depends:  get_monthly_cru
    Ref:      Bernacchi et al. (2001), Improved temperature response 
              functions for models of Rubisco-limited photosynthesis, 
              Plant, Cell and Environment, 24, 253--259.
    """
    # Define constants
    ERROR_VAL = -9999.0
    gs25 = 4.220  # Pa, assuming 25 deg C & 98.716 kPa)
    dha = 37830   # J/mol
    kR = 8.3145   # J/mol/K
    #
    # Get monthly mean air temperature, deg C
    tair_data = get_monthly_cru(tair_file, my_month, 'tmp')
    #
    # Get nodata indexes:
    nodata_idx = numpy.where(tair_data == ERROR_VAL)
    #
    gs_data = dha*(tair_data - 25.0)
    gs_data /= (298.15*kR*(tair_data + 273.15))
    gs_data = numpy.exp(gs_data)
    gs_data *= gs25
    #
    gs_data[nodata_idx] *= 0.0
    gs_data[nodata_idx] += ERROR_VAL
    #
    return gs_data

def calc_iabs(evi_file, swd_file):
    """
    Name:     calc_iabs
    Inputs:   - str, MODIS EVI file name (evi_file)
              - str, WATCH WFDEI shortwave radiation file (swd_file)
    Output:   numpy.ndarray, Iabs, mol m^-2
    Features: Returns 360x720 gridded monthly absorbed light
              NODATA: -9999.0
    Depends:  - get_monthly_watch
              - get_monthly_modis
              - upscale_evi
              - swd_to_ppfd
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    # Get monthly-averaged PPFD, mol m^-2:
    swd_data = get_monthly_watch(swd_file, 'SWdown')
    ppfd_data = swd_to_ppfd(swd_data)
    #
    # Get monthly-average EVI (used as fAPAR):
    evi_temp = get_monthly_modis(evi_file, 'CMG 0.05 Deg Monthly EVI')
    evi_data = upscale_evi(evi_temp)
    #
    # Find pixels without data:
    nodata_idx = numpy.where((evi_data == ERROR_VAL) | (ppfd_data == ERROR_VAL))
    #
    # Calculate Iabs:
    iabs_data = evi_data*ppfd_data
    iabs_data[nodata_idx] *= 0.0
    iabs_data[nodata_idx] += ERROR_VAL
    #
    return iabs_data

def calc_k(elv_file, tair_file, my_month):
    """
    Name:     calc_k
    Input:    - str, CRU elevation file (elv_file)
              - str, CRU air temperature file (tair_file)
              - datetime.date, month of interest (my_month)
    Output:   numpy.ndarray, K, Pa
    Features: Returns the 360x720 gridded temperature & pressure dependent 
              Michaelis-Menten coefficient
              NODATA: -9999.0
    Depends:  - get_monthly_cru
              - elv_to_patm
    Ref:      Bernacchi et al. (2001), Improved temperature response 
              functions for models of Rubisco-limited photosynthesis, 
              Plant, Cell and Environment, 24, 253--259.
    """
    # Define constants
    ERROR_VAL = -9999.0 
    kc25 = 39.97     # Pa, assuming 25 deg C & 98.716 kPa
    ko25 = (2.748e4) # Pa, assuming 25 deg C & 98.716 kPa
    dhac = 79430     # J/mol
    dhao = 36380     # J/mol
    kR = 8.3145      # J/mol/K
    kco = 2.09476e5  # ppm, US Standard Atmosphere
    #
    # Get monthly mean air temperature, deg C
    tair_data = get_monthly_cru(tair_file, my_month, 'tmp')
    #
    # Get atmospheric pressure from elevation data:
    elv_data = get_elv(elv_file)
    patm_data = elv_to_patm(elv_data)
    #
    # Get nodata indexes:
    nodata_idx = numpy.where((tair_data == ERROR_VAL) | (elv_data == ERROR_VAL))
    #
    # Calculate CO2 coefficient:
    vc = dhac*(tair_data - 25.0)
    vc /= (298.15*kR*(tair_data + 273.15))
    vc = numpy.exp(vc)
    vc *= kc25
    #
    # Calculated O2 coefficient:
    vo = dhao*(tair_data - 25.0)
    vo /= (298.15*kR*(tair_data + 273.15))
    vo = numpy.exp(vo)
    vo *= ko25
    #
    # Calculate Michaelis-Menton K:
    k = kco*(1e-6)*patm_data
    k /= vo
    k += 1.0
    k *= vc
    #
    k[nodata_idx] *= 0.0
    k[nodata_idx] += ERROR_VAL
    #
    return k

def calc_vpd(tair_file, vap_file, my_month):
    """
    Name:     calc_vpd
    Inputs:   - str, CRU air temperature file (tair_file)
              - str, CRU vapor pressure file (vap_file)
              - datetime.date, month of interest (my_month)
    Output:   numpy.ndarray, VPD, kPa
    Features: Returns 360x720 gridded monthly CRU-based vapor pressure deficit
              NODATA: -9999.0
    Depends:  get_monthly_cru
    Ref:      Abtew, W. and A. Melesse (2013), Vapor Pressure Calculation 
              Methods in "Evaporation and Evapotranspiration: Measurements and 
              Estimations," Springer, NY, pp. 53--54.
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    # Get monthly mean air temperature and vapor pressure
    tair_data = get_monthly_cru(tair_file, my_month, 'tmp')   # deg C
    vap_data = get_monthly_cru(vap_file, my_month, 'vap')     # hPa
    #
    # Get no data indexes:
    nodata_idx = numpy.where((tair_data == ERROR_VAL) | (vap_data == ERROR_VAL))
    #
    # Calculate saturated vapor pressure, kPa (es):
    # Eq. 5.1, Abtew & Melesse (2013)
    es_data = 17.27*tair_data
    es_data /= (tair_data + 237.3)
    es_data = numpy.exp(es_data)
    es_data *= 0.611
    #
    # Calculate vapor pressure deficit, kPa
    vpd_data = es_data - 0.1*vap_data
    vpd_data[nodata_idx] *= 0.0
    vpd_data[nodata_idx] += ERROR_VAL
    #
    return vpd_data

def elv_to_patm(elv):
    """
    Name:     elv_to_patm
    Input:    numpy.ndarray, elevation, m (elv)
    Output:   numpy.ndarray, atmospheric pressure, Pa (patm)
    Features: Returns the 360x720 gridded atmospheric pressure based on the 
              barometric formula
              NODATA: -9999.0
    Ref:      Allen et al. (1998), Crop evapotranspiration - Guidelines for 
              computing crop water requirements - FAO irrigation and drainage
              paper 56.
    """
    # Define constants:
    ERROR_VAL = -9999.0
    kPo = 101325   # standard atmosphere, Pa (Allen, 1973)
    kTo = 298.15   # base temperature, K (Prentice, unpublished)
    kL = 0.0065    # temperature lapse rate, K/m (Allen, 1973)
    kG = 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
    kR = 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
    kMa = 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
    #
    # Get no data indices:
    nodata_idx = numpy.where(elv == ERROR_VAL)
    #
    # Convert elevation to pressure, Pa:
    p_exp = (kG*kMa/(kR*kL))
    patm = -1.0*kL*elv
    patm /= kTo
    patm += 1.0
    patm = numpy.power(patm, p_exp)
    patm *= kPo
    #
    patm[nodata_idx] *= 0.0
    patm[nodata_idx] += ERROR_VAL
    #
    return patm

def get_cru_file(paths, voi):
    """
    Name:     get_cru_file
    Input:    - str, directory paths for CRU data files (paths)
              - str, variable of interest (voi)
    Output:   str OR list of file names
    Features: Returns the CRU TS file for given variable of interest
    """
    # Read through all files within the paths for voi:
    my_files = []
    for d in paths:
        temp_files = glob.glob(d + "*" + voi + "*.*")
        for f in temp_files:
            my_files.append(f)
    #
    if my_files:
        if len(my_files) > 1:
            print "Found duplicate files!"
        else:
            print "Found file!"
            my_files = my_files[0]
    else:
        print "No files found!"
    #
    return my_files

def get_elv(my_file):
    """
    Name:     get_elv
    Input:    str, input file (my_file)
    Output:   numpy.ndarray, elevation, m
    Features: Returns 360x720 gridded elevation data
              NODATA: -9999.0
              [0,0] = (-89.75, -179.75)
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    f = numpy.loadtxt(my_file)
    noval_idx = numpy.where(f == -999.0)
    f[noval_idx] = ERROR_VAL + 0.0*f[noval_idx]
    #
    return f

def get_foh_evi(lon, lat, data):
    """
    Name:     get_foh_evi
    Input:    - float, 0.5 degree pixel longitude, degrees (lon)
              - float, 0.5 degree pixel latitude, degrees (lat)
              - numpy nd.array, 0.5 degree data (data)
    Output:   numpy nd.array (my_grid_evi)
    Features: Returns array of one hundred EVI values at 0.05 degree resolution 
              for a single 0.5 degree resolution pixel
    Depends:  - get_foh_grid
              - grid_to_index
    """
    my_grid_evi = numpy.array([])
    my_grid_pnts = get_foh_grid(lon, lat)
    my_grid_indx = grid_to_index(my_grid_pnts)
    for indx_pair in my_grid_indx:
        x,y = indx_pair
        zval = data[y][x]
        # Note: MODIS missing data is -0.3
        # use NaN for invalid data (easy to remove for averaging):
        if zval == -0.3:
            zval = numpy.NaN
        my_grid_evi = numpy.append(my_grid_evi, [zval])
    #
    return my_grid_evi

def get_foh_grid(hdg_lon, hdg_lat):
    """
    Name:     get_foh_grid
    Input:    None.
    Output:   list, list of tuples (foh_grid)
    Features: Returns a list of one hundred 0.05 degree lon-lat pairs based 
              on the class's half degree (HDG) coordinates (i.e., hdg_lon, 
              hdg_lat)
    """
    # Initialize five one-hundreths grid:
    foh_grid = []
    #
    # Define half-degree and five one-hundredths resolutions:
    hdg_res = 0.5
    foh_res = 0.05
    #
    # Calculate the binding box at 0.5 deg:
    westing = hdg_lon - 0.5*hdg_res
    northing = hdg_lat + 0.5*hdg_res
    #
    # Initialize centroid offsetting for foh_grid:
    foh_lon = westing + 0.5*foh_res
    foh_lat = northing - 0.5*foh_res
    #
    # Iterate over the 10x10 box:
    for y in xrange(10):
        lat = foh_lat - y*foh_res
        for x in xrange(10):
            lon = foh_lon + x*foh_res
            foh_grid.append((lon, lat))
    #
    return foh_grid

def get_lat_lon(x,y,r):
    """
    Name:     get_lat_lon
    Input:    - int, x-index, e.g., 0--319 (x)
              - int, y-index, e.g., 0--719 (y)
              - float, resolution, e.g., 0.5 (r)
    Output:   tuple, latitude-longitude pair (lat, lon)
    Features: Returns lon-lat pair for an x-y index pair and pixel resolution, 
              where [0,0] = (-89.75, -179.75) 
    """
    # Offset lat, lon to pixel centroid
    lon = -180.0 + (0.5*r)
    lat = -90.0 + (0.5*r)
    #
    # Offset lat, lon based on pixel index
    lon = lon + (x*r)
    lat = lat + (y*r)
    #
    return (lat, lon)

def get_modis_file(paths, m):
    """
    Name:     get_modis_file
    Inputs:   - str, directory paths for MODIS files (paths)
              - datetime.date, month of interest (m)
    Output:   str OR list of file names
    Features: Returns MODIS MOD13 or MYD13 file for a given month
    """
    # Read through all files within the paths for voi:
    my_date = "A%d%03d" % (m.year, m.timetuple().tm_yday)
    my_files = []
    for d in paths:
        temp_files = glob.glob(d + "*" + my_date + "*.*")
        for f in temp_files:
            my_files.append(f)
    #
    if my_files:
        if len(my_files) > 1:
            print "Found duplicate files!"
        else:
            print "Found file!"
            my_files = my_files[0]
    else:
        print "No files found!"
    #
    return my_files

def get_modis_xy(lon, lat):
    """
    Name:     get_modis_xy
    Input:    - float, longitude (lon)
              - float, latitude (lat)
    Output:   tuple, x-y indices
    Features: Returns x and y indices for a given MODIS lon-lat pair, i.e.
              [0,0] = (89.975, -179.975) 
    """
    x = (lon + 180.0)/0.05 - 0.5
    y = (90.0 - lat)/0.05 - 0.5
    #
    return (int(x), int(y))

def get_monthly_cru(my_file, m, v):
    """
    Name:     get_monthly_cru
    Input:    - str, CRU netcdf file (my_file)
              - datetime.date, month of interest (m)
              - str, variable of interest (v)
    Output:   numpy.ndarray
    Features: Returns 360x720 monthly CRU TS dataset for a given month and 
              variable of interest (e.g., tmp and vap)
              NODATA: -9999.0
              [0,0] = (-89.75, -179.75)
    Depends:  get_time_index
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    # Open netCDF file for reading:
    f = netcdf.NetCDFFile(my_file, "r")
    #
    # Save the base time stamp:
    bt = datetime.date(1900,1,1)
    #
    # Read the time data as array:
    f_time = f.variables['time'].data
    #
    # Find the time index for the current date:
    ti = get_time_index(bt, m, f_time)
    #
    # Get the spatial data for current time:
    f_var = f.variables[v]
    f_noval = f_var.missing_value
    f_temp = f_var.data[ti]
    f_data = numpy.copy(f_temp)
    f_var = None
    f.close()
    #
    noval_idx = numpy.where(f_data == f_noval)
    f_data[noval_idx] *= 0.0
    f_data[noval_idx] += ERROR_VAL
    #
    return f_data

def get_monthly_modis(my_file, v):
    """
    Name:     get_monthly_modis
    Input:    - str, MODIS HDF file (my_file)
              - str, variable of interest (v)
    Output:   numpy.ndarray
    Features: Returns 360x720 monthly MODIS dataset for a given variable of 
              interest (i.e., )
              NODATA: -0.3
              [0,0] = ()
    """
    # CMG 0.05 Deg Monthly EVI
    #   FillValue: -3000
    #   scale_factor: 10000.0
    #   valid_range: [-2000, 10000]
    #   dimensions: 3600, 7200
    try:
        # Try to open file for reading:
        f = SD.SD(my_file)
    except:
        print "Unexpected error opening file", my_file
        f_data = None
    else:
        # Check that VOI is in the dataset:
        if v in f.datasets().keys():
            # Pull data from dataset into array:
            f_select = f.select(v)
            #f_fillval = f_select.getfillvalue()
            #f_range = f_select.getrange()
            f_calval = f_select.getcal()[0]
            f_offset = f_select.getcal()[2]
            f_temp = f.select(v).get()
            f_data = f_temp.astype(float)/f_calval + f_offset
        else:
            print "Could not open dataset", v
            f_data = None
            #
        # Close HDF file and return data:
        f.end()
    finally:
        return f_data

def get_monthly_watch(my_file, v):
    """
    Name:     get_monthly_watch
    Inputs:   - str, file name (my_file)
              - str, variable of interest (v)
    Output:   numpy.ndarray, monthly shortwave, W m^-2
    Features: Returns 360x720 monthly average of WATCH WFDEI daily data
              NODATA: -9999.0
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    # Open netCDF file as a Python object:
    f = netcdf.NetCDFFile(my_file, "r")
    if v in f.variables.keys():
        f_var = f.variables[v]
        f_max = f_var.actual_max
        f_min = f_var.actual_min
        f_temp = f_var.data
        #
        f_data = numpy.copy(f_temp[0])
        num_days = f_temp.shape[0]
        for i in xrange(1, num_days, 1):
            f_data += f_temp[i]
        f_data /= float(num_days)
        #
        # Find valid data:
        noval_idx = numpy.where((f_data < f_min) | (f_data > f_max))
        f_data[noval_idx] = ERROR_VAL + 0.0*f_data[noval_idx]
        #
        # Set very low values to zero:
        lowval_idx = numpy.where((f_data < 1e-6) & (f_data > f_min))
        f_data[lowval_idx] = 0.0*f_data[lowval_idx]
    #
    f.close()
    #
    return f_data    

def get_time_index(bt, ct, aot):
    """
    Name:     get_time_index
    Input:    - datetime.date, base timestamp (bt)
              - datetime.date, current timestamp (ct)
              - numpy.ndarray, array of days since base timestamp (aot)
    Output:   int, time index for given month
    Features: Finds the index in an array of CRU TS days for a given timestamp 
    """
    # For CRU TS 3.2, the aot is indexed for mid-month days, e.g. 15--16th
    # therefore, to make certain that ct index preceeds the index for the
    # correct month in aot, make the day of the current month less than
    # the 15th or 16th (i.e., replace day with '1'):
    ct = ct.replace(day=1)
    #
    # Calculate the time difference between ct and bt:
    dt = (ct - bt).days
    #
    # Find the first index of where dt would be in the sorted array:
    try:
        idx = numpy.where(aot > dt)[0][0]
    except IndexError:
        print "Month searched in CRU file is out of bounds!"
        idx = None
    else:
        if dt < 0:
            print "Month searched in CRU file is out of bounds!"
            idx = None
    finally:
        return idx

def get_watch_file(paths, m):
    """
    Name:     get_watch_file
    Inputs:   - str, directory paths for WATCH files (paths)
              - datetime.date, month of interest (m)
    Output:   str OR list of file names
    Features: Returns WATCH WFDEI file for a given month
    """
    my_date = "%d%02d" % (m.year, m.month)
    my_files = []
    for d in paths:
        temp_files = glob.glob(d + "*" + my_date + "*.*")
        for f in temp_files:
            my_files.append(f)
    #
    if my_files:
        if len(my_files) > 1:
            print "Found duplicate files!"
        else:
            print "Found file!"
            my_files = my_files[0]
    else:
        print "No files found!"
    #
    return my_files

def grid_to_index(grid):
    """
    Name:     grid_to_index
    Input:    list, list of lon-lat tuples (grid)
    Output:   list, list of tuples (foh_indices)
    Features: Returns a list of x-y indices based for a given list of 0.05 
              degree lon-lat pairs
    Depends:  get_modis_xy
    """
    foh_indices = []
    for grid_pair in grid:
        lon, lat = grid_pair
        x, y = get_modis_xy(lon, lat)
        foh_indices.append((x,y))
    #
    return foh_indices

def swd_to_ppfd(swd):
    """
    Name:     swd_to_ppfd
    Input:    numpy.ndarray, shortwave radiation (W m^-2)
    Output:   numpy.ndarray, PPFD (mol m^-2)
    Features: Convertes shortwave radiation to photosynthetic photon flux 
              density; preserved error value
    """
    # Define constants:
    ERROR_VAL = -9999.0
    kfFEC = 2.04   # from flux to energy conversion, umol/J (Meek et al., 1984)
    ksec = 8.64e4  # seconds in a day
    #
    ppfd  = numpy.copy(swd)
    good_idx = numpy.where(swd != ERROR_VAL)
    ppfd[good_idx] *= ksec    # W m^-2 to J m^-2
    ppfd[good_idx] *= kfFEC   # J m^-2 to umol m^-2
    ppfd[good_idx] *= (1e-6)  # umol m^-2 to mol m^-2
    #
    return ppfd

def upscale_evi(d):
    """
    Name:     upscale_evi
    Input:    numpy.ndarray, 0.05 deg MODIS EVI (d)
    Output:   numpy.ndarray (hdg_evi)
    Features: Returns 360x720 gridded EVI resampled from 0.05 deg MODIS EVI
              NODATA: -9999.0
    Depends:  - get_lon_lat
              - get_foh_evi
    """
    # Define error value:
    ERROR_VAL = -9999.0
    #
    # Initialize data array of floating points values:
    hdg_evi = numpy.zeros(shape=(360,720))
    #
    # Iterate through 0.5 deg lat-lon pairs:
    for y in xrange(360):
        for x in xrange(720):
            (lat, lon) = get_lat_lon(x, y, 0.5)
            #
            # Get 0.05 EVI values within this 0.5 cell:
            my_evi_data = get_foh_evi(lon, lat, d)
            good_idx = numpy.where(~numpy.isnan(my_evi_data))[0]
            #
            # Check that there's data in the array:
            if good_idx:
                # Calculate the average EVI for 0.5 pixel:
                ave_evi = my_evi_data[good_idx].mean()
                #
                # Adjust points outside threshold:
                if ave_evi < 0.0:
                    ave_evi = 0.0
                elif ave_evi > 1.0:
                    ave_evi = 1.0
                #
            else:
                # Assign global error value
                ave_evi = ERROR_VAL
                #
            hdg_evi[y,x] = ave_evi
    #
    return hdg_evi

###############################################################################
## MAIN PROGRAM
###############################################################################
# Define the directory paths to input datasets:
mac = False
if mac:
    cru_paths = ['/Users/twdavis/Projects/data/cru/cru_ts_3_22/',
                 '/Users/twdavis/Projects/data/cru/cru_ts_3_00/']
    modis_paths = ['/Users/twdavis/Projects/data/modis/vi_cgm_monthly/aqua/',
                   '/Users/twdavis/Projects/data/modis/vi_cgm_monthly/terra/']
    watch_paths = ['/Users/twdavis/Projects/data/watch/SWdown_daily_WFDEI/',]
else:
    cru_paths = ['/usr/local/share/database/cru/',]
    modis_paths = ['/usr/local/share/database/modis/evi/aqua/',
                   '/usr/local/share/database/modis/evi/terra/']
    watch_paths = ['/usr/local/share/database/watch/netcdf/swdown/',]

# Initialize month
my_month = datetime.date(2002, 7, 1)

# Find the necessary input files:
elv_file = get_cru_file(cru_paths, 'elv')
tair_file = get_cru_file(cru_paths, 'tmp')
vap_file = get_cru_file(cru_paths, 'vap')
evi_file = get_modis_file(modis_paths, my_month)
swd_file = get_watch_file(watch_paths, my_month)
#
# CHECK TO MAKE SURE YOU HAVE THE RIGHT FILE FOR EACH VARIABLE!
#

# Extract raw observation data:
elv_data = get_elv(elv_file)
tair_data = get_monthly_cru(tair_file, my_month, 'tmp')

# Calculate variables:
iabs_data = calc_iabs(evi_file, swd_file)           # mol m^-2
vpd_data = calc_vpd(tair_file, vap_file, my_month)  # kPa
gs_data = calc_gstar(tair_file, my_month)           # Pa
k_data = calc_k(elv_file, tair_file, my_month)      # Pa
