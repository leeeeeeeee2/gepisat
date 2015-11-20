# README
---------------
* LAST UPDATED: 2015-11-20
* TEAM: labprentice
* REPO: gepisat (private)

## Contents
--------------------
### check/
This directory holds consistency checks for the various P-model and GePiSaT versions.

* __check_r_version_beni.R__

### data/
This directory holds standard input data used for consistency checking.
The monthly data is for one of the grass FLUXNET sites in Oensingen, Switzerland for the year 2002.

* __mfsun_CH-Oe1_2002.txt__
    * Monthly sunshine fraction (unitless)
* __mprec_CH-Oe1_2002.txt__
    * Monthly precipitation (mm/month)
* __mtemp_CH-Oe1_2002.txt__
    * Monthly air temperature (Celcius)
* __mvapr_CH-Oe1_2002.txt__
    * Monthly vapor pressure (hPa)
* __mvpd_CH-Oe1_2002.txt__
    * Monthly vapor pressure deficit (kPa)

### doc/
This directory holds the documentation for the GePiSaT model.

* __gepisat_doc.pdf__
    * The current PDF format of the model documentation
* __gepisat_doc.tex__
    * The main LaTeX document
* __gepisat.bib__
    * The BibLaTeX file for documentation references
* __img/__
    * Contains EPS figures
* __tex/__
    * Contains LaTeX section files for documentation

### main_py/
This directory contains the main GePiSaT model Python code.

* __table_maker.py__
    * Data manager script for PostgreSQL database
    * Reads multiple input file types
    * Produces variable lists and data sets for db_setup.py
* __db_setup.py__
    * Setup script for PostgreSQL database
    * Creates tables with appropriate schemas
    * Populates tables with data
* __model.py__
    * Main model code
    * Performs monthly flux partitioning of NEE and PPFD observation pairs
    * Gap-fills PPFD observations
    * Calculates GPP based on flux partitioning parameters
    * Integrates PPFD and GPP to monthly totals
    * Reads monthly meteorological data from database
    * Estimates LUE

## main_r/
This directory contains the R code function definitions for the production model (P-model), which is based on some of the methods in GePiSaT.

* __pmodel.R__
    * Function definitions

## tools/
This directory holds a variety of analyzing, plotting and processing tools.

* __analysis/__
    * __lue_analyzer.R__
        * Analysis functions for GePiSaT light-use efficiency outputs
    * __summary_stats.R__
        * Analysis for optimizing GePiSaT dynamic parameterization
* __misc/__
    * __map_country.py__
        * Utility function to map continents with their associated country names and abbreviations
* __plotting/__
    * __plot_gapfill.R__
        * This script plots the monthly PPFD observations and the gap-filling product
    * __plot_gpp.R__
        * This script creates plots of GPP based on the LUE station text files
            * Three climate regions / forest regions for each climate region
            * Box and whisker plot of monthly GPP for each region
            * Calculates monthly averages of GPP +/- st deviation
            * Writes out results
    * __plot_lue.R__
        * This script processes the LUE data from monthly LUE files (output from __model.py__) and produces plots
    * __plot_outliers.R__
        * This script reads the observation and outlier-free datasets output by __model.py__, plots them and highlights the observation pairs identified as outliers
    * __plot_partitioning.R__
        * This script reads the observation and outlier-free datasets output by model.py, and plots the linear and hyperbolic partitioning
    * __plot_ts.R__
        * This script reads the data output from __timeseries.py__ to plot the gridded data associated with flux tower locations.
* __processing/__
    * __catfiles_lue.pl__
        * This script will concatenate all the files listed in the working directory that end with ".txt" into a single file
    * __file_handler-osx.pl__
        * This script is meant to speed up file handling in Mac OSX (or Linux).
            * Creates subdirectories if they do not already exist
            * Moves files into their corresponding subdirectory
    * __file_handler-win.pl__
        * This script is meant to speed up file handling in Windows.
            * Creates subdirectories if they do not already exist
            * Moves files into their corresponding subdirectory
    * __annual_gpp.py__
        * This script processes the GPP output from GePiSaT into monthly and annual totals
        * Interpolation of monthly GPP is performed for cases where one or two consecutive months are missing in an attempt to complete the time series; only complete time series of monthly GPP are integrated to annual totals
    * __cru_dat.py__
        * This script reads CRU TS 3.00 data file (i.e., .dat) containing the half-degree land-surface elevation data
    * __cru_netcdf.py__
        * This script reads NetCDF files (e.g., CRU TS 3.21), extracts variables of interest (e.g., monthly average max and min air temperature, monthly average vapor pressure), and writes the variables of interest to ASCII raster format
    * __gepisat_gpp.py__
        * This script reads the output from GePiSaT for gapfilling PPFD (aka ST-NAME-GF_YYYY-MM-01.txt), calculates the daily PPFD totals, finds the associated flux partitioning parameters (via summary_statistics.txt), and calculates daily GPP
    * __gepisat_nlsr.py__
        * This script performs non-linear least squared regression for the next-gen LUE model
    * __gepisat_ts.py__
        * This script connects to the GePiSaT database, reads observations, and saves them to file as a time series for analysis
    * __glas_netcdf.py__
        * This script reads the NetCDF file for the 2005 GLAS-derived canopy height based on the work of Simard et al., 2011.
    * __glas_tiff.py__
        * This script opens a TIFF raster file for processing
        * NOTE: pixel indexing starts in the top-left (NW) corner (similar to MODIS)
    * __grid_centroid.py__
        * Given a point location, i.e., lon and lat, this script finds the four grid centroids surrounding the point, calculates the distances between the point and centroids, and returns the closest centroid.
        * If calculated distances are equal, the nearest centroid defaults north and east
    * __lue_vars_grid.py__
        * This script creates netCDF files of gridded datasets required for the calculation of GePiSaT's predicted GPP, i.e.: ```GPP = phi_o * m' * fa * Iabs```
        * where:
            * ```phi_o``` is the intrinsic quantum efficiency
            * ```m'``` is the CO2 substrate limitation term
            * ```fa``` is the index of plant-available moisture
            * ```Iabs``` is absorbed light [```mol m^-2```]
        * The required gridded datasets are:
            * ```Iabs```: absorbed light, [```mol m^-2```]
            * ```VPD```: vapor pressure deficit, [```kPa```]
            * ```Gs```: photorespiratory compensation point, [```Pa```]
            * ```K```: Michaelis-Menten coefficient, [```Pa```]
            * ```ns```: relative viscosity of water, unitless
            * ```fa```: index of plant-available moisture, unitless
            * ```beta```: ratio of unit costs carboxylation:transpiration, unitless
        * The required gridded observations (inputs) are:
            * ```elv```: elevation AMSV, meters
            * ```evi```: monthly enhanced vegatation index
            * ```swdown```: daily downwelling shortwave radiation, W m^-2
            * ```tmp```: monthly mean daily air temperature, degrees C
            * ```vap```: monthly mean atmospheric vapor pressure, hPa
        * The sources of inputs are:
            * ```elv```: CRU TS3.00 [netCDF] (0.5 deg x 0.5 deg)
            * ```evi```: MODIS MYD13C2 & MOD13C2 [HDF] (0.05 deg x 0.05 deg)
            * ```swdown```: WATCH WFDEI [netCDF] (0.5 deg x 0.5 deg)
            * ```tmp```: CRU TS3.22 [netCDF] (0.5 deg x 0.5 deg)
            * ```vap```: CRU TS3.22 [netCDF] (0.5 deg x 0.5 deg)
    * __modis_hdf.py__
        * This script reads MODIS land-surface products (in HDF file format) and produces a CSV file for importing geo-referenced data into GIS.
            * includes assigning longitude and latitude to measurements based on the gridded data's resolution
        * This script can also resample the MODIS CGM (0.05 deg) grid to a lower resolution (e.g., 0.5 deg) and process EVI to file
    * __modis_hdf2ncdf.py__
        * This script reads MODIS land-surface products (in HDF file format), upscales the data to 0.5 degree resolution, and saves the data to netCDF file format
    * __noaa_netcdf.py__
        * This script reads the netCDF file of 0.5 degree resolution mean monthly soil moisture data based on [van den Dool, Huang & Fan (2003)](http://www.esrl.noaa.gov/psd/data/gridded/data.cpcsoil.html)
    * __raster_to_ts.py__
        * This script creates a time series based on a single pixel (or list of pixels) read from monthly ASCII raster at 0.5 degree resolution
    * __watch_netcdf.py__
        * This script reads a netCDF file (WATCH Forcing Data) and extracts data for the variable of interest for each grid point. The data extraction is then formatted for ASCII raster file output
        * __Note:__ WFDEI contains 67,209 pixels of terrestrial *observations* plus an additional 27,533 pixels over Antarctica for a total of 94,742 pixels
        
## Global ecosystem Production in Space and Time (GePiSaT)
------------------------------------------------------------------------------------
This project is aimed to develop a modeling system for global hindcasting and analysis of spatial and temporal patterns in terrestrial gross primary production (GPP). This system takes a simplistic approach which makes the best use of observational data (e.g., flux towers, meteorological stations, and remote-sensing satellites) while defensibly representing the principal ecophysiological processes that govern GPP.

The model is developed over three stages:

### Stage 1: Partitioning carbon flux data
The first modeling stage consists of partitioning high time-resolution CO2 flux data (i.e., net ecosystem exchange), F, using *in situ* photosynthetic photon flux density (PPFD) measurements, Q. Both F and Q observation pairs are obtainable from the free fair-use [FLUXNET](http://fluxnet.ornl.gov/) archive of eddy covariance flux towers around the world.

### Stage 2: Light-use efficiency modeling
The second modeling stage estimates the light-use efficiency (LUE) base on monthly aggregated GPP and monthly aggregated, gap-filled PPFD. An analysis of the empirical dependencies of LUE on vegetational and environmental factors are investigated in order to yield a simple predictive model of LUE.

### Stage 3: Global GPP prediction
The third and final modeling stage generates spatial fields of GPP based on remotely sensed reflectances and predictive LUE (obtained from stage 2).

### Model inputs
Stage 1 requires only observational pairs of net ecosystem exchange (NEE) and PPFD from flux tower data sets. These data are available in convenient pre-processed format from the [FLUXNET Synthesis Dataset](http://www.fluxdata.org/) listed under the [Free Fair-Use](http://www.fluxdata.org/Shared%20Documents/Policy_Free_Final.pdf) data policy.

Stage 2 begins with the gap-filling of flux tower PPFD observations. The current methodology utilizes [WATCH](http://www.eu-watch.org/) 0.5° x 0.5° gridded forcing data for the ERA Interim ([WFDEI](http://www.eu-watch.org/gfx_content/documents/README-WFDEI(1).pdf)). Daily shortwave solar radiation flux (Rs), W/m^2, is available in netCDF file format.

For the LUE modeling part of stage 2, a variety of additional variables are needed. For the basic LUE model: __LUE = GPP/(fAPAR PPFD)__, MODIS-based monthly 0.5° x 0.5° gridded enhanced vegetation index (EVI) data is used as a proxy for fAPAR (fractionally absorbed photosynthetically active radiation). There are two [MODIS](http://modis.gsfc.nasa.gov/) satellites ([Aqua and Terra](http://nsidc.org/data/modis/terra_aqua_differences/)) that have associated monthly 0.05° x 0.05° gridded data in HDF4 format, which can then be up-scaled (e.g., arithmetically averaged) to 0.5° resolution.

The next-generation formulation is __LUE = φo *m*__ where φo is the intrinsic quantum efficiency and *m* is the water and light use compensation coefficient.  In order to calculate part of *m*, additional variables are needed, including:

* Michaelis-Menten coefficient for Rubisco-limited photosynthesis (K), Pa
    * air temperature (Tc), °C
    * elevation (Elv), meters
* photorespiratory compensation point (Γ\*), Pa
    * air temperature (Tc), °C
* vapor pressure deficit (VPD), Pa
    * ambient vapor pressure (ea), Pa
    * air temperature (Tc), °C
* ambient CO2 concentration (*ca*), Pa
* viscosity of water (η), mPa s
    * air temperature (Tc), °C
    * elevation (Elv), meters

The temperature and pressure dependencies of K and Γ\* are well known (Bernacchi et al., 2001) as well as the temperature and pressure dependency of η (Huber et al., 2009). Elevation (Elv) can used to approximate the atmospheric pressure (Cavcar, 2000).  VPD may be calculated based on mean air temperature (to calculate saturated vapor pressure) and ambient vapor pressure (Abtew & Melesse, 2013). A global estimate of *ca*  may be based on the [globally averaged marine surface annual mean data](http://www.esrl.noaa.gov/gmd/ccgg/trends/global.html#global_data).

### Code structure
The main Python code has three files: __db_setup.py__, __table_maker.py__ and __model.py__ and makes use of a [PostgreSQL](http://www.postgresql.org/) database for storing and organizing all the observation data (i.e., NEE, PPFD, SWdown, fAPAR, Tc, VPD, *ca*,  Elv).

#### db_setup.py
This script interfaces with the GePiSaT PostgreSQL database to initialize the database tables and populate them with data.

#### table_maker.py
The purpose of this file is to convert the format of any source data (e.g., MODIS, WATCH, NOAA) to that which conforms to the GePiSaT database design. The GePiSaT database consists of three tables: met_data (station meta data), var_list (observation variable information), and data_set (time series of observations). Each station in met_data has an associated list of available variables (from var_list), each of which has an associated set of observation data (i.e., data_set).

#### model.py
This script performs the stages of the GePiSaT model (currently stages 1 and 2).


## References
-------------------
* Abtew, W. and A. Melesse (2013) "Vapor Pressure Calculation Methods," Evaporation and Evapotranspiration: Measurements and Estimations, Springer, New York.
* Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, Jr. and P. Long (2001) Improved temperature response functions for models of Rubisco-limited photosynthesis, *Plant, Cell and Environment*, vol. 24, pp. 253-259.
* Cavcar, M. (2000) *The International Standard Atmosphere (ISA)*, Anadolu University, Turkey.
* Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. Sengers, M. J. Assael, ... K. Miyagawa (2009) *J. Phys. Chem. Ref. Data*, vol. 38 (2), pp. 101-125.
