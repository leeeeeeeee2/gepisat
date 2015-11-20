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
        * This script reads the observation and outlier-free datasets output by __model.py__, and plots the linear and hyperbolic partitioning
    * __plot_ts.R__
        * This script reads the data output from __gepisat_ts.py__ to plot the gridded data associated with flux tower locations.
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
