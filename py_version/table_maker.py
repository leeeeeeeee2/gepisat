#!/usr/bin/python
#
# table_maker.py
#
# VERSION 3.0
# LAST UPDATED 2017-01-13
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2017 Prentice Lab
#
# This file is part of the GePiSaT (Global ecosystem Production in Space and
# Time) model.
#
# GePiSaT is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# GePiSaT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GePiSaT.  If not, see <http://www.gnu.org/licenses/>.
#
# ---------
# citation:
# ---------
# Davis, T.W., B.D. Stocker, X.M.P. Gilbert, T.F. Keenan, H. Wang, B.J. Evans,
# and I.C. Prentice. The Global ecosystem Production in Space and Time
# (GePiSaT) Model of the terrestrial biosphere: Part 1 - Flux partitioning
# and gap-filling gross primary production. Geosci. Model Dev.
#
# ------------
# description:
# ------------
# This script reads raw observation data, e.g., flux tower and gridded
# satellite data, and outputs CSV files formatted for the PostgreSQL database
# tables.  Currently, fluxdata.org CSV files, WATCH (WFDEI) 0.5 degree
# satellite netCDF files, and MODIS 0.05 HDF files, CRU netcdf, GLAS geoTIFF,
# and CRU dat files are supported.
#
# Only two of the three tables are output for flux tower data
# (var_list and data_set).
#
# ----------
# changelog:
# ----------
# 01. Updated input file path and output file name for bitnami [13.06.17]
# 02. Fixed IOError output (filename is FILE, not f) [13.06.17]
# 03. Reading data to make data_set table [13.06.17]
# 04. Changed output for both var_list and data_set to be restricted to
#     core-vars only [13.06.18]
# 05. Changed DATA class to FLUXDATA, because it's specific to fluxdata.org
#     station data [13.06.18]
# 06. Abbreviated the list of core variables to those with qc flags [13.07.04]
# 07. Added core variable qc flag dictionary [13.07.04]
# 08. Implemented output for data-set based on only observations [13.07.04]
# 09. Added "re" module for string searching [13.07.04]
# 10. Added scipy.io for netcdf handling [13.08.09]
# 11. Changed class "LINE" to "VAR" [13.08.09]
# 12. Added "SWdown" to coreVars and variableUnits hashes [13.08.09]
# 13. Moved main to "process_flux" function call to give room for gridded data
#     processing [13.08.09]
# 14. Added process_grid function for WATCH WFDEI [13.08.19]
# 15. Updated VAR class getStation function
#     separated it based on variable type 'flux' or 'grid' [13.08.19]
# 16. Changed WFD grid station naming scheme from lat+lon combo to
#     3 char + 6 digits (e.g., "WFD000123") [13.08.27]
# --> changed "WFD" to "HDF" [13.09.23]
# 17. Reordered processing of grids; spatial first (row-major) to get var_list,
#     then temporal to get data_set [13.08.27]
# 18. Updated get_flux_station to search filename for station ID [13.08.27]
# 19. Individualized data_set output files [13.08.27]
# 20. Added varlist_flag (to run varlist output once only) [13.08.27]
# 21. Deleted lathash and lonhash (not used) [13.08.30]
# 22. Fixed station parts error (moved from var_flag condition) [13.08.30]
# 23. Updated directory names (i.e., Projects) [13.09.09]
# 24. Added pxl_lat in process_watch function [13.09.23]
# 25. Added FAPAR to VAR class [13.09.27]
# 26. Added numpy to modules list [13.09.27]
# 27. Added pyhdf.SD to modules list [13.09.27]
# 28. Created MODISDATA class [13.09.27]
# 29. Created process_modis() function [13.09.27]
# 30. Changed vars in FLUXDATA class to vari [13.09.27]
# --> vars is a built-in function name
# 31. New output handling in process_flux for dat and var files [13.10.14]
# 32. Abbreviated core vars dict but kept the numbering [13.10.16]
# 33. Added core variable check to FLUXDATA class [13.10.16]
# 34. Added VPD to core variables list in units of kPa [13.11.06]
# 35. Added add_one_month() function [13.11.06]
# 36. Created process_cru() function [13.11.06]
# 37. Created CRUDATA class [13.11.06]
# 38. Added RH100 (canopy height, m) to core variables list [13.11.20]
# 39. Created process_glas() function for canopy height [13.11.20]
# 40. Corrected errors with stationid naming for GLAS gridded data [13.11.21]
# --> Created GLASDATA class to perform the calculations
# 41. Changed GLAS canopy height source data to TIF raster [13.11.28]
# 42. Added Image to module list [13.11.28]
# 43. Renamed process_cru to process_cru_vpd [14.01.17]
# 44. New process_cru for single netCDF file [14.01.17]
# 45. Added new variables to VAR class [14.01.17]
# --> Tc (monthly mean temp, deg. C)
# --> Pre (monthly precip, mm)
# 46. Fixed calc_vpd in CRUDATA class [14.01.17]
# --> skip 89 years (not 90)
# 47. Added new variable to VAR class [14.01.22]
# --> Cld (monthly cloudiness, % [0.0-100.0])
# 48. Added new variable to VAR class [14.01.23]
# --> Elv (CRU TS3.00 elevation data)
# 49. Created process_cru_elv function [14.01.23]
# 50. Deleted RH100 from GePiSaT database [14.02.04]
# --> renumbered CO2 at varid 21
# 51. Created process_alpha() [14.02.04]
# 52. Fixed error in VPD calc [14.02.18]
# --> 237.3 not 273.3
# 53. Added get_monthly_cru & get_time_index functions [14.02.19]
# 54. Deleted CRUDATA class [14.02.19]
# 55. Added calculate_vpd function [14.02.19]
# 56. Updated process_cru_vpd function [14.02.19]
# 57. Updated process_cru function [14.02.19]
# 58. General updates to style and inline comments [14.04.16]
# 59. More housekeeping [14.09.03]
# 60. Improved initializationo of vpd in calculate_vpd() [14.09.03]
# 61. PEP8 style fixes [15.11.18]
# 62. Style fixes for Python 3 compatibility [17.01.08]
# 63. Updated calculate_vpd function to allow mean air temp [17.01.08]
# 64. Updated process cru elevation w/ os.path.joins [17.01.09]
# 65. Added logging statements [17.01.13]
# 66. Separated classes and functions into their own modules [17.01.13]
#
###############################################################################
# IMPORT MODULES
###############################################################################
import logging
import os

# from database.crudata import process_cru
from database.crudata import process_cru_elv
# from database.crudata import process_cru_vpd
from database.fluxdata import process_flux_2012
# from database.glasdata import process_glas
# from database.modisdata import process_modis
# from database.splashdata import process_alpha
from database.watchdata import process_watch


###############################################################################
# MAIN PROGRAM
###############################################################################
if __name__ == '__main__':
    # Create a root logger:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Instantiating logging handler and record format:
    root_handler = logging.StreamHandler()
    rec_format = "%(asctime)s:%(levelname)s:%(name)s:%(funcName)s:%(message)s"
    formatter = logging.Formatter(rec_format, datefmt="%Y-%m-%d %H:%M:%S")
    root_handler.setFormatter(formatter)

    # Send logging handler to root logger:
    root_logger.addHandler(root_handler)

    # Define data directories:
    home_dir = os.path.expanduser("~")
    data_dir = os.path.join(home_dir, "Data")
    flux_dir = os.path.join(data_dir, "flux_data")
    watch_dir = os.path.join(data_dir, "watch")
    modis_dir = os.path.join(data_dir, "modis")
    cru_dir = os.path.join(data_dir, "cru")
    glas_dir = os.path.join(data_dir, "glas")
    alpha_dir = os.path.join(data_dir, 'splash')
    out_dir = os.path.join(data_dir, "out")

    if False:
        # Process flux data:
        process_flux_2012(flux_dir)

        # Process WATCH data:
        watch_voi = 'SWdown'
        process_watch(watch_dir, watch_voi)

        # Process MODIS data:
        # modis_voi = "CMG 0.05 Deg Monthly EVI"
        # process_modis(modis_dir, modis_voi)

        # Process CRU TS data:
        process_cru_elv(cru_dir)
        # process_cru_vpd(cru_dir, out_dir)
        # process_cru("cld", cru_dir, out_dir)
