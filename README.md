# README
---------------
* LAST UPDATED: 2014-10-28
* TEAM: labprentice
* REPO: gepisat (private)

## Contents
--------------------
### doc/
This directory holds the current PDF documentation for the GePiSaT model.

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

### Model Inputs
Stage 1 requires only observational pairs of net ecosystem exchange (NEE) and PPFD from flux tower data sets. These data are available in convenient pre-processed format from the [FLUXNET Synthesis Dataset](http://www.fluxdata.org/) listed under the [Free Fair-Use](http://www.fluxdata.org/Shared%20Documents/Policy_Free_Final.pdf) data policy. 

Stage 2 begins with the gap-filling of flux tower PPFD observations. The current methodology utilizes [WATCH](http://www.eu-watch.org/) 0.5° x 0.5° gridded forcing data for the ERA Interim ([WFDEI](http://www.eu-watch.org/gfx_content/documents/README-WFDEI(1).pdf)). Daily shortwave solar radiation flux (Rs), W/m^2, is available in netCDF file format.

For the LUE modeling part of stage 2, a variety of additional variables are needed. For the basic LUE model: __LUE = GPP/(fAPAR PPFD)__, MODIS-based monthly 0.5° x 0.5° gridded enhanced vegetation index (EVI) data is used as a proxy for fAPAR (fractionally absorbed photosynthetically active radiation). There are two [MODIS](http://modis.gsfc.nasa.gov/) satellites ([Aqua and Terra](http://nsidc.org/data/modis/terra_aqua_differences/)) that have associated monthly 0.05° x 0.05° gridded data in HDF4 format, which can then be up-scaled (e.g., arithmetically averaged) to 0.5° resolution. 

The next-generation formulation is __LUE = φo *m*__ where φo is the intrinsic quantum efficiency and *m* is the water and light use compensation coefficient.  In order to calculate part of *m*, additional variables are needed, including: 

* Michaelis-Menten coefficient for Rubisco-limited photosynthesis (K), Pa
* photorespiratory compensation point (Γ\*), Pa
* vapor pressure deficit (VPD), Pa
* ambient CO2 concentration (*ca*), Pa

The temperature dependencies of K and Γ\* are well known (Bernacchi et al., 2001). VPD may be calculated based on mean air temperature (to calculate saturated vapor pressure) and ambient vapor pressure (Abtew & Melesse, 2013). A global estimate of *ca*  may be based on the [globally averaged marine surface annual mean data](http://www.esrl.noaa.gov/gmd/ccgg/trends/global.html#global_data). 

## References
-------------------
* Abtew, W. and A. Melesse (2013) "Vapor Pressure Calculation Methods," Evaporation and Evapotranspiration: Measurements and Estimations, Springer, New York.
* Bernacchi, C. J., E. L. Singsaas, C. Pimentel, A. R. Portis, Jr. and P. Long (2001) Improved temperature response functions for models of Rubisco-limited photosynthesis, *Plant, Cell and Environment*, vol. 24, pp. 253-259.


