# RStudio 0.98.495 / 0.98.501
#
# lue_analyzer.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-02-06 -- created
# 2014-11-16 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script is for analyzing the new light use efficiency results for 
# GePiSaT v 2.0:
#
# A_i,j = phi * m_i,j * PPFD_i,j
# where:
#   A = monthly GPP
#   phi = the intrinsic quantum efficiency
#   m = the new LUE efficiency compensation coefficients
#   PPFD = flux tower observation
#   the subsript i is the time index and j is the station index
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 00. created [14.02.06]
# 01. created plots for basic lue, lue w/ alpha, and new lue [14.02.07]
# 02. added ablines to plots [14.02.07]
#     NO GOOD -- REMOVED
# 03. added legend for LUE [14.02.07]
# 04. fixed error in estimate.sqba function [14.02.09]
# 05. sensitivity plots for new lue model [14.02.11]
# 06. changed estimate.sqba function [14.02.09]
# --> range from 0--2000 Pa for VPD (instead of 4000 Pa)
# 07. added color to plots for vpd bins [14.02.12]
# --> plot has been extracted to function, plot_vpd()
# 08. extended range in estimate.sqba for chi [14.02.12]
# --> now between 0.5 and 0.9
# 09. boxplot comparisons of content variables [14.02.17]
# 10. created sigmoid function for temperature [14.02.17]
# 11. updated calc_gstar & calc_k [14.02.27]
# --> kc, ko, gstar are constant; Oi changes with elevation
# 12. updated calc_m [14.02.27]
# 13. added lue_conversions() function [14.02.27]
# 14. added geom.mean() function [14.02.27]
# 15. removed unused functions and cleaned up some of the code [14.02.27]
# 16. created plot_mbD() function [14.02.27]
# 17. added plot_temp() function [14.03.07]
# 18. added plot_alpha() functin [14.03.07]
# 19. added plot_temp_alpha() function [14.03.07]
# 20. fixed temperature colors for 20-25 deg. [14.03.20]
# 21. replaced plot with errbar (Hmisc package) [14.03.24]
# 22. added plot_r2 function [14.03.30]
# 23. separated content data filter from lue_conversions [14.03.30]
# --> it's its own function filter_contents [14.03.30]
# 24. added lue.analyzer & plot.lue.analyzer [14.04.18]
# 25. added flag to plot.lue.analyzer [14.04.19]
# 26. added monthly lue analyzer [14.05.14]
# 27. added filter_stations function [14.05.22]
# 28. added vegetation class to monthly.lue.analyzer function [14.05.22]
# 29. added bin fields (veg height & temp) to monthly.lue.analyzer [14.05.28]
# 30. updated plot vpd [14.06.03]
# 31. updated plot alpha [14.06.10]
# 32. added get_meta function [14.06.10]
# 33. updated Michaelis-Menten constants from Bernacchi et al., 2001 [14.06.27]
# 34. plot_xy [14.07.24]
# 35. updates to plot_xy function [14.07.28]
# 36. added atmospheric pressure to lue.analyzer [14.08.04]
# 37. created plot_lue_xy() function [14.08.04]
# 38. added station average and monthly LUE results output [14.10.07]
# 39. created plot_xy_vpd and plot_xy_temp functions [14.10.16]
# 40. updated plot_alphs [14.10.22]
# --> added weights to linear regression (GPP_err)
# --> replaced r-squared with adj.r.squared
# 41. updated plot_alpha [14.10.27]
# --> new alpha ranges
# --> added abline / bg color for legend
# 41. modified palette for color plots [14.10.27]
# 42. updated filter_stations function [14.11.03]
# --> removes crops and wetlands
# 43. created climate_class function [14.11.03]
# 44. created modeled_lue & bp_modeled_lue functions [14.11.03]
# 45. created bp_lue_dr2 and modeled_lue_alt functions [14.11.16]
#
# To crop ps to pdf:
# 1. ps2pdf -dPDFSETTINGS=/prepress [file_name].ps
# 2. pdfcrop --margins=25 [file_name].pdf
#
# ~~~~~
# todo:
# ~~~~~
# CHECK calc_m, calc_mp, lueResid, binned.lue.regression FUNCTIONS
# WHAT AM I DOING?!
# -- CHECK WHAT PARAMETERS ARE BEING ESTIMATED (B/A OR M)?
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### LIBRARIES ################################################################
# /////////////////////////////////////////////////////////////////////////////
library(Hmisc)         # for errorbars (e.g., plot_alpha)
library(lattice)       # for nice histograms (e.g., hist_modeled_lue)
library(minpack.lm)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: filter_contents
# *
# * Input: data frame
#
# * Return: data frame
#
# * Features: This function filters rows from the input data frame
#             where fAPAR, Tc, or VPD are below a threshold (~0)
# ************************************************************************
filter_contents <- function(my.obj){
  my.filt <- matrix(ncol=(dim(my.obj)[2]),nrow=1)
  my.filt <- as.data.frame(my.filt)
  names(my.filt) <- names(my.obj)
  my.filt <- my.obj[
    which(my.obj$fAPAR >= 0 & my.obj$Tc.deg_C >= 0 & my.obj$VPD.kPa > 0),
    ]
  my.filt
}

# ************************************************************************
# * Name: calc_gstar
# *
# * Input: air temperature (tc), degrees C
# *
# * Return: Photorespiratory compensation point (Gs), Pa
# *
# * Features: This function calculates the temperature-dependent 
# *           Gamma_star based on Bernacchi et al., 2001
# *
# * Note: The Bernacchi et al., 2001 constant was converted from mole
# *       fraction to partial pressure based on the assumed location of
# *       experimentation at the University of Illinois.
# * 
# ************************************************************************
calc_gstar <- function(tc){
  # Define constants:
  Gs25.ppm <- 42.75  # Bernacchi et al., 2001
  dHa <- 37830       # J/mol
  R <- 8.3145        # J/mol/K
  #
  # Convert ppm to partial pressure:
  Gs25.pa <- 4.220     # Bernacchi et al., 2001 (assuming 25 C & 98.716 kPa)
  #
  # Calculate Gs, Pa:
  Gs <- Gs25.pa * exp(dHa * (tc-25.0)/(298.15*R*(tc+273.15)))
  Gs
}

# ************************************************************************
# * Name: calc_k
# *
# * Input: air temperature (tc), degrees C
# *        atm. pressure (patm), Pa
# *
# * Return: Michaelis-Menten photosynthetic rate (K), Pa
# *
# * Features: This function calculates the temperature and pressure 
# *           dependent Michaelis-Menten coefficient, K, based on the
# *           constants defined in Bernacchi et al., 2001. The US 
# *           Standard Atmosphere is referenced for the atmospheric
# *           oxygen concentration (ppm).
# *
# * Note: The Bernacchi et al., 2001 constants were converted from mole
# *       fraction to partial pressure based on the assumed location of
# *       experimentation at the University of Illinois.
# *
# ************************************************************************
calc_k <- function(tc, patm){
  # Michaelis-Menten photosynthetic rate
  # Kubien et al., 2008
  #
  # Required functions:
  # Kc :: Michaelis-Menton CO2 constant, Pa
  # Ko :: Michaelis-Menton O2 constant, Pa
  Kc <- function(Tc){
    # Tc :: air temperature, deg. C
    # P  :: atmospheric pressure, Pa
    Kc25.ppm <- 404.9; # Bernacchi et al., 2001
    Kc25.pa <- 39.97   # Bernacchi et al., 2001 (assuming 25 C & 98.716 kPa)
    dHa <- 79.43; # kJ/mol
    R <- 0.0083145; # kJ/mol/K
    Kc <- Kc25.pa * exp(dHa * (Tc-25) / (298.15*R*(Tc+273.15)))
    Kc
  }
  Ko <- function(Tc){
    # Tc :: air temperature, deg. C
    # P  :: atmospheric pressure, Pa
    Ko25.ppm <- 278400; # Bernacchi et al., 2001
    Ko25.pa <- 27480;   # Bernacchi et al., 2001 (assuming 25 C & 98.716 kPa)
    dHa <- 36.38; # kJ/mol
    R <- 0.0083145; # kJ/mol/K
    Ko <- Ko25.pa * exp(dHa * (Tc-25) / (298.15*R*(Tc+273.15)))
    Ko
  }
  # O :: leaf inter-cellular oxygen concentration, Pa
  O.ppm <- 209476; # US Standard Atmosphere
  O.pa <- O.ppm * 10^(-6) * patm; # Pa
  K <- Kc(tc) * (1 + O.pa/Ko(tc))
  K
}

# ************************************************************************
# * Name: calc_n
# *
# * Input: air temperature (tc), degrees C
# *
# * Return: Viscosity of water, mPa s
# *
# * Features: This function calculates the temperature dependent 
# *           viscosity of water by means of a simple exponential 
# *
# * Ref: Vogel equation
# ************************************************************************
calc_n <- function(tc){
  # Temperature needs to be in Kelvin
  0.024263*exp(578.919/((tc+273.15)-137.546))
}

# ************************************************************************
# * Name: temp_sigmoid
# *
# * Input: air temperature (t), degrees C
# *
# * Return: coefficient (ts), (between 0 and 1)
# *
# * Features: This script returns a coefficient (0,1) for a given air 
# *           temperature based on a logistic function (sigmoid) that 
# *           ranges from 0 when temperature is below 0 to 1 when air
# *           temperature is greater than 10
# ************************************************************************
temp_sigmoid <- function(t){
  # Logistic function:
  ts <- 1*(1+1*exp(-0.5*(t-5)))^(-1)
  ts
}

# ************************************************************************
# * Name: lue_conversions
# *
# * Input: data content object
# *
# * Return: data content object
# *
# * Features: This script applies the calculations for unit conversions
# *           (e.g., VPD from kPa to Pa, [CO2] from ppm to Pa), 
# *           calculates Gamma_star, Michaelis-Menten K, and low-
# *           temperature coefficient (temp sigmoid)
# *
# * Depends: - calc_gstar()
# *          - calc_k()
# *          - calc_n()
# *          - temp_sigmoid()
# ************************************************************************
lue_conversions <- function(my.obj){
  # Performs the conversions & calculations on LUE data frame
  # (includes filtering for Tc and fAPAR)
  # Convert VPD to Pa:
  my.obj$VPD.Pa <- (1e3)*my.obj$VPD.kPa
  #my.obj$VPD.Pa[which(my.obj$VPD.Pa < 0)] <- 0
  #
  # Convert CO2 from ppm to Pa:
  my.obj$CO2.Pa <- (1e-6)*my.obj$CO2.ppm*my.obj$Patm.Pa
  #
  # Calculate G* and K:
  my.obj$Gstar.Pa <- calc_gstar(my.obj$Tc.deg_C)
  my.obj$K.Pa <- calc_k(my.obj$Tc.deg_C, my.obj$Patm.Pa)
  #
  # Calculate temperature sigmoid:
  my.obj$kTc <- temp_sigmoid(my.obj$Tc.deg_C)
  #
  # Calculate viscosity of water, mPa s:
  my.obj$eta.mPa_s <- calc_n(my.obj$Tc.deg_C)
  #
  my.obj
}

# ************************************************************************
# * Name: get_meta
# *
# * Input: meta file directory (meta.dir)
# *
# * Return: meta data, data frame (my.params)
# *         station name, char array (st.name)
# *
# * Features: This function retrieves meta data (vegetation & climate) 
# *           for a specific flux tower
# ************************************************************************
get_meta <- function(meta.dir, st.name){
  # Read meta data for specific flux tower
  meta.file <- list.files(path=meta.dir, pattern="Fluxdata_Met-Data*")
  meta.data <- read.csv(
    paste(meta.dir, meta.file, sep=""),
    header=T
  )
  #
  # Find row for flux tower:
  meta.row <- which(meta.data$stationid == st.name)
  #
  # Get climate and vegetation types:
  if (length(meta.row) == 0){
    # Could not find flux tower in meta data?
    meta.clim <- NA()
    meta.veg <- NA()
  } else {
    meta.clim <- as.character(meta.data$climateid[meta.row])
    meta.veg <- as.character(meta.data$classid[meta.row])
  }
  #
  # Send back meta data as data frame:
  meta.frame <- data.frame(
    climate=meta.clim,
    class=meta.veg
  )
  meta.frame
}

# ************************************************************************
# * Name: filter_stations
# *
# * Input: - data frame, LUE contents (my.obj)
# *        - data frame, climate/class IDs (my.ids)
# *
# * Return: data frame (my.filt)
# *
# * Features: This function filters rows from the input data frame
# *           for stations that are not crops or wetlands
# ************************************************************************
filter_stations <- function(my.obj, my.ids){
  # Initialize filter data frame:
  my.filt <- matrix(ncol=(dim(my.obj)[2]),nrow=1)
  my.filt <- as.data.frame(my.filt)
  names(my.filt) <- names(my.obj)
  #
  my.stations <- levels(unique(my.obj$Station))
  #
  # Apply filter:
  num.found <- 0
  for (station in my.stations){
    my.class <- my.ids[station, 'class']
    if(my.class %in% c('CRO', 'WET')){
      num.found <- num.found + 1
    } else {
      my.filt <- rbind(my.filt, my.obj[which(my.obj$Station == station),])
    }
  }
  cat(paste('Filtered', as.character(num.found), 'stations'))
  #
  # Note that row 1 is NA by default (don't return it)
  my.filt[-1,]
}

# ************************************************************************
# * Name: read_sumstats
# *
# * Input: char array <- directory name
# *        char array <- flux station name + timestamp
# *
# * Return: table
# *
# * Features: This script reads the GePiSaT model summary statistics 
# *           file located in the directory (input) and returns a table
# *           of R-squared values corresponding to the various model 
# *           fits made in GePiSaT for a specific flux tower and time
# *           (input)
# ************************************************************************
read_sumstats <- function(my.dir, my.names){
  # Reads summary statistics file in directory
  stats.file = list.files(
    path = my.dir, 
    pattern = "summary_statistics*"
  )
  stats.content <- read.csv(
    paste(my.dir,stats.file,sep=""), 
    header=T
  )
  #
  # Get the R2 values from appropriate model selection.
  # "model_select" field contains numbers 0--4:
  #    -----   ---------
  #    Model :   Field
  #    -----   ---------
  #      0   :  [no fit]
  #      1   :  r2_obs_h
  #      2   :  r2_ro_h
  #      3   :  r2_obs_l
  #      4   :  r2_ro_l
  #
  # Separate indexes for each model type:
  mod.0 <- which(stats.content$model_select == 0)
  mod.1 <- which(stats.content$model_select == 1)
  mod.2 <- which(stats.content$model_select == 2)
  mod.3 <- which(stats.content$model_select == 3)
  mod.4 <- which(stats.content$model_select == 4)
  #
  # Get the right R2 values:
  r2.mod.1 <- stats.content$r2_obs_h[mod.1]
  r2.mod.2 <- stats.content$r2_ro_h[mod.2]
  r2.mod.3 <- stats.content$r2_obs_l[mod.3]
  r2.mod.4 <- stats.content$r2_ro_l[mod.4]
  #
  # Create a searchable list of station name + timestamp:
  my.list <- paste(
    as.character(stats.content$name),
    as.character(stats.content$month),
    sep="_"
  )
  #
  # Create a list for saving the R2 values:
  my.results <- rep(-9999.0, length(my.list))
  my.results[mod.1]<-r2.mod.1
  my.results[mod.2]<-r2.mod.2
  my.results[mod.3]<-r2.mod.3
  my.results[mod.4]<-r2.mod.4
  #
  my.results <- as.table(my.results)
  row.names(my.results) <- my.list
  #
  # Return those that match my.names:
  my.results[my.names]
}

# ************************************************************************
# * Name: plot_vpd
# *
# * Input: data content object (my.obj)
# *        x-axis dataset (my.x)
# *        x-axis title (my.label)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of VPD
# ************************************************************************
plot_vpd <- function(my.obj, my.x, my.label, stid=""){
  # Set plotting extents:
  max.x <- (
    floor(max(my.x)+600) - 
      (floor(max(my.x)+600) %% 200)
  )
  max.y <- (
    floor(max(my.obj$GPP.mol_m2 + my.obj$GPP_err)+5) -
      (floor(max(my.obj$GPP.mol_m2 + my.obj$GPP_err)+5) %% 5)
  )
  #
  # Highlight VPD (coerse with boolean operations)
  vpd.old.cols <- (
    1 + 
      1*(my.obj$VPD.kPa > 0 & my.obj$VPD.kPa < 0.5) + 
      2*(my.obj$VPD.kPa >= 0.5 & my.obj$VPD.kPa < 1.0) +
      3*(my.obj$VPD.kPa >= 1.0 & my.obj$VPD.kPa < 1.5) +
      4*(my.obj$VPD.kPa >= 1.5 & my.obj$VPD.kPa < 2.0) +
      5*(my.obj$VPD.kPa >= 2.0 & my.obj$VPD.kPa < 2.5)
  )
  vpd.cols <- (
    1 + 
      1*(my.obj$VPD.kPa > 0.1 & my.obj$VPD.kPa < 0.2) + 
      2*(my.obj$VPD.kPa >= 0.2 & my.obj$VPD.kPa < 0.3) +
      3*(my.obj$VPD.kPa >= 0.3 & my.obj$VPD.kPa < 0.4) +
      4*(my.obj$VPD.kPa >= 0.4 & my.obj$VPD.kPa < 0.5) +
      5*(my.obj$VPD.kPa >= 0.5)
  )
  #
  # Pull the statistics from the model fit:
  fit <- lm(my.obj$GPP.mol_m2 ~ 0 + my.x)
  fit.coef <- summary(fit)$coefficients[1,1]
  r.sq <- summary(fit)$r.squared
  #
  rp = vector('expression',2)
  rp[1] = substitute(
    expression(epsilon == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r.sq,dig=3)))[2]
  #
  par(mar=c(4.5,4.5,1,1))
  errbar(
    my.x, 
    my.obj$GPP.mol_m2,
    (my.obj$GPP.mol_m2+my.obj$GPP_err),
    (my.obj$GPP.mol_m2-my.obj$GPP_err),
    add=FALSE,
    pch=21,
    errbar.col=vpd.cols,
    col=vpd.cols,
    xlim=c(0,max.x),
    ylim=c(0,max.y),
    xlab=NA, 
    ylab=NA,
    axes=F
  )
  box(lwd=2)
  # Add x-axis tick markers and labels:
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  # Add y-axis tick markers and labels:
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  #
  mtext(side=1, as.expression(my.label), line=2)
  mtext(side=2, expression(GPP~(mol~CO[2]%.%m^{-2})), line=2)
  mtext(side=1, as.expression(stid), line=3)
  legend('bottomright', legend = rp, bty = 'n')
  legend('topright',
         legend=c("0.0-0.1 kPa",
                  "0.1-0.2 kPa",
                  "0.2-0.3 kPa",
                  "0.3-0.4 kPa",
                  "0.4-0.5 kPa",
                  "    >0.5 kPa"
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_xy_vpd
# *
# * Input: data content object (my.obj)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots x verus y, and color-codes the data 
# *           points based on binned ranges of VPD
# ************************************************************************
plot_xy_vpd <- function(my.obj, stid=""){
  # Highlight VPD (coerse with boolean operations)
  vpd.old.cols <- (
    1 + 
      1*(my.obj$VPD.kPa > 0 & my.obj$VPD.kPa < 0.5) + 
      2*(my.obj$VPD.kPa >= 0.5 & my.obj$VPD.kPa < 1.0) +
      3*(my.obj$VPD.kPa >= 1.0 & my.obj$VPD.kPa < 1.5) +
      4*(my.obj$VPD.kPa >= 1.5 & my.obj$VPD.kPa < 2.0) +
      5*(my.obj$VPD.kPa >= 2.0 & my.obj$VPD.kPa < 2.5)
  )
  vpd.cols <- (
    1 + 
      1*(my.obj$VPD.kPa > 0.1 & my.obj$VPD.kPa < 0.2) + 
      2*(my.obj$VPD.kPa >= 0.2 & my.obj$VPD.kPa < 0.3) +
      3*(my.obj$VPD.kPa >= 0.3 & my.obj$VPD.kPa < 0.4) +
      4*(my.obj$VPD.kPa >= 0.4 & my.obj$VPD.kPa < 0.5) +
      5*(my.obj$VPD.kPa >= 0.5)
  )
  #
  # Get monthly LUE values:
  phi <- my.obj$GPP.mol_m2/(my.obj$fAPAR*my.obj$PPFD.mol_m2)
  #
  # Save x and y values:
  my.y <- (my.obj$CO2.Pa-my.obj$Gstar.Pa)/(phi*(my.obj$CO2.Pa+2*my.obj$Gstar.Pa))
  my.x <- 3*my.obj$Gstar.Pa*sqrt(
    1.6*my.obj$eta.mPa_s*my.obj$VPD.Pa/(my.obj$K.Pa+my.obj$Gstar.Pa)
  )/(my.obj$CO2.Pa+2*my.obj$Gstar.Pa)
  #
  # Make regression:
  fit <- lm(my.y ~ my.x)
  #
  # Get regression coefficients (y = A + Bx)
  A <- summary(fit)$coefficients['(Intercept)','Estimate']
  B <- summary(fit)$coefficients['my.x','Estimate']
  A.err <- summary(fit)$coefficients['(Intercept)','Std. Error']
  B.err <- summary(fit)$coefficients['my.x','Std. Error']
  #
  # Get regression fitness:
  r2 <- summary(fit)$adj.r.squared
  #
  # Calculate phi_o and beta:
  phi_o <- 1.0/A
  phi_o.err <- A.err/A^2
  beta <- (A/B)^2
  beta.err <- sqrt((2*A/B^2)^2*A.err^2 + (-2*A^2/B^3)^2*B.err^2)
  #
  # Plot expression:
  rp = vector('expression',3)
  rp[1] = substitute(
    expression(phi[o] == MYLUE %+-% LUEERR),
    list(MYLUE = format(phi_o,dig=3), LUEERR = format(phi_o.err, dig=3)))[2]
  rp[2] = substitute(
    expression(beta == MYBETA %+-% BETAERR),
    list(MYBETA = format(beta,dig=3), BETAERR = format(beta.err, dig=3)))[2]
  rp[3] = substitute(
    expression(italic(R)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  #
  par(mar=c(4.5,4.5,1,1))
  plot(my.x, my.y, xlab=NA, ylab=NA, axes=F, pch=1, col=vpd.cols, 
       xlim=c(0,max(my.x)))
  box(lwd=2)
  # Add x-axis tick markers and labels:
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  # Add y-axis tick markers and labels:
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  #
  mtext(side=1, as.expression(x~(mPa%.%s)^0.5), line=2)
  mtext(side=2, as.expression(y~(unitless)), line=3)
  mtext(side=1, as.expression(stid), line=3)
  legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c("0.0-0.1 kPa",
                  "0.1-0.2 kPa",
                  "0.2-0.3 kPa",
                  "0.3-0.4 kPa",
                  "0.4-0.5 kPa",
                  "    >0.5 kPa"
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_alpha
# *
# * Input: data content object (my.obj)
# *        x-axis dataset (my.x)
# *        x-axis title, expression (my.label)
# *        meta data directory, char array (mdir)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of Cramer-Prentice alpha
# *
# * Depends: get_meta
# ************************************************************************
plot_alpha <- function(my.obj, my.x, my.label, mdir, stid=""){
  # Set plotting extents:
  max.x <- (
    floor(max(my.x)+600) - 
      (floor(max(my.x)+600) %% 200)
  )
  max.y <- (
    floor(max(my.obj$GPP.mol_m2 + my.obj$GPP_err)+5) -
      (floor(max(my.obj$GPP.mol_m2 + my.obj$GPP_err)+5) %% 5)
  )
  #
  # Get climate & vegetation meta data:
  my.meta <- get_meta(mdir, stid)
  #
  # Highlight alpha (coerse with boolean operations)
  alpha.cols <- (
    1 + 
      1*(my.obj$ALPHA >= 0.3 & my.obj$ALPHA < 0.5) + 
      2*(my.obj$ALPHA >= 0.5 & my.obj$ALPHA < 0.7) +
      3*(my.obj$ALPHA >= 0.7 & my.obj$ALPHA < 0.9) +
      4*(my.obj$ALPHA >= 0.9 & my.obj$ALPHA < 1.1) +
      5*(my.obj$ALPHA >= 1.1 & my.obj$ALPHA < 1.259) +
      6*(my.obj$ALPHA >= 1.259)
  )
  #
  # Pull the statistics from the model fit:
  fit <- lm(my.obj$GPP.mol_m2 ~ 0 + my.x, weights=(1/my.obj$GPP_err))
  fit.coef <- summary(fit)$coefficients['my.x', 'Estimate']
  r.sq <- summary(fit)$adj.r.squared
  #
  rp = vector('expression',2)
  rp[1] = substitute(
    expression(epsilon == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r.sq,dig=3)))[2]
  #
  par(mar=c(4.5,4.5,1,1))
  errbar(
    my.x, 
    my.obj$GPP.mol_m2,
    (my.obj$GPP.mol_m2+my.obj$GPP_err),
    (my.obj$GPP.mol_m2-my.obj$GPP_err),
    add=FALSE,
    pch=19,
    xlab=NA, 
    ylab=NA,
    axes=F,
    errbar.col=alpha.cols,
    col=alpha.cols,
    xlim=c(0,max.x),
    ylim=c(0,max.y)
  )
  abline(fit, lty=2, lwd=2, col='gray70')
  box(lwd=2)
  # Add x-axis tick markers and labels:
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  # Add y-axis tick markers and labels:
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  #
  mtext(side=1, as.expression(my.label), line=2)
  mtext(side=2, expression(GPP~(mol~CO[2]%.%m^{-2})), line=2)
  mtext(side=1, as.expression(stid), line=3)
  #
  mtext(paste("Climate ID: ", as.character(my.meta$climate), "\n",
              "Class ID: ", as.character(my.meta$class)),
        side=1, line=3, adj=1, cex=0.8)
  #
  legend('bottomright', legend = rp, bty = 'n')
  #legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c(expression(alpha < 0.3),
                  expression(alpha:~0.3-0.5),
                  expression(alpha:~0.5-0.7),
                  expression(alpha:~0.7-0.9),
                  expression(alpha:~0.9-1.1),
                  expression(alpha:~1.1-1.25),
                  expression(alpha == 1.26)
         ),
         pch=c(19,19,19,19,19,19,19), 
         col=c(1,2,3,4,5,6,7),
         box.lwd=0,
         box.col='white',
         bg='white',
         inset=0.01
  )
}

# ************************************************************************
# * Name: plot_temp
# *
# * Input: data content object (my.obj)
# *        x-axis dataset (my.x)
# *        x-axis title (my.label)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of air temperature
# ************************************************************************
plot_temp <- function(my.obj, my.x, my.label, stid=""){
  # Set maxPPFD:
  max.PPFD <- max(my.obj$PPFD.mol_m2)
  max.GPP <- max(my.obj$GPP.mol_m2 + my.obj$GPP_err)
  #
  # Highlight temperature (coerse with boolean operations)
  temp.cols <- (
    1 + 
      1*(my.obj$Tc.deg_C > 5 & my.obj$Tc.deg_C < 10) + 
      2*(my.obj$Tc.deg_C >= 10 & my.obj$Tc.deg_C < 15) +
      3*(my.obj$Tc.deg_C >= 15 & my.obj$Tc.deg_C < 20) +
      4*(my.obj$Tc.deg_C >= 20 & my.obj$Tc.deg_C < 25) +
      5*(my.obj$Tc.deg_C >= 25)
  )
  #
  # Pull the statistics from the model fit:
  fit <- lm(my.obj$GPP.mol_m2 ~ 0 + my.x)
  fit.coef <- summary(fit)$coefficients[1,1]
  r.sq <- summary(fit)$r.squared
  #
  rp = vector('expression',2)
  rp[1] = substitute(
    expression(epsilon == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r.sq,dig=3)))[2]
  #
  par(mar=c(5.5,5.5,1,1))
  errbar(
    my.x, 
    my.obj$GPP.mol_m2,
    (my.obj$GPP.mol_m2+my.obj$GPP_err),
    (my.obj$GPP.mol_m2-my.obj$GPP_err),
    add=FALSE,
    pch=21,
    xlab=my.label, 
    ylab=expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")),
    errbar.col=temp.cols,
    col=temp.cols,
    xlim=c(0,max.PPFD),
    ylim=c(0,max.GPP)
  )
  box(lwd=2)
  title(sub = stid)
  #abline(a=0,b=fit.coef,lty=2)
  legend('bottomright', legend = rp, bty = 'n')
  #legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c(expression(0-5~degree*C),
                  expression(5-10~degree*C),
                  expression(10-15~degree*C),
                  expression(15-20~degree*C),
                  expression(20-25~degree*C),
                  expression(25-30~degree*C)
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_xy_temp
# *
# * Input: data content object (my.obj)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of air temperature
# ************************************************************************
plot_xy_temp <- function(my.obj, stid=""){
  # Highlight temperature (coerse with boolean operations)
  temp.cols <- (
    1 + 
      1*(my.obj$Tc.deg_C >= 5 & my.obj$Tc.deg_C < 10) + 
      2*(my.obj$Tc.deg_C >= 10 & my.obj$Tc.deg_C < 15) +
      3*(my.obj$Tc.deg_C >= 15 & my.obj$Tc.deg_C < 20) +
      4*(my.obj$Tc.deg_C >= 20 & my.obj$Tc.deg_C < 25) +
      5*(my.obj$Tc.deg_C >= 25)
  )
  #
  # Get monthly LUE values:
  phi <- my.obj$GPP.mol_m2/(my.obj$fAPAR*my.obj$PPFD.mol_m2)
  #
  # Save x and y values:
  my.y <- (my.obj$CO2.Pa-my.obj$Gstar.Pa)/(phi*(my.obj$CO2.Pa+2*my.obj$Gstar.Pa))
  my.x <- 3*my.obj$Gstar.Pa*sqrt(
    1.6*my.obj$eta.mPa_s*my.obj$VPD.Pa/(my.obj$K.Pa+my.obj$Gstar.Pa)
  )/(my.obj$CO2.Pa+2*my.obj$Gstar.Pa)
  #
  # Make regression:
  fit <- lm(my.y ~ my.x)
  #
  # Get regression coefficients (y = A + Bx)
  A <- summary(fit)$coefficients['(Intercept)','Estimate']
  B <- summary(fit)$coefficients['my.x','Estimate']
  A.err <- summary(fit)$coefficients['(Intercept)','Std. Error']
  B.err <- summary(fit)$coefficients['my.x','Std. Error']
  #
  # Get regression fitness:
  r2 <- summary(fit)$adj.r.squared
  #
  # Calculate phi_o and beta:
  phi_o <- 1.0/A
  phi_o.err <- A.err/A^2
  beta <- (A/B)^2
  beta.err <- sqrt((2*A/B^2)^2*A.err^2 + (-2*A^2/B^3)^2*B.err^2)
  #
  # Plot expression:
  rp = vector('expression',3)
  rp[1] = substitute(
    expression(phi[o] == MYLUE %+-% LUEERR),
    list(MYLUE = format(phi_o,dig=3), LUEERR = format(phi_o.err, dig=3)))[2]
  rp[2] = substitute(
    expression(beta == MYBETA %+-% BETAERR),
    list(MYBETA = format(beta,dig=3), BETAERR = format(beta.err, dig=3)))[2]
  rp[3] = substitute(
    expression(italic(R)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  #
  par(mar=c(4.5,4.5,1,1))
  plot(my.x, my.y, xlab=NA, ylab=NA, axes=F, pch=1, col=temp.cols, 
       xlim=c(0,max(my.x)))
  box(lwd=2)
  # Add x-axis tick markers and labels:
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  # Add y-axis tick markers and labels:
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  #
  mtext(side=1, as.expression(x~(mPa%.%s)^0.5), line=2)
  mtext(side=2, as.expression(y~(unitless)), line=3)
  mtext(side=1, as.expression(stid), line=3)
  legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c(expression(0-5~degree*C),
                  expression(5-10~degree*C),
                  expression(10-15~degree*C),
                  expression(15-20~degree*C),
                  expression(20-25~degree*C),
                  expression(25-30~degree*C)
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_temp_alpha
# *
# * Input: data content object (my.obj)
# *        x-axis dataset (my.x)
# *        x-axis title (my.label)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of air temperature and 
# *           scales the data point size (circle diameter) based on 
# *           binned ranges of Cramer-Prentice alpha
# ************************************************************************
plot_temp_alpha <- function(my.obj, my.x, my.label, stid=""){
  # Set maxPPFD:
  max.PPFD <- max(my.obj$PPFD.mol_m2)
  max.GPP <- max(my.obj$GPP.mol_m2)
  #
  # Highlight temperature (coerse with boolean operations)
  temp.cols <- (
    1 + 
      1*(my.obj$Tc.deg_C > 5 & my.obj$Tc.deg_C < 10) + 
      2*(my.obj$Tc.deg_C >= 10 & my.obj$Tc.deg_C < 15) +
      3*(my.obj$Tc.deg_C >= 15 & my.obj$Tc.deg_C < 20) +
      4*(my.obj$Tc.deg_C >= 20 & my.obj$Tc.deg_C < 25) +
      5*(my.obj$Tc.deg_C >= 25)
  )
  alpha.symbol <- (
    1 + 
      1*(my.obj$ALPHA > 0.3 & my.obj$ALPHA <= 0.4) + 
      2*(my.obj$ALPHA > 0.2 & my.obj$ALPHA <= 0.3) +
      3*(my.obj$ALPHA > 0.1 & my.obj$ALPHA <= 0.2) +
      4*(my.obj$ALPHA >= 0.0 & my.obj$ALPHA <= 0.1)
  )
  #
  # Pull the statistics from the model fit:
  fit <- lm(my.obj$GPP.mol_m2 ~ 0 + my.x)
  fit.coef <- summary(fit)$coefficients[1,1]
  r.sq <- summary(fit)$r.squared
  #
  rp = vector('expression',2)
  rp[1] = substitute(
    expression(epsilon == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r.sq,dig=3)))[2]
  #
  par(mar=c(5.5,5.5,1,1))
  symbols(
    my.x, 
    my.obj$GPP.mol_m2,
    xlab=my.label, 
    ylab=expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")),
    sub = stid,
    circles=alpha.symbol,
    inches=1/3,
    fg=temp.cols,
    xlim=c(0,max.PPFD),
    ylim=c(0,max.GPP)
  )
  box(lwd=2)
  #abline(a=0,b=fit.coef,lty=2)
  legend('bottomright', legend = rp, bty = 'n')
  #legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c(expression(0-5~degree*C),
                  expression(5-10~degree*C),
                  expression(10-15~degree*C),
                  expression(15-20~degree*C),
                  expression(20-25~degree*C),
                  expression(25-30~degree*C)
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: geom.mean
# *
# * Input: data array (my.vals)
# *
# * Return: value (my.product)
# *
# * Features: This script calculates the geometric mean of a list of 
# *           values.
# ************************************************************************
geom.mean <- function(my.vals){
  # Computes the geometric mean of values
  n <- length(my.vals)
  my.product = 1.0
  for (i in my.vals){
    my.product = my.product * i
  }
  my.product = my.product^(1.0/n)
  my.product
}

# ************************************************************************
# * Name: estimate.beta
# *
# * Input: - vapor pressure deficit, Pa (vpd)
# *        - Michaelis-Menten K, Pa (k)
# *        - photorespiratory compensation point (gs)
# *        - viscosity of water, mPa s (n)
# *
# * Return: value (beta)
# *
# * Features: This script estimates the magnitude of the parameter 
# *           beta (unknown parameter in next-gen LUE model) based on
# *           a reformulation of the Chi equation and 'squiggle' 
# *           equation, based on values of n, VPD, K and Γ*.
# *
# * UPDATED EQUATION: b/a  =  1.6 [D/(K + Γ*)] [(ci − Γ*)/(ca − ci)]^2
# *                   b/a  =  1.6 [D/(K + Γ*)] [(χ − Γ*)/(1 − χ)]^2
# *                   beta = n*(b/a)
# ************************************************************************
estimate.beta <- function(vpd, k, gs, n){
  # Vary chi depending on VPD (assumed linear)
  #   0.5 :: very dry (VPD = 2500 Pa); 
  #   0.9 :: very wet (VPD = 0 Pa)
  #
  chi <- 0.5 + (0.5-0.9)*(2500-vpd)/(-2500)
  ba <- 1.6*(vpd/(k+gs))*((chi - gs)/(1.0-chi))^2
  beta <- n*(ba)
  beta
}

# ************************************************************************
# * Name: calc_m
# *
# * Input: data content object (my.obj)
# *
# * Return: list of values
# *
# * Features: This script calculates the coefficient m in the next-gen
# *           LUE equation based on an estimate of the unknown parameter
# *           beta, n*(b/a).
# *
# *           m = (pa-Gs)/(pa+2Gs+3Gs*(1.6n*VPD/(beta*(K+Gs))^0.5))
# *               denom1 = beta*(K+Gs)
# *               denom2 = pa+2Gs+3Gs*(1.6n*VPD/denom1)^0.5
# *               m = (pa-Gs)/denom2
# * 
# ************************************************************************
calc_m <- function(my.obj, beta){
  # Calculate m:
  numer <- (my.obj$CO2.Pa - my.obj$Gstar.Pa)
  denom1 <- beta*(my.obj$K.Pa + my.obj$Gstar.Pa)
  if(any(denom1 == 0)){
    denom1 = denom1 + 1.0e-6
  }
  denom2 <- (
    my.obj$CO2.Pa + 2.0*my.obj$Gstar.Pa + 
      3.0*my.obj$Gstar.Pa*(1.6*my.obj$eta.mPa_s*my.obj$VPD.Pa/denom1)^0.5
  )
  if(any(denom2 == 0)){
    denom2 = denom2 + 1.0e-6
  }
  # Return:
  (numer/denom2)
}

# ************************************************************************
# * Name: calc_mp
# *
# * Input: data content object (my.obj)
# *
# * Return: list of values
# *
# * Features: This script calculates the coefficient m' in the next-gen
# *           LUE equation.
# *
# *           m' = m*[1 - (0.41/m)^(2/3)]
# ************************************************************************
calc_mp <- function(my.obj, m){
  # Calculate m':
  # NOTE: 2/3rds power cannot be applied to negative numbers
  m.temp <- m
  m.temp[which(m.temp<0)]<-NA
  m.temp <- 1-(0.41/m.temp)^(2/3)
  m.temp[which(m.temp<0)]<-NA
  #mprime <- my.obj$m*sqrt(m.temp)
  mprime <- m*sqrt(m.temp)
  mprime
}

# ************************************************************************
# * Name: plot_mbD
# *
# * Input: data content object (my.obj)
# *
# * Return: None.
# *
# * Features: This script plots b versus m, highlighting levels of VPD.
# ************************************************************************
plot_mbD <- function(my.obj){
  # Highlight VPD (coerse with boolean operations)
  vpd.cols <- (
    1 + 
      1*(my.obj$D > 0.2 & my.obj$D < 0.4) + 
      2*(my.obj$D >= 0.4 & my.obj$D < 0.6) +
      3*(my.obj$D >= 0.6 & my.obj$D < 0.8) +
      4*(my.obj$D >= 0.8 & my.obj$D < 1.0) +
      5*(my.obj$D >= 1.0)
  )
  #
  par(mar=c(4.75,4.75,1,1))
  plot(
    my.obj$b, 
    my.obj$m,
    xlab=expression(sqrt(b/a)), 
    ylab="m",
    col=vpd.cols,
  )
  box(lwd=2)
  legend('bottomright',
         legend=c("0.0-0.2 kPa",
                  "0.2-0.4 kPa",
                  "0.4-0.6 kPa",
                  "0.6-0.8 kPa",
                  "0.8-1.0 kPa",
                  ">1.0 kPa"
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(1,2,3,4,5,6),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_r2
# *
# * Input: data content object (my.obj)
# *        x-axis dataset (my.x)
# *        x-axis title (my.label)
# *        station name, optional for subtitle (stid)
# *
# * Return: None.
# *
# * Features: This script plots GPP versus x-axis data (my.x), applies
# *           error bars for GPP (y-axis), and color-codes the data 
# *           points based on binned ranges of R-squared (fitness)
# ************************************************************************
plot_r2<-function(my.obj, my.x, my.label, my.dir, stid = ""){
  # Plots basic LUE with R2 highlighted:
  #
  # Set maxPPFD:
  max.PPFD <- max(my.obj$PPFD.mol_m2)
  max.GPP <- max(my.obj$GPP.mol_m2 + my.obj$GPP_err)
  #
  # Get R2 values:
  st.titles <- paste(
    as.character(my.obj$Station), 
    as.character(my.obj$Timestamp), 
    sep="_"
  )
  my.obj$r2 <- read_sumstats(my.dir, st.titles)
  #
  # Highlight R2 (coerse with boolean operations)
  r2.cols <- (
    1 + 
      1*(my.obj$r2 >= 0.2 & my.obj$r2 < 0.3) + 
      2*(my.obj$r2 >= 0.3 & my.obj$r2 < 0.4) +
      3*(my.obj$r2 >= 0.4 & my.obj$r2 < 0.5) +
      4*(my.obj$r2 >= 0.5 & my.obj$r2 < 0.6) +
      5*(my.obj$r2 >= 0.6 & my.obj$r2 < 0.8) +
      7*(my.obj$r2 >= 0.8)
  )
  #
  # Pull the statistics from the model fit:
  fit <- lm(my.obj$GPP.mol_m2 ~ 0 + my.x)
  fit.coef <- summary(fit)$coefficients[1,1]
  r.sq <- summary(fit)$r.squared
  #
  rp = vector('expression',2)
  rp[1] = substitute(
    expression(epsilon == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r.sq,dig=3)))[2]
  #
  par(mar=c(5.5,5.5,1,1))
  errbar(
    my.x, 
    my.obj$GPP.mol_m2,
    (my.obj$GPP.mol_m2+my.obj$GPP_err),
    (my.obj$GPP.mol_m2-my.obj$GPP_err),
    add=FALSE,
    pch=21,
    xlab=my.label, 
    ylab=expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")),
    errbar.col=r2.cols,
    col=r2.cols,
    xlim=c(0,max.PPFD),
    ylim=c(0,max.GPP)
  )
  box(lwd=2)
  title(sub = stid)
  #abline(a=0,b=fit.coef,lty=2)
  legend('bottomright', legend = rp, bty = 'n')
  #legend('topleft', legend = rp, bty = 'n')
  legend('topright',
         legend=c(expression(0.2-0.3),
                  expression(0.3-0.4),
                  expression(0.4-0.5),
                  expression(0.5-0.6),
                  expression(0.6-0.8),
                  expression(0.8-1.0)
         ),
         pch=c(1,1,1,1,1,1), 
         col=c(2,3,4,5,6,8),
         bty = 'n'
  )
}

# ************************************************************************
# * Name: plot_xy
# *
# * Input: data content object (my.obj)
# *        boolean for using temperature sigmoid (tau)
# *        boolean for using Cramer-Prentice alpha (alpha)
# *        numeric for alpha function (k)
# *        char for legend position (pos)
# *        char for station id name
# *
# * Return: None.
# *
# * Features: This script plots y = A + Bx where y and x are defined
# *           based on Colin's method such that phi_o and beta can be 
# *           extracted from the fitted regression parameters
# *
# * Note: In Colin's method, to add tau (the temperature sigmoid) or 
# *       alpha, simply multiply it to y, i.e, y = y*tau
# ************************************************************************
plot_xy <- function(my.obj, tau, alpha, k=0, pos, stid=""){
  # Get monthly LUE values:
  phi <- my.obj$GPP.mol_m2/(my.obj$fAPAR*my.obj$PPFD.mol_m2)
  #
  # Save x and y values:
  my.y <- (my.obj$CO2.Pa-my.obj$Gstar.Pa)/(phi*(my.obj$CO2.Pa+2*my.obj$Gstar.Pa))
  if (tau){
    my.y <- (my.obj$kTc*my.y)
  }
  if (alpha){
    my.y <- (my.obj$ALPHA)^k*my.y
  }
  my.x <- 3*my.obj$Gstar.Pa*sqrt(
    1.6*my.obj$eta.mPa_s*my.obj$VPD.Pa/(my.obj$K.Pa+my.obj$Gstar.Pa)
  )/(my.obj$CO2.Pa+2*my.obj$Gstar.Pa)
  #
  # Make regression:
  fit <- lm(my.y ~ my.x)
  #
  # Get regression coefficients (y = A + Bx)
  A <- summary(fit)$coefficients['(Intercept)','Estimate']
  B <- summary(fit)$coefficients['my.x','Estimate']
  A.err <- summary(fit)$coefficients['(Intercept)','Std. Error']
  B.err <- summary(fit)$coefficients['my.x','Std. Error']
  #
  # Get regression fitness:
  r2 <- summary(fit)$adj.r.squared
  #
  # Calculate phi_o and beta:
  phi_o <- 1.0/A
  phi_o.err <- A.err/A^2
  beta <- (A/B)^2
  beta.err <- sqrt((2*A/B^2)^2*A.err^2 + (-2*A^2/B^3)^2*B.err^2)
  #
  # Plot expression:
  rp = vector('expression',3)
  rp[1] = substitute(
    expression(phi[o] == MYLUE %+-% LUEERR),
    list(MYLUE = format(phi_o,dig=3), LUEERR = format(phi_o.err, dig=3)))[2]
  rp[2] = substitute(
    expression(beta == MYBETA %+-% BETAERR),
    list(MYBETA = format(beta,dig=3), BETAERR = format(beta.err, dig=3)))[2]
  rp[3] = substitute(
    expression(italic(R)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  #
  rq = vector('expression',3)
  rq[1] = substitute(
    expression(tau == MYTAU),
    list(MYTAU = as.character(tau)))[2]
  rq[2] = substitute(
    expression(alpha == MYALPHA),
    list(MYALPHA = as.character(alpha)))[2]
  rq[3] = substitute(
    expression(k == MYK),
    list(MYK = format(k,dig=2)))[2]
  # Plot
  par(mar=c(5.5,4.5,1,1))
  plot(my.x, my.y, xlab=NA, ylab=NA, axes=F, pch=19, col='darkgray', xlim=c(0,max(my.x)))
  abline(fit, col='red', lwd=2, lty=2)
  box(lwd=2)
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  mtext(side=1, as.expression(x~(mPa%.%s)^0.5), line=2)
  mtext(side=1, stid, line=4)
  mtext(side=2, as.expression(y~(unitless)), line=3)
  mtext(rq[1], side=1, line=2, adj=0, cex=0.8)
  mtext(rq[2], side=1, line=3, adj=0, cex=0.8)
  mtext(rq[3], side=1, line=4, adj=0, cex=0.8)
  legend(pos, legend=rp, bty='n')
}

# ************************************************************************
# * Name: plot_lue_xy
# *
# * Input: lue analyzer object (my.lue)
# *        char for legend position (pos)
# *        char for station id name
# *
# * Return: None.
# *
# * Features: This script plots y = A + Bx where y and x are defined in
# *           such a way that phi_o and beta can be extracted from the
# *           fitted regression parameters
# *
# ************************************************************************
plot_lue_xy <- function(my.lue, pos, stid=""){
  # Calculate missing variables for LUE:
  my.lue$gstar <- calc_gstar(my.lue$mgs.tc)
  my.lue$eta <- calc_n(my.lue$mgs.tc)
  my.lue$k <- calc_k(my.lue$mgs.tc, my.lue$patm)
  #
  # Save x and y values:
  # (Colin's method)
  my.y <- (my.lue$co2-my.lue$gstar)/(my.lue$lue*(my.lue$co2+2*my.lue$gstar))
  my.x <- 3*my.lue$gstar*sqrt(
    1.6*my.lue$eta*my.lue$mgs.vpd/(my.lue$k+my.lue$gstar)
  )/(my.lue$co2+2*my.lue$gstar)
  #
  # NOTE: there is uncertainty in the y variable due to the uncertainty
  #       in the LUE parameter (i.e., lue_err), which is, as of now, 
  #       accounted for in the regression parameters by means of a
  #       weighted regression
  #
  # Define weights for linear model:
  #  assume NaNs have no error
  my.lue$lue_err[which(is.nan(my.lue$lue_err))]<-NA  # set nans to na
  my.weights <- 1/my.lue$lue_err                     # take reciprocal for weights
  my.weights <- my.weights/max(my.weights, na.rm=T)  # divide by max (range 0--1)
  my.weights[which(is.na(my.weights))]<- 1           # set weight=1 for error-free
  #
  # Make regression:
  fit <- lm(my.y ~ my.x, weights = my.weights)
  #
  # Get regression coefficients (y = A + Bx)
  A <- summary(fit)$coefficients['(Intercept)','Estimate']
  B <- summary(fit)$coefficients['my.x','Estimate']
  A.err <- summary(fit)$coefficients['(Intercept)','Std. Error']
  B.err <- summary(fit)$coefficients['my.x','Std. Error']
  #
  # Get regression fitness:
  r2 <- summary(fit)$adj.r.squared
  #
  # Calculate phi_o and beta:
  # (Colin's method)
  phi_o <- 1.0/A
  phi_o.err <- A.err/A^2
  beta <- (A/B)^2
  beta.err <- sqrt((2*A/B^2)^2*A.err^2 + (-2*A^2/B^3)^2*B.err^2)
  #
  # Plot expression:
  rp = vector('expression',3)
  rp[1] = substitute(
    expression(phi[o] == MYLUE %+-% LUEERR),
    list(MYLUE = format(phi_o,dig=3), LUEERR = format(phi_o.err, dig=3)))[2]
  rp[2] = substitute(
    expression(beta == MYBETA %+-% BETAERR),
    list(MYBETA = format(beta,dig=3), BETAERR = format(beta.err, dig=3)))[2]
  rp[3] = substitute(
    expression(italic(R)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  #
  # Plot
  par(mar=c(5.5,4.5,1,1))
  plot(my.x, my.y, xlab=NA, ylab=NA, axes=F, pch=19, col='darkgray',xlim=c(0,max(my.x)))
  abline(fit, col='red', lwd=2, lty=2)
  box(lwd=2)
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  mtext(side=1, as.expression(x~(mPa%.%s)^0.5), line=2)
  mtext(side=1, stid, line=4)
  mtext(side=2, as.expression(y~(unitless)), line=3)
  legend(pos, legend=rp, bty='n')
}

# ************************************************************************
# * Name: lue_analyzer
# *
# * Input: data frame with lue analyzer results for each station (my.frame)
# *        data content object (my.obj)
# *        location where summary stats file is located (my.dir)
# *        flux tower station name (st.name)
# *
# * Return: data frame
# *
# * Features: This script finds the climate and elevation data for a
# *           flux tower and calculates the mean growing season air
# *           temperature and VPD, basic LUE (with associated error)
# *           atmospheric pressure and CO2
# ************************************************************************
lue.analyzer <- function(my.frame, my.obj, my.dir, st.name, mac){
  # ~~~~~~~~~~~~ GET CLIMATE ~~~~~~~~~~~~~
  # Open meta data file and find climate variable:
  met.file = list.files(
    path = my.dir, 
    pattern = "Fluxdata_Met-Data_2002-06.csv"
  )
  met.content <- read.csv(
    paste(my.dir,met.file,sep=""), 
    header=T
  )
  #
  # Climate variable: "climateid"
  my.climate <- as.character(
    met.content[which(met.content$stationid == st.name),'climateid']
  )
  #
  # ~~~~~ GET MEAN GROWING SEASON & ELEVATION ~~~~~
  if (mac){
    mgs.path <- paste(
      "/Users/twdavis/Dropbox/", # Mac
      "Work/Imperial/flux/results/2002-06/time-series_analysis/",
      sep=""
    )
  } else {
    mgs.path <- paste(
      "/home/user/Dropbox/",     # Linux
      "Work/Imperial/flux/results/2002-06/time-series_analysis/",
      sep=""
    )
  }
  mgs.file <- list.files(
    path = mgs.path,
    pattern = "*.csv"
  )
  mgs.content <- read.csv(
    paste(mgs.path, mgs.file, sep=""),
    header=T
  )
  # Get elevation (meters)
  mgs.ele <- mgs.content[which(mgs.content$site == st.name), 'ele.m']
  #
  # Get growing season mean temperatures (deg. C)
  mgs.tc <- 0*seq(from=1, to=5, by=1)
  mgs.tc[1] <- mgs.content[which(mgs.content$site == st.name), 'tc.2002']
  mgs.tc[2] <- mgs.content[which(mgs.content$site == st.name), 'tc.2003']
  mgs.tc[3] <- mgs.content[which(mgs.content$site == st.name), 'tc.2004']
  mgs.tc[4] <- mgs.content[which(mgs.content$site == st.name), 'tc.2005']
  mgs.tc[5] <- mgs.content[which(mgs.content$site == st.name), 'tc.2006']
  #
  # Get growing season mean VPD (kPa)
  mgs.vpd <- 0*seq(from=1, to=5, by=1)
  mgs.vpd[1] <- mgs.content[which(mgs.content$site == st.name), 'vpd.2002']
  mgs.vpd[2] <- mgs.content[which(mgs.content$site == st.name), 'vpd.2003']
  mgs.vpd[3] <- mgs.content[which(mgs.content$site == st.name), 'vpd.2004']
  mgs.vpd[4] <- mgs.content[which(mgs.content$site == st.name), 'vpd.2005']
  mgs.vpd[5] <- mgs.content[which(mgs.content$site == st.name), 'vpd.2006']
  #
  # ~~~~~~~~~~~~ GET LUE ~~~~~~~~~~~~~
  # Fit the basic LUE regression:
  my.x <- my.obj$fAPAR * my.obj$PPFD.mol_m2
  fit <- lm(
    my.obj$GPP.mol_m2 ~ 0 + my.x
  )
  my.lue <- summary(fit)$coefficients[1,'Estimate']
  my.err <- summary(fit)$coefficients[1,'Std. Error']
  #
  # ~~~~~~~~~~~~ GET Patm & CO2 ~~~~~~~~~~~~~
  my.patm <- mean(my.obj$Patm.Pa)
  my.co2 <- mean(my.obj$CO2.Pa)
  #
  # ~~~~~~~~~~~~  CALC MAX X & Y   ~~~~~~~~~~~~~
  max.x <- (
    floor(max(my.x)+100) - 
      (floor(max(my.x)+100) %% 100)
  )
  if (max.x == 900){
    max.x <- 1000
  }
  max.y <- (
    floor(max(my.obj$GPP.mol_m2)+10) - 
      (floor(max(my.obj$GPP.mol_m2)+10) %% 10)
  )
  # ~~~~~~~~~~~~ UPDATE DATA FRAME ~~~~~~~~~~~~~
  my.frame[st.name, 'lue'] <- my.lue
  my.frame[st.name, 'lue_err'] <- my.err
  my.frame[st.name, 'maxX'] <- max.x
  my.frame[st.name, 'maxY'] <- max.y
  my.frame[st.name, 'clim'] <- my.climate
  my.frame[st.name, 'elv'] <- mgs.ele
  my.frame[st.name, 'mgs.tc'] <- mean(mgs.tc)
  my.frame[st.name, 'mgs.vpd'] <- mean(mgs.vpd)
  my.frame[st.name, 'patm'] <- my.patm
  my.frame[st.name, 'co2'] <- my.co2
  #
  # Return data frame:
  my.frame
}

# ************************************************************************
# * Name: plot.lue.analyzer
# *
# * Input: data frame, based on lue.analyzer function (lue.frame)
# *        plotting flag (flag)
# *
# * Return: None.
# *
# * Features: This script plots LUE by climate types (flag=1) or 
# *           precipitation zones (flag=2)
# ************************************************************************
plot.lue.analyzer <- function(lue.frame, flag){
  # Plots LUE by climate type
  # Definition of input variables:
  #   lue.frame :: data frame from lue.analyzer()
  #   flag      :: tells plotter which type of plot to plot
  #                1 - Koppen-Geiger main climate regions
  #
  # Get max x and y values:
  max.x <- max(lue.frame$maxX)
  max.y <- max(lue.frame$maxY)
  #
  if(flag == 1){
    # Plot main Koppen-Geiger climate types:
    #   A :: equitorial
    #   B :: arid
    #   C :: warm temperature
    #   D :: cold
    all.climates <- substr(
      lue.frame$clim,
      start=1,
      stop=1
    )
    #
    # Initialize plot preferences:
    par(mfrow=c(2,2))
    par(mar=c(4.5,4.5,1,1))
    # -----------
    # Climate 'A'
    # -----------
    my.clim.a <- which(all.climates == 'A')
    x.a <- seq(from=0, to=lue.frame$maxX[my.clim.a[1]])
    y.a <- x.a*lue.frame$lue[my.clim.a[1]]
    mean.lue.a <- mean(lue.frame$lue[my.clim.a])
    plot(
      x.a, 
      y.a, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=1
    )
    box(lwd=2)
    #
    # Add lines for the remaining climate 'A' plots:
    for(i in seq(from=2, to=length(my.clim.a), by=1)){
      x.a <- seq(from=0, to=lue.frame$maxX[my.clim.a[i]])
      y.a <- x.a*lue.frame$lue[my.clim.a[i]]
      lines(x.a,y.a,lwd=2,col=1)
    }
    #
    # Add mean LUE abline for climate 'A'
    abline(0, mean.lue.a, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Equitorial'),
      lty=1,
      lwd=2,
      col=1,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.a,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # -----------
    # Climate 'B'
    # -----------
    my.clim.b <- which(all.climates == 'B')
    x.b <- seq(from=0, to=lue.frame$maxX[my.clim.b[1]])
    y.b <- x.b*lue.frame$lue[my.clim.b[1]]
    mean.lue.b <- mean(lue.frame$lue[my.clim.b])
    plot(
      x.b, 
      y.b, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=2
    )
    box(lwd=2)
    for(i in seq(from=2, to=length(my.clim.b), by=1)){
      x.b <- seq(from=0, to=lue.frame$maxX[my.clim.b[i]])
      y.b <- x.b*lue.frame$lue[my.clim.b[i]]
      lines(x.b,y.b,lwd=2,col=2)
    }
    abline(0, mean.lue.b, col="grey", lty=2, lwd=2)
    legend(
      'topleft',
      c('Arid'),
      lty=1,
      lwd=2,
      col=2,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.b,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # -----------
    # Climate 'C'
    # -----------
    my.clim.c <- which(all.climates == 'C')
    x.c <- seq(from=0, to=lue.frame$maxX[my.clim.c[1]])
    y.c <- x.c*lue.frame$lue[my.clim.c[1]]  
    mean.lue.c <- mean(lue.frame$lue[my.clim.c])
    plot(
      x.c, 
      y.c, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=3
    )
    box(lwd=2)
    for(i in seq(from=2, to=length(my.clim.c), by=1)){
      x.c <- seq(from=0, to=lue.frame$maxX[my.clim.c[i]])
      y.c <- x.c*lue.frame$lue[my.clim.c[i]]
      lines(x.c,y.c,lwd=2,col=3)
    }
    abline(0, mean.lue.c, col="grey", lty=2, lwd=2)
    legend(
      'topleft',
      c('Temperate'),
      lty=1,
      lwd=2,
      col=3,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.c,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # -----------
    # Climate 'D'
    # -----------
    my.clim.d <- which(all.climates == 'D')
    x.d <- seq(from=0, to=lue.frame$maxX[my.clim.d[1]])
    y.d <- x.d*lue.frame$lue[my.clim.d[1]]
    mean.lue.d <- mean(lue.frame$lue[my.clim.d])
    plot(
      x.d, 
      y.d, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=4
    )
    box(lwd=2)
    for(i in seq(from=2, to=length(my.clim.d), by=1)){
      x.d <- seq(from=0, to=lue.frame$maxX[my.clim.d[i]])
      y.d <- x.d*lue.frame$lue[my.clim.d[i]]
      lines(x.d,y.d,lwd=2,col=4)
    }
    abline(0, mean.lue.d, col="grey", lty=2, lwd=2)
    legend(
      'topleft',
      c('Cold'),
      lty=1,
      lwd=2,
      col=4,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.d,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
  } else if(flag==2){
    # Plot Koppen-Geiger precipitation types:
    #   f :: fully humid
    #   w :: dry winter
    #   m :: monsoonal
    #   S :: steppe
    #   s :: dry summer
    #
    # Get list of precipitation zone codes:
    all.climates <- substr(
      lue.frame$clim,
      start=2,
      stop=2
    )
    #
    # Initialize plot preferences:
    par(mfrow=c(3,2))
    par(mar=c(4.5,4.5,1,1))
    #
    # ----------------------
    # Precipitation zone 'f'
    # ----------------------
    my.precip.f <- which(all.climates == 'f')
    x.f <- seq(from=0, to=lue.frame$maxX[my.precip.f[1]])
    y.f <- x.f*lue.frame$lue[my.precip.f[1]]
    mean.lue.f <- mean(lue.frame$lue[my.precip.f])
    plot(
      x.f, 
      y.f, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=1
    )
    box(lwd=2)
    #
    # Add lines for the remaining climate 'A' plots:
    for(i in seq(from=2, to=length(my.precip.f), by=1)){
      x <- seq(from=0, to=lue.frame$maxX[my.precip.f[i]])
      y <- x*lue.frame$lue[my.precip.f[i]]
      lines(x,y,lwd=2,col=1)
    }
    #
    # Add mean LUE abline for climate 'A'
    abline(0, mean.lue.f, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Humid'),
      lty=1,
      lwd=2,
      col=1,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.f,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # ----------------------
    # Precipitation zone 'w'
    # ----------------------
    my.precip.w <- which(all.climates == 'w')
    x.w <- seq(from=0, to=lue.frame$maxX[my.precip.w[1]])
    y.w <- x.w*lue.frame$lue[my.precip.w[1]]
    mean.lue.w <- mean(lue.frame$lue[my.precip.w])
    plot(
      x.w, 
      y.w, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=2
    )
    box(lwd=2)
    #
    # Add lines for the remaining precip 'w' plots:
    for(i in seq(from=2, to=length(my.precip.w), by=1)){
      x <- seq(from=0, to=lue.frame$maxX[my.precip.w[i]])
      y <- x*lue.frame$lue[my.precip.w[i]]
      lines(x,y,lwd=2,col=2)
    }
    #
    # Add mean LUE abline for precip 'w'
    abline(0, mean.lue.w, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Dry winter'),
      lty=1,
      lwd=2,
      col=2,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.w,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # ----------------------
    # Precipitation zone 'm'
    # ----------------------
    my.precip.m <- which(all.climates == 'm')
    x.m <- seq(from=0, to=lue.frame$maxX[my.precip.m[1]])
    y.m <- x.m*lue.frame$lue[my.precip.m[1]]
    mean.lue.m <- mean(lue.frame$lue[my.precip.m])
    plot(
      x.m, 
      y.m, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=3
    )
    box(lwd=2)
    #
    # Add lines for the remaining precip 'm' plots:
    # NOTE: there's only one 'm'
    #
    # Add mean LUE abline for precip 'm'
    abline(0, mean.lue.m, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Monsoon'),
      lty=1,
      lwd=2,
      col=3,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.m,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # ----------------------
    # Precipitation zone 'S'
    # ----------------------
    my.precip.S <- which(all.climates == 'S')
    x.S <- seq(from=0, to=lue.frame$maxX[my.precip.S[1]])
    y.S <- x.S*lue.frame$lue[my.precip.S[1]]
    mean.lue.S <- mean(lue.frame$lue[my.precip.S])
    plot(
      x.S, 
      y.S, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=4
    )
    box(lwd=2)
    #
    # Add lines for the remaining precip 'S' plots:
    for(i in seq(from=2, to=length(my.precip.S), by=1)){
      x <- seq(from=0, to=lue.frame$maxX[my.precip.S[i]])
      y <- x*lue.frame$lue[my.precip.S[i]]
      lines(x,y,lwd=2,col=4)
    }
    #
    # Add mean LUE abline for precip 'S'
    abline(0, mean.lue.S, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Steppe'),
      lty=1,
      lwd=2,
      col=4,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.S,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    # ----------------------
    # Precipitation zone 's'
    # ----------------------
    my.precip.s <- which(all.climates == 's')
    x.s <- seq(from=0, to=lue.frame$maxX[my.precip.s[1]])
    y.s <- x.s*lue.frame$lue[my.precip.s[1]]
    mean.lue.s <- mean(lue.frame$lue[my.precip.s])
    plot(
      x.s, 
      y.s, 
      type='l', 
      lwd=2,
      xlab=expression(fPAR~x~PPFD~(mol%.%m^{-2})),
      ylab=expression(GPP~(mol%.%m^{-2})),
      xlim=c(0,max.x),
      ylim=c(0,max.y),
      col=5
    )
    box(lwd=2)
    #
    # Add lines for the remaining precip 's' plots:
    for(i in seq(from=2, to=length(my.precip.s), by=1)){
      x <- seq(from=0, to=lue.frame$maxX[my.precip.s[i]])
      y <- x*lue.frame$lue[my.precip.s[i]]
      lines(x,y,lwd=2,col=5)
    }
    #
    # Add mean LUE abline for precip 's'
    abline(0, mean.lue.s, col="grey", lty=2, lwd=2)
    #
    # Add legend for climate name and mean LUE:
    legend(
      'topleft',
      c('Dry summer'),
      lty=1,
      lwd=2,
      col=5,
      bty='n'
    )
    legend(
      'topleft',
      c('',format(mean.lue.s,digits=3)),
      lty=c(0,2),
      lwd=c(0,2),
      col=c('grey','grey'),
      bty='n'
    )
    #
  }
    #
}

# ************************************************************************
# * Name: monthly.lue.analyzer
# *
# * Input: data frame (my.obj)
# *        boolean for directory (mac)
# *
# * Return: data frame (my.obj)
# *
# * Features: This script adds additional columns to the input data
# *           frame for monthly LUEmo=(GPP/(fAPAR*PPFD)), station 
# *           elevation, mean growing season temperature, mean growing
# *           season VPD, year
# ************************************************************************
monthly.lue.analyzer <- function(my.obj, mac){
  # Get monthly LUE, MGS Tc, MGS VPD, and elevation for analysis
  # * monthly LUE is the slope between the origin and monthly points of 
  #   (GPP, fAPAR*PPFD)
  # * MGS (mean growing season) data can be found in an external file
  #
  if (mac){
    met.path <- paste(
      "/Users/twdavis/Dropbox/", # Mac
      "Work/Imperial/GePPiSaT/gepisat/psql_data/met_data/flux/",
      sep=""
    )
    mgs.path <- paste(
      "/Users/twdavis/Dropbox/", # Mac
      "Work/Imperial/flux/results/2002-06/time-series_analysis/",
      sep=""
    )
  } else {
    met.path <- paste(
      "/home/user/Dropbox/",     # Linux
      "Work/Imperial/GePPiSaT/gepisat/psql_data/met_data/flux/",
      sep=""
    )
    mgs.path <- paste(
      "/home/user/Dropbox/",     # Linux
      "Work/Imperial/flux/results/2002-06/time-series_analysis/",
      sep=""
    )
  }
  #
  # Get vegetation type from meta data:
  met.file = list.files(
    path = met.path, 
    pattern = "Fluxdata_Meta-Data.csv"
  )
  met.content <- read.csv(
    paste(met.path,met.file,sep=""), 
    header=T
  )
  #
  # Get mean growing season data
  mgs.file <- list.files(
    path = mgs.path,
    pattern = "*.csv"
  )
  mgs.content <- read.csv(
    paste(mgs.path, mgs.file, sep=""),
    header=T
  )
  #
  # Calculate monthly LUE:
  # GPP = LUE * fAPAR * PPFD  =>  LUE = GPP / (fAPAR * PPFD)
  my.obj$LUEmo <- my.obj$GPP.mol_m2 / (my.obj$fAPAR * my.obj$PPFD.mol_m2)
  #
  # Add year column to my.obj for correct assignment of MGS
  my.obj$year <- as.numeric(substr(as.character(my.obj$Timestamp),1, 4))
  #
  # Add MGS and elevation data to my.obj
  all.sites <- as.character(mgs.content$site)
  all.years <- seq(from=2002, to=2006, by=1)
  my.obj$elv <- 0*seq(from=1, to=dim(my.obj)[1], by=1)
  my.obj$veg <- 0*seq(from=1, to=dim(my.obj)[1], by=1)
  my.obj$mgs.tc <- 0*seq(from=1, to=dim(my.obj)[1], by=1)
  my.obj$mgs.vpd <- 0*seq(from=1, to=dim(my.obj)[1], by=1)
  #
  for (site in all.sites){
    # Assign vegetation type:
    my.obj$veg[which(my.obj$Station == site)] <- as.character(
      met.content[which(met.content$stationid == site),'classid']
    )
    # Assign elevation:
    my.obj$elv[which(my.obj$Station == site)] <- mgs.content$ele.m[which(mgs.content$site == site)]
    #
    # Assign annual Tc and VPD:
    for (year in all.years){
      if (year == 2002){
        if (length(my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)]) > 0){
          my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$tc.2002[which(mgs.content$site == site)]
          my.obj$mgs.vpd[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$vpd.2002[which(mgs.content$site == site)]
        }
      } else if (year == 2003){
        if (length(my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)]) > 0){
          my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$tc.2003[which(mgs.content$site == site)]
          my.obj$mgs.vpd[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$vpd.2003[which(mgs.content$site == site)]
        }
      } else if (year == 2004){
        if (length(my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)]) > 0){
          my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$tc.2004[which(mgs.content$site == site)]
          my.obj$mgs.vpd[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$vpd.2004[which(mgs.content$site == site)]
        }
      } else if (year == 2005){
        if (length(my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)]) > 0){
          my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$tc.2005[which(mgs.content$site == site)]
          my.obj$mgs.vpd[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$vpd.2005[which(mgs.content$site == site)]
        }
      } else if (year == 2006){
        if (length(my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)]) > 0){
          my.obj$mgs.tc[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$tc.2006[which(mgs.content$site == site)]
          my.obj$mgs.vpd[which(my.obj$Station == site & my.obj$year == year)] <- 
            mgs.content$vpd.2006[which(mgs.content$site == site)]
        }
      }
    } #/ END YEAR LOOP /#
  } #/ END SITE LOOP /#
  #
  # Bin vegetation classes together (based on assumed height)
  my.obj$vbin <- 0 + 
    1*(my.obj$veg == 'EBF') + 
    1*(my.obj$veg == 'MF') + 
    1*(my.obj$veg == 'ENF') + 
    1*(my.obj$veg == 'DBF')
  #
  # Bin air temperature (5-degree breaks):
  my.obj$tbin <- 0 +
    1*(my.obj$Tc.deg_C >= 5 & my.obj$Tc.deg_C < 10) +
    2*(my.obj$Tc.deg_C >= 10 & my.obj$Tc.deg_C < 15) +
    3*(my.obj$Tc.deg_C >= 15 & my.obj$Tc.deg_C < 20) +
    4*(my.obj$Tc.deg_C >= 20 & my.obj$Tc.deg_C < 25) +
    5*(my.obj$Tc.deg_C >= 25)
  #
  # Return object:
  my.obj
}

# ************************************************************************
# * Name: f.of.alpha
# *
# * Input: Cramer-Prentice alpha
# *        exponent k (optional)
# *
# * Return: f_alpha
# *
# * Features: This script processes the Cramer-Prentice alpha (0,1) 
# *           using a power function with an exponent between 0 and 1
# ************************************************************************
f.of.alpha <- function(alpha, k=0.1){
  f_alpha <- alpha^(k)
  f_alpha
}

# ************************************************************************
# * Name: monthly.lue.regression
# *
# * Input: x <- regressor (e.g., Tc, VPD, Elv)
# *        y <- monthly LUE
# *        x.text <- expression (x axis label)
# *        y.text <- expression (y axis label)
# *
# * Return: None.
# *
# * Features: This function plots monthly LUE versus a x-axis variable
# ************************************************************************
monthly.lue.regression <- function(x, y, x.text, y.text){
  par(mar=c(4.5,4.5,1,1))
  plot(
    x,
    y,
    type='p',
    pch=20,
    col='gray',
    xlab=NA,
    ylab=NA,
    axes=F
  )
  box(lwd=2)
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  mtext(side=1, x.text, line=2)
  mtext(side=2, y.text, line=3)
  #
  fit <- lm(y ~ x)
  abline(fit, col='red',lwd=1.5)
  #
  fit.coef <- as.numeric(fit$coefficients[2])
  n <- length(fit$residuals)
  r2 <- summary(fit)$adj.r.squared
  my.p <- summary(fit)$coefficients[2,4]
  #
  rp = vector('expression',4)
  rp[1] = substitute(
    expression(italic(m) == MYLUE),
    list(MYLUE = format(fit.coef,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(n) == MYN),
    list(MYN = format(n,dig=2)))[2]
  rp[3] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  rp[4] = substitute(
    expression(italic(p) == MYOTHERVALUE), 
    list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  legend('topleft', legend = rp, bty = 'n')
}

# ************************************************************************
# * Name: binned.lue.regression
# *
# * Input: my.obj <- content object (requires lue_conversions)
# *        bin_v <- vegetation bin number (0:short, 1:tall)
# *        bin_t <- temperature bin number 
# *                 0 :  0- 5 deg C
# *                 1 :  5-10 deg C
# *                 2 : 10-15 deg C
# *                 3 : 15-20 deg C
# *                 4 : 20-25 deg C
# *                 5 : 25-30 deg C
# *        pars <- minpack.lm object
# *
# * Return: None.
# *
# * Features: This function plots GPP versus PHI_O * m' * Iabs
# ************************************************************************
binned.lue.regression <- function(my.obj, bin_v, bin_t, pars){
  # Fit coefficients
  #PHI_O <- pars$par$a
  B_OVER_A <- pars$par$b
  #
  # Get binned x and y points:
  my.points <- which(
    my.obj$vbin == bin_v &
      my.obj$tbin == bin_t
  )
  #
  # Reduce object to selected bins:
  my.obj <- my.obj[my.points,]
  #
  # Calculate x & y:
  x <- PHI_O * 
    my.obj$PPFD.mol_m2 * 
    my.obj$fAPAR * 
    calc_m(my.obj, sqrt(B_OVER_A))
  y <- my.obj$GPP.mol_m2
  #
  # Get max plotting extents:
  max.y <- (
    floor(max(y[which(!is.na(x))])+5) - 
      (floor(max(y[which(!is.na(x))])+5) %% 5)
  )
  max.x <- (
    floor(max(x, na.rm=T)+2) -
      (floor(max(x, na.rm=T)+2) %% 2)
  )
  #
  # Set bin options to subtitle expression:
  veg_height = ""
  if (bin_v == 0){
    veg_height = "short"
  } else if (bin_v == 1) {
    veg_height = "tall"
  }
  tc_bin = list(a=NA,b=NA)
  if (bin_t == 0){
    tc_bin$a = 0
    tc_bin$b = 5
  } else if (bin_t == 1) {
    tc_bin$a = 5
    tc_bin$b = 10
  } else if (bin_t == 2) {
    tc_bin$a = 10
    tc_bin$b = 15
  } else if (bin_t == 3) {
    tc_bin$a = 15
    tc_bin$b = 20
  } else if (bin_t == 4) {
    tc_bin$a = 20
    tc_bin$b = 25
  } else if (bin_t == 5) {
    tc_bin$a = 25
    tc_bin$b = 30
  }
  (subtitle <- substitute(
    list(H[veg] == MYHEIGHT,T[c] == MYA-MYB~degree*C),
    list(MYHEIGHT = veg_height, MYA = tc_bin$a, MYB = tc_bin$b)
  ))
  # Calc stats and create plot legend:
  n <- length(pars$fvec)
  #lue.p <- summary(pars)$coefficients[1,4] # \ These were available with minpack.lm
  #boa.p <- summary(pars)$coefficients[2,4] # / nls.lm object.
  sxy <- sum((y[which(!is.na(x))]-mean(y[which(!is.na(x))]))*(x[which(!is.na(x))]-mean(x[which(!is.na(x))])))
  sxx <- sum((x[which(!is.na(x))]-mean(x[which(!is.na(x))]))^2)
  syy <- sum((y[which(!is.na(x))]-mean(y[which(!is.na(x))]))^2)
  r <- sxy / sqrt(sxx) / sqrt(syy)
  r2 <- r^2
  rp = vector('expression',4)
  rp[1] = substitute(
    expression(italic(phi[o]) == MYLUE),
    list(MYLUE = format(PHI_O,dig=3)))[2]
  rp[2] = substitute(
    expression(italic(b/a) == MYBA),
    list(MYBA = format(B_OVER_A,dig=2)))[2]
  rp[3] = substitute(
    expression(italic(n) == MYN),
    list(MYN = format(n,dig=2)))[2]
  rp[4] = substitute(
    expression(italic(r)^2 == MYVALUE), 
    list(MYVALUE = format(r2,dig=3)))[2]
  #
  # Create the plot:
  par(mar=c(5.25,4.5,1,1))
  plot(
    x, 
    y, 
    xlim=c(0,max.x),
    ylim=c(0,max.y),
    type='p',
    pch=20,
    col='gray',
    xlab=NA,
    ylab=NA,
    axes=F
  )
  box(lwd=2)
  abline(a=0,b=1, col='black',lwd=1.5, lty=2) # 1:1 line
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, las=1, lwd=0, line=-0.4)
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, las=1, lwd=0, line=-0.4)
  #mtext(side=1, expression(phi[o]%.%italic(m*minute)%.%I[abs]~(mol%.%m^{-2})), line=2)
  mtext(side=1, expression(phi[o]%.%italic(m)%.%I[abs]~(mol%.%m^{-2})), line=2)
  mtext(side=2, expression(GPP~(mol%.%m^{-2})), line=2)
  mtext(side=1, as.expression(subtitle), line=3)
  legend('topleft', legend = rp, bty = 'n')
}

# ************************************************************************
# * Name: lueResid
# *
# * Input: par <- list of parameters (a: PHI_O, b: b/a)
# *        my.obj <- content object (requires lue_conversions)
# *
# * Return: vector of residuals
# *
# * Features: This function calculates GPP - (PHI_O * m' * Iabs), where
# *           NA values (due to the calculation of m') are set equal to
# *           zero (should not cause issue with minpack.lm, which 
# *           minimizes the sum of the squared residuals)
# *
# * Depends: calc_m
# * 
# * @TODO: change the minimization function to 1 - r, where r is the 
# *        Pearson correlation coefficient (sign of linearity)
# ************************************************************************
lueResid <- function(par, my.obj){
  # Calculate Pearson's correlation coefficient between GPP and 
  # PHI_O*fAPAR*m*PPFD
  y <- my.obj$GPP.mol_m2
  x <- PHI_O * my.obj$PPFD.mol_m2 * my.obj$fAPAR * calc_m(my.obj, par$b)
  sxy <- sum((y[which(!is.na(x))]-mean(y[which(!is.na(x))]))*(x[which(!is.na(x))]-mean(x[which(!is.na(x))])))
  sxx <- sum((x[which(!is.na(x))]-mean(x[which(!is.na(x))]))^2)
  syy <- sum((y[which(!is.na(x))]-mean(y[which(!is.na(x))]))^2)
  r <- sxy / sqrt(sxx) / sqrt(syy)
  #
  1-r
}

# ************************************************************************
# * Name: climate_class
# *
# * Input: - char, list of stations (my.stations)
# *        - char, directory for meta data
# *
# * Return: data frame, station climate and classes
# *
# * Features: This function maps station names to their climate and 
# *           vegetation class
# *
# * Depends: get_meta
# * 
# ************************************************************************
climate_class <- function(my.stations, mdir){
  # Create data frame for climate and class variables:
  num.stations <- length(my.stations)
  my.data <- matrix(nrow = num.stations, ncol = 2)
  my.data <- as.data.frame(my.data)
  names(my.data) <- c('climate', 'class')
  rownames(my.data) <- my.stations
  #
  for (station in my.stations){
    my.meta <- get_meta(mdir, station)
    my.data[station, 'climate'] <- as.character(my.meta$climate)
    my.data[station, 'class'] <- as.character(my.meta$class)
  }
  #
  return(my.data)
}

# ************************************************************************
# * Name: modeled_lue
# *
# * Input: - data frame, LUE contents (my.obj)
# *        - data frame, climate/class IDs (my.ids)
# *
# * Return: data frame, station LUE w/ climate and classes
# *
# * Features: This function calculate LUE via linear regression and 
# *           saves the output, along with climate and class info
# * 
# ************************************************************************
modeled_lue <- function(my.obj, my.ids){
  my.stations <- levels(unique(my.obj$Station))
  my.data <- matrix(nrow = length(my.stations), ncol = 5)
  my.data <- as.data.frame(my.data)
  names(my.data) <- c('lue', 'lue_err', 'r2', 'climate', 'class')
  rownames(my.data) <- my.stations
  #
  for(station in my.stations){
    # Linear regression: GPP v. fAPAR x PPFD
    my.idx <- which(as.character(my.obj$Station) == station)
    my.y <- my.obj$GPP.mol_m2[my.idx]
    my.var <- (1/my.obj$GPP_err[my.idx])
    my.x <- (my.obj$fAPAR[my.idx])*(my.obj$PPFD.mol_m2[my.idx])
    my.fit <- lm(my.y ~ 0 + my.x, weights = my.var)
    #
    my.lue <- summary(my.fit)$coefficients['my.x', 'Estimate']
    my.lue.err <- summary(my.fit)$coefficients['my.x', 'Std. Error']
    my.r2 <- summary(my.fit)$adj.r.squared
    #
    my.data[station, 'lue'] <- my.lue
    my.data[station, 'lue_err'] <- my.lue.err
    my.data[station, 'r2'] <- my.r2
    my.data[station, 'climate'] <- my.ids[station, 'climate']
    my.data[station, 'class'] <- my.ids[station, 'class']
  }
  return(my.data)
}

# ************************************************************************
# * Name: modeled_lue_alt
# *
# * Input: - data frame, LUE contents (my.obj)
# *        - data frame, climate/class IDs (my.ids)
# *
# * Return: data frame, station LUE w/ CPA correction, climate and 
# *         classes
# *
# * Features: This function calculates LUE with CPA correction via linear 
# *           regression and saves the output, along with climate and 
# *           class info
# * 
# ************************************************************************
modeled_lue_alt <- function(my.obj, my.ids){
  my.stations <- levels(my.obj$Station)
  my.data <- matrix(nrow = length(my.stations), ncol = 5)
  my.data <- as.data.frame(my.data)
  names(my.data) <- c('lue', 'lue_err', 'r2', 'climate', 'class')
  rownames(my.data) <- my.stations
  #
  for(station in my.stations){
    # Linear regression: GPP v. fAPAR x PPFD
    my.idx <- which(as.character(my.obj$Station) == station)
    my.y <- my.obj$GPP.mol_m2[my.idx]
    my.var <- (1/my.obj$GPP_err[my.idx])
    my.x <- (my.obj$fAPAR[my.idx])*(my.obj$PPFD.mol_m2[my.idx])*(
      my.obj$ALPHA[my.idx]/1.26)^(1/4)
    my.fit <- lm(my.y ~ 0 + my.x, weights = my.var)
    #
    my.lue <- summary(my.fit)$coefficients['my.x', 'Estimate']
    my.lue.err <- summary(my.fit)$coefficients['my.x', 'Std. Error']
    my.r2 <- summary(my.fit)$adj.r.squared
    #
    my.data[station, 'lue'] <- my.lue
    my.data[station, 'lue_err'] <- my.lue.err
    my.data[station, 'r2'] <- my.r2
    my.data[station, 'climate'] <- my.ids[station, 'climate']
    my.data[station, 'class'] <- my.ids[station, 'class']
  }
  return(my.data)
}

# ************************************************************************
# * Name: bp_modeled_lue
# *
# * Input: - data frame, LUE contents (my.obj)
# *        - char, path to meta data file (mdir)
# *
# * Return: None.
# *
# * Features: This function creates two box plots of modeled LUE
# *           (1) with respect to climate groups
# *           (2) with respect to vegetation classes
# *
# * Depends: - climate_class, function
# *          - modeled_lue, function
# *          - lattice, library
# * 
# ************************************************************************
bp_modeled_lue <- function(my.obj, mdir){
  my.stations <- levels(my.obj$Station)
  my.ids <- climate_class(my.stations, mdir)
  my.lue <- modeled_lue(my.obj, my.ids)
  my.climates <- unique(my.lue$climate)
  my.classes <- unique(my.lue$class)
  #
  bwplot(
    lue ~ climate,
    data = my.lue,
    xlab = list(label = 'Climate', cex = 1.5),
    ylab = list(label = 'LUE', cex=1.5),
    scales = list(cex = 1.5, font = 1)
  )
  #
  bwplot(
    lue ~ class,
    data = my.lue,
    xlab = list(label = 'Vegetation Class', cex = 1.5),
    ylab = list(label = 'LUE', cex = 1.5),
    scales = list(cex = 1.5)
  )
}

# ************************************************************************
# * Name: bp_lue_dr2
# *
# * Input: - data frame, LUE contents (my.obj)
# *        - char, path to meta data file (mdir)
# *
# * Return: None.
# *
# * Features: This function creates two box plots of LUE fits 
# *           (w/ and w/o CPA correction)
# *           (1) with respect to climate groups
# *           (2) with respect to vegetation classes
# *
# * Depends: - climate_class, function
# *          - modeled_lue, function
# *          - modeled_lue_alt, function
# *          - lattice, library
# * 
# ************************************************************************
bp_lue_dr2 <- function(my.obj, mdir, opt){
  my.stations <- levels(my.obj$Station)
  my.ids <- climate_class(my.stations, mdir)
  my.lue <- modeled_lue(my.obj, my.ids)
  my.lue.alt <- modeled_lue_alt(my.obj, my.ids)
  my.climates <- unique(my.lue$climate)
  my.classes <- unique(my.lue$class)
  #
  my.data <- matrix(nrow = length(my.stations), ncol = 5)
  my.data <- as.data.frame(my.data)
  names(my.data) <- c('r2_orig', 'r2_cpa', 'diff_r2', 'climate', 'class')
  row.names(my.data) <- my.stations
  for (station in my.stations){
    my.data[station, 'r2_orig'] <- my.lue[station, 'r2']
    my.data[station, 'r2_cpa'] <- my.lue.alt[station, 'r2']
    my.data[station, 'class'] <- my.lue[station, 'class']
    my.data[station, 'climate'] <- my.lue[station, 'climate']
  }
  my.data$diff_r2 <- (my.data$r2_cpa - my.data$r2_orig)
  #
  if (opt == 1){
    bwplot(
      diff_r2 ~ climate,
      data = my.data,
      xlab = list(label = 'Climate', cex = 1.5),
      ylab = list(label = expression(Delta*R^2), cex=1.5),
      scales = list(cex = 1.5, font = 1)
    )
  } else {
    bwplot(
      diff_r2 ~ class,
      data = my.data,
      xlab = list(label = 'Vegetation Class', cex = 1.5),
      ylab = list(label = expression(Delta*R^2), cex = 1.5),
      scales = list(cex = 1.5)
    )
  }
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### CONSTANTS ################################################################
# /////////////////////////////////////////////////////////////////////////////
mac = FALSE

# Read directories:
if(mac){
  meta.path <- paste(
    "/Users/twdavis/Dropbox/",     # Mac
    "Work/Imperial/flux/data/psql-data/flux/",
    sep=""  
  )
  lue.path <- paste(
    "/Users/twdavis/Dropbox/",       # Mac
    "Work/Imperial/flux/results/2002-06/lue/",
    sep=""
  )
  ts.path <- paste(
    "/Users/twdavis/Dropbox/", 
    "Work/Imperial/flux/results/2002-06/time-series_analysis/",
    sep=""
  )
} else {
  meta.path <- paste(
    "/home/user/Dropbox/",         # Linux
    "Work/Imperial/flux/data/psql-data/flux/",
    sep=""  
  )
  lue.path <- paste(
    "/home/user/Dropbox/",           # Linux
    "Work/Imperial/flux/results/2002-06/lue/",
    sep =""
  )  
  ts.path <- paste(
    "/home/user/Dropbox/",
    "Work/Imperial/flux/results/2002-06/time-series_analysis/",
    sep=""
    )
}

PHI_O <- 0.093

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### MAIN #####################################################################
# /////////////////////////////////////////////////////////////////////////////
lue.all.file <- list.files(path=lue.path, pattern = "*_All_Data-26.txt")
content <- read.csv(paste(lue.path, lue.all.file, sep=""), header=T)

# Backup content and filter & convert:
content.bak <- content
content <- filter_contents(content)
content <- lue_conversions(content)

# If you are partitioning the All Data file:
all.stations <- levels(unique(content$Station))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SITE AVERAGE LUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# written 2014-10-07
# station average LUE; weighted linear regression
# prepared for Wang Han
content <- content.bak
content$Iabs.mol_m2 <- content$fAPAR*content$PPFD.mol_m2

ts.file <- list.files(path=ts.path, pattern = "^timeseries_v26.*txt")
ts.data <- read.csv(paste(ts.path, ts.file, sep=""), header=T)
junk.ts <- do.call(rbind, sapply(as.character(ts.data$Month), strsplit, split="-"))
ts.data$year <- as.numeric(junk.ts[,1])
ts.data$month <- as.numeric(junk.ts[,2])
rm(junk.ts)

num_rows <- length(all.stations)
output.data <- matrix(nrow=num_rows, ncol=5)
output.data <- as.data.frame(output.data)
names(output.data) <- c('Station', 'St_lon', 'St_lat', 'LUE', 'LUE_err')
output.data$Station <- all.stations

for (i in seq(num_rows)){
  my.station <- output.data[i, 'Station']
  j <- which(ts.data[,'Station'] == my.station)[1]
  output.data[i, 'St_lon'] <- ts.data[j, 'St_Lon']
  output.data[i, 'St_lat'] <- ts.data[j, 'St_Lat']
  #
  k <- which(as.character(content[,'Station']) == my.station)
  if (length(k) > 2){
    y <- content[k, 'GPP.mol_m2']
    x <- content[k, 'Iabs.mol_m2']
    my.weights <- 1/(content[k, 'GPP_err']*max(1/content[k, 'GPP_err']))
    fit <- lm(y ~ 0 + x, weights=my.weights)
    output.data[i, 'LUE'] <- summary(fit)$coefficients[1, 'Estimate']
    output.data[i, 'LUE_err'] <- summary(fit)$coefficients[1, 'Std. Error']
    rm(x, y, my.weights, fit)
  } else {
    output.data[i, 'LUE'] <- -9999
    output.data[i, 'LUE_err'] <- -9999
  }
  rm(my.station, j, k)
}
rm(ts.data)
rm(i, num_rows)

output.file <- paste(lue.path, 'LUE_station_ave-24.txt', sep="")
write.table(x=output.data,
            file=output.file,
            sep=",",
            col.names=T,
            row.names=F,
            quote=F,
            eol='\n')
rm(output.file, output.data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MONTHLY LUE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# written 2014-10-07
# monthly LUE (i.e., GPP/Iabs) and climatology
# prepared for Wang Han
content <- content.bak
content$LUE <- content$GPP.mol_m2/(content$fAPAR*content$PPFD.mol_m2)
content$LUE_err <- content$GPP_err/(content$fAPAR*content$PPFD.mol_m2)
junk.ts <- do.call(rbind, sapply(as.character(content$Timestamp), strsplit, split="-"))
content$year = as.numeric(junk.ts[,1])
content$month = as.numeric(junk.ts[,2])
rm(junk.ts)
as.data.frame(cbind(as.character(content$Timestamp),content$years,content$months))

ts.file <- list.files(path=ts.path, pattern = "*26.txt")
ts.data <- read.csv(paste(ts.path, ts.file, sep=""), header=T)
junk.ts <- do.call(rbind, sapply(as.character(ts.data$Month), strsplit, split="-"))
ts.data$year <- as.numeric(junk.ts[,1])
ts.data$month <- as.numeric(junk.ts[,2])
rm(junk.ts)

num_rows <- dim(content)[1]
output.data <- matrix(nrow=num_rows, ncol=13)
output.data <- as.data.frame(output.data)
names(output.data) <- c('Station', 'Year', 'Month', 'St_lon', 'St_lat', 
                        'Grid_lon', 'Grid_lat', 'LUE', 'LUE_err', 'Tair_C',
                        'VPD_Pa', 'SF', 'Pre_mm')
output.data$Station <- as.character(content$Station)
output.data$Year <- content$year
output.data$Month <- content$month
output.data$LUE <- content$LUE
output.data$LUE_err <- content$LUE_err

for (i in seq(num_rows)){
  my.station <- as.character(output.data[i, 'Station'])
  my.year <- output.data[i, 'Year']
  my.month <- output.data[i, 'Month']
  #
  j <- which(as.character(ts.data$Station) == my.station & 
                    ts.data$year == my.year & 
                    ts.data$month == my.month)
  #
  if (any(j) & length(j) == 1){
    output.data[i, 'St_lon'] <- ts.data[j, 'St_Lon']
    output.data[i, 'St_lat'] <- ts.data[j, 'St_Lat']
    output.data[i, 'Grid_lon'] <- ts.data[j, 'Grid_Lon']
    output.data[i, 'Grid_lat'] <- ts.data[j, 'Grid_Lat']
    output.data[i, 'Tair_C'] <- ts.data[j, 'Tair_C']
    output.data[i, 'VPD_Pa'] <- (1e3)*ts.data[j, 'VPD_kPa']
    output.data[i, 'SF'] <- ts.data[j, 'SF']
    output.data[i, 'Pre_mm'] <- ts.data[j, 'Pre_mm']
  }
}
rm(i, j, num_rows)
rm(my.month, my.year, my.station)

output.file <- paste(lue.path, 'LUE_monthly_ts.txt', sep="")
write.table(x=output.data,
            file=output.file,
            sep=",",
            col.names=T,
            row.names=F,
            quote=F,
            eol='\n')
rm(ts.data, ts.file)
rm(output.file, output.data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DELTA R SQUARED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# written 2014-11-16
# difference in model fits of station averaged LUE; weighted linear regression
content <- content.bak
content$Iabs.mol_m2 <- content$fAPAR*content$PPFD.mol_m2
my_clims <- climate_class(all.stations, meta.path)

num_rows <- length(all.stations)
output.data <- matrix(nrow=num_rows, ncol=5)
output.data <- as.data.frame(output.data)
names(output.data) <- c('Station', 'R2_orig', 'R2_cpa', 'climate', 'class')
output.data$Station <- all.stations
output.data$climate <- my_clims[, 'climate']
output.data$class <- my_clims[, 'class']
rm(my_clims)

for (i in seq(num_rows)){
  my.station <- output.data[i, 'Station']
  k <- which(as.character(content[,'Station']) == my.station)
  if (length(k) > 2){
    y <- content[k, 'GPP.mol_m2']
    x_a <- content[k, 'Iabs.mol_m2']
    x_b <- x_a*(content[k, 'ALPHA']^(1/4))
    #
    my.weights <- 1/(content[k, 'GPP_err'])
    #
    fit_a <- lm(y ~ 0 + x_a, weights=my.weights)
    fit_b <- lm(y ~ 0 + x_b, weights=my.weights)
    #
    r2_a <- summary(fit_a)$adj.r.squared
    r2_b <- summary(fit_b)$adj.r.squared
    #
    output.data[i, 'R2_orig'] <- r2_a
    output.data[i, 'R2_cpa'] <- r2_b
    #
    rm(x_a, x_b, y, my.weights, fit_a, fit_b, r2_a, r2_b)
  } else {
    output.data[i, 'R2_orig'] <- -9999
    output.data[i, 'R2_cpa'] <- -9999
  }
  rm(my.station, k)
}
rm(i, num_rows)

output.data$DeltaR2 <- (output.data$R2_cpa - output.data$R2_orig)

output.file <- paste(lue.path, 'LUE_R2_all_stations-26.txt', sep="")
write.table(x=output.data,
            file=output.file,
            sep=",",
            col.names=T,
            row.names=F,
            quote=F,
            eol='\n')
rm(ts.data, ts.file)
rm(output.file, output.data)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ALPHA HISTOGRAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
content <- content.bak
my.breaks <- seq(from=0.1, to=1.3, by=0.2)
par(mar=c(4,4.5,1,4))
hist(content$ALPHA, my.breaks, xlim=c(0, 1.4), ylim=c(0, 2000), main="", 
     xlab='', col='gray50', labels=T, axes=F)
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, expression(Cramer~Prentice~alpha), line=2)
mtext(side=2, 'Frequency', line=3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INDIVIDUAL STATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# For storing m and b estimates:
# mb.est <- matrix(nrow = length(all.stations), ncol=3)
# mb.est <- as.data.frame(mb.est)
# names(mb.est) <- c("m","b","D")
# rownames(mb.est) <- all.stations

# Define output file:
ps.out.name <- "LUE_Analysis_Ind-Stations_VPD.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_VPD-sig.ps"

ps.out.name <- "LUE_Analysis_Ind-Stations_T.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_T-sig.ps"

ps.out.name <- "LUE_Analysis_Ind-Stations_ALPHA-26.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_ALPHA_alt-26.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_ALPHA-sig.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_ALPHA-pow.ps"

ps.out.name <- "LUE_Analysis_Ind-Stations_T-ALPHA.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_T-ALPHA-sig.ps"

ps.out.name <- "LUE_Analysis_Ind-Stations_R2.ps"
ps.out.name <- "LUE_Analysis_Ind-Stations_R2-sig.ps"

ps.out.name <- "LUE_Analysis_Ind_Stations_X-Y_af.ps"

out.path <- lue.path
postscript(
  file = paste(out.path,ps.out.name, sep=""), 
  width = 6, 
  height = 4, 
  title = ps.out.name,
  paper = "special", 
  horizontal = FALSE,
  onefile = TRUE)

# Iterate through files:
for (stationid in all.stations){
  # Reset content:
  content <- lue_conversions(content.bak)
  #
  # Get indexes for current station id (if using All Data file):
  station.idx <- which(content$Station == stationid)
  content <- content[station.idx,]
  #
  # Filter data:
  content <- filter_contents(content)
  #
  # Process if enough data is available:
  if (dim(content)[1] > 2){
  #if (dim(content)[1] > 2 & climate_class(stationid, meta.path)$class %in% c('DBF', 'EBF', 'ENF', 'GRA', 'MF')){
    # ~~~~~~~~~~
    # Basic LUE:
    # ~~~~~~~~~~
    my.x <- (content$PPFD.mol_m2)*(content$fAPAR)*(content$ALPHA/1.26)^(0.25)
    x.label <- expression(list(fAPAR%*%PPFD%*%(alpha/1.26)^0.25, mol%.%m^{-2}))
    #
    # ~~~~~
    # Tsig:
    # ~~~~~
    #my.x <- content$PPFD.mol_m2 * content$fAPAR * content$kTc
    #x.label <- expression(paste(italic(T)[italic(sig)], 
    #                      " x fAPAR x PPFD (", 
    #                      mol%.%m^{-2},
    #                      ")"))
    #
    # ~~~~~
    # Apow:
    # ~~~~~
    #my.x <- content$PPFD.mol_m2 * content$fAPAR * content$ALPHA^{0.25}
    #x.label <- expression(alpha^{0.25}~x~fAPAR~x~PPFD~(mol%.%m^{-2}))
    #
    # ~~~~~~~~~~~~~~~~~
    # New x & y method:
    # ~~~~~~~~~~~~~~~~~
    #plot_xy(content, method='a', basic=F, tau=F, alpha=F, k=0, pos='topleft', 
    #        stid=stationid)
    #
    #plot_vpd(content, my.x, x.label, stationid)
    #plot_temp(content, my.x, x.label, stationid)
    palette(c('firebrick1', 'darkorange', 'gold', 'chartreuse', 'cyan', 'deepskyblue', 'blue'))
    plot_alpha(content, my.x, x.label, meta.path, stationid)
    palette('default')
    #plot_temp_alpha(content, my.x, x.label, stationid)
    #plot_r2(content, my.x, x.label, lue.path, stationid)
    #
  }
}
dev.off()

palette(c('firebrick1', 'darkorange', 'gold', 'chartreuse', 'cyan', 'deepskyblue', 'blue'))
plot(seq(7), pch=19, col=seq(7))
plot_mbD(mb.est)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
 #### LUE ANALYZER ###########################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Create new blank lue analyzer data frame:
my.matrix <- matrix(nrow = length(all.stations), ncol=10)
my.lue.analyzer <- as.data.frame(my.matrix, row.names=all.stations)
names(my.lue.analyzer) <- c(
  "lue",  "lue_err", "maxX",   "maxY",
  "clim", "elv",     "mgs.tc", "mgs.vpd", 
  "patm", "co2"
)
rm(my.matrix)

for (stationid in all.stations){
  # Reset content:
  content <- lue_conversions(content.bak)
  #
  # Get indexes for current station id (if using All Data file):
  station.idx <- which(content$Station == stationid)
  content <- content[station.idx,]
  #
  # Filter data:
  content <- filter_contents(content)
  #
  # Update lue.analyzer
  my.lue.analyzer <- lue.analyzer(my.lue.analyzer, 
                                  content, 
                                  meta.path, 
                                  stationid, 
                                  mac)
  #
  # Be tidy:)
  rm(station.idx)
}
plot.lue.analyzer(my.lue.analyzer, flag=2)
plot_lue_xy(my.lue.analyzer, pos='topleft', stid='Mean Growing Season LUE')

# X-Y analysis:
# @TODO
# - calc_gstar(my.lue.analyzer$mgs.tc)
# - calc_k
# - calc_n


# ~~~~~~~~~~~~ MONTHLY LUE AS PREDICTED BY ELV, TC, AND VPD ~~~~~~~~~~~~~
# Monthly LUE is determined by the slope from the origin (0,0) to the 
# point (GPP, fAPAR*PPFD), i.e., LUE = A/Iabs
content <- filter_contents(content.bak)
content <- filter_stations(content)
content <- lue_conversions(content)
my.lue.mo.analyzer <- monthly.lue.analyzer(content, mac)
my.lue.mo.analyzer$elv[which(my.lue.mo.analyzer$elv < -5000)] <- NA

# ************* Binned LUE *************
# Create an output file for holding plots:
ps.out.name <- "Binned_LUE_XY-ta05_Prediction_v24.ps"
postscript(
  file = paste(out.path,ps.out.name, sep=""), 
  width = 6, 
  height = 4, 
  title = ps.out.name,
  paper = "special", 
  horizontal = FALSE,
  onefile = TRUE
)
veg.bin.all <- seq(from=0, to=1, by=1)
tc.bin.all <- seq(from=0, to=5, by=1)
for (veg.bin in veg.bin.all){
  for (tc.bin in tc.bin.all){
    # Save subgroup of monthly LUE analyzer content:
    my.sub <- my.lue.mo.analyzer[which(
      my.lue.mo.analyzer$vbin == veg.bin & my.lue.mo.analyzer$tbin == tc.bin
      ),]
    #
    plot_xy(my.sub, method='a', basic=F, tau=T, alpha=T, k=0.5, pos='topleft', 
            stid=paste("veg.bin=",veg.bin," temp.bin=", tc.bin, sep=""))
    #
    if(FALSE){
      # Estimate a starting value of beta:
      beta.est <- estimate.beta(
        mean(my.sub$VPD.Pa),
        mean(my.sub$K.Pa),
        mean(my.sub$Gstar.Pa),
        mean(my.sub$eta.mPa_s)
      )
      #
      attempt <- try(
        nls.lm(
          par=list(b=beta.est),
          lower=c(0),
          upper=c(10000),
          fn = lueResid,
          my.obj = my.sub
        )
      )
      binned.lue.regression(my.lue.mo.analyzer,veg.bin,tc.bin,attempt)
    }
  }
}
dev.off()

# ************* (alpha) Tc *************
k.a <- 0.5
x.a <- f.of.alpha(my.lue.mo.analyzer$ALPHA, k.a) * my.lue.mo.analyzer$mgs.tc
tx.a <- expression(italic(f)(alpha)%*%Temperature~(degree*C))
y.a <- my.lue.mo.analyzer$LUEmo
ty.a <- expression(Monthly~LUE)
monthly.lue.regression(x.a, y.a, tx.a, ty.a)
# ************* VPD *************
k.b <- 0.15
x.b <- f.of.alpha(my.lue.mo.analyzer$ALPHA, k.b) * my.lue.mo.analyzer$mgs.vpd
tx.b <- expression(italic(f)(alpha)%*%Vapor~Pressure~Deficit~(kPa))
y.b <- my.lue.mo.analyzer$LUEmo
ty.b <- expression(Monthly~LUE)
monthly.lue.regression(x.b, y.b, tx.b, ty.b)
# ************* Elv *************
k.c <- 0.001
x.c <- f.of.alpha(my.lue.mo.analyzer$ALPHA, k.c) * my.lue.mo.analyzer$elv
tx.c <- expression(italic(f)*(alpha)%*%Elevation~(m))
y.c <- my.lue.mo.analyzer$LUEmo
ty.c <- expression(Monthly~LUE)
monthly.lue.regression(x.c, y.c, tx.c, ty.c)

# ~~~~~~~~~~~~ LUE AS PREDICTED BY ELV, TC, AND VPD ~~~~~~~~~~~~~
par(mfrow=c(1,3))
plot(
  my.lue.analyzer$mgs.tc, 
  my.lue.analyzer$lue,
  log(my.lue.analyzer$mgs.tc), 
  log(my.lue.analyzer$lue),
  type='p',
  cex.lab = 1.5,
  cex.axis = 1.5,
  xlab="log Temperature (C)",
  ylab="log LUE"
)
fit1 <- lm(
  my.lue.analyzer$lue ~ my.lue.analyzer$mgs.tc
)
fit1l <- lm(
  log(my.lue.analyzer$lue) ~ log(my.lue.analyzer$mgs.tc)
)

plot(
  my.lue.analyzer$mgs.vpd, 
  my.lue.analyzer$lue,
  log(my.lue.analyzer$mgs.vpd), 
  log(my.lue.analyzer$lue),
  type='p',
  cex.lab = 1.5,
  cex.axis = 1.5,
  xlab="VPD (kPa)",
  ylab="LUE"
)
fit2 <- lm(
  my.lue.analyzer$lue ~ my.lue.analyzer$mgs.vpd
)
fit2l <- lm(
  log(my.lue.analyzer$lue) ~ log(my.lue.analyzer$mgs.vpd)
)

plot(
  my.lue.analyzer$elv[which(my.lue.analyzer$elv > -100)], 
  my.lue.analyzer$lue[which(my.lue.analyzer$elv > -100)],
  log(my.lue.analyzer$elv[which(my.lue.analyzer$elv > 0)]), 
  log(my.lue.analyzer$lue[which(my.lue.analyzer$elv > 0)]),
  type='p',
  cex.lab = 1.5,
  cex.axis = 1.5,
  xlab="log Elevation (meters)",
  ylab="log LUE"
)
fit3 <- lm(
  my.lue.analyzer$lue[which(my.lue.analyzer$elv > -100)] ~ my.lue.analyzer$elv[which(my.lue.analyzer$elv > -100)]
  )
fit3l <- lm(
  log(my.lue.analyzer$lue[which(my.lue.analyzer$elv > 0)]) ~ log(my.lue.analyzer$elv[which(my.lue.analyzer$elv > 0)])
)
y3 <- my.lue.analyzer$lue[which(my.lue.analyzer$elv > -100)]
x3 <- 1/my.lue.analyzer$elv[which(my.lue.analyzer$elv > -100)]
fit3s <- lm(y3 ~ x3)

par(mfrow=c(1,1))
plot(
  1/sqrt(my.lue.analyzer$mgs.vpd), 
  my.lue.analyzer$lue,
  type='p',
#  cex.lab = 1.5,
#  cex.axis = 1.5,
  xlab="VPD (kPa)",
  ylab="LUE"
)
x <- 1/(my.lue.analyzer$mgs.vpd)
fit2s <- lm(
  my.lue.analyzer$lue ~ x
)
summary(fit2s)

x4 <- (1/sqrt(my.lue.analyzer$mgs.vpd) - my.lue.analyzer$mgs.tc)/(1/sqrt(my.lue.analyzer$mgs.vpd) + 2*my.lue.analyzer$mgs.tc)
fit4s <- lm(my.lue.analyzer$lue ~ x4)

y <- my.lue.analyzer$lue[which(my.lue.analyzer$elv > -100)]
x1 <- my.lue.analyzer$elv[which(my.lue.analyzer$elv > -100)]
x2 <- my.lue.analyzer$mgs.vpd[which(my.lue.analyzer$elv > -100)]
x3 <- my.lue.analyzer$mgs.tc[which(my.lue.analyzer$elv > -100)]
fit4 <- lm(y ~ x1 * x2 * x3)

x4 <- x1 * x2
x5 <- x1 * x3
x6 <- x1 * x2 * x3
fit5 <- lm(y ~ x1 + x4 + x5 + x6)
fit6 <- lm(y ~ x1 + x4 + x6)
fit7 <- lm(y ~ x4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ALL STATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Comparison of ALPHA & fAPAR
alfa <- cbind(content$ALPHA, content$fAPAR)
alfa <- as.data.frame(alfa)
names(alfa) <- c("ALPHA", "fAPAR")
boxplot(alfa)

# \\\\\\\\\\\\ PLOTTING ///////////// #
# Straight GPP v PPFD
my.x1 <- content$PPFD.mol_m2
x.label1 <- expression(paste("PPFD (", mol%.%m^{-2},")"))

# Basic LUE (GPP v fAPAR x PPFD)
my.x2 <- content$PPFD.mol_m2 * content$fAPAR
x.label2 <- expression(paste("fAPAR x PPFD (", mol%.%m^{-2},")"))

# (GPP v ALPHA x PPFD)
my.x3 <- content$PPFD.mol_m2 * content$ALPHA
x.label3 <- expression(paste(alpha," x PPFD (", mol%.%m^{-2},")"))

# (GPP v ALPHA x fAPAR x PPFD)
my.x4 <- content$PPFD.mol_m2 * content$ALPHA * content$fAPAR
x.label4 <- expression(paste(alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

# Calculate m-values:
my.m1 <- calc_m(content, b=0.5)
my.m2 <- calc_m(content, b=1)
my.m3 <- calc_m(content, b=5)
my.m4 <- calc_m(content, b=10)
my.m5 <- calc_m(content, b=20)
my.m7 <- calc_m(content, b=2.5)
my.m8 <- calc_m(content, b=1.25)
my.m9 <- calc_m(content, b=0.75)

my.m10 <- calc_m(content, b=10)
my.m15 <- calc_m(content, b=15)
my.m20 <- calc_m(content, b=20)
my.m25 <- calc_m(content, b=25)
my.m30 <- calc_m(content, b=30)

# See how b influences m:
mvals <- cbind(my.m10, my.m15, my.m20, my.m25, my.m30)
mvals <- as.data.frame(mvals)
names(mvals)<-c("10", "15", "20", "25", "30")
par(mar=c(4.75,4.75,1,1))
m.box <- boxplot(
  mvals,
  xlab=expression(sqrt(b/a)),
  ylab=expression(italic(m))
)
# Save the outliers from the boxplot:
mm1<-match(m.box$out[which(m.box$group == 1)],my.m10)
mm2<-match(m.box$out[which(m.box$group == 2)],my.m15)
mm3<-match(m.box$out[which(m.box$group == 3)],my.m20)
mm4<-match(m.box$out[which(m.box$group == 4)],my.m25)
mm5<-match(m.box$out[which(m.box$group == 5)],my.m30)
mmu<-unique(c(mm1,mm2,mm3,mm4,mm5))
# Check VPD of outliers:
content$VPD.Pa[mmu]
hist(content$VPD.Pa, 
     freq=T,
     main="",
     xlab="VPD (Pa)",
     col="grey"
)
m.box$stats[1,1]

which(my.m10 < m.box$stats[2,1])

plot(
  c(0.5, 0.75, 1, 1.25, 2.5, 5, 10, 20),
  c(mean(my.m1), mean(my.m9), mean(my.m2), mean(my.m8), mean(my.m7), 
    mean(my.m3), mean(my.m4), mean(my.m5)),
  type="b",
  xlab = expression(sqrt(b/a)),
  ylab = expression(mean(italic(m))),
  ylim = c(0,1)
)

#  GPP v m(assuming b=0.5) x PPFD
my.x5 <- content$PPFD.mol_m2 * my.m1
x.label5 <- expression(paste(italic(m)[b==0.5]," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1) x PPFD
my.x6 <- content$PPFD.mol_m2 * my.m2
x.label6 <- expression(paste(italic(m)[b==1]," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=5) x PPFD
my.x7 <- content$PPFD.mol_m2 * my.m3
x.label7 <- expression(paste(italic(m)[b==5]," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=10) x PPFD
my.x8 <- content$PPFD.mol_m2 * my.m4
x.label8 <- expression(paste(italic(m)[b==10]," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=20) x PPFD
my.x9 <- content$PPFD.mol_m2 * my.m5
x.label9 <- expression(paste(italic(m)[b==20]," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=0.5) x fAPAR x PPFD
my.x10 <- content$PPFD.mol_m2 * my.m1 * content$fAPAR
x.label10 <- expression(paste(italic(m)[b==0.5]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=0.5) x ALPHA x fAPAR x PPFD
my.x11 <- content$PPFD.mol_m2 * my.m1 * content$fAPAR * content$ALPHA
x.label11 <- expression(paste(italic(m)[b==0.5]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1) x fAPAR x PPFD
my.x12 <- content$PPFD.mol_m2 * my.m2 * content$fAPAR
x.label12 <- expression(paste(italic(m)[b==1]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1) x ALPHA x fAPAR x PPFD
my.x13 <- content$PPFD.mol_m2 * my.m2 * content$fAPAR * content$ALPHA
x.label13 <- expression(paste(italic(m)[b==1]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=5) x fAPAR x PPFD
my.x14 <- content$PPFD.mol_m2 * my.m3 * content$fAPAR
x.label14 <- expression(paste(italic(m)[b==5]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=5) x ALPHA x fAPAR x PPFD
my.x15 <- content$PPFD.mol_m2 * my.m3 * content$fAPAR * content$ALPHA
x.label15 <- expression(paste(italic(m)[b==5]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=10) x fAPAR x PPFD
my.x16 <- content$PPFD.mol_m2 * my.m4 * content$fAPAR
x.label16 <- expression(paste(italic(m)[b==10]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=10) x ALPHA x fAPAR x PPFD
my.x17 <- content$PPFD.mol_m2 * my.m4 * content$fAPAR * content$ALPHA
x.label17 <- expression(paste(italic(m)[b==10]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=20) x fAPAR x PPFD
my.x18 <- content$PPFD.mol_m2 * my.m5 * content$fAPAR
x.label18 <- expression(paste(italic(m)[b==20]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=20) x ALPHA x fAPAR x PPFD
my.x19 <- content$PPFD.mol_m2 * my.m5 * content$fAPAR * content$ALPHA
x.label19 <- expression(paste(italic(m)[b==20]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=0.5) x ALPHA x PPFD
my.x20 <- content$PPFD.mol_m2 * my.m1 * content$ALPHA
x.label20 <- expression(paste(italic(m)[b==0.5]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1) x ALPHA x PPFD
my.x21 <- content$PPFD.mol_m2 * my.m2 * content$ALPHA
x.label21 <- expression(paste(italic(m)[b==1]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=5) x ALPHA x PPFD
my.x22 <- content$PPFD.mol_m2 * my.m3 * content$ALPHA
x.label22 <- expression(paste(italic(m)[b==5]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=10) x ALPHA x PPFD
my.x23 <- content$PPFD.mol_m2 * my.m4 * content$ALPHA
x.label23 <- expression(paste(italic(m)[b==10]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=20) x ALPHA x PPFD
my.x24 <- content$PPFD.mol_m2 * my.m5 * content$ALPHA
x.label24 <- expression(paste(italic(m)[b==20]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

# GPP v m(estimate b pointwise) x PPFD
b.est.all <- estimate.sqba(content$VPD.Pa,content$K.Pa)
my.m6 <- calc_m(content, b=b.est.all)
my.x25 <- content$PPFD.mol_m2 * my.m6
x.label25 <- expression(paste(italic(m)[b==opt]," x PPFD (", mol%.%m^{-2},")"))

# GPP v m(estimate b pointwise) x fAPAR x PPFD
my.x26 <- content$PPFD.mol_m2 * my.m6 * content$fAPAR
x.label26 <- expression(paste(italic(m)[b==opt]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

# GPP v m(estimate b pointwise) x fAPAR x PPFD
my.x27 <- content$PPFD.mol_m2 * my.m6 * content$ALPHA
x.label27 <- expression(paste(italic(m)[b==opt]," x ",alpha," x PPFD (", mol%.%m^{-2},")"))

# GPP v m(estimate b pointwise) x fAPAR x PPFD
my.x28 <- content$PPFD.mol_m2 * my.m6 * content$ALPHA * content$fAPAR
x.label28 <- expression(paste(italic(m)[b==opt]," x ",alpha," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=2.5) x fAPAR x PPFD
my.x29 <- content$PPFD.mol_m2 * my.m7 * content$fAPAR
x.label29 <- expression(paste(italic(m)[b==2.5]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1.5) x fAPAR x PPFD
my.x30 <- content$PPFD.mol_m2 * my.m8 * content$fAPAR
x.label30 <- expression(paste(italic(m)[b==1.25]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

#  GPP v m(assuming b=1.5) x fAPAR x PPFD
my.x31 <- content$PPFD.mol_m2 * my.m9 * content$fAPAR
x.label31 <- expression(paste(italic(m)[b==0.75]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

my.x99 <- content$PPFD.mol_m2 * content$fAPAR * my.m10
x.label99 <- expression(paste(italic(m)[b==0]," x fAPAR x PPFD (", mol%.%m^{-2},")"))

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
# ~~~~~~~~~~ Select model you want to plot: ~~~~~~~~~~ #
# //////////////////////////////////////////////////// #
# the fAPAR x PPFD's : 2, 10, 12, 14, 16, 18, 29
my.x <- my.x99
x.label <- x.label99
plot_vpd(content, my.x, x.label)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #
# ~~~~~~~~~~ Plot beta limits ~~~~~~~~~~ #
# ////////////////////////////////////// #
# 0 < beta < 85000
#   0.0 < VPD.Pa < 2200  => 
#      0.0 < beta (k=5) < 9500
#      0.0 < beta (k=15) < 4000
#      0.0 < beta (k=25) < 2600
#      0.0 < beta (k=50) < 1400
#      0.0 < beta (k=100) < 680
#   1.0 < Gstar.Pa < 5.5; mean = 2.3 * large positive influence
#   5.0 < K.Pa < 100
#   0.5 < eta.mPa_s < 2.0; mean = 1.2  * large influence
my.params <- list(vpd=500, k=5, gs=5.5, n=2, pa=37)
my.betas <- 0*seq(101)
my.ms <- 0*seq(101)
for(i in seq(101)){
  my.betas[i]<- estimate.beta(
    vpd=(22*(i-1)),
    k=my.params$k,
    gs=my.params$gs,
    n=my.params$n
  )
  my.ms[i] <- (my.params$pa-my.params$gs)/(
    my.params$pa+
      2*my.params$gs+
      3*my.params$gs*sqrt(1.6*my.params$n*22*(i-1)/(my.betas[i]*(my.params$k+my.params$gs))))
}
par(mar=c(4,4.5,1,1))
plot(
  seq(from=0, to=2200, by=22), 
  my.betas,
  type='l',
  col='red',
  lwd=2,
  xlab=expression(VPD~(Pa)),
  ylab=expression(beta),
  ylim=c(0, 110000)
)
for(j in seq(95)){
  for(i in seq(101)){
    my.betas[i]<- estimate.beta(
      vpd=(22*(i-1)),
      k=(5+j),
      gs=my.params$gs,
      n=my.params$n
    )
  }
  lines(seq(from=0,to=2200,by=22),my.betas,lwd=2,col='darkgray')
}
lines(seq(from=0,to=2200,by=22),my.betas,lwd=2,col='black')


my.betas <- 0*seq(96)
my.ms <- 0*seq(96)
for(i in seq(96)){
  my.betas[i]<- estimate.beta(
    vpd=22,
    k=(4+i),
    gs=my.params$gs,
    n=my.params$n
  )
  my.ms[i] <- (my.params$pa-my.params$gs)/(
    my.params$pa+
      2*my.params$gs+
      3*my.params$gs*sqrt(1.6*my.params$n*22*(i-1)/(my.betas[i]*(my.params$k+my.params$gs))))
}
plot(
  seq(from=5, to=100, by=1), 
  my.betas,
  type='l',
  col='black',
  lwd=2,
  xlab=expression(K~(Pa)),
  ylab=expression(beta),
  ylim=c(0, 10000)
)
for(j in seq(2,100,2)){
  for(i in seq(96)){
    my.betas[i]<- estimate.beta(
      vpd=(22*(j)),
      k=(4+i),
      gs=my.params$gs,
      n=my.params$n
    )
  }
  lines(seq(from=5,to=100,by=1),my.betas,lwd=2,col='darkgray')
}
lines(seq(from=5,to=100,by=1),my.betas,lwd=2,col='red')
