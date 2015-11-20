## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
##  Test: Calculate GPP for monthly input data for CH-Oe1
##  year: 2000, elv: 450, lon: 7.73, lat: 47.3
##  using standard input data, read from files
## ////////////////////////////////////////////////////////////////////
source('../main_r/pmodel.R')

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
nmonth <- 12

##---------------------------------------------------------------------
## STANDARD PARAMETERS
##---------------------------------------------------------------------
elv     <- 450.0                          # m.a.s.l.
patm    <- calc_patm(elv)                 # using function defined in pmodel.R
co2     <- 376.0                          # (ppm)
ca      <- co2_to_ca( co2, patm )         # atms. CO2 (Pa); using function defined in pmodel.R; 
fpar    <- 1.0                            # fraction of absorbed photosynthetically active radiation (unitless)
cpalpha <- 1.26                           # Cramer-Prentice alpha (unitless)
beta    <- 244.033                        # Unit cost of carboxylation
kPo     <- 101325.0                       # standard atmosphere, Pa (Allen, 1973)
kTo     <- 25.0                           # base temperature, deg C (Prentice, unpublished)

##---------------------------------------------------------------------
## GET MONTHLY CLIMATE
##---------------------------------------------------------------------
## temperature: read monthly file (deg C)
filn_temp <- "../data/mtemp_CH-Oe1_2002.txt"
df.mtemp <- read.csv( filn_temp )
mtemp <- df.mtemp$temp

## precipitation, convert from daily to monthly values (mm/month)
filn_prec <- "../data/mprec_CH-Oe1_2002.txt"
df.mprec <- read.csv( filn_prec )
mprec <- df.mprec$prec

## sunshine fraction (unitless)
filn_fsun <- "../data/mfsun_CH-Oe1_2002.txt"
df.mfsun <- read.csv(filn_fsun)
mfsun <- df.mfsun$fsun

## vapour pressure (mvapr, hPa)
filn_vapr <- "../data/mvapr_CH-Oe1_2002.txt"
df.mvapr <- read.csv(filn_vapr)
mvapr <- df.mvapr$vapr

## Calculate vapour pressure deficit from vapour pressure and temperature (mvpd, kPa)
mvpd <- rep(NA,nmonth)
for (moy in 1:nmonth){
  mvpd[moy] <- calc_vpd( mtemp[moy], mvapr[moy] )
}

## PPFD from file produced as output from STASH (sofun implementation, '*.m.qm.out')
filn_ppfd <- "../data/CH-Oe1_2002.m.qm.out"
df.ppfd <- read.table( filn_ppfd, col.names=c("year","ppfd") )
istart <- which.min( abs(df.ppfd$year-2002.0) )  
df.ppfd <- df.ppfd[ istart:(istart+nmonth-1), ]  ## take subset of year 2002


##------------------------------------------------------------
## RUN P-MODEL
## ... and several "components" (K, viscosity, Gamma-star, 
## Chi...) using functions defined within pmodel.R
## with monthly input for PPFD, air temperature, and VPD
##------------------------------------------------------------
mgpp <- rep( NA, nmonth )
mrd  <- rep( NA, nmonth )
mvcmax    <- rep( NA, nmonth )
mnrubisco <- rep( NA, nmonth )
mluenet   <- rep( NA, nmonth )
for (moy in 1:nmonth){
  out <- pmodel( fpar=fpar, ppfd=df.ppfd$ppfd[moy], co2=co2, tc=mtemp[moy], cpalpha=cpalpha, vpd=mvpd[moy], elv=elv )
  mgpp[moy]      <- out$gpp
  mrd[moy]       <- out$rd
  mvcmax[moy]    <- out$vcmax
  mnrubisco[moy] <- out$actnv
  mluenet[moy]   <- out$lue
}

## CALCULATE K
## Michaelis-Menten coefficient (Pa)
kmm       <- sapply( mtemp, FUN = function(x) calc_k(x, patm) )
kmm_colin <- sapply( mtemp, FUN = calc_k_colin )

## CALCULATE VISCOSITY
visc25    <- viscosity_h2o( kTo, kPo )
visc      <- mapply( viscosity_h2o, mtemp, patm )
visc_star <- visc/visc25
visc_vogel<- sapply( mtemp, FUN = viscosity_h2o_vogel )
visc_out  <- visc * 1e3  # convert from Pa s to mPa s
visc_vogel_out <- visc_vogel * 1e3  # convert from Pa s to mPa s

## CALCULATE GAMMA-STAR
gstar       <- sapply( mtemp, FUN = calc_gstar_gepisat )
gstar_colin <- sapply( mtemp, FUN = calc_gstar_colin )

## CALCULATE CHI 
chi_wh     <- rep( NA, nmonth )
chi_simpl  <- rep( NA, nmonth )
chi_full   <- rep( NA, nmonth )
for (moy in 1:nmonth){
  chi_wh[moy]    <- lue_approx( mtemp[moy], mvpd[moy], elv, ca, gstar )$chi
  chi_simpl[moy] <- lue_vpd_simpl( kmm[moy], gstar[moy], visc_star[moy], ca, mvpd[moy], beta )$chi
  chi_full[moy]  <- lue_vpd_full(  kmm[moy], gstar[moy], visc_star[moy], ca, mvpd[moy], beta )$chi
}

##------------------------------------------------------------
## PRINT TO SCREEN
##------------------------------------------------------------
# zzz <- file("check_pmodel_r_version_beni.txt","w")
cat( "------------------------------------------------------------", "\n"   ) #, file=zzz )
cat( "PMODEL, R-VERSION", "\n"   ) #, file=zzz )
cat( "------------------------------------------------------------", "\n"   ) #, file=zzz )
cat( "Using as input/parameters:", "\n"  ) #, file=zzz )
cat( "Elevation (m.asl.):               ", elv, "\n"   ) #, file=zzz )
cat( "Atmospheric CO2 (ppm):            ", co2, "\n"   ) #, file=zzz )
cat( "fAPAR (unitless):                 ", fpar, "\n"   ) #, file=zzz )
cat( "Cramer-Prentice alpha:            ", cpalpha, "\n"   ) #, file=zzz )
cat( "beta (unit cost of carboxylation):", beta, "\n"   ) #, file=zzz )
cat( "Standard atmosphere (Pa):         ", kPo, "\n"  ) #, file=zzz )
cat( "Standard temperature (deg C):     ", kTo, "\n"  ) #, file=zzz )
cat( "------------------------------------------------------------", "\n"   ) #, file=zzz )
cat( "This is converted by P-model to:  ", "\n"  ) #, file=zzz )
cat( "Atmospheric pressure (Pa):        ", patm, "\n"  ) #, file=zzz )
cat( "Ambient CO2 partial pressure (Pa):", ca, "\n"  ) #, file=zzz )
cat( "------------------------------------------------------------", "\n"  ) #, file=zzz )
cat( "Input files:", "\n"   ) #, file=zzz )
cat( "temperature (given in deg C):                      ", "\n" ) #, file=zzz )
cat( "   ", filn_temp, "\n" ) #, file=zzz )
cat( "precipitation (given in mm/month):                      ", "\n" ) #, file=zzz )
cat( "   ", filn_prec, "\n" ) #, file=zzz )
cat( "sunshine fraction (given in unitless):                      ", "\n" ) #, file=zzz )
cat( "   ", filn_fsun, "\n" ) #, file=zzz )
cat( "vapour pressure (given in hPa):                      ", "\n" ) #, file=zzz )
cat( "   ", filn_vapr, "\n" ) #, file=zzz )
cat( "PPFD (given in mol/m2):                           ", "\n" ) #, file=zzz )
cat( "   ", filn_ppfd, "\n" ) #, file=zzz )
cat( "------------------------------------------------------------", "\n"  ) #, file=zzz )
cat( "Outputs:", "\n"   ) #, file=zzz )
cat( "Michaelis-Menten K (Pa):                      ", "\n" ) #, file=zzz )
cat( "   ", kmm, "\n" ) #, file=zzz )
cat( "Michaelis-Menten K, using 'calc_k_colin' (Pa):", "\n" ) #, file=zzz )
cat( "   ", kmm_colin, "\n" ) #, file=zzz )
cat( "Viscosity (mPa s):", "\n" ) #, file=zzz )
cat( "   ", visc_out, "\n" ) #, file=zzz )
cat( "Gamma-star, using 'calc_gstar_gepisat' (Pa):", "\n" ) #, file=zzz )
cat( "   ", gstar, "\n" ) #, file=zzz )
cat( "Gamma-star, using 'calc_gstar_colin' (Pa):", "\n" ) #, file=zzz )
cat( "   ", gstar_colin, "\n" ) #, file=zzz )
cat( "Chi, using full method 'lue_vpd_full' (unitless):", "\n" ) #, file=zzz )
cat( "   ", chi_full, "\n" ) #, file=zzz )
cat( "Chi, using simplified method 'lue_vpd_simpl' (unitless):", "\n" ) #, file=zzz )
cat( "   ", chi_simpl, "\n" ) #, file=zzz )
cat( "Chi, using Wang-Han method 'lue_approx' (unitless):", "\n" ) #, file=zzz )
cat( "   ", chi_wh, "\n" ) #, file=zzz )
cat( "GPP (mol C m-2 s-1):", "\n" ) #, file=zzz )
cat( "   ", mgpp, "\n" ) #, file=zzz )
cat( "Rd (mol C m-2 s-1):", "\n" ) #, file=zzz )
cat( "   ", mrd, "\n" ) #, file=zzz )
cat( "Vcmax (mol C m-2 s-1):", "\n" ) #, file=zzz )
cat( "   ", mvcmax, "\n" ) #, file=zzz )
cat( "------------------------------------------------------------", "\n" ) #, file=zzz  )
# close(zzz)


##---------------------------------------------------------------------
## USE DATA FOR DIRECT COMPARISON WITH DATA PROVIDED BY WANG HAN AND 
## TYLER AND FROM OBSERVATIONS
##---------------------------------------------------------------------
## Wang Han data
## read input data for GPP calculation
indata <- read.csv( "/alphadata01/bstocker/sofun/trunk/input_raw/Wang_Han_beta_estimates.csv", header=TRUE )

## get day of year (doy) and year
indata$doy  <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$yday+1
indata$year <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$year+1900
indata$moy  <- as.POSIXlt( indata$Timestamp, format="%d/%m/%Y" )$mon+1

## take sub-set of this for station CH-Oe1
indata <- indata[ indata$station == "CH-Oe1", ]
indata <- indata[ indata$year == 2002, ]

## Tyler data
mgpp_tyler <- read.csv( "/alphadata01/bstocker/sofun/trunk/components/mgpp_gepisat_tyler.csv", header=TRUE )

## GPP data from observations
filnam <- "/alphadata01/bstocker/sofun/trunk/input/CH-Oe1_daily_gpp_med_STANDARD.txt"
df.gpp <- read.table( filnam, header=FALSE, col.names="gpp" )
source("/alphadata01/bstocker/utilities/daily2monthly.R") ## this is in our 'utilities' repository
mgpp_obs <- daily2monthly( df.gpp$gpp, method="sum" )


## SAVE DATA FOR USE BY TYLER AND WANG-HAN
benioutput <- data.frame( year=rep(2002,nmonth), moy=1:nmonth, kpp_Pa=kmm
  , visc_Pa_s=visc, gstar_Pa=gstar, chi_WH_method=chi_wh, chi_simpl=chi_simpl, chi_full=chi_full 
  )
save( benioutput, file="benioutput.Rdata" )
write.csv( benioutput, file="benioutput.txt", row.names=FALSE )

# COMPARE DATA
# temperature: ok
plot( 1:nmonth, mtemp, type="l", col="red" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$Tair_degC[idx] )
}

## VPD xxx not ok xxx: is WH's data really in kPa (not Pa)?
# pdf("vpd_comparison.pdf")
plot( 1:nmonth, mvpd, type="l", col="red", xlab="MOY"  )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$D_kPa[idx] )
}
# dev.off()

## K - Michaelis-Menten coefficient: ok
# pdf("Kc_comparison.pdf")
plot( 1:nmonth, kmm, type="l", col="red", xlab="MOY"  )
lines(1:nmonth, kmm_colin, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$K_Pa[idx] )
}
legend( "topleft", c("GePiSaT method", "Colin method"), lty=1, bty="n", col=c("red","green") )
# dev.off()

## Viscosity: systematically lower values by my code, but tested to be identical with python implementation in gepisat
# pdf("viscosity_comparison.pdf")
plot( 1:nmonth, visc_out, type="l", col="red", xlab="MOY", ylab="viscosity [mPa s]" )
lines(1:nmonth, visc_vogel_out, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$ns[idx] )
}
legend( "bottomleft", c("Huber method", "Vogel method"), lty=1, bty="n", col=c("red","green") )
# dev.off()

## Gamma-star: ok
# pdf("gammastar_comparison.pdf")
plot( 1:nmonth, gstar, type="l", col="red", xlab="MOY"  )
lines(1:nmonth, gstar_colin, col="green" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$Gs_Pa[idx] )
}
legend( "topleft", c("GePiSaT method", "Colin method"), lty=1, bty="n", col=c("red","green") )
# dev.off()

## CHI - WH method
# pdf("chi_comparison.pdf")
plot( 1:nmonth, chi_wh, type="l", col="red", xlab="MOY"  )
lines( 1:nmonth, chi_full, col="blue" )
for (idx in 1:dim(indata)[1]){
  points( indata$moy[idx], indata$chi_wang_han[idx] )
  points( indata$moy[idx], indata$chi_vpd[idx], col="green" )
}
legend( "topleft", c("Wang-Han method", "theoretical full method"), lty=1, bty="n", col=c("red","blue") )
# dev.off()

# ## GPP 
# pdf("gpp_comparison.pdf")
plot(  1:nmonth, mgpp, type="l", col="red", xlab="MOY", ylab="GPP (mol C m-2 month-1)"  )
lines( 1:nmonth, mgpp-mrd, lty=2, col="red", xlab="MOY"  )
lines( 1:nmonth, mgpp_tyler$gpp, col="blue")
lines( 1:nmonth, mluenet*df.ppfd$ppfd, col="magenta")
points( 1:nmonth, mgpp_obs )
legend( "topleft", c("simulted GPP", "simulated GPP, GePiSaT", "simulated GPP-Rd"), lty=1, bty="n", col=c("red","blue","magenta"))
legend( "topleft", c("","","","FLUXNET data"), bty="n", lty=0, pch=c(NA,NA,NA,1) )
title("CH-Oe1, fAPAR=1, year 2002")
# dev.off()

# ## VCMAX 
# # pdf("vcmax_comparison.pdf")
# plot( 1:nmonth, mvcmax, type="l", col="red", xlab="MOY"  )
# # dev.off()

# ## LEAF N 
# # pdf("nrubisco_comparison.pdf")
# plot( 1:nmonth, mnrubisco, type="l", col="red", xlab="MOY"  )
# # dev.off()

