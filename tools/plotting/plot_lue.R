# RStudio v. 0.97.551
#
# plot_lue.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-01-09 -- created
# 2014-05-08 -- last updated
# * based in-part on query_plots.R, created 2013-07-12
#
# ------------
# description:
# ------------
# This script processes the LUE data from monthly LUE files (output from 
# model.py) and produces plots.
#
# To convert the multi-page PS file to PDF:
# 1. convert PS to PDF:
#    ps2pdf -dPDFSETTINGS=/prepress ALL_LUE.ps
# 2. crop whitespace from PDF:
#    pdfcrop --margins 25 ALL_LUE.pdf
#
# ----------
# changelog:
# ----------
# 01. copied relevant info from query_plots.R [14.01.09]
# 02. created an outlier_list for storing station and dates of outliers found
# 03. updated for the latest GePiSaT_All_Data.txt files [14.05.08]
# 04. improved plotting parameters [14.05.08]
#
# --------
# sources:
# --------
# 01. peirce_dev.R (if you are identifying outliers)
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////
filter_contents <- function(my.obj){
  # Filter out negative temperature and fAPAR:
  my.filt <- matrix(ncol=(dim(my.obj)[2]),nrow=1)
  my.filt <- as.data.frame(my.filt)
  names(my.filt) <- names(my.obj)
  my.filt <- my.obj[
    which(my.obj$fAPAR >= 0 & my.obj$Tc.deg_C >= 0 & my.obj$VPD.kPa > 0),
    ]
  my.filt
}

find_outliers <- function(model.fit, N){
  # Returns a list of indexes that breach Peirce's deviation
  # model.fit :: R lm-object
  # N :: number of observations
  n_index <- integer(0)
  if(N>2){
    sq.err <- as.numeric(model.fit$residuals)^2       # squared-error
    mse <- sum(sq.err)/(N-1)                          # mean-squared-error
    x2 <- peirce_dev(N, 1, 1)                         # Peirce's deviation
    d2 <- mse*x2                                      # threshold squared-error
    n_index <- which(sq.err >= d2)                    # list of outlier indices
    n_found <- length(n_index)                        # number of outliers found
    if (n_found == 0){
      x2 = peirce_dev(N, 2, 1)
      d2 = mse*x2
      n_index = which(sq.err >= d2)
      n_found = length(n_index)
    }
    n = 1
    while(n <= n_found){
      n = n+1
      x2 = peirce_dev(N,n,1)
      d2 = mse*x2
      n_index = which(sq.err >= d2)
      n_found = length(n_index)
    }
  }
  n_index
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FIND LUE FILES ###########################################################
# /////////////////////////////////////////////////////////////////////////////
lue.file.path=paste(
  "/home/user/Dropbox/",
  "Work/Imperial/flux/results/2002-06/lue/",
  sep="")
lue.files = list.files(path = lue.file.path, pattern = "^GePiSaT.*txt$");

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### PROCESS LUE DATA #########################################################
# /////////////////////////////////////////////////////////////////////////////
# Get LUE data:
lue.all.file <- lue.files[9]
content <- read.csv(paste(lue.file.path, lue.all.file, sep=""), header=T)
content.bak <- content

# Get the station names:
all.stations <- levels(unique(content$Station))

#ps.out.name <- "ALL_LUE_out-no-legend.ps"
ps.out.name <- "Basic_LUE_labeled_v24.ps"
#
# Create a blank character array for storing outlier information
#outlier_list = character()
#
postscript(
  file = paste(lue.file.path,ps.out.name, sep=""), 
  title = ps.out.name,
  width = 4, 
  height = 4, 
  paper = "special", 
  horizontal = FALSE,
  onefile = TRUE
)

for (stationid in all.stations){
  # Reset content:
  content <- content.bak
  #
  # Get indexes for current station id and filter content for current station
  station.idx <- which(content$Station == stationid)
  content <- content[station.idx,]
  content <- filter_contents(content)
  #
  # Create new column for fitting parameter (i.e., PPFD*fAPAR):
  content$PARAM <- content$fAPAR * content$PPFD.mol_m2
  #
  # Fit linear model:
  fit <- lm(content$GPP.mol_m2 ~ 0 + content$PARAM)
  #hist(fit$residuals)
  #
  # Save fitting parameters:
  n <- length(content$GPP.mol_m2)
  fit.coef <- summary(fit)$coefficients[1,1]
  #
  # <<<< NEW R2 CALC >>>>
  #mse <- mean(as.numeric(fit$residuals)^2)
  #sxy <- sum(
  #  (content$GPP.mol_m2 - mean(content$GPP.mol_m2))*
  #    (content$PARAM - mean(content$PARAM))
  #)
  #sxx <- sum((content$GPP.mol_m2 - mean(content$GPP.mol_m2))^2)
  #syy <- sum((content$PARAM - mean(content$PARAM))^2)
  #r <- sxy / sqrt(sxx) / sqrt(syy)
  #r2 <- r^2
  #
  r2 <- summary(fit)$adj.r.squared
  my.p <- summary(fit)$coefficients[1,4]
  #
  # <----- Find outliers ----->
  #my_outliers <- find_outliers(fit, n)
  #
  # Plot parameters:
  rp = vector('expression',4)
  rp[1] = substitute(
    expression(LUE == MYLUE),
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
  #
  fit.coef.str <- sprintf("%0.4f", fit.coef)
  main.text <- stationid
  #
  # <----- Save outlier meta to list: ----->
  #outlier_finds <- paste(
  #  rep(stationid,length(my_outliers)),
  #  levels(content$Timestamp)[my_outliers],
  #  sep=" "
  #  )
  #outlier_list <- c(outlier_list, outlier_finds)
  #
  # Calculate plotting extremes:
  max.x <- (
    floor(max(content$PARAM)+100) - 
      (floor(max(content$PARAM)+100) %% 100)
  )
  if (max.x == 900){
    max.x <- 1000
  }
  max.y <- (
    floor(max(content$GPP.mol_m2)+10) - 
      (floor(max(content$GPP.mol_m2)+10) %% 10)
  )
  # ------------------------
  # Save plot as postscript:
  # ------------------------
  par(mar=c(5.5,4.5,1,1))
  plot(
    content$PARAM, 
    content$GPP.mol_m2, 
    sub = NULL, 
    main = NULL,
    xlab = NA,
    ylab = NA,
    axes = F,
    xlim=c(0,max.x), 
    ylim=c(0,max.y)
  )
  box(lwd=2)
  #
  # Add Y axis ticks and labels:
  axis(side=2, las=1, tck=-0.02, labels=NA)
  axis(side=2, lwd=0, line=-0.4, las=1)
  #
  # Add X axis ticks and labels:
  axis(side=1, las=1, tck=-0.02, labels=NA)
  axis(side=1, lwd=0, line=-0.4)
  #
  mtext(
    side=2, 
    expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")), 
    line=2
  )
  mtext(
    side=1, 
    expression(paste("fPAR x PPFD (", mol%.%m^{-2},")")), 
    line=2
  )
  mtext(side=1, main.text, line=3)
  #
  legend('topleft', legend = rp, bty = 'n')
  abline(a = 0, b = fit.coef, lty=2, lwd=2)
  # <----- Highlight outliers: ----->
  #points(
  #  content$PARAM[my_outliers],
  #  content$GPP[my_outliers],
  #  col="red",
  #  pch=20
  #)
  #
  # -----------------
  # Save plot as bmp: 
  # -----------------
  #bmp(
  #  filename = paste(lue.file.path,bmp.out.name, sep=""), 
  #  width = 2400, 
  #  height = 2400, 
  #  units = "px",
  #  bg = "white", 
  #  res = 800
  #  )
  #par(mar=c(1,1,1,1))
  #plot(
  #  content$PARAM, 
  #  content$GPP,  
  #  main = NULL,
  #  xlab = "",
  #  ylab = "",
  #  axes = FALSE,
  #  frame.plot = TRUE,
  #  xlim=c(0,max(content$PARAM)), 
  #  ylim=c(0,max(content$GPP))
  #  )
  #box(lwd=2)
  #Axis(side=1, labels=FALSE, tck = -0.02)
  #Axis(side=2, labels=FALSE, tck = -0.02)
  #legend('topleft', legend = rp, bty = 'n')
  #abline(a = 0, b = fit.coef, lty=2, lwd=2)
  #dev.off()
  #
  # <----- reset variables ----->
  #rm(content,fit,my_outliers,outlier_finds)
  rm(fit, fit.coef, fit.coef.str, main.text, max.x, max.y, my.p, r2, rp, n, station.idx)
}
dev.off()

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### UNIVERSAL LUE ############################################################
# /////////////////////////////////////////////////////////////////////////////
# Find GePiSaT All Data file:
u.lue.path <- paste("/Users/twdavis",
                    "/Dropbox/Work/Imperial/flux/data/model_results/2002-06/",
                    sep=""
)
u.lue.file <- list.files(path = u.lue.path, pattern="*-17.txt")
#
# Read data & create fitting parameter column:
u.lue.data <- read.csv(paste(u.lue.path, u.lue.file, sep=""), header=T)
u.lue.data$PARAM <- u.lue.data$fAPAR * u.lue.data$PPFD
#
# Linear model fit:
u.lue.fit <- lm(
  u.lue.data$GPP[which(u.lue.data$PARAM > 0)] ~ 
    0 + u.lue.data$PARAM[which(u.lue.data$PARAM > 0)]
  )
n <- length(u.lue.data$GPP[which(u.lue.data$PARAM > 0)])
fit.coef <- summary(u.lue.fit)$coefficients[1,1]
p <- summary(u.lue.fit)$coefficients[1,4]
#
# Regression goodness-of-fit:
sxy <- sum(
  (
    u.lue.data$GPP[which(u.lue.data$PARAM > 0)] - 
      mean(u.lue.data$GPP[which(u.lue.data$PARAM > 0)])
    )*(
      u.lue.data$PARAM[which(u.lue.data$PARAM > 0)] - 
        mean(u.lue.data$PARAM[which(u.lue.data$PARAM >0)])
      )
)
sxx <- sum(
  (u.lue.data$GPP[which(u.lue.data$PARAM > 0)] - 
     mean(u.lue.data$GPP[which(u.lue.data$PARAM > 0)]))^2
)
syy <- sum(
  (u.lue.data$PARAM[which(u.lue.data$PARAM > 0)] - 
     mean(u.lue.data$PARAM[which(u.lue.data$PARAM > 0)]))^2
)
r <- sxy / sqrt(sxx) / sqrt(syy)
r2 <- r^2
#
# Plot parameters:
rp = vector('expression',4)
rp[1] = substitute(
  expression(LUE == MYLUE),
  list(MYLUE = format(fit.coef,dig=3)))[2]
rp[2] = substitute(
  expression(italic(n) == MYN),
  list(MYN = format(n,dig=2)))[2]
rp[3] = substitute(
  expression(italic(r)^2 == MYVALUE), 
  list(MYVALUE = format(r2,dig=3)))[2]
rp[4] = substitute(
  expression(italic(p) == MYOTHERVALUE), 
  list(MYOTHERVALUE = format(p, digits = 2)))[2]
#
# Plot:
par(mar=c(5.5,4.75,1,1))
plot(
  u.lue.data$PARAM, 
  u.lue.data$GPP, 
  sub = "Universal LUE", 
  main = NULL,
  xlab=expression(paste("fPAR x PPFD (", mol%.%m^{-2},")")), 
  ylab=expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")), 
  xlim=c(0,max(u.lue.data$PARAM)), 
  ylim=c(0,max(u.lue.data$GPP))
)
box(lwd=2)
legend('topleft', legend = rp, bty = 'n')
abline(a = 0, b = fit.coef, lty=2, lwd=2, col="red")
#
# Color-code points based on VPD:
vpd.list <- u.lue.data$VPD[which(u.lue.data$PARAM > 0)]
vpd.max <- max(vpd.list)
vpd.min <- min(vpd.list)
par(mar=c(4.75,4.75,1,1))
vpd.hist <- hist(
  vpd.list, 
  breaks=12, 
  freq=FALSE, 
  main="",
  xlab="VPD", 
  col="lightgreen",
  xlim=c(-0.5,2.0),
  ylim=c(0,2.5)
  )
vpd.blue <- which(u.lue.data$PARAM > 0 
                  & u.lue.data$VPD < vpd.hist$breaks[2] 
                  & u.lue.data$VPD >= vpd.hist$breaks[1])
vpd.lightblue <- which(u.lue.data$PARAM > 0
                       & u.lue.data$VPD < vpd.hist$breaks[3]
                       & u.lue.data$VPD >= vpd.hist$breaks[2])
vpd.cyan <- which(u.lue.data$PARAM > 0
                       & u.lue.data$VPD < vpd.hist$breaks[4]
                       & u.lue.data$VPD >= vpd.hist$breaks[3])
vpd.lightcyan <- which(u.lue.data$PARAM > 0
                       & u.lue.data$VPD < vpd.hist$breaks[5]
                       & u.lue.data$VPD >= vpd.hist$breaks[4])
vpd.turquoise <- which(u.lue.data$PARAM > 0
                        & u.lue.data$VPD < vpd.hist$breaks[6]
                        & u.lue.data$VPD >= vpd.hist$breaks[5])
vpd.lightgreen <- which(u.lue.data$PARAM > 0
                   & u.lue.data$VPD < vpd.hist$breaks[7]
                   & u.lue.data$VPD >= vpd.hist$breaks[6])
vpd.green <- which(u.lue.data$PARAM > 0
                    & u.lue.data$VPD < vpd.hist$breaks[8]
                    & u.lue.data$VPD >= vpd.hist$breaks[7])
vpd.yellow <- which(u.lue.data$PARAM > 0
                    & u.lue.data$VPD < vpd.hist$breaks[9]
                    & u.lue.data$VPD >= vpd.hist$breaks[8])
vpd.orange <- which(u.lue.data$PARAM > 0
                    & u.lue.data$VPD < vpd.hist$breaks[10]
                    & u.lue.data$VPD >= vpd.hist$breaks[9])
vpd.red <- which(u.lue.data$PARAM > 0
                 & u.lue.data$VPD < vpd.hist$breaks[11]
                 & u.lue.data$VPD >= vpd.hist$breaks[10])
vpd.pink <- which(u.lue.data$PARAM > 0
                 & u.lue.data$VPD < vpd.hist$breaks[12]
                 & u.lue.data$VPD >= vpd.hist$breaks[11])
vpd.magenta <- which(u.lue.data$PARAM > 0
                     & u.lue.data$VPD < vpd.hist$breaks[13]
                     & u.lue.data$VPD >= vpd.hist$breaks[12])
# Cyan -- 0.0<VPD<0.2:
points(u.lue.data$PARAM[vpd.cyan], u.lue.data$GPP[vpd.cyan], pch=20, col="cyan")
# Lightcyan -- 0.2<VPD<0.4:
points(u.lue.data$PARAM[vpd.lightcyan], u.lue.data$GPP[vpd.lightcyan], pch=20, col="lightcyan")
# Turquoise -- 0.4<VPD<0.6:
points(u.lue.data$PARAM[vpd.turquoise], u.lue.data$GPP[vpd.turquoise], pch=20, col="turquoise")
# Lightgreen -- 0.6<VPD<0.8:
points(u.lue.data$PARAM[vpd.lightgreen], u.lue.data$GPP[vpd.lightgreen], pch=20, col="lightgreen")
# Green -- 0.8<VPD<1.0:
points(u.lue.data$PARAM[vpd.green], u.lue.data$GPP[vpd.green], pch=20, col="green")
# Lightblue -- -0.2<VPD<0.0:
points(u.lue.data$PARAM[vpd.lightblue], u.lue.data$GPP[vpd.lightblue], pch=20, col="lightblue")
# Yellow -- 1.0<VPD<1.2:
points(u.lue.data$PARAM[vpd.yellow], u.lue.data$GPP[vpd.yellow], pch=20, col="yellow")
# Blue -- -0.4<VPD<-0.2:
points(u.lue.data$PARAM[vpd.blue], u.lue.data$GPP[vpd.blue], pch=20, col="blue")
# Orange -- 1.2<VPD<1.4:
points(u.lue.data$PARAM[vpd.orange], u.lue.data$GPP[vpd.orange], pch=20, col="orange")
# Pink:
points(u.lue.data$PARAM[vpd.pink], u.lue.data$GPP[vpd.pink], pch=20, col="pink")
# Red -- 1.6<VPD<1.8:
points(u.lue.data$PARAM[vpd.red], u.lue.data$GPP[vpd.red], pch=20, col="red")
# Magenta -- 1.8<VPD<2.0
points(u.lue.data$PARAM[vpd.magenta], u.lue.data$GPP[vpd.magenta], pch=20, col="magenta")
#
# ~~~~~~ Look at individual VPD bins ~~~~~~ #
my.points <- vpd.lightcyan
my.info <- c()
my.info$sub <- "LUE: 0.2<VPD<0.4"
my.info$col <- "lightcyan"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
my.fit <- lm(
  u.lue.data$GPP[my.points] ~ 0 + u.lue.data$PARAM[my.points])
my.n <- length(u.lue.data$GPP[my.points])
my.coef <- summary(my.fit)$coefficients[1,1]
my.p <- summary(my.fit)$coefficients[1,4]
#
#
my.sxy <- sum(
  (u.lue.data$GPP[my.points] - mean(u.lue.data$GPP[my.points]))*(
    u.lue.data$PARAM[my.points] - mean(u.lue.data$PARAM[my.points]))
)
my.sxx <- sum(
  (u.lue.data$GPP[my.points] - mean(u.lue.data$GPP[my.points]))^2
)
my.syy <- sum(
  (u.lue.data$PARAM[my.points] - mean(u.lue.data$PARAM[my.points]))^2
)
my.r <- my.sxy / sqrt(my.sxx) / sqrt(my.syy)
my.r2 <- my.r^2
#
#
my.rp = vector('expression',4)
my.rp[1] = substitute(
  expression(LUE == MYLUE),
  list(MYLUE = format(my.coef,dig=3)))[2]
my.rp[2] = substitute(
  expression(italic(n) == MYN),
  list(MYN = format(my.n,dig=2)))[2]
my.rp[3] = substitute(
  expression(italic(r)^2 == MYVALUE), 
  list(MYVALUE = format(my.r2,dig=3)))[2]
my.rp[4] = substitute(
  expression(italic(p) == MYOTHERVALUE), 
  list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
#
#
par(mar=c(5.5,4.75,1,1))
plot(
  u.lue.data$PARAM[my.points], 
  u.lue.data$GPP[my.points], 
  sub = my.info$sub, 
  main = NULL,
  xlab=expression(paste("fPAR x PPFD (", mol%.%m^{-2},")")), 
  ylab=expression(paste("GPP (", mol, " ", CO[2]%.%m^{-2},")")), 
  xlim=c(0,max(u.lue.data$PARAM)), 
  ylim=c(0,max(u.lue.data$GPP)),
  pch = 20,
  col = my.info$col
)
box(lwd=2)
legend('topleft', legend = my.rp, bty = 'n')
points(
  u.lue.data$PARAM[my.points],
  u.lue.data$GPP[my.points],
  pch=1,
  col="grey"
  )
abline(a = 0,b = my.coef, lty=2, lwd=2)

