# RStudio
#
# summary_stats.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-07-18 -- created
# 2014-11-19 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script reads the summary statistics file that is output by 
# model.py and analyzes the relationship between the observation data
# and the optimization parameters (both hyperbolic and linear fits)
# for the purposes of improving the dynamic parameterization.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 01. processing test7 results [13.07.22]
# 02. processing test9 results [13.10.01]
# 03. renamed "summary_stats.R" [13.10.11]
# 04. separating obs/ro_h/ro_l [14.03.11]
# 05. started looking at r2 [14.03.27]
# 06. added add_fapar function [14.06.03]
# 07. updated plotting: [14.06.03]
# --> removed estimate (crosses)
# --> divided Foo & Alpha by fAPAR
# 08. added add_meta, get_fapar_params, & fapar_plot functions [14.06.05] 
# 09. added LUE stats [14.11.19]
#
# ~~~~
# todo
# ~~~~
# 1. finish updating get_params for the new fAPAR augmented datasets
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### LIBRARIES ################################################################
# /////////////////////////////////////////////////////////////////////////////

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: get_params
# *
# * Input: data content (my.obj)
# *        field name indicator (my.val)
# *        flux tower station name (my.st)
# *
# * Return: list
# *
# * Features: This function returns the necessary plotting params 
# *           required for make_plot function
# ************************************************************************
get_params <- function(my.obj, my.val, my.st){
  # Return list:
  my.return.list <- list()
  #
  # Save station name:
  my.return.list$station <- my.st
  #
  # Optimization & estimation params:
  all.opt.fields <- c(
    "foo_opt_obs_h",   "foo_opt_ro_h",
    "alpha_opt_obs_h", "alpha_opt_ro_h",
    "alpha_opt_obs_l", "alpha_opt_ro_l",
    "r_opt_obs_h",     "r_opt_ro_h",
    "r_opt_obs_l",     "r_opt_ro_l"
  )
  all.est.fields <- c(
    "foo_est_obs_h",   "foo_est_ro_h",
    "alpha_est_obs_h", "alpha_est_ro_h",
    "alpha_est_obs_l", "alpha_est_ro_l",
    "r_est_obs_h",     "r_est_ro_h",
    "r_est_obs_l",     "r_est_ro_l"
  )
  #
  # Get field names:
  my.opt.field <- all.opt.fields[my.val]
  my.est.field <- all.est.fields[my.val]
  #
  my.return.list$opt.field <- my.opt.field
  my.return.list$est.field <- my.est.field
  #
  opt.data <- my.obj[,my.opt.field]
  est.data <- my.obj[,my.est.field]
  #
  # Save plotting data to return list:
  my.return.list$opt.data <- opt.data
  my.return.list$est.data <- est.data
  my.return.list$x.data <- my.obj$yrfr
  #
  # Get plotting limits:
  my.lims <- c(0,0)
  my.dy <- 0
  if (my.val == 1 | my.val == 2){
    # ~~~~~~~~~~~~~~~
    # FLUX PARAM: Foo
    # ~~~~~~~~~~~~~~~
    my.lims <- c(0,100)
    my.dy <- 10
    if (my.val == 1){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(italic(F[infinity])))
    } else {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(italic(F[infinity])))
    }
    # 
    # Look for max between limits:
    opt.max <- max(opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])])
    est.max <- max(est.data[which(est.data >= my.lims[1] & est.data <= my.lims[2])])
    # switch comment for plotting estimates:
    #my.max <- max(opt.max, est.max)
    my.max <- opt.max
    #
    # Adjust dy for small max values:
    if (my.max < 5){
      my.dy <- 1
    } else if (my.max < 18){
      my.dy <- 2
    } else if (my.max < 35){
      my.dy <- 5
    } else if (my.max > 80){
      my.dy <- 20
    } 
    #
    # Set max & min plotting values:
    min.y <- 0 - 0.2*my.dy
    max.y <- (
      floor(my.max+my.dy) -
        (floor(my.max+my.dy) %% (my.dy)) + 
        0.2*my.dy
    )
  } else if (my.val > 2 & my.val < 7) {
    # ~~~~~~~~~~~~~~~~~
    # FLUX PARAM: Alpha
    # ~~~~~~~~~~~~~~~~~
    my.lims <- c(0,1)
    my.dy <- 0.05
    if (my.val == 3){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(alpha))
    } else if (my.val == 4) {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(alpha))
    } else if (my.val == 5) {
      my.return.list$model <- 3
      my.return.list$expr <- expression(Linear~Flux~Param~Obs~(alpha))
    } else {
      my.return.list$model <- 4
      my.return.list$expr <- expression(Linear~Flux~Param~RO~(alpha))
    }
    #
    # Look for max between limits:
    opt.max <- max(opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])])
    est.max <- max(est.data[which(est.data >= my.lims[1] & est.data <= my.lims[2])])
    # switch comment for plotting estimates:
    #my.max <- max(opt.max, est.max)
    my.max <- opt.max
    #
    #
    # Adjust dy for small max values:
    if (my.max < 0.05){
      my.dy <- 0.01
    }
    #
    # Set plotting max & min:
    min.y <- 0 - 0.2*my.dy
    max.y <- (
      floor(100*my.max+100*my.dy) -
        (floor(100*my.max+100*my.dy) %% (100*my.dy))
    )
    max.y <- max.y / 100 + 0.2*my.dy
  } else if (my.val > 6 & my.val < 11) {
    # ~~~~~~~~~~~~~
    # FLUX PARAM: R
    # ~~~~~~~~~~~~~
    my.lims <- c(-10, 50)
    my.dy <- 5
    if (my.val == 7){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(italic(R)))
    } else if (my.val == 8) {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(italic(R)))
    } else if (my.val == 9){
      my.return.list$model <- 3
      my.return.list$expr <- expression(Linear~Flux~Param~Obs~(italic(R)))
    } else {
      my.return.list$model <- 4
      my.return.list$expr <- expression(Linear~Flux~Param~RO~(italic(R)))
    }
    #
    # Look for max between limits:
    opt.max <- max(opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])])
    est.max <- max(est.data[which(est.data >= my.lims[1] & est.data <= my.lims[2])])
    #
    opt.min <- min(opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])])
    est.min <- min(est.data[which(est.data >= my.lims[1] & est.data <= my.lims[2])])
    # Switch comments for plotting estimates:
    #my.max <- max(opt.max, est.max)
    #my.min <- min(opt.min, est.min)
    my.max <- opt.max
    my.min <- opt.min
    #
    # Adjust dy for small max values:
    if (my.max < 4){
      my.dy <- 1
    } else if (my.max < 10){
      my.dy <- 2
    }
    #
    # Set plotting max & min:
    min.y <- floor(my.min) - 
      (floor(my.min) %% (my.dy)) - 
      0.2*my.dy
    max.y <- (
      floor(my.max+my.dy) -
        (floor(my.max+my.dy) %% (my.dy)) +
        0.2*my.dy
    )
  }
  #
  my.return.list$plot.lims <- c(min.y,max.y)
  #
  # Sequency for ticks and labels:
  my.return.list$at <- seq((min.y + 0.2*my.dy), (max.y - 0.2*my.dy), my.dy)
  #
  # Save plotting character for optimized/model selected data:
  #   pch = 1 (open circles) 
  #   pch = 16 (filled circles)
  # switch comments for highlighting model selection
  #my.return.list$pch.vals <- 1 + 15*(my.obj$model_select == my.return.list$model)
  my.return.list$pch.vals <- 16
  #
  my.return.list
}

# ************************************************************************
# * Name: get_fapar_params
# *
# * Input: data content (my.obj)
# *        field name indicator (my.val)
# *        flux tower station name (my.st)
# *
# * Return: list of plotting params
# *           station :: flux tower name (character)
# *           opt.field :: parameter field name (character)
# *           opt.data :: list of parameter values
# *           fapar.data :: list of parameter values divided by fapar
# *           x.data :: list of x values (fractions of years)
# *           y.line :: mtext line location for Y-axis label
# *           model :: identifier of model selection
# *           expr :: Y-axis plotting label (expression)
# *           plot.lims :: list of Y-axis plotting limits
# *           at :: list of Y-axis tick mark locations
# *           pch.vals :: list of plotting character codes
# *           climate :: climate ID for flux tower (character)
# *           class :: vegetation class ID for flux tower (character)
# *
# * Features: This function returns the necessary plotting params 
# *           divided by fAPAR that are required for make_plot function
# ************************************************************************
get_fapar_params <- function(my.obj, my.val, my.st){
  # Return list:
  my.return.list <- list()
  #
  # Save station name:
  my.return.list$station <- my.st
  #
  # Optimization & estimation params:
  all.opt.fields <- c(
    "foo_opt_obs_h",   "foo_opt_ro_h",
    "alpha_opt_obs_h", "alpha_opt_ro_h",
    "alpha_opt_obs_l", "alpha_opt_ro_l",
    "r_opt_obs_h",     "r_opt_ro_h",
    "r_opt_obs_l",     "r_opt_ro_l"
  )
  #
  # Get field name:
  my.opt.field <- all.opt.fields[my.val]
  #
  my.return.list$opt.field <- my.opt.field
  #
  opt.data <- my.obj[,my.opt.field]
  opt.fapar.data <- opt.data / my.obj$fAPAR
  #
  # Save plotting data:
  my.return.list$opt.data <- opt.data
  my.return.list$fapar.data <- opt.fapar.data
  my.return.list$x.data <- my.obj$yrfr
  #
  # Get plotting limits:
  my.lims <- c(0,0)
  my.dy <- 0
  if (my.val == 1 | my.val == 2){
    # ~~~~~~~~~~~~~~~
    # FLUX PARAM: Foo
    # ~~~~~~~~~~~~~~~
    my.return.list$y.line <- 2
    my.lims <- c(0,100)
    my.dy <- 10
    if (my.val == 1){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(italic(F[infinity])))
    } else {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(italic(F[infinity])))
    }
    # 
    # Look for max between limits:
    opt.max <- max(
      opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    fapar.max <- max(
      opt.fapar.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    my.max <- max(opt.max,fapar.max)
    #
    # Adjust dy for small max values:
    if (my.max < 5){
      my.dy <- 1
    } else if (my.max < 18){
      my.dy <- 2
    } else if (my.max < 35){
      my.dy <- 5
    } else if (my.max > 500) {
      my.dy <- 100
    } else if (my.max > 300) {
      my.dy <- 50
    } else if (my.max > 180){
      my.dy <- 40
    } else if (my.max > 80){
      my.dy <- 20
    } 
    #
    # Set max & min plotting values:
    min.y <- 0 - 0.2*my.dy
    max.y <- (
      floor(my.max+my.dy) -
        (floor(my.max+my.dy) %% (my.dy)) + 
        0.2*my.dy
    )
  } else if (my.val > 2 & my.val < 7) {
    # ~~~~~~~~~~~~~~~~~
    # FLUX PARAM: Alpha
    # ~~~~~~~~~~~~~~~~~
    my.return.list$y.line <- 3
    my.lims <- c(0,1)
    my.dy <- 0.05
    if (my.val == 3){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(alpha))
    } else if (my.val == 4) {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(alpha))
    } else if (my.val == 5) {
      my.return.list$model <- 3
      my.return.list$expr <- expression(Linear~Flux~Param~Obs~(alpha))
    } else {
      my.return.list$model <- 4
      my.return.list$expr <- expression(Linear~Flux~Param~RO~(alpha))
    }
    #
    # Look for max between limits:
    opt.max <- max(
      opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    fapar.max <- max(
      opt.fapar.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    my.max <- max(opt.max, fapar.max)
    #
    # Adjust dy for small max values:
    if (my.max < 0.05){
      my.dy <- 0.01
    } else if (my.max < 0.1) {
      my.dy <- 0.02
    } else if (my.max > 10) {
      my.dy <- 2
    } else if (my.max > 4.5) {
      my.dy <- 1
    } else if (my.max > 1.5) {
      my.dy <- 0.5
    } else if (my.max > 0.9) {
      my.dy <- 0.2
    } else if (my.max > 0.4){
      my.dy <- 0.1
    }
    #
    # Set plotting max & min:
    min.y <- 0 - 0.2*my.dy
    max.y <- (
      floor(100*my.max+100*my.dy) -
        (floor(100*my.max+100*my.dy) %% (100*my.dy))
    )
    max.y <- max.y / 100 + 0.2*my.dy
  } else if (my.val > 6 & my.val < 11) {
    # ~~~~~~~~~~~~~
    # FLUX PARAM: R
    # ~~~~~~~~~~~~~
    my.return.list$y.line <- 2
    my.lims <- c(-10, 50)
    my.dy <- 5
    if (my.val == 7){
      my.return.list$model <- 1
      my.return.list$expr <- expression(Hyper~Flux~Param~Obs~(italic(R)))
    } else if (my.val == 8) {
      my.return.list$model <- 2
      my.return.list$expr <- expression(Hyper~Flux~Param~RO~(italic(R)))
    } else if (my.val == 9){
      my.return.list$model <- 3
      my.return.list$expr <- expression(Linear~Flux~Param~Obs~(italic(R)))
    } else {
      my.return.list$model <- 4
      my.return.list$expr <- expression(Linear~Flux~Param~RO~(italic(R)/fAPAR))
    }
    #
    # Look for max between limits:
    opt.max <- max(
      opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    fapar.max <- max(
      opt.fapar.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    my.max <- max(opt.max, fapar.max)
    #
    opt.min <- min(
      opt.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    fapar.min <- min(
      opt.fapar.data[which(opt.data >= my.lims[1] & opt.data <= my.lims[2])]
    )
    my.min <- min(opt.min, fapar.min)
    #
    # Adjust dy for small max values:
    if (my.max < 2){
      my.dy <- 1
    } else if (my.max < 10){
      my.dy <- 2
    }
    #
    # Set plotting max & min:
    min.y <- floor(my.min) - 
      (floor(my.min) %% (my.dy)) - 
      0.2*my.dy
    max.y <- (
      floor(my.max+my.dy) -
        (floor(my.max+my.dy) %% (my.dy)) +
        0.2*my.dy
    )
  }
  #
  my.return.list$plot.lims <- c(min.y,max.y)
  #
  # Sequence for ticks and labels:
  my.return.list$at <- seq((min.y + 0.2*my.dy), (max.y - 0.2*my.dy), my.dy)
  #
  # Save plotting characters 
  #   par=1 (open circles)
  #   par=3 (crosses)
  #   par=16 (filled circles)
  #   par=8 (asterisk)
  my.return.list$pch.vals <- c(1,8)
  #
  my.return.list
}

# ************************************************************************
# * Name: make_plot
# *
# * Input: parameters list (my.params)
# *        plot color (my.color)
# *
# * Return: None.
# *
# * Features: This function plots optimization and estimation parameters
# *           with respect to time
# ************************************************************************
make_plot <- function(my.params, my.color='gray'){
  # Extract optimized and estimated y axis points:
  y.opt <- my.params$opt.data
  y.est <- my.params$est.data
  #
  # Set out-of-bounds data to NA
  y.opt[which(y.opt < my.params$plot.lims[1])] <- NA
  y.opt[which(y.opt > my.params$plot.lims[2])] <- NA
  #
  # Define plotting character for optimization points
  #   If model select -> pch=16
  #   else -> pch=1
  my.pch.vals <- my.params$pch.vals
  #
  # Setup plotting space:
  par(mar=c(4.5,4.5,1,1))
  plot(
    my.params$x.data,
    y.opt,
    xlim=c(2001.867,2007.133),
    ylim=my.params$plot.lims,
    type='n',
    xlab = NA,
    ylab = NA,
    xaxs = 'i',
    yaxs = 'i',
    axes=F
  )
  #
  # Add optimization data:
  points(
    my.params$x.data,
    y.opt,
    col = my.color,
    pch = my.pch.vals
  )
#   lines(
#     my.params$x.data,
#     y.opt,
#     col = my.color,
#     lwd=1.5,
#     lty=1
#   )
  #
  # Add estimation data (optional):
  #points(
  #  my.obj$yrfr,
  #  y.est,
  #  col = 'black',
  #  pch = 3
  #)
  #
  # Dress-up plot with axes and labels
  box(lwd=2)
  axis(
    at = seq(2002, 2007, 1),
    side = 1,
    labels = seq(2002, 2007, 1),
    pos = my.params$plot.lims[1],
    tck = -0.02,
    mgp = c(0,0.5,0)
  )
  axis(
    at = seq(2002, 2007, 1/12), 
    side = 1, 
    labels = FALSE, 
    pos = my.params$plot.lims[1], 
    tck = -0.01, 
  )
  axis(
    at = my.params$at,
    side=2, 
    labels = my.params$at,
    pos = 2001.867,
    las=1, 
    tck=-0.01,
    mgp = c(0,0.5,0)
  )
  #
  # Added axis title & legend:
  mtext(side=1, my.params$station, line=3)
  mtext(side=2, as.expression(my.params$expr), line=2)
#   legend(
#     'top',
#     legend=c('Optimized', 'Selected'),
#     col=c(my.color, my.color),
#     pch=c(1, 16),
#     bty='o',
#     bg = 'white',
#     box.col = 'white',
#     horiz=TRUE,
#     inset = 0.01
#   )
}

# ************************************************************************
# * Name: fapar_plot
# *
# * Input: parameters list (my.params)
# *        plot colors (my.color)
# *
# * Return: None.
# *
# * Features: This function plots optimization parameters and those
# *           divided by fAPAR with respect to time
# ************************************************************************
fapar_plot <- function(my.params, my.color=c('gray','black')){
  # Extract optimized and estimated y axis points:
  y.opt <- my.params$opt.data
  y.fapar <- my.params$fapar.data
  #
  # Set out-of-bounds data to NA
  y.opt[which(y.opt < my.params$plot.lims[1])] <- NA
  y.opt[which(y.opt > my.params$plot.lims[2])] <- NA
  #
  # Get plotting characters
  my.pch.vals <- my.params$pch.vals
  #
  # Setup plotting space:
  par(mar=c(4.5,4.5,1,1))
  plot(
    my.params$x.data,
    y.opt,
    xlim=c(2001.867,2007.133),
    ylim=my.params$plot.lims,
    type='n',
    xlab = NA,
    ylab = NA,
    xaxs = 'i',
    yaxs = 'i',
    axes=F
  )
  #
  # Add optimization data:
  points(
    my.params$x.data,
    y.opt,
    col = my.color[1],
    pch = my.pch.vals[1]
  )
  #
  # Add fapar data:
  points(
    my.params$x.data,
    y.fapar,
    col = my.color[2],
    pch = my.pch.vals[2]
  )
  #
  # Dress-up plot with axes and labels
  box(lwd=2)
  axis(
    at = seq(2002, 2007, 1),
    side = 1,
    labels = seq(2002, 2007, 1),
    pos = my.params$plot.lims[1],
    tck = -0.02,
    mgp = c(0,0.5,0)
  )
  axis(
    at = seq(2002, 2007, 1/12), 
    side = 1, 
    labels = FALSE, 
    pos = my.params$plot.lims[1], 
    tck = -0.01, 
  )
  axis(
    at = my.params$at,
    side=2, 
    labels = my.params$at,
    pos = 2001.867,
    las=1, 
    tck=-0.01,
    mgp = c(0,0.5,0)
  )
  #
  # Added axis title & legend:
  mtext(side=1, my.params$station, line=3)
  mtext(side=2, as.expression(my.params$expr), line=my.params$y.line)
  mtext(
    paste(
      "Climate ID: ", as.character(my.params$climate), "\n",
      "Class ID: ", as.character(my.params$class)
    ),
    side=1,
    line=3,
    adj=1,
    cex=0.8
  )
  #   legend(
  #     'top',
  #     legend=c('Optimized', 'Selected'),
  #     col=c(my.color, my.color),
  #     pch=c(1, 16),
  #     bty='o',
  #     bg = 'white',
  #     box.col = 'white',
  #     horiz=TRUE,
  #     inset = 0.01
  #   )
}

# ************************************************************************
# * Name: add_fapar
# *
# * Input: data frame with content (my.obj)
# *        time series data directory (ts.path)
# *
# * Return: data frame with content + fapar
# *
# * Features: This function reads the LUE output and appends the 
# *           associated monthly fAPAR to a new column in the data frame
# ************************************************************************
add_fapar <- function(my.obj, ts.path){
  # Add column for fAPAR:
  my.obj$fAPAR <- 0*seq(from=1, to=dim(my.obj)[1], by=1)
  my.obj$fAPAR <- NA
  #
  # Get LUE file:
  ts.file <- list.files(path=ts.path, pattern = "timeseries.txt$")
  ts.content <- read.csv(paste(ts.path, ts.file, sep=""), header=T)
  #
  # Filter for station:
  my.station <- as.character(my.obj$name)[1]
  ts.content <- ts.content[which(ts.content$station == my.station),]
  #
  # Find months & add fAPAR data:
  for (my.month in as.character(my.obj$Timestamp)){
    my.obj$fAPAR[which(as.character(my.obj$Timestamp) == my.month)] <- 
      ts.content$fpar[which(as.character(ts.content$timestamp) == my.month)]
  }
  #
  my.obj
}

# ************************************************************************
# * Name: add_meta
# *
# * Input: plotting parameters (my.params)
# *        meta file directory (meta.dir)
# *
# * Return: updated plotting parameters (my.params)
# *
# * Features: This function retrieves meta data (vegetation & climate) 
# *           for a specific flux tower and adds it to plotting 
# *           parameters
# ************************************************************************
add_meta <- function(my.params, meta.dir){
  # Read meta data for specific flux tower
  meta.file <- list.files(path=meta.dir, pattern="Fluxdata_Met-Data*")
  meta.data <- read.csv(
    paste(meta.dir, meta.file, sep=""),
    header=T
  )
  #
  # Find row for flux tower:
  meta.row <- which(meta.data$stationid == my.params$station)
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
  # Add meta data:
  my.params$climate <- meta.clim
  my.params$class <- meta.veg
  my.params
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### DEFINITIONS ##############################################################
# /////////////////////////////////////////////////////////////////////////////
mac = F
if (mac){
  in.file.path=paste("/Users/twdavis/Dropbox/Work/Imperial/flux/results/",
                     "2002-06/summary_stats/", sep="")
} else {
  in.file.path=paste("/home/user/Dropbox/Work/Imperial/flux/results/",
                     "2002-06/summary_stats/", sep="")
}
stats.file = list.files(path = in.file.path, 
                        pattern = "summary_statistics-v24*")
content <- read.csv(paste(in.file.path,stats.file,sep=""), header=T);
content.bak <- content

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### MAIN #####################################################################
# /////////////////////////////////////////////////////////////////////////////
## TEST 6
# BEST FIT: fh = 3.64 * std.n (r2 = 0.84)
# BEST FIT: ah = 0.00491* std.n 
#                + 0.00235 * krt.p 
#                - 0.00233 * krt.n (r2 = 0.76)
# BEST FIT: rh = 0.0224 * ave.n + 0.661 * std.n (r2 = 0.82)
# BEST FIT: al = 0.000941 * krt.p 
#                + 0.00292 * std.n 
#                - 0.000675 * max.n (r2 = 0.88)
# BEST FIT: rl = 0.818 * std.n 
#                + 0.71 * ave.n 
#                + 0.00537 * ave.p 
#                - 0.00689 * std.p (r2 = 0.89)
#
## TEST 7
# BEST FIT: fh = 3.83 * std.n (r2 = 0.86)
# BEST FIT: ah = 1.96*(max.n.ah - min.n.ah)/(max.p.ah - min.p.ah) (r2 = 0.71)
# BEST FIT: rh = 0.69 * std.n (r2 = 0.87)
# BEST FIT: al = 0.672*(max.n.al - min.n.al)/(max.p.al - min.p.al) (r2 = 0.93)
# BEST FIT: rl = 0.899 * std.n 
#                + 0.827 * ave.n 
#                + 0.00628 * ave.p 
#                - 0.008 * std.p (r2 = 0.90)
#
## TEST 9
# BEST FIT: fh = 3.8367 * std.n (r2adj = 0.863)
# BEST FIT: ah = 1.9935*(max.n.ah - min.n.ah)
#                /(max.p.ah - min.p.ah) (r2adj = 0.718)
# BEST FIT: rh = 0.6931 * std.n (r2adj = 0.875)
# BEST FIT: al = 0.6710*(max.n.al - min.n.al)
#                /(max.p.al - min.p.al) (r2adj = 0.93)
# BEST FIT: rl = 0.8974 * std.n 
#                + 0.8400 * ave.n 
#                + 0.006344 * ave.p 
#                - 0.007969 * std.p (r2adj = 0.907)
#
# TEST 10 (146 stations / 5 years)
# BEST FIT: fh = 3.7113 (+/- 0.0224) * std.n (r2adj = 0.874)
# BEST FIT: ah = 1.8301 (+/- 0.01807)*(max.n.ah - min.n.ah)
#                /(max.p.ah - min.p.ah) (r2adj = 0.712)
# BEST FIT: rh = 0.5804 (+/- 0.00394) * std.n (r2adj = 0.827)
# BEST FIT: al = 0.6114 (+/- 0.00286)*(max.n.al - min.n.al)
#                /(max.p.al - min.p.al) (r2adj = 0.945)
# BEST FIT: rl = 0.7491 (+/- 0.0108) * std.n 
#                + 0.6159 (+/- 0.0100) * ave.n 
#                + 0.003756 (+/- 0.00017) * ave.p 
#                - 0.005416 (+/- 0.00024) * std.p (r2adj = 0.784)


#### Model H - Foo ####
plot(content$foo_opt_ro_h)
plot(content$foo_opt_obs_h)
fh.obs.list <- which(content$foo_opt_obs_h <= 100 & content$foo_opt_obs_h > 0)
fh.ro.list <- which(content$foo_opt_ro_h <= 100 & content$foo_opt_ro_h > 0)
fh.obs <- content$foo_opt_obs_h[fh.obs.list]
fh.ro <- content$foo_opt_ro_h[fh.ro.list]
plot(fh.obs)
plot(fh.ro)

krt.p.fh.obs <- content$krt_ppfd_obs[fh.obs.list]
krt.n.fh.obs <- content$krt_nee_obs[fh.obs.list]
max.p.fh.obs <- content$max_ppfd_obs[fh.obs.list]
max.n.fh.obs <- content$max_nee_obs[fh.obs.list]
min.p.fh.obs <- content$min_ppfd_obs[fh.obs.list]
min.n.fh.obs <- content$min_nee_obs[fh.obs.list]
skw.p.fh.obs <- content$skw_ppfd_obs[fh.obs.list]
skw.n.fh.obs <- content$skw_nee_obs[fh.obs.list]
ave.n.fh.obs <- content$ave_nee_obs[fh.obs.list]
ave.p.fh.obs <- content$ave_ppfd_obs[fh.obs.list]
std.n.fh.obs <- content$std_nee_obs[fh.obs.list]
std.p.fh.obs <- content$std_ppfd_obs[fh.obs.list]

krt.p.fh.ro <- content$krt_ppfd_ro_h[fh.ro.list]
krt.n.fh.ro <- content$krt_nee_ro_h[fh.ro.list]
max.p.fh.ro <- content$max_ppfd_ro_h[fh.ro.list]
max.n.fh.ro <- content$max_nee_ro_h[fh.ro.list]
min.p.fh.ro <- content$min_ppfd_ro_h[fh.ro.list]
min.n.fh.ro <- content$min_nee_ro_h[fh.ro.list]
skw.p.fh.ro <- content$skw_ppfd_ro_h[fh.ro.list]
skw.n.fh.ro <- content$skw_nee_ro_h[fh.ro.list]
ave.n.fh.ro <- content$ave_nee_ro_h[fh.ro.list]
ave.p.fh.ro <- content$ave_ppfd_ro_h[fh.ro.list]
std.n.fh.ro <- content$std_nee_ro_h[fh.ro.list]
std.p.fh.ro <- content$std_ppfd_ro_h[fh.ro.list]

cbind(
  krt.p.fh, krt.n.fh, 
  max.p.fh, max.n.fh, 
  min.p.fh, min.n.fh, 
  skw.p.fh, skw.n.fh, 
  ave.n.fh, ave.p.fh, 
  std.n.fh, std.p.fh)

plot(std.n.fh.obs, fh.obs, xlab='Std. Dev. NEE', ylab='Foo')
plot(std.n.fh.ro, fh.ro, xlab='Std. Dev. NEE', ylab='Foo')

fit_fh_obs <- lm(fh.obs ~ 0 + std.n.fh.obs)
summary(fit_fh_obs)
hist(fit_fh_obs$residuals)

fit_fh_ro <- lm(fh.ro ~ 0 + std.n.fh.ro)
summary(fit_fh_ro)
hist(fit_fh_ro$residuals)

fh <- c(fh.obs, fh.ro)
std.n.fh <- c(std.n.fh.obs, std.n.fh.ro)
fit_fh <- lm(fh ~ 0 + std.n.fh)
summary(fit_fh)
hist(fit_fh$residuals)
#
# BEST FIT: fh = 3.5883 (+/- 0.01643) * std.n (r2adj = 0.858)

#### Model H - ALPHA ####
plot(content$alpha_opt_ro_h)
plot(content$alpha_opt_obs_h)

ah.obs.list <- which(content$alpha_opt_obs_h < 1 & content$alpha_opt_obs_h > 0)
ah.ro.list <- which(content$alpha_opt_ro_h < 1 & content$alpha_opt_ro_h > 0)

ah.obs <- content$alpha_opt_obs_h[ah.obs.list]
ah.ro <- content$alpha_opt_ro_h[ah.ro.list]

plot(ah.obs)
plot(ah.ro)

krt.p.ah.obs <- content$krt_ppfd_obs[ah.obs.list]
krt.n.ah.obs <- content$krt_nee_obs[ah.obs.list]
skw.p.ah.obs <- content$skw_ppfd_obs[ah.obs.list]
skw.n.ah.obs <- content$skw_nee_obs[ah.obs.list]
max.p.ah.obs <- content$max_ppfd_obs[ah.obs.list]
min.p.ah.obs <- content$min_ppfd_obs[ah.obs.list]
max.n.ah.obs <- content$max_nee_obs[ah.obs.list]
min.n.ah.obs <- content$min_nee_obs[ah.obs.list]
ave.n.ah.obs <- content$ave_nee_obs[ah.obs.list]
ave.p.ah.obs <- content$ave_ppfd_obs[ah.obs.list]
std.n.ah.obs <- content$std_nee_obs[ah.obs.list]
std.p.ah.obs <- content$std_ppfd_obs[ah.obs.list]

krt.p.ah.ro <- content$krt_ppfd_ro_h[ah.ro.list]
krt.n.ah.ro <- content$krt_nee_ro_h[ah.ro.list]
skw.p.ah.ro <- content$skw_ppfd_ro_h[ah.ro.list]
skw.n.ah.ro <- content$skw_nee_ro_h[ah.ro.list]
max.p.ah.ro <- content$max_ppfd_ro_h[ah.ro.list]
min.p.ah.ro <- content$min_ppfd_ro_h[ah.ro.list]
max.n.ah.ro <- content$max_nee_ro_h[ah.ro.list]
min.n.ah.ro <- content$min_nee_ro_h[ah.ro.list]
ave.n.ah.ro <- content$ave_nee_ro_h[ah.ro.list]
ave.p.ah.ro <- content$ave_ppfd_ro_h[ah.ro.list]
std.n.ah.ro <- content$std_nee_ro_h[ah.ro.list]
std.p.ah.ro <- content$std_ppfd_ro_h[ah.ro.list]

range.n.ah.obs <- max.n.ah.obs - min.n.ah.obs
range.p.ah.obs <- max.p.ah.obs - min.p.ah.obs
ratio.np.ah.obs <- range.n.ah.obs / range.p.ah.obs

range.n.ah.ro <- max.n.ah.ro - min.n.ah.ro
range.p.ah.ro <- max.p.ah.ro - min.p.ah.ro
ratio.np.ah.ro <- range.n.ah.ro / range.p.ah.ro

# Secondary outlier removal
ah.obs[which(ratio.np.ah.obs > 0.1)] <- NA
ratio.np.ah.obs[which(ratio.np.ah.obs > 0.1)] <- NA
ah.obs[which(ratio.np.ah.ro > 0.1)] <- NA
ratio.np.ah.ro[which(ratio.np.ah.ro > 0.1)] <- NA

scatterplot3d(ave.n.ah, skw.p.ah, ah)
plot(
  ratio.np.ah.obs, 
  ah.obs, 
  xlab="Max(NEE)-Min(NEE) / Max(PPFD)-Min(PPFD)", 
  ylab='alpha')

fit_ah_obs <- lm(ah.obs ~ 0 + ratio.np.ah.obs)
summary(fit_ah_obs)
hist(fit_ah_obs$residuals)

plot(
  ratio.np.ah.ro, 
  ah.ro, 
  xlab="Max(NEE)-Min(NEE) / Max(PPFD)-Min(PPFD)", 
  ylab='alpha')

fit_ah_ro <- lm(ah.ro ~ 0 + ratio.np.ah.ro)
summary(fit_ah_ro)
hist(fit_ah_ro$residuals)

ratio.np.ah <- c(ratio.np.ah.obs, ratio.np.ah.ro)
ah <- c(ah.obs, ah.ro)
fit_ah <- lm(ah ~ 0 + ratio.np.ah)
summary(fit_ah)
hist(fit_ah$residuals)
#
# BEST FIT: ah = 1.4506 (+/- 0.01249)*(max.n.ah - min.n.ah)/(max.p.ah - min.p.ah) (r2adj = 0.6188)

#### Model H - R ####
plot(content$r_opt_obs_h)
plot(content$r_opt_ro_h)

rh.obs.list <- which(
  content$r_opt_obs_h > min(content$min_nee_obs) & content$r_opt_obs_h < max(content$max_nee_obs)
  )
rh.ro.list <- which(
  content$r_opt_ro_h > 0 & content$r_opt_ro_h < max(content$max_nee_ro_h)
  )

rh.obs <- content$r_opt_obs_h[rh.obs.list]
rh.ro <- content$r_opt_ro_h[rh.ro.list]

plot(rh.obs)
plot(rh.ro)

krt.p.rh <- content$krt_ppfd[rh.list]
krt.n.rh <- content$krt_nee[rh.list]
skw.p.rh <- content$skw_ppfd[rh.list]
skw.n.rh <- content$skw_nee[rh.list]
ave.p.rh <- content$ave_ppfd[rh.list]
ave.n.rh <- content$ave_nee[rh.list]
std.p.rh <- content$std_ppfd[rh.list]
std.n.rh <- content$std_nee[rh.list]
max.p.rh <- content$max_ppfd[rh.list]
max.n.rh <- content$max_nee[rh.list]
min.p.rh <- content$min_ppfd[rh.list]
min.n.rh <- content$min_nee[rh.list]
cbind(
  krt.p.rh, krt.n.rh, 
  max.p.rh, max.n.rh, 
  min.p.rh, min.n.rh, 
  skw.p.rh, skw.n.rh, 
  ave.n.rh, ave.p.rh, 
  std.n.rh, std.p.rh)

scatterplot3d(ave.n.rh, std.n.rh, rh)
plot(std.n.rh, rh, xlab='Std. Dev. NEE', ylab='R')

fit_rh <- lm(rh ~ 0 + std.n.rh)
summary(fit_rh)
hist(fit_rh$residuals)
#
# BEST FIT: rh = 0.5804 (+/- 0.00394) * std.n (r2adj = 0.827)

#### Model L - ALPHA ####
plot(content$alpha_opt_obs_l)
plot(content$alpha_opt_ro_l)

al.obs.list <- which(content$alpha_opt_obs_l > -1 & content$alpha_opt_obs_l < 1)
al.ro.list <- which(content$alpha_opt_ro_l > -1 & content$alpha_opt_ro_l < 1)

al.obs <- content$alpha_opt_obs_l[al.obs.list]
al.ro <- content$alpha_opt_ro_l[al.ro.list]

plot(al.obs)
plot(al.ro)

krt.p.al <- content$krt_ppfd[al.list]
krt.n.al <- content$krt_nee[al.list]
skw.p.al <- content$skw_ppfd[al.list]
skw.n.al <- content$skw_nee[al.list]
ave.p.al <- content$ave_ppfd[al.list]
ave.n.al <- content$ave_nee[al.list]
std.p.al <- content$std_ppfd[al.list]
std.n.al <- content$std_nee[al.list]
max.p.al <- content$max_ppfd[al.list]
max.n.al <- content$max_nee[al.list]
min.p.al <- content$min_ppfd[al.list]
min.n.al <- content$min_nee[al.list]

ratio.np.al <- (max.n.al - min.n.al)/(max.p.al - min.p.al)

scatterplot3d(krt.p.al, krt.n.al, al)
plot(
  ratio.np.al, 
  al, 
  xlab="Max(NEE)-Min(NEE) / Max(PPFD)-Min(PPFD)", 
  ylab='alpha')

fit_al <- lm(al ~ 0 + ratio.np.al)
summary(fit_al)
hist(fit_al$residuals)
#
# BEST FIT: al = 0.6114 (+/- 0.00286)*(max.n.al - min.n.al)/(max.p.al - min.p.al) (r2adj = 0.945)


#### Model L - R ####
plot(content$r_opt_obs_l)
plot(content$r_opt_ro_l)

rl.obs.list <- which(
  content$r_opt_obs_l > min(content$min_nee_obs) 
  & content$r_opt_obs_l < max(content$max_nee_obs)
  )
rl.ro.list <- which(content$r2_ro_l >= 0.4)

rl.obs <- content$r_opt_obs_l[rl.obs.list]
rl <- content$r_opt_ro_l[rl.list]

plot(rl.obs)

krt.p.rl <- content$krt_ppfd[rl.list]
krt.n.rl <- content$krt_nee[rl.list]
skw.p.rl <- content$skw_ppfd[rl.list]
skw.n.rl <- content$skw_nee[rl.list]
ave.p.rl <- content$ave_ppfd[rl.list]
ave.n.rl <- content$ave_nee[rl.list]
std.p.rl <- content$std_ppfd[rl.list]
std.n.rl <- content$std_nee[rl.list]
max.p.rl <- content$max_ppfd[rl.list]
max.n.rl <- content$max_nee[rl.list]
min.p.rl <- content$min_ppfd[rl.list]
min.n.rl <- content$min_nee[rl.list]

scatterplot3d(krt.p.rl, krt.n.rl, rl)
fit_rl <- lm(rl ~ 0 + std.n.rl + ave.n.rl + ave.p.rl + std.p.rl)
summary(fit_rl)
fit.rl.obs = 0.749*std.n.rl + 0.6159*ave.n.rl + 0.003757*ave.p.rl - 0.005416*std.p.rl
plot(
  fit_rl$fitted, 
  fit.rl.obs, 
  xlab="Modeled", 
  ylab='Observed')
hist(fit_rl$residuals)
#
# BEST FIT: rl = 0.7491 (+/- 0.0108) * std.n 
#                + 0.6159 (+/- 0.0100) * ave.n 
#                + 0.003756 (+/- 0.00017) * ave.p 
#                - 0.005416 (+/- 0.00024) * std.p (r2adj = 0.784)

#### MODEL SELECTION ####
# Indexes of sites with problems:
# ~~~~ BE-Lon
belon.months <- c(
  "2006-10-01",
  "2005-10-01",
  "2005-03-01",
  "2006-09-01",
  "2005-09-01",
  "2004-04-01",
  "2006-04-01",
  "2006-05-01",
  "2005-08-01",
  "2004-05-01"
)
belon.indexes <- match(belon.months, as.character(content[which(content$name == "BE-Lon"),]$month))
belon.data <- content[which(content$name == "BE-Lon"),][belon.indexes,]

# ~~~~ CH-Oe2
choe2.months <- c(
  "2005-03-01",
  "2005-07-01",
  "2005-08-01"
)
choe2.indexes <- match(choe2.months, as.character(content[which(content$name == "CH-Oe2"),]$month))
choe2.data <- content[which(content$name == "CH-Oe2"),][choe2.indexes,]

# ~~~~ IE-Ca1
ieca1.months <- c(
  "2005-02-01",
  "2005-10-01",
  "2006-10-01",
  "2004-10-01",
  "2006-09-01",
  "2005-09-01",
  "2004-09-01",
  "2005-04-01",
  "2004-08-01"
)
ieca1.indexes <- match(ieca1.months, as.character(content[which(content$name == "IE-Ca1"),]$month))
ieca1.data <- content[which(content$name == "IE-Ca1"),][ieca1.indexes,]

# ~~~~ IT-BCi
itbci.months <- c(
  "2006-03-01",
  "2004-09-01",
  "2005-09-01",
  "2006-04-01",
  "2006-09-01"
)
itbci.indexes <- match(itbci.months, as.character(content[which(content$name == "IT-BCi"),]$month))
itbci.data <- content[which(content$name == "IT-BCi"),][itbci.indexes,]

# ~~~~ IT-LMa
itlma.months <- c(
  "2006-06-01",
  "2006-07-01"
)
itlma.indexes <- match(itlma.months, as.character(content[which(content$name == "IT-LMa"),]$month))
itlma.data <- content[which(content$name == "IT-LMa"),][itlma.indexes,]

# ~~~~ IT-Ro2
itro2.months <- c(
  "2002-03-01",
  "2002-04-01",
  "2003-04-01",
  "2004-04-01",
  "2005-04-01",
  "2004-09-01",
  "2006-04-01"
)
itro2.indexes <- match(itro2.months, as.character(content[which(content$name == "IT-Ro2"),]$month))
itro2.data <- content[which(content$name == "IT-Ro2"),][itro2.indexes,]

# ~~~~ US-FPe
usfpe.months <- c(
  "2002-05-01",
  "2003-09-01",
  "2004-03-01",
  "2002-04-01",
  "2002-06-01",
  "2006-07-01",
  "2006-08-01"
)
usfpe.indexes <- match(usfpe.months, as.character(content[which(content$name == "US-FPe"),]$month))
usfpe.data <- content[which(content$name == "US-FPe"),][usfpe.indexes,]
#

# Separate data indexes based on model selection:
mod.0 <- which(content$model_select == 0)
mod.1 <- which(content$model_select == 1)
mod.2 <- which(content$model_select == 2)
mod.3 <- which(content$model_select == 3)
mod.4 <- which(content$model_select == 4)

# Histograms of R-squared & Pearson's R
par(mfrow=c(1,2))
model.breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5)
ms.ht <- hist(content$model_select, main="", xlab="Model Selection v24",breaks=model.breaks)

r2.breaks <- c(0.1,0.3,0.5,0.7,0.9,1.1)
m1.ht <- hist(content$r2_obs_h[mod.1], main="", xlab="Model 1 R-squared",breaks=r2.breaks)
m2.ht <- hist(content$r2_ro_h[mod.2], main="", xlab="Model 2 R-squared",breaks=r2.breaks)
m3.ht <- hist(content$r2_obs_l[mod.3], main="", xlab="Model 3 R-squared",breaks=r2.breaks)
m4.ht <- hist(content$r2_ro_l[mod.4], main="", xlab="Model 4 R-squared",breaks=r2.breaks)

# Percent fits with R2 > 0
sum(m1.ht$counts[3:5])/sum(m1.ht$counts)
sum(m2.ht$counts[3:5])/sum(m2.ht$counts)
sum(m3.ht$counts[3:5])/sum(m3.ht$counts)
sum(m4.ht$counts[3:5])/sum(m4.ht$counts)
sum(c(m1.ht$counts[3:5],m2.ht$counts[3:5],m3.ht$counts[3:5],m4.ht$counts[3:5]))/sum(c(m1.ht$counts,m2.ht$counts,m3.ht$counts,m4.ht$counts))
sum(c(m1.ht$counts[3:5],m2.ht$counts[3:5],m3.ht$counts[3:5],m4.ht$counts[3:5]))/length(content$model_select)


#### RSQUARED ANALYSIS ####
length(which(content$r2_obs_h < 0.5 & content$r2_obs_h > 0.2))
length(which(content$r2_ro_h < 0.5 & content$r2_ro_h > 0.2))
hist(content$r2_obs_h[which(content$r2_obs_h >= 0)])
hist(content$r2_ro_h[which(content$r2_ro_h >= 0)])

length(which(content$r2_obs_l < 0.5 & content$r2_obs_l > 0.2))
length(which(content$r2_ro_l < 0.5 & content$r2_ro_l > 0.2))
hist(content$r2_obs_l[which(content$r2_obs_l >= 0)])
hist(content$r2_ro_l[which(content$r2_ro_l >= 0)])

length(which(content$r2_obs_h > 0.2 & content$r2_obs_h < 0.5 
      & content$r2_ro_h > 0.2 & content$r2_ro_h < 0.5
      & content$r2_obs_l > 0.2 & content$r2_obs_l < 0.5
      & content$r2_ro_l > 0.2 & content$r2_ro_l < 0.5
))

# The 176 months that completely failed:
failed.months <- paste(
  as.character(
    content$name[(
        which(content$r2_obs_h < 0 
              & content$r2_ro_h < 0 
              & content$r2_obs_l < 0 
              & content$r2_ro_l < 0
        )
      )]
  ),
  as.character(
    content$month[(
      which(content$r2_obs_h < 0 
            & content$r2_ro_h < 0 
            & content$r2_obs_l < 0 
            & content$r2_ro_l < 0
      )
    )]
  ),
  sep=": "
)

 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### VARIABILITY OF PARAMETERS ####
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Get unique station names:
all.stations <- unique(as.character(content$name))

# Save timestamp, year and month data:
content$Timestamp <- content$month
content$year <- as.numeric(substr(as.character(content$Timestamp),1, 4))
content$month <- as.numeric(substr(as.character(content$Timestamp),6,7))
content$yrfr <- content$year + (content$month - 1)/12

# fAPAR directory
fapar.dir <- paste(
  "/home/user/Dropbox/",       # Linux
  "/Users/twdavis/Dropbox/",   # Mac
  "Work/Imperial/flux/results/2002-06/time-series_analysis/",
  sep=""
)

# meta data file
meta.dir <- paste(
  "/home/user/Dropbox/",       # Linux
  "/Users/twdavis/Dropbox/",   # Mac
  "Work/Imperial/flux/data/psql-data/flux/",
  sep=""  
)

# Optimized and estimated parameters:
out.names <- c(
  "Flux-Partitioning_Hobs-Foo_v24.ps",
  "Flux-Partitioning_Hro-Foo_v24.ps",
  "Flux-Partitioning_Hobs-Alpha_v24.ps",
  "Flux-Partitioning_Hro-Alpha_v24.ps",
  "Flux-Partitioning_Lobs-Alpha_v24.ps",
  "Flux-Partitioning_Lro-Alpha_v24.ps",
  "Flux-Partitioning_Hobs-R_v24.ps",
  "Flux-Partitioning_Hro-R_v24.ps",
  "Flux-Partitioning_Lobs-R_v24.ps",
  "Flux-Partitioning_Lro-R_v24.ps"
)
out.fapar.names <- c(
  "Flux-Partitioning_Hobs-Foo-fPAR_v24.ps",
  "Flux-Partitioning_Hro-Foo-fPAR_v24.ps",
  "Flux-Partitioning_Hobs-Alpha-fPAR_v24.ps",
  "Flux-Partitioning_Hro-Alpha-fPAR_v24.ps",
  "Flux-Partitioning_Lobs-Alpha-fPAR_v24.ps",
  "Flux-Partitioning_Lro-Alpha-fPAR_v24.ps",
  "Flux-Partitioning_Hobs-R-fPAR_v24.ps",
  "Flux-Partitioning_Hro-R-fPAR_v24.ps",
  "Flux-Partitioning_Lobs-R-fPAR_v24.ps",
  "Flux-Partitioning_Lro-R-fPAR_v24.ps"
)
out.path <- in.file.path

p <- 6
f <- TRUE

if (f){
  postscript(
    file = paste(out.path,out.fapar.names[p], sep=""),
    width = 6, 
    height = 4, 
    title = out.names[p],
    paper = "special", 
    horizontal = FALSE,
    onefile = TRUE
  )
} else {
  postscript(
    file = paste(out.path,out.names[p], sep=""),
    width = 6, 
    height = 4, 
    title = out.names[p],
    paper = "special", 
    horizontal = FALSE,
    onefile = TRUE
  )
} 

for (station in all.stations){
  # Get content for specific station:
  st.content <- content[which(content$name == station),]
  #
  if (f){
    st.content <- add_fapar(st.content, fapar.dir)
    st.params <- get_fapar_params(st.content, p, station)
    st.params <- add_meta(st.params, meta.dir)
    fapar_plot(st.params, c('black','red'))
  } else {
    st.params <- get_params(st.content, p, station)
    make_plot(st.params, 'darkgray')
  }
}
dev.off()

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### LUE STATS ################################################################
# /////////////////////////////////////////////////////////////////////////////
#
# Summary:
# phio ~ (1.76e-1)*Gs_min + 
#        (1.70e-2)*GPP_std + 
#        -(6.89e-3)*K_min + 
#        -(3.53e-4)*Iabs_ave
#
# beta ~ -1.25 +
#        (2.52e2)*eta_max +
#        (4.84e3)*eta_ave +
#        -(7.27e2)*eta_std +
#        -(2.05e1)*eta_krt +
#        (3.06e-1)*Iabs_max +
#        -(7.54e-1)*Iabs_std +
#        -(3.65e1)*Iabs_skw +
#        -(5.53e0)*GPP_ave +
#        (4.64e1)*ca_max +
#        (3.29e3)*Gs_ave +
#        (3.77e1)*Gs_skw +
#        -(1.05e2)*K_ave +
#        -(3.33e-1)*D_max +
#        -(3.25e-1)*D_min +
#        (5.49e-1)*D_ave + 
#        (6.90e1)*D_skw
#
lue.file.path <- '/home/user/Dropbox/Work/Imperial/flux/results/2002-06/lue/'
lue.file <- list.files(path=lue.file.path, pattern="GePiSaT_nxgn*")[3]
lue.content <- read.csv(paste(lue.file.path, lue.file, sep=""), na.strings='-9999.000000')

# Filter bad fits:
lue.content <- lue.content[which(lue.content$r_sq > 0.4),]
#lue.content <- lue.content[which(lue.content$beta < 2e2),]

par(mar=c(4,4.5,1,1))
plot(sort(lue.content$phio_opt), pch=20, xlab=NA, ylab=NA, axes=F)
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=2, expression(italic(phi[o])~(modeled)), line=2)

par(mar=c(4.5,4.5,1,1))
hist(lue.content$phio_opt)

par(mar=c(4,4.5,1,1))
plot(sort(log(lue.content$beta)), pch=20, xlab=NA, ylab=NA, axes=F)
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=2, expression(log~italic(beta)~(modeled)), line=2)

hist(lue.content$beta_opt)

sort(lue.content$beta)[floor(length(lue.content$beta)/2)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict beta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta.content <- lue.content[which(lue.content$beta_p < 0.05),]
beta.content <- beta.content[which(beta.content$beta < 100),]

plot(lue.content$beta_est, lue.content$beta_opt, xlab=NA, 
     ylab=NA, axes=F, ylim=c(0.75,7.25), xlim=c(1.75,7.25))
abline(0, 1, lty=2, col='gray40')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, expression(italic(beta)~(estimated)), line=2)
mtext(side=2, expression(italic(beta)~(modelled)), line=2)

fit_beta <- lm(beta ~ eta_max + eta_min + eta_ave + eta_std + eta_skw + eta_krt + 
                 Iabs_max + Iabs_min + Iabs_ave + Iabs_std + Iabs_skw + Iabs_krt + 
                 GPP_max + GPP_min + GPP_ave + GPP_std + GPP_skw + GPP_krt + 
                 ca_max + ca_min + ca_ave + ca_std + ca_skw + ca_krt + 
                 Gs_max + Gs_min + Gs_ave + Gs_std + Gs_skw + Gs_krt + 
                 K_max + K_min + K_ave + K_std + K_skw + K_krt + 
                 D_max + D_min + D_ave + D_std + D_skw + D_krt,
               data = lue.content)

# rsqr ~ 0.47
lue.content <- lue.content[which(lue.content$beta < 2e2),]
fit_beta <- lm(beta ~ 0 + GPP_ave + GPP_skw + Iabs_min + eta_max + eta_ave + eta_skw,
               data = lue.content)

# rsqr ~ 0.44
lue.content <- lue.content[which(lue.content$beta < 2e2),]
fit_beta <- lm(beta ~ eta_ave + 
                 Iabs_max + Iabs_std + Iabs_skw +  
                 GPP_ave +
                 ca_max + 
                 Gs_ave + Gs_skw + 
                 K_ave +
                 D_max,
               data = lue.content)

# rsqr ~ 0.60
lue.content <- lue.content[which(lue.content$beta < 2e2),]
fit_beta <- lm(beta ~ eta_max + eta_ave + eta_std + eta_krt + 
                 Iabs_max + Iabs_std + Iabs_skw +  
                 GPP_ave +
                 ca_max + 
                 Gs_ave + Gs_skw + 
                 K_ave +
                 D_max + D_min + D_ave + D_skw,
               data = lue.content)


my.b <- beta.content$beta_opt
my.d <- (beta.content$D_max + beta.content$D_min)/2
my.k <- beta.content$K_ave
my.g <- (beta.content$Gs_max + beta.content$Gs_min)/2
my.n <- (beta.content$eta_max + beta.content$eta_min)/2
my.chi <- 0.5 + (0.5 - 0.9)*(2500 - my.d)/(-2500)
my.ba <- 1.6*my.d/(my.k + my.g)
my.ba <- my.ba*((my.chi - my.g)/(1.0 - my.chi))^2
my.ba <- my.ba*my.n

fit_beta <- lm(my.b ~ 0 + my.ba)
summary(fit_beta)
row.names(summary(fit_beta)$coefficients)[which(as.numeric(summary(fit_beta)$coefficients[,'Pr(>|t|)']) == max(as.numeric(summary(fit_beta)$coefficients[,'Pr(>|t|)'])))]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict phio
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(lue.content$phio_est, lue.content$phio_opt, xlab=NA, 
     ylab=NA, axes=F)
abline(0, 1, lty=2, col='gray40')
box(lwd=2)
axis(side=1, las=1, tck=-0.02, labels=NA)
axis(side=1, las=1, lwd=0, line=-0.4)
axis(side=2, las=1, tck=-0.02, labels=NA)
axis(side=2, las=1, lwd=0, line=-0.4)
mtext(side=1, expression(italic(phi[o])~(estimated)), line=2)
mtext(side=2, expression(italic(phi[o])~(modelled)), line=2)


fit_phio <- lm(phio ~ K_max + K_min + K_ave + K_std + K_skw + K_krt + 
                 Gs_max + Gs_min + Gs_ave + Gs_std + Gs_skw + Gs_krt + 
                 eta_max + eta_min + eta_ave + eta_std + eta_skw + eta_krt + 
                 Iabs_max + Iabs_min + Iabs_ave + Iabs_std + Iabs_skw + Iabs_krt + 
                 GPP_max + GPP_min + GPP_ave + GPP_std + GPP_skw + GPP_krt + 
                 ca_max + ca_min + ca_ave + ca_std + ca_skw + ca_krt + 
                 D_max + D_min + D_ave + D_std + D_skw + D_krt,
               data = lue.content)

# rsqr ~ 0.41
fit_phio <- lm(phio ~ Iabs_max + Iabs_skw + GPP_ave,
               data = lue.content)

# rsqr ~ 0.75
fit_phio <- lm(phio ~ 0 + Gs_min + Gs_std,
               data = lue.content)

# rsqr ~ 0.89
fit_phio <- lm(phio ~ 0 + K_min + Gs_min + Iabs_ave + GPP_std,
               data = lue.content)

summary(fit_phio)

# 
row.names(summary(fit_phio)$coefficients)[which(as.numeric(summary(fit_phio)$coefficients[,'Pr(>|t|)']) == max(as.numeric(summary(fit_phio)$coefficients[,'Pr(>|t|)'])))]
