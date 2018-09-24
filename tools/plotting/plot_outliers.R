# plot_outliers.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2012-10-10 -- created
# 2018-09-24 -- last updated
#
# ------------
# description:
# ------------
# This script reads the observation and outlier-free datasets output by
# model.py, plots them and highlights the observation pairs identified
# as outliers.
#
# NOTE: observation and outlier files are organized into their own separate
# station folders by using a script (e.g., file_handler-osx.pl,
# file_handler-win.pl, or file_handler-any.py) located in the GePiSaT
# Bitbucket repository (/toos/processing).
#
# ----------
# changelog:
# ----------
# - Added fix for files with no observations [18.09.24]
# - Updated folder hierarchy [18.09.24]
# - Added optim_params() function [14.03.18]
# - Added list.dirs function [14.03.18]
# - Addressed NAs in function calls (na.rm=TRUE)
# - Addressed the "None" in read.csv by replacing them with NAs
# - Created separate functions for plotting linear and hyperbolic plots


#### FUNCTIONS ####


# ************************************************************************
# Name: modelH
#
# Input:
#
# Return:
#
# Features:
# ************************************************************************
modelH <- function(x, foo, alpha, r) {
    # Define linear arguments:
    a <- (1.0*alpha*r - 1.0*alpha*foo);
    b <- 1.0*foo*r;
    c <- alpha;
    d <- foo;

    y <- (a*x + b)/(c*x + d)

    return(y)
}

# ************************************************************************
# Name: modelL
#
# Input:
#
# Return:
#
# Features: Returns the linearly modeled NEE based on PPFD
# ************************************************************************
modelL <- function(x, alpha, r) {
    # Define linear arguments:
    a <- -1.0*alpha;
    b <- r;

    y = a*x + b

    return(y)
}

# ************************************************************************
# Name: plot_outlier_modelL
#
# Input:
#
# Return:
#
# Features:
# ************************************************************************
plot_outlier_modelL <- function(obs, ro, alpha, r, year.month){
  # Plotting extents:
  max.y <- max(obs$nee_obs, ro$nee_obs_l, na.rm = TRUE)
  min.y <- min(obs$nee_obs, ro$nee_obs_l, na.rm = TRUE)
  max.x <- max(obs$ppfd_obs, ro$ppfd_obs_l, na.rm = TRUE)
  min.x <- 0

  # Save fitting line:
  x <- seq(
    min(obs$ppfd_obs, na.rm = TRUE),
    max(obs$ppfd_obs, na.rm = TRUE),
    length = 100
    )
  y <- modelL(x, alpha, r)

  # Plot
  par(mar = c(5.5,4.75,1,1))
  plot(
    obs$ppfd_obs,
    obs$nee_obs,
    type = "p",
    pch = 20,
    col = "red",
    main = NULL,
    sub = year.month,
    xlab = expression(paste("PPFD (", mu, mol %.% m^{-2} %.% s^{-1},")")),
    ylab = expression(paste("NEE (", mu, mol %.% m^{-2} %.% s^{-1},")")),
    #xlab = "",
    #ylab = "",
    #axes = FALSE,
    #frame.plot = TRUE,
    xlim = c(min.x,max.x),
    ylim = c(min.y,max.y)
  )
  box(lwd = 2)
  points(
    ro$ppfd_obs_l,
    ro$nee_obs_l,
    pch = 20,
    col = "gray"
  )
  #points(
  #  obs$ppfd_obs,
  #  obs$nee_obs,
  #  pch = 1,
  #  col = "black"
  #)
  lines(
    x,
    y,
    lty = 2,
    col = "blue",
    lwd = 3
  )
}

# ************************************************************************
# Name: plot_outlier_modelH
#
# Input:
#
# Return:
#
# Features:
# ************************************************************************
plot_outlier_modelH <- function(obs, ro, alpha, r, foo, year.month){
  # Save plotting extents:
  max.y <- max(obs$nee_obs, ro$nee_obs_h, na.rm = TRUE)
  min.y <- min(obs$nee_obs, ro$nee_obs_h, na.rm = TRUE)
  max.x <- max(obs$ppfd_obs, ro$ppfd_obs_h, na.rm = TRUE)
  min.x <- 0

  # Save fitting line:
  x <- seq(
    min(obs$ppfd_obs, na.rm = TRUE),
    max(obs$ppfd_obs, na.rm = TRUE),
    length = 100
    )
  y <- modelH(x, foo, alpha, r)

  # Plot:
  par(mar = c(5.5,4.75,1,1))
  plot(
    obs$ppfd_obs,
    obs$nee_obs,
    type = "p",
    pch = 20,
    col = "red",
    main = NULL,
    sub = year.month,
    xlab = expression(paste("PPFD (", mu, mol %.% m^{-2} %.% s^{-1},")")),
    ylab = expression(paste("NEE (", mu, mol %.% m^{-2} %.% s^{-1},")")),
    #xlab = "",
    #ylab = "",
    #axes = FALSE,
    #frame.plot = TRUE,
    xlim = c(min.x,max.x),
    ylim = c(min.y,max.y)
  )
  box(lwd = 2)
  points(
    ro$ppfd_obs_h,
    ro$nee_obs_h,
    pch = 20,
    col = "gray"
  )
  #points(
  #  obs$ppfd_obs,
  #  obs$nee_obs,
  #  pch = 1,
  #  col = "black"
  #)
  lines(
    x,
    y,
    lty = 2,
    col = "blue",
    lwd = 3
  )
}

# ************************************************************************
# Name: list.dirs
#
# Input:
#
# Return:
#
# Features:
#
# Ref:      J. Ulrich (2011), How to obtain a list of directories
#              within a directory, like list.files(), but instead
#              "list.dirs()", Stack Overflow,
#              http://stackoverflow.com/questions/4749783/how-to-obtain-
#              a-list-of-directories-within-a-directory-like-list-files-
#              but-i
# ************************************************************************
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {

    all <- list.files(path, pattern, all.dirs,
                      full.names, recursive = FALSE, ignore.case)
    return(all)
}

# ************************************************************************
# Name: optim_params
#
# Input:
#
# Return:
#
# Features: Returns optimization params table
# ************************************************************************
optim_params <- function(my.path, my.file){
    headerparts <- read.csv(
        paste(my.path, my.file, sep = ""), header = F)[1:4,1:6]

    colnames(headerparts) <- c("names", "alpha", "R", "Foo", "RMSE", "R2")
    myvars = setNames(data.frame(t(headerparts[,-1])), headerparts[,1])

    # Retrieve optimization parameters:
    mh.alpha <- as.numeric(
        levels(myvars$MH_opt)[as.numeric(myvars$MH_opt['alpha'])])
    mh.guess.alpha <- as.numeric(
        levels(myvars$MH_guess)[as.numeric(myvars$MH_guess['alpha'])])

    mh.guess.r <- as.numeric(
        levels(myvars$MH_guess)[as.numeric(myvars$MH_guess['R'])])
    mh.r <- as.numeric(
        levels(myvars$MH_opt)[as.numeric(myvars$MH_opt['R'])])

    mh.guess.foo <- as.numeric(
        levels(myvars$MH_guess)[as.numeric(myvars$MH_guess['Foo'])])
    mh.foo <- as.numeric(
        levels(myvars$MH_opt)[as.numeric(myvars$MH_opt['Foo'])])

    ml.guess.alpha <- as.numeric(
        levels(myvars$ML_guess)[as.numeric(myvars$ML_guess['alpha'])])
    ml.alpha <- as.numeric(
        levels(myvars$ML_opt)[as.numeric(myvars$ML_opt['alpha'])])

    ml.guess.r <- as.numeric(
        levels(myvars$ML_guess)[as.numeric(myvars$ML_guess['R'])])
    ml.r <- as.numeric(
        levels(myvars$ML_opt)[as.numeric(myvars$ML_opt['R'])])

    mh.r2 <- as.numeric(
        levels(myvars$MH_opt)[as.numeric(myvars$MH_opt['R2'])])
    ml.r2 <- as.numeric(
        levels(myvars$ML_opt)[as.numeric(myvars$ML_opt['R2'])])

    mh.rmse <- as.numeric(
        levels(myvars$MH_opt)[as.numeric(myvars$MH_opt['RMSE'])])
    ml.rmse <- as.numeric(
        levels(myvars$ML_opt)[as.numeric(myvars$ML_opt['RMSE'])])

    # Make a table:
    optimization.params <- matrix(c(mh.alpha,       ml.alpha,
                                    mh.guess.alpha, ml.guess.alpha,
                                    mh.r,           ml.r,
                                    mh.guess.r,     ml.guess.r,
                                    mh.foo,         NA,
                                    mh.guess.foo,   NA,
                                    mh.rmse,        ml.rmse,
                                    mh.r2,          ml.r2), ncol = 2,byrow = T)
    optimization.params <- as.data.frame(optimization.params)
    colnames(optimization.params) <- c("model_h", "model_l")
    rownames(optimization.params) <- c(
        "alpha","alpha_est",
        "r","r_est",
        "foo","foo_est",
        "rmse","r2"
    )

    return(optimization.params)
}


#### DEFINITIONS ####

base_dir <- "C:/Workspace/GePiSat/results/2003-14/"
part_dir <- paste(base_dir, "partition/", sep = "")
stations <- list.dirs(path = part_dir, pattern = "*-*")
out_path <- paste(base_dir, "plots/outlier_plots/", sep = "")


#### MAIN ####


# -------------------------------------------------- #
# -------------- LINEAR OUTLIERS ------------------- #
# -------------------------------------------------- #
# Step through each directory:
for (st.id in stations) {
  file.path <- paste(part_dir, st.id, "/", sep = "")

  # Find files in subdirectory:
  obs.files.all <- list.files(path = file.path, pattern = "*[:1:]_obs.txt");
  ro.files.all <- list.files(path = file.path, pattern = "*_ro.txt$");
  obs.num.files <- length(obs.files.all)
  ro.num.files <- length(ro.files.all)

  if (obs.num.files == ro.num.files) {
        num.files <- obs.num.files
    } else {
        cat(paste("Error with number of files in station", st.id))
        num.files <- 0
    }

  # Create postscript outfile:
  ps.file.name.l <- paste(st.id, "_modL-outliers.ps", sep = "")
  postscript(
    file = paste(out_path, ps.file.name.l, sep = ""),
    width = 6,
    height = 6,
    paper = "special",
    horizontal = FALSE,
    onefile = TRUE
  )

  # Step through each file in subdirectory:
  for (i in seq(num.files)) {
    # Obs and ro file names:
    obs.file <- obs.files.all[i]
    ro.file <- ro.files.all[i]

    obs.file.len <- length(count.fields(paste(file.path, obs.file, sep = "")))
    ro.file.len <- length(count.fields(paste(file.path, ro.file, sep = "")))

    if (obs.file.len > 5 & ro.file.len > 5) {
      # Save station info:
      station.name <- gsub("_.*$", "", obs.file)
      year.month <- gsub("(^.*_)(.*)(_obs\\.txt$)", "\\2", obs.file)

      # Get obs/ro params from file headers:
      obs.params <- optim_params(file.path, obs.file)
      ro.params <- optim_params(file.path, ro.file)

      # Read contents:
      obs.contents <- read.csv(
        paste(file.path, obs.file, sep = ""),
        header = T,
        skip = 4,
        na.strings = "None"
      )
      ro.contents <- read.csv(
        paste(file.path, ro.file, sep = ""),
        header = T,
        skip = 4,
        na.strings = "None"
      )


      # Plot LINEAR:
      plot_outlier_modelL(
        obs.contents,
        ro.contents,
        ro.params['alpha','model_l'],
        ro.params['r','model_l'],
        year.month
      )

      rm(
        obs.contents, ro.contents,
        obs.params,   ro.params
      )
    }
  }

  # Close postscript file:
  dev.off()
}

# -------------------------------------------------- #
# ------------ HYPERBOLIC OUTLIERS ----------------- #
# -------------------------------------------------- #
# Step through each directory:
for (st.id in stations) {
  file.path <- paste(part_dir, st.id, "/", sep = "")

  # Find files in subdirectory:
  obs.files.all <- list.files(path = file.path, pattern = "*[:1:]_obs.txt");
  ro.files.all <- list.files(path = file.path, pattern = "*_ro.txt$");
  obs.num.files <- length(obs.files.all)
  ro.num.files <- length(ro.files.all)
  if (obs.num.files == ro.num.files) {
    num.files <- obs.num.files
  } else {
    cat(paste("Error with number of files in station", st.id))
    num.files <- 0
  }

  # Create postscript outfile:
  ps.file.name.h <- paste(st.id, "_modH-outliers.ps", sep = "")
  postscript(
    file = paste(out_path, ps.file.name.h, sep = ""),
    width = 6,
    height = 6,
    paper = "special",
    horizontal = FALSE,
    onefile = TRUE
  )

  # Step through each file in subdirectory:
  for (i in seq(num.files)) {
    obs.file <- obs.files.all[i]
    ro.file <- ro.files.all[i]

    obs.file.len <- length(count.fields(paste(file.path, obs.file, sep = "")))
    ro.file.len <- length(count.fields(paste(file.path, ro.file, sep = "")))

    if (obs.file.len > 5 & ro.file.len > 5) {
      # Save station info:
      station.name <- gsub("_.*$", "", obs.file)
      year.month <- gsub("(^.*_)(.*)(_obs\\.txt$)", "\\2", obs.file)

      # Get obs/ro params from file headers:
      obs.params <- optim_params(file.path, obs.file)
      ro.params <- optim_params(file.path, ro.file)

      # Read contents:
      obs.contents <- read.csv(
        paste(file.path, obs.file, sep = ""),
        header = T,
        skip = 4,
        na.strings = "None"
      )
      ro.contents <- read.csv(
        paste(file.path, ro.file, sep = ""),
        header = T,
        skip = 4,
        na.strings = "None"
      )

      # Plot HYPER:
      plot_outlier_modelH(
        obs.contents,
        ro.contents,
        ro.params['alpha','model_h'],
        ro.params['r','model_h'],
        ro.params['foo','model_h'],
        year.month
      )

      rm(
        obs.contents, ro.contents,
        obs.params,   ro.params
      )
    }
  }

  # Close postscript file:
  dev.off()
}
