# plot_partitioning.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-03-12 -- created
# 2018-09-24 -- last updated
#
# ~~~~~~
# notes:
# ~~~~~~
# based on plot_outliers.R
# replaces query_plots.R
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script reads the observation and outlier-free datasets output by
# model.py, and plots the linear and hyperbolic partitioning.
#
# NOTE: observation and outlier files are organized into their own separate
# station folders by using a script (e.g., file_handler-osx.pl,
# file_handler-win.pl, or file_handler-any.py) located in the GePiSaT
# Bitbucket repository (/toos/processing).
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - Fixed issue with partition files with no observations [18.09.24]
# - Amended folder hierarchy for obs/ro files in same directory [18.09.24]
# - Updated plotting functions [14.11.09]
# - Added function headers [14.11.09]
# - Created optim_params() function [14.03.15]
# - Separated ro_h from ro_l in plots [14.03.15]
# - Started to separated ro_h from ro_l in plots [14.03.13]
# - Make stations read directory names [14.03.13]
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

  y <- a*x + b

  return(y)
}

# ************************************************************************
# * Name: optim_params
# *
# * Input:
# *
# * Return:
# *
# * Features: Returns optimization params table
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
                                  mh.r2,          ml.r2),ncol = 2,byrow = T)
  optimization.params <- as.data.frame(optimization.params)
  colnames(optimization.params) <- c("model_h","model_l")
  rownames(optimization.params) <- c(
    "alpha","alpha_est",
    "r","r_est",
    "foo","foo_est",
    "rmse","r2"
  )

  return(optimization.params)
}

# ************************************************************************
# Name: plot_obs
#
# Input: - str, observation file path (obs.file.path)
#        - str, current observation file name (cur.file)
#
# Return:
#
# Features: Produces a plot of NEE:PPFD partitioning (observations)
#
# Depends:  - optim_params
#           - modelH
#           - modelL
# ************************************************************************
plot_obs <- function(obs.file.path, cur.file){
  # Get station name and year:
  station.name <- gsub("_.*$", "", cur.file)
  station.year <- gsub("(^.*_)(.*)(_obs\\.txt$)", "\\2", cur.file)

  # Get optimization params:
  optimization.params <- optim_params(obs.file.path, cur.file)
  op_aL <- optimization.params['alpha','model_l']
  op_rL <- optimization.params['r','model_l']
  op_aH <- optimization.params['alpha','model_h']
  op_rH <- optimization.params['r','model_h']
  op_fH <- optimization.params['foo','model_h']
  op_R2L <- optimization.params['r2','model_l']
  op_R2H <- optimization.params['r2','model_h']

  # Expressions for plot:
  rpa = vector('expression',2)
  rpa[1] = substitute(
    expression(italic(F)[lin] == LINR-LINA~italic(Q)),
    list(LINA = format(op_aL,dig = 3),
         LINR = format(op_rL,dig = 3))
  )[2]
  rpa[2] = substitute(
    expression(
      italic(F)[hyp] == HYPR - (HYPAF~italic(Q))/(HYPA~italic(Q)+HYPF)
    ),
    list(HYPR = format(op_rH,dig = 3),
         HYPA = format(op_aH,dig = 3),
         HYPF = format(op_fH,dig = 3),
         HYPAF = format((op_aH*op_fH), dig = 3))
  )[2]

  rpb = vector('expression',2)
  rpb[1] = substitute(
    expression(italic(R)[lin]^2 == LR2),
    list(LR2 = format(op_R2L,dig = 3)))[2]
  rpb[2] = substitute(
    expression(italic(R)[hyp]^2 == HR2),
    list(HR2 = format(op_R2H,dig = 3)))[2]

  # Read in data:
  content <- read.csv(paste(obs.file.path,cur.file,sep = ""), header = T,
                      skip = 4, na.strings = "None")

  if (nrow(content) > 0) {
    # Get lines for fits:
    max_x <- max(content$ppfd_obs, na.rm = TRUE)
    mh.x <- seq(0, max_x, length = 100)
    mh.y <- modelH(mh.x, op_fH, op_aH, op_rH)
    ml.y <- modelL(mh.x, op_aL, op_rL)
    
    # min_y <- 1.1*min(content$nee_obs)  # -10
    # max_y <- 1.5*max(content$nee_obs)  # 6
    
    # Plot observations and fits:
    par(mar = c(4.5,4.5,1,1));
    plot(content$ppfd_obs,
         content$nee_obs,
         type = "p",
         pch = 20,
         col = "gray",
         xlab = NA,
         ylab = NA,
         # ylim = c(min_y, max_y),
         axes = F)
    box(lwd = 2)
    lines(mh.x,mh.y, lty = 1, col = "red", lwd = 3)
    lines(mh.x,ml.y, lty = 2, col = "blue", lwd = 3)
    
    axis(side = 1, las = 1, tck = -0.02, labels = NA)
    axis(side = 1, las = 1, lwd = 0, line = -0.4)
    axis(side = 2, las = 1, tck = -0.02, labels = NA)
    axis(side = 2, las = 1, lwd = 0, line = -0.4)
    mtext(side = 1,
          expression(paste("PPFD (", mu, mol %.% m^{-2} %.% s^{-1},")")),
          line = 2)
    mtext(side = 2,
          expression(paste("NEE (", mu, mol %.% m^{-2} %.% s^{-1},")")),
          line = 2)
    mtext(side = 1, station.year, line = 3)
    
    legend('bottomleft',
           legend = rpa,
           bty = 'n',
           inset = 0.02,
           y.intersp = 1.5)
    legend('topright',
           legend = rpb,
           bty = 'n',
           inset = 0.02,
           y.intersp = 1.5)
    legend("top",
           inset = 0.02,
           #y.intersp = 1.5,
           col = c("red","blue"),
           lty = c(1,2),
           lwd = c(3,3),
           legend = c("Model H","Model L"),
           bty = 'n',
           horiz = T)
  }
}

# ************************************************************************
# Name: plot_ro_h
#
# Input:
#
# Return:
#
# Features: Creates a plot of NEE:PPFD partitions
#           observations w/o outliers using hyperbolic model
# ************************************************************************
plot_ro_h <- function(ro.file.path, cur.file){
  # Get station name and year:
  station.name <- gsub("_.*$", "", cur.file)
  station.year <- gsub("(^.*_)(.*)(_ro\\.txt$)", "\\2", cur.file)

  # Get optimization parts
  optimization.params <- optim_params(ro.file.path, cur.file)
  op_fH <- optimization.params['foo','model_h']
  op_aH <- optimization.params['alpha','model_h']
  op_rH <- optimization.params['r','model_h']

  # Read in data:
  content <- read.csv(paste(ro.file.path, cur.file, sep = ""),
                      header = T, skip = 4, na.strings = "None")

  # Get lines for fits:
  max_x <- max(content$ppfd_obs_h, na.rm = TRUE)
  mh.x <- seq(0, max_x, length = 100)
  mh.y <- modelH(mh.x, op_fH, op_aH, op_rH)

  # min_y <- -10
  # max_y <- 8

  # Plot expressions:
  rpa = vector('expression', 1)
  rpa[1] = substitute(
    expression(italic(F) == HYPR - (HYPFA~italic(Q))/(HYPA~italic(Q)+HYPF)),
    list(
      HYPR = format(op_rH, dig = 3),
      HYPA = format(op_aH, dig = 3),
      HYPF = format(op_fH, dig = 3),
      HYPFA = format((op_fH*op_aH), dig = 3))
  )[2]
  rpb = vector('expression', 1)
  rpb[1] = substitute(
    expression(italic(R)^2 == HR2),
    list(HR2 = format(optimization.params['r2','model_h'], dig = 3))
  )[2]

  # Plot observations and fits:
  par(mar = c(4.5,4.5,1,1))
  plot(content$ppfd_obs_h,
       content$nee_obs_h,
       type = "p",
       pch = 20,
       col = "gray",
       xlab = NA,
       ylab = NA,
       # ylim = c(min_y, max_y),
       axes = F)
  box(lwd = 2)
  lines(mh.x,mh.y, lty = 1, col = "red", lwd = 3)

  axis(side = 1, las = 1, tck = -0.02, labels = NA)
  axis(side = 1, las = 1, lwd = 0, line = -0.4)
  axis(side = 2, las = 1, tck = -0.02, labels = NA)
  axis(side = 2, las = 1, lwd = 0, line = -0.4)
  mtext(side = 1,
        expression(paste("PPFD (", mu, mol %.% m^{-2} %.% s^{-1},")")),
        line = 2)
  mtext(side = 2,
        expression(paste("NEE (", mu, mol %.% m^{-2} %.% s^{-1},")")),
        line = 2)
  mtext(side = 1, station.year, line = 3)

  legend('bottomleft', legend = rpa, inset = 0.02, bty = 'n')
  legend('bottomright', legend = rpb, inset = 0.02, bty = 'n')
  legend("topright",
         inset = 0.02,
         y.intersp = 1.5,
         col = c("red"),
         lty = c(1),
         lwd = c(3),
         legend = c("Model H (ro)"),
         bty = 'n')
}

# ************************************************************************
# Name: plot_ro_l
#
# Input:
#
# Return:
#
# Features: Creates a plot of NEE:PPFD partitions
#           observations w/o outliers using linear model
# ************************************************************************
plot_ro_l <- function(ro.file.path, cur.file){
  # Get station name and year:
  station.name <- gsub("_.*$", "", cur.file)
  station.year <- gsub("(^.*_)(.*)(_ro\\.txt$)", "\\2", cur.file)

  # Get optimization params:
  optimization.params <- optim_params(ro.file.path, cur.file)
  op_aL <- optimization.params['alpha','model_l']
  op_rL <- optimization.params['r','model_l']

  # Read in data:
  content <- read.csv(paste(ro.file.path,cur.file,sep = ""), header = T,
                      skip = 4, na.strings = "None")

  # Get lines for fits:
  max_x <- max(content$ppfd_obs_l, na.rm = TRUE)
  ml.x <- seq(0, max_x, length = 100)
  ml.y <- modelL(ml.x, op_aL, op_rL)

  # min_y <- -10
  # max_y <- 8

  # Plot expressions:
  rpa = vector('expression', 1)
  rpa[1] = substitute(
    expression(italic(F) == LINR-LINA~italic(Q)),
    list(LINA = format(op_aL, dig = 3),
         LINR = format(op_rL, dig = 3))
  )[2]

  rpb = vector('expression', 1)
  rpb[1] = substitute(
    expression(italic(R)^2 == LR2),
    list(LR2 = format(optimization.params['r2','model_l'],dig = 3)))[2]

  # Plot observations and fits:
  par(mar = c(4.5,4.5,1,1))
  plot(content$ppfd_obs_l,
       content$nee_obs_l,
       type = "p",
       pch = 20,
       col = "gray",
       xlab = NA,
       ylab = NA,
       # ylim = c(min_y, max_y),
       axes = F)
  box(lwd = 2)
  lines(ml.x,ml.y, lty = 2, col = "blue", lwd = 3)
  axis(side = 1, las = 1, tck = -0.02, labels = NA)
  axis(side = 1, las = 1, lwd = 0, line = -0.4)
  axis(side = 2, las = 1, tck = -0.02, labels = NA)
  axis(side = 2, las = 1, lwd = 0, line = -0.4)
  mtext(side = 1,
        expression(paste("PPFD (", mu, mol %.% m^{-2} %.% s^{-1},")")),
        line = 2)
  mtext(side = 2,
        expression(paste("NEE (", mu, mol %.% m^{-2} %.% s^{-1},")")),
        line = 2)
  mtext(side = 1, station.year, line = 3)

  legend('bottomleft', legend = rpa, inset = 0.02, bty = 'n')
  legend('bottomright', legend = rpb, inset = 0.02, bty = 'n')
  legend("topright",
         inset = 0.02,
         y.intersp = 1.5,
         col = c("blue"),
         lty = c(2),
         lwd = c(3),
         legend = c("Model L (ro)"),
         bty = 'n')
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
#           within a directory, like list.files(), but instead
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


#### DEFINITIONS ####

base_dir <- "C:/Workspace/GePiSat/results/2003-14/"
part_dir <- paste(base_dir, "partition/", sep = "")
stations <- list.dirs(path = part_dir, pattern = "*-*")
out_path <- paste(base_dir, "plots/partition_plots/", sep = "")


#### MAIN ####


# -------------------------------------------------- #
# ---------------- OBSERVATIONS -------------------- #
# -------------------------------------------------- #
# Step through each directory:
for (st.id in stations) {
  obs.file.path <- paste(part_dir, st.id, "/", sep = "")

  # Find files in subdirectory:
  obs.files.all <- list.files(path = obs.file.path, pattern = "*[:1:]_obs.txt")
  num.files <- length(obs.files.all)

  # Create postscript outfile:
  ps.file.name <- paste(st.id, "_partitioning_obs.ps", sep = "")

  postscript(
    file = paste(out_path, ps.file.name, sep = ""),
    width = 8,
    height = 4,
    paper = "special",
    horizontal = FALSE,
    onefile = TRUE
  )

  # Step through each file in subdirectory:
  for (obs.file in obs.files.all) {
    file_len <- length(count.fields(paste(obs.file.path, obs.file, sep = "")))
    if (file_len > 5) {
      plot_obs(obs.file.path, obs.file)
    }
  }

  # Close postscript file:
  dev.off()
}


# -------------------------------------------------- #
# ------------------ OUTLIERS ---------------------- #
# ------------------  linear  ---------------------- #
# -------------------------------------------------- #
# Step through each directory:
for (st.id in stations) {
  ro.file.path <- paste(part_dir, st.id, "/", sep = "")

  # Find files in subdirectory:
  ro.files.all = list.files(path = ro.file.path, pattern = "*_ro.txt")
  num.files <- length(ro.files.all)

  # Create postscript outfile:
  ps.file.name <- paste(st.id, "_partitioning_ro_l.ps", sep = "")
  postscript(
    file = paste(out_path, ps.file.name, sep = ""),
    width = 8,
    height = 4,
    paper = "special",
    horizontal = FALSE,
    onefile = TRUE
  )

  # Step through each file in subdirectory:
  for (ro.file in ro.files.all) {
    file_len <- length(count.fields(paste(ro.file.path, ro.file, sep = "")))
    if (file_len > 5) {
      plot_ro_l(ro.file.path, ro.file)
    }
  }

  # Close postscript file:
  dev.off()
}

# -------------------------------------------------- #
# ------------------ OUTLIERS ---------------------- #
# ------------------ hyperbol ---------------------- #
# -------------------------------------------------- #
# Step through each directory:
for (st.id in stations) {
  ro.file.path <- paste(part_dir, st.id, "/", sep = "")

  # Find files in subdirectory:
  ro.files.all = list.files(path = ro.file.path, pattern = "*_ro.txt")
  num.files <- length(ro.files.all)

  # Create postscript outfile:
  ps.file.name <- paste(st.id, "_partitioning_ro_h.ps", sep = "")
  postscript(
    file = paste(out_path, ps.file.name, sep = ""),
    width = 8,
    height = 4,
    paper = "special",
    horizontal = FALSE,
    onefile = TRUE
  )

  # Step through each file in subdirectory:
  for (ro.file in ro.files.all) {
    file_len <- length(count.fields(paste(ro.file.path, ro.file, sep = "")))
    if (file_len > 5) {
      plot_ro_h(ro.file.path, ro.file)
    }
  }

  # Close postscript file:
  dev.off()
}
