# RStudio 0.97
#
# plot_gapfill.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-03-13 -- created
# 2014-10-30 -- last updated
#
# ------------
# description:
# ------------
# This script plots the monthly PPFD observations and the gap-filling product
#
# Note: the gapfill data is assumed to be in station-specific sub-directories
# (see file_handler.pl)
#
# ----------
# changelog:
# ----------
# 00. created [14.03.13]
# 01. added get_max_ppfd() [14.03.18]
# --> incorporated max PPFD in plots
# 02. added num2month() [14.10.30]
#
# -----
# todo:
# -----
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### LIBRARIES ################################################################
# /////////////////////////////////////////////////////////////////////////////
library(chron)

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////

# ************************************************************************
# * Name: list.dirs
# *
# * Input: - char (path)
# *        - char (pattern)
# *        - boolean (all.dirs)
# *        - boolean (full.names)
# *        - boolean (ignore.case)
# *
# * Return: char vector
# *
# * Features: This reads all the directory names for a given path
# *
# *
# * Ref: J. Ulrich (2011) How to obtain a list of directories within a
# *        directory, like list.files(), but instead "list.dirs()", Stack
# *        Overflow, url: http://stackoverflow.com/questions/4749783/how-
# *        to-obtain-a-list-of-directories-within-a-directory-like-list-
# *        files-but-i
# ************************************************************************
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  all <- list.files(path, pattern, all.dirs,
                    full.names, recursive=FALSE, ignore.case)
  return(all)
}


# ************************************************************************
# * Name: get_max_ppfd
# *
# * Input: - char (my.dir)
# *        - char vector (my.all.files)
# *
# * Return: numeric vector
# *
# * Features: Returns pair of max values (actual and rounded up for 
# *           plotting)
# *
# ************************************************************************
get_max_ppfd <- function(my.dir, my.all.files){
  # This function reads through PPFD in all files and returns max value
  #
  # Initialize max PPFD:
  PPFD.max <- 0
  #
  # Read through all files:
  for (my.file in my.all.files){
    # Open file and read contents:
    contents <- read.csv(
      paste(my.dir,my.file,sep=""), 
      header=T, 
      na.strings="-9999.000"
    )
    #
    # Get max values:
    max.ppfd.gf <- max(contents$PPFDgf, na.rm=TRUE)
    max.ppfd.obs <- max(contents$PPFDobs, na.rm=TRUE)
    #
    # Do comparisons:
    if (max.ppfd.gf > PPFD.max){
      PPFD.max <- max.ppfd.gf
    }
    #
    if (max.ppfd.obs > PPFD.max){
      PPFD.max <- max.ppfd.obs
    }
  }
  #
  # Perform some minor rounding:
  if (PPFD.max < 1000){
    PPFD.ret <- signif(PPFD.max, 1) # round to the nearest hundred
  } else {
    PPFD.ret <- signif(PPFD.max, 2) # round to the nearest hundred
  }
  if (PPFD.ret < PPFD.max){
    PPFD.ret <- PPFD.ret + 100.0 # bump up to next 100
  }
  # Return max PPFD:
  c(PPFD.max, PPFD.ret)
}

# ************************************************************************
# * Name: plot.gpp
# *
# * Input: - char, file directory (my.dir)
# *        - char, file name (my.file)
# *        - double (y.max)
# *
# * Return: None.
# *
# * Depends: num2month
# *
# * Features: Plots observed and gap-filled PPFD
# *
# ************************************************************************
plot.gpp <- function(my.dir, my.file, y.max){
  # Read data from file:
  content <- read.csv(paste(my.dir,my.file,sep=""), 
                       header=T, 
                       na.strings="-9999.000")
  #
  # Grab month and year from file name:
  my.stat <- gsub("(^.*-.*)(-.*_.*\\.txt$)", "\\1", my.file)
  my.year <- gsub("(^.*_)(.*)(-.*-.*\\.txt$)", "\\2", my.file)
  my.month <- gsub("(^.*_)(.*-)(.*)(-.*\\.txt$)", "\\3", my.file)
  my.month <- num2month(my.month)
  #
  # Create time object:
  my.datetime <- do.call(
    rbind, sapply(as.character(content$Timestamp), strsplit, split=" ")
  )
  my.dates <- as.character(my.datetime[,1])
  my.times <- as.character(my.datetime[,2])
  my.tstamp <- chron(my.dates,my.times,format=c(date="y-m-d", time="h:m:s"))
  #
  par(mar=c(5,4.75,1,1.5))
  plot(
    my.tstamp,
    content$PPFDgf,
    ylab = expression(PPFD~(mu*mol%.%m^{-2}%.%s^{-1})),
    xlab = paste(my.stat," (", my.month, " ", my.year, ")", sep=""),
    ylim = c(0,y.max),
    type="l",
    col="red",
    lwd=1,
    lty=1
  )
  box(lwd=2)
  lines(
    my.tstamp,
    content$PPFDobs,
    lty=1,
    col="blue",
    lwd=2
  )
  legend(
    "top", 
    inset=0.01, 
    #y.intersp = 1.25, 
    col=c("blue","red"), 
    lty=c(1,1), 
    lwd=c(2,1), 
    c("Obs","GF"), 
    bg="white",
    box.col="white",
    horiz=TRUE
  )
}

# ************************************************************************
# * Name: num2month
# *
# * Input: char, month two-digit number (m)
# *
# * Return: char, three character month
# *
# * Features: Maps month number to character
# *
# * Ref: Sharpie (2009) "How to add variable key/value pair to list 
# *        object?," Stack Overflow, url: 
#*         http://stackoverflow.com/questions/1105659/how-to-add-variable-
#*         key-value-pair-to-list-object
# ************************************************************************
num2month <- function(m) {
  n = list()
  n[['01']] <- 'Jan'
  n[['02']] <- 'Feb'
  n[['03']] <- 'Mar'
  n[['04']] <- 'Apr'
  n[['05']] <- 'May'
  n[['06']] <- 'Jun'
  n[['07']] <- 'Jul'
  n[['08']] <- 'Aug'
  n[['09']] <- 'Sep'
  n[['10']] <- 'Oct'
  n[['11']] <- 'Nov'
  n[['12']] <- 'Dec'
  return(n[m])
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### DEFINITIONS ##############################################################
# /////////////////////////////////////////////////////////////////////////////
gf.file.path <- paste("/home/user/Projects/gepisat/results/2002-06/",
                      "attempt25/gapfill_v25/",
                      sep = "")

stations <- list.dirs(path=gf.file.path, pattern="*-*")

out.file.path <- paste("/home/user/Projects/gepisat/results/2002-06/",
                       "attempt25/gf_plots/", 
                       sep="")

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### MAIN #####################################################################
# /////////////////////////////////////////////////////////////////////////////
for (st.id in stations){
  # Assign directory for current station:
  file.dir <- paste(gf.file.path, st.id, "/", sep="")
  #
  # Read all files in directory:
  files.all = list.files(path = file.dir, pattern = "*[:1:].txt")
  num.files <- length(files.all)
  #
  # Find max PPFD:
  max.vals = get_max_ppfd(file.dir, files.all)
  #
  # Create postscript outfile:
  ps.file.name <- paste(toupper(st.id), "_gapfill.ps", sep="")
  postscript(file = paste(out.file.path, ps.file.name, sep=""), 
             width = 6, 
             height = 5, 
             paper = "special", 
             horizontal = FALSE,
             title = ps.file.name,
             onefile = TRUE)
  #
  # Step through each file in subdirectory:
  for (cur.file in files.all){
    plot.gpp(file.dir, cur.file, max.vals[2])
  }
  #
  # Close postscript file:
  dev.off()
}
