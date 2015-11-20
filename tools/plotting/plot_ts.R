# RStudio v 0.98.501
#
# plot_ts.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2014-03-20 -- created
# 2014-10-22 -- last updated
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script reads the data output from timeseries.py to plot the gridded
# data associated with flux tower locations.
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# 00. created [14.03.30]
# 01. added tower-specific monthly gpp [14.04.09]
#     -> there are now six plots, two columns by three rows
# 02. added climate and vegetation meta data [14.04.09]
# 03. added mean growing season data frame, processing, and writeout [14.05.08]
# 04. updated get.gpp and get.meta functions [14.06.05]
# 05. updated for latest (v.26) timeseries data [14.10.22]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### LIBRARIES ################################################################
# /////////////////////////////////////////////////////////////////////////////

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////
# ************************************************************************
# * Name: make.ts
# *
# * Input: matrix of time series data (my.data)
# *
# * Return: time series object
# *
# * Features: This function converts a matrix of time series data to a
# *           time series object (monthly data from 2002-2006)
# ************************************************************************
make.ts <- function(my.data){
  my.ts <- ts(
    my.data,
    frequency = 12,
    start = c(2002,1)
  )
  my.ts
}

# ************************************************************************
# * Name: get.gpp
# *
# * Input: flux tower name (my.site)
# *        GPP file directory (my.dir)
# *
# * Return: time series object
# *
# * Features: This function retrieves GPP data for a specific flux tower
# *           and returns a time series object (monthly GPP between 
# *           2002-2006)
# *
# * Notes:
# *     Methodology is based on code retrieved from:
# *     http://stackoverflow.com/questions/6058677/how-to-create-na-for-
# *     missing-data-in-a-time-series
# ************************************************************************
get.gpp <- function(my.site, my.dir){
  # Read GPP data from file:
  gpp.dir <- my.dir
  gpp.file <- list.files(path=gpp.dir, pattern="GePiSaT_All_Data*")
  gpp.data <- read.csv(
    paste(gpp.dir, gpp.file, sep=""),
    header=T
  )
  #
  # Create a full timeseries sequence:
  my.start <- as.Date("2002/01/01")
  my.full <- seq(my.start, by='1 month', length=60)
  #
  # Find rows for specific site:
  gpp.rows <- which(gpp.data$Station == my.site)
  #
  if (length(gpp.rows) == 0){
    # Flux tower does not have any GPP data, fill with zeros
    site.data <- rep(0, 60)
    site.dates <- my.full
  } else {
    # Retrieve site-specific monthly GPP and associated timestamps:
    site.data <- gpp.data$GPP.mol_m2[gpp.rows]
    site.dates <- as.Date(as.character(gpp.data$Timestamp[gpp.rows]))
  }
  #
  # Make site data into a data frame:
  site.gpp <- data.frame(
    date=site.dates,
    gpp=site.data
  )
  #
  # Fill in station GPP to full time series:
  full.gpp <- data.frame(
    Date=my.full,
    GPP=with(site.gpp, gpp[match(my.full, date)])
  )
  #
  # Return data frame with *gappy* GPP:
  full.gpp
}

# ************************************************************************
# * Name: get.meta
# *
# * Input: flux tower name (my.site)
# *        meta file directory (my.dir)
# *
# * Return: data frame
# *
# * Features: This function retrieves meta data (vegetation & climate) 
# *           for a specific flux tower and returns it in a data frame
# ************************************************************************
get.meta <- function(my.site, my.dir){
  # Read meta data for specific flux tower
  meta.dir <- my.dir
  meta.file <- list.files(path=meta.dir, pattern="Fluxdata_Meta-Data*")
  meta.data <- read.csv(
    paste(meta.dir, meta.file, sep=""),
    header=T
  )
  #
  # Find row for flux tower:
  meta.row <- which(meta.data$stationid == my.site)
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
# * Name: get.mean.growing
# *
# * Input: data content object (my.obj)
# *        mean growing season object (gs.obj)
# *
# * Return: mean growing season object (gs.obj)
# *
# * Features: This function fills an empty mean growing season object
# *           with mean growing season temperature and VPD.  Mean 
# *           growing season is assumed to be any month with mean 
# *           air temperature greater than 5 deg. C
# ************************************************************************
get.mean.growing <- function(my.obj, gs.obj){
  # Check for the right amount of months (5 years, 60 months):
  if (dim(my.obj)[1] == 60){
    # Get station name:
    st.name <- as.character(my.obj$station[1])
    #
    # Pull years from timestamp and cast them to numeric:
    all.years <- as.numeric(substr(as.character(my.obj$timestamp), 1, 4))
    #
    # Iterate through each year:
    for(yr in seq(from=2002, to=2006, by=1)){
      yr.rows <- which(all.years == yr)
      #
      # Get growing season temperatures (Tc > 5 deg C):
      tc.col <- paste('tc', yr, sep='.')
      gs.tc <- my.obj$tc[yr.rows][which(my.obj$tc[yr.rows] > 5)]
      if (length(gs.tc) > 0){
        gs.obj[st.name, tc.col] <- mean(gs.tc)
      } else {
        gs.obj[st.name, tc.col] <- NA
      }
      #
      # Get growing season VPD (Tc > 5 deg C):
      vpd.col <- paste('vpd', yr, sep='.')
      gs.vpd <- my.obj$vpd[yr.rows][which(my.obj$tc[yr.rows] > 5)]
      if (length(gs.vpd) > 0){
        gs.obj[st.name, vpd.col] <- mean(gs.vpd)
      } else {
        gs.obj[st.name, vpd.col] <- NA
      }
    }
  }
  # Return:
  gs.obj
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### CONSTANTS ################################################################
# /////////////////////////////////////////////////////////////////////////////
file.path <- paste("/home/user/Projects/gepisat/results/2002-06/",
                   "attempt26/", sep="")

gpp.path <- paste("/home/user/Projects/gepisat/results/2002-06/",
                  "attempt26/", sep="")

meta.path <- paste("/home/user/Projects/gepisat/",
                   "data/fluxdata.org/", sep="")

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### MAIN #####################################################################
# /////////////////////////////////////////////////////////////////////////////
my.file <- list.files(path=file.path, pattern = "^timeseries.*txt$")
content <- read.csv(paste(file.path, my.file, sep=""), 
                    header=T,
                    na.strings="-9999.000000")

# Save station names:
my.stations <- levels(content$Station)

# Create a data frame for storing mean growing season temperature & vpd:
dummy.vect <- 0.0*seq(from=1, to=length(my.stations), by=1)
growing.season <- matrix(
  cbind(
    dummy.vect,  # 2002 Tc
    dummy.vect,  # 2002 VPD
    dummy.vect,  # 2003 Tc
    dummy.vect,  # 2003 VPD
    dummy.vect,  # 2004 Tc
    dummy.vect,  # 2004 VPD
    dummy.vect,  # 2005 Tc
    dummy.vect,  # 2005 VPD
    dummy.vect,  # 2006 Tc
    dummy.vect   # 2006 VPD
  ),
  ncol=10,
  nrow=length(my.stations)
)
colnames(growing.season) <- c(
  "tc.2002","vpd.2002",
  "tc.2003","vpd.2003",
  "tc.2004","vpd.2004",
  "tc.2005","vpd.2005",
  "tc.2006","vpd.2006"
)
growing.season <- as.data.frame(
  growing.season,
  row.names = my.stations
)

# Postscript:
ps.out <- paste("GePiSaT","time-series.ps", sep="_")
postscript(
  file = paste(file.path,ps.out, sep=""), 
  title = "GePiSaT Time Series",
  width = 6, 
  height = 4, 
  paper = "special", 
  horizontal = FALSE,
  onefile = TRUE
)

for (station in my.stations){
  # Save row numbers for current station:
  my.rows <- which(content$Station == station)
  #
  # Get climate and class from flux tower meta data:
  my.meta <- get.meta(station, meta.path)
  #
  # Retrieve station GPP data:
  #my.gpp <- get.gpp(station, gpp.path)
  #
  # Get growing season means:
  #growing.season <- get.mean.growing(content[my.rows,], growing.season)
  #
  my.data <- matrix(
    cbind(
      content$Tair_C[my.rows],
      content$VPD_kPa[my.rows],
      content$alpha[my.rows],
      content$fAPAR[my.rows],
      content$SF[my.rows],
      content$Pre_mm[my.rows]
    ),
    ncol = 6,
    nrow = 60
  )
  colnames(my.data) <- c(
    "Ta","VPD, kPa","CPA","fPAR","SF","PPT, mm"
  )
  my.ts <- make.ts(my.data)
  #
  plot(my.ts,
       main=station,
       xlab="",
       cex.lab=0.8,
       cex.main=1,
       yax.flip=T,
       col="blue"
  )
  #
  mtext(paste("Climate ID: ", as.character(my.meta$climate), "\n",
              "Class ID: ", as.character(my.meta$class)),
        side=3, line=1, adj=1, cex=0.8)
  #
  # Delete unnecessary variables:
  rm(my.rows, my.data, my.meta, my.ts)
}
dev.off()

# Write growing season results to file:
out.file <- paste(file.path, "GePiSaT_mean_growing_season.txt", sep="")
write.table(
  growing.season, 
  file=out.file,
  sep=",",
  row.names=TRUE,
  col.names=TRUE
)
