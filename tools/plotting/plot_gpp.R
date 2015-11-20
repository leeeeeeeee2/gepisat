# RStudio v. 0.97.551
#
# plot_gpp.R
#
# written by Tyler W. Davis
# Imperial College London
#
# 2013-10-30 -- created
# 2013-11-08 -- last upated
#
# ------------
# description:
# ------------
# This script creates plots of GPP based on the LUE station text files.
# 1. Three climate regions / forest regions for each climate region
# 2. Box and whisker plot of monthly GPP for each region
# 3. Calculates monthly averages of GPP +/- st deviation
# 4. Writes out results
#
# ----------
# changelog:
# ----------
# 01. Added box and whisker plot function [13.10.31]
# 02. Created process_gpp function [13.11.08]
# 03. Updated monthly stats [13.11.08]
# --> zero fill months with no data
# --> added monthly max and min
# 04. Added "name" field to list objects [13.11.08]
#
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### FUNCTIONS ################################################################
# /////////////////////////////////////////////////////////////////////////////
monthly_append <- function(climate.object, data){
  # This function appends contents to lists in the climate object
  # Add monthly values to lists
  jan <- data$GPP[which(data$Month == 1)]
  if(length(jan) != 0){
    climate.object$jan <- c(climate.object$jan, jan)
  }
  feb <- data$GPP[which(data$Month == 2)]
  if (length(feb) != 0){
    climate.object$feb <- c(climate.object$feb, feb)
  }
  mar <- data$GPP[which(data$Month == 3)]
  if (length(mar) != 0){
    climate.object$mar <- c(climate.object$mar, mar)
  }
  apr <- data$GPP[which(data$Month == 4)]
  if (length(apr) != 0){
    climate.object$apr <- c(climate.object$apr, apr)
  }
  may <- data$GPP[which(data$Month == 5)]
  if (length(may) != 0){
    climate.object$may <- c(climate.object$may, may)
  }
  jun <- data$GPP[which(data$Month == 6)]
  if (length(jun) != 0){
    climate.object$jun <- c(climate.object$jun, jun)
  }
  jul <- data$GPP[which(data$Month == 7)]
  if (length(jul) != 0){
    climate.object$jul <- c(climate.object$jul, jul)
  }
  aug <- data$GPP[which(data$Month == 8)]
  if (length(aug) != 0){
    climate.object$aug <- c(climate.object$aug, aug)
  }
  sep <- data$GPP[which(data$Month == 9)]
  if (length(sep) != 0){
    climate.object$sep <- c(climate.object$sep, sep)
  }
  oct <- data$GPP[which(data$Month == 10)]
  if (length(oct) != 0){
    climate.object$oct <- c(climate.object$oct, oct)
  }
  nov <- data$GPP[which(data$Month == 11)]
  if (length(nov) != 0){
    climate.object$nov <- c(climate.object$nov, nov)
  }
  dec <- data$GPP[which(data$Month == 12)]
  if (length(dec) != 0){
    climate.object$dec <- c(climate.object$dec, dec)
  }
  climate.object
}

monthly_stats <- function(climate.object){
  # This function appends two new lists for mean/sd for 12 months
  if (length(climate.object$jan) != 0){
    jan.ave <- mean(climate.object$jan)
    jan.std <- sd(climate.object$jan)
    jan.max <- max(climate.object$jan)
    jan.min <- min(climate.object$jan)
  } else {
    climate.object$jan <- 0
    jan.ave <- 0
    jan.std <- 0
    jan.max <- 0
    jan.min <- 0
  }
  if (length(climate.object$feb) != 0){
    feb.ave <- mean(climate.object$feb)
    feb.std <- sd(climate.object$feb)
    feb.max <- max(climate.object$feb)
    feb.min <- min(climate.object$feb)
  } else {
    climate.object$feb <- 0
    feb.ave <- 0
    feb.std <- 0
    feb.max <- 0
    feb.min <- 0
  }
  if (length(climate.object$mar) != 0){
    mar.ave <- mean(climate.object$mar)
    mar.std <- sd(climate.object$mar)
    mar.max <- max(climate.object$mar)
    mar.min <- min(climate.object$mar)
  } else {
    climate.object$mar <- 0
    mar.ave <- 0
    mar.std <- 0
    mar.max <- 0
    mar.min <- 0
  }
  if (length(climate.object$apr) != 0){
    apr.ave <- mean(climate.object$apr)
    apr.std <- sd(climate.object$apr)
    apr.max <- max(climate.object$apr)
    apr.min <- min(climate.object$apr)
  } else {
    climate.object$apr <- 0
    apr.ave <- 0
    apr.std <- 0
    apr.max <- 0
    apr.min <- 0
  }
  if (length(climate.object$may) != 0){
    may.ave <- mean(climate.object$may)
    may.std <- sd(climate.object$may)
    may.max <- max(climate.object$may)
    may.min <- min(climate.object$may)
  } else {
    climate.object$may <- 0
    may.ave <- 0
    may.std <- 0
    may.max <- 0
    may.min <- 0
  }
  if (length(climate.object$jun) != 0){
    jun.ave <- mean(climate.object$jun)
    jun.std <- sd(climate.object$jun)
    jun.max <- max(climate.object$jun)
    jun.min <- min(climate.object$jun)
  } else {
    climate.object$jun <- 0
    jun.ave <- 0
    jun.std <- 0
    jun.max <- 0
    jun.min <- 0
  }
  if (length(climate.object$jul) != 0){
    jul.ave <- mean(climate.object$jul)
    jul.std <- sd(climate.object$jul)
    jul.max <- max(climate.object$jul)
    jul.min <- min(climate.object$jul)
  } else {
    climate.object$jul <- 0
    jul.ave <- 0
    jul.std <- 0
    jul.max <- 0
    jul.min <- 0
  }
  if (length(climate.object$aug) != 0){
    aug.ave <- mean(climate.object$aug)
    aug.std <- sd(climate.object$aug)
    aug.max <- max(climate.object$aug)
    aug.min <- min(climate.object$aug)
  } else {
    climate.object$aug <- 0
    aug.ave <- 0
    aug.std <- 0
    aug.max <- 0
    aug.min <- 0
  }
  if (length(climate.object$sep) != 0){
    sep.ave <- mean(climate.object$sep)
    sep.std <- sd(climate.object$sep)
    sep.max <- max(climate.object$sep)
    sep.min <- min(climate.object$sep)
  } else {
    climate.object$sep <- 0
    sep.ave <- 0
    sep.std <- 0
    sep.max <- 0
    sep.min <- 0
  }
  if (length(climate.object$oct) != 0){
    oct.ave <- mean(climate.object$oct)
    oct.std <- sd(climate.object$oct)
    oct.max <- max(climate.object$oct)
    oct.min <- min(climate.object$oct)
  } else {
    climate.object$oct <- 0
    oct.ave <- 0
    oct.std <- 0
    oct.max <- 0
    oct.min <- 0
  }
  if (length(climate.object$nov) != 0){
    nov.ave <- mean(climate.object$nov)
    nov.std <- sd(climate.object$nov)
    nov.max <- max(climate.object$nov)
    nov.min <- min(climate.object$nov)
  } else {
    climate.object$nov <- 0
    nov.ave <- 0
    nov.std <- 0
    nov.max <- 0
    nov.min <- 0
  }
  if (length(climate.object$dec) != 0){
    dec.ave <- mean(climate.object$dec)
    dec.std <- sd(climate.object$dec)
    dec.max <- max(climate.object$dec)
    dec.min <- min(climate.object$dec)
  } else {
    climate.object$dec <- 0
    dec.ave <- 0
    dec.std <- 0
    dec.max <- 0
    dec.min <- 0
  }
  climate.object$mean <- c(
    jan.ave, feb.ave, mar.ave, apr.ave,
    may.ave, jun.ave, jul.ave, aug.ave,
    sep.ave, oct.ave, nov.ave, dec.ave
    )
  climate.object$sdev <- c(
    jan.std, feb.std, mar.std, apr.std,
    may.std, jun.std, jul.std, aug.std,
    sep.std, oct.std, nov.std, dec.std
    )
  climate.object$max <- c(
    jan.max, feb.max, mar.max, apr.max,
    may.max, jun.max, jul.max, aug.max,
    sep.max, oct.max, nov.max, dec.max
    )
  climate.object$min <- c(
    jan.min, feb.min, mar.min, apr.min,
    may.min, jun.min, jul.min, aug.min,
    sep.min, oct.min, nov.min, dec.min
    )
  climate.object
}

box_and_whisker <- function(climate.object, write.out, file.dir){
  # This function takes the climate.object and converts it to a data frame
  # and produces a box plot.
  # Length of each month:
  jan.N <- length(climate.object$jan)
  feb.N <- length(climate.object$feb)
  mar.N <- length(climate.object$mar)
  apr.N <- length(climate.object$apr)
  may.N <- length(climate.object$may)
  jun.N <- length(climate.object$jun)
  jul.N <- length(climate.object$jul)
  aug.N <- length(climate.object$aug)
  sep.N <- length(climate.object$sep)
  oct.N <- length(climate.object$oct)
  nov.N <- length(climate.object$nov)
  dec.N <- length(climate.object$dec)
  # 
  # Lists of month names:
  jan.L <- rep(1, jan.N)
  feb.L <- rep(2, feb.N)
  mar.L <- rep(3, mar.N)
  apr.L <- rep(4, apr.N)
  may.L <- rep(5, may.N)
  jun.L <- rep(6, jun.N)
  jul.L <- rep(7, jul.N)
  aug.L <- rep(8, aug.N)
  sep.L <- rep(9, sep.N)
  oct.L <- rep(10, oct.N)
  nov.L <- rep(11, nov.N)
  dec.L <- rep(12, dec.N)
  #
  # Month List:
  months <- c(
    jan.L, feb.L, mar.L, apr.L,
    may.L, jun.L, jul.L, aug.L,
    sep.L, oct.L, nov.L, dec.L)
  #
  # GPP list:
  gpp <- c(
    climate.object$jan, climate.object$feb, climate.object$mar,
    climate.object$apr, climate.object$may, climate.object$jun,
    climate.object$jul, climate.object$aug, climate.object$sep,
    climate.object$oct, climate.object$nov, climate.object$dec)
  #
  # Create dataframe:
  headers <- c("Month", "GPP")
  my.df <- as.data.frame(cbind(months, gpp))
  names(my.df) <- headers
  my.df$Month = factor(my.df$Month, labels = c(
    "Jan", "Feb", "Mar", "Apr",
    "May", "Jun", "Jul", "Aug",
    "Sep", "Oct", "Nov", "Dec"))
  if (write.out == 1){
    postscript(
      file = paste(file.dir, climate.object$name, ".ps", sep=""), 
      width = 8,
      height = 6,
      paper = "special", 
      horizontal = FALSE,
      onefile = TRUE
    )
    par(mar=c(5.5,5.25,1,1))
    boxplot(
      GPP ~ Month,
      data = my.df,
      main = "",
      xlab = climate.object$name,
      ylab = expression(paste("GPP, ", mol~CO[2]%.%m^{-2}))
      )
    box(lwd=2)
    dev.off()
  } else {
    par(mar=c(4.25,4.75,1,1))
    boxplot(
      GPP ~ Month,
      data = my.df,
      main = "",
      xlab = climate.object$name,
      ylab = expression(paste("GPP, ", mol~CO[2]%.%m^{-2})),
      varwidth=TRUE,
      col="lightgreen"
      )
    box(lwd=2)
  }
}

bp_annual <- function(file.dir){
  annual.gpp <- read.table(
    paste(file.dir, "Annual_GPP-All_Stations.txt", sep=""),
    header = TRUE,
    na.strings = "0"
    )
  #
  # Retrieve annual GPP lists:
  gpp.2002 <- annual.gpp$X2002[!is.na(annual.gpp$X2002)]
  gpp.2003 <- annual.gpp$X2003[!is.na(annual.gpp$X2003)]
  gpp.2004 <- annual.gpp$X2004[!is.na(annual.gpp$X2004)]
  gpp.2005 <- annual.gpp$X2005[!is.na(annual.gpp$X2005)]
  gpp.2006 <- annual.gpp$X2006[!is.na(annual.gpp$X2006)]
  #
  # Lists of year names:
  gpp.2002.L <- rep(2002, length(gpp.2002))
  gpp.2003.L <- rep(2003, length(gpp.2003))
  gpp.2004.L <- rep(2004, length(gpp.2004))
  gpp.2005.L <- rep(2005, length(gpp.2005))
  gpp.2006.L <- rep(2006, length(gpp.2006))
  #
  # Create year and GPP lists:
  year.list <- c(gpp.2002.L, gpp.2003.L, gpp.2004.L, gpp.2005.L, gpp.2006.L)
  gpp.list <- c(gpp.2002, gpp.2003, gpp.2004, gpp.2005, gpp.2006)
  #
  # Create dataframe:
  headers <- c("Year", "GPP")
  my.df <- as.data.frame(cbind(year.list, gpp.list))
  names(my.df) <- headers
  my.df$Year = factor(
    my.df$Year, 
    labels = c("2002", "2003", "2004", "2005", "2006")
    )
  #
  # Box plot:
  par(mar=c(3.75,4.75,1,1))
  boxplot(
    GPP ~ Year,
    data = my.df,
    main = "",
    xlab = "",
    ylab = expression(paste("GPP, ", mol~CO[2]%.%m^{-2})),
    varwidth=TRUE
  )
  box(lwd=2)
}

process_gpp <- function(file.dir, climate.object){
  # Define all years:
  all.years <- c("2002", "2003", "2004", "2005", "2006")
  #
  # Initialize zero place holders:
  zero.vals <- rep(0, length(climate.object$stations))
  climate.gpp <- as.data.frame(
    cbind(zero.vals, zero.vals, zero.vals, zero.vals, zero.vals),
    row.names = climate.object$stations)
  names(climate.gpp) <- all.years
  #
  # Open each file for processing
  for (station in climate.object$stations){
    # Open file and read contents:
    my.file <- paste(file.dir, station, "_LUE.txt", sep="")
    content <- read.csv(my.file, header=T)
    #
    # Add month column:
    content$Month <- as.numeric(substr(content$Timestamp, 6, 7))
    content$Year <- as.numeric(substr(content$Timestamp, 1, 4))
    #
    # Unique years:
    years <- unique(content$Year)
    for (yr in years){
      # Reset annual accumulated GPP:
      gpp.annual <- 0
      months <- length(content$GPP[which(content$Year == yr)])
      if (months >= 10){
        gpp.ave <- mean(content$GPP[which(content$Year == yr)])
        gpp.ann <- gpp.ave * 12
        climate.gpp[station, as.character(yr)] <- gpp.ann
      }
    }
    #
    # Update monthly GPP:
    climate.object <- monthly_append(climate.object, content)
    rm(content)
  }
  #
  # Write monthly GPP results to file:
  write.table(
    climate.gpp,
    file = paste(file.dir, "Annual_GPP-", climate.object$name, ".txt", sep=""),
    row.names = TRUE,
    col.names = TRUE)
  #
  # Calculate monthly stats:
  climate.object <- monthly_stats(climate.object)
  climate.object
}

plot_annual <- function(file.dir){
  # Find files:
  arid.file <- paste(file.dir, "Annual_GPP-Arid_All.txt", sep="")
  boreal.file <- paste(file.dir, "Annual_GPP-Boreal_All.txt", sep="")
  subtropical.file <- paste(file.dir, "Annual_GPP-Subtropical_All.txt", sep="")
  temperate.file <- paste(file.dir, "Annual_GPP-Temperate_All.txt", sep="")
  tropical.file <- paste(file.dir, "Annual_GPP-Tropical_All.txt", sep="")
  #
  # Open contents:
  arid <- read.table(arid.file, header=T, na.strings="0")
  boreal <- read.table(boreal.file, header=T, na.strings="0")
  subtropical <- read.table(subtropical.file, header=T, na.strings="0")
  temperate <- read.table(temperate.file, header=T, na.strings="0")
  tropical <- read.table(tropical.file, header=T, na.strings="0")
  #
  # Find max GPP (for plotting boundary)
  gpp.max <- 0
  gpp.max <- max(
    max(arid), max(boreal), max(subtropical), max(temperate), max(tropical) 
    )
  #
  # Length of each climate class:
  arid.n <- length(arid[,1])
  boreal.n <- length(boreal[,1])
  subtropical.n <- length(subtropical[,1])
  temperate.n <- length(temperate[,1])
  tropical.n <- length(tropical[,1])
  #
  # Plot:
  par(mar=c(2.5, 4.75, 1, 1))
  plot(
    rep(2002, temperate.n), 
    temperate$X2002,
    xlim=c(2001, 2006),
    ylim=c(0, gpp.max),
    pch = 20,
    col="forestgreen",
    main = "",
    xlab = "",
    ylab = expression(paste("GPP (", mol%.%m^{-2},")")),
    axes = FALSE
    )
  axis(side=1, lwd=2)
  axis(side=2, lwd=2)
  points(
    rep(2002, temperate.n),
    temperate$X2002,
    pch = 1,
    col = "black"
    )
  points(
    rep(2003, temperate.n),
    temperate$X2003,
    pch = 20,
    col = "forestgreen"
    )
  points(
    rep(2003, temperate.n),
    temperate$X2003,
    pch = 1,
    col = "black"
    )
  points(
    rep(2004, temperate.n),
    temperate$X2004,
    pch = 20,
    col = "forestgreen"
    )
  points(
    rep(2004, temperate.n),
    temperate$X2004,
    pch = 1,
    col = "black"
    )
  points(
    rep(2005, temperate.n),
    temperate$X2005,
    pch = 20,
    col = "forestgreen"
    )
  points(
    rep(2005, temperate.n),
    temperate$X2005,
    pch = 1,
    col = "black"
    )
  points(
    rep(2006, temperate.n),
    temperate$X2006,
    pch = 20,
    col = "forestgreen"
    )
  points(
    rep(2006, temperate.n),
    temperate$X2006,
    pch = 1,
    col = "black"
  )
  #
  # Subtropical
  points(
    rep(2002, subtropical.n),
    subtropical$X2002,
    pch = 20,
    col = "gold"
    )
  points(
    rep(2002, subtropical.n),
    subtropical$X2002,
    pch = 1,
    col = "black"
    )
  points(
    rep(2003, subtropical.n),
    subtropical$X2003,
    pch = 20,
    col = "gold"
    )
  points(
    rep(2003, subtropical.n),
    subtropical$X2003,
    pch = 1,
    col = "black"
    )
  points(
    rep(2004, subtropical.n),
    subtropical$X2004,
    pch = 20,
    col = "gold"
  )
  points(
    rep(2004, subtropical.n),
    subtropical$X2004,
    pch = 1,
    col = "black"
  )
  points(
    rep(2005, subtropical.n),
    subtropical$X2005,
    pch = 20,
    col = "gold"
  )
  points(
    rep(2005, subtropical.n),
    subtropical$X2005,
    pch = 1,
    col = "black"
  )
  points(
    rep(2006, subtropical.n),
    subtropical$X2006,
    pch = 20,
    col = "gold"
  )
  points(
    rep(2006, subtropical.n),
    subtropical$X2006,
    pch = 1,
    col = "black"
  )
  #
  # Boreal
}

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### DEFINITIONS ##############################################################
# /////////////////////////////////////////////////////////////////////////////
# All stations:
all.stations <- c()
all.stations$name <- "All_Stations"
all.stations$stations <- c(
  "AT-Neu", 
  "AU-Fog", 
  "AU-How", 
  "AU-Tum", 
  "AU-Wac", 
  "BE-Bra", 
  "BE-Jal", 
  "BE-Lon", 
  "BE-Vie", 
  "BR-Sa3", 
  "BW-Ghg", 
  "BW-Ghm", 
  "CA-Man", 
  "CA-Mer", 
  "CA-NS1", 
  "CA-NS2", 
  "CA-NS3", 
  "CA-NS4", 
  "CA-NS5", 
  "CA-NS6", 
  "CA-NS7", 
  "CA-Qcu", 
  "CA-Qfo", 
  "CA-SF1", 
  "CA-SF2", 
  "CA-SF3", 
  "CH-Oe1", 
  "CH-Oe2", 
  "CZ-BK1", 
  "CZ-wet", 
  "DE-Geb", 
  "DE-Gri", 
  "DE-Hai", 
  "DE-Kli", 
  "DE-Meh", 
  "DE-Tha", 
  "DE-Wet", 
  "DK-Fou", 
  "DK-Lva", 
  "DK-Ris", 
  "DK-Sor", 
  "ES-ES1", 
  "ES-ES2", 
  "ES-LMa", 
  "ES-VDA", 
  "FI-Hyy", 
  "FI-Kaa", 
  "FI-Sod", 
  "FR-Fon", 
  "FR-Gri", 
  "FR-Hes", 
  "FR-LBr", 
  "FR-Lq1", 
  "FR-Lq2", 
  "FR-Pue", 
  "HU-Bug", 
  "HU-Mat", 
  "ID-Pag", 
  "IE-Ca1", 
  "IE-Dri", 
  "IL-Yat", 
  "IT-Amp", 
  "IT-BCi", 
  "IT-Cas", 
  "IT-Col", 
  "IT-Cpz", 
  "IT-LMa", 
  "IT-Lav", 
  "IT-Lec", 
  "IT-MBo", 
  "IT-Mal", 
  "IT-Non", 
  "IT-PT1", 
  "IT-Pia", 
  "IT-Ren", 
  "IT-Ro1", 
  "IT-Ro2", 
  "IT-SRo", 
  "NL-Ca1", 
  "NL-Haa", 
  "NL-Hor", 
  "NL-Lan", 
  "NL-Loo", 
  "NL-Lut", 
  "NL-Mol", 
  "PL-wet", 
  "PT-Esp", 
  "PT-Mi1", 
  "PT-Mi2", 
  "RU-Cok", 
  "RU-Fyo", 
  "RU-Ha1", 
  "RU-Ha2", 
  "RU-Ha3", 
  "RU-Zot", 
  "SE-Deg", 
  "SE-Faj", 
  "SE-Fla", 
  "SE-Nor", 
  "SE-Sk1", 
  "SE-Sk2", 
  "SK-Tat", 
  "UK-AMo", 
  "UK-EBu", 
  "UK-ESa", 
  "UK-Gri", 
  "UK-Ham", 
  "UK-Her", 
  "UK-PL3", 
  "US-ARM", 
  "US-Aud", 
  "US-Bar", 
  "US-Bkg", 
  "US-Blo", 
  "US-Bo1", 
  "US-FPe", 
  "US-Goo", 
  "US-Ha1", 
  "US-Ho1", 
  "US-Ho2", 
  "US-Los", 
  "US-MMS", 
  "US-MOz", 
  "US-Ne1", 
  "US-Ne2", 
  "US-Ne3", 
  "US-Oho", 
  "US-PFa", 
  "US-SP1", 
  "US-SP2", 
  "US-SP3", 
  "US-Syv", 
  "US-Ton", 
  "US-UMB", 
  "US-Var", 
  "US-WCr", 
  "US-Wi0", 
  "US-Wi1", 
  "US-Wi2", 
  "US-Wi4", 
  "US-Wi5", 
  "US-Wi6", 
  "US-Wi7", 
  "US-Wi8", 
  "US-Wi9", 
  "ZA-Kru"
  )

# Station IDs for boreal, temperate, and tropical forest:
boreal.enf <- c()
boreal.enf$name <- "Boreal_ENF"
boreal.enf$stations <- c(
  "CA-Man",
  "CA-NS1",
  "CA-NS2",
  "CA-NS3",
  "CA-NS4",
  "CA-NS5",
  "CA-Qcu",
  "CA-Qfo",
  "CA-SF1",
  "CA-SF2",
  "FI-Hyy",
  "FI-Sod",
  "RU-Zot",
  "SE-Fla")
temperate.dbf <- c()
temperate.dbf$name <- "Temperate_DBF"
temperate.dbf$stations <- c(
  "DK-Sor",
  "FR-Fon",
  "FR-Hes",
  "DE-Hai",
  "UK-Ham",
  "UK-PL3")
temperate.enf <- c()
temperate.enf$name <- "Temperate_ENF"
temperate.enf$stations <- c(
  "AU-Tum",
  "AU-Wac",
  "FR-LBr",
  "DE-Tha",
  "DE-Wet",
  "IT-Lav",
  "IT-Ren",
  "NL-Loo",
  "UK-Gri")
temperate.mf <- c()
temperate.mf$name <- "Temperate_MF"
temperate.mf$stations <- c(
  "BE-Jal",
  "BE-Vie")
tropical.ebf <- c()
tropical.ebf$name <- "Tropical_EBF"
tropical.ebf$stations <- c(
  "BR-Sa3",
  "ID-Pag")
temperate.cro <- c()
temperate.cro$name <- "Temperate_CRO"
temperate.cro$stations <- c(
  "BE-Lon",
  "DK-Fou",
  "DK-Ris",
  "FR-Gri",
  "DE-Geb",
  "DE-Kli",
  "IE-Ca1",
  "NL-Lan",
  "NL-Lut",
  "NL-Mol",
  "CH-Oe2",
  "UK-ESa",
  "UK-Her"
  )
subtropical.cro <- c()
subtropical.cro$name <- "Subtropical_CRO"
subtropical.cro$stations <- c(
  "IT-Cas",
  "US-ARM",
  "IT-BCi",
  "ES-ES2"
  )
temperate.csh <- c()
temperate.csh$name <- "Temperate_CSH"
temperate.csh$stations <- c("US-Los")
subtropical.dbf <- c()
subtropical.dbf$name <- "Subtropical_DBF"
subtropical.dbf$stations <- c(
  "IT-Col",
  "IT-Non",
  "IT-PT1",
  "US-MMS",
  "US-MOz",
  "IT-Ro1",
  "IT-Ro2"
  )
subtropical.ebf <- c()
subtropical.ebf$name <- "Subtropical_EBF"
subtropical.ebf$stations <- c(
  "IT-Lec",
  "FR-Pue",
  "IT-Cpz",
  "PT-Esp",
  "PT-Mi1"
  )
subtropical.enf <- c()
subtropical.enf$name <- "Subtropical_ENF"
subtropical.enf$stations <- c(
  "US-SP1",
  "US-SP2",
  "US-SP3",
  "IT-SRo",
  "ES-ES1",
  "US-Blo"
  )
boreal.gra <- c()
boreal.gra$name <- "Boreal_GRA"
boreal.gra$stations <- c(
  "RU-Ha1",
  "RU-Ha2",
  "RU-Ha3"
  )
arid.gra <- c()
arid.gra$name <- "Arid_GRA"
arid.gra$stations <- c(
  "US-Aud",
  "US-FPe"
  )
subtropical.gra <- c()
subtropical.gra$name <- "Subtropical_GRA"
subtropical.gra$stations <- c(
  "US-Goo",
  "PT-Mi2",
  "US-Var"
  )
temperate.gra <- c()
temperate.gra$name <- "Temperate_GRA"
temperate.gra$stations <- c(
  "DK-Lva",
  "FR-Lq1",
  "FR-Lq2",
  "DE-Gri",
  "DE-Meh",
  "HU-Bug",
  "HU-Mat",
  "IE-Dri",
  "IT-LMa",
  "IT-Mal",
  "IT-MBo",
  "NL-Ca1",
  "NL-Haa",
  "NL-Hor",
  "ES-VDA",
  "CH-Oe1",
  "UK-EBu",
  "US-Bkg"
  )
boreal.osh <- c()
boreal.osh$name <- "Boreal_OSH"
boreal.osh$stations <- c(
  "CA-NS6",
  "CA-NS7",
  "CA-SF3",
  "RU-Cok"
  )
boreal.all <- c()
boreal.all$name <- "Boreal_All"
boreal.all$stations <- c(
  "CA-Man",
  "CA-NS1",
  "CA-NS2",
  "CA-NS3",
  "CA-NS4",
  "CA-NS5",
  "CA-Qcu",
  "CA-Qfo",
  "CA-SF1",
  "CA-SF2",
  "FI-Hyy",
  "FI-Sod",
  "RU-Zot",
  "SE-Fla",
  "RU-Ha1",
  "RU-Ha2",
  "RU-Ha3",
  "CA-NS6",
  "CA-NS7",
  "CA-SF3",
  "RU-Cok",
  "FI-Kaa",
  "SE-Deg"
  )

arid.all <- c()
arid.all$name <- "Arid_All"
arid.all$stations <- c(
  "US-Aud",
  "US-FPe",
  "BW-Ghg",
  "BW-Ghm"
  )

subtropical.all <- c()
subtropical.all$name <- "Subtropical_All"
subtropical.all$stations <- c(
  "US-ARM",
  "IT-Amp",
  "ZA-Kru",
  "IL-Yat",
  "IT-Cas",
  "IT-BCi",
  "ES-ES2",
  "IT-Col",
  "IT-Non",
  "IT-PT1",
  "US-MMS",
  "US-MOz",
  "IT-Ro1",
  "IT-Ro2",
  "IT-Lec",
  "FR-Pue",
  "IT-Cpz",
  "PT-Esp",
  "PT-Mi1",
  "US-SP1",
  "US-SP2",
  "US-SP3",
  "IT-SRo",
  "ES-ES1",
  "US-Blo",
  "US-Goo",
  "PT-Mi2",
  "US-Var",
  "IT-Pia",
  "ES-LMa",
  "US-Ton"
  )

temperate.all <- c()
temperate.all$name <- "Temperate_All"
temperate.all$stations <- c(
  "BE-Lon",
  "DK-Fou",
  "DK-Ris",
  "FR-Gri",
  "DE-Geb",
  "DE-Kli",
  "IE-Ca1",
  "NL-Lan",
  "NL-Lut",
  "NL-Mol",
  "CH-Oe2",
  "UK-ESa",
  "UK-Her",
  "DK-Sor",
  "FR-Fon",
  "FR-Hes",
  "DE-Hai",
  "UK-Ham",
  "UK-PL3",
  "AU-Tum",
  "AU-Wac",
  "FR-LBr",
  "DE-Tha",
  "DE-Wet",
  "IT-Lav",
  "NL-Loo",
  "UK-Gri",
  "IT-Ren",
  "DK-Lva",
  "FR-Lq1",
  "FR-Lq2",
  "DE-Gri",
  "DE-Meh",
  "HU-Bug",
  "HU-Mat",
  "IE-Dri",
  "IT-LMa",
  "IT-Mal",
  "IT-MBo",
  "NL-Ca1",
  "NL-Haa",
  "NL-Hor",
  "ES-VDA",
  "CH-Oe1",
  "UK-EBu",
  "BE-Jal",
  "BE-Vie",
  "CZ-wet",
  "SE-Faj",
  "UK-AMo",
  "US-Bo1",
  "US-Ne1",
  "US-Ne2",
  "US-Ne3",
  "US-Los",
  "US-Oho",
  "US-Bar",
  "US-Ha1",
  "US-UMB",
  "US-WCr",
  "US-Wi1",
  "US-Wi8",
  "RU-Fyo",
  "SE-Nor",
  "SE-Sk1",
  "SE-Sk2",
  "US-Ho1",
  "US-Ho2",
  "US-Wi0",
  "US-Wi2",
  "US-Wi4",
  "US-Wi5",
  "US-Wi9",
  "US-Bkg",
  "US-PFa",
  "US-Syv",
  "US-Wi6",
  "US-Wi7",
  "CA-Mer"
  )

tropical.all <- c()
tropical.all$name <- "Tropical_All"
tropical.all$stations <- c(
  "AU-How",
  "BR-Sa3",
  "AU-Fog",
  "ID-Pag"
  )

# Directory Paths:
lue.file.path <- "/home/user/Projects/gepisat/results/2002-06/attempt17/lue/"
#lue.file.path <- "/Users/twdavis/Projects/gepisat/results/2002-06/lue/"
#lue.file.path <- "C:\\Users\\Tyler\\Projects\\gepisat\\results\\2002-06\\lue\\"

# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### MAIN #####################################################################
# /////////////////////////////////////////////////////////////////////////////
# -------------------
# Process annual GPP:
# -------------------
all.stations <- process_gpp(lue.file.path, all.stations)
box_and_whisker(all.stations, write.out=0, lue.file.path)

# -------------------
# Process Boreal ENF:
# -------------------
boreal.enf <- process_gpp(lue.file.path, boreal.enf)
box_and_whisker(boreal.enf, write.out=0, lue.file.path)

# ----------------------
# Process Temperate DBF:
# ----------------------
temperate.dbf <- process_gpp(lue.file.path, temperate.dbf)
box_and_whisker(temperate.dbf, write.out=0, lue.file.path)

# ----------------------
# Process Temperate ENF:
# ----------------------
temperate.enf <- process_gpp(lue.file.path, temperate.enf)
box_and_whisker(temperate.enf, write.out=0, lue.file.path)

# ---------------------
# Process Temperate MF:
# ---------------------
temperate.mf <- process_gpp(lue.file.path, temperate.mf)
box_and_whisker(temperate.mf, write.out=0, lue.file.path)

# ---------------------
# Process Tropical EBF:
# ---------------------
tropical.ebf <- process_gpp(lue.file.path, tropical.ebf)
box_and_whisker(tropical.ebf, write.out=0, lue.file.path)

# ----------------------
# Process Temperate CRO:
# ----------------------
temperate.cro <- process_gpp(lue.file.path, temperate.cro)
box_and_whisker(temperate.cro, write.out=1, lue.file.path)

# ------------------------
# Process Subtropical CRO:
# ------------------------
subtropical.cro <- process_gpp(lue.file.path, subtropical.cro)
box_and_whisker(subtropical.cro, write.out=1, lue.file.path)

# ------------------------
# Process Subtropical DBF:
# ------------------------
subtropical.dbf <- process_gpp(lue.file.path, subtropical.dbf)
box_and_whisker(subtropical.dbf, write.out=1, lue.file.path)

# ------------------------
# Process Subtropical EBF:
# ------------------------
subtropical.ebf <- process_gpp(lue.file.path, subtropical.ebf)
box_and_whisker(subtropical.ebf, write.out=1, lue.file.path)

# ------------------------
# Process Subtropical ENF:
# ------------------------
subtropical.enf <- process_gpp(lue.file.path, subtropical.enf)
box_and_whisker(subtropical.enf, write.out=1, lue.file.path)

# -------------------
# Process Boreal GRA:
# -------------------
boreal.gra <- process_gpp(lue.file.path, boreal.gra)
box_and_whisker(boreal.gra, write.out=1, lue.file.path)

# -------------------
# Process Arid GRA:
# -------------------
arid.gra <- process_gpp(lue.file.path, arid.gra)
box_and_whisker(arid.gra, write.out=1, lue.file.path)

# ------------------------
# Process Subtropical GRA:
# ------------------------
subtropical.gra <- process_gpp(lue.file.path, subtropical.gra)
box_and_whisker(subtropical.gra, write.out=1, lue.file.path)

# ----------------------
# Process Temperate GRA:
# ----------------------
temperate.gra <- process_gpp(lue.file.path, temperate.gra)
box_and_whisker(temperate.gra, write.out=1, lue.file.path)

# -------------------
# Process Boreal OSH:
# -------------------
boreal.osh <- process_gpp(lue.file.path, boreal.osh)
box_and_whisker(boreal.osh, write.out=1, lue.file.path)

# ------------------------
# Process Climate regions:
# ------------------------
arid.all <- process_gpp(lue.file.path, arid.all)
boreal.all <- process_gpp(lue.file.path, boreal.all)
subtropical.all <- process_gpp(lue.file.path, subtropical.all)
temperate.all <- process_gpp(lue.file.path, temperate.all)
tropical.all <- process_gpp(lue.file.path, tropical.all)
box_and_whisker(tropical.all, write.out=0, lue.file.path)

# --------------------
# Save output to file:
# --------------------
results <- cbind(
  all.stations$mean, all.stations$sdev, all.stations$max, all.stations.min,
  arid.gra$mean, arid.gra$sdev, arid.gra$max, arid.gra$min,
  boreal.enf$mean, boreal.enf$sdev, boreal.enf$max, boreal.enf$min,
  boreal.gra$mean, boreal.gra$sdev, boreal.gra$max, boreal.gra$min,
  boreal.osh$mean, boreal.osh$sdev, boreal.osh$max, boreal.osh$min,
  subtropical.cro$mean, subtropical.cro$sdev, subtropical.cro$max, subtropical.cro$min,
  subtropical.dbf$mean, subtropical.dbf$sdev, subtropical.dbf$max, subtropical.dbf$min,
  subtropical.ebf$mean, subtropical.ebf$sdev, subtropical.ebf$max, subtropical.ebf$min,
  subtropical.enf$mean, subtropical.enf$sdev, subtropical.enf$max, subtropical.enf$min,
  subtropical.gra$mean, subtropical.gra$sdev, subtropical.gra$max, subtropical.gra$min,
  temperate.cro$mean, temperate.cro$sdev, temperate.cro$max, temperate.cro$min,
  temperate.dbf$mean, temperate.dbf$sdev, temperate.dbf$max, temperate.dbf$min,
  temperate.enf$mean, temperate.enf$sdev, temperate.enf$max, temperate.enf$min,
  temperate.gra$mean, temperate.gra$sdev, temperate.gra$max, temperate.gra$min,
  temperate.mf$mean, temperate.mf$sdev, temperate.mf$max, temperate.mf$min,
  tropical.ebf$mean, tropical.ebf$sdev, tropical.ebf$max, tropical.ebf$min
  )
col.heads <- c(
  "All_Stations_ave", "All_Stations_std", "All_Stations_max", "All_Stations_min",
  "Arid_GRA_ave", "Arid_GRA_std", "Arid_GRA_max", "Arid_GRA_min",
  "Boreal_ENF_ave", "Boreal_ENF_std", "Boreal_ENF_max", "Boreal_ENF_min",
  "Boreal_GRA_ave", "Boreal_GRA_std", "Boreal_GRA_max", "Boreal_GRA_min",
  "Boreal_OSH_ave", "Boreal_OSH_std", "Boreal_OSH_max", "Boreal_OSH_min",
  "Subtropical_CRO_ave", "Subtropical_CRO_std", "Subtropical_CRO_max", "Subtropical_CRO_min",
  "Subtropical_DBF_ave", "Subtropical_DBF_std", "Subtropical_DBF_max", "Subtropical_DBF_min",
  "Subtropical_EBF_ave", "Subtropical_EBF_std", "Subtropical_EBF_max", "Subtropical_EBF_min",
  "Subtropical_ENF_ave", "Subtropical_ENF_std", "Subtropical_ENF_max", "Subtropical_ENF_min",
  "Subtropical_GRA_ave", "Subtropical_GRA_std", "Subtropical_GRA_max", "Subtropical_GRA_min",
  "Temperate_CRO_ave", "Temperate_CRO_std", "Temperate_CRO_max", "Temperate_CRO_min",
  "Temperate_DBF_ave", "Temperate_DBF_std", "Temperate_DBF_max", "Temperate_DBF_min",
  "Temperate_ENF_ave", "Temperate_ENF_std", "Temperate_ENF_max", "Temperate_ENF_min",
  "Temperate_GRA_ave", "Temperate_GRA_std", "Temperate_GRA_max", "Temperate_GRA_min",
  "Temperate_MF_ave", "Temperate_MF_std", "Temperate_MF_max", "Temperate_MF_min",
  "Tropical_EBF_ave", "Tropical_EBF_std", "Tropical_EBF_max", "Tropical_EBF_min"
  )


results <- cbind(
  all.stations$mean, all.stations$sdev, all.stations$max, all.stations$min,
  arid.all$mean, arid.all$sdev, arid.all$max, arid.all$min,
  boreal.all$mean, boreal.all$sdev, boreal.all$max, boreal.all$min,
  subtropical.all$mean, subtropical.all$sdev, subtropical.all$max, subtropical.all$min,
  temperate.all$mean, temperate.all$sdev, temperate.all$max, temperate.all$min,
  tropical.all$mean, tropical.all$sdev, tropical.all$max, tropical.all$min
)
col.heads <- c(
  "All_Stations_ave", "All_Stations_std", "All_Stations_max", "All_Stations_min",
  "Arid_ave", "Arid_std", "Arid_max", "Arid_min",
  "Boreal_ave", "Boreal_std", "Boreal_max", "Boreal_min",
  "Subtropical_ave", "Subtropical_std", "Subtropical_max", "Subtropical_min",
  "Temperate_ave", "Temperate_std", "Temperate_max", "Temperate_min",
  "Tropical_ave", "Tropical_std", "Tropical_max", "Tropical_min"
  )


row.heads <- c(
  "Jan", "Feb", "Mar", "Apr", 
  "May", "Jun", "Jul", "Aug",
  "Sep", "Oct", "Nov", "Dec"
  )

r.df <- as.data.frame(results, row.names=row.heads)
names(r.df) <- col.heads
write.table(
  r.df, 
  file=paste(lue.file.path, "GPP_Timeseries.txt", sep=""),
  sep=",",
  row.names=TRUE,
  col.names=TRUE
  )
