etsrad <- function(user.lon, user.lat, user.day){
  # This script calculates the half-hourly extraterrestrial solar radiation
  # based on the user's coordinates and the time of year
  #
  # User-provided variable definitions:
  #   user.lon :: longitude [decimal degrees]; range: -180 to 180
  #   user.lat :: latitude [decimal degrees]; range: -90 to 90
  #   user.day :: Julian day of the year [day]; range: 1 to 366
  if (user.lon > 180 || user.lon < -180){
    stop("Warning: Longitude outside range of validity (-180 to 180)!")
  }
  if (user.lat > 90 || user.lat < -90){
    stop("Warning: Latitude outside range of validity (-90 to 90)!")
  }
  if (user.day < 1 || user.day > 366){
    stop("Warning: Day outside range of validity (1 to 366)!")
  }
  #
  ########################### FUNCTION VARIABLE ###############################
  et.srad <- list()
  #
  ########################### FUNCTION CONSTANTS ##############################
  # Half-hourly time series:
  local.hh <- (seq(48) - 1.0)*0.5
  #
  # One-second time series (for calculating total daylight hours):
  local.sec <- (seq(86400) - 1.0)/3600.0
  #
  # Solar constant [W m-2]:
  # ref: Duffie & Beckman (2013); Wehrli (1985)
  kGsc = 1367.0
  #
  # Conversion factor from solar radiation to PPFD [umol J-1]:
  # ref: Chen et al. (1993); Meek et al. (1984)
  kSRadToPPFD = 2.04
  #
  # Total number of days per year [day]:
  kDaysPerYear = 365.242
  #
  ######################## STEP 1. TIME ZONE HOUR #############################
  # Correction for local time zone; based on one hour per 15 degrees longitude:
  # * positive or negative deviation from UTC/GMT
  if (user.lon < 0){
    temp.lon <- -1.0 * user.lon
    temp.tzh <- floor(temp.lon / 15)
    tz.hour <- -1.0 * temp.tzh
  } else {
    tz.hour <- floor(user.lon / 15)
  }
  #
  ######################## STEP 2. ECCENTRICITY ###############################
  # Correction for the eccentricity of Earth's orbit around the sun
  # ref: Spencer (1971)
  dr <- 1.0 + 0.033*cos(2.0 * pi * user.day / kDaysPerYear)
  #
  ######################## STEP 3. DECLINATION ANGLE ##########################
  # Declination angle of the Earth [degrees]
  # * ranges from -23.45 in the winter to 23.45 in the summer
  # ref: Cooper (1969)
  delta <- 23.45 * sin(2.0 * pi / kDaysPerYear * (284 + user.day))
  #
  ######################## STEP 4. EQUATION OF TIME ###########################
  # The equation of time calculates the difference between apparent solar time 
  # and mean solar time [min]
  # ref: Woolf (1968); Stine and Geyer (2001)
  B <- (user.day - 1) * 2 * pi / kDaysPerYear
  EOT <- 0.258*cos(B) - 7.416*sin(B) - 3.648*cos(2*B) - 9.228*sin(2*B)
  #
  ######################## STEP 5. LONGITUDE ##################################
  # Correction factor for longitude [hr]
  # ref: Stine and Geyer (2001)
  LC <- (15 * tz.hour - user.lon)/15
  #
  ######################## STEP 6. SOLAR TIME #################################
  # Calculate the solar time based on the local time [hr]
  # ref: Stine and Geyer (2001)
  ts.hh <- local.hh + EOT/60 - LC
  ts.sec <- local.sec + EOT/60 - LC
  #
  ######################## STEP 7. HOUR ANGLE #################################
  # Solar zenith angle [degrees]
  # * solar noon is at 0 degrees (i.e., vertical line)
  # * based once again on 15 deg movement per hour
  # ref: Stine and Geyer (2001)
  w.hh <- 15*(ts.hh - 12)
  w.sec <- 15*(ts.sec - 12)
  #
  ######################## STEP 8. ET SRAD ####################################
  # Extraterrestrial solar radiation w.r.t. a flat horizon
  # ref: Eq. 1.10.2 in Duffie and Beckman (2013)
  Ra.hh <- kGsc * dr * (
    cos(user.lat*pi/180) * cos(delta*pi/180) * cos(w.hh*pi/180)
    + sin(user.lat*pi/180) * sin(delta*pi/180)
    )
  Ra.sec <- kGsc * dr * (
    cos(user.lat*pi/180) * cos(delta*pi/180) * cos(w.sec*pi/180)
    + sin(user.lat*pi/180) * sin(delta*pi/180)
  )
  Ra.hh[Ra.hh < 0 ] <- 0
  et.srad$srad <- Ra.hh
  et.srad$ppfd <- Ra.hh * kSRadToPPFD
  #
  ######################## STEP 9. DAYLIGHT HOURS #############################
  # Convert one-second solar radiation to PPFD (zero negative values):
  ppfd.sec <- Ra.sec * kSRadToPPFD
  ppfd.sec[ppfd.sec < 0] <- 0
  #
  # Hourly time factor (i.e., dt) [hr]
  dt <- 24.0/86400.0
  #
  # Count the number of positive radiation values:
  num.pos <- length(ppfd.sec[ppfd.sec > 0])
  #
  # Multiply number of positive occurrences by time factor [hr]:
  daylight.hours <- num.pos * dt
  et.srad$dayhours <- daylight.hours
  #
  ######################## STEP 10. PPFD INTEGRATION ##########################
  # Integrate PPFD over the entire day [mol m-2]
  # * Simpson's rule for integrating one-second PPFD over the day:
  #   int(fx) <- (h/3) * [f(x1) + 4f(x2) + 2f(x3) ... f(xn)]
  h <- 1.0
  s <- ppfd.sec[1] + ppfd.sec[length(ppfd.sec)]
  for (i in seq(from=2, to=length(ppfd.sec)+1, by=2)){
    s <- s + 4 * ppfd.sec[i]
  }
  for (j in seq(from=3, to=length(ppfd.sec), by=2)){
    s <- s + 2 * ppfd.sec[j]
  }
  s <- s * h / 3
  ppfd.integral <- s / 1000000.0
  et.srad$dayppfd <- ppfd.integral
  #
  ############################### REFERENCES ##################################
  # Chen, J. M., T. A. Black, D. T. Price, and R. E. Carter (1993). “Model for
  #   calculating photosynthetic photon flux densities in forest openings on
  #   slopes”. In: Journal of Applied Meteorology 32.10, pp. 1656–1665
  # Cooper, P. I. (1969). “The absorption of radiation in solar stills”. In: 
  #   Solar Energy 12.3, pp. 333–346 
  # Duffie, J. A. and W. A. Beckman (2013). Solar engineering of thermal pro-
  #   cesses. 4th ed. New Jersey: John Wiley and Sons
  # Meek, D. W., J. L. Hatfield, T. A. Howell, S. B. Idso, and R. J. Reginato
  #   (1984). “A generalized relationship between photosynthetically active
  #   radiation and solar radiation”. In: Agronomy Journal 76, pp. 939–945
  # Spencer, J. W. (1971). “Fourier series representation of the position of
  #   the sun”. In: Search 2, p. 172.
  # Stine, W. B. and M. Geyer (2001). “Power from the Sun”. In: Available on-
  #   line: http://www.powerfromthesun.net/Book/chapter03/chapter03.
  # Wehrli, C. (1985). Extraterrestrial solar spectrum. Tech. rep. 615. Davos
  #   Dorf, Switzerland: Physikalisch-Meteorologisches Observatorium/World
  #   Radiation Center (PMO/WRC)
  # Woolf, H. M. (1968). On the computation of solar evaluation angles and the
  #   determination of sunrise and sunset times. Tech. rep. NASA-TM-X-164.
  #   National Aeronautics and Space Administration (NASA).
  #
  ############################## RETURN VALUES ################################
  et.srad
}