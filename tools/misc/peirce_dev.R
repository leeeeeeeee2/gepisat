# ************************************************************************
# * Name: peirce_dev
# *
# * Input: - double, total number of observations (N)
# *        - double, number of outliers to be removed (n)
# *        - double, number of unknowns in model (m)
# *
# * Return: double, 
# *
# * Features: Calculates the squared-error tolerance for outlier 
# *           identification
# *
# * Ref: Peirce, B. (1852) Criterion for the rejection of doubtful
# *        observations, Astronomical Journal, vol. 2 (21), pp. 161-163.
# ************************************************************************
peirce_dev <- function(N, n, m){
  # Check number of observations:
  if (N > 1){
    # Calculate Q (Nth root of Gould's equation B):
    # Note: 1/N exponent is factored to each individual term to prevent
    # OverflowError with large N (e.g., >142)
    Q <- (n^(n/N)*(N - n)^((N - n)/N))/N
    #
    # Initialize R values (as floats):
    Rnew <- 1.0
    Rold <- 0.0  # <- Necessary to prompt while loop
    #
    while(abs(Rnew-Rold) > (N*2.0e-16)){
      # Calculate Lamda (1/(N-n)th root of Gould's equation A'):
      ldiv <- Rnew^n
      if (ldiv == 0){
        ldiv <- 1.0e-6
      }
      Lamda <- ((Q^N)/(ldiv))^(1.0/(N - n))
      #
      # Calculate x-squared (straight-forward Gould's equation C):
      x2 <- 1.0 + (N - m - n)/n * (1.0 - Lamda^2.0)
      #
      # If x2 goes negative, return as zero:
      if (x2 < 0){
        x2 <- 0
        Rold <- Rnew
      } else {
        #
        # Use x-squared to update R (Gould's equation D):
        # NOTE: error function (erfc) is replaced with pnorm:
        # source: 
        # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/Normal.html
        Rold <- Rnew
        Rnew <- exp((x2 - 1)/2.0)*(2*pnorm(
          sqrt(x2)/sqrt(2)*sqrt(2), 
          lower=FALSE))
      }
    }
  } else {
    x2 <- 0
  }
  x2
}