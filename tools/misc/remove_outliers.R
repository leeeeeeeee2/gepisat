# ************************************************************************
# * Name: remove_outliers
# *
# * Input: - numeric list (fc.data)
# *        - numeric list (ppfd.data)
# *        - numeric list (fitted.data)
# *
# * Return: numeric list
# *
# * Features: Returns outlier free data pairs
# *
# *           Depends on:
# *           - peirce_dev()
# *
# * Ref: Peirce, B. (1852) Criterion for the rejection of doubtful
# *        observations, Astronomical Journal, vol. 2 (21), pp. 161-163.
# ************************************************************************
remove_outliers <- function(fc.data, ppfd.data, fitted.data){
  # Identify and remove ouliers from fc:ppfd data pairs
  # Define Peirce's parameters:
  peirce.cap.n <- length(fc.data)
  peirce.lc.n <- 1.0
  peirce.m <- 3.0
  #
  # Calculate the mean squared error of fit:
  se <- (fc.data - fitted.data)^2
  mse <- sum(se)/(length(fc.data) - peirce.m)
  #
  # Calculate Peirce's tolerance:
  peirce.x.sqr <- peirce_dev(peirce.cap.n, peirce.lc.n, peirce.m)
  peirce.delta.sqr <- mse * peirce.x.sqr
  #
  # Count number of ouliers identified:
  outlier.index <- se > peirce.delta.sqr
  outliers.found <- length(which(outlier.index == "TRUE"))
  #
  # Run a second time in-case no outliers were found:
  if (outliers.found == 0){
    peirce.lc.n <- 2.0
    peirce.x.sqr <- peirce_dev(peirce.cap.n, peirce.lc.n, peirce.m)
    peirce.delta.sqr <- mse * peirce.x.sqr
    outlier.index <- se > peirce.delta.sqr
    outliers.found <- length(which(outlier.index == "TRUE"))
    # Reset Peirce's n:
    peirce.lc.n <- 1.0
  }
  # Loop Peirce's "n" until all outliers are accounted for:
  while (peirce.lc.n <= outliers.found){
    peirce.lc.n <- peirce.lc.n + 1
    peirce.x.sqr <- peirce_dev(peirce.cap.n, peirce.lc.n, peirce.m)
    peirce.delta.sqr <- mse * peirce.x.sqr
    outlier.index <- se > peirce.delta.sqr
    outliers.found <- length(which(outlier.index == "TRUE"))
  }
  #
  # Find where exceedance occurs and remove outliers:
  non.outlier.index <- se < peirce.delta.sqr
  fc.data.no <- fc.data[non.outlier.index]
  ppfd.data.no <- ppfd.data[non.outlier.index]
  outlier.free.data <- cbind(fc.data.no, ppfd.data.no)
}