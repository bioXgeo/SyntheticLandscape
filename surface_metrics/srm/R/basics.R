#' Calculates the average roughness of a surface.
#'
#' Finds the average roughness of a surface (Sa) as the absolute
#' deviation of surface heights from the mean surface height.
#' Height is measured as the value of a raster and may not
#' necessarily represent actual height.
#'
#' @param x A raster.
#' @return A value of average roughness in the units of the
#'   original raster.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the surface roughness
#' roughness <- sa(normforest)
sa <- function(x) {
  z <- getValues(x)
  zbar <- mean(z, na.rm = TRUE)
  N <- length(z)

  val <- sum(abs(z - zbar)) / N

  return(val)
}

# root mean square roughness = standard deviation of surface heights
sq <- function(x) {
  z <- getValues(x)

  val <- sd(z, na.rm = TRUE)

  return(val)
}

# surface skewness = asymmetry of surface height distribution
ssk <- function(x, adj = TRUE) {
  z <- getValues(x)
  zbar <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  N <- length(z)

  val_unadj <- (sum((z - zbar) ^ 3) / N) / (s ^ 3)

  if (adj == TRUE) {
    val <- (sqrt((N * (N - 1))) / (N - 2)) * val_unadj # adjusted Fisher-Pearson coefficient of skewness
  } else {
    val <- val_unadj # Fisher-Pearson coefficient of skewness
  }

  return(val)
}

# surface kurtosis = peaked-ness of surface distribution
sku <- function(x, excess = TRUE) {
  z <- getValues(x)
  zbar <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  N <- length(z)

  val_unadj <- (sum((z - zbar) ^ 4) / N) / (s ^ 4)

  if (excess == TRUE) {
    val <- val_unadj - 3 # excess kurtosis (i.e., diff from normal distribution kurtosis)
  } else {
    val <- val_unadj
  }

  return(val)
}

# maximum valley depth = lowest value in the landscape
sv <- function(x) {
  z <- getValues(x)

  val <- abs(min(z, na.rm = TRUE))

  return(val)
}

# maximum peak height = highest value in the landscape
sp <- function(x) {
  z <- getValues(x)

  val <- abs(max(z, na.rm = TRUE))

  return(val)
}

# mean peak height = average peak height
smean <- function(x) {
  z <- getValues(x)

  val <- mean(z[z > 0], na.rm = TRUE)

  return(val)
}
