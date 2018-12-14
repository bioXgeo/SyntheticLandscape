#' Calculates the rotated Bearing Area curve.
#'
#' Finds a rotated version of the Bearing Area (Abbott-Firestone)
#' curve from a raster. The resulting function should be
#' rotated 90 degrees clockwise to get the actual Bearing
#' Area curve.
#'
#' @param x A raster.
#' @return A function describing the rotated Bearing Area curve.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the rotated Bearing Area curve.
#' ba_func <- bearing_area(normforest)
#'
#' # rotate the values and re-plot
#' xval <- environment(ba_func)$y
#' yval <- (1 - environment(ba_func)$x)
#' plot(yval ~ xval)
bearing_area <- function(x) {
  z <- getValues(x)

  # basic values
  N <- length(z)
  s <- sd(z)
  zbar <- mean(z, na.rm = TRUE)

  f <- ecdf(1 - z)

  return(f)
}

#' Plots the Bearing Area curve.
#'
#' Calculates and plots the Bearing Area curve for a raster
#' using the \code{bearing_area()} function (with correctly
#' rotated results).
#'
#' If \code{divisions = TRUE}, the lines representing the
#' best fit line to the flattest 40 percent of the curve will be
#' shown, as well as both the x and y interception points
#' of that line.
#'
#' @param x A raster.
#' @param divisions Logical, defaults to \code{FALSE}. If
#'   \code{TRUE}, divisions of the curve will be plotted.
#'   See details section for more information.
#' @return Plots the Bearing Area curve for a raster.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # plot the bearing area curve
#' plot_ba_curve(normforest, divisions = TRUE)
plot_ba_curve <- function(x, divisions = FALSE) {
  f <- bearing_area(x)

  xval <- environment(f)$y
  yval <- (1 - environment(f)$x)

  plot(yval ~ xval)

  if (divisions == TRUE) {
    line_fit <- find_flat(x, perc = 0.4)

    abline(line_fit[[1]]$coefficients[[1]], line_fit[[1]]$coefficients[[2]], col = 'blue')
    abline(h = line_fit[[3]], col = 'red')
    abline(h = line_fit[[4]], col = 'red')
    abline(v = line_fit[[5]], col = 'green')
    abline(v = line_fit[[6]], col = 'green')
  }
}

#' Finds the flattest part of the Bearing Area curve.
#'
#' Locates the flattest x percentage of the Bearing Area
#' curve. Meant to locate the flattest 40 percent of the
#' Bearing Area curve as used in several roughness parameter
#' calculations.
#'
#' @param x A raster.
#' @param perc Numeric between 0 and 1. The percentage of
#'   the curve over which to fit the line.
#' @return A list containing the equation for the best fit
#'   line, the predicted values from that line, the high
#'   and low y-intercept values for the intersection points
#'   of the line with the Bearing Area curve, and the high
#'   and low x-intercept values for the intersection points
#'   of the line with the Bearing Area curve.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate the flattest 40% of the bearing area curve
#' line_data <- find_flat(normforest, perc = 0.4)
#'
#' # extract the equation of the line
#' bf_line <- line_data[[1]]
find_flat <- function(x, perc = 0.4) {
  f <- bearing_area(x)

  xval <- environment(f)$y
  yval <- (1 - environment(f)$x)

  # find 40% of curve with least decline
  # use symmetric difference quotient to estimate the derivative at evenly spaced points
  # then find 40% consecutive section with lowest mean slope
  even_x <- seq(0, 1, length.out = 100000)
  even_y <- (1 - quantile(f, probs = even_x))
  forty_length <- perc * length(even_x)
  h <- 0.001
  slopes <- slopecalc(even_x, h, f = f) # calculate slope at every point
  means <- slopemeans(slopes) # calculate averages for each 40% segment

  # x value of start of 40% section with smallest decline
  slope_min <- means[means$slope == min(means$slope),]

  # calculate least-squares line for 40% of curve with smallest decline (lowest slope)
  lm_data <- data.frame(x = xval[xval >= slope_min$xstart & xval <= slope_min$xend],
                        y = yval[xval >= slope_min$xstart & xval <= slope_min$xend])
  ls_line <- lm(y ~ x, data = lm_data)

  # get value of ls line between 0 and 1
  pred_data <- remove_rownames(data.frame(x = even_x, y = even_y))
  pred_data$y <- predict(ls_line, newdata = pred_data)

  # what is the ls line y-value at x = 0, x = 1?
  ls_int_high <- pred_data$y[pred_data$x == 0]
  ls_int_low <- pred_data$y[pred_data$x == 1]

  # Smr1/Smr2 = x values that correspond to cdf y values at ls_int_high/low
  Smr1 <- f(1 - ls_int_high)
  Smr2 <- f(1 - ls_int_low)

  return(list(ls_line, pred_data, ls_int_high, ls_int_low, Smr1, Smr2))
}

height_ba <- function(x, xval) {
  f <- bearing_area(x)

  val <- (1 - quantile(f, probs = c(xval))[[1]])

  return(val)
}

sdc <- function(x, low, high) {
  val_low <- height_ba(x, low)
  val_high <- height_ba(x, high)

  val <- val_low - val_high

  return(val)
}

# surface bearing index = ratio of Sq to height from top of surface to height at
# 5% of bearing area (z05)
sbi <- function(x) {
  Sq <- sq(x)
  z05 <- height_ba(x, 0.05)

  val <- Sq / z05

  return(val)
}

# valley fluid retention index = void volume (area under Abbott curve) in 'valley' zone
# see fig.2a from Kedron et al. (2018)
svi <- function(x) {
  f <- bearing_area(x)

  val <- simpsons_above(model = f, b = 1, a = 0.8, n = 500)

  return(val)
}

# core fluid retention index = void volume (area above Abbott curve) in the core zone
# see fig.2a from Kedron et al. (2018)
sci <- function(x) {
  f <- bearing_area(x)

  core_above <- simpsons_above(model = f, b = 1, a = 0.05, n = 1000)
  Svi <- svi(x)
  val <- core_above - Svi

  return(val)
}

# core roughness depth = height difference between the intersection points of the found least
# mean square line in the Abbott curve (see fig. 2b from Kedron et al. 2018)
sk <- function(x) {
  line_info <- find_flat(x, perc = 0.4)

  ls_int_high <- line_info[[3]]
  ls_int_low <- line_info[[4]]

  val <- abs(ls_int_high - ls_int_low)

  return(val)
}

# reduced valley depth = height of triangle drawn at 100% on Abbott curve
svk <- function(x) {
  f <- bearing_area(x)

  line_info <- find_flat(x, perc = 0.4)

  smr2 <- line_info[[5]]

  val <- abs((1 - quantile(f, probs = 1)) - (1 - quantile(f, probs = smr2)))[[1]]

  return(val)
}

# reduced peak height = height of upper left triangle in abbott curve
spk <- function(x) {
  f <- bearing_area(x)

  line_info <- find_flat(x, perc = 0.4)

  smr1 <- line_info[[4]]

  val <- abs((1 - quantile(f, probs = 0)) - (1 - quantile(f, probs = smr1)))[[1]]

  return(val)
}
