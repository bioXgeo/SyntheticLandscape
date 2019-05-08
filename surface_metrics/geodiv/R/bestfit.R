# functions to find the best fit line for forty percent of slope with lowest slope

#' @export
#' Determines the slopes along a the bearing area curve.
#'
#' Calculates the slopes along the bearing area curve
#' of a raster. Slopes are determined at points x,
#' from point x - h to x + h.
#'
#' @param x A raster.
#' @param h Spacing before and after each point.
#' 2h is the distance over which slopes are calculated.
#' @param f Bearing area function as calculated with
#' bearing_area.
#' @return A dataframe with the slope for each segment
#' with centerpoint x.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the slopes along the bearing area curve
#' ba <- bearing_area(normforest)
#' x <- seq(0, 1, length.out = 100000)
#' slopes <- slopecalc(x = x, h = 0.01, f = ba)
slopecalc <- function(x, h, f) {
  xplus <- x + h
  xminus <- x - h
  x <- x

  # figure out ends (need distance on both sides of x)
  space <- 1 / length(x)
  end_length <- h / space
  if (end_length < 1) { end_length = 1}
  begin <- ceiling(1 + end_length)
  end <- floor(length(x) - end_length)
  begin1 <- begin - 1
  end1 <- end + 1

  fxh_pos <- seq(1, length(x))
  fxh_neg <- seq(1, length(x))
  fxh_pos[1:begin1] <- (1 - stats::quantile(f, probs = xplus[1:begin1]))
  fxh_neg[1:begin1] <- (1 - stats::quantile(f, probs = x[1:begin1]))
  # variation on newton's difference quotient at far end
  fxh_pos[end1:length(x)] <- (1 - stats::quantile(f, probs = x[end1:length(x)]))
  fxh_neg[end1:length(x)] <- (1 - stats::quantile(f, probs = xminus[end1:length(x)]))
  # symmetric difference quotient everywhere else (99800 points)
  fxh_pos[begin:end] <- (1 - stats::quantile(f, probs = xplus[begin:end]))
  fxh_neg[begin:end] <- (1 - stats::quantile(f, probs = xminus[begin:end]))

  diff_quo <- (fxh_pos - fxh_neg) / (2 * h)
  slopes <- data.frame(slope = diff_quo, x = x)

  return(slopes)
}

#' @export
#' Determines the average slope along larger segments of
#' the bearing area curve of a raster.
#'
#' Calculates the average slope over every segment
#' of a specified percentage length of the total bearing
#' area curve.
#'
#' @param slopes A dataframe containing all slopes along
#' the bearing area curve, calculated using the slopecalc
#' function.
#' @param l Percentage of the curve over which to calculate
#' mean slope.
#' @return A dataframe with the average slope over segments
#' beginning at specified x locations along the bearing area
#' curve. 'slope' represents the mean slope over the segment,
#' 'xstart' is the beginning x location of the segment, and
#' 'xend' is the concluding x location of the segment.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # find the average slope of segments of the bearing area
#' # curve.
#' ba <- bearing_area(normforest)
#' x <- seq(0, 1, length.out = 100000)
#' slopes <- slopecalc(x = x, h = 0.01, f = ba)
#' slopes_forty <- slopemeans(slopes = slopes, l = 0.4)
slopemeans <- function(slopes, l = 0.4) {
  x <- slopes$x
  slope <- slopes$slope
  length <- length(x) * l
  end_ind <- length(x) - length
  begin_ind <- 1 + length
  xstart <- x[1:end_ind]
  istart <- seq(1, end_ind)
  xend <- x[begin_ind:length(x)]
  iend <- seq(begin_ind, length(x))

  # create matrix of slopes
  xmat <- matrix(c(xstart, xend), ncol = 2)
  mean <- sapply(seq(1, end_ind), function(i) {mean(abs(slope[istart[i]:iend[i]]))})

  data <- data.frame(slope = mean, xstart = xstart, xend = xend)
  return(data)
}
