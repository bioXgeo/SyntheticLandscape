#' Estimate the areal autocorrelation function.
#'
#' Calculates the areal autocorrelation function (AACF) as the
#' inverse of the Fourier power spectrum. \code{aacf(x)} returns
#' the AACF in both matrix and image format.
#'
#' @param x An n x n raster object.
#' @return A list containing matrix and image representations
#'   of the AACF. Both matrix and image values are normalized
#'   so that the maximum is equal to 1.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(normforest)
#'
#' # plot resulting aacf image
#' plot(aacf_list[[2]])
aacf <- function(x) {

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  # convert raster values to matrix
  zmat <- matrix(getValues(x), ncol = M, nrow = N)

  # create windows to prevent leakage
  wc <- hanning(M)
  wr <- hanning(N)

  # create matrix of weights
  w <- meshgrid(wr, wc)
  w <- w[[1]] * w[[2]]

  # apply window weights
  zmatw <- zmat * w

  # perform Fourier transform
  ft <- fft(zmatw)

  # calculate power spectrum
  ps <- (abs(ft) ^ 2) / (M * N)

  # autocorrelation function
  af <- Re(fft(ps, inverse = TRUE) / (M * N))
  af_shift <- fftshift(af) # this should be symmetric!

  # normalize to max 1
  af_norm <- af_shift / max(as.numeric(af_shift), na.rm = TRUE)

  # set values of new raster
  af_img <- setValues(x, af_norm)

  return (list(af_norm, af_img))
}


#' Calculate correlation length.
#'
#' Calculates the smallest and largest distances to specified autocorrelation
#' values (e.g., 0.2) of the areal autocorrelation function (AACF). All 180
#' degrees from the origin of the AACF image are considered for the calculation.
#'
#' @param aacf An n x n raster of the areal autocorrelation function,
#'   produced by the function \code{aacf(x)}.
#' @param threshold A numeric vector containing values between 0 and 1. Indicates
#'   the autocorrelation values to which the rates of decline are measured.
#' @param plot Logical. Defaults to \code{FALSE}. If \code{TRUE}, the AACF and
#'   lines showing the considered directions of autocorrelation from the origin
#'   will be plotted.
#' @return A list containing the minimum and maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation values < 1.
#'   Distances are in the units of the x, y coordinates of the raster image. If more
#'   than one threshold value is specified, the order of this list will be
#'   \code{[minval(t1), maxval(t1), minval(t2), maxval(t2)]}.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(normforest)
#'
#' # estimate the fastest/slowest declines to 0.20 and 0.37 (1/e) autocorrelation
#' sclvals <- scl(aacf_list[[2]])
#'
#' # calculate Scl20, the minimum distance to an autocorrelation value of 0.2 in the AACF
#' Scl20 <- sclvals[[1]]
scl <- function(x, threshold = c(0.20, 1 / exp(1)), plot = FALSE) {

  # take amplitude image, cut in half (y direction)
  half_dist <- (ymax(x) - ymin(x)) / 2
  ymin <- ymax(x) - half_dist
  x <- crop(x, c(xmin(x), xmax(x), ymin, ymax(x)))

  # get origin of image (actually bottom center)
  origin <- c(mean(coordinates(x)[, 1]), ymin(x))

  ### line calculations are taken from the plotrix function draw.radial.line
  # calculate rays extending from origin
  M <- 180
  j <- seq(0, (M - 1))
  alpha <- (pi * j) / M # angles
  px <- c(0, half_dist) # line length
  linex <- unlist(lapply(seq(1, length(alpha)), function(x) origin[1] + px * cos(alpha[x])))
  liney <- unlist(lapply(seq(1, length(alpha)), function(x) origin[2] + px * sin(alpha[x])))
  linelist <- lapply(seq(1, length(linex), 2),
                     FUN = function(i) Lines(Line(cbind(linex[i:(i + 1)], liney[i:(i + 1)])),
                                             ID = paste('l', i, sep = '')))
  lines <- SpatialLines(linelist, proj4string = CRS(proj4string(x)))

  # plot and calculate amplitude sums along rays
  if(plot == TRUE) {
    plot(x)
    lines(lines)
  }

  # get values for all points along line
  Aalpha <- list()
  for (i in 1:length(lines)) {
    Aalpha[[i]] <- extract(x, lines[i], along = TRUE, cellnumbers = TRUE)
  }

  # each line has length = half_dist, with each point approx. 1 pixel apart
  fast_dists <- list()
  for (i in 1:length(threshold)) {
    fast_dists[[i]] <- suppressWarnings(mindist(threshold[i], Aalpha, x))
  }

  slow_dists <- list()
  for (i in 1:length(threshold)) {
    slow_dists[[i]] <- suppressWarnings(maxdist(threshold[i], Aalpha, x))
  }

  return(c(fast_dists, slow_dists))
}

#' Estimate minimum correlation length.
#'
#' Internal function to calculates the minimum distances to specified
#' autocorrelation values (e.g., 0.2) of the areal autocorrelation
#' function (AACF). All 180 degrees from the origin of the AACF image
#' are considered for the calculation.
#'
#' @param threshold A number with a value between 0 and 1. Indicates
#'   the autocorrelation value to which the rates of decline are measured.
#' @param Aalpha An list of dataframes produced by \code{scl()} that contain
#'   the AACF values along lines extending in multiple directions from the
#'   AACF origin (autocorrelation = 1).
#' @param aacf An 0.5n x n raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @return A list containing the minimum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are in the units of the x, y coordinates of the raster image.
mindist <- function(threshold, Aalpha, aacfimg) {
  # get index of minimum <= threshold value
  decay_ind <- list()
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= threshold), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- coordinates(aacfimg)[, 1]
  y <- coordinates(aacfimg)[, 2]
  origin <- c(mean(coordinates(aacfimg)[, 1]), ymin(aacfimg))
  decay_celln <- list()
  for (j in 1:length(Aalpha)) {
    decay_celln[[j]] <- unlist(lapply(Aalpha[[j]], FUN = function(x) x[decay_ind[j], 1]))
  }
  decay_celln <- unlist(decay_celln)
  decay_coords <- data.frame(x = x[decay_celln], y = y[decay_celln])
  decay_coords <- decay_coords[decay_ind != Inf,]
  decay_dist <- min(crossdist.default(X = decay_coords$x[!is.na(decay_coords$x)], Y = decay_coords$y[!is.na(decay_coords$x)],
                                      x2 = origin[1], y2 = origin[2]))

  return(decay_dist)
}

#' Estimate maximum correlation length.
#'
#' Internal function to calculates the maximum distances to specified
#' autocorrelation values (e.g., 0.2) of the areal autocorrelation
#' function (AACF). All 180 degrees from the origin of the AACF image
#' are considered for the calculation.
#'
#' @param threshold A number with a value between 0 and 1. Indicates
#'   the autocorrelation value to which the rates of decline are measured.
#' @param Aalpha An list of dataframes produced by \code{scl()} that contain
#'   the AACF values along lines extending in multiple directions from the
#'   AACF origin (autocorrelation = 1).
#' @param aacf An 0.5n x n raster of the areal autocorrelation function. This
#'   is the AACF raster split in two in terms of height.
#' @return A list containing the maximum distances from an
#'   autocorrelation value of 1 to the specified autocorrelation value < 1.
#'   Distances are in the units of the x, y coordinates of the raster image.
maxdist <- function(threshold, Aalpha, aacfimg) {
  # get index of minimum <= threshold value
  decay_ind <- list()
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= threshold), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- coordinates(aacfimg)[, 1]
  y <- coordinates(aacfimg)[, 2]
  origin <- c(mean(coordinates(aacfimg)[, 1]), ymin(aacfimg))
  decay_celln <- list()
  for (j in 1:length(Aalpha)) {
    decay_celln[[j]] <- unlist(lapply(Aalpha[[j]], FUN = function(x) x[decay_ind[j], 1]))
  }
  decay_celln <- unlist(decay_celln)
  decay_coords <- data.frame(x = x[decay_celln], y = y[decay_celln])
  decay_coords <- decay_coords[decay_ind != Inf,]
  decay_dist <- max(crossdist.default(X = decay_coords$x[!is.na(decay_coords$x)], Y = decay_coords$y[!is.na(decay_coords$x)],
                                      x2 = origin[1], y2 = origin[2]))

  return(decay_dist)
}
