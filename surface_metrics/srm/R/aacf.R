#' Areal autocorrelation function.
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
#' data(orforest)
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(orforest)
#'
#' # plot resulting aacf image
#' plot(aacf_list[[2]])
aacf <- function(x) {

  # get raster dimensions
  M <- ncol(x)
  N <- nrow(x)

  zmat <- matrix(getValues(x), ncol = M, nrow = N)

  # create windows to prevent leakage
  wc <- hanning(M)
  wr <- hanning(N)

  w <- meshgrid(wr, wc)
  w <- w[[1]] * w[[2]]

  # apply window weights
  zmatw <- zmat * w

  # Fourier transform
  ft <- fft(zmatw)

  # power spectrum
  ps <- (abs(ft) ^ 2) / (M * N)

  # autocorrelation function
  af <- Re(fft(ps, inverse = TRUE) / (M * N))
  af_shift <- fftshift(af) # this should be symmetric!

  # normalize to max 1
  af_norm <- af_shift / max(as.numeric(af_shift), na.rm = TRUE)

  af_img <- setValues(rast, af_norm)

  if (plot == TRUE) {
    plot(af_img)
  }

  return (list(af_norm, af_img))
}


#' Correlation length.
#'
#' Calculates the distances to specified autocorrelation
#' values (e.g., 20%) of the areal autocorrelation function (AACF).
#'
#' @param aacf An n x n raster of the areal autocorrelation function,
#'   produced by the function \code{aacf(x)}.
#' @return A list containing matrix and image representations
#'   of the AACF. Both matrix and image values are normalized
#'   so that the maximum is equal to 1.
#' @examples
#' # import raster image
#' data(orforest)
#'
#' # calculate aacf img and matrix
#' aacf_list <- aacf(orforest)
#'
#' # plot resulting aacf image
#' plot(aacf_list[[2]])
scl <- function(aacf, threshold = c(0.20, 1 / exp(1))) {
  # borrow most from fourier functions

  # take amplitude image, cut in half (y direction)
  half_dist <- (ymax(aacf_img) - ymin(aacf_img)) / 2
  ymin <- ymax(aacf_img) - half_dist
  aacf_img <- crop(aacf_img, c(xmin(aacf_img), xmax(aacf_img), ymin, ymax(aacf_img)))

  # get origin of image (actually bottom center)
  origin <- c(mean(coordinates(aacf_img)[,1]), ymin(aacf_img))

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
  lines <- SpatialLines(linelist, proj4string = CRS(proj4string(aacf_img)))

  # plot and calculate amplitude sums along rays
  if(plot == TRUE) {
    plot(aacf_img)
    lines(lines)
  }

  # get values for all points along line
  Aalpha <- list()
  for (i in 1:length(lines)) {
    Aalpha[[i]] <- extract(aacf_img, lines[i], along = TRUE, cellnumbers = TRUE)
  }

  # each line has length = half_dist, with each point approx. 1 pixel apart
  fast_dists <- list()
  for (i in 1:length(threshold)) {
    fast_dists[[i]] <- min_dist(threshold[i], Aalpha)
  }

  slow_dists <- list()
  for (i in 1:length(threshold)) {
    slow_dists[[i]] <- max_dist(threshold[i], Aalpha)
  }

  return(c(fast_dists, slow_dists))
}

# minimum  distance to threshold values
min_dist <- function(t, Aalpha) {
  # get index of minimum <= threshold value
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= t), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- coordinates(aacf_img)[, 1]
  y <- coordinates(aacf_img)[, 2]
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

max_dist <- function(t, Aalpha) {
  # get index of minimum <= threshold value
  for (j in 1:length(Aalpha)) {
    decay_ind[[j]] <- lapply(Aalpha[[j]], FUN = function(x)
      min(which(x[, 2] <= t), na.rm = TRUE))
  }
  decay_ind <- unlist(decay_ind)

  # get distance to minimum index below threshold value
  x <- coordinates(aacf_img)[, 1]
  y <- coordinates(aacf_img)[, 2]
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
