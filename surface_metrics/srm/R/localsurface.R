# figure out peaks and valleys of surface

# from Image Metrology:
# local minimum = points where all 8 surrounding points are higher & below zero
# local maximum = points where all 8 surrounding points are lower & above zero


#' Find local peaks.
#'
#' Locates local peaks on a raster. A peak is defined as any pixel where
#' all 8 surrounding pixels have lower values, and the center pixel
#' has a positive value.
#'
#' @param x An n x n raster object.
#' @return A dataframe of local peak locations (\code{x, y}) and
#'   values (\code{val}). The raster location index (\code{ind}),
#'   row (\code{row}), and column (\code{col}) are also listed.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate peaks
#' peaks <- findpeaks(normforest)
#'
#' # calculate the summit density (# peaks/area)
#' N <- ncol(normforest)
#' M <- nrow(normforest)
#' Sds <- nrow(peaks) / ((N - 1) * (M - 1))
findpeaks <- function(x) {
  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)

  # center values, indices, and coordinates
  centers <- getValues(x)
  xcoords <- coordinates(x)[, 1]
  ycoords <- coordinates(x)[, 2]

  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M, byrow = TRUE)

  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)

  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  xcoords <- xcoords[-rm_inds]
  ycoords <- ycoords[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]

  # gather surrounding points
  xmin <- rows - 1
  xmax <- rows + 1
  ymin <- cols - 1
  ymax <- cols + 1
  ind <- seq(1, length(centers))
  surrounding <- lapply(ind, function(i) {zmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})

  # check for peak requirements
  check <- lapply(ind, function(i) {((sum(centers[i] > surrounding[[i]]) ==
                                        length(surrounding[[i]])) & centers[i] > 0)})

  # create dataframe, limit to actual peaks
  peaks <- data.frame(x = xcoords, y = ycoords, val = centers, ind = ind,
                      row = rows, col = cols, check = unlist(check))
  peaks <- peaks[peaks$check ==  TRUE,]
  peaks <- peaks[, 1:6]

  return(peaks)
}

#' Find local valleys.
#'
#' Locates local valleys on a raster. A valley is defined as any pixel where
#' all 8 surrounding pixels have higher values, and the center pixel
#' has a negative value.
#'
#' @param x An n x n raster object.
#' @return A dataframe of local valley locations (\code{x, y}) and
#'   values (\code{val}). The raster location index (\code{ind}),
#'   row (\code{row}), and column (\code{col}) are also listed.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # locate peaks and valleys
#' peaks <- findpeaks(normforest)
#' valleys <- findvalleys(normforest)
#'
#' # find top 5 peaks, valleys
#' top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
#' bottom_valleys <- valleys[order(valleys$val)[1:5],]
#'
#' # calculate the ten-point height
#' S10z <- (sum(top_peaks$val) + sum(abs(bottom_valleys$val))) / 5
findvalleys <- function(x) {
  N <- dim(x)[1] # rows
  M <- dim(x)[2] # cols

  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)

  # center values, indices, and coordinates
  centers <- getValues(x)
  ind <- seq(1, length(centers))
  xcoords <- coordinates(x)[, 1]
  ycoords <- coordinates(x)[, 2]

  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M)

  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)

  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  xcoords <- xcoords[-rm_inds]
  ycoords <- ycoords[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]

  # gather surrounding points
  xmin <- rows - 1
  xmax <- rows + 1
  ymin <- cols - 1
  ymax <- cols + 1
  ind <- seq(1, length(centers))
  surrounding <- lapply(ind, function(i) {zmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})

  # check for valley requirements
  check <- lapply(ind, function(i) {((sum(centers[i] < surrounding[[i]]) ==
                                        length(surrounding[[i]])) & centers[i] < 0)})

  # create dataframe, limit to actual valleys
  valleys <- data.frame(x = xcoords, y = ycoords, val = centers, ind = ind,
                      row = rows, col = cols, check = unlist(check))
  valleys <- valleys[valleys$check ==  TRUE,]
  valleys <- valleys[, 1:6]

  return(valleys)
}

#' Mean summit curvature.
#'
#' Calculates the mean summit curvature of a raster. Mean summit
#' curvature is the average principle curvature of local maximas
#' on the surface.
#'
#' @param x An n x n raster object.
#' @return A numeric value representing the average curvature of
#'   surface peaks.
#' @examples
#' # import raster image
#' data(normforest)
#'
#' # calculate mean summit curvature
#' Ssc <- ssc(normforest)
ssc <- function(x) {
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(x)
  xcoords <- coordinates(x)[, 1]
  ycoords <- coordinates(x)[, 2]
  deltax <- res(x)[1] / mean(res(x)[1], res(x)[2])
  deltay <- res(x)[2] / mean(res(x)[1], res(x)[2])

  # matrix of values
  zmat <- matrix(((z - min(z)) / (max(z) - min(z))), nrow = N, ncol = M, byrow = TRUE)

  # number of peaks
  peaks <- findpeaks(x)
  n <- nrow(peaks)

  # zshift of 1
  z_xpl <- zshift(x, xdist = 1, ydist = 0, scale = TRUE)
  z_xpl <- matrix(z_xpl, nrow = nrow(x), ncol = ncol(x) - 1, byrow = TRUE)

  # yshift of 1
  z_ypl <- zshift(x, xdist = 0, ydist = 1, scale = TRUE)
  z_ypl <- matrix(z_ypl, nrow = nrow(x) - 1, ncol = ncol(x), byrow = TRUE)

  # zshift of 1
  z_xmn <- zshift(x, xdist = -1, ydist = 0, scale = TRUE)
  z_xmn <- matrix(z_xmn, nrow = nrow(x), ncol = ncol(x))

  # zshift of 1
  z_ymn <- zshift(x, xdist = 0, ydist = -1, scale = TRUE)
  z_ymn <- matrix(z_ymn, nrow = nrow(x), ncol = ncol(x))

  # get z_xpl, z_ypl at peaks (add to df)
  peaks$val_xpl <- unlist(lapply(seq(1, nrow(peaks)),
                                 FUN = function(i) z_xpl[peaks$row[i], peaks$col[i]]))
  peaks$val_ypl <- unlist(lapply(seq(1, nrow(peaks)),
                                 FUN = function(i) z_ypl[peaks$row[i], peaks$col[i]]))
  peaks$val_xmn <- unlist(lapply(seq(1, nrow(peaks)),
                                 FUN = function(i) z_xmn[peaks$row[i], peaks$col[i]]))
  peaks$val_ymn <- unlist(lapply(seq(1, nrow(peaks)),
                                 FUN = function(i) z_ymn[peaks$row[i], peaks$col[i]]))

  # new center val
  peaks$val <- unlist(lapply(seq(1, nrow(peaks)),
                             FUN = function(i) zmat[peaks$row[i], peaks$col[i]]))

  # calculate curvature with normalized values
  ssc <- -(1 / (2 * n)) * sum(((((peaks$val - peaks$val_xpl) + (peaks$val - peaks$val_xmn)) / 2) / deltax),
             ((((peaks$val - peaks$val_ypl) + (peaks$val - peaks$val_ymn)) / 2) / deltay))

  return(ssc)
}

# summit density = number of local peaks per area
sds <- function(x) {
  M <- nrow(x)
  N <- ncol(x)

  peaks <- findpeaks(x)

  val <- nrow(peaks) / ((N - 1) * (M - 1))

  return(val)
}

# ten-point height = avg. height above mean surface for five highest local maxima plus avg.
# height below for five lowest local minima
s10z <- function(x) {
  peaks <- findpeaks(x)
  valleys <- findvalleys(x)

  # find top 5 peaks, valleys
  top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
  bottom_valleys <- valleys[order(valleys$val)[1:5],]

  val <- (sum(top_peaks$val, na.rm = TRUE) + sum(abs(bottom_valleys$val, na.rm = TRUE))) / 5

  return(val)
}
