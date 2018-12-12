#' Radian to degree conversion.
#'
#' Converts radian value(s) to degrees.
#'
#' @param x Numeric. Radian value(s).
#' @return Numeric of degree value(s).
rad2deg <- function(x) {(x * 180) / (pi)}

#' Degree to radian conversion.
#'
#' Converts degree value(s) to radians.
#'
#' @param x Numeric. Degree value(s).
#' @return Numeric of degree value(s).
deg2rad <- function(x) {(x * pi) / (180)}

#' Calculate a least squares polynomial plane.
#'
#' Fits a polynomial plane of order \code{n} to a raster
#' image.
#'
#' @param x A raster.
#' @param order Numeric. Indicates the polynomial order to be fit.
#' @return A matrix of the same size as the raster with values
#'   predicted from the polynomial fit.
#' @examples
#' # import raster image
#' data(orforest)
#'
#' # find the 2nd order least squares polynomial plane
#' polyfit <- fitplane(orforest, order = 2)
#'
#' # create raster of polyfit
#' x <- setValues(orforest, polyfit)
#'
#' # plot the fit
#' plot(x)
fitplane <- function(x, order) {
  # extract coordinates and values
  xcoord <- coordinates(x)[, 1]
  ycoord <- coordinates(x)[, 2]
  z <- getValues(x)

  # fit least squares polynomial with order = order
  surfmod <- surf.ls(np = order, xcoord, ycoord, z)

  # predict polynomial model over raster
  surfvals <- matrix(predict(surfmod, xcoord, ycoord), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  return(surfvals)
}

#' Finds the best fit polynomial plane.
#'
#' Finds the best fit polynomial plane for a raster image. This
#' function tests least squares polynomial fits with orders of
#' 1 - 3 and determines which order minimizes the error when the
#' fit is subtracted from the original image.
#'
#' @param x A raster.
#' @return A raster of the same size as the input with values
#'   predicted from the best polynomial fit.
#' @examples
#' # import raster image
#' data(orforest)
#'
#' # find the least squares polynomial plane
#' poly <- bestfitplane(orforest)
#'
#' # plot the fit
#' plot(x)
bestfitplane <- function(x) {
  # fit least squares plane for polynomials from orders 1 - 3
  mods <- lapply(seq(1, 3), FUN = function(i) fitplane(x, order = i))

  # convert raster to matrix
  xmat <- matrix(x, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  # calculate raster errors from best fit planes for each order tested
  errlist <- lapply(seq(1, 3), FUN = function(i) xmat - mods[[i]])
  meanerr <- lapply(seq(1, 3), FUN = function(i) mean(errlist[[i]], na.rm = TRUE))

  # determine which order is best
  bestfit <- which(as.numeric(meanerr) == min(as.numeric(meanerr), na.rm = TRUE))

  # if more than one, get first
  if (length(bestfit) > 1) {
    bestfit <- min(bestfit, na.rm = TRUE)
  }

  # fill in raster with best fit values
  bfx <- setValues(x, errlist[[bestfit]])

  print(paste('Order of polynomial that minimizes errors: ', bestfit, sep = ''))
  return(bfx)
}