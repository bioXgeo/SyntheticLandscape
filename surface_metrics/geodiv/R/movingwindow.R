#### !!!! srw and the like should not be done on very small windows (e.g., size = 3)

#' Calculate Texture Metrics per Pixel
#'
#' Calculates the various texture metrics over windows centered
#' on individual pixels. This creates a continuous surface of the
#' texture metric.
#'
#' @param x A raster.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels on each
#' side for square windows (must be an odd value), or distance from
#' center (in meters) for circular windows.
#' @param epsg_proj Numeric. Appropriate equal area EPSG code used to
#' crop raster to each circular window.
#' @param metric Character. Metric to calculate for each window. Metrics
#' are listed below.
#' @param parallel Logical. Option to run the calculations in
#' parallel on available cores.
#' @param ncores Numeric. If parallel is TRUE, number of cores on which to
#' run the calculations. Defaults to all available, minus 1.
#' @return A raster with pixel values representative of the metric
#' value for the window surrounding that pixel.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # get a surface of root mean square roughness
#' sq_img <- texture_image(x = normforest, window = 'square',
#' size = 11, epsg_proj = 5070, metric = 'sq', parallel = TRUE)
#'
#' # plot the result
#' plot(sq_img)
#' @export
texture_image <- function(x, window_type = 'square', size = 11, epsg_proj = 5070, metric, threshold = NULL,
                          low = NULL, high = NULL, parallel = TRUE, ncores = NULL){

  # data frame of x, y locations
  coords <- data.frame(xyFromCell(x, 1:ncell(x)))

  # get list of # total pixels
  pixlist <- seq(1, length(x), 1)

  # output raster
  out <- x

  if (parallel == FALSE) {
    for (i in pixlist) {
      cat('Pixel #: ', i, '\n')
      pt_coords <- coords[i, ]
      rownum <- rowFromCell(x, i)
      colnum <- colFromCell(x, i)

      if (window_type == 'square') {
        out[rownum, colnum] <- window_metric(x, 'square', size, epsg_proj = epsg_proj, rownum, colnum, metric, threshold, low, high)
      } else {
        out[rownum, colnum] <- window_metric(x, 'circle', size, epsg_proj = epsg_proj, rownum, colnum, metric, threshold, low, high)
      }
    }
  } else {
    if(missing(ncores)) {ncores <- parallel::detectCores() - 1}

    # make and start cluster
    try(snow::stopCluster(cl), silent = TRUE)
    cl <- snow::makeCluster(ncores)
    snow::makeSOCKcluster(cl)
    snow::clusterExport(cl = cl, list = list('x', 'out', 'coords', 'size',
                                             'window_type', 'epsg_proj',
                                             'metric', 'threshold', 'low', 'high'))
    invisible(snow::clusterEvalQ(cl, {
      library(raster)
      library(devtools)
      library(sf)
      devtools::load_all()}))
    result <- snow::parLapply(cl, pixlist[1:50], function(i) {
      pt_coords <- coords[i, ]
      rownum <- rowFromCell(x, i)
      colnum <- colFromCell(x, i)

      if (window_type == 'square') {
        outval <- window_metric(x, 'square', size, epsg_proj = epsg_proj, rownum, colnum, metric, threshold, low, high)
      } else {
        outval <- window_metric(x, 'circle', size, epsg_proj = epsg_proj, rownum, colnum, metric, threshold, low, high)
      }
      return(outval)
    })
    stopCluster(cl)
  }

  out <- setValues(out, unlist(result))

  return(out)
}

#' Calculate Texture Metric for Single Pixel
#'
#' Calculates the various texture metrics over a window centered
#' on an individual pixel.
#'
#' @param x A raster.
#' @param window_type Character. Type of window, either circular or square.
#' @param size Numeric. Size of window, in number of pixels on each
#' side for square windows (must be an odd value), or distance from
#' center (in meters) for circular windows.
#' @param epsg_proj Numeric. Appropriate equal area EPSG code used to
#' crop raster to each circular window.
#' @param rownum Numeric. Row number of pixel.
#' @param colnum Numeric. Column number of pixel.
#' @param metric Character. Metric to calculate for each window. Metrics
#' are listed below.
#' @param threshold Numeric. Value of autocorrelation distance (0 - 1), if
#' calculating \code{scl} or \code{str}.
#' @param low Numeric. Low value (0 - 1) if calculating \code{sdc}.
#' @param high Numeric. High value (0 - 1) if calculating \code{sdc}.
#' @return A raster with pixel values representative of the metric
#' value for the window surrounding that pixel.
#' @examples
#' library(raster)
#'
#' # import raster image
#' data(normforest)
#'
#' # get a surface of root mean square roughness
#' sq_img <- texture_image(x = normforest, window = 'circle',
#' size = 90, epsg_proj = 5070, metric = 'sq')
#'
#' # plot the result
#' plot(sq_img)
#' @export
window_metric <- function(x, window_type = 'square', size = 11, epsg_proj = 5070,
                          rownum, colnum, metric, threshold = NULL, low = NULL, high = NULL) {
  if(class(x) != 'RasterLayer') {stop('x must be a raster.')}
  if(class(window_type) != 'character') {stop('window_type must be a string.')}
  if(class(size) != 'numeric') {stop('size must be numeric.')}
  if(class(epsg_proj) != 'numeric') {stop('epsg_proj must be numeric.')}
  if(class(rownum) != 'numeric' & class(rownum) != 'integer') {stop('rownum must be numeric or integer.')}
  if(class(colnum) != 'numeric' & class(colnum) != 'integer') {stop('colnum must be numeric or integer.')}
  if(class(metric) != 'character') {stop('metric must be a character.')}
  if(!is.null(threshold) & class(threshold) != 'numeric') {stop('threshold must be numeric.')}
  if(!is.null(low) & class(low) != 'numeric') {stop('low must be numeric.')}
  if(!is.null(high) & class(high) != 'numeric') {stop('high must be numeric.')}

  if(!is.null(low) & is.null(high)) {stop('high value is required if low value is given.')}
  if(!is.null(high) & is.null(low)) {stop('high value is required if low value is given.')}
  if(!is.null(threshold) & length(threshold) > 1) {stop('too many values provided to threshold.')}
  if(length(rownum) > 1) {stop('too many values provided to rownum.')}
  if(length(colnum) > 1) {stop('too many values provided to colnum.')}
  if(length(metric) > 1) {stop('too many values provided for metric.')}
  if(length(epsg_proj) > 1) {stop('too many values provided to epsg_proj.')}
  if(length(size) > 1) {stop('too many values provided to size.')}
  if(length(window_type) > 1) {stop('too many values provided to window_type.')}

  if(!(metric %in% c('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk_adj',
                     'ssk', 'sku_exc', 'sku', 'sds', 'sfd', 'srw', 'srwi', 'shw', 'std',
                     'stdi', 'svi', 'str', 'ssc', 'sv', 'sph', 'sk', 'smean', 'spk', 'svk',
                     'scl', 'sdc'))) {stop('invalid metric argument.')}

  if (window_type == 'square') {
    # change size to distance out from center
    size <- floor(size / 2)

    # continue values to edges to account for edge effect (# pixels radius/edge)

    # first, get edge values that will be extended
    firstrow_vals <- x[1, ]
    firstcol_vals <- x[, 1]
    lastrow_vals <- x[nrow(x), ]
    lastcol_vals <- x[, ncol(x)]

    # add pixels on all sides, increasing the extent of the raster as well
    ext_x <- x
    dim(ext_x) <- c(nrow(x) + (2 * size), ncol(x) + (2 * size))
    extra_x <- size * res(x)[2]
    extra_y <- size * res(x)[1]
    extent(ext_x) <- extent(c(xmin(ext_x) - extra_x, xmax(ext_x) + extra_x,
                              ymin(ext_x) - extra_y, ymax(ext_x) + extra_y))

    # fill in top rows
    ext_x[1:size, (size + 1):(ncol(ext_x) - size)] <- rep(firstrow_vals, size)
    # fill in bottom rows
    ext_x[(nrow(ext_x) - (size - 1)):nrow(ext_x), (size + 1):(ncol(ext_x) - size)] <- rep(lastrow_vals, size)
    # fill in left columns
    ext_x[(size + 1):(nrow(ext_x) - size), 1:size] <- t(firstcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in right columns
    ext_x[(size + 1):(nrow(ext_x) - size), (ncol(ext_x) - (size - 1)):ncol(ext_x)] <- t(lastcol_vals * matrix(1, nrow = nrow(x), ncol = size))
    # fill in middle
    ext_x[(size + 1):(nrow(ext_x) - size), (size + 1):(ncol(ext_x) - size)] <- getValues(x)
    # fill in corners with nearest point value (always the same)
    ext_x_mat <- zoo::na.approx(matrix(ext_x, nrow = nrow(ext_x), ncol = ncol(ext_x)), rule = 2)
    ext_x <- setValues(ext_x, t(ext_x_mat))

    # crop to square
    newrow <- rownum + size
    newcol <- colnum + size
    cropped_x <- crop(ext_x, extent(ext_x, newrow - size, newrow + size, newcol - size, newcol + size))
  } else {

    # convert to new proj
    projx <- projectRaster(x, crs = sp::CRS(sf::st_crs(epsg_proj)$proj4string))

    # get equivalent # pixels of size
    pixeq_size <- ceiling(size / res(projx))[1]

    # extend...
    # continue values to edges to account for edge effect (# pixels radius/edge)
    firstrow_vals <- x[1, ]
    firstcol_vals <- x[, 1]
    lastrow_vals <- x[nrow(x), ]
    lastcol_vals <- x[, ncol(x)]

    # add pixels on all sides, increasing the extent of the raster as well
    ext_x <- x
    dim(ext_x) <- c(nrow(x) + (2 * pixeq_size), ncol(x) + (2 * pixeq_size))
    extra_x <- pixeq_size * res(x)[2]
    extra_y <- pixeq_size * res(x)[1]
    extent(ext_x) <- extent(c(xmin(ext_x) - extra_x, xmax(ext_x) + extra_x,
                              ymin(ext_x) - extra_y, ymax(ext_x) + extra_y))

    # fill in top rows
    ext_x[1:pixeq_size, (pixeq_size + 1):(ncol(ext_x) - pixeq_size)] <- rep(firstrow_vals, pixeq_size)
    # fill in bottom rows
    ext_x[(nrow(ext_x) - (pixeq_size - 1)):nrow(ext_x), (pixeq_size + 1):(ncol(ext_x) - pixeq_size)] <- rep(lastrow_vals, pixeq_size)
    # fill in left columns
    ext_x[(pixeq_size + 1):(nrow(ext_x) - pixeq_size), 1:pixeq_size] <- t(firstcol_vals * matrix(1, nrow = nrow(x), ncol = pixeq_size))
    # fill in right columns
    ext_x[(pixeq_size + 1):(nrow(ext_x) - pixeq_size), (ncol(ext_x) - (pixeq_size - 1)):ncol(ext_x)] <- t(lastcol_vals * matrix(1, nrow = nrow(x), ncol = pixeq_size))
    # fill in middle
    ext_x[(pixeq_size + 1):(nrow(ext_x) - pixeq_size), (pixeq_size + 1):(ncol(ext_x) - pixeq_size)] <- getValues(x)
    # fill in corners with nearest point value (always the same)
    ext_x_mat <- zoo::na.approx(matrix(ext_x, nrow = nrow(ext_x), ncol = ncol(ext_x)), rule = 2)
    ext_x <- setValues(ext_x, t(ext_x_mat))

    # get point coordinates
    pt_ind <- cellFromRowCol(x, rownum, colnum)
    coords <- data.frame(xyFromCell(x, 1:ncell(x)))
    pt_coords <- coords[pt_ind, ]

    # crop to circle
    pt_sf <- sf::st_as_sf(pt_coords, coords = c("x", "y"), crs = sf::st_crs(x)) %>%
      sf::st_transform(epsg_proj)
    poly_circ <- sf::st_buffer(pt_sf, size)
    poly_circ <- sf::st_transform(poly_circ, sf::st_crs(x))
    poly_circ <- sf::as_Spatial(poly_circ)
    cropped_x <- crop(ext_x, poly_circ)
    cropped_x <- mask(cropped_x, poly_circ)
  }

  # calculate metric
  if (metric == 'sa') {
    outval <- sa(cropped_x)
  }
  if (metric == 'sq') {
    outval <- sq(cropped_x)
  }
  if (metric == 's10z') {
    outval <- s10z(cropped_x)
  }
  if (metric == 'sdq') {
    outval <- sdq(cropped_x)
  }
  if (metric == 'sdq6') {
    outval <- sdq6(cropped_x)
  }
  if (metric == 'sdr') {
    outval <- sdr(cropped_x)
  }
  if (metric == 'sbi') {
    outval <- sbi(cropped_x)
  }
  if (metric == 'sci') {
    outval <- sci(cropped_x)
  }
  if (metric == 'ssk_adj') {
    outval <- ssk(cropped_x, adj = TRUE)
  }
  if (metric == 'ssk') {
    outval <- ssk(cropped_x, adj = FALSE)
  }
  if (metric == 'sku_exc') {
    outval <- sku(cropped_x, excess = TRUE)
  }
  if (metric == 'sku') {
    outval <- sku(cropped_x, excess = FALSE)
  }
  if (metric == 'sds') {
    outval <- sds(cropped_x)
  }
  if (metric == 'sfd') {
    outval <- sfd(cropped_x)
  }
  if (metric == 'srw') {
    outval <- srw(cropped_x, plot = FALSE)[[1]]
  }
  if (metric == 'srwi') {
    outval <- srw(cropped_x, plot = FALSE)[[2]]
  }
  if (metric == 'shw') {
    outval <- srw(cropped_x, plot = FALSE)[[3]]
  }
  if (metric == 'std') {
    outval <- mean(std(cropped_x, plot = FALSE)[[1]], na.rm = TRUE)
  }
  if (metric == 'stdi') {
    outval <- std(cropped_x, plot = FALSE)[[2]]
  }
  if (metric == 'svi') {
    outval <- svi(cropped_x)
  }
  if (metric == 'str') {
    outval <- str(cropped_x, threshold = threshold)
  }
  if (metric == 'ssc') {
    outval <- ssc(cropped_x)
  }
  if (metric == 'sv') {
    outval <- sv(cropped_x)
  }
  if (metric == 'sph') {
    outval <- sph(cropped_x)
  }
  if (metric == 'sk') {
    outval <- sk(cropped_x)
  }
  if (metric == 'smean') {
    outval <- smean(cropped_x)
  }
  if (metric == 'spk') {
    outval <- spk(cropped_x)
  }
  if (metric == 'svk') {
    outval <- svk(cropped_x)
  }
  if (metric == 'scl') {
    outval <- scl(cropped_x, threshold = threshold, plot = FALSE)
  }
  if (metric == 'sdc') {
    outval <- sdc(cropped_x, low = low, high = high)
  }
  return(outval)
}
