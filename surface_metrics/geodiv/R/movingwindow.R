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
#' size = 11, epsg_proj = 5070, metric = 'sq')
#'
#' # plot the result
#' plot(sq_img)
#' @export
texture_image <- function(x, window_type = 'square', size = 3, epsg_proj = 5070, metric,
                          parallel = TRUE, ncores = NULL){

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
        out[rownum, colnum] <- window_metric(x, 'square', size, epsg_proj = 5070, rownum, colnum, metric)
      } else {
        out[rownum, colnum] <- window_metric(x, 'circle', size, epsg_proj = 5070, rownum, colnum, metric)
      }
    }
  } else {
    if(missing(ncores)) {ncores <- detectCores() - 1}

    # make and start cluster
    stopCluster(cl)
    cl <- makeCluster(ncores)
    registerDoSNOW(cl)
    snow::clusterExport(cl = cl, list = list('x', 'out', 'coords', 'size',
                                             'window_type', 'epsg_proj',
                                             'metric', 'window_metric'))
    clusterEvalQ(cl, {
      library(raster)
      library(devtools)
      library(sf)
      devtools::load_all()})
    result <- parLapply(cl, pixlist[1:500], function(i) {
      pt_coords <- coords[i, ]
      rownum <- rowFromCell(x, i)
      colnum <- colFromCell(x, i)

      if (window_type == 'square') {
        outval <- window_metric(x, 'square', size, epsg_proj = 5070, rownum, colnum, metric)
      } else {
        outval <- window_metric(x, 'circle', size, epsg_proj = 5070, rownum, colnum, metric)
      }
      return(outval)
    })
    stopCluster(cl)
  }
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
window_metric <- function(x, window = 'square', size = 3, epsg_proj = NULL,
                          rownum, colnum, metric) {
  if (window == 'square') {
    # change size to distance out from center
    size <- floor(size / 2)

    # extend...
    # continue values to edges to account for edge effect (# pixels radius/edge)
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
    projx <- projectRaster(x, crs = CRS(st_crs(epsg_proj)$proj4string))

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
    pt_sf <- st_as_sf(pt_coords, coords = c("x", "y"), crs = st_crs(x)) %>%
      st_transform(epsg_proj)
    poly_circ <- st_buffer(pt_sf, size)
    poly_circ <- st_transform(poly_circ, st_crs(x))
    poly_circ <- as(poly_circ, 'Spatial')
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
  if (metric == 'srw') {
    outval <- srw(cropped_x, plot = FALSE)[[1]]
  }
  return(outval)
}
