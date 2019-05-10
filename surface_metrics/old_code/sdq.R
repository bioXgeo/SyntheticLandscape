# functions to find Sdq and Sdq6 using two-point and seven-point slopes

# both functions from the equations here:
# https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf

sdq <- function(rast) {
  # convert to equal area for resolution
  aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' 
  rast <- projectRaster(rast, crs = aea_crs) # note that this changes the # pixels (b/c rotation)
  
  # get dimensions
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  deltax <- res(rast)[1]
  deltay <- res(rast)[2]
  
  # create matrix of centers to get surrounding from
  zmat <- matrix(z, nrow = N, ncol = M, byrow = TRUE)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)
  
  # get rid of edge points
  rm_inds <- which(rows == max(rows) | cols == max(cols))
  z <- z[-rm_inds]
  x <- x[-rm_inds]
  y <- y[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]
  
  # for every point, get z of x + 1, z of y + 1
  xmax <- rows + 1
  ymax <- cols + 1
  ind <- seq(1, length(z))
  z_xplus <- unlist(lapply(ind, function(i) {zmat[rows[i], ymax[i]]}))
  z_yplus <- unlist(lapply(ind, function(i) {zmat[xmax[i], cols[i]]}))
  
  # calculate two-point slope
  sdq <- sqrt(sum((((z - z_xplus) / deltax) ^ 2) + 
                    (((z - z_yplus) / deltay) ^ 2), 
                  na.rm = TRUE) / length(z))
  
  return(sdq)
}

sdq6 <- function(rast) {
  # convert to equal area for resolution
  aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' 
  rast <- projectRaster(rast, crs = aea_crs) # note that this changes the # pixels (b/c rotation)
  
  # get dimensions
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  deltax <- res(rast)[1]
  deltay <- res(rast)[2]
  
  # create matrix of centers to get surrounding from
  zmat <- matrix(z, nrow = N, ncol = M, byrow = TRUE)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)
  
  # get rid of edge points
  rm_inds <- which(rows < 4 | rows > (max(rows) - 3) | cols < 4 | cols > (max(cols) - 3))
  z <- z[-rm_inds]
  x <- x[-rm_inds]
  y <- y[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]
  
  # for every point, get z of x + 1, x + 2, x + 3, x - 1, x - 2, x - 3
  yps1 <- cols + 1
  yps2 <- cols + 2
  yps3 <- cols + 3
  ymn1 <- cols - 1
  ymn2 <- cols - 2
  ymn3 <- cols - 3
  ind <- seq(1, length(z))
  z_xps1 <- unlist(lapply(ind, function(i) {zmat[rows[i], yps1[i]]}))
  z_xps2 <- unlist(lapply(ind, function(i) {zmat[rows[i], yps2[i]]}))
  z_xps3 <- unlist(lapply(ind, function(i) {zmat[rows[i], yps3[i]]}))
  z_xmn1 <- unlist(lapply(ind, function(i) {zmat[rows[i], ymn1[i]]}))
  z_xmn2 <- unlist(lapply(ind, function(i) {zmat[rows[i], ymn2[i]]}))
  z_xmn3 <- unlist(lapply(ind, function(i) {zmat[rows[i], ymn3[i]]}))
  
  # same as above, but for y
  xps1 <- rows + 1
  xps2 <- rows + 2
  xps3 <- rows + 3
  xmn1 <- rows - 1
  xmn2 <- rows - 2
  xmn3 <- rows - 3
  ind <- seq(1, length(z))
  z_yps1 <- unlist(lapply(ind, function(i) {zmat[xps1[i], cols[i]]}))
  z_yps2 <- unlist(lapply(ind, function(i) {zmat[xps2[i], cols[i]]}))
  z_yps3 <- unlist(lapply(ind, function(i) {zmat[xps3[i], cols[i]]}))
  z_ymn1 <- unlist(lapply(ind, function(i) {zmat[xmn1[i], cols[i]]}))
  z_ymn2 <- unlist(lapply(ind, function(i) {zmat[xmn2[i], cols[i]]}))
  z_ymn3 <- unlist(lapply(ind, function(i) {zmat[xmn3[i], cols[i]]}))
  
  # calculate two-point slope
  pklx <- (1 / (60 * deltax)) * (z_xps3 - (9 * z_xps2) + (45 * z_xps1) - (45 * z_xmn1) -
                                   (9 * z_xps2) - z_xmn3)
  pkly <- (1 / (60 * deltay)) * (z_yps3 - (9 * z_yps2) + (45 * z_yps1) - (45 * z_ymn1) -
                                   (9 * z_yps2) - z_ymn3)
  
  # calculate sdq6
  sdq6 <- sqrt((1 / length(z)) * sum(sum(pklx ^ 2, na.rm = TRUE), 
                                sum(pkly ^ 2, na.rm = TRUE), 
                                na.rm = TRUE))
  
  return(sdq6)
}