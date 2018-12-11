# figure out peaks and valleys of surface

# from Image Metrology:
# local minimum = points where all 8 surrounding points are higher & below zero
# local maximum = points where all 8 surrounding points are lower & above zero

findpeaks <- function(rast) {
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)
  
  # center values, indices, and coordinates
  centers <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  
  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M, byrow = TRUE)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)
  
  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  x <- x[-rm_inds]
  y <- y[-rm_inds]
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
  peaks <- data.frame(x = x, y = y, val = centers, ind = ind, 
                      row = rows, col = cols, check = unlist(check))
  peaks <- peaks[peaks$check ==  TRUE,]
  
  return(peaks)
}

findvalleys <- function(rast) {
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)
  
  # center values, indices, and coordinates
  centers <- getValues(rast)
  ind <- seq(1, length(centers))
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  
  # create matrix of centers to get surrounding from
  zmat <- matrix(centers, nrow = N, ncol = M)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)
  
  # get rid of edge points
  rm_inds <- which(rows < 2 | rows == max(rows) | cols < 2 | cols == max(cols))
  centers <- centers[-rm_inds]
  x <- x[-rm_inds]
  y <- y[-rm_inds]
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
  check <- lapply(ind, function(i) {((sum(centers[i] < surrounding[[i]]) == 
                                        length(surrounding[[i]])) & centers[i] < 0)})
  
  # create dataframe, limit to actual peaks
  valleys <- data.frame(x = x, y = y, val = centers, ind = ind, 
                      row = rows, col = cols, check = unlist(check))
  valleys <- valleys[valleys$check ==  TRUE,]
  
  return(valleys)
}

ssc <- function(rast, peaks) {
  # calculates mean summit curvature of peaks
  
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  deltax <- res(rast)[1] / mean(res(rast)[1], res(rast)[2])
  deltay <- res(rast)[2] / mean(res(rast)[1], res(rast)[2])
  
  # matrix of values
  zmat <- matrix(((z - min(z)) / (max(z) - min(z))), nrow = N, ncol = M, byrow = TRUE)
  
  # number of peaks
  n <- nrow(peaks)
  
  # zshift of 1
  z_xpl <- zshift(zmat, xdist = 1, ydist = 0)
  z_xpl <- matrix(z_xpl, nrow = nrow(rast), ncol = ncol(rast) - 1, byrow = TRUE)
  
  # yshift of 1
  z_ypl <- zshift(zmat, xdist = 0, ydist = 1)
  z_ypl <- matrix(z_ypl, nrow = nrow(rast) - 1, ncol = ncol(rast), byrow = TRUE)
  
  # zshift of 1
  z_xmn <- zshift(zmat, xdist = -1, ydist = 0)
  z_xmn <- matrix(z_xmn, nrow = nrow(rast), ncol = ncol(rast))
  
  # zshift of 1
  z_ymn <- zshift(zmat, xdist = 0, ydist = -1)
  z_ymn <- matrix(z_ymn, nrow = nrow(rast), ncol = ncol(rast))
  
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