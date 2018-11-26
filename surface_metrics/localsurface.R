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
  xmat <- matrix(centers, nrow = N, ncol = M, byrow = TRUE)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(1:M, each = N)
  
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
  surrounding <- lapply(ind, function(i) {xmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})
  
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
  xmat <- matrix(centers, nrow = N, ncol = M)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(1:M, each = N)
  
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
  surrounding <- lapply(ind, function(i) {xmat[xmin[i]:xmax[i], ymin[i]:ymax[i]][-5]})
  
  # check for peak requirements
  check <- lapply(ind, function(i) {((sum(centers[i] < surrounding[[i]]) == 
                                        length(surrounding[[i]])) & centers[i] < 0)})
  
  # create dataframe, limit to actual peaks
  valleys <- data.frame(x = x, y = y, val = centers, ind = ind, 
                      row = rows, col = cols, check = unlist(check))
  valleys <- valleys[valleys$check ==  TRUE,]
  
  return(valleys)
}
