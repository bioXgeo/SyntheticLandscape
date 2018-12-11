### calculates different offsets of z matrix, negative or positive
zshift <- function(zmat, xdist = 0, ydist = 0, xrm, yrm) {
  # xdist is distance away in x direction
  # ydist is distance away in y direction
  # xrm, yrm can be specified if you want matrix clipped more than xdist/ydist
  # xrm, yrm cannot be less than corresponding dist, or of a different sign
  # need to provide at least one of xdist or ydist
  try(if(missing(xrm)) (xrm = xdist))
  try(if(missing(yrm)) (yrm = ydist))
  
  # get dimensions
  N <- dim(zmat)[1] # rows
  M <- dim(zmat)[2] # cols
  
  # define z
  z <- as.numeric(zmat)
  
  # row/col of each center
  rows <- rep(1:N, each = M)
  cols <- rep(rep(1:M), N)
  
  # get rid of edge points x distance away
  rm_inds <- numeric(0)
  if (xrm > 0) {
    posx_rm <- which(cols > (max(cols) - xrm))
    rm_inds <- c(rm_inds, posx_rm)
  } else if (xrm < 0) {
    negx_rm <- which(cols < (abs(xrm) + 1))
    rm_inds <- c(rm_inds, negx_rm)
  } else {
    posx_rm <- NULL
    negx_rm <- NULL
    rm_inds <- c(rm_inds, posx_rm, negx_rm)
  }
  
  if (yrm > 0) {
    posy_rm <- which(rows > (max(rows) - yrm))
    rm_inds <- c(rm_inds, posy_rm)
  } else if (yrm < 0) {
    negy_rm <- which(rows < (abs(yrm) + 1))
    rm_inds <- c(rm_inds, negy_rm)
  } else {
    posy_rm <- NULL
    negy_rm <- NULL
    rm_inds <- c(rm_inds, posy_rm, negy_rm)
  }
  
  z <- z[-rm_inds]
  x <- x[-rm_inds]
  y <- y[-rm_inds]
  rows <- rows[-rm_inds]
  cols <- cols[-rm_inds]
  
  # for every point, get z of x + 1, x + 2, x + 3, x - 1, x - 2, x - 3
  yshift <- rows + ydist
  xshift <- cols + xdist
  ind <- seq(1, length(z))
  if (xdist != 0 & ydist != 0) {
    z_shift <- unlist(lapply(ind, function(i) {zmat[yshift[i], xshift[i]]}))
  } else if (xdist != 0 & ydist == 0){
    z_shift <- unlist(lapply(ind, function(i) {zmat[rows[i], xshift[i]]}))
  } else if (xdist == 0 & ydist != 0) {
    z_shift <- unlist(lapply(ind, function(i) {zmat[yshift[i], cols[i]]}))
  } else {
    z_shift = unlist(lapply(ind, function(i) {zmat[rows[i], cols[i]]}))
  }
  
  if (xdist < 0) {
    z_shift <- matrix(z_shift, nrow = length(unique(yshift)), 
                      ncol = length(unique(xshift)), byrow = TRUE)
    z_shift <- cbind(rep(rep(NA, abs(xdist)), nrow(z_shift)), z_shift)
    z_shift <- as.numeric(z_shift)
  } 
  if (ydist < 0) {
    z_shift <- matrix(z_shift, nrow = length(unique(yshift)), 
                      ncol = length(unique(xshift)), byrow = TRUE)
    z_shift <- rbind(rep(rep(NA, abs(ydist)), ncol(z_shift)), z_shift)
    z_shift <- as.numeric(z_shift)
  } 
  
  return(z_shift)
}