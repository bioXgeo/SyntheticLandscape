### fourier functions
sfd <- function(fftshift, origin, x, y) {
  # this has been checked against the matlab version and produces the same results
  # now need to confirm that it is the correct value...
  
  # from matlab fdsurfft function
  num_dir <- 24 # number of directions that the frequency space is uniformally divided
  num_rad <- 30 # number of points that the radius is uniformally divided

  M <- nrow(fftshift)
  N <- ncol(fftshift)
  xdist <- matrix(x, nrow = M, ncol = N, byrow = TRUE)[1, ] - origin[1]
  ydist <- matrix(y, nrow = M, ncol = N, byrow = TRUE)[, 1] - origin[2]
  xctr <- which(xdist == min(abs(xdist))) # should actually be index of x that = origin x
  yctr <- which(ydist == min(abs(ydist))) # same as xctr
  fim <- fftshift
  
  # calculate power spectrum
  mag <- Re(log(fim * fim + 10 ^ (-6)))
  sumBrite <- zeros(num_dir, num_rad) # accumulation magnitude for each direction and radius
  nCount <- zeros(num_dir, num_rad) # number of magnitude
  radius <- zeros(2 * num_rad, 1) # accumulation magnitude for all directions 
  radCount <- zeros(2 * num_rad, 1) # number of magnitude for all directions
  
  # Compute phase image and phase histogram
  phaseim <- zeros(M, N)
  phase <- zeros(180)
  for (j  in 1:M) {
    for (i in 1:N) {
      realv <- Re(fim[j, i])
      imagv <- Im(fim[j, i])
      if (realv == 0) {
        value = pi / 2
      } else {
        value <- atan((imagv / realv))
        phaseim[j, i] <- value
        ang <- floor(180 * (pi / 2 + value) / pi)
      }
      if (ang < 0) {
        ang <- 0 
      }
      if (ang > 179) {
        ang <- 179
      }
      phase[ang + 1] <- phase[ang + 1] + 1
    }
  }
  
  maxphase <- max(phase)
  plot(phase / maxphase, type = 'l')
  
  # accumulation of magnitude for each direction and radius
  rmax <- log(min(M, N) / 2) # maximum radius (log scale)
  for (j in 1:M) {
    if (j != yctr) {
      yval <- yctr - j
      y2 <- yval * yval
    for (i in 1:N) {
      if (i != xctr) {
        xval <- i - xctr
        rho <- log(sqrt(y2 + xval * xval))
        if (rho > 0 & rho <= rmax) {
          mval <- mag[j, i]
          temp <- yval /xval
          theta <- atan(temp)
          if (xval < 0) {
            theta <- theta + pi
          }
          if (theta < 0) {
            theta <- theta + 2 * pi
          }
          ang <- floor(num_dir * theta / (2 * pi))
          if (ang > num_dir - 1 | ang < 0) {
            ang <- num_dir - 1 
          }
          k <- floor(2 * num_rad * rho / rmax)
          h <- floor(k / 2) 
          if (k > 2 * num_rad - 1) {
            h <- num_rad - 1
            k <- 2 * num_rad - 1
          }
          if (h >= 5) {
            sumBrite[ang + 1, h + 1] <- sumBrite[ang + 1, h + 1] + mval
            nCount[ang + 1, h + 1] <- nCount[ang + 1, h + 1] + 1
          }
          if (k >= 5) {
            radius[k + 1] <- radius[k + 1] + mval
            radCount[k + 1] <- radCount[k + 1] + 1
          }
        }
      }
    }
    }
  }

  # linear regression
  slope <- as.numeric()
  intercept <- as.numeric()
  for (ang in 1:num_dir) {
    sumx <- 0
    sumy <- 0
    sumx2 <- 0
    sumxy <- 0
    sumn <- 0
    for (range in 6:num_rad) {
      if (nCount[ang, range] > 0) {
        yval <- sumBrite[ang, range] / nCount[ang, range]
        xval <- (range -1) * rmax / num_rad 
        sumx <- sumx + xval
        sumy <- sumy + yval
        sumx2 <- sumx2 + xval * xval
        sumxy <- sumxy + xval * yval
        sumn <- sumn + 1
      }
    }
    slope[ang] <- (sumn * sumxy - sumx * sumy) / (sumn * sumx2 - sumx * sumx)
    intercept[ang] <- (sumy - slope[ang] * sumx) / sumn
  }
  
  # compute average slope over all directions and scales
  sumn <- 0
  yval <- as.numeric()
  tempr <- as.numeric()
  for (k in 6:(2 * num_rad)) {
    if (radCount[k] > 0) {
      sumn <- sumn + 1
      yval[sumn] <- radius[k] / radCount[k]
      tempr[sumn] <- (k - 1) * rmax / (2 * num_rad) # might be c(k, -1)
    }
  }
  model_data <- data.frame(yval = yval, tempr = tempr)
  p <- lm(yval ~ tempr, data = model_data)
  averslope <- p$coefficients[2]
  averIC <- p$coefficients[1]
  
  fitln <- predict(p, newdata = model_data)
  plot(tempr, yval, type = 'l')
  lines(tempr, fitln, col = 'blue')
  slope[num_dir + 1] <- slope[1]
  intercept[num_dir + 1] <- intercept[1]
  
  # draw rose plot of slope and intercept
  ang <- seq(1, (num_dir + 1))
  polar.plot(intercept[ang], rad2deg(pi / num_dir + (ang - 1) * 2 * pi / num_dir))
  polar.plot(abs(slope[ang]), rad2deg(pi / num_dir + (ang - 1) * 2 * pi / num_dir))

  return(abs(averslope)[[1]])
}