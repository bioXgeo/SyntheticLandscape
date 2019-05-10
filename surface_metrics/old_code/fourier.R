### fourier functions
std <- function(z, rast, plot = TRUE) {
  
  # get raster dimensions
  M <- ncol(rast)
  N <- nrow(rast)
  
  # get fourier transform
  zmat <- matrix(z, ncol = M, nrow = N, byrow = TRUE)
  # complex spectrum from fast fourier transform
  ft <- fft(zmat)
  ft_shift <- fftshift(ft)
  
  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))
  
  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(rast, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))
  
  # get origin of image (actually bottom center)
  origin <- c(mean(coordinates(amp_img)[,1]), ymin(amp_img))
  
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
  lines <- SpatialLines(linelist, proj4string = CRS(proj4string(amp_img)))
  
  # plot and calculate amplitude sums along rays
  if(plot == TRUE) {
    plot(amp_img)
    lines(lines)
  }
  
  Aalpha <- list()
  for (i in 1:length(lines)) {
    Aalpha[i] <- extract(amp_img, lines[i], fun = sum)
  }
  
  std <- rad2deg(alpha[which(unlist(Aalpha) == max(unlist(Aalpha), na.rm = TRUE))])
  stdi <- mean(unlist(Aalpha), na.rm = TRUE) / max(unlist(Aalpha), na.rm = TRUE)
  
  return(list(std, stdi))
}

srw <- function(z, rast, plot = TRUE) {
  
  # get raster dimensions
  M <- ncol(rast)
  N <- nrow(rast)
  
  # get fourier transform
  zmat <- matrix(z, ncol = M, nrow = N, byrow = TRUE)
  # complex spectrum from fast fourier transform
  ft <- fft(zmat)
  ft_shift <- fftshift(ft)
  
  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))
  
  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(rast, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))
  
  # get origin of image (actually bottom center)
  origin <- c(mean(coordinates(amp_img)[,1]), ymin(amp_img))
  
  # calculate half circles extending from origin
  nv <- 100
  angle.inc <- 2 * pi / nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  radius <- seq(0, half_dist, half_dist / 80)
  linex <- unlist(lapply(seq(1, length(radius)), function(x) origin[1] + radius[x] * cos(angles)))
  liney <- unlist(lapply(seq(1, length(radius)), function(x) origin[2] + radius[x] * sin(angles)))
  linelist <- lapply(seq(1, length(linex), 100), 
                     FUN = function(i) Lines(list(Line(cbind(linex[i:(i + 99)], liney[i:(i + 99)]))), ID = paste('p', i, sep = '')))
  lines <- SpatialLines(linelist, proj4string = CRS(proj4string(amp_img)))
  
  # plot and get amplitude sums within each radius
  if (plot == TRUE){
    plot(amp_img)
    plot(lines, add = TRUE)
  }
  Br <- list()
  Br[1] <- 0
  for (i in 2:length(lines)) {
    templine <- crop(lines[i], extent(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))
    Br[i] <- extract(amp_img, templine, fun = sum)
  }
  if (plot == TRUE) {
    plot(unlist(Br) ~ (1 / radius), type = 'l')
  }
  Srw <- radius[which(unlist(Br) == max(unlist(Br), na.rm = TRUE))]
  Srwi <- mean(unlist(Br), na.rm = TRUE) / max(unlist(Br), na.rm = TRUE)
  
  half_val <- 0.5 * sum(unlist(Br), na.rm = TRUE)
  vals <- as.numeric()
  for (i in 1:length(unlist(Br))) {
    vals[i] <- sum(unlist(Br)[1:i], na.rm = TRUE)
  }
  shw_ind <- vals - half_val
  Shw <- radius[which(abs(shw_ind) == min(abs(shw_ind), na.rm = TRUE))]
  
  return(list(Srw, Srwi, Shw))
}

sfd <- function(z, rast, x, y) {
  # this has been checked against the matlab version and produces the same results
  # now need to confirm that it is the correct value...
  
  # get raster dimensions
  M <- ncol(rast)
  N <- nrow(rast)
  
  # get fourier transform
  zmat <- matrix(z, ncol = M, nrow = N, byrow = TRUE)
  # complex spectrum from fast fourier transform
  ft <- fft(zmat)
  ft_shift <- fftshift(ft)
  
  # amplitude spectrum
  amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))
  
  # take amplitude image, cut in half (y direction)
  amp_img <- setValues(rast, amplitude)
  half_dist <- (ymax(amp_img) - ymin(amp_img)) / 2
  ymin <- ymax(amp_img) - half_dist
  amp_img <- crop(amp_img, c(xmin(amp_img), xmax(amp_img), ymin, ymax(amp_img)))
  
  # get origin of image (actually bottom center)
  origin <- c(mean(coordinates(amp_img)[,1]), ymin(amp_img))
  
  # from matlab fdsurfft function
  num_dir <- 24 # number of directions that the frequency space is uniformally divided
  num_rad <- 30 # number of points that the radius is uniformally divided

  M <- nrow(ft_shift)
  N <- ncol(ft_shift)
  xdist <- matrix(x, nrow = M, ncol = N, byrow = TRUE)[1, ] - origin[1]
  ydist <- matrix(y, nrow = M, ncol = N, byrow = TRUE)[, 1] - origin[2]
  xctr <- which(xdist == min(abs(xdist))) # should actually be index of x that = origin x
  yctr <- which(ydist == min(abs(ydist))) # same as xctr
  fim <- ft_shift
  
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