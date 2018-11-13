### Surface metric calculations for synthetic geo paper

# Written by ACS 3 Nov 2018
# Last edited by ACS 12 Nov 2018

### helpful resources
# https://www.keyence.com/ss/products/microscope/roughness/surface/parameters.jsp
# https://link-springer-com.proxy2.cl.msu.edu/content/pdf/10.1007%2F978-1-84800-297-5_22.pdf
# https://www.usna.edu/Users/physics/vanhoy/_files/SP425/LabDocs/STM/SPIP%204_4_3%20R/SPIP_Manual_4_2.pdf

# also see below for other abbott-firestone curve definition
# ftp://ftp.astmtmc.cmu.edu/docs/diesel/mack/minutes/2013/Mack.2013-02-07.Meeting/21763%20pdf.pdf

# probably need to flip min/max values to get correct curve -- done 11/12/18

# see https://trac.osgeo.org/grass/browser/grass/trunk/raster/r.surf.area/area.c 
# for surface area (from GRASS GIS)

# https://stackoverflow.com/questions/40376299/r-fft-fourier-spectrum-of-image

# load packages -----------------------------------------------------------

library(raster)
library(tibble)
library(spatialEco)

# functions ---------------------------------------------------------------

# simpson's rule area under curve for empirical functions
simpsons_above <- function(model, a, b, n = 100) {
  # determines area above curve (just within rect defined by area below)
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)
  y <- quantile(model, probs = x) # get y-values of inverse cdf function
  s <- ((b - a) / (3 * n)) * (y[[1]] + 
                                sum(4 * y[seq(2, n - 1, by = 2)]) + 
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])
  # get inverse of s for actual area above curve
  area_above <- ((max(y) - min(y)) * (max(x) - min(x))) - s
  return(area_above)
}

simpsons_below <- function(model, a, b, n = 100) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)
  y <- quantile(model, probs = x) # get y-values of inverse cdf function
  s <- ((b - a) / (3 * n)) * (y[[1]] + 
                                sum(4 * y[seq(2, n - 1, by = 2)]) + 
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])
  return(s)
}

# import testing data -----------------------------------------------------

# ndvi raster over oregon
# any small raster will do for now
rast <- raster('/home/annie/Documents/satellite-field-comp/data/images/p45_r30/summer_ndvi_p45_r30_2000_2016_30m.tif')

# crop for example
newext <- c(-123, -122.9, 43, 43.1)
newrast <- crop(rast, newext)

# get values of raster
vals <- getValues(newrast)

# calculate best-fit plane and variables ------------------------------------

# everything is relative to the best-fitting plane, so calculate the formula
# for that first.
# ax + by + c = z ---> least squares fitting of plane

# fitting a plane only changes things very slightly for most rasters
A <- as.matrix(data.frame(coordinates(newrast)[, 1],
                    coordinates(newrast)[, 2],
                    rep(1, nrow(coordinates(newrast)))))
b <- as.matrix(vals, ncol = 1)
fit <- (solve((t(A) %*% A)) %*% t(A)) %*% b
errors <- b - (A %*% fit)

# plot new raster
newrast2 <- newrast
newrast2 <- setValues(newrast2, errors)
plot(newrast2)

# basic values
N <- length(errors)
s <- sd(errors)
z <- errors
zbar <- mean(errors)
x <- coordinates(newrast2)[, 1]
y <- coordinates(newrast2)[, 2]

# calculate cdf and other necessary components ----------------------------

# calculate cumulative probability density function of surface 'height' (= ndvi)
mod <- ecdf((1 - z))

# get cumulative values
newy <- (1 - environment(mod)$x)
newx <- environment(mod)$y

# flip axes so that value is on the y-axis and frequency is on the x-axis
plot(newy ~ newx)#, ylim = (rev(range(newy)))) # approximation of the bearing area curve

# find 40% of curve with least decline
# use symmetric difference quotient to estimate the derivative at evenly spaced points
# then find 40% consecutive section with lowest mean slope
even_x <- seq(0, 1, length.out = 100000)
even_y <- (1 - quantile(mod, probs = even_x))
forty_length <- 0.4 * length(even_x)
h <- 0.001
i <- 1 # can't do difference quotient at very end
while (i <= length(even_x)) {
  if (i < 101) {
    # newton's difference quotient at near end
    fxh_pos <- (1 - quantile(mod, probs = (even_x[i] + h))[[1]])
    fxh_neg <- (1 - quantile(mod, probs = (even_x[i]))[[1]])
  } else if (i > 99900) {
    # variation on newton's difference quotient at far end
    fxh_pos <- (1 - quantile(mod, probs = (even_x[i]))[[1]])
    fxh_neg <- (1 - quantile(mod, probs = (even_x[i] - h))[[1]])
  } else {
    # symmetric difference quotient everywhere else (99800 points)
    fxh_pos <- (1 - quantile(mod, probs = (even_x[i] + h))[[1]])
    fxh_neg <- (1 - quantile(mod, probs = (even_x[i] - h))[[1]])
  }
  diff_quo <- (fxh_pos - fxh_neg) / (2 * h)
  if (i == 1) {
    slopes <- data.frame(slope = diff_quo, x = even_x[i])
  } else {
    slopes <- rbind(slopes, data.frame(slope = diff_quo, x = even_x[i]))
  }
  
  i <- i + 1
}

# calculate averages for each 40% segment
for (i in 1:(length(even_x) - forty_length)) {
  if (i == 1) {
    means <- data.frame(slope = mean(slopes$slope[i:(i + forty_length)]), 
                        xstart = slopes$x[i], xend = slopes$x[(i + forty_length)])
  } else {
    means <- rbind(means, data.frame(slope = mean(slopes$slope[i:(i + forty_length)]), 
                        xstart = slopes$x[i], xend = slopes$x[(i + forty_length)]))
  }
}

# x value of start of 40% section with smallest decline
slope_min <- means[means$slope == min(means$slope),]

# calculate least-squares line for 40% of curve with smallest decline (lowest slope)
lm_data <- data.frame(x = newx[newx >= slope_min$xstart & newx <= slope_min$xend],
                      y = (1 - newy[newx >= slope_min$xstart & newx <= slope_min$xend]))
ls_line <- lm(y ~ x, data = lm_data)
plot(newy ~ newx) # approximation of the bearing area curve
abline(ls_line$coefficients[[1]], ls_line$coefficients[[2]], col = 'blue')

# get value of ls line between 0 and 1
pred_data <- remove_rownames(data.frame(x = even_x, y = (1 - even_y)))
pred_data$y <- predict(ls_line, newdata = pred_data)

# what is the ls line y-value at x = 0, x = 1?
ls_int_high <- pred_data$y[pred_data$x == 0]
ls_int_low <- pred_data$y[pred_data$x == 1]
abline(h = ls_int_high, col = 'red')
abline(h = ls_int_low, col = 'red')

# Smr1/Smr2 = x values that correspond to cdf y values at ls_int_high/low 
Smr1 <- mod(ls_int_high)
Smr2 <- mod(ls_int_low)
abline(v = Smr1, col = 'green')
abline(v = Smr2, col = 'green')

# calculate height at various percentiles of z (see fig.2 in Kedron et al. 2018 paper)
z05 <- (1 - quantile(mod, probs = c(0.05))[[1]])
z80 <- (1 - quantile(mod, probs = c(0.80))[[1]])

# basic surface metrics ---------------------------------------------------

# average roughness = average absolute deviation of surface heights from mean
Sa <- sum(abs(z - zbar)) / N

# root mean square roughness = standard deviation of surface heights
Sq <- s

# surface skewness = asymmetry of surface height distribution
Ssk <- (sum((z - zbar)^3) / N) / (s^3) # Fisher-Pearson coefficient of skewness
Ssk_adj <- (sqrt((N * (N - 1))) / (N - 2)) * Ssk # adjusted Fisher-Pearson coefficient of skewness

# surface kurtosis = peaked-ness of surface distribution
Sku <- (sum((z - zbar)^4) / N) / (s^4) 
Sku_exc <- Sku - 3 # excess kurtosis (i.e., diff from normal distribution kurtosis)

# maximum valley depth = lowest value in the landscape
Sv <- abs(min(z))

# maximum peak height = highest value in the landscape
Sp <- abs(max(z))

# mean peak height = average peak height
Smean <- mean(z[z > 0])

# abbott curve variables --------------------------------------------------

# surface bearing index = ratio of Sq to height from top of surface to height at 
# 5% of bearing area (z05)
Sbi <- Sq / z05

# valley fluid retention index = void volume (area under Abbott curve) in 'valley' zone
# see fig.2a from Kedron et al. (2018)
Svi <- simpsons_above(model = mod, b = 1, a = 0.8, n = 500)

# core fluid retention index = void volume (area above Abbott curve) in the core zone
# see fig.2a from Kedron et al. (2018)
core_above <- simpsons_above(model = mod, b = 1, a = 0.05, n = 1000)
Sci <- core_above - Svi

# core roughness depth = height difference between the intersection points of the found least
# mean square line in the Abbott curve (see fig. 2b from Kedron et al. 2018)
Sk <- abs(ls_int_high - ls_int_low)

# reduced valley depth = height of triangle drawn at 100% on Abbott curve 
Svk <- abs((1 - quantile(mod, probs = 1)) - (1 - quantile(mod, probs = Smr2)))[[1]]

# reduced peak height = height of upper left triangle in abbott curve
Spk <- abs((1 - quantile(mod, probs = 0)) - (1 - quantile(mod, probs = Smr1)))[[1]]

# calculate local variables -----------------------------------------------
# from Image Metrology:
# local minimum = points where all 8 surrounding points are higher & below zero
# local maximum = points where all 8 surrounding points are lower & above zero

# local peaks (maxima) and valleys (minima) -- coordinates and values
# this takes a while (should parallelize!)
N <- dim(newrast2)[1]
M <- dim(newrast2)[2]
peaks <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)
valleys <- data.frame(x = NA, y = NA, val = NA, ind = NA, row = NA, col = NA)
for (i in 2:(N - 1)) {
  for (j in 2:(M - 1)) {
    center <- newrast2[i, j]
    surrounding <- newrast2[(i - 1):(i + 1), (j - 1):(j + 1)]
    surrounding <- surrounding[-5]
    peak <- ((sum(center > surrounding) == length(surrounding)) & center > 0) # has to be positive
    valley <- ((sum(center < surrounding) == length(surrounding)) & center < 0) # has to be negative
    
    # get raster index
    ind <- (((i - 1) * M) + j)
    
    if (peak == TRUE) {
      ind <- (((i - 1) * M) + j)
      peaks <- rbind(peaks, data.frame(x = x[ind], y = y[ind], val = center, ind = ind, row = i, col = j))
    } 
    if (valley == TRUE) {
      valleys <- rbind(valleys, data.frame(x = x[ind], y = y[ind], val = center, ind = ind, row = i, col = j))
    }
  }
}

# clean up peak/valley dfs
peaks <- peaks[-1,]
valleys <- valleys[-1,]

# summit density = number of local peaks per area
Sds <- nrow(peaks) / ((N - 1) * (M - 1))

# ten-point height = avg. height above mean surface for five highest local maxima plus avg. 
# height below for five lowest local minima
top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
bottom_valleys <- valleys[order(valleys$val)[1:5],]
S10z <- (sum(top_peaks$val) + sum(abs(bottom_valleys$val))) / 5

# fourier variables -------------------------------------------------------

# get fourier transform
zmat <- matrix(z, ncol = M, nrow = N, byrow = TRUE)
ft <- fft(zmat)
