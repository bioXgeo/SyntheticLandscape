### Surface metric calculations for synthetic geo paper

# Written by ACS 3 Nov 2018
# Last edited by ACS 3 Nov 2018

### helpful resources
# https://www.keyence.com/ss/products/microscope/roughness/surface/parameters.jsp
# https://link-springer-com.proxy2.cl.msu.edu/content/pdf/10.1007%2F978-1-84800-297-5_22.pdf
# https://www.usna.edu/Users/physics/vanhoy/_files/SP425/LabDocs/STM/SPIP%204_4_3%20R/SPIP_Manual_4_2.pdf

# load packages -----------------------------------------------------------

library(raster)
library(tibble)

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
mod <- ecdf(z)

# get cumulative values
newy <- environment(mod)$x
newx <- environment(mod)$y

# flip axes so that value is on the y-axis and frequency is on the x-axis
plot(newy ~ newx, ylim = (rev(range(newy)))) # approximation of the bearing area curve

# find 40% of curve with least decline
# get each 40% segment, pick one with smallest slope
even_x <- seq(0, 1, length.out = 100000)
even_y <- quantile(mod, probs = even_x)
forty_length <- 0.4 * length(even_x)
i <- 1
while (i <= (length(even_x) - forty_length)) {
  y1 <- even_y[[i]]
  y2 <- even_y[[i + forty_length]]
  x1 <- even_x[[i]]
  x2 <- even_x[[i + forty_length]]
  if (i == 1) {
    slopes <- data.frame(slope = (y2 - y1) / (x2 - x1), xstart = x1, xend = x2)
  } else {
    slopes <- rbind(slopes, data.frame(slope = (y2 - y1) / (x2 - x1), xstart = x1, xend = x2))
  }
  
  i <- i + 1
}
# x value of start of 40% section with smallest decline
slope_min <- slopes[slopes$slope == min(slopes$slope),]

# calculate least-squares line for 40% of curve with smallest decline (lowest slope)
lm_data <- data.frame(x = newx[newx >= slope_min$xstart & newx <= slope_min$xend],
                      y = newy[newx >= slope_min$xstart & newx <= slope_min$xend])
ls_line <- lm(y ~ x, data = lm_data)
plot(newy ~ newx, ylim = (rev(range(newy)))) # approximation of the bearing area curve
abline(ls_line$coefficients[[1]], ls_line$coefficients[[2]], col = 'blue')

# get value of ls line between 0 and 1
pred_data <- remove_rownames(data.frame(x = even_x, y = even_y))
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
z05 <- quantile(mod, probs = c(0.05))[[1]]
z80 <- quantile(mod, probs = c(0.80))[[1]]

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
Sbi <- Sq/z05

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
Svk <- abs(quantile(mod, probs = 1) - quantile(mod, probs = Smr2))[[1]]

# reduced peak height = height of upper left triangle in abbott curve
Spk <- abs(quantile(mod, probs = 0) - quantile(mod, probs = Smr1))[[1]]

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
