### Surface metric calculations for synthetic geo paper

# Written by ACS 3 Nov 2018
# Last edited by ACS 26 Nov 2018

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

# best resource for fourier and surface area/slope variables
# https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf

# changes needed ----------------------------------------------------------

# should really create a function to create various shifts in x and y directions,
# since this is used a lot in several of the other functions

# load packages -----------------------------------------------------------

library(raster)
library(tibble)
library(spatialEco)
library(dplyr)
library(plotly)
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/zshift.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/fftshift.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/simpsons.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/bestfit.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/localsurface.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/sdq.R')
source('/home/annie/Documents/SyntheticLandscape/surface_metrics/surfacearea.R')

# functions ---------------------------------------------------------------

# radian/degree conversions
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

# import testing data -----------------------------------------------------

# ndvi raster over oregon
# any small raster will do for now
rast <- raster('/home/annie/Documents/satellite-field-comp/data/images/p45_r30/summer_ndvi_p45_r30_2000_2016_30m.tif')

# crop for example
newext <- c(-123, -122.9, 43, 43.1)
newrast <- crop(rast, newext)

# get values of raster
vals <- getValues(newrast)

# set scale (meters here)
xscale <- 30
yscale <- 30

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
slopes <- slopecalc(even_x, h) # calculate slope at every point
means <- slopemeans(slopes) # calculate averages for each 40% segment

# x value of start of 40% section with smallest decline
slope_min <- means[means$slope == min(means$slope),]

# calculate least-squares line for 40% of curve with smallest decline (lowest slope)
lm_data <- data.frame(x = newx[newx >= slope_min$xstart & newx <= slope_min$xend],
                      y = newy[newx >= slope_min$xstart & newx <= slope_min$xend])
ls_line <- lm(y ~ x, data = lm_data)
plot(newy ~ newx) # approximation of the bearing area curve
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
Smr1 <- mod(1 - ls_int_high)
Smr2 <- mod(1 - ls_int_low)
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

# re-define N/M to raster
N <- nrow(newrast2)
M <- ncol(newrast2)

# local peaks (maxima) and valleys (minima) -- coordinates and values
peaks <- findpeaks(newrast2)
valleys <- findvalleys(newrast2)

# find top 5 peaks, valleys
top_peaks <- peaks[order(peaks$val, decreasing = TRUE)[1:5],]
bottom_valleys <- valleys[order(valleys$val)[1:5],]

# summit density = number of local peaks per area
Sds <- nrow(peaks) / ((N - 1) * (M - 1))

# ten-point height = avg. height above mean surface for five highest local maxima plus avg. 
# height below for five lowest local minima
S10z <- (sum(top_peaks$val) + sum(abs(bottom_valleys$val))) / 5

# root mean square slope = variance in local slope across the surface
# note: different than terrain result
Sdq <- sdq(newrast2)

# area root mean square slope = similar to Sdq but includes more neighbor pixels in slope computation
Sdq6 <- sdq6(newrast2)

# surface area ratio = ratio between surface area to area of flat plane with same x,y dimensions
Sdr <- sdr(newrast2)

# mean summit curvature = average principal curvature of local maximas on the surface
Ssc <- ssc(newrast2, peaks)

# fourier variables -------------------------------------------------------

# get fourier transform
zmat <- matrix(z, ncol = M, nrow = N, byrow = TRUE)
# complex spectrum from fast fourier transform
ft <- fft(zmat)
ft_shift <- fftshift(ft)

# amplitude and phase spectrum
amplitude <- sqrt((Re(ft_shift) ^ 2) + (Im(ft_shift) ^ 2))
phase <- atan(Im(ft_shift) / Re(ft_shift))
power <- (Im(ft_shift) ^ 2) + (Im(ft_shift) ^ 2)

amplitude <- (sqrt(Re(ft) ^ 2) + (Im(ft) ^ 2))

# get log-log frequency-amplitude
loglog_data <- data.frame(amp = log(matrix(amplitude, ncol = 1)), freq = log(seq(1, length(amplitude))))
loglog <- lm(amp ~ freq, data = loglog_data)

### TESTING - HERE DOWN

# plot amplit
test <- newrast2
test <- setValues(test, log(amplitude)) # plot amplitude
test2 <- setValues(test, phase)
# Each pixel in the Fourier transform has a coordinate (h,k) representing 
# the contribution of the sine wave with x-frequency h, and y-frequency k 
# in the Fourier transform.

# view amplitude
p <- plot_ly(showscale = FALSE) %>%
  add_surface(z = ~log(amplitude))
p

# check if all angles exist...
summary(matrix(rad2deg(phase), ncol = 1)) # symmmetric 90 to 90, so yes

# dominant texture direction = angle of dominating texture in image calculated from
# fourier spectrum
Std <- rad2deg(phase[which(amplitude == max(amplitude))])

# texture direction index = relative dominance of Std over other directions of texture, 
# defined as avg. amplitude sum over all directions divided by amplitude sum of 
# dominating direction
Stdi <- mean(amplitude) / sum(amplitude[rad2deg(phase) == Std])

# fractal dimension = calculated for different angles of angular spectrum by analyzing 
# fourier amplitude spectrum
Sfd <- 2 - loglog$coefficients[2]
