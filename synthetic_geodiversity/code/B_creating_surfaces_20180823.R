### making synthetic geodiversity rasters ####
library(raster)
library(gstat)
library(plotrix)
library(NLMR)
library(RColorBrewer)

setwd("~/Desktop/synthetic_geodiversity")

today <- format(Sys.Date(), "%Y%m%d")

# goal is to make rasters that are 144 x 144 km with 30 m pixels
# meaning 4800 pixels on a side. For simplicity values range from 0 to 100.

total.area <- matrix(1, nrow = 4800, ncol = 4800)
total.raster <- raster(total.area)
extent(total.raster) <- c(0, 144, 0, 144)

## CONTINUOUS ###
cont.values <- seq(from = 0, to = 100, length.out = 4800)
continuous.grid <- total.area * cont.values
continuous.surface <- raster(continuous.grid)
extent(continuous.surface) <- c(0, 144, 0, 144)

## STRIPES ##
stripe.values <- rep(c(rep(100, 100), rep(50, 100), rep(0, 100)), 16)
stripe.grid <- total.area * stripe.values
stripe.surface <- raster(t(stripe.grid))
extent(stripe.surface) <- c(0, 144, 0, 144)

## CHECKERBOARD ##

check1 <- rep(c(rep(100, 100), rep(50, 100), rep(0, 100)), 16)
check2 <- rep(c(rep(50, 100), rep(0, 100), rep(100, 100)), 16)
check3 <- rep(c(rep(0, 100), rep(100, 100), rep(50, 100)), 16)

short.area <- matrix(1, nrow = 4800, ncol = 100)

check1.area <- (short.area) * check1
check2.area <- (short.area) * check2
check3.area <- (short.area) * check3

check.grid <- cbind(check1.area, check2.area, check3.area)

check.grid.all <- cbind(check.grid, check.grid, check.grid, check.grid, check.grid,
                        check.grid, check.grid, check.grid, check.grid, check.grid,
                        check.grid, check.grid, check.grid, check.grid, check.grid,
                        check.grid)
check.surface <- raster(check.grid.all)
extent(check.surface) <- c(0, 144, 0, 144)

## HILL ## (Cribbed from Quentin)

sd_max <- 100
# sd gets smaller as you go away from origin
sd_distance <- function(d, sd0) sd0 - .0294 * d
# Create tiled surface of values with that sd

hill.area <- total.area

for (i in 1:4800) {
  for (j in 1:4800) {
    hill.area[i, j] <- sd_distance(sqrt((i - 2400) ^ 2 + (j - 2400) ^ 2), sd_max)
  }
  #print(i)
}

hill.surface <- raster(hill.area)
extent(hill.surface) <- c(0, 144, 0, 144)

### now make a small clusters gaussian field
# note I'm setting a seed here so this will be repeatable
date()
gf.nug20.range100 <- nlm_gaussianfield(ncol = ncol(total.raster),
                        nrow = nrow(total.raster),
                        res = res(total.raster)[1],
                        autocorr_range = 100,
                        nug = 20,
                        mag_var = 100,
                        mean = 50,
                        rescale = F,
                        user_seed = 7)
date()

#### now a smoother gaussian field
date()
gf.nug5.range500 <- nlm_gaussianfield(ncol = ncol(total.raster),
                                      nrow = nrow(total.raster),
                                      res = res(total.raster)[1],
                                      autocorr_range = 500,
                                      nug = 5,
                                      mag_var = 100,
                                      mean = 50,
                                      rescale = F,
                                      user_seed = 7)
date()


## now stack the surfaces and save!

all.surf <- stack(stripe.surface, check.surface, continuous.surface, 
                  hill.surface, gf.nug20.range100, gf.nug5.range500)

names(all.surf) <- c("stripe", "checker", "continuous", "hill", 
                     "gf_nug20_range100", "gf_nug5_range500")

writeRaster(all.surf, 
            paste0("data/synthetic_surfaces_30m_", today, ".tif"))

## now aggregate to 1000 m grain and save

empty.coarse.surf <- raster(matrix(1, nrow = 144, ncol = 144))
extent(empty.coarse.surf) <- c(0, 144, 0, 144)

all.surf.ag <- aggregate(all.surf, 
                         fact = 33.333, 
                         fun = mean, 
                         expand = TRUE, 
                         na.rm = TRUE)

all.surf.ag <- resample(all.surf.ag, 
                        empty.coarse.surf, 
                        method = "ngb")

writeRaster(all.surf.ag, 
            paste0("data/synthetic_surfaces_1000m_", today, ".tif"))

# this is figuring out how to pack in as many sampling circles as possible

sample.coords <- as.data.frame(matrix(NA, nrow = 12, ncol = 2))
names(sample.coords) <- c("x", "y")
sample.coords[1:3, 1] <- c(21, 62, 103)
sample.coords[1:3, 2] <- 20.5
sample.coords[4:6, 1] <- c(41.5, 82.5, 123.5)
sample.coords[4:6, 2] <- 54.9
sample.coords[7:9, 1] <- c(21, 62, 103)
sample.coords[7:9, 2] <- 89.25
sample.coords[10:12, 1] <- c(41.5, 82.5, 123.5)
sample.coords[10:12, 2] <- 123.6

write.csv(sample.coords, 
          paste0("data/sampling_centroids_x12_", today, ".csv"), 
          row.names = F)

plot(all.surf$continuous)
points(sample.coords, pch = 4)
for(i in 1:12) {
  draw.circle(sample.coords[i,1], sample.coords[i,2], 5, lty = 1, lwd = 0.5)
  draw.circle(sample.coords[i,1], sample.coords[i,2], 20, lty = 2, lwd = 1.5)
}



