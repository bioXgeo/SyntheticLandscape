### extracting data and summarizing for synthetic geodiversity ####
library(raster)
library(rgdal)
library(plotrix)
library(e1071)
library(ggplot2)
library(reshape)
library(devtools)
devtools::load_all('/home/annie/Documents/SyntheticLandscape/surface_metrics/geodiv/')

# setwd and today's date for file writing
setwd("/home/annie/Documents/synthetic_geodiversity")
today <- format(Sys.Date(), "%Y%m%d")

# set data directory
data.dir <- "/home/annie/Documents/synthetic_geodiversity/data/"

# set layer names
geo.names <- c("Stripe", "Check", "Slope", "Hill", "G.Patchy", "G.Smooth",
               "AZ.therm", "KS.therm")
geo.names.short <- c("St", "Ch","Sl", "Hi", "GP", "GS", "AZ.TH", "KS.TH")

# read in synthetic data from step A
in.fine <- stack(paste0(data.dir, 
                        "synthetic_surfaces_30m_20190521.tif"))

# realized the gaussian patchy layer has some negative values, which could mess
# up the math later, so setting all those neg values to zero
#in.fine[[5]][in.fine[[5]] < 0] <- 0

# from annie: convert anything that isn't exactly 0-100 to that range
for (i in seq(4, 6)) {
  temp <- in.fine[[i]]
  
  gauss.min <- cellStats(temp, stat = 'min', na.rm = TRUE)
  gauss.max <- cellStats(temp, stat = 'max', na.rm = TRUE)
  gauss.diff <- gauss.max - gauss.min
  
  in.fine[[i]] <- ((temp - gauss.min) / gauss.diff) * 100
}

# read in real thermal data and remove crs and extent so it can be stacked with
# synthetic data
in.AZ.landsat <- raster(paste0(data.dir, 
                                "arizona/arizona_landsat_square_20190521.tif"))
crs(in.AZ.landsat) <- NA
extent(in.AZ.landsat) <- extent(in.fine)

in.KS.landsat <- raster(paste0(data.dir, 
                               "kansas/kansas_landsat_square_20190521.tif"))
crs(in.KS.landsat) <- NA
extent(in.KS.landsat) <- extent(in.fine)

all.fine <- stack(in.fine, in.AZ.landsat, in.KS.landsat)
names(all.fine) <- geo.names
unique.ids.fine <- paste0(geo.names, "30")

# now give it all a projection because some of the fxns below don't run without
# a projection
crs(all.fine) <- "+proj=utm +zone=12 +datum=WGS84 +units=km +no_defs"

rm(in.fine)

# write this so we can use again as needed
writeRaster(all.fine, paste0(data.dir, 
                             "all_fine_res_layers_",
                             today,
                             ".tif"),
            overwrite = TRUE)


# now the coarse data!
in.coarse <- stack(paste0(data.dir, 
                          "synthetic_surfaces_1000m_20190521.tif"))

# modis
in.AZ.modis <- raster(paste0(data.dir, 
                               "arizona/arizona_modis_square_20190521.tif"))
crs(in.AZ.modis) <- NA

# aaannnnnddd these don't line up perfectly, so cropping. :\ would be better to
# go back and fix code in step B.
#extent(in.AZ.modis) <- extent(c(0, 145, 0, 145))
#in.AZ.modis <- crop(in.AZ.modis, extent(in.coarse))

# alt from annie: grab inner-most square
in.AZ.modis <- crop(in.AZ.modis, extent(in.AZ.modis, 8, 151, 8, 151))
extent(in.AZ.modis) <- extent(c(0, 144, 0, 144))

in.KS.modis <- raster(paste0(data.dir, 
                             "kansas/kansas_modis_square_20190521.tif"))
crs(in.KS.modis) <- NA

#extent(in.KS.modis) <- extent(c(0, 145, 0, 145))
#in.KS.modis <- crop(in.KS.modis, extent(in.coarse))

# alt from annie: grab inner-most square
in.KS.modis <- crop(in.KS.modis, extent(in.KS.modis, 8, 151, 8, 151))
extent(in.KS.modis) <- extent(c(0, 144, 0, 144))

# now all together
all.coarse <- stack(in.coarse, in.AZ.modis, in.KS.modis)

# from annie: convert anything that isn't exactly 0-100 to that range
for (i in seq(3, 7)) {
  temp <- all.coarse[[i]]
  
  gauss.min <- cellStats(temp, stat = 'min', na.rm = TRUE)
  gauss.max <- cellStats(temp, stat = 'max', na.rm = TRUE)
  gauss.diff <- gauss.max - gauss.min
  
  all.coarse[[i]] <- ((temp - gauss.min) / gauss.diff) * 100
}

names(all.coarse) <- geo.names
unique.ids.coarse <- paste0(geo.names, "1k")
crs(all.coarse) <- "+proj=utm +zone=12 +datum=WGS84 +units=km +no_defs"

writeRaster(all.coarse, paste0(data.dir, 
                             "all_coarse_res_layers_",
                             today,
                             ".tif"),
            overwrite = TRUE)


# read in centroids csv
in.centers <- read.csv("./data/sampling_centroids_x12_20190521.csv",
                       stringsAsFactors = FALSE)

## this code generates circle shapefiles but only needs to run once, so is
## commented out now.

# # in order for draw.circles to work right you need to plot the image they go
# # on to first (this is a bit hacky and makes this run super slow - like a minute
# # on kyla's desktop)
# plot(in.fine[[1]])
# circ.5 <- draw.circle(in.centers[1,1], in.centers[1,2], 5, lty = 1, lwd = 0.5)
# circ.5.mat <- cbind(circ.5$x, circ.5$y)
# circ.5.poly <- Polygon(circ.5.mat)
# circ.5.polys <- Polygons(list(circ.5.poly), 1)
# circ.5.list <- list(circ.5.polys)
# 
# circ.20 <- draw.circle(in.centers[1,1], in.centers[1,2], 20, lty = 2, lwd = 1.5)
# circ.20.mat <- cbind(circ.20$x, circ.20$y)
# circ.20.poly <- Polygon(circ.20.mat)
# circ.20.polys <- Polygons(list(circ.20.poly), 1)
# circ.20.list <- list(circ.20.polys)
# 
# for(i in 2:12) {
#   #plot(in.fine[[1]]) # hopefully don't need this if prev lines run
#   circ.5 <- draw.circle(in.centers[i,1], in.centers[i,2], 5, lty = 1, lwd = 0.5)
#   circ.5.mat <- cbind(circ.5$x, circ.5$y)
#   circ.5.poly <- Polygon(circ.5.mat)
#   circ.5.polys <- Polygons(list(circ.5.poly), i)
#   circ.5.list <- c(circ.5.list, circ.5.polys)
#   
#   circ.20 <- draw.circle(in.centers[i,1], in.centers[i,2], 20, lty = 2, lwd = 1.5)
#   circ.20.mat <- cbind(circ.20$x, circ.20$y)
#   circ.20.poly <- Polygon(circ.20.mat)
#   circ.20.polys <- Polygons(list(circ.20.poly), i)
#   circ.20.list <- c(circ.20.list, circ.20.polys)
#   
#   print(i)
# }
# 
# # make spatial polygons data frames and write to shapefiles for fun
# circ.5.SP <- SpatialPolygons(circ.5.list)
# data.forSPDF <- as.data.frame(matrix(1:12, nrow = 12, ncol = 1))
# names(data.forSPDF) <- "plot.id"
# circ.5.SPDF <- SpatialPolygonsDataFrame(circ.5.SP, data = data.forSPDF)
# 
# writeOGR(circ.5.SPDF, 
#          dsn = "./data",
#          layer = paste0("sample_circles_5km_", today),
#          driver = "ESRI Shapefile")
# 
# circ.20.SP <- SpatialPolygons(circ.20.list)
# circ.20.SPDF <- SpatialPolygonsDataFrame(circ.20.SP, data = data.forSPDF)
# 
# writeOGR(circ.20.SPDF, 
#          dsn = "./data",
#          layer = paste0("sample_circles_20km_", today),
#          driver = "ESRI Shapefile")

# instead of re-running that we can just read those in now!
circ.5.SPDF <- readOGR(paste0(data.dir, 
                              "circles/sample_circles_5km_20181031.shp"))

circ.20.SPDF <- readOGR(paste0(data.dir, 
                              "circles/sample_circles_20km_20181031.shp"))



#### Notes: 
# 2) add in gradient metrics and write out
# 3) check out the extra metrics required
# 4) figure out how to implement and write out
# 5) clean code, summarize, write to github, hpc, and dropbox
# 6) send kyla link to dropbox with revised code and results



### OK! now let's extract some data and write it! ## 
##### 5 km radius first ############################
variables <- c("st.deviation",
               "kurtosis", 
               "skewness",
               "range",
               "sbi",
               "prof.convexity.mean",
               "prof.convexity.sd",
               "bolstad.curvature.mean",
               "bolstad.curvature.sd",
               "slope.mean",
               "slope.sd",
               "TPI.mean",
               "TPI.sd",
               "TRI.mean",
               "TRI.sd",
               "surf.area.sum",
               "aspect.sine.mean",
               "aspect.sine.sd",
               "roughness.mean",
               "roughness.sd",
               "text.dir",
               "text.dir.ind",
               "RWI.mean",
               "RWI.sd",
               "variogram.diff",
               "summit.density"
               )
out.fine <- as.data.frame(matrix(NA, 
                                 nrow = 12 * 8, 
                                 ncol = length(variables) + 4))
names(out.fine) <- c("surf.name", "surf.type", "fine.coarse", "poly", variables)

out.coarse <- out.fine

# first let's do the simple metrics (summary of original data) 
# this is super slow -> methinks its the extraction step
date()
for (y in 1:12) {
  in.circ <- extract(all.fine, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$surf.name[i] <- unique.ids.fine[x]
    out.fine$surf.type[i] <- names(all.fine)[x]
    out.fine$fine.coarse[i] <- "fine"
    out.fine$poly[i] <- y

    out.fine$st.deviation[i] <- sd(in.circ[,x], na.rm = TRUE)
    out.fine$kurtosis[i] <- kurtosis(in.circ[,x], na.rm = TRUE)
    out.fine$skewness[i] <- skewness(in.circ[,x], na.rm = TRUE)
    out.fine$range[i] <- max(in.circ[,x], na.rm = TRUE) - 
      min(in.circ[,x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ TEXTURE DIRECTION ###############################

# this takes a while.
for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.fine[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    stdvals <- std(rast)
    out.fine$text.dir[i] <- mean(stdvals[[1]], na.rm = TRUE)
    out.fine$text.dir.ind[i] <- stdvals[[2]]
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ SURFACE BEARING INDEX ############################

for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.fine[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    out.fine$sbi[i] <- sbi(rast)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ SUMMIT DENSITY ############################

for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.fine[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    out.fine$sds[i] <- sds(rast)
    cat('Done with layer ', x, 'polygon ', y, '.', '\n')
  }
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ CURVATURE ###############################

# note the spatialEco curvature functions only work on a single raster layer, hence
# the looping here...
all.fine.curv <- curvature(all.fine[[1]], 
                          s = 3, 
                          type = "bolstad")
for (a in 2:8) {
  gradient <- curvature(all.fine[[a]], 
                        s = 3, 
                        type = "bolstad")
  all.fine.curv <- stack(all.fine.curv, gradient)
  print(a)
}
names(all.fine.curv) <- paste0("curvature.", geo.names.short)

## CURVATURE to table##
for (y in 1:12) {
  in.circ <- extract(all.fine.curv, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$bolstad.curvature.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$bolstad.curvature.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ CONVEXITY ###############################

# note the spatialEco curvature functions only work on a single raster layer, hence
# the looping here...
all.fine.conv <- curvature(all.fine[[1]], 
                           s = 3, 
                           type = "profile")
for (a in 2:8) {
  gradient <- curvature(all.fine[[a]], 
                        s = 3, 
                        type = "profile")
  all.fine.conv <- stack(all.fine.conv, gradient)
  print(a)
}
names(all.fine.conv) <- paste0("profile.convexity.", geo.names.short)

## CURVATURE to table##
for (y in 1:12) {
  in.circ <- extract(all.fine.conv, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$prof.convexity.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$prof.convexity.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE SLOPE LAYERS ###############################
# note the gdal terrain functions only work on a single raster layer, hence
# the looping here...
all.fine.slope <- terrain(all.fine[[1]], 
                         opt = "slope", 
                         unit = "radians", 
                         neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.fine[[a]], 
                      opt = "slope", 
                      unit = "radians", 
                      neighbors = 8)
  all.fine.slope <- stack(all.fine.slope, gradient)
  print(a)
}
names(all.fine.slope) <- paste0("slope.", geo.names.short)

## SLOPE to table##
for (y in 1:12) {
  in.circ <- extract(all.fine.slope, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$slope.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$slope.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################# SURFACE AREA #####################################

all.fine.area <- (0.03 * 0.03) / cos(all.fine.slope)
names(all.fine.area) <- paste0("area.", geo.names.short)

## SURF AREA ##
for (y in 1:12) {
  in.circ <- extract(all.fine.area, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$surf.area.sum[i] <- sum(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE TPI LAYERS ###############################

all.fine.tpi <- terrain(all.fine[[1]], 
                         opt = "TPI", 
                         neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.fine[[a]], 
                      opt = "TPI", 
                      neighbors = 8)
  all.fine.tpi <- stack(all.fine.tpi, gradient)
  print(a)
}
names(all.fine.tpi) <- paste0("tpi.", geo.names.short)

## TPI ##
for (y in 1:12) {
  in.circ <- extract(all.fine.tpi, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$TPI.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$TPI.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE TRI LAYERS ###############################

all.fine.tri <- terrain(all.fine[[1]], 
                       opt = "TRI", 
                       neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.fine[[a]], 
                      opt = "TRI", 
                      neighbors = 8)
  all.fine.tri <- stack(all.fine.tri, gradient)
  print(a)
}
names(all.fine.tri) <- paste0("tri.", geo.names.short)

## TRI ##
for (y in 1:12) {
  in.circ <- extract(all.fine.tri, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.fine$TRI.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$TRI.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE ASPECT LAYERS ###############################

all.fine.asp <- terrain(all.fine[[1]], 
                       opt = "aspect", 
                       unit = "radians",
                       neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.fine[[a]], 
                      opt = "aspect", 
                      unit = "radians",
                      neighbors = 8)
  all.fine.asp <- stack(all.fine.asp, gradient)
  print(a)
}

all.fine.sine.asp <- sin(all.fine.asp)
names(all.fine.sine.asp) <- paste0("asp.", geo.names.short)

## ASPECT ##
for (y in 1:12) {
  in.circ <- extract(all.fine.sine.asp, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$aspect.sine.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$aspect.sine.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE ROUGHNESS LAYERS ###############################

all.fine.rough <- terrain(all.fine[[1]], 
                       opt = "roughness", 
                       neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.fine[[a]], 
                      opt = "roughness", 
                      neighbors = 8)
  all.fine.rough <- stack(all.fine.rough, gradient)
  print(a)
}
names(all.fine.rough) <- paste0("rough.", geo.names.short)


## ROUGHNESS ##
for (y in 1:12) {
  in.circ <- extract(all.fine.rough, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$roughness.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.fine$roughness.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

write.csv(out.fine, paste0("./data/fine_surf_metrics_ALL_", today, ".csv"))

################################################################################
############################## NOW COARSE ######################################
################################################################################

# first let's do the simple metrics (summary of original data) 
# this is super slow -> methinks its the extraction step
date()
for (y in 1:12) {
  in.circ <- extract(all.coarse, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$surf.name[i] <- unique.ids.coarse[x]
    out.coarse$surf.type[i] <- names(all.coarse)[x]
    out.coarse$fine.coarse[i] <- "coarse"
    out.coarse$poly[i] <- y
    
    out.coarse$st.deviation[i] <- sd(in.circ[, x], na.rm = TRUE)
    out.coarse$kurtosis[i] <- kurtosis(in.circ[, x], na.rm = TRUE)
    out.coarse$skewness[i] <- skewness(in.circ[, x], na.rm = TRUE)
    out.coarse$range[i] <- max(in.circ[, x], na.rm = TRUE) - 
      min(in.circ[,x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ TEXTURE DIRECTION ###############################

# this takes a while.
for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.coarse[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    stdvals <- std(rast)
    out.coarse$text.dir[i] <- mean(stdvals[[1]], na.rm = TRUE)
    out.coarse$text.dir.ind[i] <- stdvals[[2]]
    cat('Done with layer ', x, 'polygon ', y, '.', '\n')
  }
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ SURFACE BEARING INDEX ############################

for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.coarse[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    out.coarse$sbi[i] <- sbi(rast)
    cat('Done with layer ', x, 'polygon ', y, '.', '\n')
  }
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ SUMMIT DENSITY ############################

for (y in 1:12) {
  for (x in 1:8) {  
    rast <- all.coarse[[x]]
    rast <- crop(rast, circ.5.SPDF[y, ])
    rast <- mask(rast, circ.5.SPDF[y, ])
    i <- (y - 1) * 8 + x
    out.coarse$sds[i] <- sds(rast)
    cat('Done with layer ', x, 'polygon ', y, '.', '\n')
  }
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ CURVATURE ###############################

# note the spatialEco curvature functions only work on a single raster layer, hence
# the looping here...
all.coarse.curv <- curvature(all.coarse[[1]], 
                           s = 3, 
                           type = "bolstad")
for (a in 2:8) {
  gradient <- curvature(all.coarse[[a]], 
                        s = 3, 
                        type = "bolstad")
  all.coarse.curv <- stack(all.coarse.curv, gradient)
  print(a)
}
names(all.coarse.curv) <- paste0("curvature.", geo.names.short)

## CURVATURE to table##
for (y in 1:12) {
  in.circ <- extract(all.coarse.curv, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$bolstad.curvature.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$bolstad.curvature.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ CONVEXITY ###############################

# note the spatialEco curvature functions only work on a single raster layer, hence
# the looping here...
all.coarse.conv <- curvature(all.coarse[[1]], 
                           s = 3, 
                           type = "profile")
for (a in 2:8) {
  gradient <- curvature(all.coarse[[a]], 
                        s = 3, 
                        type = "profile")
  all.coarse.conv <- stack(all.coarse.conv, gradient)
  print(a)
}
names(all.coarse.conv) <- paste0("profile.convexity.", geo.names.short)

## CURVATURE to table##
for (y in 1:12) {
  in.circ <- extract(all.coarse.conv, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$prof.convexity.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$prof.convexity.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE SLOPE LAYERS ###############################
# note the gdal terrain functions only work on a single raster layer, hence
# the looping here...
all.coarse.slope <- terrain(all.coarse[[1]], 
                          opt = "slope", 
                          unit = "radians", 
                          neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.coarse[[a]], 
                      opt = "slope", 
                      unit = "radians", 
                      neighbors = 8)
  all.coarse.slope <- stack(all.coarse.slope, gradient)
  print(a)
}
names(all.coarse.slope) <- paste0("slope.", geo.names.short)

## SLOPE to table##
for (y in 1:12) {
  in.circ <- extract(all.coarse.slope, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$slope.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$slope.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################# SURFACE AREA #####################################

all.coarse.area <- (1*1)/cos(all.coarse.slope)
names(all.coarse.area) <- paste0("area.", geo.names.short)

## SURF AREA ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.area, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$surf.area.sum[i] <- sum(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE TPI LAYERS ###############################

all.coarse.tpi <- terrain(all.coarse[[1]], 
                        opt = "TPI", 
                        neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.coarse[[a]], 
                      opt = "TPI", 
                      neighbors = 8)
  all.coarse.tpi <- stack(all.coarse.tpi, gradient)
  print(a)
}
names(all.coarse.tpi) <- paste0("tpi.", geo.names.short)

## TPI ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.tpi, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$TPI.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$TPI.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE TRI LAYERS ###############################

all.coarse.tri <- terrain(all.coarse[[1]], 
                        opt = "TRI", 
                        neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.coarse[[a]], 
                      opt = "TRI", 
                      neighbors = 8)
  all.coarse.tri <- stack(all.coarse.tri, gradient)
  print(a)
}
names(all.coarse.tri) <- paste0("tri.", geo.names.short)

## TRI ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.tri, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$TRI.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$TRI.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE ASPECT LAYERS ###############################

all.coarse.asp <- terrain(all.coarse[[1]], 
                        opt = "aspect", 
                        unit = "radians",
                        neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.coarse[[a]], 
                      opt = "aspect", 
                      unit = "radians",
                      neighbors = 8)
  all.coarse.asp <- stack(all.coarse.asp, gradient)
  print(a)
}

all.coarse.sine.asp <- sin(all.coarse.asp)
names(all.coarse.sine.asp) <- paste0("asp.", geo.names.short)

## ASPECT ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.sine.asp, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$aspect.sine.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$aspect.sine.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################ MAKE ROUGHNESS LAYERS ###############################

all.coarse.rough <- terrain(all.coarse[[1]], 
                          opt = "roughness", 
                          neighbors = 8)
for (a in 2:8) {
  gradient <- terrain(all.coarse[[a]], 
                      opt = "roughness", 
                      neighbors = 8)
  all.coarse.rough <- stack(all.coarse.rough, gradient)
  print(a)
}
names(all.coarse.rough) <- paste0("rough.", geo.names.short)


## ROUGHNESS ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.rough, circ.5.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y - 1) * 8 + x
    out.coarse$roughness.mean[i] <- mean(in.circ[, x], na.rm = TRUE)
    out.coarse$roughness.sd[i] <- sd(in.circ[, x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

write.csv(out.coarse, paste0("./data/coarse_surf_metrics_ALL_", today, ".csv"))

# combine the two dataframes
out.all <- rbind.data.frame(out.fine, out.coarse)
write.csv(out.all, 
          paste0("./data/synth_data_summarized_5km_radius_", today, ".csv"), 
          row.names = FALSE)


################################################################################
################### PLOTTING TIME 5 km #########################################
################################################################################


# # read in the data (if you don't want to re-run above)
# out.all <- read.csv("./data/synth_data_summarized_5km_radius_20181102.csv", 
#                     stringsAsFactors = FALSE)

# remove na columns
trash <- which(colSums(is.na(out.all)) == nrow(out.all))
out.all <- out.all[, -trash]

# now melt the data into a ggplot-friendly 'tidy' matrix
out.gg <- melt(out.all, id = c("surf.name", "surf.type", "fine.coarse", "poly"))
out.gg$fine.coarse <- factor(out.gg$fine.coarse, 
                             levels = c("fine", "coarse"))
out.gg$surf.type <- factor(out.gg$surf.type, 
                           levels = unique(out.all$surf.type))

# out.gg$variable <- factor(out.gg$variable, 
#                           levels = c(names(in.fine), names(in.coarse)))

# trying to plot with ggplot
png(filename = paste0("./figures/geodiv_metrics_5km_radius_", today, ".png"),
    width = 7, height = 15, units = "in", res = 300)

  ggplot(data = out.gg, mapping = aes(x = surf.name, 
                                      y = value, 
                                      col = fine.coarse,
                                      fill = fine.coarse)) +
    geom_dotplot(binaxis = "y", 
                 stackdir = "center", 
                 dotsize = 0.8, 
                 aes(surf.type)) +
    facet_grid(variable ~ fine.coarse, 
               scales = "free") +
    labs(x = "Geo Surface", y = "", 
         title = "Surface Geodiversity Metrics (5 km radius)") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 7))
  
dev.off()

##### now same for the 20 km radius ############################
variables <- c("st.deviation",
               "kurtosis", 
               "skewness",
               "range")
out.fine <- as.data.frame(matrix(NA, 
                                 nrow = 12 * 6, 
                                 ncol = length(variables) + 4))
names(out.fine) <- c("surf.name", "surf.type", "fine.coarse", "poly", variables)

out.coarse <- out.fine

date()
for (y in 1:12) {
  in.circ <- extract(in.fine, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:6) {  
    i <- (y - 1) * 6 + x
    out.fine$surf.name[i] <- unique.ids.fine[x]
    out.fine$surf.type[i] <- names(in.fine)[x]
    out.fine$fine.coarse[i] <- "fine"
    out.fine$poly[i] <- y
    
    out.fine$st.deviation[i] <- sd(in.circ[, x])
    out.fine$kurtosis[i] <- kurtosis(in.circ[, x])
    out.fine$skewness[i] <- skewness(in.circ[, x])
    out.fine$range[i] <- max(in.circ[, x]) - min(in.circ[, x])
  }
  print(paste("done with polygon #", y, date()))
}

### NOW COARSE ###

date()
for (y in 1:12) {
  in.circ <- extract(in.coarse, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:6) {  
    i <- (y - 1) * 6 + x
    out.coarse$surf.name[i] <- unique.ids.coarse[x]
    out.coarse$surf.type[i] <- names(in.coarse)[x]
    out.coarse$fine.coarse[i] <- "coarse"
    out.coarse$poly[i] <- y
    
    out.coarse$st.deviation[i] <- sd(in.circ[, x])
    out.coarse$kurtosis[i] <- kurtosis(in.circ[, x])
    out.coarse$skewness[i] <- skewness(in.circ[, x])
    out.coarse$range[i] <- max(in.circ[, x]) - min(in.circ[, x])
  }
  print(paste("done with polygon #", y, date()))
}

################################################################################
################### PLOTTING TIME ##############################################
################################################################################

# combine the two dataframes
out.all <- rbind.data.frame(out.fine, out.coarse)
write.csv(out.all, 
          paste0("synth_data_summarized_20km_radius_", today, ".csv"), 
          row.names = FALSE)

# now melt the data into a ggplot-friendly 'tidy' matrix
out.gg <- melt(out.all, id = c("surf.name", "surf.type", "fine.coarse", "poly"))
out.gg$fine.coarse <- factor(out.gg$fine.coarse, 
                             levels = c("fine", "coarse"))
out.gg$surf.type <- factor(out.gg$surf.type, 
                           levels = geo.names.short)

# out.gg$variable <- factor(out.gg$variable, 
#                           levels = c(names(in.fine), names(in.coarse)))

# trying to plot with ggplot
ggplot(data = out.gg, mapping = aes(x = surf.name, 
                                    y = value, 
                                    col = fine.coarse,
                                    fill = fine.coarse)) +
  geom_dotplot(binaxis = "y", 
               stackdir = "center", 
               dotsize = 1.2, 
               aes(surf.type)) +
  facet_grid(variable ~ fine.coarse, 
             scales = "free") +
  labs(x = "Geo Surface", y = "", 
       title = "Surface Geodiversity Metrics (20 km radius)") +
  theme_bw() +
  theme(legend.position = "none")



