### extracting data and summarizing for synthetic geodiversity ####
library(raster)
library(rgdal)
library(plotrix)
library(e1071)
library(ggplot2)
library(reshape)

# setwd and today's date for file writing
setwd("C:/Users/kdahlin/Dropbox/NASA_biodiversity/synthetic_geodiversity")
today <- format(Sys.Date(), "%Y%m%d")

# set data directory
data.dir <- "Y:/data/synthetic_geodiversity/"

# set layer names
geo.names <- c("Stripe", "Check", "Slope", "Hill", "G.Patchy", "G.Smooth",
               "AZ.therm", "KS.therm")
geo.names.short <- c("St", "Ch","Sl", "Hi", "GP", "GS", "AZ.TH", "KS.TH")

# read in synthetic data from step A
in.fine <- stack(paste0(data.dir, 
                        "synthetic_surfaces_30m_20181031.tif"))

# realized the gaussian patchy layer has some negative values, which could mess
# up the math later, so setting all those neg values to zero
in.fine[[5]][in.fine[[5]] < 0] <- 0

# read in real thermal data and remove crs and extent so it can be stacked with
# synthetic data
in.AZ.landsat <- raster(paste0(data.dir, 
                                "arizona/arizona_landsat_square.tif"))
crs(in.AZ.landsat) <- NA
extent(in.AZ.landsat) <- extent(in.fine)

in.KS.landsat <- raster(paste0(data.dir, 
                               "kansas/kansas_landsat_square.tif"))
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
                          "synthetic_surfaces_1000m_20181031.tif"))

# modis
in.AZ.modis <- raster(paste0(data.dir, 
                               "arizona/arizona_modis_square.tif"))
crs(in.AZ.modis) <- NA

# aaannnnnddd these don't line up perfectly, so cropping. :\ would be better to
# go back and fix code in step B.
extent(in.AZ.modis) <- extent(c(0,145,0,145))
in.AZ.modis <- crop(in.AZ.modis, extent(in.coarse))

in.KS.modis <- raster(paste0(data.dir, 
                             "kansas/kansas_modis_square.tif"))
crs(in.KS.modis) <- NA

extent(in.KS.modis) <- extent(c(0,145,0,145))
in.KS.modis <- crop(in.KS.modis, extent(in.coarse))

# now all together
all.coarse <- stack(in.coarse, in.AZ.modis, in.KS.modis)

names(all.coarse) <- geo.names
unique.ids.coarse <- paste0(geo.names, "1k")
crs(all.coarse) <- "+proj=utm +zone=12 +datum=WGS84 +units=km +no_defs"

writeRaster(all.coarse, paste0(data.dir, 
                             "all_coarse_res_layers_",
                             today,
                             ".tif"),
            overwrite = TRUE)


# read in centroids csv
in.centers <- read.csv("./data/sampling_centroids_x12_20181031.csv",
                       stringsAsFactors = FALSE)

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


### OK! now let's extract some data and write it! ## 
##### 5 km radius first ############################
variables <- c("st.deviation",
               "kurtosis", 
               "skewness",
               "range",
               "SBI.mean",
               "SBI.sd",
               "prof.convexity.mean",
               "prof.convexity.sd",
               "max.curvature.mean",
               "max.curvature.sd",
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
               "text.dir.mean",
               "text.dir.sd",
               "RWI.mean",
               "RWI.sd",
               "variogram.diff"
               )
out.fine <- as.data.frame(matrix(NA, 
                                 nrow = 12*8, 
                                 ncol = length(variables)+4))
names(out.fine) <- c("surf.name", "surf.type", "fine.coarse", "poly", variables)

out.coarse <- out.fine

# first let's do the simple metrics (summary of original data) 
# this is super slow -> methinks its the extraction step
date()
for (y in 1:12) {
  in.circ <- extract(all.fine, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
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
  in.circ <- extract(all.fine.slope, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$slope.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.fine$slope.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.fine, paste0("./data/intermed_to_delete_", today, ".csv"))

############################# SURFACE AREA #####################################

all.fine.area <- (0.03*0.03)/cos(all.fine.slope)
names(all.fine.area) <- paste0("area.", geo.names.short)

## SURF AREA ##
for (y in 1:12) {
  in.circ <- extract(all.fine.area, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$surf.area.sum[i] <- sum(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.fine.tpi, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$TPI.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.fine$TPI.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.fine.tri, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$TRI.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.fine$TRI.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.fine.sine.asp, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$aspect.sine.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.fine$aspect.sine.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.fine.rough, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.fine$roughness.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.fine$roughness.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$surf.name[i] <- unique.ids.coarse[x]
    out.coarse$surf.type[i] <- names(all.coarse)[x]
    out.coarse$fine.coarse[i] <- "coarse"
    out.coarse$poly[i] <- y
    
    out.coarse$st.deviation[i] <- sd(in.circ[,x], na.rm = TRUE)
    out.coarse$kurtosis[i] <- kurtosis(in.circ[,x], na.rm = TRUE)
    out.coarse$skewness[i] <- skewness(in.circ[,x], na.rm = TRUE)
    out.coarse$range[i] <- max(in.circ[,x], na.rm = TRUE) - 
      min(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse.slope, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$slope.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.coarse$slope.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

############################# SURFACE AREA #####################################

all.coarse.area <- (1*1)/cos(all.coarse.slope)
names(all.coarse.area) <- paste0("area.", geo.names.short)

## SURF AREA ##
for (y in 1:12) {
  in.circ <- extract(all.coarse.area, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$surf.area.sum[i] <- sum(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse.tpi, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$TPI.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.coarse$TPI.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse.tri, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$TRI.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.coarse$TRI.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse.sine.asp, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$aspect.sine.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.coarse$aspect.sine.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
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
  in.circ <- extract(all.coarse.rough, circ.20.SPDF[y,])
  in.circ <- as.data.frame(in.circ[[1]])
  for (x in 1:8) {  
    i <- (y-1)*8 + x
    out.coarse$roughness.mean[i] <- mean(in.circ[,x], na.rm = TRUE)
    out.coarse$roughness.sd[i] <- sd(in.circ[,x], na.rm = TRUE)
  }
  print(paste("done with polygon #", y, date()))
}

write.csv(out.coarse, paste0("./data/intermed_to_delete_", today, ".csv"))

write.csv(out.coarse, paste0("./data/coarse_surf_metrics_ALL_", today, ".csv"))

# combine the two dataframes
out.all <- rbind.data.frame(out.fine, out.coarse)
write.csv(out.all, 
          paste0("./data/synth_data_summarized_20km_radius_", today, ".csv"), 
          row.names = FALSE)


################################################################################
################### PLOTTING TIME 20 km ########################################
################################################################################


# # read in the data (if you don't want to re-run above)
# out.all <- read.csv("./data/synth_data_summarized_5km_radius_20181102.csv", 
#                     stringsAsFactors = FALSE)

# remove na columns
trash <- which(colSums(is.na(out.all)) == nrow(out.all))
out.all <- out.all[,-trash]

# now melt the data into a ggplot-friendly 'tidy' matrix
out.gg <- melt(out.all, id = c("surf.name", "surf.type", "fine.coarse", "poly"))
out.gg$fine.coarse <- factor(out.gg$fine.coarse, 
                             levels = c("fine", "coarse"))
out.gg$surf.type <- factor(out.gg$surf.type, 
                           levels = unique(out.all$surf.type))

# out.gg$variable <- factor(out.gg$variable, 
#                           levels = c(names(in.fine), names(in.coarse)))

# trying to plot with ggplot
png(filename = paste0("./figures/geodiv_metrics_20km_radius_", today, ".png"),
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
         title = "Surface Geodiversity Metrics (20 km radius)") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = 7))
  
dev.off()

