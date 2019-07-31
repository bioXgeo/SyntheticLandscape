### reading in thermal landsat and MODIS data, formatting, and clipping ####

library(raster)
library(rgdal)

setwd("~/Desktop/synthetic_geodiversity/data/")

### Arizona ###
# set location for AZ data (HPCC nasabio directory)
az.loc <- "~/Desktop/synthetic_geodiversity/data/arizona/"

# read in landsat band 11
az.in.landsat <- raster(paste0(az.loc, 
                               "LC08_L1TP_037036_20180702_20180717_01_T1/",
                               "LC08_L1TP_037036_20180702_20180717_01_T1_B11.TIF"))
# get the projection
az.proj <- as.character(crs(az.in.landsat))

# making a square polygon (note, I figured out the coordinates for these
# by looking at them in ArcGIS) - starting in LL, going counter clockwise
az.sq.x <- c(324800, 324800 + 144000, 324800 + 144000, 324800)
az.sq.y <- c(3752000, 3752000, 3752000 + 144000, 3752000 + 144000)
az.sq.mat <- cbind(az.sq.x, az.sq.y)
az.sq.poly <- Polygon(az.sq.mat)
az.sq.poly2 <- Polygons(list(az.sq.poly), 1)
az.sq <- SpatialPolygons(list(az.sq.poly2), proj4string = CRS(az.proj))
data.forSPDF <- as.data.frame(matrix(1, nrow = 1, ncol = 1))
names(data.forSPDF) <- "plot.id"
az.sq.out <- SpatialPolygonsDataFrame(az.sq, data = data.forSPDF)

# write this as a shapefile so you can extract same data via AppEEARS for MODIS
writeOGR(az.sq.out, dsn = ".", layer = "arizona_square_20190521", driver = "ESRI Shapefile")

# now back to Landsat data
# set zeros to NAs
az.in.landsat[az.in.landsat == 0] <- NA

# take a look!
plot(az.in.landsat)
plot(az.sq.out, add = TRUE)

# read in QA band and mask out bad pix
az.in.QA <- raster(paste0(az.loc, 
                          "LC08_L1TP_037036_20180702_20180717_01_T1/",
                          "LC08_L1TP_037036_20180702_20180717_01_T1_BQA.TIF"))

# set 1s (fill values) to NA
az.in.QA[az.in.QA == 1] <- NA

# now set all the keep values to 1
# values from https://landsat.usgs.gov/collectionqualityband
az.in.QA[az.in.QA == 2720] <- 1
az.in.QA[az.in.QA == 2724] <- 1
az.in.QA[az.in.QA == 2728] <- 1
az.in.QA[az.in.QA == 2732] <- 1

# now set the remaining values to NAs
az.in.QA[az.in.QA > 1] <- NA

# apply this mask (just 1s and NAs) to the landsat radiance layer
az.landsat.masked <- az.in.landsat * az.in.QA

# convert to TOA radiance 
# (from https://landsat.usgs.gov/using-usgs-landsat-8-product)
# numbers from the .txt file for this scene (in same folder)

az.toa.rad.landsat <- az.landsat.masked * 3.3420E-04 + 0.10000

# convert to brightness temperature 
# (from https://landsat.usgs.gov/using-usgs-landsat-8-product)
# numbers from the .txt file for this scene (in same folder)

az.BT.landsat <- 1201.1442 / log(480.8883 / (az.toa.rad.landsat + 1))

az.landsat.sq <- crop(az.BT.landsat, az.sq.out)

# take a look
plot(az.landsat.sq)

# now scale this so it ranges from 0 to 100, like the synthetic landscapes
min.az <- cellStats(az.landsat.sq, min)
max.az <- cellStats(az.landsat.sq, max)
new.max.az <- max.az - min.az

az.landsat.sq.sc <- ((az.landsat.sq - min.az) / new.max.az) * 100

# write this Raster!
writeRaster(az.landsat.sq.sc, paste0(az.loc, "arizona_landsat_square_20190521.tif"))

cellStats(az.landsat.sq.sc, mean)
cellStats(az.landsat.sq.sc, sd)

#### NOW MODIS 
#### NOTE THIS IS WHAT YOU WANT FOR MODIS: 
## https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/modtbga_v006
## ! still need to deal with QA for MODIS. 
az.in.modis <- raster(paste0(az.loc, 
                             "MODIS/",
                             "MODTBGA.006_BAND32_1_doy2018183_aid0001.tif"))

# reproject MODIS data to match Landsat (with 1000 m resolution)
az.modis.out <- projectRaster(az.in.modis, res = 1000, crs = crs(az.landsat.sq.sc),
                              method = "bilinear", alignOnly = FALSE)

min.az <- cellStats(az.modis.out, min)
max.az <- cellStats(az.modis.out, max)
new.max.az <- max.az - min.az

az.modis.sq.sc <- ((az.modis.out - min.az)/new.max.az)*100

# write this Raster!
writeRaster(az.modis.sq.sc, paste0(az.loc, "arizona_modis_square_20190521.tif"))

cellStats(az.modis.sq.sc, mean)
cellStats(az.modis.sq.sc, sd)


### Kansas ###
# set location for KS data
ks.loc <- "~/Desktop/synthetic_geodiversity/data/kansas/"

# read in landsat band 11
ks.in.landsat <- raster(paste0(ks.loc, 
                               "LC08_L1TP_028033_20180601_20180614_01_T1/",
                               "LC08_L1TP_028033_20180601_20180614_01_T1_B11.TIF"))
# get the projection
ks.proj <- as.character(crs(ks.in.landsat))

# making a square polygon (note, I figured out the coordinates for these
# by looking at them in ArcGIS) - starting in LL, going counter clockwise
ks.sq.x <- c(608000, 608000 + 144000, 608000 + 144000, 608000)
ks.sq.y <- c(4240000, 4240000, 4240000 + 144000, 4240000 + 144000)
ks.sq.mat <- cbind(ks.sq.x, ks.sq.y)
ks.sq.poly <- Polygon(ks.sq.mat)
ks.sq.poly2 <- Polygons(list(ks.sq.poly), 1)
ks.sq <- SpatialPolygons(list(ks.sq.poly2), proj4string = CRS(ks.proj))
data.forSPDF <- as.data.frame(matrix(1, nrow = 1, ncol = 1))
names(data.forSPDF) <- "plot.id"
ks.sq.out <- SpatialPolygonsDataFrame(ks.sq, data = data.forSPDF)

# write this as a shapefile so you can extract same data via AppEEARS for MODIS
writeOGR(ks.sq.out, dsn = ".", layer = "kansas_square_20190521", driver = "ESRI Shapefile")

# now back to Landsat data
# set zeros to NAs
ks.in.landsat[ks.in.landsat == 0] <- NA

# take a look!
plot(ks.in.landsat)
plot(ks.sq.out, add = TRUE)

# read in QA band and mask out bad pix
ks.in.QA <- raster(paste0(ks.loc, 
                          "LC08_L1TP_028033_20180601_20180614_01_T1/",
                          "LC08_L1TP_028033_20180601_20180614_01_T1_BQA.TIF"))

# set 1s (fill values) to NA
ks.in.QA[ks.in.QA == 1] <- NA

# now set all the keep values to 1
# values from https://landsat.usgs.gov/collectionqualityband
ks.in.QA[ks.in.QA == 2720] <- 1
ks.in.QA[ks.in.QA == 2724] <- 1
ks.in.QA[ks.in.QA == 2728] <- 1
ks.in.QA[ks.in.QA == 2732] <- 1

# now set the remaining values to NAs
ks.in.QA[ks.in.QA > 1] <- NA

# apply this mask (just 1s and NAs) to the landsat radiance layer
ks.landsat.masked <- ks.in.landsat * ks.in.QA

# convert to TOA radiance 
# (from https://landsat.usgs.gov/using-usgs-landsat-8-product)
# numbers from the .txt file for this scene (in same folder)

ks.toa.rad.landsat <- ks.landsat.masked * 3.3420E-04 + 0.10000

# convert to brightness temperature 
# (from https://landsat.usgs.gov/using-usgs-landsat-8-product)
# numbers from the .txt file for this scene (in same folder)

ks.BT.landsat <- 1201.1442 / log(480.8883 / (ks.toa.rad.landsat + 1))

ks.landsat.sq <- crop(ks.BT.landsat, ks.sq.out)

# take a look
plot(ks.landsat.sq)

# now scale this so it ranges from 0 to 100, like the synthetic landscapes
min.ks <- cellStats(ks.landsat.sq, min)
max.ks <- cellStats(ks.landsat.sq, max)
new.max.ks <- max.ks - min.ks

ks.landsat.sq.sc <- ((ks.landsat.sq - min.ks) / new.max.ks) * 100

# write this Raster!
writeRaster(ks.landsat.sq.sc, paste0(ks.loc, "kansas_landsat_square_20190521.tif"))

cellStats(ks.landsat.sq.sc, mean)
cellStats(ks.landsat.sq.sc, sd)

#### NOW MODIS 
## NOTE THIS IS WHAT YOU WANT FOR MODIS: 
## https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/modtbga_v006
## ! still need to deal with QA for MODIS. 
ks.in.modis <- raster(paste0(ks.loc, 
                             "MODIS/",
                             "MODTBGA.006_BAND32_1_doy2018152_aid0001.tif"))

# reproject MODIS data to match Landsat (with 1000 m resolution)
ks.modis.out <- projectRaster(ks.in.modis, res = 1000, crs = crs(ks.landsat.sq.sc),
                              method = "bilinear", alignOnly = FALSE)

min.ks <- cellStats(ks.modis.out, min)
max.ks <- cellStats(ks.modis.out, max)
new.max.ks <- max.ks - min.ks

ks.modis.sq.sc <- ((ks.modis.out - min.ks) / new.max.ks) * 100

# write this Raster!
writeRaster(ks.modis.sq.sc, paste0(ks.loc, "kansas_modis_square_20190521.tif"))

cellStats(ks.modis.sq.sc, mean)
cellStats(ks.modis.sq.sc, sd)










