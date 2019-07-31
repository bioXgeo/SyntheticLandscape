### assessing clumpiness ###

library(raster)
library(automap)

setwd("C:/Users/kdahlin/Dropbox/NASA_biodiversity/clumpiness")

in.half <- as.matrix(read.csv("fake_data_half.csv"))
in.half <- raster(in.half)

in.quarters <- as.matrix(read.csv("fake_data_quarters.csv"))
in.quarters <- raster(in.quarters)

in.grid <- as.matrix(read.csv("fake_data_grid.csv"))
in.grid <- raster(in.grid)

in.patchy <- as.matrix(read.csv("fake_data_patchy.csv"))
in.patchy <- raster(in.patchy)

in.all <- stack(in.half, in.quarters, in.grid, in.patchy)
extent(in.all) <- c(0,1,0,1)
crs(in.all) <- "+proj=longlat"
names(in.all) <- c("in.half", "in.quarters", "in.grid", "in.patchy")

patchy.points <- rasterToPoints(in.all$in.patchy, spatial = TRUE)
patchy.var <- autofitVariogram(in.patchy ~ 1, patchy.points, 
                               model = c("Sph", "Exp", "Gau"))

writeRaster(in.all, "C:/Users/kdahlin/Dropbox/NASA_biodiversity/synthetic_geodiversity/data/circles.tif")

x11()
par(mfrow = c(1,4))
plot(in.all)

out.table <- as.data.frame(matrix(NA, nrow = 4, ncol = 10))
out.table[,1] <- names(in.all)
names(out.table) <- c("data", "mean", "sd", "skew", "sill", "nugget", "range",
                      "nd.sill.nugget.ratio", "TPI.mean", "TRI.mean")

for (i in 1:4) {
  in.data <- in.all[[i]]
  out.mean <- cellStats(in.data, stat = "mean")
  out.table$mean[i] <- out.mean
  out.sd <- cellStats(in.data, stat = "sd")
  out.table$sd[i] <- out.sd
  out.skew <- cellStats(in.data, stat = "skew")
  out.table$skew[i] <- out.skew

  data.points <- rasterToPoints(in.data, spatial = TRUE)
  data.info <- rasterToPoints(in.data, spatial = FALSE)
  data.var <- autofitVariogram(data.info[,3] ~ 1, data.points, 
                                 model = c("Sph", "Exp", "Gau"))
  out.sill <- data.var$var_model$psill[2] + data.var$var_model$psill[1]
  out.table$sill[i] <- out.sill
  out.nugget <- data.var$var_model$psill[1]
  out.table$nugget[i] <- out.nugget
  out.range <- data.var$var_model$range[2]
  out.table$range[i] <- out.range
  out.nd.sill.nug <- (out.sill - out.nugget)/(out.sill + out.nugget)
  out.table$nd.sill.nugget.ratio[i] <- out.nd.sill.nug
  
  out.TPI <- terrain(in.data, opt = 'TPI')
  TPI.mean <- cellStats(out.TPI, stat = "mean")
  out.table$TPI.mean[i] <- TPI.mean
  
  out.TRI <- terrain(in.data, opt = 'TRI')
  TRI.mean <- cellStats(out.TRI, stat = "mean")
  out.table$TRI.mean[i] <- TRI.mean
  
}

x11()
par(mfrow = c(1,9))
barplot(as.matrix(out.table[,2]), beside = TRUE, legend.text = as.character(out.table$data),
        main = names(out.table)[2])
for(i in 3:10) {
  barplot(as.matrix(out.table[,i]), beside = TRUE, main = names(out.table)[i])
}
  

### summary stats & clumpiness estimate ideas ####
# mean
# sd
# skew
# variogram sill
# variogram range
# variogram sill:nugget ratio (larger ratio == more spatial structure)
# mean TPI
# mean TRI
