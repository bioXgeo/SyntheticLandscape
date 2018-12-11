# function to calculate surface area, surface area ratio
# from https://www.ntmdt-si.ru/data/media/files/manuals/image_analisys_p9_nov12.e.pdf
# diagram for actual equations used: file:///home/annie/Downloads/9783642364570-c2%20(1).pdf, page 17/30
flatsa <- function(rast) {
  # calculates surface area of flat plane covering sample area (N-1, M-1)
  
  # convert to equal area for resolution
  # aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' 
  # rast <- projectRaster(rast, crs = aea_crs) # note that this changes the # pixels (b/c rotation)
  
  # get dimensions
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  deltax <- res(rast)[1]
  deltay <- res(rast)[2]
  
  # normalize deltax, y
  divide <- mean(deltax, deltay)
  deltax <- deltax / divide
  deltay <- deltay / divide
  
  sa <- ((N - 1) * (M - 1)) * deltax * deltay
  
  return(sa)
}

sdr <- function(rast) {
  # calculates surface area of sample area (N-1, M-1)
  # because most satellite data has units that are not of the same scale as x,y,
  # this is not super meaningful
  # i decided to call x, y = 1, and then scaled z to between 0 and 1 to best match units
  
  # get area of flat plane
  flat_area <- flatsa(rast)
  
  # get dimensions
  N <- dim(rast)[1] # rows
  M <- dim(rast)[2] # cols
  
  # z values, coordinates, and resolution (change in x, y)
  z <- getValues(rast)
  x <- coordinates(rast)[, 1]
  y <- coordinates(rast)[, 2]
  deltax <- res(rast)[1]
  deltay <- res(rast)[2]
  
  # create matrix of centers to get surrounding from
  zmat <- matrix(((z - min(z)) / (max(z) - min(z))), nrow = N, ncol = M, byrow = TRUE)
  
  z <- zshift(zmat, xdist = 0, ydist = 0, xrm = 1, yrm = 1)
  z_ypl <- zshift(zmat, xdist = 0, ydist = 1, xrm = 1)
  z_xpl <- zshift(zmat, xdist = 1, ydist = 0, yrm = 1)
  z_xplypl <- zshift(zmat, xdist = 1, ydist = 1)
  
  # normalize deltas
  divide <- mean(deltax, deltay)
  deltax <- deltax / divide
  deltay <- deltay / divide
  
  # calculate area - problem with ndvi not being equal to x, y units
  # NEED TO SCALE TO MATCH X, Y UNITS -- this should be checked!
  akl1 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2))))
    
  akl2 <- ((1 / 2) *
             (sqrt((deltay ^ 2) + ((z_ypl - z) ^ 2)) *
                sqrt((deltax ^ 2) + ((z_xplypl - z_ypl) ^ 2)))) +
    ((1 / 2) *
       (sqrt((deltay ^ 2) + ((z_xplypl - z_xpl) ^ 2)) *
          sqrt((deltax ^ 2) + ((z_xpl - z) ^ 2))))
  
  akl <- (1 / 2) * (akl1 + akl2)
  
  # calculate area ratio
  adr <- ((sum(akl, na.rm = TRUE) - flat_area) / flat_area) * 100
  
  return(adr)
}