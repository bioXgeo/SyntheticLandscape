# Test script for geodiv R package ----------------------------------------

# This tests all functions associated with the geodiv package and writes
# errors to a log file.

# Instructions (all steps are coded out below):
# 1: Load geodiv: clone git repository onto local machine
# 2: Run geodiv using both the 'normforest' and 'orforest' rasters included
#    in the package.
# 3: Run geodiv with a test raster of your choosing (make sure that it is
#    fairly small).
# 4: Add resulting log files to the geodiv Google Drive directory.
# 5: Send any ideas that people might want to do extra.

# load geodiv -------------------------------------------------------------

sink('/home/annie/Desktop/geodiv_logfile.txt', split = TRUE)

library(devtools)
library(roxygen2)

devtools::load_all()

# function to run all geodiv functions ------------------------------------

run_all <- function(r){

  # run all functions
  cat('Sa: ', '\n')
  sa <- sa(r)
  print(sa)
  cat('\n')

  cat('Sq: ', '\n')
  sq <- sq(r)
  print(sq)
  cat('\n')

  cat('S10z: ', '\n')
  s10z <-s10z(r)
  print(s10z)
  cat('\n')

  cat('Sdq: ', '\n')
  sdq <- sdq(r)
  print(sdq)
  cat('\n')

  cat('Sdq6: ', '\n')
  sdq6 <- sdq6(r)
  print(sdq6)
  cat('\n')

  cat('Sdr: ', '\n')
  sdr <- sdr(r)
  print(sdr)
  cat('\n')

  cat('Sbi: ', '\n')
  sbi <- sbi(r)
  print(sbi)
  cat('\n')

  cat('Sci: ', '\n')
  sci <- sci(r)
  print(sci)
  cat('\n')

  cat('Ssk (Adjacent): ', '\n')
  ssk_adj <- ssk(r, adj = TRUE)
  print(ssk_adj)
  cat('\n')

  cat('Ssk: ', '\n')
  ssk_nonadj <- ssk(r, adj = FALSE)
  print(ssk_nonadj)
  cat('\n')

  cat('Sku (excess): ', '\n')
  sku_excess <- sku(r, excess = TRUE)
  print(sku_excess)
  cat('\n')

  cat('Sku: ', '\n')
  sku_nonexcess <- sku(r, excess = FALSE)
  print(sku_nonexcess)
  cat('\n')

  cat('Sds: ', '\n')
  sds <- sds(r)
  print(sds)
  cat('\n')

  cat('Sfd: ', '\n')
  sfd <- sfd(r)
  print(sfd)
  cat('\n')

  cat('Srw function - no plotting: ', '\n')
  srwvals <- srw(r, plot = FALSE) # this takes a while
  print(srwvals)
  cat('\n')

  cat('Srw function - plotting: ', '\n')
  srwvals <- srw(r, plot = TRUE)
  print(srwvals)
  cat('\n')

  cat('Srw: ', '\n')
  srw <- srwvals[[1]]
  print(srw)
  cat('\n')

  cat('Srwi: ', '\n')
  srwi <- srwvals[[2]]
  print(srwi)
  cat('\n')

  cat('Shw: ', '\n')
  shw <- srwvals[[3]]
  print(shw)
  cat('\n')

  cat('Std function - no plotting: ', '\n')
  stdvals <- std(r, plot = FALSE)
  print(stdvals)
  cat('\n')

  cat('Std function - plotting: ', '\n')
  stdvals <- std(r, plot = TRUE)
  print(stdvals)
  cat('\n')

  cat('Std: ', '\n')
  std <- stdvals[[1]]
  print(std)
  cat('\n')

  cat('Stdi: ', '\n')
  stdi <- stdvals[[2]]
  print(stdi)
  cat('\n')

  cat('Svi: ', '\n')
  svi <- svi(r)
  print(svi)
  cat('\n')

  cat('Str function: ', '\n')
  strvals <- str(r, threshold = c(0.2, 1 / exp(1)))
  print(strvals)
  cat('\n')

  cat('Str (20%): ', '\n')
  str20 <- strvals[[1]]
  print(str20)
  cat('\n')

  cat('Str (37%): ', '\n')
  str37 <- strvals[[2]]
  print(str37)
  cat('\n')

  cat('Ssc: ', '\n')
  ssc <- ssc(r)
  print(ssc)
  cat('\n')

  cat('Sv: ', '\n')
  sv <- sv(r)
  print(sv)
  cat('\n')

  cat('Sp: ', '\n')
  sp <- sp(r)
  print(sp)
  cat('\n')

  cat('Sk: ', '\n')
  sk <- sk(r)
  print(sk)
  cat('\n')

  cat('Smean: ', '\n')
  smean <- smean(r)
  print(smean)
  cat('\n')

  cat('Spk: ', '\n')
  spk <- spk(r)
  print(spk)
  cat('\n')

  cat('Svk: ', '\n')
  svk <- svk(r)
  print(svk)
  cat('\n')

  cat('Scl function - plotting: ', '\n')
  sclvals <- scl(r, threshold = c(0.2, 1 / exp(1)), plot = TRUE)
  print(sclvals)
  cat('\n')

  cat('Scl function - no plotting: ', '\n')
  sclvals <- scl(r, threshold = c(0.2, 1 / exp(1)), plot = FALSE)
  print(sclvals)
  cat('\n')

  cat('Scl (20%): ', '\n')
  scl20 <- sclvals[[1]] # never gets to 0.2
  print(scl20)
  cat('\n')

  cat('Scl (37%): ', '\n')
  scl37 <- sclvals[[2]]
  print(scl37)
  cat('\n')

  cat('Sdc 0-5%: ', '\n')
  sdc0_5 <- sdc(r, 0, 0.05)
  print(sdc0_5)
  cat('\n')

  cat('Sdc 50-55%: ', '\n')
  sdc50_55 <- sdc(r, 0.50, 0.55)
  print(sdc50_55)
  cat('\n')

  cat('Sdc 80-85%: ', '\n')
  sdc80_85 <- sdc(r, 0.80, 0.85)
  print(sdc80_85)
  cat('\n')

  return(cat('Finished with raster.', '\n'))
}

# run geodiv using orforest data ------------------------------------------

data(orforest)
cat('orforest raster: ', '\n')
print(orforest)

# run all functions
run_all(orforest)

# run everything for normforest -------------------------------------------

data(normforest)
cat('normforest raster: ', '\n')
print(normforest)

# run all functions
run_all(normforest)

# test for another raster -------------------------------------------------

rast <- raster('filename') # add file here
cat('other raster: ', '\n')
print(rast)

# remove any trend
r <- remove_plane(rast)

# run all functions
run_all(r)

sink()