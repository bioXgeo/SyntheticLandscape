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

# run geodiv using orforest data ------------------------------------------

data(orforest)
cat('orforest raster: ', '\n')
print(orforest)

# run all functions
sa <- sa(orforest)
sq <- sq(orforest)
s10z <-s10z(orforest)
sdq <- sdq(orforest)
sdq6 <- sdq6(orforest)
sdr <- sdr(orforest)
sbi <- sbi(orforest)
sci <- sci(orforest)
ssk_adj <- ssk(orforest, adj = TRUE)
ssk_nonadj <- ssk(orforest, adj = FALSE)
sku_excess <- sku(orforest, excess = TRUE)
sku_nonexcess <- sku(orforest, excess = FALSE)
sds <- sds(orforest)
sfd <- sfd(orforest)
srwvals <- srw(orforest, plot = FALSE) # this takes a while
srwvals <- srw(orforest, plot = TRUE)
srw <- srwvals[[1]]
srwi <- srwvals[[2]]
shw <- srwvals[[3]]
stdvals <- std(orforest, plot = FALSE)
stdvals <- std(orforest, plot = TRUE)
std <- stdvals[[1]]
stdi <- stdvals[[2]]
svi <- svi(orforest)
strvals <- str(orforest, threshold = c(0.2, 1/exp(1)))
str20 <- strvals[[1]]
str37 <- strvals[[2]]
ssc <- ssc(orforest)
sv <- sv(orforest)
sp <- sp(orforest)
sk <- sk(orforest)
smean <- smean(orforest)
spk <- spk(orforest)
svk <- svk(orforest)
sclvals <- scl(orforest, threshold = c(0.2, 1/exp(1)), plot = TRUE)
sclvals <- scl(orforest, threshold = c(0.2, 1/exp(1)), plot = FALSE)
scl20 <- sclvals[[1]]
scl37 <- sclvals[[3]]
sdc0_5
sdc50_55
sdc80_85
