
bearing_area <- function(x) {
  z <- getValues(x)

  # basic values
  N <- length(errors)
  s <- sd(errors)
  z <- errors
  zbar <- mean(errors)

  f <- ecdf(1 - z)

  yval <- (1 - environment(f)$x)
  xval <- environment(mod)$y

  return(f)
}

plot_ba_curve <- function(f, divisions = FALSE) {
  xval <- environment(f)$y
  yval <- (1 - environment(f)$x)

  plot(yval ~ xval)

  if (divisions == TRUE) {
    line_fit <- find_flat(f, perc = 0.4)

    abline(line_fit[[1]]$coefficients[[1]], line_fit[[1]]$coefficients[[2]], col = 'blue')
    abline(h = line_fit[[2]], col = 'red')
    abline(h = line_fit[[3]], col = 'red')
    abline(v = line_fit[[4]], col = 'green')
    abline(v = line_fit[[5]], col = 'green')
  }
}

find_flat <- function(x, perc = 0.4) {
  f <- bearing_area(x)

  # find 40% of curve with least decline
  # use symmetric difference quotient to estimate the derivative at evenly spaced points
  # then find 40% consecutive section with lowest mean slope
  even_x <- seq(0, 1, length.out = 100000)
  even_y <- (1 - quantile(f, probs = even_x))
  forty_length <- perc * length(even_x)
  h <- 0.001
  slopes <- slopecalc(even_x, h) # calculate slope at every point
  means <- slopemeans(slopes) # calculate averages for each 40% segment

  # x value of start of 40% section with smallest decline
  slope_min <- means[means$slope == min(means$slope),]

  # calculate least-squares line for 40% of curve with smallest decline (lowest slope)
  lm_data <- data.frame(x = newx[newx >= slope_min$xstart & newx <= slope_min$xend],
                        y = newy[newx >= slope_min$xstart & newx <= slope_min$xend])
  ls_line <- lm(y ~ x, data = lm_data)

  # get value of ls line between 0 and 1
  pred_data <- remove_rownames(data.frame(x = even_x, y = even_y))
  pred_data$y <- predict(ls_line, newdata = pred_data)

  # what is the ls line y-value at x = 0, x = 1?
  ls_int_high <- pred_data$y[pred_data$x == 0]
  ls_int_low <- pred_data$y[pred_data$x == 1]

  # Smr1/Smr2 = x values that correspond to cdf y values at ls_int_high/low
  Smr1 <- f(1 - ls_int_high)
  Smr2 <- f(1 - ls_int_low)

  return(list(ls_line, pred_data, ls_int_high, ls_int_low, Smr1, Smr2))
}

height_ba <- function(x, xval) {
  f <- bearing_area(x)

  val <- (1 - quantile(f, probs = c(xval))[[1]])

  return(val)
}

sdc <- function(x, low, high) {
  val_low <- height_ba(x, low)
  val_high <- height_ba(x, high)

  val <- val_low - val_high

  return(val)
}

# surface bearing index = ratio of Sq to height from top of surface to height at
# 5% of bearing area (z05)
sbi <- function(x) {
  Sq <- sq(x)
  z05 <- height_ba(x, 0.05)

  val <- Sq / z05

  return(val)
}

# valley fluid retention index = void volume (area under Abbott curve) in 'valley' zone
# see fig.2a from Kedron et al. (2018)
svi <- function(x) {
  f <- bearing_area(x)

  val <- simpsons_above(model = f, b = 1, a = 0.8, n = 500)

  return(val)
}

# core fluid retention index = void volume (area above Abbott curve) in the core zone
# see fig.2a from Kedron et al. (2018)
sci <- function(x) {
  f <- bearing_area(x)

  core_above <- simpsons_above(model = f, b = 1, a = 0.05, n = 1000)
  Svi <- svi(x)
  val <- core_above - Svi

  return(val)
}

# core roughness depth = height difference between the intersection points of the found least
# mean square line in the Abbott curve (see fig. 2b from Kedron et al. 2018)
sk <- function(x) {
  line_info <- find_flat(x, perc = 0.4)

  ls_int_high <- line_info[[3]]
  ls_int_low <- line_info[[4]]

  val <- abs(ls_int_high - ls_int_low)

  return(val)
}

# reduced valley depth = height of triangle drawn at 100% on Abbott curve
svk <- function(x) {
  f <- bearing_area(x)

  line_info <- find_flat(x, perc = 0.4)

  smr2 <- line_info[[5]]

  val <- abs((1 - quantile(f, probs = 1)) - (1 - quantile(f, probs = smr2)))[[1]]

  return(val)
}

# reduced peak height = height of upper left triangle in abbott curve
spk <- function(x) {
  f <- bearing_area(x)

  line_info <- find_flat(x, perc = 0.4)

  smr1 <- line_info[[4]]

  val <- abs((1 - quantile(f, probs = 0)) - (1 - quantile(f, probs = smr1)))[[1]]

  return(val)
}
