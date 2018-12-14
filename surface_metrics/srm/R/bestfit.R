# functions to find the best fit line for forty percent of slope with lowest slope
slopecalc <- function(x, h, f) {
  xplus <- even_x + h
  xminus <- even_x - h
  x <- even_x

  # figure out ends
  space <- 1 / length(x)
  end_length <- h / space
  if (end_length < 1) { end_length = 1}
  begin <- 1 + end_length
  end <- length(x) - end_length
  begin1 <- begin - 1
  end1 <- end + 1

  fxh_pos <- seq(1, length(x))
  fxh_neg <- seq(1, length(x))
  fxh_pos[1:begin1] <- (1 - quantile(f, probs = xplus[1:begin1]))
  fxh_neg[1:begin1] <- (1 - quantile(f, probs = x[1:begin1]))
  # variation on newton's difference quotient at far end
  fxh_pos[end1:length(x)] <- (1 - quantile(f, probs = x[end1:length(x)]))
  fxh_neg[end1:length(x)] <- (1 - quantile(f, probs = xminus[end1:length(x)]))
  # symmetric difference quotient everywhere else (99800 points)
  fxh_pos[begin:end] <- (1 - quantile(f, probs = xplus[begin:end]))
  fxh_neg[begin:end] <- (1 - quantile(f, probs = xminus[begin:end]))

  diff_quo <- (fxh_pos - fxh_neg) / (2 * h)
  slopes <- data.frame(slope = diff_quo, x = x)

  return(slopes)
}

# calculate averages for each 40% segment
slopemeans <- function(slope_data, segment_perc = 0.4) {
  x <- slope_data$x
  slope <- slope_data$slope
  length <- length(x) * segment_perc
  end_ind <- length(x) - length
  begin_ind <- 1 + length
  xstart <- x[1:end_ind]
  istart <- seq(1, end_ind)
  xend <- x[begin_ind:length(x)]
  iend <- seq(begin_ind, length(x))

  # create matrix of slopes
  xmat <- matrix(c(xstart, xend), ncol = 2)
  mean <- sapply(seq(1, end_ind), function(i) {mean(abs(slope[istart[i]:iend[i]]))})

  data <- data.frame(slope = mean, xstart = xstart, xend = xend)
  return(data)
}
