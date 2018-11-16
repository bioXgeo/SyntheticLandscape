# simpson's rule area under curve for empirical functions

# modified for area above curve
simpsons_above <- function(model, a, b, n = 100) {
  # determines area above curve (just within rect defined by area below)
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)
  y <- quantile(model, probs = x) # get y-values of inverse cdf function
  s <- ((b - a) / (3 * n)) * (y[[1]] + 
                                sum(4 * y[seq(2, n - 1, by = 2)]) + 
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])
  # get inverse of s for actual area above curve
  area_above <- ((max(y) - min(y)) * (max(x) - min(x))) - s
  return(area_above)
}

# area under curve
simpsons_below <- function(model, a, b, n = 100) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b - a) / n # sub-interval width
  x <- seq(a, b, by = h)
  y <- quantile(model, probs = x) # get y-values of inverse cdf function
  s <- ((b - a) / (3 * n)) * (y[[1]] + 
                                sum(4 * y[seq(2, n - 1, by = 2)]) + 
                                sum(2 * y[seq(3, n - 1, by = 2)]) +
                                y[[n]])
  return(s)
}