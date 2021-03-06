% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localsurface.R
\name{s10z}
\alias{s10z}
\title{Ten-Point Height}
\usage{
s10z(x)
}
\arguments{
\item{x}{A raster object.}
}
\value{
A numeric value representing the ten-point height.
}
\description{
Calculates the average height abover the mean surface for the five highest local maxima
plus the average height below the mean surface for the five lowest local minima.
}
\examples{
# import raster image
data(normforest)

# calculate ten-point height.
S10z <- s10z(normforest)
}
