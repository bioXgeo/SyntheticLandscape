% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{remove_plane}
\alias{remove_plane}
\title{Removes the Best Fit Polynomial Plane from a Raster}
\usage{
remove_plane(x)
}
\arguments{
\item{x}{A raster.}
}
\value{
A raster of the same size as the input with values
  equal to the difference between the original and bestfit
  plane rasters.
}
\description{
Finds the best fit polynomial plane for a raster image and
subtracts it from the actual raster values. The remaining
raster has positive values where the actual values are higher
than the plane and negative values where the actual value
are lower than the plane.
}
\examples{
library(raster)

# import raster image
data(orforest)

# remove the least squares polynomial plane
new_rast <- remove_plane(orforest)

# plot
plot(new_rast)
}
