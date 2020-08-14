\name{dMix}
\alias{dMix}

\title{
Mixture density
}

\description{
Evaluate the (log) density function of a mixture specified as \code{Mix} object.
}

\usage{
dMix(x, obj, log = FALSE)
}

\arguments{
  \item{x}{vector of quantiles.}

  \item{obj}{object of class \code{Mix}.}
  
  \item{log}{logical; if \code{TRUE}, probabilities/densities \eqn{f} are returned as \eqn{log(f)}.}
}

\value{
\code{dMix(x)} returns the numeric vector of probability values \eqn{f(x)}, logged if \code{log} is \code{TRUE}.
}

\seealso{
\code{\link{Mix}} for the construction of \code{Mix} objects, 
\code{\link{rMix}} for random number generation (and construction of an \code{rMix} object) and 
\code{\link{plot.Mix}} which makes use of \code{dMix}.
}

\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))

# evaluate density at points x
x <- seq(7, 20, length = 501)
dens <- dMix(x, normLocMix)
plot(x, dens, type = "l")

# compare to plot.Mix
plot(normLocMix)
}

\keyword{cluster}