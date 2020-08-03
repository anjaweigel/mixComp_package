\name{Mix}
\alias{Mix}
\alias{print.Mix}
\alias{is.Mix}

\title{
Mixtures of Univariate Distributions
}

\description{
Objects of class \code{Mix} represent finite mixtures of any univariate distribution. Methods for construction, printing and plotting are provided.
}

\usage{
Mix(dist, w = NULL, theta.list = NULL, name = NULL, \dots)
  
is.Mix(x)
  
\method{print}{Mix}(x, \dots)  

}

\arguments{
  \item{dist}{a character string giving the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers. For example, to create a gaussian mixture, \code{dist} has to be specified as \code{norm} instead of \code{normal}, \code{gaussian} etc. for the package to find the functions \code{dnorm} and \code{rnorm}.}
  
  \item{w}{numeric vector of length \eqn{p}, specifying the mixture weights \eqn{w[i]} of the components, \eqn{i = 1,\dots,p}. If the weights don't add up to 1, they will be scaled accordingly. Uses equal weights for all components by default.}
  
  \item{theta.list}{named list specifying the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly. For a gaussian mixture, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm}. Alternatively, the component parameters can be supplied directly as named vectors of length \eqn{p} via \dots}
    
  \item{name}{optional name tag of the result (used for printing and plotting).}  

  \item{x}{
    \describe{
      \item{in \code{is.Mix()}:}{R object.}
      \item{in \code{print.Mix()}:}{object of class \code{Mix}.}
  }}

  \item{\dots}{further arguments passed to the print method.}
}  

\details{
Note that the \code{Mix} function will the random number generator (RNG) state.
}

\value{
An object of class \code{Mix} (implemented as a matrix) with the following attributes:
  \item{dim}{dimensions of the matrix.}
  \item{dimnames}{a \code{\link[base]{dimnames}} attribute for the matrix.}
  \item{name}{as entered via \code{name}.}
  \item{dist}{as entered via \code{dist}.}
  \item{discrete}{logical indicating whether the mixture distribution is discrete.}
  \item{theta.list}{as entered via \code{theta.list}.}
}

\seealso{
\code{\link{dMix}} for the density, 
\code{\link{rMix}} for random numbers (and construction of an \code{rMix} object) and
\code{\link{plot.Mix}} for the plot method.
}

\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))

# plot 'Mix' object
plot(normLocMix)
plot(poisMix)
}

\keyword{cluster}
