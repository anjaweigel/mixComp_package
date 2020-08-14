\name{rMix}
\alias{rMix}
\alias{is.rMix}
\alias{print.rMix}

\title{
Generate a Random Sample from a Mixture Distribution
}

\description{
Generate a random sample of size \code{n}, distributed according to a mixture specified as \code{Mix} object. Returns an object of class \code{rMix}.
}

\usage{
rMix(n, obj)

is.rMix(x)

\method{print}{rMix}(x, \dots)  
}

\arguments{
  \item{n}{integer specifying the number of observations.}
    
  \item{obj}{object of class \code{Mix}.}
  
  \item{x}{
  \describe{
      \item{in \code{is.rMix()}:}{R object.}
      \item{in \code{print.rMix()}:}{object of class \code{rMix}.}
  }}
      
  \item{\dots}{further arguments passed to the print method.}
}

\details{
For a mixture of \eqn{p} components, generate the number of observations in each component as multinomial, and then use an implemented random variate generation function for each component. The integer (multinomial) numbers are generated via \code{\link[base]{sample}}.
}

\value{
An object of class \code{rMix} with the following attributes (for further explanations
see \code{\link{Mix}}):
  \item{name}{name of the \code{Mix} object that was given as input.}
  \item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}
  \item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}
  \item{theta.list}{named list specifying the parameter values of the \eqn{p} components.}
  \item{w}{numeric vector of length \eqn{p}, specifying the mixture weights \eqn{w[i]} of the components, \eqn{i = 1,\dots,p}.}
  \item{indices}{numeric vector of length \code{n} containing integers between \eqn{1} and \eqn{p} specifying which mixture component each observation belongs to.}
}

\seealso{
\code{\link{dMix}} for the density, 
\code{\link{Mix}} for the construction of \code{Mix} objects and
\code{\link{plot.rMix}} for the plot method.
}

\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))

# generate n random samples
set.seed(1)
x <- rMix(1000, normLocMix)
hist(x)
}

\keyword{cluster}
