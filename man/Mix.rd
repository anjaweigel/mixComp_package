\name{Mix}
\alias{Mix}
\alias{print.Mix}
\alias{is.Mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixtures of Univariate Distributions
}
\description{
Objects of class \code{Mix} represent finite mixtures of any univariate distribution. Methods for construction, printing and plotting are provided.
}
\usage{
Mix(dist, w = NULL, param.list = NULL, name = NULL, \dots)
  
is.Mix(x)
  
\method{print}{Mix}(x, \dots)  

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{dist}{a character string giving the (abbreviated) name of the component distribution.
    Adding "d" or "r" to the beginning of the string has to give the name of the
    density function and the random number generation function corresponding to the
    component distribution.}
  
  \item{w}{numeric vector of length \eqn{p}, say, specifying the mixture weights \eqn{w[i]}
  of the components, \eqn{i = 1,\dots,p}. If the weights don't add up to 1, they will be scaled
  accordingly. Defaults to equal proportions.}
  
  \item{param.list}{named list specifying the parameter values of the \eqn{p} components. 
    The names have to match the formals of the density function and the random number
    generation function exactly. Alternatively, the parameter values can be supplied directly 
    as named vectors of length \eqn{p} via \dots.}
    
  \item{name}{optional name tag of the result (used for printing and plotting).}  

  \item{x}{
    \describe{
      \item{in \code{is.Mix()}:}{R object.}
      \item{in \code{print.Mix()}:}{object of class \code{Mix}.}
  }}

  \item{\dots}{further arguments passed to the print method.}
}  
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{Mix} (implemented as a matrix) with the following attributes:
  \item{dim}{dimensions of the matrix.}
  \item{dimnames}{a \code{\link[base]{dimnames}} attribute for the matrix.}
  \item{name}{\code{name} (see above)}
  \item{dist}{\code{dist} (see above)}
  \item{discrete}{logical indicating whether the mixture distribution is discrete.}
  \item{param.list}{\code{param.list} (see above)}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dMix}} for the density, \code{\link{rMix}} for random numbers (and construction of an \code{rMix} object) and \code{\link{plot.Mix}} for the plot method.
}
\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))

# plot 'Mix' object
plot(normLocMix)
plot(poisMix)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
