\name{dMix}
\alias{dMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture density
}
\description{
Evaluate the (log) density function of a mixture specified as \code{Mix} object.
}
\usage{
dMix(x, obj, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles.}

  \item{obj}{object of class \code{Mix}.}
  
  \item{log}{logical; if \code{TRUE}, probabilities/densities \eqn{p} are returned as 
    \eqn{log(p)}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{dMix(x)} returns the numeric vector of density/probability values \eqn{f(x)}, logged if \code{log} is \code{TRUE}.
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
\code{\link{Mix}} for the construction of \code{Mix} objects, \code{\link{rMix}} for 
random number generation (and construction of an \code{rMix} object) and \code{\link{plot.Mix}} which makes use of \code{dMix}.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}