\name{rMix}
\alias{rMix}
\alias{is.rMix}
\alias{print.rMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Random Numbers Following a Mixture Distribution
}
\description{
Generate \code{n} random numbers, distributed according to a mixture specified as \code{Mix} object. Gives back an object of class \code{rMix}.
}
\usage{
rMix(n, obj)

is.rMix(x)

\method{print}{rMix}(x, \dots)  
}
%- maybe also 'usage' for other objects documented here.
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
For a mixture of \eqn{p} components, generate the number of observations in each component as multinomial, and then use 
a standard random number generation function for each component.
The integer (multinomial) numbers are generated via \code{\link[base]{sample}}.
}

\value{
An object of class \code{rMix} with the following attributes (for further explanations
see \code{\link{Mix}}):
  \item{name}{name of the \code{Mix} object that was given as input.}
  \item{dist}{(abbreviated) name of the  component distribution.}
  \item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}
  \item{param.list}{named list specifying the parameter values of the \eqn{p} components.}
  \item{w}{numeric vector of length \eqn{p}, specifying the mixture weights \eqn{w[i]}
    of the components, \eqn{i = 1,\dots,p}.}
  \item{indices}{numeric vector of length \code{n} with only integers between 
    \eqn{1} and \eqn{p} specifying which mixture component each observation belongs to.}
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
\code{\link{dMix}} for the density and \code{\link{Mix}} for the construction of \code{Mix} objects.
}
\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))

# generate n random samples
x <- rMix(1000, normLocMix)
hist(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
