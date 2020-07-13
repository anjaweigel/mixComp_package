\name{plot.Mix}
\alias{plot.Mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{\code{plot} Method for \code{\link{Mix}} Objects}

\description{
  \code{plot} method for \code{\link{Mix}} objects drawing the mixture density, optionally additonally showing the component densities.
}
\usage{
\method{plot}{Mix}(x, ylim, xlim = NULL, xout = NULL, n = 511, type = NULL,
     xlab = "x", ylab = "f(x)", main = attr(obj, "name"), 
     lwd = 1.4, log = FALSE, components = TRUE, h0 = FALSE, 
     parComp = list(col= NULL, lty = 3, lwd = 1),
     parH0   = list(col= NULL, lty = 3, lwd = 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{Mix}.}
  
  \item{ylim}{range of y values to use; if not specified (or
    containing \code{NA}), the function tries to construct reasonable 
    default values.}
  
  \item{xlim}{range of x values to use; particularly important if
    \code{xout} is not specified. If both are left unspecified, the
    function tries to construct reasonable default values.}
  
  \item{xout}{numeric or \code{NULL} giving the abscissae at which to
    draw the density.}
  
  \item{n}{number of points to generate if \code{xout} is unspecified (for continous
           distributions).}
    
  \item{type}{character denoting type of plot, see, e.g. \code{\link[graphics]{lines}}. 
    Defaults to \code{"l"} if the mixture distribution is continuous and to 
    \code{"h"} otherwise.}
    
  \item{xlab,ylab}{labels for the x and y axis with defaults.}
  
  \item{main}{main title of plot, defaulting to the \code{\link{Mix}}
    name.}

  \item{lwd}{line width for plotting with a non-standard default.}
  
  \item{log}{logical; if \code{TRUE}, probabilities/densities \eqn{p} are plotted
     as \eqn{log(p)}. Only works if \code{components} is false.}
  
  \item{h0}{logical indicating whether the line \eqn{y = 0} should be drawn.}
  
  \item{components}{logical indicating whether the individual mixture components should
    be plotted, defaulting to \code{TRUE}.}
    
  \item{parH0}{graphical parameters for drawing the line \eqn{y = 0} if
    \code{h0} is true.}
    
  \item{parComp}{graphical parameters for drawing the individual components if
    \code{components} is true.}
    
  \item{\dots}{further arguments passed to the function plotting the mixture
    density.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{Mix}} for the construction of \code{Mix} objects, \code{\link{dMix}} for 
the density of a mixture.
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
