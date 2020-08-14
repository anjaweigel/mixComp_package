\name{plot.rMix}
\alias{plot.rMix}

\title{
  \code{plot} Method for \code{\link{rMix}} Objects
}

\description{
plot method for \code{rMix} objects, plotting the histogram of the random sample, with the option of additionally plotting the components (stacked or plotted over one another).
}

\usage{
\method{plot}{rMix}(x, xlab = attr(obj, "name"), ylim = NULL,
     main = paste("Histogram of", attr(obj, "name")), 
     breaks = NULL, col = "grey", components = TRUE, stacked = FALSE, 
     component.colors = NULL, freq = TRUE, plot = TRUE, ...)
}

\arguments{
  \item{x}{object of class \code{rMix}.}
  
  \item{xlab}{label for the x axis with default.}
  
  \item{ylim}{range of y values to use; if not specified (or containing \code{NA}), default values are used.}
  
  \item{main}{main title of the plot with default.}
  
  \item{breaks}{see \code{\link[graphics]{hist}}. If left unspecified the function tries to construct reasonable default values.}
    
  \item{col}{a colour to be used to fill the bars of the histogram evaluated on the whole data.}
  
  \item{components}{logical indicating whether the plot should show to which component the observations belong (either by plotting individual histograms or by overlaying a stacked barplot), defaulting to \code{TRUE}. Ignored if \code{plot} is \code{FALSE}.}
    
  \item{stacked}{logical indicating whether the component plots should be stacked or plotted over one another, defaulting to \code{FALSE}. Ignored if \code{components} is \code{FALSE} or ignored itself.}
    
  \item{component.colors}{the colors for the component plots. If left unspecified default colors are used.}
  
  \item{freq}{logical, if \code{TRUE}, the histogram graphic is a representation of frequencies, if \code{FALSE}, probability densities. See \code{\link[graphics]{hist}}.}
    
  \item{plot}{logical, if \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned. See \code{\link[graphics]{hist}}.}
    
  \item{\dots}{further arguments passed to the histogram function evaluated on the whole data as well as the component data (if \code{components} is \code{TRUE} and \code{stacked} is \code{FALSE}).}
}

\seealso{
\code{\link{rMix}} for the creation of \code{rMix} objects.
}

\examples{
# define 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))

# generate n random samples
set.seed(1)
x <- rMix(1000, normLocMix)
plot(x)
}

\keyword{cluster}
