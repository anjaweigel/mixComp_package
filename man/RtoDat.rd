\name{RtoDat}
\alias{RtoDat}

\title{
Converting \code{rMix} to \code{datMix} Objects
}

\description{
Converting an object of class \code{rMix} to an object of class \code{datMix}, so that it can be passed to functions estimating the mixture complexity.
}

\usage{
RtoDat(obj, theta.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL, 
       Hankel.function = NULL)
}

\arguments{
  \item{obj}{object of class \code{rMix}.}
    
  \item{theta.bound.list}{a named list specifying the upper and the lower bound for the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly (the \code{dist} attribute is specified when creating the \code{Mix} object \code{obj}). For a gaussian mixture with \code{dist = norm}, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm}. Has to be supplied if a method that estimates the component weights and parameters is to be used.}
    
  \item{MLE.function}{function (or list of functions) which takes as input the data and gives as output the maximum likelihood estimator for the parameter(s) of a one component mixture (i.e. the standard MLE of the component distribution \code{dist}).  If the component distribution has more than one parameter, a list of functions has to be supplied and the order of the MLE functions has to match the order of the component parameters in \code{theta.bound.list} (e.g. for a normal mixture, if the first entry of \code{theta.bound.list} is the bounds of the mean, then then first entry of \code{MLE.function} has to be the MLE of the mean). If this argument is supplied and the \code{datMix} object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the \code{MLE.function} attribute will be used for the single component case. In case the objective function is either not a likelihood or corresponds to a mixture with more than 1 component, numerical optimization will be used based on \code{\link{Rsolnp}}'s function \code{solnp}, but \code{MLE.function} will be used to calculate the initial values passed to \code{solnp}. Specifying \code{MLE.function} is optional and if it is not, for example because the MLE solution does not exists in closed form, numerical optimization is used to find the relevant MLE's.}

  \item{Hankel.method}{character string in \code{c("explicit", "translation", "scale")},  specifying the method of estimating the moments of the mixing distribution used to calculate the relevant Hankel matrix. Has to be specified when using \code{nonparamHankel}, \code{paramHankel} or \code{paramHankel.scaled}. For further details see the details section of \code{\link{datMix}}.}
  
  \item{Hankel.function}{function needed for the moment estimation via \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as \code{dist}. For further details see the details section of \code{\link{datMix}}.}
}

\value{
An object of class \code{datMix} with the following attributes (for further explanations
see \code{\link{datMix}}):
  \item{dist}{}
  \item{discrete}{}
  \item{theta.bound.list}{}
  \item{MLE.function}{}
  \item{Hankel.method}{}
  \item{Hankel.function}{}
}

\seealso{
\code{\link{datMix}} for direct generation of a \code{datMix} object from a vector of observations.
}

\examples{
### generating 'Mix' object
normLocMix <- Mix("norm", w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))

### generating 'rMix' from 'Mix' object (with 1000 observations)
set.seed(0)
normLocRMix <- rMix(1000, normLocMix)

### generating 'datMix' from 'R' object

## generate list of parameter bounds

norm.bound.list <- vector(mode = "list", length = 2)
names(norm.bound.list) <- c("mean", "sd")
norm.bound.list$mean <- c(-Inf, Inf)
norm.bound.list$sd <- c(0, Inf)

## generate MLE functions

# for "mean"
MLE.norm.mean <- function(dat) mean(dat)
# for "sd" (not using the sd function as it uses (n-1) as denominator)
MLE.norm.sd <- function(dat){
n <- length(dat)
var_hat <- (1/n)*sum((dat-mean(dat))^2)
sqrt(var_hat)
} 
# combining the functions to a list
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                      "MLE.norm.sd" = MLE.norm.sd)

## function giving the j^th raw moment of the standard normal distribution,
## needed for calculation of the Hankel matrix via the "translation" method
## (assuming gaussian components with variance 1)

mom.std.norm <- function(j){
  if(j \%\% 2 == 0){
    prod(seq(1, j-1, by = 2))
  }
  else 0
}


normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)
                      
### using 'datMix' object to estimate the mixture

set.seed(0)
res <- paramHankel.scaled(normLoc.dM)
plot(res)
}

\keyword{cluster}
