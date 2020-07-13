\name{RtoDat}
\alias{RtoDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converting \code{rMix} to \code{datMix} Objects
}
\description{
Converting an object of class \code{rMix} to an object of class \code{datMix}, so that it can be passed to functions estimating the mixture complexity (in some cases along with the mixture parameters).
}
\usage{
RtoDat(obj, param.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL, 
       Hankel.function = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{object of class \code{rMix}.}
    
  \item{param.bound.list}{a named list containing the upper and the lower bound of every
    parameter of the component distribution. The names have to match the formals of the 
    density function and the random number generation function exactly. Has to be supplied
    if methods that estimate the parameter values are to be used.}
    
  \item{MLE.function}{function (or list of functions) which takes as input the data and gives 
    as output the maximum likelihood estimator for the parameter(s) of a one component mixture     
    (i.e. the standard MLE for the component distribution \code{dist}). Passed to all 
    functions which optimize over the likelihood to get analytical instead of numeric results 
    for \eqn{j = 1} (\eqn{j} being the current estimate of the model complexity),
    and used for parameter initialization in functions estimating distribtions parameters (i.e.
    all but \code{\link{nonparamHankel}}). In a list, the order of the MLE 
    functions has to match the order of the parameters in \code{param.bound.list}. Numerical
    optimization is used if this is not supplied.}
    
  \item{Hankel.method}{character string in \code{c("natural", "explicit", "translation", 
    "scale")}. Has to be specified if the \code{datMix} object is to be passed to a function 
    that calculates the Hankel matrix to determine how to estimate the moments of the mixing distribution.
    For further details, see the details section of \code{\link{datMix}}.}
  
  \item{Hankel.function}{function (or list of functions) needed for the moment estimation via 
    \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as 
    \code{dist}. For further details, see the details section of \code{\link{datMix}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{datMix} with the following attributes (for further explanations
see \code{\link{datMix}}):
  \item{dist}{}
  \item{discrete}{}
  \item{param.bound.list}{}
  \item{MLE.function}{}
  \item{Hankel.method}{}
  \item{Hankel.function}{}
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


normLoc.dM <- RtoDat(normLocRMix, param.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)
                      
### using 'datMix' object to estimate the mixture

set.seed(0)
res <- paramHankel.scaled(normLoc.dM)
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
