\name{datMix}
\alias{datMix}
\alias{is.datMix}
\alias{print.datMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Object for Which to Estimate the Mixture Complexity
}
\description{
Function to generate a \code{datMix} object to be passed to other \code{mixcomp} functions estimating the mixture complexity (in some cases along with the mixture parameters).
}
\usage{
datMix(dat, dist, param.bound.list = NULL, MLE.function = NULL, 
       Hankel.method = NULL, Hankel.function = NULL)

is.datMix(x)
       
\method{print}{datMix}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{a numeric vector containing the observations from the mixture model.}

  \item{dist}{a character string giving the (abbreviated) name of the component distribution.
    Adding "d" or "r" to the beginning of the string has to give the name of the
    density function and the random number generation function corresponding to the
    component distribution.}
    
  \item{param.bound.list}{a named list containing the upper and the lower bound for every
    parameter of the component distribution. The names have to match the argument names of the  
    density function and the random number generation function exactly. Has to be supplied
    if methods that estimate the parameter values are to be used.}
    
  \item{MLE.function}{function (or list of functions) which takes as input the data and gives 
    as output the maximum likelihood estimator for the parameter(s) of a one component 
    mixture (i.e. the standard MLE for the component distribution \code{dist}). Passed to all
    functions which optimize over the likelihood to get analytical instead of numeric results 
    for \eqn{j = 1} (\eqn{j} being the current estimate of the model complexity), and used for 
    parameter initialization in functions estimating distribution parameters (i.e. all but \code{nonparamHankel}). In a list, 
    the order of the MLE functions has to match the order of the parameters in 
    \code{param.bound.list}. Numerical optimization is used if this is not supplied.}
    
  \item{Hankel.method}{character string in \code{c("natural", "explicit", "translation", 
    "scale")}. Has to be specified if the \code{datMix} object is to be passed to a function 
    that calculates the Hankel matrix to determine how to estimate the moments of the 
    mixing distribution. For further details see below.}
  
  \item{Hankel.function}{function (or list of functions) needed for the moment estimation via 
    \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as 
    \code{dist}. For further details see below.}
    
  \item{x}{
    \describe{
      \item{in \code{is.datMix()}:}{R object.}
      \item{in \code{print.datMix()}:}{object of class \code{datMix}.}
  }}
  
  \item{\dots}{further arguments passed to the print method.}
}
\details{

If the \code{datMix} object is supposed to be passed to a function that calculates the Hankel matrix 
(i.e. \code{\link{nonparamHankel}}, \code{\link{paramHankel}} or \code{\link{paramHankel.scaled}}), 
the arguments \code{Hankel.method} and \code{Hankel.function} have to be specified. The 
\code{Hankel.method}s that can be used to generate the estimate of the (raw) moments of the mixing
distribution and the corresponding \code{Hankel.function}s are the following (here the subscript \eqn{j}
is used (instead of \eqn{p} as in the original paper) to be consistent in the package documentation):

\describe{
\item{\code{"natural"}} {see Dacunha-Castelle & Gassiat (1997), page 283 equation (3). For this method, 
                        the functions \eqn{\psi_j} and \eqn{f_j} have to be supplied as a list to \code{Hankel.function},
                        with \eqn{\psi_j} as first element and \eqn{f_j} as second element. The function \eqn{\psi_j} contains the 
                        data vector as first argument and \eqn{j} as second, and gives back the vector
                        \eqn{\psi_j(X_i), 1 <= i <= n}. \eqn{f_j} contains the average of \eqn{\psi_j(X_i)} as first argument 
                        and \eqn{j} as second (even if it is unused in the function body). If a single function is given as
                        input this will be taken as \eqn{\psi_j}, and \eqn{f_j} will be taken to be the identity function.}


\item{\code{"explicit"}} {For this method, \code{Hankel.function} contains a function which explicitly estimates the moments of 
                         the mixing distribution. Note that \code{"natural"} is equivalent to using \code{"explicit"} with
                         \code{Hankel.function} 
                         \eqn{f_j((1/n) * sum_i(\psi_j(X_i)))} (i.e. \code{f(mean(psi(dat, j)), j)}).}


\item{\code{"translation"}} {see Dacunha-Castelle & Gassiat (1997), page 284 example 3.1. For this method, \code{Hankel.function}
                            contains the function giving back \eqn{E[Y^j]} (i.e the raw moments of the random variable Y) as a 
                            function of \code{j}, where Y is defined as in 3.1. (4).}

\item{\code{"scale"}} {see Dacunha-Castelle & Gassiat (1997), page 285 example 3.2. For this method, \code{Hankel.function} contains
                      the function giving back \eqn{E[Y^j]} (i.e the raw moments of the random variable Y) as a function of 
                      \code{j}, where Y is defined as in the second equation of 3.2 (\eqn{X_t = \sigma_t * Y_t}).}
}
                                               
If the \code{datMix} object is supposed to be passed to a function that estimates the distribution parameters 
(i.e. all but \code{\link{nonparamHankel}}), the argument \code{param.bound.list} has to be specified,
and \code{MLE.function} is used instead of numerical calculation if supplied.
                                               
}
\value{
An object of class \code{datMix} with the following attributes (for further explanations
see above):
  \item{dist}{}
  \item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}
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
\code{\link{RtoDat}} for the conversion of \code{RMix} to \code{datMix} objects.
}
\examples{
## observations from a (presumed) mixture model
obs <- faithful$waiting

## generate list of parameter bounds (assuming gaussian components)
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

        
## generate 'datMix' object
faithful.dM <- datMix(obs, dist = "norm", param.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = "translation",
                      Hankel.function = mom.std.norm)
                      
## using 'datMix' object to estimate the mixture complexity
set.seed(1)
res <- paramHankel.scaled(faithful.dM)
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
