\name{datMix}
\alias{datMix}
\alias{is.datMix}
\alias{print.datMix}

\title{
Create Object for Which to Estimate the Mixture Complexity
}

\description{
Function to generate a \code{datMix} object to be passed to other \code{mixComp} functions used for estimating the mixture complexity.
}

\usage{
datMix(dat, dist, theta.bound.list = NULL, MLE.function = NULL, 
       Hankel.method = NULL, Hankel.function = NULL)

is.datMix(x)
       
\method{print}{datMix}(x, \dots)
}

\arguments{
  \item{dat}{a numeric vector containing the observations from the mixture model.}

  \item{dist}{a character string giving the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers. For example, to create a gaussian mixture, \code{dist} has to be specified as \code{norm} instead of \code{normal}, \code{gaussian} etc. for the package to find the functions \code{dnorm} and \code{rnorm}.}
    
  \item{theta.bound.list}{a named list specifying the upper and the lower bound for the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly. For a gaussian mixture, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm}. Has to be supplied if a method that estimates the component weights and parameters is to be used.}
  
  \item{MLE.function}{function (or list of functions) which takes as input the data and gives as output the maximum likelihood estimator for the parameter(s) of a one component mixture (i.e. the standard MLE of the component distribution \code{dist}).  If the component distribution has more than one parameter, a list of functions has to be supplied and the order of the MLE functions has to match the order of the component parameters in \code{theta.bound.list} (e.g. for a normal mixture, if the first entry of \code{theta.bound.list} is the bounds of the mean, then then first entry of \code{MLE.function} has to be the MLE of the mean). If this argument is supplied and the \code{datMix} object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the \code{MLE.function} attribute will be used for the single component case. In case the objective function is either not a likelihood or corresponds to a mixture with more than 1 components, numerical optimization will be used based on \code{\link{Rsolnp}}'s function \code{solnp}, but \code{MLE.function} will be used to calculate the initial values passed to \code{solnp}. Specifying \code{MLE.function} is optional and if it is not, for example because the MLE solution does not exists in closed form, numerical optimization is used to find the relevant MLE's.}

  \item{Hankel.method}{character string in \code{c("explicit", "translation", "scale")},  specifying the method of estimating the moments of the mixing distribution used to calculate the relevant Hankel matrix. Has to be specified when using \code{nonparamHankel}, \code{paramHankel} or \code{paramHankel.scaled}. For further details see below.}
  
  \item{Hankel.function}{function needed for the moment estimation via \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as \code{dist}. For further details see below.}
    
  \item{x}{
    \describe{
      \item{in \code{is.datMix()}:}{R object.}
      \item{in \code{print.datMix()}:}{object of class \code{datMix}.}
  }}
  
  \item{\dots}{further arguments passed to the print method.}
}

\details{
If the \code{datMix} object is supposed to be passed to a function that calculates the Hankel matrix of the moments of the mixing distribution (i.e. \code{\link{nonparamHankel}}, \code{\link{paramHankel}} or \code{\link{paramHankel.scaled}}), the arguments \code{Hankel.method} and \code{Hankel.function} have to be specified. The \code{Hankel.method}s that can be used to generate the estimate of the (raw) moments of the mixing distribution and the corresponding \code{Hankel.function}s are the following, where \eqn{j} specifies an estimate of the number of components:

\describe{
\item{\code{"explicit"}}{For this method, \code{Hankel.function} contains a function with arguments called \code{dat} and \code{j}, explicitly estimating the moments of the mixing distribution from the data and the currently assumed mixture complexity. Note that what Dacunha-Castelle & Gassiat (1997) called the "natural" estimator in their original paper is equivalent to using \code{"explicit"} with \code{Hankel.function} \eqn{f_j((1/n) * sum_i(\psi_j(X_i)))}.}

\item{\code{"translation"}}{This method corresponds to Dacunha-Castelle & Gassiat's (1997) example 3.1. It is applicable if the family of component distributions \eqn{(G_\theta)} is given by \eqn{dG_\theta(x) = dG(x-\theta)}, where \eqn{G} is a known probability distribution whose moments can be given explicitly. \code{Hankel.function} contains a function of \eqn{j} returning the \eqn{j}th (raw) moment of \eqn{G}.}

\item{\code{"scale"}}{This method corresponds to Dacunha-Castelle & Gassiat's (1997) example 3.2. It is applicable if the family of component distributions \eqn{(G_\theta)} is given by \eqn{dG_\theta(x) = dG(x\\theta)}, where \eqn{G} is a known probability distribution whose moments can be given explicitly. \code{Hankel.function} contains a function of \eqn{j} returning the \eqn{j}th (raw) moment of \eqn{G}.}
}
                                               
If the \code{datMix} object is supposed to be passed to a function that estimates the component weights and parameters (i.e. all but \code{\link{nonparamHankel}}), the argument \code{theta.bound.list} has to be specified, and \code{MLE.function} will be used in the estimation process if it is supplied (otherwise the MLE is found numerically).

Note that the \code{datMix} function will change the random number generator (RNG) state.
}

\value{
An object of class \code{datMix} with the following attributes (for further explanations
see above):
  \item{dist}{}
  \item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}
  \item{theta.bound.list}{}
  \item{MLE.function}{}
  \item{Hankel.method}{}
  \item{Hankel.function}{}
}

\seealso{
\code{\link{RtoDat}} for the conversion of \code{rMix} to \code{datMix} objects.
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
# for "sd" (the sd function uses (n-1) as denominator)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
# combining the functions to a list
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                      "MLE.norm.sd" = MLE.norm.sd)

## function giving the j^th raw moment of the standard normal distribution,
## needed for calculation of the Hankel matrix via the "translation" method
## (assuming gaussian components with variance 1)

mom.std.norm <- function(j){
  ifelse(j \%\% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

        
## generate 'datMix' object
faithful.dM <- datMix(obs, dist = "norm", theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = "translation",
                      Hankel.function = mom.std.norm)
                      
## using 'datMix' object to estimate the mixture complexity
set.seed(1)
res <- paramHankel.scaled(faithful.dM)
plot(res)
}

\keyword{cluster}
