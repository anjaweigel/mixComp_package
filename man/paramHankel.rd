\name{paramHankel}
\alias{paramHankel}
\alias{paramHankel.scaled}
\alias{print.paramEst}
\alias{plot.paramEst}

\title{
Estimate a Mixture's Complexity (and Component Weights/Parameters) Based on Hankel Matrix
}

\description{
Estimation of a mixture's complexity as well as component weights and parameters based on estimating the determinant of the Hankel matrix of the moments of the mixing distribution and comparing it to determinant values generated via parametric bootstrap.
}

\usage{
paramHankel(obj, j.max = 10, B = 1000, ql = 0.025, 
            qu = 0.975, control = c(trace = 0), \dots)

paramHankel.scaled(obj, j.max = 10, B = 100, ql = 0.025, 
                   qu = 0.975, control = c(trace = 0), \dots)
                   
\method{print}{paramEst}(x, \dots)
                   
\method{plot}{paramEst}(x, mixture = TRUE, components = TRUE, ylim = NULL, 
     cex.main = 0.9, \dots)
}

\arguments{
  \item{obj}{object of class \code{datMix}.}
  
  \item{j.max}{integer stating the maximal number of components to be considered.}  
  
  \item{B}{integer specifying the number of bootstrap replicates.}
  
  \item{ql}{numeric between \eqn{0} and \eqn{1} specifying the lower bootstrap quantile to which the observed determinant value will be compared.}
    
  \item{qu}{numeric between \eqn{0} and \eqn{1} specifying the upper bootstrap quantile to which the observed determinant value will be compared.}  
    
  \item{control}{control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.}  
  \item{x}{object of class \code{paramEst}.}
  
  \item{mixture}{logical indicating whether the estimated mixture density should be drawn, set to \code{TRUE} by default.}

  \item{components}{logical indicating whether the individual mixture components should be drawn, set to \code{TRUE} by default.}
  
  \item{ylim}{range of y values to use; if not specified (or containing \code{NA}), the function tries to construct reasonable default values itself.}
    
  \item{cex.main}{The magnification to be used for main titles relative to the current setting of \code{cex}, see \code{\link[graphics]{par}}.} 
    
  \item{\dots}{
    \describe{
      \item{in \code{paramHankel()} and \code{paramHankel.scaled()}:}{further arguments passed to the \code{\link[boot]{boot}} function.}
      \item{in \code{plot.hankDet()}:}{further arguments passed to the \code{\link[graphics]{hist}} function plotting the data.}
            \item{in \code{print.hankDet()}:}{further arguments passed to the \code{\link[stats]{printCoefmat}} function.}
    }}
}

\details{
Define \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as

\eqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p)}.

The \code{paramHankel} procedure initially assumes the mixture to only contain a single component, setting \eqn{j = 1}, and then sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}, until the algorithm terminates. To do so, it determines the MLE for a \eqn{j}-component mixture, generates \code{B} parametric bootstrap samples of size \eqn{n} from the distribution the MLE corresponds to and calculates \code{B} determinants of the corresponding \eqn{(j+1)x(j+1)} Hankel matrices of the first \eqn{2j} raw moments of the mixing distribution (for details see \code{\link{nonparamHankel}}). The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the determinant value based on the original data lies outside of the interval \eqn{[ql, qu]}, a range specified by the \code{ql} and \code{qu} empirical quantiles of the bootstrapped determinants. Otherwise, \eqn{j} is returned as the complexity estimate. 

\code{paramHankel.scaled} functions similarly to \code{paramHankel} with the exception that the bootstrapped determinants are scaled by the empirical standard deviation of the bootstrap sample. To scale the original determinant, \code{B} nonparametric bootstrap samples of size \eqn{n} are generated from the data, the corresponding determinants are calculated and their empirical standard deviation is used.

The MLEs are calculated via the \code{MLE.function} attribute (of the \code{datMix} object \code{obj}) for \eqn{j = 1}, if it is supplied. For all other \eqn{j} (and also for \eqn{j = 1} in case \code{MLE.function = NULL}) the solver \code{\link[Rsolnp]{solnp}} is used to calculate the minimum of the negative log likelihood. The initial values supplied to the solver are calculated as follows: the data is clustered into \eqn{j} groups by the function \code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function} (if supplied to the \code{datMix} object, otherwise numerical optimization is used here as well). The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates.
}

\value{
Object of class \code{paramEst} with the following attributes

\item{dat}{data based on which the complexity is estimated.}

\item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}

\item{ndistparams}{integer specifying the number of parameters identifying the component distribution, i.e. if \eqn{\theta \subseteq R^d} then \code{ndistparams}\eqn{ = d}.}

\item{formals.dist}{string vector specifying the names of the formal arguments identifying the distribution \code{dist} and used in \code{ddist} and \code{rdist}, e.g. for a gaussian mixture (\code{dist = norm}) amounts to \code{mean} and \code{sd}, as these are the formal arguments used by \code{dnorm} and \code{rnorm}.}

\item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}

\item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}

\item{pars}{Say the complexity estimate is equal to some \eqn{j}. Then \code{pars} is a numeric vector of size \eqn{(d+1)*j-1} specifying the component weight and parameter estimates, given as 

\eqn{(w_1, ... w_{j-1}, \theta 1_1, ... \theta 1_j, \theta 2_1, ... \theta d_j)}.}

\item{values}{numeric vector of function values gone through during optimization at iteration \eqn{j}, the last entry being the value at the optimum.}

\item{convergence}{indicates whether the solver has converged (0) or not (1 or 2) at iteration \eqn{j}.}
}

\seealso{
\code{\link{nonparamHankel}} for estimation of the mixture complexity based on the Hankel
matrix without parameter estimation, 
\code{\link[Rsolnp]{solnp}} for the solver, 
\code{\link{datMix}} for creation of the \code{datMix} object.
}

\examples{
## create 'Mix' object
poisMix <- Mix("pois", w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))


## create random data based on 'Mix' object (gives back 'rMix' object)
set.seed(1)
poisRMix <- rMix(1000, obj = poisMix)


## create 'datMix' object for estimation

# generate list of parameter bounds
poisList <- vector(mode = "list", length = 1)
names(poisList) <- "lambda"
poisList$lambda <- c(0, Inf)

# generate MLE function
MLE.pois <- function(dat){
  mean(dat)
}

# generate function needed for estimating the j^th moment of the 
# mixing distribution via Hankel.method "explicit"

explicit.pois <- function(dat, j){
  res <- 1
  for (i in 0:(j-1)){
    res <- res*(dat-i)
  }
  return(mean(res))
}
        
# generating 'datMix' object
pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois,
                  Hankel.method = "explicit", Hankel.function = explicit.pois)


## complexity and parameter estimation
set.seed(1)
res <- paramHankel(pois.dM)
plot(res)
}

\keyword{cluster}
