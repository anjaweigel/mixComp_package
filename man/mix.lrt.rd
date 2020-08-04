\name{mix.lrt}
\alias{mix.lrt}

\title{
Estimate a Mixture's Complexity Based on Likelihood Ratio Test Statistics
}

\description{
Estimation of a mixture's the complexity as well as its component weights and parameters based on comparing the likelihood ratio test statistic (LRTS) to a bootstrapped quantile.
}

\usage{
mix.lrt(obj, j.max = 10, B = 100, quantile = 0.95,
        control = c(trace = 0), ...)
}

\arguments{
  \item{obj}{object of class \code{datMix}.}
  
  \item{j.max}{integer giving the maximal complexity to be considered.}
  
  \item{B}{integer specifying the number of bootstrap replicates.}
    
  \item{quantile}{numeric between \eqn{0} and \eqn{1} specifying the bootstrap quantile to which the observed LRTS will be compared.}
  
  \item{control}{control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.}    
  \item{\dots}{further arguments passed to the \code{\link[boot]{boot}} function.}
}

\details{

Define the \eqn{complexity} of a finite mixture \eqn{F} as the smallest 
integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as

\eqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p)}.

To estimate \eqn{p}, \code{mix.lrt} sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}, by finding the maximum likelihood estimator (MLE) for the density of a mixture with \eqn{j} and \eqn{j+1} components and calculating the corresponding likelihood ratio test statistic (LRTS). Next, a parametric bootstrap procedure is used to generate \code{B} samples of size \eqn{n} from a \eqn{j}-component mixture given the previously calculated MLE. For each of the bootstrap samples, the MLEs corresponding to densities of mixtures with \eqn{j} and \eqn{j+1} components are calculated, as well as the LRTS. The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the LRTS based on the original data is larger than the chosen \code{quantile} of its bootstrapped counterparts. Otherwise, \eqn{j} is returned as the complexity estimate. 

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
\code{\link[Rsolnp]{solnp}} for the solver,
\code{\link{datMix}} for the creation of the \code{datMix} object.
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
# for "sd" (the sd function uses (n-1) as denominator)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
}
# combining the functions to a list
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                      "MLE.norm.sd" = MLE.norm.sd)

## generating 'datMix' object
normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list)
                
                      
### complexity and parameter estimation 
set.seed(0)
res <- mix.lrt(normLoc.dM, B = 30)
plot(res)
}

\keyword{cluster}

