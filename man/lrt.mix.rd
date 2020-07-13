\name{mix.lrt}
\alias{mix.lrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Mixture Complexity Based on Likelihood Ratio Test Statistics
}
\description{
Estimation of the complexity and the component parameters of an underlying mixture based on comparing the likelihood ratio test statistic (LRTS) to a bootstrapped quantile.
}
\usage{
mix.lrt(obj, j.max = 10, B = 100, quantile = 0.95,
        control = c(trace = 0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{object of class \code{datMix}.}
  
  \item{j.max}{integer giving the maximal complexity to be considered.}
  
  \item{B}{integer specifying the number of bootstrap replicates.}
    
  \item{quantile}{numeric between \eqn{0} and \eqn{1} specifying the bootstrap quantile to which
    the observed value will be compared.}
  
  \item{control}{control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.}    
    
  \item{\dots}{further arguments passed to the \code{\link[boot]{boot}} function.}
}
\details{

Define the \eqn{order} or \eqn{complexity} of a finite mixture \eqn{F} as the smallest 
integer \eqn{p}, such that its density \eqn{f} can be written as

\eqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p)}.

To estimate \eqn{p}, \code{mix.lrt} iteratively increases the assumed order \eqn{j}, finds 
the maximum likelihood estimator (MLE) for both, the density of a mixture with \eqn{j} and \eqn{j+1}
components, and calculates the corresponding likelihood ratio test statistic (LRTS). Then, a 
parametric bootstrap is used to generate \code{B} samples of size \eqn{n} from a \eqn{j} component
mixture (given the previously calculated MLE). For each of the samples, again the MLEs corresponding
to densities with \eqn{j} and \eqn{j+1} components are calculated, as well as the LRTS. The original
LRTS is then compared to the \code{quantile} quantile of the bootstrapped counterparts; if it lies 
within this range, \eqn{j} is returned as the order estimate, otherwise \eqn{j} is increased by 1 
and the procedure is started over.

The MLEs are calculated via the \code{MLE.function} attribute for \eqn{j = 1}, if it is supplied. 
For all other \eqn{j} (and also for \eqn{j = 1} in case \code{MLE.function = NULL}) the solver 
\code{\link[Rsolnp]{solnp}} is used to calculate the minimum of the negative log likelihood. As 
initial values (for \code{\link[Rsolnp]{solnp}}), the data is clustered into \eqn{j} groups via 
\code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function}
(if supplied, otherwise numerical optimization is used here as well). The size of the groups is 
taken as initial component weights and the MLEs are taken as initial parameter estimates.
}
\value{
Object of class \code{paramEst} with the following attributes

\item{dat}{numeric of underlying data of \code{obj}.}

\item{dist}{character string stating the (abbreviated) name of the component distribution.}

\item{ndistparams}{integer specifying the number of parameters identifying the distribution       
  \code{dist}.}

\item{formals.dist}{string vector specifying the names of the formals identifying the distribution   \code{dist}.}

\item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}

\item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}

\item{pars}{Say \eqn{d = } \code{ndistparams}. Then \code{pars} is a numeric of parameter estimates of size \eqn{(d+1)*p-1}, given as 

\eqn{(w_1, ... w_{p-1}, \theta 1_1, ... \theta 1_p, \theta 2_1, ... \theta d_p)}.}

\item{values}{numeric of function values gone through during optimization, where the last one is the value at the optimum.}

\item{convergence}{indicates whether the solver has converged (0) or not (1 or 2).}
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
\code{\link[Rsolnp]{solnp}} for the solver, \code{\link{datMix}} for the creation of the \code{datMix} object.
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

normLoc.dM <- RtoDat(normLocRMix, param.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list)
                
                      
### using 'datMix' object to estimate the mixture

set.seed(0)
res <- mix.lrt(normLoc.dM, B = 30)
plot(res)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}% use one of  RShowDoc("KEYWORDS")

