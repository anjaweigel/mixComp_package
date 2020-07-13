\name{paramHankel}
\alias{paramHankel}
\alias{paramHankel.scaled}
\alias{print.paramEst}
\alias{plot.paramEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Mixture Complexity and Parameters Based on Hankel Matrix
}
\description{
Estimation of the complexity and the component parameters of an underlying mixture based on 
estimating the determinant of the Hankel matrix made up of the moments of the mixing distribution
and comparing it to determinant values generated via parametric bootstrap.
}
\usage{
paramHankel(obj, j.max = 10, B = 1000, ql = 0.025, 
            qu = 0.975, control = c(trace = 0), \dots)

paramHankel.scaled(obj, j.max = 10, B = 100, ql = 0.025, 
                   qu = 0.975, control = c(trace = 0), \dots)
                   
\method{plot}{paramEst}(x, mixture = TRUE, components = TRUE, ylim = NULL, 
     cex.main = 0.9, \dots)

\method{print}{paramEst}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{object of class \code{datMix}.}
  
  \item{j.max}{integer stating the maximal number of components to be considered.}  
  
  \item{B}{integer specifying the number of bootstrap replicates.}
  
  \item{ql}{numeric between \eqn{0} and \eqn{1} specifying the lower bootstrap quantile to which
    the observed value will be compared.}
    
  \item{qu}{numeric between \eqn{0} and \eqn{1} specifying the upper bootstrap quantile to which
    the observed value will be compared.}  
    
  \item{control}{control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.}  
  
  \item{x}{object of class \code{paramEst}.}
  
  \item{mixture}{logical indicating whether the estimated mixture density should
    be drawn, defaulting to \code{TRUE}.}

  \item{components}{logical indicating whether the individual mixture components should
    be drawn, defaulting to \code{TRUE}.}
  
  \item{ylim}{range of y values to use; if not specified (or
    containing \code{NA}), the function tries to construct reasonable 
    default values itself.}
    
  \item{cex.main}{The magnification to be used for main titles relative to the current setting 
      of \code{cex}, see \code{\link[graphics]{par}}.} 
    
  \item{\dots}{
    \describe{
      \item{in \code{paramHankel()} and \code{paramHankel.scaled()}:}{further arguments passed 
      to the \code{\link[boot]{boot}} 
      function.}
      \item{in \code{plot.hankDet()}:}{further arguments passed to the 
      \code{\link[graphics]{hist}} function plotting the data.}
    }}
  
}
\details{

Define the \eqn{order} or \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer 
\eqn{p}, such that its pdf/pmf \eqn{f} can be written as

\eqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p)}.

\code{paramHankel} estimates \eqn{p} by iteratively increasing the assumed order \eqn{j} 
and calculating the determinant of the \eqn{(j+1)}x\eqn{(j+1)} Hankel matrix made up of the 
first \eqn{2j} raw moments of the mixing distribution (for details see \code{\link{nonparamHankel}}).
Then, for a given \eqn{j}, the MLE for a \eqn{j} component mixture is calculated and \code{B}
parametric bootstrap samples of size \eqn{n} (size of the data) are generated from the distribution
corresponding to the MLE. For each of the \code{B} samples the determinant of the resulting 
\eqn{(j+1)}x\eqn{(j+1)} Hankel matrix is calculated, and the original determinant value is compared
to the bootstrap quantiles \code{ql} and \code{qu}. If the original determinant lies within this range,
\eqn{j} is returned as the order estimate, otherwise \eqn{j} is increased by 1 and the procedure is 
started over. 

\code{paramHankel.scaled} does the same as \code{paramHankel} with the exception that the determinants
are devided by their standard deviation. For the bootstrapped determinants, this denominator is simply
calculated as the empirical standard deviation of the bootstrap sample. For the original determinant, 
\code{B} nonparametric bootstrap samples of size \eqn{n} are generated from the data, the corresponding
determinants are calculated and their empirical standard deviation is used.

The MLEs are calculated via the \code{MLE.function} attribute for \eqn{j = 1}, if it is supplied. For 
all other \eqn{j} (and also for \eqn{j = 1} in case \code{MLE.function = NULL}) the solver 
\code{\link[Rsolnp]{solnp}} is used to calculate the minimum of the negative log likelihood. As initial
values (for \code{\link[Rsolnp]{solnp}}), the data is clustered into \eqn{j} groups via 
\code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function}
(if supplied, otherwise numerical optimization is used here as well). The size of the groups is taken 
as initial component weights and the MLEs corresponding to each group are taken as initial parameter estimates.


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

\item{pars}{Say \eqn{d = } \code{ndistparams}. Then \code{pars} is a numeric of parameter estimates of size
            \eqn{(d+1)*p-1}, given as 

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
\code{\link{nonparamHankel}} for estimation of the mixture complexity based on the Hankel
matrix without parameter estimation, \code{\link[Rsolnp]{solnp}} for the solver, \code{\link{datMix}}
for creation of the \code{datMix} object.
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
# mixing distribution via Hankel.method "natural"

psi.pois <- function(dat, j){
  res <- 1
  for (i in 0:(j-1)){
    res <- res*(dat-i)
  }
  res
}
        
# generating 'datMix' object
pois.dM <- RtoDat(poisRMix, param.bound.list = poisList, MLE.function = MLE.pois,
                  Hankel.method = "natural", Hankel.function = psi.pois)


## complexity and parameter estimation
set.seed(1)
res <- paramHankel(pois.dM)
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
