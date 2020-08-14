\name{L2.disc}
\alias{L2.disc}
\alias{L2.boot.disc}

\title{
Estimate a Discrete Mixture's Complexity Based on L2 Distance
}

\description{
Estimation of a discrete mixture's complexity as well as its component weights and parameters by minimizing the squared L2 distance to the empirical probability mass function. 
}

\usage{
L2.disc(obj, j.max = 10, n.inf = 1000, threshold = "SBC", control = c(trace = 0))

L2.boot.disc(obj, j.max = 10, n.inf = 1000, B = 100, 
             ql = 0.025, qu = 0.975, control = c(trace = 0), ...)

}

\arguments{
  \item{obj}{object of class \code{datMix}.}
    
  \item{j.max}{integer stating the maximal number of components to be considered.}
  
  \item{n.inf}{integer; the L2 distance contains an infinite sum, which will be approximated by a sum ranging from 0 to \code{n.inf}.}
    
  \item{threshold}{function or character string in \code{c("LIC", "SBC")} specifying which threshold should be used to compare two mixture estimates of complexities \eqn{j} and \eqn{j+1}. If the difference in minimized squared distances is smaller than the relevant threshold, \eqn{j} will be returned as complexity estimate.}
    
  \item{B}{integer specifying the number of bootstrap replicates.}
    
  \item{ql}{numeric between \eqn{0} and \eqn{1} specifying the lower quantile to which the observed difference in minimized squared distances will be compared.}
    
  \item{qu}{numeric between \eqn{0} and \eqn{1} specifying the upper quantile to which the observed difference in minimized squared distances will be compared.}  
  
  \item{control}{control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.}
    
  \item{\dots}{further arguments passed to the \code{\link{boot}} function.}
}

\details{

Define the \eqn{complexity} of a finite discrete mixture \eqn{F} as the smallest integer \eqn{p}, such that its probability mass function (pmf) \eqn{f} can be written as

\eqn{f(x) = w_1*g(x;\theta_1) + \dots + w_p*g(x;\theta_p)}.

Further, let \eqn{g, f} be two probability mass functions. The squared L2 distance between \eqn{g} and \eqn{f} is given by

\eqn{L_2^2(g,f) = \sum(g(x)-f(x))^2}.

To estimate \eqn{p}, \code{L2.disc} iteratively increases the assumed complexity \eqn{j} and finds the ``best'' estimate for both, the pmf of a mixture with \eqn{j} and \eqn{j+1} components, by calculating the parameters that minimize the squared L2 distances to the empirical probability mass function. The infinite sum contained in the objective function will be approximated by a sum ranging from 0 to \code{n.inf}, set to 1000 by default. Once the ``best'' parameters have been obtained, the difference in squared distances is compared to a predefined \code{threshold}. If this difference is smaller than the threshold, the algorithm terminates and the true complexity is estimated as \eqn{j}, otherwise \eqn{j} is increased by 1 and the procedure is started over. The predefined thresholds are the \code{"LIC"} given by \eqn{(0.6*log((j+1)/j))/n} and the \code{"SBC"} given by \eqn{(0.6*log(n)*log((j+1)/j))/n}, \eqn{n} being the sample size. Note that, if a customized function is to be used, it may only take the arguments \code{j} and \code{n}.

\code{L2.boot.disc} works similarly to \code{L2.disc} with the exception that the difference in squared distances is not compared to a predefined threshold but a value generated by a bootstrap procedure. At every iteration (of \eqn{j}), the function sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}, using a parametric bootstrap to generate \code{B} samples of size \eqn{n} from a \eqn{j}-component mixture given the previously calculated ``best'' parameter values. For each of the bootstrap samples, again the ``best'' estimates corresponding to pmfs with \eqn{j} and \eqn{j+1} components are calculated, as well as their difference in squared L2 distances from the empirical probability mass function. The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the original difference in squared distances lies outside of the interval \eqn{[ql, qu]}, specified by the \code{ql} and \code{qu} empirical quantiles of the bootstrapped differences. Otherwise, \eqn{j} is returned as the complexity estimate. 

To calculate the minimum of the L2 distance (and the corresponding parameter values), the solver \code{\link[Rsolnp]{solnp}} is used. The initial values supplied to the solver are calculated as follows: the data is clustered into \eqn{j} groups by the function \code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function} (if supplied to the \code{datMix} object \code{obj}, otherwise numerical optimization is used here as well). The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates.
}

\value{
Object of class \code{paramEst} with the following attributes:

\item{dat}{data based on which the complexity is estimated.}

\item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}

\item{ndistparams}{integer specifying the number of parameters identifying the component distribution, i.e. if \eqn{\theta \subseteq R^d} then \code{ndistparams}\eqn{ = d}.}

\item{formals.dist}{string vector specifying the names of the formal arguments identifying the distribution \code{dist} and used in \code{ddist} and \code{rdist}, e.g. for a gaussian mixture (\code{dist = norm}) amounts to \code{mean} and \code{sd}, as these are the formal arguments used by \code{dnorm} and \code{rnorm}.}

\item{discrete}{logical indicating whether the underlying mixture distribution is discrete. Will always be \code{TRUE} in this case.}

\item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}

\item{pars}{Say the complexity estimate is equal to some \eqn{j}. Then \code{pars} is a numeric vector of size \eqn{(d+1)*j-1} specifying the component weight and parameter estimates, given as 

\eqn{(w_1, ... w_{j-1}, \theta 1_1, ... \theta 1_j, \theta 2_1, ... \theta d_j)}.}

\item{values}{numeric vector of function values gone through during optimization at iteration \eqn{j}, the last entry being the value at the optimum.}

\item{convergence}{indicates whether the solver has converged (0) or not (1 or 2) at iteration \eqn{j}.}
}

\references{
T. Umashanger and T. Sriram, "L2E estimation of mixture complexity for count data", Computational Statistics and Data Analysis 51, 4379-4392, 2007.
}

\seealso{
\code{\link{hellinger.disc}} for the same estimation method using the Hellinger distance, \code{\link[Rsolnp]{solnp}} for the solver,
\code{\link{datMix}} for the creation of the \code{datMix} object.
}

\examples{
## create 'Mix' object
poisMix <- Mix("pois", w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))


## create random data based on 'Mix' object (gives back 'rMix' object)
set.seed(1)
poisRMix <- rMix(10000, obj = poisMix)


## create 'datMix' object for estimation

# generate list of parameter bounds
poisList <- vector(mode = "list", length = 1)
names(poisList) <- "lambda"
poisList$lambda <- c(0, Inf)

# generate MLE function
MLE.pois <- function(dat){
  mean(dat)
}

# generating 'datMix' object
pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois)


## complexity and parameter estimation 
set.seed(1)
res <- L2.disc(pois.dM)
plot(res)
}

\keyword{cluster}
