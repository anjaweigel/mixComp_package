\name{nonparamHankel}
\alias{nonparamHankel}
\alias{print.hankDet}
\alias{plot.hankDet}

\title{
Estimate Mixture Complexity Based on Hankel Matrix
}

\description{
Estimation of a mixture's complexity based on estimating the determinant of the Hankel matrix of the moments of the mixing distribution. The estimated determinants can be scaled and/or penalized.
}

\usage{
nonparamHankel(obj, j.max = 10, pen.function = NULL, scaled = FALSE, 
               B = 1000, ...)

\method{print}{hankDet}(x, ...)

\method{plot}{hankDet}(x, type = "b", xlab = "j", ylab = NULL, mar = NULL, 
     ylim = c(min(0, min(obj)), max(obj)), ...)
}

\arguments{
    \item{obj}{object of class \code{datMix}.}
    
    \item{j.max}{integer specifying the maximal number of components to be considered.}
    
    \item{pen.function}{a function with arguments \code{j} and \code{n} specifying the penalty added to the determinant value given sample size \eqn{n} and the currently assumed complexity \eqn{j}. If left empty no penalty will be added. If non-empty and \code{scaled} is \code{TRUE}, the penalty function will be added after the determinants are scaled.}
    
    \item{scaled}{logical specifying whether the vector of estimated determinants should be scaled.}
                  
    \item{B}{integer specifying the number of bootstrap replicates used for scaling of the determinants. Ignored if \code{scaled} is \code{FALSE}.}
      
  \item{x}{object of class \code{hankDet}.}    
      
  \item{type}{character denoting type of plot, see, e.g. \code{\link[graphics]{lines}}. Defaults to \code{"b"}.}
    
  \item{xlab,ylab}{labels for the x and y axis with defaults (the default for \code{ylab} is created within the function, if no value is supplied).}
  
  \item{mar}{numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot, see \code{\link[graphics]{par}}.}
    
  \item{ylim}{range of y values to use.}  
  
  \item{\dots}{
    \describe{
      \item{in \code{nonparamHankel()}:}{further arguments passed to the \code{\link[boot]{boot}} function if \code{scaled} is \code{TRUE}.}
      \item{in \code{plot.hankDet()}:}{further arguments passed to \code{\link[base]{plot}}.}
      \item{in \code{print.hankDet()}:}{further arguments passed to \code{\link[base]{print}}.}
    }}
}

\details{
Define the \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as

\deqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p).}

\code{nonparamHankel} estimates \eqn{p} by iteratively increasing the assumed complexity \eqn{j} and calculating the determinant of the \eqn{(j+1)}x\eqn{(j+1)} Hankel matrix made up of the first \eqn{2j} raw moments of the mixing distribution. As shown by Dacunha-Castelle & Gassiat (1997), once the correct complexity is reached (i.e. for all \eqn{j >= p}), this determinant is zero. 

This suggests an estimation procedure for \eqn{p} based on initially finding a consistent estimator of the moments of the mixing distribution and then choosing the estimator \eqn{estim_p} as the value of \eqn{j} which yields a sufficiently small value of the 
determinant. Since the estimated determinant is close to 0 for all \eqn{j >= p}, this could lead to choosing \eqn{estim_p} rather larger than the true value. The function therefore returns all estimated determinant values corresponding to complexities up to \code{j.max},
so that the user can pick the lowest \eqn{j} generating a sufficiently small determinant. In addition, the function allows the inclusion of a penalty term as a function of the sample size \code{n} and the currently assumed complexity \code{j} which will be added to the determinant value (by supplying \code{pen.function}), and/or scaling of the determinants (by setting \code{scaled  = TRUE}). For scaling, a nonparametric bootstrap is used to calculate the covariance of the estimated determinants, with \code{B} being the size of the bootstrap sample. The inverse of the square root of this covariance matrix (i.e. the matrix \eqn{S^(-1)} such that $A = SS$, where A is the covariance matrix) is then multiplied with the estimated determinant vector to get the scaled determinant vector.

For a thorough discussion of the methods that can be used for the estimation of the moments see the details section of \code{\link{datMix}}.
}

\value{
The vector of estimated determinants (optionally scaled and/or penalized), given back as an object of class \code{hankDet} with the following attributes:
\item{scaled}{logical indicating whether the determinants are scaled.}
\item{pen }{logical indicating whether a penalty was added to the determinants.}
\item{dist }{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}
}

\references{
D. Dacunha-Castelle and E. Gassiat, "The estimation of the order of a mixture model", Bernoulli, Volume 3, Number 3, 279-299, 1997.
}

\seealso{
\code{\link{paramHankel}} for a similar approach which estimates the component weights and parameters on top of the complexity, 
\code{\link{datMix}} for the creation of the \code{datMix} object.
}

\examples{
## create 'Mix' object
geomMix <- Mix("geom", w = c(0.1, 0.6, 0.3), prob = c(0.8, 0.2, 0.4))

## create random data based on 'Mix' object (gives back 'rMix' object)
set.seed(1)
geomRMix <- rMix(1000, obj = geomMix)

## create 'datMix' object for estimation

# explicit function giving the estimate for the j^th moment of the 
# mixing distribution, needed for Hankel.method "explicit"

explicit.fct.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}
        
## generating 'datMix' object
geom.dM <- RtoDat(geomRMix, Hankel.method = "explicit", 
                  Hankel.function = explicit.fct.geom)

## function for penalization
pen <- function(j, n){
  (j*log(n))/(sqrt(n))
}

## estimate determinants
set.seed(1)
geomdets_pen <- nonparamHankel(geom.dM, pen.function = pen, j.max = 5)
plot(geomdets_pen, main = "Three component geometric mixture")
}

\keyword{cluster}
