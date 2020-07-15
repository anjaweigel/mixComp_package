

## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "explicit" 

.deterHankel.explicit <- function(dat, inds = 1:length(dat), Hankel.function, j.max = NULL,
                                  j = NULL){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel
  
  dat <- dat[inds] # for bootstrap
  
  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max) 
    mn <- 1
    mx <- j.max
  } 
  
  Hankel.functionV <- Vectorize(Hankel.function, vectorize.args = "j")
  cp_moments <- Hankel.functionV(dat, 1:(2*mx))
  cp_add1 <- c(1, cp_moments)
  
  for(i in mn:mx){    
    
    H <- hankel.matrix((i + 1), cp_add1[1:((2*i)+1)])
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants
    
  }
  
  return(det.vec)
}


## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "translation" 

.deterHankel.translation <- function(dat, inds = 1:length(dat), Hankel.function, j.max = NULL,
                                     j = NULL){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel
  
  dat <- dat[inds] # for bootstrap
  n <- length(dat)
  
  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max)
    mn <- 1
    mx <- j.max
  } 
  
  # construct elements of triagular linear system
  # empirical moments of mixture distribution
  X_hat <- (1/n) * mapply(function(x, y){sum(x^y)}, rep(list(dat), 2*mx), 1:(2*mx))
  Hankel.functionV <- Vectorize(Hankel.function)
  # theoretical moments of G
  EY <- Hankel.functionV(1:(2*mx))
  
  b <- X_hat - EY # vector b
  EY_A <- c(1, EY)

  ## generate matrix A
  A <- matrix(0, nrow = 2*mx, ncol = 2*mx)
  for (i_row in 1:(2*mx)){
    for (i_col in 1:i_row){
      A[i_row, i_col] <- choose(i_row, i_row-i_col)*EY_A[i_row-i_col+1]  
    }
  }   
  
  # solve system
  cp_moments <- solve(a = A, b = b)
  
  cp_add1 <- c(1, cp_moments)
  for(i in mn:mx){
    
    H <- hankel.matrix((i + 1), cp_add1)
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants
    
  }
  
  return(det.vec)
}



## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "scale" 

.deterHankel.scale <- function(dat, Hankel.function, inds = 1:length(dat), j.max = NULL, j = NULL,  
                               message = FALSE){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel
  
  dat <- dat[inds] # for bootstrap
  n <- length(dat)
  
  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max)
    mn <- 1
    mx <- j.max
  } 
  
  # construct elements of triagular linear system
  # empirical moments of mixture distribution
  X_hat <- (1/n) * mapply(function(x, y){sum(x^y)}, rep(list(dat), 2*mx), 1:(2*mx))
  Hankel.functionV <- Vectorize(Hankel.function)
  # theoretical moments of G
  EY <- Hankel.functionV(1:(2*mx))
  
  if (sum(EY == 0) != 0){  # cannot devide by zero, in that case take squares of everything
    
    if(message == TRUE) message("Moments of mixing distribution M^k are 0 for some k. Calculating M^(2k) instead.")
    
    X_hat <- (1/n) * mapply(function(x,y){sum(x^y)}, rep(list(dat),2*mx), 
                            seq(from = 2, to = (4*mx), by = 2))
    EY <- Hankel.functionV(seq(from = 2, to = (4*mx), by = 2))
    
  }  
  
  cp_moments <- X_hat/EY
  cp_add1 <- c(1, cp_moments)
  
  for(i in mn:mx){  
    
    H <- hankel.matrix((i + 1), cp_add1)
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants
    
  }
  
  return(det.vec)
}


## Purpose: return the correct function for calculating the determinant based on
##          "Hankel.method"

.moments.map <- function(Hankel.method) {
  
  if (Hankel.method == "explicit"){
    fun <- .deterHankel.explicit
  }
  else if (Hankel.method == "translation"){
    fun <- .deterHankel.translation
  }
  else if (Hankel.method == "scale"){
    fun <- .deterHankel.scale
  }
  
  fun
}



## Purpose: check whether the inputs to any of the functions estimating mixtuer complexity
##          are of the correct form

.input.checks.functions <- function(obj, B, j.max, ql, qu, quantile, n.inf, thrshL2, 
                                    thrshHel, discrete, continuous, pen.function, scaled, Hankel,
                                    param, bandwidth, sample.n, sample.plot){
 
  if(Hankel == TRUE){
    # if the function uses the Hankel matrix of the moments of the mixing distribution
    
    if(is.null(attr(obj, "Hankel.method"))) 
      stop("'Hankel.method' has to be defined as a datMix object attribute!")
    if(!attr(obj, "Hankel.method") %in% c("explicit", "translation", "scale")){
      stop("Hankel.method is not one of the implemented methods. 
           Please enter \"explicit\", \"translation\" or \"scale\"")
    }

    if(is.null(attr(obj, "Hankel.function"))) 
      stop("'Hankel.function' has to be defined as a datMix object attribute!")
    if(!is.function(attr(obj, "Hankel.function")))
      stop("'Hankel.function' has to be a function!")
  }
  
  if(param == TRUE){ 
    # if the functions estimates the weights as well as the component parameters 
    if(is.null(attr(obj, "theta.bound.list")))
      stop("'theta.bound.list' has to be specified as a datMix object attribute!")
  }
  
  if(is.null(obj) || !is.datMix(obj)) stop("'obj' has to be a 'datMix' object!")
  
  if(is.null(j.max) || !is.numeric(j.max) || !(as.integer(j.max) == j.max)) 
    stop("'j.max' has to be an integer!")

  if(missing(pen.function)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(!is.null(pen.function) && !is.function(pen.function))
    stop("'pen.function' has to be NULL or a function!")
  else if(!is.null(pen.function) && (!all(names(formals(pen.function)) %in% c("n", "j"))
                                     || !all(c("n", "j") %in% names(formals(pen.function)))))
    stop("'pen.function' must contain the arguments \"n\" and \"j\".")
  
  if(missing(scaled)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(scaled) || !(is.logical(scaled))) stop("'scaled' has to be logical!")
  
  if(missing(B)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if (is.null(B) || !is.numeric(B) || !(as.integer(B) == B)) stop("'B' has to be an integer!")
  
  if(missing(ql) & missing(qu)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(ql) ||  is.null(qu) || !is.numeric(ql) || !is.numeric(qu) ||
          !(ql >= 0 & ql <= 1) || !(qu >= 0 & qu <= 1) || ql > qu)
    stop("'ql' and 'qu' have to be numerics between 0 and 1 with qu > ql!")
  
  if(missing(quantile)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(!is.numeric(quantile) || !(quantile >= 0 & quantile <= 1)) 
    stop("'quantile' has to be a numeric between 0 and 1!")
  
  if(missing(n.inf)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(n.inf) || !is.numeric(n.inf) || !(as.integer(n.inf) == n.inf)) 
    stop("'n.inf' has to be an integer!")
  
  if(missing(thrshL2)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(thrshL2) || !(is.function(thrshL2)|| thrshL2 %in% c("LIC", "SBC"))) 
    stop("'threshold' has to be a user-entered function or an element of c(\"LIC\", \"SBC\")")
  else if(is.function(thrshL2) && (!all(names(formals(thrshL2)) %in% c("n", "j"))
                                   || !all(c("n", "j") %in% names(formals(thrshL2)))))
    stop("The function 'threshold' needs to have arguments \"n\" and \"j\".")
  else if(!is.function(thrshL2) && thrshL2 == "LIC") 
    warning("While being used in Umashanger & Sriram's original paper,  asymptotically, the threshold 'LIC' does not go to 0 slower than the difference in squared L2 distances once the correct order p is reached and the estimator is therefore not consistent.",
            call. = FALSE)
  
  if(missing(thrshHel)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(thrshHel) || !(is.function(thrshHel)|| thrshHel %in% c("AIC", "SBC"))) 
    stop("'threshold' has to be a user-entered function or an element of c(\"AIC\", \"SBC\")")
  else if(is.function(thrshHel) && (!all(names(formals(thrshHel)) %in% c("n", "j"))
                                    || !all(c("n", "j") %in% names(formals(thrshHel)))))
    stop("The function 'threshold' needs to have arguments \"n\" and \"j\".")
  else if(!is.function(thrshHel) && thrshHel == "AIC") 
    warning("While being used in Woo & Sriram's original paper, asymptotically, the threshold 'AIC' does not go to 0 slower than the difference in squared Hellinger distances once the correct order p is reached and the estimator is therefore not consistent.",
            call. = FALSE)

  if(missing(discrete)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(discrete) || discrete != TRUE) 
    stop("The 'discrete' attribute of the Rdat object is not set to TRUE, however
         this function only works for discrete data.")
  
  if(missing(continuous)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(continuous) || continuous != TRUE) 
    stop("The 'discrete' attribute of the Rdat object is set to TRUE, however
         this function only works for continuous data.")
 
  if(missing(bandwidth)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(bandwidth) || 
          (bandwidth != "adaptive" && (!is.numeric(bandwidth) || length(bandwidth) != 1)))
    stop("'bandwidth' has to be either \"adaptive\" or a single numeric!")
  
  if(missing(sample.n)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(sample.n) || !is.numeric(sample.n) || !(as.integer(sample.n) == sample.n)) 
    stop("'sample.n' has to be an integer!")
  
  if(missing(sample.plot)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(sample.plot) || !(is.logical(sample.plot))) 
    stop("'sample.plot' has to be logical!")
  
} 


## Purpose: calculate the scaled vector of determinants of the Hankel matrix of the 
##          moments of the mixing distribution

.deterHankel.scaled <- function(dat, Hankel.method, Hankel.function, j.max = 10, B = 1000, ...){
  
  fun <- .moments.map(Hankel.method = Hankel.method)
  
  if(Hankel.method == "scale") # get warning if squared moments have to be calculated
    D_hat <- fun(dat = dat, Hankel.function = Hankel.function, j.max = j.max, message = TRUE)
  else D_hat <- fun(dat = dat, Hankel.function = Hankel.function, j.max = j.max) 
  
  # bootstrapped vector of determinants
  bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
             j.max = j.max, ...) 
  D_boot <- bt$t
  
  cov_boot <- cov(D_boot)
  cov_boot_sqrt <- sqrtm(cov_boot)
  cov_boot_sqrt_inv <- try(solve(cov_boot_sqrt), silent = TRUE)
  if(inherits(cov_boot_sqrt_inv, "try-error")) {
    stop("The matrix squareroot of the estimated covariance matrix could not be inverted: singular system.")
  }
  Y_boot <- cov_boot_sqrt_inv %*% D_hat  
  
  det.vec <- abs(Y_boot)
  return(as.vector(det.vec))
}


## Purpose: function returning the vector of estimated determinants for orders up to
##          j.max; returns an object of class "hankDet"

nonparamHankel <- function(obj, j.max = 10, pen.function = NULL, scaled = FALSE, B = 1000,
                           ...){

  .input.checks.functions(obj, j.max = j.max, pen.function = pen.function, scaled = scaled,
                          B = B, Hankel = TRUE, param = FALSE)
  
  Hankel.method <- attr(obj, "Hankel.method")
  Hankel.function <- attr(obj, "Hankel.function")
  
  dat <- as.numeric(obj)
  
  if(!is.null(pen.function)){
    pen <- TRUE
    n <- length(dat)
  }
  else pen <- FALSE
  
  if (scaled == TRUE){
    res <- .deterHankel.scaled(dat, Hankel.function = Hankel.function, 
                               Hankel.method = Hankel.method, j.max = j.max, B = B, ...)
    
  } else {
    
    fun <- .moments.map(Hankel.method = Hankel.method)
    res <- fun(dat = dat, j.max = j.max, Hankel.function = Hankel.function) 
  }  
  
  if (pen == TRUE) res <- res + pen.function(j = 1:j.max, n = n)
  
  class(res) <- "hankDet"
  attr(res, "scaled") <- scaled
  attr(res, "pen") <- pen
  attr(res, "dist") <- attr(obj, "dist")
  return(res)
}


## Purpose: print method for "hankDet" objects

print.hankDet <- function(x){
  
  obj <- x
  scaled <- attr(obj, "scaled")
  pen <- attr(obj, "pen")
  dist <- attr(obj, "dist")
  
  if(scaled == TRUE & pen == TRUE){
    header <- paste("\nEstimation of the scaled and penalized determinants for a '", dist, 
                    "' mixture model:\n", "\n", sep = "")
  } else if(scaled == TRUE & pen == FALSE){
    header <- paste("\nEstimation of the scaled determinants for a '", dist, 
                    "' mixture model:\n", "\n", sep = "")
  } else if(scaled == FALSE & pen == TRUE){
    header <- paste("\nEstimation of the penalized determinants for a '", dist, 
                    "' mixture model:\n", "\n", sep = "")
  } else {
    header <- paste("\nEstimation of the determinants for a '", dist, "' mixture model:\n",
                    "\n", sep = "")
  } 
  
  cat(header)
  cmat <- matrix(c(1:length(obj), obj), nrow = length(obj), ncol = 2)
  colnames(cmat) <- c("Number of components", "Determinant")
  rownames(cmat) <- rep("", length(obj))
  print(cmat)
  
}



## Purpose: plot method for "hankDet" objects

plot.hankDet <- function(x, type = "b", xlab = "Complexity Estimate",  ylab = NULL, mar = NULL,
                         ylim = c(min(0, min(obj)), max(obj)), ...){
  
  obj <- x
  scaled <- attr(obj, "scaled")
  pen <- attr(obj, "pen")
  
  if(is.null(ylab)){
    if(scaled == TRUE & pen == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"), 
                                   textstyle("(scaled and penalized)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    } 
    else if (scaled == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"), 
                                   textstyle("(scaled)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    } 
    else if (pen == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"), 
                                   textstyle("(penalized)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    } else {
      ylab <- "Determinant value"
      if(is.null(mar)) mar <- c(5,5,4,2) + 0.1 # ylab needs less space on the left
    } 
  }
  
  par(mar = mar)
  plot(x = 1:length(obj), y = obj, type = type, ylab = ylab, xlab = xlab, ylim = ylim, ...)
  
}
