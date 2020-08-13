
## Purpose: returns the negative log likelihood function of parameters x
##          (length of x determines the mixture complexity)

.get.negloglik.dist <- function(dat, dist, formals.dist, ndistparams, dist_call){
  
  function(x){
   
    j <- (length(x) + 1)/(ndistparams + 1) # current complexity estimate
    n <- length(dat)
    w <- x[1:(j-1)]

    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j-1)], nrow = n, ncol = j, byrow = TRUE)
    }
    theta.list.long$x <- dat
    
    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j-1)]))
    if(any(w < 0)) return(0)
    res <- as.matrix(mat) %*% w
    if(any(is.na(res))) return(0)
    
    # avoid log(0)
    res[res == 0] <- ifelse(suppressWarnings(min(res[res != 0])) < .Machine$double.xmin, 
                            min(res[res != 0]), .Machine$double.xmin)
    
    return(-sum(log(res)))
  }
}



## Purpose: returns the negative log likelihood function of parameters x
##          when the mixture consists of only a single component

.get.negloglik.dist.0 <- function(dat, dist, formals.dist, ndistparams, dist_call){
  
  function(x){
  
    n <- length(dat)
    dist_call <- get(paste("d", dist, sep = ""))
    
    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist  
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- rep(x[i], n)
    }
    theta.list.long$x <- dat
    
    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    res <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    if(any(is.na(res))) return(0)
    
    # avoid log(0)
    res[res == 0] <- ifelse(suppressWarnings(min(res[res != 0])) < .Machine$double.xmin, 
                            min(res[res != 0]), .Machine$double.xmin)
    
    return(-sum(log(res)))
  }
}


## Purpose: returns the MLE function contained in a 'datMix' object as a list

.get.MLE.function <- function(obj){
  
  if(!is.null(attr(obj, "MLE.function"))){
    if(is.list(attr(obj, "MLE.function"))) MLE.function <- attr(obj, "MLE.function")
    else MLE.function <- list(attr(obj, "MLE.function"))
  } 
  else MLE.function <- NULL
  return(MLE.function)
  
}


## Purpose: returns the restrictions set when optimizing the negative log likelihood
##          function

.get.restrictions <- function(j, ndistparams, lower, upper){
  
  # first j-1 weights need to sum up to less than 1
  ineq <- function(x){
    return(sum(x[1:(j - 1)]))
  }
  
  # lower bound for all parameters
  lx <- rep(0, j - 1) # lower bound for weights
  for(i in 1:ndistparams){lx <- c(lx, rep(lower[i], j))} # lower bound for component parameters
  
  # upper bound for all parameters
  ux <- rep(1, j - 1) # upper bound for weights
  for(i in 1:ndistparams){ux <- c(ux, rep(upper[i], j))} # upper bound for component parameters
  
  return(list("ineq" = ineq, "lx" = lx, "ux" = ux))
}  


## Purpose: returns the parameters (weights and component parameters)
##          used in parametric bootstrap

.get.bootstrapparams <- function(formals.dist, ndistparams, mle.est, j)  {
  
  theta.list <- vector(mode = "list", length = ndistparams)
  names(theta.list) <- formals.dist
  for (i in 1:ndistparams){
    theta.list[[i]] <-  mle.est[(i*j):((1 + i)*j - 1)]
  }
  
  if(length(mle.est) != ndistparams){
    w <- c(mle.est[1:(j - 1)], 1 - sum(mle.est[1:(j - 1)]))
  } else w <- 1
  
  return(list("theta.list" = theta.list, "w" = w))  
}



## Purpose: returns initial values for the calculation of the initial values
##          passed to solnp (in case numerical optimization has to be used
##          to determine the initial values)

.get.initialvals.init <- function(j, ndistparams, lower, upper){
  
  lower[lower == -Inf] <- -100
  upper[upper == Inf] <- 100
  lower.new <- lower + (upper-lower)/4 # such that we are in the midde 50% of the range
  upper.new <- upper - (upper-lower)/4 # such that we are in the midde 50% of the range
  
  if(j != 1) initial <- rep(1/j, (j - 1)) # uniform weights
  else initial <- NULL
  
  initial <- c(initial, runif(j*ndistparams, rep(lower.new, each = j), 
                              rep(upper.new, each = j)))
  initial
}



## Purpose: returns initial values to be passed to solnp

.get.initialvals <- function(dat, j, ndistparams, MLE.function = NULL, lower, upper, dist,
                             formals.dist, dist_call){
  
  # cluster data into j groups
  cluster.vec <- clara(dat, j, rngR = TRUE, pamLike = TRUE, medoids.x = FALSE)$clustering
  groups <- mapply(function(i){sum(abs(cluster.vec - i) <= .Machine$double.eps)}, 1:j)
  
  # check that each group has at least 2 observations
  if(any(groups == 1)){
    ind <- which(groups == 1)
    for(i in ind){
      distance <- abs(dat[cluster.vec == i] - dat)
      distance[distance == 0] <- max(distance)
      closest <- which(distance == min(distance))
      cluster.vec[closest] <- i
    }
  }
  
  # initial weights proportional to number of observations in each cluster
  initial <- (mapply(function(i){sum(abs(cluster.vec - i) <= .Machine$double.eps)}, 1:j)/length(cluster.vec))[-j]
  
  if(!is.null(MLE.function)){
    
    # apply MLE function to each of the clusters
    init.mat <- mapply(function(i){sapply(MLE.function, 
                                          function(fun) fun(dat[cluster.vec == i]))}, 1:j)
    
  } else {
    
    # use numerical optimization to find the MLE for each of the clusters
    init <- .get.initialvals.init(j = 1, ndistparams = ndistparams, lower = lower,
                                  upper = upper)
    likelihood0list <- mapply(function(i){.get.negloglik.dist.0(dat[cluster.vec == i], 
                                                                dist, formals.dist, ndistparams, 
                                                                dist_call)}, 1:j)
    init.mat <- mapply(function(i){solnp(init, likelihood0list[[i]], LB = lower, UB = upper, 
                                         control= c(trace = 0))$pars}, 1:j)
    
  }
  
  initial.theta <- as.vector(t(init.mat))
  
  # don't want starting values on the boarder of the feasible region
  initial.theta[initial.theta == lower] <- 1.05*initial.theta[initial.theta == lower]
  initial.theta[initial.theta == lower] <- 0.05 # if lower == 0
  initial.theta[initial.theta == upper] <- 0.95*initial.theta[initial.theta == upper]
  initial.theta[initial.theta == upper] <- -0.05 # if upper == 0
  
  initial <- c(initial, initial.theta)
  initial
}



## Purpose: brings weights back into feasible region in case solnp goes outside of it
##          (this may happen for convergence and conditioning reasons)

.augment.pars <- function(pars, j){
  
  if(any(pars[1:(j-1)] < 0)) pars[1:(j - 1)][pars[1:(j - 1)] < 0] <- 1e-08
  if(sum(pars[1:(j-1)]) >= 1) pars[1:(j - 1)] <- pars[1:(j - 1)]/(sum(pars[1:(j - 1)]) + 1e-08)
  pars
  
}



## Purpose: returns an object of class 'paramEst', which contains the result
##          of the function estimating the mixture complexity (alongside the
##          weights and component parameters)

.return.paramEst <- function(j, j.max, dat, pars, values, conv, dist, ndistparams, formals.dist,
                             discrete, MLE.function){
  
  if(j < j.max) cat(paste("\nThe estimated order is ", j, ".", sep = ""))
  else cat(paste("\nThe estimation procedure reached the value 'j.max' equal to ", j.max,
                 ".\nThe result showcases the last parameter estimates.", sep = ""))
  
  paramEst <- j
  attr(paramEst, "dat") <- as.vector(dat)
  attr(paramEst, "pars") <- pars
  attr(paramEst, "values") <- values
  attr(paramEst, "convergence") <- conv
  attr(paramEst, "dist") <- dist
  attr(paramEst, "ndistparams") <- ndistparams
  attr(paramEst, "formals.dist") <- formals.dist
  attr(paramEst, "discrete") <- discrete
  attr(paramEst, "mle.fct") <- MLE.function
  class(paramEst) <- "paramEst"
  invisible(paramEst)
  
}


## Purpose: prints the results of the numerical optimization procedure

.printresults <- function(res, j, dist, formals.dist, ndistparams){
  
  niter <- length(res$values)
  funv <- res$values[niter] # function value
  conv <- res$convergence
  if(conv == 0) conv <- TRUE
  else conv <- FALSE
  
  pars <- res$pars
  pars_complete <- numeric(j*(ndistparams + 1))
  if(j != 1){ 
    # add last weight completely determined by the others as they sum to 1
    pars_complete[1:(j - 1)] <- pars[1:(j - 1)]
    pars_complete[j] <- 1-sum(pars[1:(j - 1)])
    pars_complete[(j + 1):length(pars_complete)] <- pars[j:length(pars)]
  } else {
    pars_complete[1] <- 1
    pars_complete[2:length(pars_complete)] <- pars[1:length(pars)]
  }
  
  cmat <- matrix(pars_complete, nrow = j, ncol = ndistparams + 1)
  colnames(cmat) <- c("w", formals.dist)
  rownames(cmat) <- paste("Component ", 1:j, ":", sep = "")
  
  cat(paste("\nParameter estimation for a ", j, " component '", dist, "' mixture model:",
            "\nFunction value: ",
            format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
            "\n", "\n", sep = ""))
  printCoefmat(cmat)
  
  if(conv == TRUE){
    cat(paste( "Converged in ", niter, " iterations.\n", strrep("-", 70), sep = ""))
  } 
  else cat("Optimization algorithm did not converge.\n", strrep("-", 70), sep = "")
}



## Purpose: prints the results obtained by applying the MLE function

.printresultsMLE <- function(res, dist, formals.dist, ndistparams, likelihood0){
  
  pars <- c(1, res) # parameters
  funv <- likelihood0(res) # function value
  
  cmat <- matrix(pars, nrow = 1, ncol = ndistparams + 1)
  colnames(cmat) <- c("w", formals.dist)
  rownames(cmat) <- paste("Component ", 1, ":", sep = "")
  
  cat(paste("\nParameter estimation for a 1 component '", dist, "' mixture model:",
            "\nFunction value: ",
            format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
            "\n", "\n", sep = ""))
  printCoefmat(cmat)
  cat("Optimization via user entered MLE-function.\n", strrep("-", 70), sep = "")

}


## Purpose: returs a list of commonly needed variables; called at the beginning of
##          most function estimating the mixture complexity

.get.list <- function(obj){

  theta.bound.list <- attr(obj, "theta.bound.list")
  formals.dist <- names(theta.bound.list)
  ndistparams <- length(formals.dist)
  bounds <- as.vector(matrix(unlist(theta.bound.list, use.names = FALSE), 
                             nrow = ndistparams, byrow = TRUE))
  dat <- as.numeric(obj)
  dist <-  attr(obj, "dist")
  dist_call <- get(paste("d", dist, sep = ""))
  
  
  return(list(
    
    dist = dist,
    theta.bound.list = theta.bound.list,
    formals.dist = formals.dist,
    ndistparams = ndistparams,
    dist_call = dist_call,
    
    bounds = bounds,
    lower = bounds[1:ndistparams],
    upper = bounds[(ndistparams + 1):(2*ndistparams)],
    
    dat = dat,
    N = length(dat),
    n.max = max(dat),
    discrete = attr(obj, "discrete"),
    continuous = !attr(obj, "discrete"),
    
    Hankel.method = attr(obj, "Hankel.method"),
    Hankel.function = attr(obj, "Hankel.function"),
    MLE.function = .get.MLE.function(obj)))

}



## Purpose: method of estimating the mixture complexity (as well as the weights
##          and component parameters) returning a 'paramEst' object

paramHankel <- function(obj, j.max = 10, B = 1000, ql = 0.025,  qu = 0.975, 
                        control = c(trace = 0), ...){
  
  # check relevant inputs
  .input.checks.functions(obj, B = B, j.max = j.max, ql = ql, qu = qu,
                          Hankel = TRUE, param = TRUE)
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())
  
  fun <- .moments.map(Hankel.method = Hankel.method)
  likelihood <- .get.negloglik.dist(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  j <- 0 
  
  repeat{
    
    j <- j + 1 # current complexity estimate
    
    if (j == 1 && !is.null(MLE.function)){ # Calculate MLE via the MLE function
      # only possible for j equal to 1 and if a MLE function was supplied
      
      mle.est <- sapply(MLE.function, function(fun) fun(dat))
      .printresultsMLE(mle.est, dist, formals.dist, ndistparams, likelihood0)
      values <- likelihood0(mle.est)
      conv <- NULL
   
    } else { # Calculate MLE of a j component mixture numerically
      
      restrictions <- .get.restrictions(j = j, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      ineq <- restrictions$ineq
      lx <- restrictions$lx
      ux <- restrictions$ux
  
      initial <- .get.initialvals(dat, j, ndistparams, MLE.function, lower, upper, 
                                  dist, formals.dist, dist_call)
      
      if(j != 1){ # need to include weight restrictions in optimization
        
        opt <- solnp(initial, fun = likelihood, ineqfun = ineq, ineqLB = 0, ineqUB = 1,
                     LB = lx, UB = ux, control = control)
        opt$pars <- .augment.pars(opt$pars, j)
        opt$values[length(opt$values)] <- likelihood(opt$pars)
        
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)
        
      } else { # already know w = 1 (single component mixture)
        
        opt <- solnp(initial, fun = likelihood0, LB = lx, UB = ux, control = control)
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)
        
      }
    }

    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams, 
                                            mle.est = mle.est, j = j)
    Mix.mle <- Mix(dist = dist, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                   name = "Mix-boot")
    
    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }
    
    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, sim = "parametric", ran.gen = ran.gen, mle = Mix.mle, ...)
    det0 <- bt$t0 # original determinant value
    det.boot <- bt$t # bootstrapped determinant values
   
    q_lower <- quantile(det.boot, probs = ql)
    q_upper <- quantile(det.boot, probs = qu)
    
    if(det0 >= q_lower && det0 <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j == j.max){
      break
    }
  }
  
  .return.paramEst(j, j.max, dat, mle.est, values, conv, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}


## Purpose: method of estimating the mixture complexity (as well as the weights
##          and component parameters) returning a 'paramEst' object

paramHankel.scaled <- function(obj, j.max = 10, B = 100, ql = 0.025, qu = 0.975,
                               control = c(trace = 0), ...){
  
  # check relevant inputs
  .input.checks.functions(obj, B = B, j.max = j.max, ql = ql, qu = qu, Hankel = TRUE, 
                          param = TRUE)
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  fun <- .moments.map(Hankel.method = Hankel.method)
  likelihood <- .get.negloglik.dist(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  j <- 0 
  
  repeat{
    
    j <- j + 1 # current complexity estimate
    
    if (j == 1 && !is.null(MLE.function)){ # Calculate MLE via the MLE function
      # only possible for j equal to 1 and if a MLE function was supplied
      
      mle.est <- sapply(MLE.function, function(fun) fun(dat))
      .printresultsMLE(mle.est, dist, formals.dist, ndistparams, likelihood0)
      values <- likelihood0(mle.est)
      conv <- NULL
      
    } else { # Calculate MLE of a j component mixture numerically
      
      restrictions <- .get.restrictions(j = j, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      ineq <- restrictions$ineq
      lx <- restrictions$lx
      ux <- restrictions$ux
      initial <- .get.initialvals(dat, j, ndistparams, MLE.function, lower, upper, 
                                  dist, formals.dist, dist_call)
      
      if(j != 1){ # need to include weight restrictions in optimization
        
        opt <- solnp(initial, fun = likelihood, ineqfun = ineq, ineqLB = 0, ineqUB = 1, 
                     LB = lx, UB = ux, control = control)
        opt$pars <- .augment.pars(opt$pars, j)
        opt$values[length(opt$values)] <- likelihood(opt$pars)
        
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)
        
      } else { # already know w = 1 (single component mixture)
        
        opt <- solnp(initial, fun = likelihood0, LB = lx, UB = ux, control = control)
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)
        
      }
    }
    
    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams, 
                                            mle.est = mle.est, j = j)
    Mix.mle <- Mix(dist = dist, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                   name = "Mix.boot")
    
    ran.gen <- function(dat, mle){ 
      rMix(n = length(dat), obj = mle)
    }
    
    # parametric bootstrap for bootstrapped scaled determinants                  
    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, sim = "parametric", ran.gen = ran.gen, mle = Mix.mle, ...)
    det0 <- bt$t0 # original determinant value
    det.boot.param <- bt$t # bootstrapped determinants
    det.boot.scaled <- det.boot.param/(sd(det.boot.param))
    
    q_lower <- quantile(det.boot.scaled, probs = ql)
    q_upper <- quantile(det.boot.scaled, probs = qu)
    
    # non parametric bootstrap for scaling of det0
    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, ...)
    det.boot.np <- bt$t
    det0.scaled <- det0/(sd(det.boot.np))
    
    if(det0.scaled >= q_lower && det0.scaled <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j == j.max){
      break
    }
    
  }
  
  .return.paramEst(j, j.max, dat, mle.est, values, conv, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}



# Purpose: plot method for 'paramEst' objects

plot.paramEst <- function(x, mixture = TRUE, components = TRUE, ylim = NULL, cex.main = 0.9,
                          ...){

  obj <- x
  dat <- attr(obj, "dat")
  class(dat) <- "rMix"
  n <- length(dat)
  pars <- attr(obj, "pars")
  dist <- attr(obj, "dist")
  dist_call <- get(paste("d", dist, sep = ""))
  ndistparams <- attr(obj, "ndistparams")
  formals.dist <- attr(obj, "formals.dist")
  discrete <- attr(obj, "discrete")
  
  # construct reasonable x values
  if(discrete == TRUE){
    x <- seq(min(dat), max(dat), by = 1)
  } else {
    x <- seq(min(dat), max(dat), by = 0.01)
  }
  
  theta.list.long <- vector(mode = "list", length = ndistparams)
  names(theta.list.long) <- formals.dist  
  for(i in 1:ndistparams){
    theta.list.long[[i]] <- matrix(pars[(i*obj):((1 + i)*obj - 1)], nrow = length(x),
                                   ncol = obj, byrow = TRUE)
  }
  theta.list.long$x <- x
  
  if(length(pars) == ndistparams){ # i.e. the estimated mixture complexity equals 1
    w <- 1
    y <- do.call(dist_call, theta.list.long)
  } else {
    w <- c(pars[1:(obj - 1)], 1-sum(pars[1:(obj-1)]))
    y <- do.call(dist_call, theta.list.long) %*% w
  } 
  
  if(is.null(ylim) || anyNA(ylim)){ # construct reasonable ordinate values
    plt.inv <- suppressWarnings(plot(dat, plot = FALSE, right = FALSE))
    mx <- max(plt.inv$density, y)
    ylim <- c(0, mx)
  }
  
  cols <- .get.colors(alpha = 0.9)[c(6, 7)]
  main <- paste("Estimated ", obj, " component '", dist, "' mixture model", sep = "")
  plot(dat, freq = FALSE, ylim = ylim, main = main, xlab = "Data", border = "darkgrey",
       components = FALSE, cex.main = cex.main, ...)
  
  if(components == TRUE){ # plot components
    
    y.comp <- matrix(w, nrow = length(x), ncol = obj, 
                byrow = TRUE) * do.call(dist_call, theta.list.long)
   
    mapply(function(i) lines(x = x, y = y.comp[, i], col = cols[1], lwd = 1, lty = 2),
           i = 1:obj)
  }
  
  if(mixture == TRUE){ # plot mixture
    
    if(discrete == TRUE) lines(theta.list.long$x, y, col = cols[2], lwd = 1, type = "b",
                               pch = 20, cex = 0.70)
    else lines(theta.list.long$x, y, col = cols[2], lwd = 1) 
    
  }
  
}



# Purpose: print method for 'paramEst' objects

print.paramEst <- function(x, ...){
  
  obj <- x
  j <- as.numeric(obj)
  pars <- attr(obj, "pars")
  values <- attr(obj, "values")
  conv <- attr(obj, "conv")
  dist <- attr(obj, "dist")
  ndistparams <- attr(obj, "ndistparams")
  formals.dist <- attr(obj, "formals.dist")
  mle.fct <- attr(obj, "mle.fct")
  niter <- length(values)
  funv <- values[niter]
  
  if(!is.null(mle.fct) & j == 1){ # parameters estimated via MLE.function
    
    pars <- c(1, pars)
    
    cmat <- matrix(pars, nrow = 1, ncol = ndistparams + 1)
    colnames(cmat) <- c("w", formals.dist)
    rownames(cmat) <- paste("Component ", 1, ":", sep = "")
    
    cat(paste("\nParameter estimation for a 1 component '", dist, "' mixture model:",
              "\nFunction value: ",
              format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
              "\n", "\n", sep = ""))
    printCoefmat(cmat, ...)
    cat("Optimization via user entered MLE-function.\n", strrep("-", 70), sep = "")
    
  } else { # parameters estimated numerically
    
    if(conv == 0) conv <- TRUE
    else conv <- FALSE
    
    # add last weight completely determined by the others as they sum to 1
    pars_complete <- numeric(j*(ndistparams + 1))
    if(j != 1){
      pars_complete[1:(j-1)] <- pars[1:(j-1)]
      pars_complete[j] <- 1-sum(pars[1:(j-1)])
      pars_complete[(j+1):length(pars_complete)] <- pars[j:length(pars)]
    } else {
      pars_complete[1] <- 1
      pars_complete[2:length(pars_complete)] <- pars[1:length(pars)]
    }
    
    cmat <- matrix(pars_complete, nrow = j, ncol = ndistparams +1)
    colnames(cmat) <- c("w", formals.dist)
    rownames(cmat) <- paste("Component ", 1:j, ":", sep = "")
    
    cat(paste("\nParameter estimation for a ", j, " component '", dist, "' mixture model:",
              "\nFunction value: ",
              format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
              "\n", "\n", sep = ""))
    printCoefmat(cmat, ...)
    
    if(conv == TRUE){
      cat(paste( "Converged in ", niter, " iterations.\n", sep = ""))
    } 
    else cat("Optimization algorithm did not converge.\n", sep = "")
  }
  cat(paste("\nThe estimated order is ", j, ".", sep = ""))
  
}