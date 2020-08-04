
## Purpose: returns the L2 distance function corresponding to parameters x

.get.fmin.L2 <- function(dat, dist, formals.dist, ndistparams, j, n.inf, N, dist_call){

  function(x){

    w <- x[1:(j - 1)]

    # first term of difference: sum over values up to "n.inf"
    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j - 1)], nrow = (n.inf + 1), ncol = j,
                                     byrow = TRUE)
    }
    theta.list.long$x <- 0:n.inf

    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
    if(any(w < 0)) return(sqrt(.Machine$integer.max))
    f.theta <- as.matrix(mat) %*% w
    if(any(is.na(f.theta))) return(sqrt(.Machine$integer.max))
    f.theta.sq <- f.theta^2
    f.1 <- sum(f.theta.sq) # first term of difference
    
    # second term of difference: sum over the data as we multiply with the empirical
    #                            distribution function
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j - 1)], nrow = N, ncol = j,
                                     byrow = TRUE)
    }
    theta.list.long$x <- dat

    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
    if(any(w < 0)) return(sqrt(.Machine$integer.max))
    f.theta.obs <- as.matrix(mat) %*% w
    if(any(is.na(f.theta.obs))) return(sqrt(.Machine$integer.max))
    f.2 <- (2/N)*sum(f.theta.obs) # second term of difference

    return(f.1 - f.2)
  }
}



## Purpose: returns the L2 distance function corresponding to parameters x
##          when the mixture consists of only a single component

.get.fmin.L2.0 <- function(dat, dist, formals.dist, ndistparams, n.inf, N,  dist_call){
  
  function(x){
    
    # first term of difference: sum over values up to "n.inf"
    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- rep(x[i], n.inf + 1)
    }
    theta.list.long$x <- 0:n.inf
    
    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    f.theta <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    if(any(is.na(f.theta))) return(sqrt(.Machine$integer.max))
    f.theta.sq <- f.theta^2
    f.1 <- sum(f.theta.sq)
    
    # second term of difference: sum over the data as we multiply with the empirical 
    #                            distribution function
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- rep(x[i], N)
    }
    theta.list.long$x <- dat
    
    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    f.components.obs <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    if (any(is.na(f.components.obs))) return(sqrt(.Machine$integer.max))
    f.theta.obs <- sum(f.components.obs)
    f.2 <- (2/N)*sum(f.theta.obs)
    
    return(f.1 - f.2)
  }
} 



## Purpose: L2 distance based method of estimating the mixture complexity of a
##          discrete mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object

L2.disc <- function(obj, j.max = 10, n.inf = 1000, threshold = "SBC", control = c(trace = 0)){
 
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  # check relevant inputs
  .input.checks.functions(obj, thrshL2 = threshold, j.max = j.max, n.inf = n.inf, 
                          discrete = discrete, Hankel = FALSE, param = TRUE)
  j0 <- 0 
  
  repeat{
    
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    
    if(is.function(threshold)){
      thresh <- threshold(n = N, j = j0)
    }
    else if(threshold == "LIC"){
      thresh <- (0.6*log((j1)/j0))/N 
    } 
    else if (threshold == "SBC"){
      thresh <- (0.6*log(N)*log((j1)/j0))/N 
    }
    
    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...
      
      theta.j0 <- theta.j1
      L2.j0 <- L2.j1 
      conv.j0 <- conv.j1
      values.j0 <- values.j1
      
    } else { # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)
      
      fmin <- .get.fmin.L2.0(dat = dat, dist = dist, formals.dist = formals.dist, 
                             ndistparams = ndistparams, n.inf = n.inf, N = N, dist_call)
      
      restrictions <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      lx <- restrictions$lx
      ux <- restrictions$ux
      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
      
      opt <- solnp(initial.j0, fun = fmin, LB = lx, UB = ux, control = control)
      .printresults(opt, j0, dist, formals.dist, ndistparams)
      theta.j0 <- opt$pars
      L2.j0 <- opt$values[length(opt$values)]
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
    }
    
    # optimization for j1. Starts from j1 = 2 so we always need to include weight 
    # restrictions in optimization
    
    fmin <- .get.fmin.L2(dat = dat, dist = dist, formals.dist = formals.dist, 
                         ndistparams = ndistparams, j = j1, n.inf = n.inf, N = N, dist_call)
    
    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq <- restrictions$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)
    
    opt <- solnp(initial.j1, fun = fmin, LB = lx.j1, UB = ux.j1, ineqfun = ineq, 
                 ineqLB = 0, ineqUB = 1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    L2.j1 <- opt$values[length(opt$values)] <- fmin(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values
    
    .printresults(opt, j1, dist, formals.dist, ndistparams)
    
    if((L2.j0 - L2.j1) < thresh){
      break
    } else if(j0 == j.max){
      break
    }
    
  }
  
  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function = NULL)
}



## Purpose: L2 distance based method of estimating the mixture complexity of a
##          discrete mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object (using bootstrap)

L2.boot.disc <- function(obj, j.max = 10, n.inf = 1000, B = 100, ql = 0.025, qu = 0.975,
                         control = c(trace = 0), ...){
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())
  
  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max, B = B, n.inf = n.inf, ql = ql, qu = qu,
                          discrete = discrete, Hankel = FALSE, param = TRUE)
  j0 <- 0 
  
  repeat{
    
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...
      
      theta.j0 <- theta.j1
      L2.j0 <- L2.j1 
      conv.j0 <- conv.j1
      values.j0 <- values.j1
      
      # also need to pass over the restrictions as they will be used in the bootstrap
      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1
      
    } else { # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)
      
      fmin.j0 <- .get.fmin.L2.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                ndistparams = ndistparams, n.inf = n.inf, N = N, dist_call)
      
      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux
      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper,
                                     dist, formals.dist)
      
      opt <- solnp(initial.j0, fun = fmin.j0, LB = lx.j0, UB = ux.j0, control = control)
      .printresults(opt, j0, dist, formals.dist, ndistparams)
      theta.j0 <- opt$pars
      L2.j0 <- opt$values[length(opt$values)]
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
    }
    
    # optimization for j1. Starts from j1 = 2 so we always need to include weight 
    # restrictions in optimization

    fmin.j1 <- .get.fmin.L2(dat = dat, dist = dist, formals.dist = formals.dist,
                            ndistparams = ndistparams, j = j1, n.inf = n.inf, N = N, dist_call)
    
    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, 
                                   dist, formals.dist)
    
    opt <- solnp(initial.j1, fun = fmin.j1, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                 LB = lx.j1, UB = ux.j1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    L2.j1 <- opt$values[length(opt$values)] <- fmin.j1(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values
    
    .printresults(opt, j1, dist, formals.dist, ndistparams)
    
    diff.0 <- L2.j0 - L2.j1
    
    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams, 
                                            mle.est = theta.j0, j = j0)
    Mix.boot <- Mix(dist = dist, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                    name = "Mix.boot")
    
    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }
    
    # counting bootstrap iterations to print progression
    bs_iter <- -1
    
    stat <- function(dat){
      
      assign("bs_iter", bs_iter + 1, inherits = TRUE)
      if(bs_iter != 0){ 
        
        # don't include first iteration as this just uses the original data
        # to calculate t0
        cat(paste("Running bootstrap iteration ", bs_iter, " testing for ", j0, 
                  " components.\n", sep = ""))
        
      } else cat(paste("\n")) 
      
      # in the bootstrap we have to calculate the values for j0 and j1 as the bootstrap
      # data changes in every iteration (cannot reuse last j1 values as j0)
      initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper,
                                        dist, formals.dist)
      
      # calculate optimal parameters for j0
      if(j0 != 1){ # need to include weight restrictions in optimization
        
        fmin.boot0 <- .get.fmin.L2(dat = dat, dist = dist, formals.dist = formals.dist,
                                   ndistparams = ndistparams, j = j0, n.inf = n.inf, N = N, dist_call)
        opt.boot0 <- solnp(initial.boot0, fun = fmin.boot0, ineqfun = ineq.j0, ineqLB = 0,
                            ineqUB = 1, LB = lx.j0, UB = ux.j0, control = control)
        opt.boot0$pars <- .augment.pars(opt.boot0$pars, j0)
        L2.boot0 <- fmin.boot0(opt.boot0$pars)
        
      } else { # already know w = 1 (single component mixture)
        
        fmin.boot0 <- .get.fmin.L2.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                     ndistparams = ndistparams, n.inf = n.inf, N = N, dist_call)
        opt.boot0 <- solnp(initial.boot0, fun = fmin.boot0, LB = lx.j0, UB = ux.j0,
                            control = control)
        L2.boot0 <- opt.boot0$values[length(opt.boot0$values)]
        
      }
      
      # calculate optimal parameters for j1 (always need weight restrictions since j1 
      # starts from 2)
      fmin.boot1 <- .get.fmin.L2(dat = dat, dist = dist, formals.dist = formals.dist,
                                 ndistparams = ndistparams, j = j1, n.inf = n.inf, N = N, dist_call)
      
      initial.boot1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper,
                                        dist, formals.dist)
      
      opt.boot1 <- solnp(initial.boot1, fun = fmin.boot1, ineqfun = ineq.j1, ineqLB = 0, 
                         ineqUB = 1, LB = lx.j1, UB = ux.j1, control = control)
      opt.boot1$pars <- .augment.pars(opt.boot1$pars, j1)
      L2.boot1 <- fmin.boot1(opt.boot1$pars)
      
      return(L2.boot0 - L2.boot1)
    }
    
    bt <- boot(dat, statistic = stat, R = B, sim = "parametric", ran.gen = ran.gen,
               mle = Mix.boot, ...)
    diff.boot <- bt$t
    
    q_lower <- quantile(diff.boot, probs = ql)
    q_upper <- quantile(diff.boot, probs = qu)
    
    if(diff.0 >= q_lower && diff.0 <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if (j0 == j.max){
      break
    }
    
  }
  
  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function = NULL)
}



## Purpose: returns the squareroot of the empirical mass function needed for hellinger
##          distance calculation

.get.f.n.sqrt <- function(dat, n.max, N){
  
  # calculating square root of the empirical mass function
  f.n <- as.numeric(table(dat)[match(0:n.max, sort(unique(dat)))]/N)
  f.n[is.na(f.n)] <- 0
  f.n.sqrt <- sqrt(f.n)
  
}



## Purpose: returns the Hellinger distance function corresponding to parameters x

.get.fmin.hellinger <- function(dat, dist, formals.dist, ndistparams, j, n.max, N, 
                                f.n.sqrt, dist_call){
  
  function(x){
    
    w <- x[1:(j-1)]
    
    # calculating square root of mixture distribution corresponding to the parameters x
    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist  
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j-1)], nrow = (n.max+1), ncol = j, byrow = TRUE)
    }
    theta.list.long$x <- 0:n.max
    
    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
    if(any(w < 0)) return(sqrt(.Machine$integer.max))
    f.theta <- as.matrix(mat) %*% w
    if(any(is.na(f.theta))) return(sqrt(.Machine$integer.max))
    f.theta.sqrt <- sqrt(f.theta)
    
    # calculate Hellinger distance to empirical mass function
    H2 <- f.n.sqrt * f.theta.sqrt
    return(2 - 2 * sum(H2))
  }
}



## Purpose: returns the Hellinger distance function corresponding to parameters x
##          when the mixture consists of only a single component

.get.fmin.hellinger.0 <- function(dat, dist, formals.dist, ndistparams, n.max, N,
                                  f.n.sqrt, dist_call){
  
  function(x){
    
    # calculating square root of mixture distribution corresponding to the parameters x
    # (single component mixture)
    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist  
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- rep(x[i], n.max+1)
    }
    theta.list.long$x <- 0:n.max
    
    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    f.components <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    if(any(is.na(f.components))) return(sqrt(.Machine$integer.max))
    f.components.sqrt <- sqrt(f.components)
    
    # calculate Hellinger distance to empirical mass function
    H2 <- f.n.sqrt*f.components.sqrt
    return(2 - 2*sum(H2))
  }
}



## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          discrete mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object 

hellinger.disc <- function(obj, j.max = 10, threshold = "SBC", control = c(trace = 0)){
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())
  
  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max, thrshHel = threshold,
                          discrete = discrete, Hankel = FALSE, param = TRUE)
  j0 <- 0
  
  if(is.character(threshold)){ 
    # otherwise it is a function and will be calculated further down
    if(threshold == "AIC") thresh <- (ndistparams + 1)/N
    if(threshold == "SBC") thresh <- ((ndistparams + 1) * log(N))/(2 * N)   
  }
  
  repeat{
    
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    if(is.function(threshold)){
      thresh <- threshold(n = N, j = j0)
    }
    
    f.n.sqrt <- .get.f.n.sqrt(dat, n.max, N)
      
    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...
      
      theta.j0 <- theta.j1
      Hellinger.j0 <- Hellinger.j1
      conv.j0 <- conv.j1
      values.j0 <- values.j1
      
    } else { # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)
      
      fmin <- .get.fmin.hellinger.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                    ndistparams = ndistparams, n.max = n.max, N = N,
                                    f.n.sqrt = f.n.sqrt, dist_call)
      
      restrictions <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      lx <- restrictions$lx
      ux <- restrictions$ux
      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper,
                                     dist, formals.dist)
      
      opt <- solnp(initial.j0, fun = fmin, LB = lx, UB = ux, control = control)
      .printresults(opt, j0, dist, formals.dist, ndistparams)
      theta.j0 <- opt$pars
      Hellinger.j0 <- opt$values[length(opt$values)]
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
    }
    
    # optimization for j1. Starts from j1 = 2 so we always need to include weight 
    # restrictions in optimization
    
    fmin <- .get.fmin.hellinger(dat = dat, dist = dist, formals.dist = formals.dist,
                                ndistparams = ndistparams, j = j1, n.max = n.max, N = N,
                                f.n.sqrt = f.n.sqrt, dist_call)
    
    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)
    
    opt <- solnp(initial.j1, fun = fmin, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                 LB = lx.j1, UB = ux.j1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    Hellinger.j1 <- opt$values[length(opt$values)] <- fmin(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values
    
    .printresults(opt, j1, dist, formals.dist, ndistparams)
    
    if((Hellinger.j0 - Hellinger.j1) < thresh){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j0 == j.max){
      break
    }
    
  }
  
  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}

## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          discrete mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object (using bootstrap)

hellinger.boot.disc <- function(obj, j.max = 10, B = 100, ql = 0.025, qu = 0.975, 
                                control = c(trace = 0), ...){
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())
  
  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max, B = B, ql = ql, qu = qu,
                          discrete = discrete, Hankel = FALSE, param = TRUE)
  j0 <- 0
  
  repeat{
    
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    f.n.sqrt <- .get.f.n.sqrt(dat, n.max, N)
    
    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...
      
      theta.j0 <- theta.j1
      Hellinger.j0 <- Hellinger.j1 
      conv.j0 <- conv.j1
      values.j0 <- values.j1
      
      # also need to pass over the restrictions as they will be used in the bootstrap
      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1
      
    } else { # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)
      
      fmin.j0 <- .get.fmin.hellinger.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                       ndistparams = ndistparams, n.max = n.max, N = N,
                                       f.n.sqrt = f.n.sqrt, dist_call)
      
      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux
      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
      
      opt <- solnp(initial.j0, fun = fmin.j0, LB = lx.j0, UB = ux.j0, control = control)
      .printresults(opt, j0, dist, formals.dist, ndistparams)
      theta.j0 <- opt$pars
      Hellinger.j0 <- opt$values[length(opt$values)] 
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
    }
    
    # optimization for j1. Starts from j1 = 2 so we always need to include weight 
    # restrictions in optimization
    
    fmin.j1 <- .get.fmin.hellinger(dat = dat, dist = dist, formals.dist = formals.dist,
                                   ndistparams = ndistparams, j = j1, n.max = n.max, N = N,
                                   f.n.sqrt = f.n.sqrt, dist_call)
    
    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist, 
                                   formals.dist)
    
    opt <- solnp(initial.j1, fun = fmin.j1, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1, 
                 LB = lx.j1, UB = ux.j1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    Hellinger.j1 <- opt$values[length(opt$values)] <- fmin.j1(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values
    
    .printresults(opt, j1, dist, formals.dist, ndistparams)
    
    diff.0 <- Hellinger.j0 - Hellinger.j1
    
    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams, 
                                            mle.est = theta.j0, j = j0)
    Mix.boot <- Mix(dist = dist, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                    name = "Mix.boot")
    
    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }
    
    # counting bootstrap iterations to print progression
    bs_iter <- -1
    
    stat <- function(dat){
      
      assign("bs_iter", bs_iter + 1, inherits = TRUE)
      if(bs_iter != 0){
        
        # don't include first iteration as this just uses the original data
        # to calculate t0
        cat(paste("Running bootstrap iteration ", bs_iter, " testing for ", j0, 
                  " components.\n", sep = ""))
        
      } else cat(paste("\n"))
      
      f.n.sqrt.boot <- .get.f.n.sqrt(dat, n.max, N)
      
      # in the bootstrap we have to calculate the values for j0 and j1 as the bootstrap
      # data changes in every iteration (cannot reuse last j1 values as j0)
      
      initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                        formals.dist)
      
      # calculate optimal parameters for j0
      if(j0 != 1){ # need to include weight restrictions in optimization
        
        fmin.boot0 <- .get.fmin.hellinger(dat = dat, dist = dist, formals.dist = formals.dist,
                                          ndistparams = ndistparams, j = j0, n.max = n.max, N = N,
                                          f.n.sqrt = f.n.sqrt.boot, dist_call)
        opt.boot0 <- solnp(initial.boot0, fun = fmin.boot0, ineqfun = ineq.j0, ineqLB = 0, 
                           ineqUB = 1, LB = lx.j0, UB = ux.j0, control = control)
        opt.boot0$pars <- .augment.pars(opt.boot0$pars, j0)
        hellinger.boot0 <- fmin.boot0(opt.boot0$pars)
        
      } else { # already know w = 1 (single component mixture)
        
        fmin.boot0 <- .get.fmin.hellinger.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                            ndistparams = ndistparams, n.max = n.max, N = N,
                                            f.n.sqrt = f.n.sqrt.boot, dist_call)
        opt.boot0 <- solnp(initial.boot0, fun = fmin.boot0, LB = lx.j0, UB = ux.j0,
                           control = control)
        hellinger.boot0 <- opt.boot0$values[length(opt.boot0$values)]
        
      }
      
      # calculate optimal parameters for j1 (always need weight restrictions since j1 
      # starts from 2)
      
      fmin.boot1 <- .get.fmin.hellinger(dat = dat, dist = dist, formals.dist = formals.dist,
                                        ndistparams = ndistparams, j = j1, n.max = n.max, N = N,
                                        f.n.sqrt = f.n.sqrt.boot, dist_call)
      
      initial.boot1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper,
                                        dist, formals.dist)
      
      opt.boot1 <- solnp(initial.boot1, fun = fmin.boot1, ineqfun = ineq.j1, ineqLB = 0, 
                         ineqUB = 1, LB = lx.j1, UB = ux.j1, control = control)
      opt.boot1$pars <- .augment.pars(opt.boot1$pars, j1)
      hellinger.boot1 <- fmin.boot1(opt.boot1$pars)
      
      return(hellinger.boot0 - hellinger.boot1)
      
    }
    
    bt <- boot(dat, statistic = stat, R = B, sim = "parametric", ran.gen = ran.gen,
               mle = Mix.boot, ...)
    diff.boot <- bt$t
    
    q_lower <- quantile(diff.boot, probs = ql)
    q_upper <- quantile(diff.boot, probs = qu)
    
    if(diff.0 >= q_lower && diff.0 <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if (j0 == j.max){
      break
    }
    
  }
  
  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}
