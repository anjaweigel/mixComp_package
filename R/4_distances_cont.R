
## Purpose: estimating the adaptive Kernel density estimator found in Cutler & 
##          Cordero-Brana and sampling from it

ADAPkde <- function(dat, ndistparams, n, j, init, dist, formals.dist, dist_call,
                    sample.n, sample.plot, bs_iter = NULL){

  if(j == 1) w <- 1 
  else w <- c(init[1:(j - 1)], 1 - sum(init[1:(j - 1)]))
  
  theta.list.long <- vector(mode = "list", length = ndistparams)
  names(theta.list.long) <- formals.dist 
  for(i in 1:ndistparams){
    theta.list.long[[i]] <- matrix(init[(i*j):((1 + i)*j - 1)], nrow = n, ncol = j, byrow = TRUE)
  }
  theta.list.long$x <- dat
  
  # get matrix of a's given an estimate of the weights w and an estimate of the component
  # parameters theta.list.long (both supplied via 'init') at points dat
  
  w.array <- matrix(w, nrow = n, ncol = j, byrow = TRUE)
  f.array <- array(do.call(dist_call, args = theta.list.long), dim = c(n, j))
  a.array <- w.array * f.array / rowSums(array(w.array * f.array, dim = c(n, j)))
  # a.array will be infinite if f.array is infinite for any point (can happen when
  # solnp does not converge). 
  a.array[which(is.na(a.array))] <- max(a.array, na.rm = TRUE)
  
  # calculate \hat{\sigma}_i, i in 1,...j i.e. for every component as the empirical 
  # standard deviation of those observations whose weighted density is highest in 
  # component i
  
  ind <- apply(f.array, 1, function(x) which(x == max(x)))
  # if two components have the same value (should not really happen) choose the
  # component with the lower index
  if(is.list(ind)) ind <- sapply(ind, min)
  
  sigma <- mapply(function(i) sd(dat[ind == i]), 1:j)
  # default value in case no observation has its highest density in component i
  sigma[is.na(sigma)] <- 1 
  bw <- matrix(rep(2.283 * sigma / n^(0.283), each = n), nrow = n, ncol = j, byrow = FALSE)
  
  theta.list.long$mean <- matrix(dat, nrow = n, ncol = j, byrow = FALSE)
  theta.list.long$sd <- bw
  
  # adaptive kernel density estimate
  kde <- function(y, tll = theta.list.long, a = a.array, N = n){
    return(1/N*sum(a*dnorm(x = y, mean = tll$mean, sd = tll$sd)))
  }
  kde <- Vectorize(kde)
  
  # sampling from the adaptive kde  
  mean.ind <- sample(1:n, size = sample.n, replace = TRUE)
  means <- dat[mean.ind]
  weights <- array(a.array[mean.ind, ]/bw[1, ], dim = c(sample.n, j))
  sd.sample <- function(i) sample(1:j, size = 1, prob = weights[i,])
  sd.sample <- Vectorize(sd.sample)
  sd.ind <- sd.sample(1:sample.n)
  sds <- bw[1, sd.ind]
  sample <- rnorm(sample.n, means, sds)
  
  if(sample.plot == TRUE){
    
    if(is.null(bs_iter)){ # not in bootstrap loop yet
      
      txt <- "Sample from kde based on original data"
      hist(sample, freq = FALSE, breaks = 100, col = "light grey", 
           main = txt, xlab = "Sample Value")
      vals <- seq(min(dat), max(dat), length.out = 100)
      kdevals <- kde(vals)
      lines(vals, kdevals) 
      
    } else if(bs_iter != 0){ 
      # bs_iter equal 0 just recomputes statistic based on original data
      
      txt <- paste("Sample from bootstrap kde: iteration ", bs_iter, sep = "")
      hist(sample, freq = FALSE, breaks = 100, col = "light grey", 
           main = txt, xlab = "Sample Value")
      vals <- seq(min(dat), max(dat), length.out = 100)
      kdevals <- kde(vals)
      lines(vals, kdevals) 
      
    }
  }
  
  return(list(kde = kde, sample = sample))
  
} 



## Purpose: returns the approximated Hellinger distance function corresponding to 
##          parameters x

.get.fmin.hellinger.c <- function(kde, dat, formals.dist, ndistparams, dist,
                                  sample, kdevals, dist_call){

   fmin <- function(x){
     
    j <- (length(x) + 1)/(ndistparams + 1)
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
    if(any(w < 0)) return(0)

    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j - 1)], nrow = length(sample), ncol = j,
                                     byrow = TRUE)
    }
    theta.list.long$x <- sample

    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    fvals <- as.matrix(mat[, ]) %*% w
    # Additional constraint of not allowing all weight being put on a single observation;
    # this may happen because we approximate the integral with a sum and this sum can be
    # maximized by increasing the value of a single summand sufficiently while letting
    # the others go to 0
    if(any(is.na(fvals)) || any(fvals > 1000)) return(0)

    return (- sum(sqrt(fvals / kdevals)))
  }

}



## Purpose: returns the approximated Hellinger distance function corresponding to 
##          parameters x when the mixture consists of only a single component

.get.fmin.hellinger.c.0 <- function(kde, dat, formals.dist, ndistparams, dist, sample,
                                    kdevals, dist_call){
  
  fmin <- function(x){
    
    theta <- setNames(x, formals.dist)
    theta.list <- split(unname(theta), names(theta))
    theta.list$x <- sample
    
    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    fvals <- suppressWarnings(do.call(dist_call, args = theta.list))
    # Additional constraint of not allowing all weight being put on a single observation;
    # this may happen because we approximate the integral with a sum and this sum can be
    # maximized by increasing the value of a single summand sufficiently while letting
    # the others go to 0
    if(any(is.na(fvals)) || any(fvals > 1000)) return(0)
    
    return (- sum(sqrt(fvals / kdevals)))
  }
  
}



## Purpose: returns the true Hellinger distance between the kde and the mixture
#           corresponding to parameters x

.get.hellingerD <- function(x, j, ndistparams, formals.dist, kde, dist, values){

  if(length(x) == ndistparams) w <- 1
  else w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
  theta <- x[j:length(x)]
  
  theta <- setNames(as.vector(theta), rep(formals.dist, each = j))
  theta.list <- split(unname(theta), names(theta))
  Mix.obj <- Mix(dist, w = w, theta.list = theta.list)
  
  objective <- try(integrate(function(x){sqrt(dMix(x, Mix.obj) * kde(x))}, -Inf, Inf,
                             subdivisions = 1000L)[[1]], silent = TRUE)
  if(inherits(objective, "try-error")){
    cat(" Error while calculating the value of the true objective function. \n Returning the value of the approximated objective function instead. \n")
    objective <- values[length(values)]
  }

  return(2 - 2 * objective)
}



## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          continuous mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object

hellinger.cont <- function(obj, bandwidth = 1, j.max = 10, threshold = "SBC", 
                           sample.n = 3000, sample.plot = TRUE, control = c(trace = 0)){
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, env = environment())
  continuous <- !discrete
  
  # check relevant inputs
  .input.checks.functions(obj, bandwidth = bandwidth, j.max = j.max, thrshHel = threshold,
                          sample.n = sample.n, sample.plot = sample.plot,
                          continuous = continuous, Hankel = FALSE, param = TRUE)
  j0 <- 0

  if(is.character(threshold)){ 
    # otherwise it is a function and will be calculated further down
    if(threshold == "AIC") thresh <- (ndistparams + 1)/N
    if(threshold == "SBC") thresh <- ((ndistparams + 1) * log(N))/(2 * N)   
  }
  
  repeat{

    j0 <- j0 + 1  # current complexity estimate
    j1 <- j0 + 1
  
    if(is.function(threshold)){
      thresh <- threshold(n = N, j = j0)
    }
    
    initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)
    restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    lx.j0 <- restrictions.j0$lx
    ux.j0 <- restrictions.j0$ux
    
    if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)
      
      # will not really be used in the computation of the adaptive kde anyways (for j0 == 1);
      # just supplying something so the same function (ADAPkde) can be used
      if(j0 == 1) theta.j1 <- initial.j0
      
      # compute the kde and draw a sample from it
      kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.j1, dist, formals.dist, dist_call,
                          sample.n, sample.plot)
      kde <- kde.list$kde
      sample <- kde.list$sample
      
    } else { # use the standard gaussian Kernel estimate with the supplied bandwidth
      
      # compute the kde and draw a sample from it
      kde <- kdensity::kdensity(dat, bw = bandwidth, kernel = "gaussian")
      rkernel <- function(n) rnorm(n, sd = bandwidth)
      sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n) 
      
      if(sample.plot == TRUE){
        hist(sample, freq = FALSE, breaks = 100)
        lines(seq(min(dat), max(dat), length.out = 100), kde(seq(min(dat), max(dat), length.out = 100)))
      }
      
    }
    
    kdevals <- kde(sample)
    
    # calculate optimal parameters for j0
    if(j0 > 1){ # need to include weight restrictions in optimization
      
      fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                   kdevals, dist_call)
      ineq.j0 <- restrictions.j0$ineq

      opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                   LB = lx.j0, UB = ux.j0, control = control)

      
    } else { # already know w = 1 (single component mixture)
      
      fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)
      
      opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
    }
    
    # if we estimate multiple components check that all weights satisfy the constraints
    if(j0 != 1) theta.j0 <- opt$pars <- .augment.pars(opt$pars, j0)
    else theta.j0 <- opt$pars
    
    Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, 
                                                                      formals.dist, kde, dist, opt$values[length(opt$values)])
    conv.j0 <- opt$convergence
    values.j0 <- opt$values
    .printresults(opt, j0, dist, formals.dist, ndistparams)
    
    # calculate optimal parameters for j1 (always need weight restrictions since j1 
    # starts from 2)
    
    fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                  kdevals, dist_call)
    
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
    Hellinger.j1 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j1, j1,  ndistparams,
                                                                      formals.dist, kde, dist, opt$values[length(opt$values)])
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
                   discrete = !continuous, MLE.function)
}



## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          continuous mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object (using bootstrap)

hellinger.boot.cont <- function(obj, bandwidth = 1, j.max = 10, B = 100, ql = 0.025,
                                qu = 0.975, sample.n = 3000, sample.plot = TRUE,
                                control = c(trace = 0), ...){

  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, env = environment())
  continuous <- !discrete
  
  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max,  B = B, ql = ql, qu = qu,
                          continuous = continuous, Hankel = FALSE, param = TRUE)
  
  j0 <- 0 
  
  repeat{
  
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)
    restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    lx.j0 <- restrictions.j0$lx
    ux.j0 <- restrictions.j0$ux
    
    if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)
      
      # will not really be used in the computation of the adaptive kde anyways (for j0 == 1);
      # just supplying something so the same function (ADAPkde) can be used
      if(j0 == 1) theta.j1 <- initial.j0
      
      # compute the kde and draw a sample from it
      kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.j1, dist, formals.dist, dist_call,
                          sample.n, sample.plot)
      kde <- kde.list$kde
      sample <- kde.list$sample
      
    } else { # use the standard gaussian Kernel estimate with the supplied bandwidth
     
      # compute the kde and draw a sample from it
      kde <- kdensity::kdensity(dat, bw = bandwidth, kernel = "gaussian")
      rkernel <- function(n) rnorm(n, sd = bandwidth)
      sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n)
      
      if(sample.plot == TRUE){
        hist(sample, freq = FALSE, breaks = 100)
        lines(seq(min(dat), max(dat), length.out = 100), kde(seq(min(dat), max(dat), length.out = 100)))
      }
      
    }
    
    kdevals <- kde(sample)
    
    # calculate optimal parameters for j0
    if(j0 > 1){ # need to include weight restrictions in optimization
      
      fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                    kdevals, dist_call)
      ineq.j0 <- restrictions.j0$ineq
      
      opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                   LB = lx.j0, UB = ux.j0, control = control)
      
      
    } else { # already know w = 1 (single component mixture)
      
      fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)
      
      opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
    }
    
    # if we estimate multiple components check that all weights satisfy the constraints
    if(j0 != 1) theta.j0 <- opt$pars <- .augment.pars(opt$pars, j0)
    else theta.j0 <- opt$pars
    
    Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, 
                                                                      formals.dist, kde, dist, opt$values[length(opt$values)])
    conv.j0 <- opt$convergence
    values.j0 <- opt$values
    .printresults(opt, j0, dist, formals.dist, ndistparams)
    
    # calculate optimal parameters for j1 (always need weight restrictions since j1 
    # starts from 2)
    
    fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                  kdevals, dist_call)
    
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
    Hellinger.j1 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j1, j1,  ndistparams,
                                                                      formals.dist, kde, dist, opt$values[length(opt$values)])
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
    bs_iter <- - 1
    
    stat <- function(dat){
      
      assign("bs_iter", bs_iter + 1, inherits = TRUE)
      if(bs_iter != 0){
        
        # don't include first iteration as this just uses the original data
        # to calculate t0
        cat(paste("Running bootstrap iteration ", bs_iter, " testing for ", j0, 
                  " components.\n", sep = ""))
        
      } else cat(paste("\n"))
      
      # calculate optimal parameters for j0
      
      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
        
      
      if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)
        
        # compute the kde and draw a sample from it
        kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.j1, dist, formals.dist, dist_call,
                            sample.n, sample.plot, bs_iter)
        kde <- kde.list$kde
        sample <- kde.list$sample
        
      } else { # use the standard gaussian Kernel estimate with the supplied bandwidth
        
        # compute the kde and draw a sample from it
        kde <- kdensity::kdensity(dat, bw = bandwidth, kernel = "gaussian")
        rkernel <- function(n) rnorm(n, sd = bandwidth)
        sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n)
        
        if(sample.plot == TRUE){
          hist(sample, freq = FALSE, breaks = 100)
          lines(seq(min(dat), max(dat), length.out = 100), kde(seq(min(dat), max(dat), length.out = 100)))
        }

      }
      
      kdevals <- kde(sample)
      
      if(j0 > 1){ # need to include weight restrictions in optimization
        
        fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)
        opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)
        
      } else { # already know w = 1 (single component mixture)
        
        fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                        kdevals, dist_call)
        opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
      }
      
      # if we estimate multiple components check that all weights satisfy the constraints
      if(j0 != 1) theta.boot0 <- .augment.pars(opt$pars, j0)
      else theta.boot0 <- opt$pars
      
      Hellinger.boot0 <- .get.hellingerD(theta.boot0, j0,  ndistparams, formals.dist, kde, dist, opt$values[length(opt$values)])
      
      # calculate optimal parameters for j1 (always need weight restrictions since j1 
      # starts from 2)
      
      fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                    kdevals, dist_call)
      
      initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
      
      opt <- solnp(initial.j1, fun = fmin, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                   LB = lx.j1, UB = ux.j1, control = control)
      
      theta.boot1 <- .augment.pars(opt$pars, j1)
      Hellinger.boot1 <- .get.hellingerD(theta.boot1, j1,  ndistparams, formals.dist, kde, dist, opt$values[length(opt$values)])

      return(Hellinger.boot0 - Hellinger.boot1)
      
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
                   discrete =!continuous, MLE.function)
}
