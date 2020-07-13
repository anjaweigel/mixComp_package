
## Purpose: LRT based method of estimating the mixture complexity (as well as 
##          the weights and component parameters) returning a 'paramEst' object
##          (using bootstrap)


mix.lrt <- function(obj, j.max = 10, B = 100, quantile = 0.95, control = c(trace = 0), ...){

  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max, quantile = quantile, Hankel = FALSE,
                          param = TRUE)
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, env = environment())
  
  
  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood <- .get.negloglik.dist(dat = dat, dist = dist, formals.dist = formals.dist, 
                                    ndistparams = ndistparams, dist_call)
  j0 <- 0
  
  repeat{
    
    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1
    
    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...
      
      mle.est0 <- mle.est1
      L0 <- L1
      conv.j0 <- conv.j1
      values.j0 <- values.j1
      
      # also need to pass over the restrictions as they will be used in the bootstrap
      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1
      
    } else {  # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)
      
      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      ineq.j0 <- restrictions.j0$ineq
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux
      
      if (!is.null(MLE.function)){ # Calculate MLE via the MLE function
        
        mle.est0 <- sapply(MLE.function, function(fun) fun(dat))
        L0 <- likelihood0(mle.est0)
        .printresultsMLE(mle.est0, dist, formals.dist, ndistparams, likelihood0)
        conv.j0 <- NULL
        values.j0 <- L0
        
      } else { # Calculate MLE of a j component mixture numerically
        
        initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                       formals.dist)
        
        opt <- solnp(initial.j0, fun = likelihood0, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)
        
        .printresults(opt, j0, dist, formals.dist, ndistparams)
        mle.est0 <- opt$pars 
        L0 <- likelihood0(mle.est0)
        conv.j0 <- opt$convergence
        values.j0 <- opt$values
        
      }
    }
    
    # optimization for j1. Starts from j1 = 2 so we always need to include weight 
    # restrictions in optimization
    
    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)
    
    opt <- solnp(initial.j1, fun = likelihood, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1, 
                 LB = lx.j1, UB = ux.j1, control = control)
    mle.est1 <- opt$pars <- .augment.pars(opt$pars, j1)
    L1 <- opt$values[length(opt$values)] <- likelihood(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values
    
    .printresults(opt, j1, dist, formals.dist, ndistparams)
    
    LRT <- 2*(L0 - L1)
    
    # parameters used for parametric bootstrap and corresponding 'Mix' object                 
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams, 
                                            mle.est = mle.est0, j = j0)
    Mix.mle <- Mix(dist = dist, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
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
      
      likelihood_boot <- .get.negloglik.dist(dat = dat, dist = dist, formals.dist = formals.dist, 
                                             ndistparams = ndistparams, dist_call)
      
      # calculate optimal parameters for j0
      
      if (j0 == 1){  # already know w = 1 (single component mixture)
        
        likelihood_boot0 <- .get.negloglik.dist.0(dat = dat, dist = dist, formals.dist = formals.dist, 
                                                  ndistparams = ndistparams, dist_call)
        
        if(!is.null(MLE.function)){ # Calculate MLE via the MLE function
          
          opt.boot0 <- sapply(MLE.function, function(fun) fun(dat))
          L0.boot <- likelihood_boot0(opt.boot0)
          
        } else {  # Calculate MLE of a j0 = 1 component mixture numerically
          
          initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                            formals.dist)
          
          opt.boot0 <- solnp(initial.boot0, fun = likelihood_boot0, LB = lx.j0, UB = ux.j0, 
                             control = control)$values
          L0.boot <- opt.boot0[length(opt.boot0)]
          
        }
        
      } else { # need to include weight restrictions in optimization for j0 != 1
          
        initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                          formals.dist)
        opt.boot0 <- solnp(initial.boot0, fun = likelihood_boot, ineqfun = ineq.j0, ineqLB = 0, 
                           ineqUB = 1, LB = lx.j0, UB = ux.j0, control = control)
        opt.boot0$pars <- .augment.pars(opt.boot0$pars, j0)
        L0.boot <- likelihood_boot(opt.boot0$pars)
          
      }
      
      # calculate optimal parameters for j1 (always need restrictions since j1 
      # starts from 2)
      
      initial.boot1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper,
                                        dist, formals.dist)
      
      opt.boot1 <- solnp(initial.boot1, fun = likelihood_boot, ineqfun = ineq.j1, ineqLB = 0, 
                         ineqUB = 1, LB = lx.j1, UB = ux.j1, control = control)
      opt.boot1$pars <- .augment.pars(opt.boot1$pars, j1)
      L1.boot <- likelihood_boot(opt.boot1$pars)
      
      return(2*(L0.boot - L1.boot))
    }
    
    
    bt <- boot(dat, statistic = stat, R = B, sim = "parametric", ran.gen = ran.gen,
               mle = Mix.mle, ...)
    LRT.boot <- bt$t

    qu <- quantile(LRT.boot, quantile)
    
    if(LRT <= qu){
      # so that the printed result reflects that the order j.max was actually estimated 
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if (j0 == j.max){
      break
    }
  }
  
  .return.paramEst(j0, j.max, dat, mle.est0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}
