# Purpose: variable that are defined globally within this package

globalVariables(c("dist", "theta.bound.list", "formals.dist", "ndistparams", "dist_call",
                  "bounds", "lower", "upper",
                  "dat", "N", "n.max", "discrete", "continuous",
                  "Hankel.method", "Hankel.function", "MLE.function"))


# Purpose: Gives colors for plots

.get.colors <- function(alpha){
  
  cols <- list(c(0, 0, 1), # blue
               c(0, 1, 0), # green
               c(0, 1, 1), # turqoise
               c(1, 0.6, 0), # orange
               c(0.57, 0.1, 0.33), # lilac
               c(0.35, 0, 0.22), # dark lilac 
               c(0, 0.3, 0.3), # dark turqoise
               c(0, 0.5, 0.5)) # medium turqoise
  cols <- sapply(cols, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
  
}


# Purpose: Constructor for 'Mix' (normal mixture) objects

Mix <- function(dist, w = NULL, theta.list = NULL, name = NULL, ...){ 
  
  if(!is.character(dist)) stop("'dist' needs to be a character string specfying the 
                               distribution of the mixture components!")
  if(!is.null(theta.list))
    if(!is.list(theta.list))
      stop("Input to theta.list has to be of class 'list'!")
    else theta.list <- theta.list
  # component parameters can also be entered via ...
  else theta.list <- list(...)
  
  if(is.null(unlist(theta.list)))
    stop("The component parameters have to be entered as 'theta.list' or as input to ...!")
  
  formals.dist <- names(theta.list)
  ndistparams <- length(formals.dist)
  
  # check if the vectors of component parameters (e.g. mean and sd for normal distribution)
  # are of equal length
  equal.length <- function(x){
    len <- sapply(x, length)
    diff(range(len)) < .Machine$double.eps
  }
  if(!equal.length(theta.list)) 
    stop("The elements of theta.list (or the inputs to ...) must be of equal length!")
  
  if(!is.numeric(unlist(theta.list)))
    stop("The elements of theta.list (or the inputs to ...) must all be numeric!")
  
  # true mixture complexity
  p <- length(theta.list[[1]])
  
  if(is.null(w)){ # construct default value for the weights
    w <- rep.int(1/p, p)
  } else {
    
    if(length(w) != p || !is.numeric(w) || any(w < 0))
      stop("'w' must be a numeric >= 0 with same length as the elements of theta.list 
           (or the inputs to ...)")
    s <- sum(w)
    if(abs(s - 1) > .Machine$double.eps) w <- w/s
    
  }
  
  if(is.null(name)) { # construct default name
    
    sformat <- function(v) sapply(v, format, digits = 1)
    pPar <- function(pp) {
      pp <- if(p >= 10) c(sformat(pp[1:9]), "....") else sformat(pp)
      paste(pp, collapse = "'")
    }
    name <- paste0(dist, "Mix", format(p))
    for(i in 1:ndistparams){
      name <- paste0(name, "_", pPar(theta.list[[i]]))
    }
    
  }
  if(!is.character(name)) name <- as.character(name)
  
  # construct matrix to be printed with weights as first column and other parameters
  # as further columns
  dat <- matrix(NA, nrow = p, ncol = ndistparams + 1)
  colnames(dat) <- as.character(1:(ndistparams + 1))
  for(i in 1:ndistparams){
    dat[, i+1] <- theta.list[[i]]
    colnames(dat)[i+1] <- formals.dist[i]
  }
  dat[, 1] <- w
  colnames(dat)[1] <- "w"
  
  # check whether there exists a function generating random numbers and returning the
  # probability denisty/mass for the distribution specified as 'dist' and if their
  # arguments are compatible with the component parameter names specified in 'theta.list'
  # (or via ...)
  dist_call_r <- try(get(paste("r", dist, sep = "")), silent = TRUE)
  dist_call_d <- try(get(paste("d", dist, sep = "")), silent = TRUE)
  if(inherits(dist_call_r, "try-error") || inherits(dist_call_d, "try-error"))
    stop("combining the string \"dist\" with \"d\" or \"r\" has to yield an existing function name!")
  
  if(!all(formals.dist %in% names(formals(dist_call_r))))
    stop(paste("The names of theta.list do not match the names of the formal arguments
               of the function r", dist, sep = ""))
  if(!all(formals.dist %in% names(formals(dist_call_d))))
    stop(paste("The names of theta.list do not match the names of the formal arguments
               of the function d", dist, sep = ""))
  
  # check whether the resulting mixture distribution is discrete
  is.int <- function(j){
    is.integer(do.call(dist_call_r, 
                       args = c(n = 1, lapply(theta.list, function (x) x[j]))))
  }
  rand <- mapply(is.int, 1:p)
  if(sum(rand) == p) discrete <- TRUE
  else discrete <- FALSE
  
  structure(name = name, dist = dist, discrete = discrete, theta.list = theta.list, class = "Mix", 
            .Data = dat)
}


## Purpose: is 'obj' a "Mix", i.e. a mixture object?

is.Mix <- function(x){
  
  obj <- x
  dist_call_d <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  dist_call_r <- try(get(paste("r", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call_d, "try-error") || inherits(dist_call_r, "try-error")) return(FALSE)
  theta.list <- attr(obj, "theta.list")
  formals.dist <- names(theta.list)
  
  inherits(obj, "Mix") && is.matrix(obj) && 
    (!is.null(w <- obj[, "w"])) && dim(obj)[1] == length(w) &&
    dim(obj)[2] == (length(theta.list) + 1) &&
    is.numeric(w) && all(w >= 0) && abs(sum(w) - 1) < 1000*.Machine$double.eps &&
    length(w) == length(theta.list[[1]]) && 
    length(unique(sapply(theta.list, function(x) length(x)))) == 1 && 
    all(formals.dist %in% names(formals(dist_call_d))) && 
    all(formals.dist %in% names(formals(dist_call_r))) &&
    is.logical(attr(obj, "discrete"))
  
}


## Purpose: density evaluation for "Mix" objects (mixtures)

dMix <- function(x, obj, log = FALSE){

  if(!is.numeric(x)) stop("'x' has to be numeric!")
  if(!is.Mix(obj)) stop("'obj' has to be a 'Mix' object!")
  if(!is.logical(log)) stop("'log' has to be logical!")
  
  w <- obj[,"w"]
  p <- length(w) # number of components
  theta.list <- attr(obj, "theta.list")
  theta.list <- lapply(theta.list, function(y) matrix(y, nrow = length(x), ncol = p,
                                                      byrow = TRUE))
  theta.list$x <- x
  
  dist_call <- get(paste("d", attr(obj, "dist"), sep = ""))
  y <- rowSums(matrix(w, nrow = length(x), 
                      ncol = p, byrow = TRUE) * do.call(dist_call, theta.list))
  if(log) log(y) else y
  
}


## Purpose: Generate random numbers according to "Mix" object;
##          simultaneously creates an "rMix" object

rMix <- function(n, obj){
  
  if(!is.Mix(obj)) stop("'obj' has to be a 'Mix' object!")
  if(!is.numeric(n) || !(as.integer(n) == n) || length(n) != 1) 
    stop("'n' has to be a single integer!")
  
  w <- obj[, "w"]
  p <- length(w)
  theta.list <- attr(obj, "theta.list")
  
  ind <- sample(1:p, prob = w, size = n, replace = TRUE)
  theta.list.expanded <- lapply(theta.list, function(x) x[ind]) 
  
  dist_call <- get(paste("r", attr(obj, "dist"), sep = ""))
  dat <- do.call(dist_call, args = c(n = n, theta.list.expanded))
  
  att <- attributes(obj)
  attributes(dat) <- att[names(att) != "class" & names(att) != "dimnames" & names(att) != "dim"]
  attr(dat, "w") <- w
  attr(dat, "indices") <- ind
  class(dat) <- "rMix"
  
  return(dat)
}


## Purpose: print method for "Mix" objects (mixtures)

print.Mix <- function(x, ...){
  
  obj <- x
  if(!is.Mix(obj)) stop("obj is not a 'Mix' object!")
  
  cat(paste("'", attr(obj, "dist"), sep = ""), "Mixture' object",
      paste("\t ``", attr(obj, "name"), "''", sep=''), "\n")
  
  att <- attributes(obj); 
  att <- att[names(att) != "dist" & names(att) != "theta.list" & names(att) != "discrete" & names(att) != "name"]
  attributes(obj) <- if(length(att) > 0) att

  class(obj) <- character(0)
  print(obj, ...)
  invisible(x)
  
}


## Purpose: plot method for "Mix" objects (mixtures)

plot.Mix <- function(x, ylim, xlim = NULL, xout = NULL, n = 511, type = NULL, 
                     xlab = "x", ylab = "f(x)", main = attr(obj, "name"), lwd = 1.4,
                     log = FALSE, components = TRUE, h0 = FALSE,
                     parComp = list(col = NULL, lty = 3, lwd = 1),
                     parH0 = list(col = NULL, lty = 3, lwd = 1), ...){

  obj <- x
  if(!is.numeric(n)) stop("'n' has to be an integer!")
  if(!is.logical(components)) stop("'components' has to be logical!")
  if(!is.list(parH0)) stop("'parH0' has to be a list!")
  if(!is.list(parComp)) stop("'parComp' has to be a list!")
  
  if(components == TRUE && is.null(parComp$col)){ # default colors
    parComp$col <- .get.colors(alpha = 0.9)[8]
  }
  
  if(h0 == TRUE && is.null(parH0$col)){ # default colors
    parH0$col <- .get.colors(alpha = 0.9)[6]
  }
  
  if(is.null(type)){ # default type
    if (attr(obj, "discrete") == FALSE) type = "l"
    else type = "h"
  }

  if(is.null(xlim) && is.null(xout)){  # construct "reasonable" abscissa values
    
    set.seed(1)
    rand <- rMix(1000, obj = obj)
    d.o <- dMix(obj, x = rand)
    if (attr(obj, "discrete") == FALSE)
      rand_trunc <- rand[d.o >= max(d.o)/200] # discard values with too low density
    else rand_trunc <- rand
    xlim <- c(min(rand_trunc), max(rand_trunc))
    
  } # here we definitely have xlim values
  
  if(!is.null(xlim) && is.null(xout)){ # construct xout from xlim
    
    if (attr(obj, "discrete") == FALSE)
      xout <- seq(xlim[1], xlim[2], length = n)
    else xout <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
    
  } # here we definitely have xout values
  
  
  if(log == TRUE & components == TRUE){
    warning("Only either 'log' or 'components' can take the value true. Setting components = FALSE")
    components <- FALSE
  }

  d.o <- dMix(obj, x = xout, log = log)
  
  if(missing(ylim) || anyNA(ylim)){ # construct "reasonable" ordinate values
    if(log == FALSE)
      ylim <- c(0, max(d.o))
    else ylim <- c(min(d.o), max(d.o))
  }
  
  # plot mixture density
  plot(y = d.o, x = xout, type = type, xlim = xlim, ylim = ylim,
       main = main, xlab = xlab, ylab = ylab, lwd = lwd, ...)
  
  if(components == TRUE) { # plot individual components
    
    w <- obj[ ,"w"]
    p <- length(w) 
    dist_call <- get(paste("d", attr(obj, "dist"), sep = ""))
    theta.list <- attr(obj, "theta.list")
    theta.list <- lapply(theta.list, function(y) matrix(y, nrow = length(xout), ncol = p,
                                                        byrow = TRUE))
    theta.list$x <- xout
    y <- matrix(w, nrow = length(xout), ncol = p, 
                byrow = TRUE) * do.call(dist_call, theta.list)
    mapply(function(i) do.call(lines, c(list(x = xout, y = y[, i]), parComp)), i = 1:p)
    
  }
  
  if(h0 == TRUE) do.call(abline, c(list(h = 0), parH0))
}


## Purpose: print method for "rMix" objects (random numbers generiated via a mixture)

print.rMix <- function(x, ...){
  print(as.vector(x), ...)
  return(x)
}


## Purpose: is 'obj' a "rMix" object?

is.rMix <- function(x){
  
  obj <- x
  dist_call_d <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  dist_call_r <- try(get(paste("r", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call_d, "try-error") || inherits(dist_call_r, "try-error")) return(FALSE)
  
  theta.list <- attr(obj, "theta.list")
  formals.dist <- names(theta.list)
  
  inherits(obj, "rMix") && is.numeric(obj) &&
    (!is.null(w <- attr(obj, "w"))) &&
    is.numeric(w) && all(w >= 0) && abs(sum(w)-1) < 1000*.Machine$double.eps &&
    length(w) == length(theta.list[[1]]) && 
    length(unique(sapply(theta.list, function(x) length(x)))) == 1 && 
    all(formals.dist %in% names(formals(dist_call_d))) && 
    all(formals.dist %in% names(formals(dist_call_r))) && 
    is.logical(attr(obj, "discrete")) && all(unique(attr(obj, "indices")) %in% 1:length(w))
}


## Purpose: plot method for "rMix" objects (random numbers generiated via a mixture)

plot.rMix <- function(x, xlab = attr(obj, "name"), ylim = NULL,
                      main = paste("Histogram of", attr(obj, "name")), breaks = NULL,
                      col = "grey", components = TRUE, stacked = FALSE, component.colors = NULL, 
                      freq = TRUE, plot = TRUE, ...){
  
  obj <- x
  if(!is.logical(components)) stop("'components' has to be logical!")
  if(!is.logical(stacked)) stop("'stacked' has to be logical!")
  
  # from which component did an observation arise?
  ind <- attr(obj, "indices")
  
  if(is.null(breaks)){ # create default breaks
    
    hist <- hist(obj, plot = FALSE)
    breaks <- hist$breaks
    diffbreaks <- unique(diff(breaks)*0.5)
    nbreaks <- length(breaks)
    breaks <- rep(breaks, each = 2) + rep(c(0, diffbreaks), nbreaks) 
    # making breaks twice as fine so that components can also be shown nicely 
  }
  
  if(is.null(component.colors) && components == TRUE){ # default colors
    if(stacked == FALSE)
      component.colors <- .get.colors(alpha = 0.4)
    else component.colors <- .get.colors(alpha = 0.3)
  }
  
  if(freq == FALSE && components == TRUE){ # if we plot probabilities and individual components
    if(stacked == TRUE){
      warning("Stacked components cannot be shown when 'freq' is FALSE. Setting freq = TRUE.")
      freq <- TRUE
    } else if(is.null(ylim) || anyNA(ylim)){ 
      # ordinate values such that component probabilities surely fit in the plot
      plt.inv <- hist(obj, breaks = breaks, plot = FALSE, ...)
      breaks <- plt.inv$breaks
      ylim <- c(0, 1/min(diff(breaks)))
    } 
  }
  
  if(plot == FALSE){
    plt.inv <- hist(obj, breaks = breaks, plot = plot, ...)
    return(plt.inv)
  } else {
    hist(obj, col = col, breaks = breaks, main = main, xlab = xlab, freq = freq, 
         ylim = ylim, plot = plot, right = FALSE, ...)
    # invisible plot to reuse breaks later for components
    plt.inv <- suppressWarnings(hist(obj, breaks = breaks, plot = FALSE, ...))
  }
  
  if(components == TRUE){
    
    breaks <- plt.inv$breaks
    p <- length(attr(obj, "theta.list")[[1]]) # number of components
    
    if(stacked == FALSE){ # just plot the individual histograms
      for (i in 1:p){
        hist(obj[ind == i], col = component.colors[i], breaks = breaks,
             add = TRUE, freq = freq, plot = plot, ...)
      }
      
    } else { # plot stacked rectangles
      
      nbreaks <- length(breaks)
      mat <- matrix(0, nrow = p + 1, ncol = nbreaks - 1)
      # cumulative version of mat
      matcum <- matrix(0, nrow = p + 1, ncol = nbreaks - 1)
      
      for(i in 2:(p + 1)){ # for every component i...
        
        # ... how many observations lie within first bin 
        # (seperate because first breaks value need to be included)
        mat[i, 1] <- length(obj[obj <= breaks[2] & obj >= breaks[1] & ind == (i - 1)])
        matcum[i, 1] <- matcum[i - 1, 1] + mat[i, 1]
        
        # ... how many observations lie within other bins
        mat[i, 2:ncol(mat)] <- mapply(function(j) length(obj[obj <= breaks[j + 1] &
                                                          obj > breaks[j] &
                                                          ind == (i - 1)]), 2:ncol(mat))
        matcum[i, 2:ncol(mat)] <- mapply(function(j) matcum[i - 1, j] + mat[i, j], 2:ncol(mat))
          
        
        rect(xleft = breaks[1:(nbreaks - 1)], ybottom = matcum[i - 1, ], xright = breaks[2:(nbreaks)], 
             ytop = matcum[i, ], col = component.colors[i - 1])
      }
      
    }
  }
}


# Purpose: Constructor for 'datMix' objects, to be passed to functions estimating
#          the mixture complexity, contains all "static" information about the data

datMix <- function(dat, dist, theta.bound.list = NULL, MLE.function = NULL,
                   Hankel.method = NULL, Hankel.function = NULL){  
  
  ndistparams <- length(theta.bound.list)
  formals.dist <- names(theta.bound.list)
  
  if(!is.null(theta.bound.list)){
    # check whether distribution is discrete
    dist_call <- get(paste("r", dist, sep = ""))
    rd.list <- theta.bound.list
    rd.list <- lapply(rd.list, function(x) replace(x, x == Inf, 100))
    rd.list <- lapply(rd.list, function(x) replace(x, x == -Inf, -100))
    rd.list <- lapply(rd.list, function(x) runif(n = 1, min = x[1], max = x[2]))
    rand <- do.call(dist_call, args = c(n = 1, rd.list))
    if(is.integer(rand)) 
      discrete <- TRUE
    else discrete <- FALSE
    
  }
  
  # return 'datMix' object with relevant attributes
  structure(class = "datMix", dist = dist, 
            theta.bound.list = theta.bound.list, discrete = discrete,
            MLE.function = MLE.function, Hankel.method = Hankel.method, 
            Hankel.function = Hankel.function, .Data = dat)
}


## Purpose: is 'obj' a "datMix" object?

is.datMix <- function(x){
  
  obj <- x
  dist_call <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call, "try-error")) return(FALSE)
  
  inherits(obj, "datMix") && is.numeric(obj)
}


## Purpose: print method for "datMix" objects (don't print attributes)

print.datMix <- function(x, ...){
  print(as.vector(x), ...)
}                                


## Purpose: Convert "RMix" object to "datMix" object

RtoDat <- function(obj, theta.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL, 
                   Hankel.function = NULL){
  
  new <- obj
  att <- attributes(obj)
  # discard other information like component parameter values as these are 
  # supposed to be estimated from the "datMix" object
  attributes(new) <- att[names(att) == "dist" | names(att) == "discrete"]
  class(new) <- "datMix"
  attr(new, "theta.bound.list") <- theta.bound.list
  attr(new, "MLE.function") <- MLE.function
  attr(new, "Hankel.method") <- Hankel.method
  attr(new, "Hankel.function") <- Hankel.function
  return(new)
  
}
