#
#======================= Dimension-Adaptive Mixture Discriminant Analysis
#
# Copyright (C) 2021 - Michael Fop

damda <- function( learn, data,
                   K = learn$K,
                   H = 0:2, 
                   regularize = FALSE,
                   control_em = damda::control_em(),
                   control_reg = damda::control_reg(),
                   verbose = TRUE )
  # Automatic model selection via BIC
  #
{
  data <- as.matrix(data)
  N <- nrow(data)
  R <- ncol(data)
  if ( !all( exists("pro", learn), exists("mu", learn), exists("sigma", learn),  exists("K", learn) ) ) {
    stop("Object 'learn' should be a list with entries 'pro', 'mu', 'sigma', and 'K'")
  }
  # TODO: check pd sigma in input
  K <- learn$K
  if ( K != length(learn$pro) ) stop("Number of classes 'K' does not match parameter dimensions")
  Xnames <- colnames(learn$sigma[,,1, drop = FALSE])
  varnames <- colnames(data)
  classnames <- colnames(learn$pro)
  if ( is.null(Xnames) | is.null(varnames) ) {
    warning("Variables from the learning phase and/or variables in the test data should have names")
  }
  if ( is.null(Xnames) ) Xnames <- paste0("V", 1:ncol(learn$sigma[,,1]))
  if ( is.null(varnames) ) varnames <- paste0("V", 1:ncol(data))
  obs <- match( Xnames, varnames ) 
  ext <- match( setdiff(varnames, Xnames), varnames )   # additional variables
  if ( length(ext) == 0 ) ext <- FALSE
  classnames <- colnames(learn$mu)
  if ( is.null(classnames) ) classnames <- as.character(1:K)
  colnames(learn$mu) <- names(learn$pro) <- classnames
  
  regpar <- if ( regularize ) {
    control_reg(gamma = control_reg$gamma, alpha = control_reg$alpha, data = data, K = K)
  } else NULL 
  
  BIC <- rep(NA, length(H))
  names(BIC) <- H
  
  learn0 <- learn
  res <- list()
  
  if ( verbose ) {
    flush.console()
    pbar <- txtProgressBar(min = 0, max = length(H), style = 3)
    on.exit( close(pbar) )
    ipbar <- 0
  }
  
  # initialize class allocation
  hc_start <- if (N <= R) mclust::hc(data, modelName = "VII", use = "SVD") else mclust::hc(data, modelName = "VVV", use = "SVD")
  #
  # an option?
  # hc_start <- if (N <= R) mclust::hc(data, modelName = "VII", use = "VARS") else mclust::hc(data, modelName = "VVV", use = "VARS")
  
  for ( h in H ) {
    i <- match(h, H)
    temp <- try( em_damda(learn, data, h, 
                          regularize = regularize,
                          regpar = regpar,
                          ctrl_em = control_em,
                          hc_start = hc_start,
                          classnames),
                 silent = !control_em$printmsg )
    if ( class(temp) == "try-error") {
      res[[i]] <- NA
      BIC[i] <- NA
    } else {
      res[[i]] <- temp
      BIC[i] <- 2*temp$loglik - temp$npar*log(temp$N)
    }
    
    if ( verbose ) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }
  
  # best model
  if ( all(is.na(BIC) ) ) {
    out <- structure( list(learn = learn0, K = K, H = H,
                           parameters = NA,  z = NA, classification = NA,
                           loglik = NA, N = N, npar = NA,
                           obs = obs, ext = ext, BIC = BIC, bic = NA) )
    attr(out, "info") <- NA
    class(out) <- "damda"
    return(out)
  } else {
    best <- which.max(BIC)
    out <- res[[ best ]]
    out$bic <- BIC[best]
    out$BIC <- BIC
    out$learn <- learn0
    out$H <- names(BIC)[best]
    out$data <- data
  }
  class(out) <- "damda"
  return(out)
}



predict.damda <- function(object, newdata, ...)
{
  if (!inherits(object, "damda")) 
    stop("object not of class \"damda\"")
  if ( missing(newdata) ) {
    newdata <- object$data
  }
  newdata <- as.matrix(newdata)
  if ( ncol(object$data ) != ncol(newdata) ) {
    stop("newdata must match ncol of object data")
  }
  object$data <- newdata
  e <- estepdamda(newdata, t(object$parameters$mu), 
                  object$parameters$sigma, object$parameters$pro)
  z <- e$z
  cl <- seq(object$K)
  colnames(z) <- cl
  cl <- max.col(z)
  out <- list(classification = cl, z = z)
  return(out)
}
