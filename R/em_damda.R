#
#======================= EM functions for damda
#
# Copyright (C) 2021 - Michael Fop

em_damda <- function(learn, data, H,
                     regularize = FALSE,
                     regpar = NULL,
                     ctrl_em = control_em(),
                     hc_start = NULL,
                     classnames)
  # Function implementing the EM algorithm for the Dimension-Adaptive Mixture Discriminant Analysis
  # NOTE: the estimation method is the INDUCTIVE approach
  #
{
  Xnames <- colnames(learn$sigma[,,1])
  varnames <- colnames(data)
  obs <- match( Xnames, varnames ) 
  ext <- match( setdiff(varnames, Xnames), varnames )   # additional variables
  if ( length(ext) == 0 ) ext <- FALSE
  
  R <- length(varnames)
  P <- length(Xnames)
  Q <- R - P
  N <- nrow(data)
  K <- length(learn$pro)
  C <- H + K
  
  # set parameters list--------------------------------------------
  parameters <- list()
  parameters$pro <- rep(NA, C)
  parameters$mu <- matrix(NA, R, C)
  parameters$sigma <- array(NA, c(R,R,C))
  # ---------------------------------------------------------------
  
  # fix parameters from training set-------------------------------
  parameters$mu[obs,1:K] <- learn$mu
  parameters$pro[1:K] <- learn$pro
  for ( l in 1:K ) parameters$sigma[obs,obs,l] <- learn$sigma[,,l]
  # ---------------------------------------------------------------
  
  
  # discovery phase------------------------------------------------
  if ( H == 0 & is.logical(ext) ) {
    # no new classes nor new variables ---> just a simple E step for prediction
    es <- estepdamda(data, t(parameters$mu), parameters$sigma, parameters$pro)
    z <- es$z
    loglik <- es$loglik
    it <- 1
    
  } else {
    # yes/no new classes --- yes/no new variables
    #
    # EM settings
    loglikPrev <- -.Machine$integer.max/2
    it <- 0
    crit <- TRUE
    
    # parameters initialization
    init <- init_em(learn, parameters, data, H, hc_start = hc_start, init = ctrl_em$init)
    if ( regularize ) {
      for ( j in 1:C ) {
        init$sigma[,,j] <- ( N*init$pro[j]*init$sigma[,,j] + regpar$scale*( (1/C)^(1/R) ) ) / (N*init$pro[j])
      }
    }
    parameters$mu <- init$mu
    parameters$sigma <- init$sigma
    parameters$pro <- init$pro
    
    es <- try( estepdamda(data, t(parameters$mu), parameters$sigma, parameters$pro),
               silent = TRUE )
    if ( regularize ) {
      Nc <- colSums(es$z)
      den <- sum( Nc + regpar$alpha )
      es$z <- es$z*(N/den) + (regpar$alpha/den)
    }
    
    if ( class(es) == "try-error") {
      out <- structure( list(learn = learn, K = K, H = H, parameters = parameters,
                             z = es$z, classification = mclust::map(es$z),
                             loglik = NA, N = N, npar = npar,
                             obs = obs, ext = ext) )
      return(out)
    }
    
    # EM ###############################################
    while ( crit ) {
      # E step
      z <- es$z
      
      # M step
      parameters <- mstep_damda(data, ext, z, K, H, learn$pro, parameters, regularize, regpar)
      
      # loglik
      es <- estepdamda(data, t(parameters$mu), parameters$sigma, parameters$pro)
      loglik <- es$loglik
      
      # check
      diff <- abs(loglik - loglikPrev) / ( 1 + abs(loglik) )
      loglikPrev <- loglik
      it <- it + 1
      #
      crit <- ( diff > ctrl_em$tol & it < ctrl_em$iter )
    }
    ####################################################
  }
  
  # count parameters
  npar <- if ( H == 0 ) {
    Q*C + (choose(Q, 2) + Q)*C + Q*P*C + C-1
  } else {
    Q*K + (choose(Q, 2) + Q)*K + Q*P*K + R*H + (choose(R, 2) + R)*H + C-1
  }
  
  cl <- factor( mclust::map(z), labels = if(H > 0) c(classnames, 1:H) else classnames )
  out <- structure( list(learn = learn, K = K, H = H, parameters = parameters,
                         z = z, classification = cl,
                         loglik = loglik, N = N, npar = npar,
                         obs = obs, ext = ext) )
  attr(out, "info") <- list(it = it, control = ctrl_em[1:2])
  return(out)
}
