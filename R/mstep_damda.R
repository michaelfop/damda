#
#======================= M-step functions
# 
# Copyright (C) 2021 - Michael Fop


mstep_damda <- function(data, ext, z, K, H, proX, parameters, regularize, regpar)
  # M step for the EM algorithm used to fit the damda model
{
  data <- as.matrix(data)
  R <- ncol(data)
  Ny <- nrow(data)
  C <- K + H
  obs <- if ( !is.logical(ext) ) setdiff(1:R, ext) else 1:R
  Nc <- colSums(z)
  
  if ( regularize ) {
    den <- sum( Nc + regpar$alpha )
    z <- z*(Ny/den) + (regpar$alpha/den)
    Nc <- colSums(z)
  }
  
  temp <- mclust::covw(data, z, normalize = FALSE)
  if ( regularize ) {
    for ( j in 1:C ) {
      temp$W[,,j] <- temp$W[,,j] + regpar$scale*( (1/C)^(1/R) )
      temp$S[,,j] <- temp$W[,,j]/Nc[j]
    }
  }
  
  if ( H == 0 ) {
    # no new classes --- yes new variables
    
    # pro
    parameters$pro <- Nc/Ny
    if ( any( Nc <= R & !regularize) ) stop("Empty component")
    
    # mu and covariance
    out <- ice(data, z, parameters$mu, parameters$sigma, temp$mean, temp$W, obs-1, ext-1, Nc, C)
    parameters$mu[ext,] <- out$mu
    parameters$sigma[obs,ext, ] <- out$cross
    parameters$sigma[ext,obs, ] <- aperm(out$cross, c(2,1,3))
    parameters$sigma[ext,ext, ] <- out$sigma
    
  } else {
    # yes new classes
    
    # # pro
    # # update by fixing the learning phase probabilities -- NOT IMPLEMENTED
    # Nh <- colSums( z[,(K+1):C, drop = FALSE] )
    # parameters$pro[(K+1):C] <- Nh/Ny
    # parameters$pro[1:K] <- ( 1 - sum(parameters$pro[(K+1):C]) ) * proX
    #
    # or simply compute them
    parameters$pro <- colMeans(z)
    #
    if ( any( parameters$pro*Ny <= R ) & !regularize ) stop("Empty component")
    
    # mu and covariance
    if ( is.logical(ext) ){
      # no new variables
      #
      eig <- apply( temp$S[,,(K+1):C, drop = FALSE], 3, function(s) eigen(s, only.values = TRUE) )
      if ( any( sapply(eig, "[[", "values") < sqrt(.Machine$double.eps) ) ) { 
        stop("Covariance not positive-definite") }
      parameters$mu[,(K+1):C] <- temp$mean[,(K+1):C, drop = FALSE]
      parameters$sigma[,,(K+1):C] <- temp$S[,,(K+1):C, drop = FALSE]
      
    } else {
      # yes new variables
      #
      # unknown classes, all variables
      eig <- apply( temp$S, 3, function(s) eigen(s, only.values = TRUE) )
      if ( any(sapply(eig, "[[", "values") < sqrt(.Machine$double.eps)) ) { 
        stop("Covariance not positive-definite") }
      parameters$mu[,(K+1):C] <- temp$mean[,(K+1):C, drop = FALSE]
      parameters$sigma[,,(K+1):C] <- temp$S[,,(K+1):C, drop = FALSE]
      #
      # known classes, extra variables
      out <- ice(data, z, parameters$mu[,1:K,drop = FALSE], parameters$sigma[,,1:K,drop = FALSE], 
                 temp$mean[,1:K, drop = FALSE], temp$W[,,1:K, drop = FALSE], 
                 obs-1, ext-1, Nc[1:K], K)
      parameters$mu[ext,1:K] <- out$mu
      parameters$sigma[obs,ext,1:K] <- out$cross
      parameters$sigma[ext,obs,1:K] <- aperm(out$cross, c(2,1,3))
      parameters$sigma[ext,ext,1:K] <- out$sigma
    }
  }
  
  return(parameters)
}
