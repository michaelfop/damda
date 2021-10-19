#
#======================= Initialize EM algorithm
#
# Copyright (C) 2021 - Michael Fop


init_em <- function(learn, parameters, data, H, hc_start, init)
  # Initialization for the EM algorithm
{
  Xnames <- colnames(learn$sigma[,,1])
  varnames <- colnames(data)
  obs <- match( Xnames, varnames ) 
  ext <- match( setdiff(varnames, Xnames), varnames )   # additional variables
  if ( length(ext) == 0 ) ext <- FALSE
  R <- length(varnames)
  P <- length(Xnames)
  Q <- R - P
  Ny <- nrow(data)
  K <- length(learn$pro)
  C <- H + K
  
  # hierarchical clustering----------------------------------------
  if ( init == "hc") {
    HC <- mclust::unmap( mclust::hclass(hc_start, G = C) )
    init <- mclust::mstep(modelName = "VVV", data, HC, prior = mclust::priorControl())
    #
    # impossible to compute M step with VVV
    tmp <- apply( init$parameters$variance$sigma, 3, function(s) eigen(s, only.values = TRUE)$values )
    if ( any( tmp < sqrt(.Machine$double.eps) ) || ( attr(init, "returnCode") < 0 ) ) {
      init <- mclust::Mclust(data, G = C, modelNames = "VVI")
    }
  } else {
    # initialization using mclust
    init <- try( suppressWarnings( mclust::Mclust(data, G = C, verbose = FALSE, initialization = list(hcPairs = hc_start),
                                                  modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "EVV", "VVV"),
                                                  prior = mclust::priorControl()) ), silent = TRUE )

    if ( class(init) == "try-error" ) {
      init <- suppressWarnings( mclust::Mclust(data, G = C, verbose = FALSE, initialization = list(hcPairs = hc_start),
                                               modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE"),
                                               prior = mclust::priorControl()) )
    }
  }
  # ---------------------------------------------------------------
  
  # matching-------------------------------------------------------
  # bhattcharyya distance between k observed class and c cluster
  kl <- bhat_dist(mu1 = parameters$mu[obs,1:K, drop = FALSE], sigma1 = parameters$sigma[obs,obs,1:K, drop = FALSE],
                  mu2 = init$parameters$mean[obs,,drop = FALSE], sigma2 = init$parameters$variance$sigma[obs,obs,, drop = FALSE])
  matching <- rep(NA, K)   # class k match to which cluster?
  names(matching) <- 1:K
  dimnames(kl) <- list(1:K, 1:C)
  KL <- kl
  for ( l in 1:K ) {
    val <- which(KL == min(KL), arr.ind = TRUE)
    i <- as.numeric( rownames( KL[val[1],val[2], drop = FALSE] ) )
    matching[i] <- as.numeric( colnames(KL[val[1],val[2], drop = FALSE]) )
    KL <- KL[-val[1], -val[2], drop = FALSE]
  }
  new <- setdiff(1:C, matching)
  # ---------------------------------------------------------------
  
  # test parameters initialization---------------------------------
  if ( H > 0 ) {
    #
    # yes new classes
    if ( is.logical(ext) ) {
      # no new variables
      #
      parameters$mu[,(K+1):C] <- init$parameters$mean[,new, drop = FALSE]
      parameters$sigma[,,(K+1):C] <- init$parameters$variance$sigma[,,new, drop = FALSE]
      parameters$pro[(K+1):C] <- init$parameters$pro[new]
      parameters$pro <- parameters$pro / sum(parameters$pro)
    }  else {
      # yes new variables
      #
      # new classes
      parameters$mu[,(K+1):C] <- init$parameters$mean[,new, drop = FALSE]
      parameters$sigma[,,(K+1):C] <- init$parameters$variance$sigma[,,new, drop = FALSE]
      parameters$pro[(K+1):C] <- init$parameters$pro[new]
      #
      # known classes
      parameters$mu[ext,1:K] <- init$parameters$mean[ext,matching]
      parameters$pro[1:K] <- init$parameters$pro[matching]
      parameters$pro <- parameters$pro / sum(parameters$pro)
      for ( l in 1:K ) {
        j <- matching[l]
        parameters$sigma[ext,ext, l] <- init$parameters$variance$sigma[ext,ext, j, drop = FALSE]
        parameters$sigma[obs,ext, l] <- init$parameters$variance$sigma[obs,ext, j, drop = FALSE]
        parameters$sigma[ext,obs, l] <- init$parameters$variance$sigma[ext,obs, j, drop = FALSE]
        tmp <- eigen(parameters$sigma[,,l], only.values = TRUE)$val
        
        if ( any( tmp < sqrt(.Machine$double.eps) ) ) {                 # ensure positive definite
          #
          parameters$sigma[obs,ext, l] <- 0
          parameters$sigma[ext,obs, l] <- 0
          #
          # # NOT TESTED
          # corr <- init$parameters$variance$sigma[,,j]*parameters$pro[l]*Ny
          # inv <- solve(parameters$sigma[obs,obs, l])
          # A <- crossprod(inv, corr[obs, obs])
          # B <- crossprod(corr[obs,ext, drop = FALSE], inv)
          # CC <- solve( tcrossprod(A, inv) )%*% t(B)
          # D <- crossprod(CC, inv)
          # E <- 1/(parameters$pro[j]*Ny) * ( tcrossprod(tcrossprod(D, corr[obs,obs]), D) - 2*B %*% CC + corr[ext,ext] )
          # parameters$sigma[ext,ext, l] <- E + D %*% CC
          # parameters$sigma[obs,ext, l] <- CC
          # parameters$sigma[ext,obs, l] <- t(C)
        }
      }
    }
  } else {
    # no new classes, yes new variables
    #
    parameters$mu[ext,1:C] <- init$parameters$mean[ext,matching]
    parameters$pro[1:C] <- init$parameters$pro[matching]
    for ( l in 1:K ) {
      j <- matching[l]
      parameters$sigma[ext,ext, l] <- init$parameters$variance$sigma[ext,ext, j, drop = FALSE]
      parameters$sigma[obs,ext, l] <- init$parameters$variance$sigma[obs,ext, j, drop = FALSE]
      parameters$sigma[ext,obs, l] <- init$parameters$variance$sigma[ext,obs, j, drop = FALSE]
      tmp <- eigen(parameters$sigma[,,l], only.values = TRUE)$val
      
      if ( any( tmp < sqrt(.Machine$double.eps) ) ) {                 # ensure positive definite
        #
        parameters$sigma[obs,ext, l] <- 0
        parameters$sigma[ext,obs, l] <- 0
        # 
        # # NOT TESTED
        # corr <- init$parameters$variance$sigma[,,j]*parameters$pro[l]*Ny
        # inv <- solve(parameters$sigma[obs,obs, l])
        # A <- crossprod(inv, corr[obs, obs])
        # B <- crossprod(corr[obs,ext, drop = FALSE], inv)
        # CC <- solve( tcrossprod(A, inv) )%*% t(B)
        # D <- crossprod(CC, inv)
        # E <- 1/(parameters$pro[j]*Ny) * ( tcrossprod(tcrossprod(D, corr[obs,obs]), D) - 2*B %*% CC + corr[ext,ext] )
        # parameters$sigma[ext,ext, l] <- E + D %*% CC
        # parameters$sigma[obs,ext, l] <- CC
        # parameters$sigma[ext,obs, l] <- t(CC)
      }
    }
  }
  
  dimnames(parameters$sigma) <- list(varnames, varnames)
  rownames(parameters$mu) <- varnames
  colnames(parameters$mu) <- 1:C
  
  attr(parameters, "kl") <- kl
  attr(parameters, "matching") <- matching
  return(parameters)
}

