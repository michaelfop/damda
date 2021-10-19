#
#======================= Utility functions
#
# Copyright (C) 2021 - Michael Fop


bhat_dist <- function(mu1, sigma1, mu2, sigma2)
  # Compute the Bhattacharyya distance
{
  if ( dim(mu1)[1] != dim(mu2)[1] ) stop("mean parameters must be column vectors with the same number of rows (variables)")
  if ( min( c(length(dim(sigma1)), length(dim(sigma1))) ) != 3 ) stop("sigma parameters must be arrays")
  K <- ncol(mu1)
  C <- ncol(mu2)
  # 1 --> parameters
  # 2 --> clusters
  
  mat <- matrix(NA, K, C)   # bhat distance between k observed class and c cluster
  for ( l in 1:K ) {
    s1 <- as.matrix( sigma1[,,l] )
    for ( c in 1:C ) {
      s2 <- as.matrix( sigma2[,,c] )
      pool <- (s1 + s2)*0.5
      inv <- solve(pool)
      mat[l,c] <- 0.125 * mahalanobis(mu1[,l], mu2[,c], inv, inverted = TRUE) +
        0.5*( determinant(pool, TRUE)$mod - 0.5 * ( determinant(s1, TRUE)$mod + determinant(s2, TRUE)$mod ) )
    }
  }
  dimnames(mat) <- list(1:K, 1:C)
  return(mat)
}


control_em <- function(tol = 1e-05, iter = 1e03, printmsg = FALSE, init = c("mclust", "hc"))
  # Set control parameters for inductive EM
{
  init <- match.arg(init, c("mclust", "hc"))
  list(tol = tol, iter = iter, printmsg = printmsg, init = init)
}



control_reg <- function(gamma = 0.01, alpha = 0, ...)
  # Set hyperparameters for regularization
{
  tmp <- list(...)
  if ( exists("data", tmp) ) {
    N <- nrow(tmp$data)
    S <- cov(tmp$data)*(N-1)/N
    R <- ncol(S)
    if (N <= R) S <- diag(diag(S))
  } else {
    S <- diag(gamma, 1)
    K <- 0
    R <- 1
    N <- 1
  }
  dt <- det(S)
  if ( dt < sqrt(.Machine$double.eps) ) dt <- sqrt( .Machine$double.eps )/10
  scale <- S/(dt^(1/R)) * gamma^(1/R) * (1/N)
  return( list(gamma = gamma, alpha = alpha, scale = scale, K = K) )
}
