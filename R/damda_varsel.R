#
#======================= Variable selection for Dimension-Adaptive Mixture Discriminant Analysis
#
# Copyright (C) 2021 - Michael Fop

bic_reg <- function(x, y)
  # Compute BIC for the regression of y on x,
  # after doing subset selection  on the predictors.
{
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- length(y)
  mod <- try( BMA::bicreg(y = y, x = x, nbest = 1), silent = TRUE )
  if ( class(mod) == "try-error" ) return( list(bic = NA) )
  subset <- which(mod$which[1,])
  mod <- lm.fit( y = y, x = cbind(1, x[,subset]) )
  
  # calculate the BIC for the regression
  sigma <- sqrt( mean(residuals(mod)^2) )
  p <- n - df.residual(mod) + 1
  bic <- list(bic = -n*log(2*pi) -2*n*log(sigma) -n -log(n)*p,
              subset = subset)
}


damda_varsel <- function( learn, data,
                          K = learn$K,
                          H = 0:1,
                          regularize = FALSE,
                          start = NULL,
                          control_em = damda::control_em(),
                          control_reg = damda::control_reg(),
                          parallel = FALSE,
                          verbose = TRUE,
                          maxit = 100)
  # Variable selection for damda using a greedy forward search
{
  if ( !all( exists("pro", learn), exists("mu", learn), exists("sigma", learn),  exists("K", learn) ) ) {
    stop("Object 'learn' should be a list with entries 'pro', 'mu', 'sigma', and 'K'")
  }
  K <- learn$K
  if ( K != length(learn$pro) ) stop("Number of classes 'K' does not match parameter dimensions")
  Xnames <- colnames(learn$sigma[,,1])
  varnames <- colnames(data)
  if ( is.null(Xnames) | is.null(varnames) ) {
    warning("Variables from the learning phase and/or variables in the test data should have names")
  }
  if ( is.null(Xnames) ) Xnames <- paste0("V", 1:ncol(learn$sigma[,,1]))
  if ( is.null(varnames) ) varnames <- paste0("V", 1:ncol(data))
  obs <- match( Xnames, varnames ) 
  ext <- match( setdiff(varnames, Xnames), varnames )   # additional variables
  if ( length(ext) == 0 ) ext <- FALSE
  # if ( max(remclass) > K - 1) stop("Cannot remove all the classes observed in the learning phase")
  classnames <- colnames(learn$mu)
  if ( is.null(classnames) ) classnames <- as.character(1:K)
  colnames(learn$mu) <- names(learn$pro) <- classnames
  
  i <- NULL # avoid 'no visible binding for global variable'
  
  Obs <- varnames[obs]        # observed variables
  Ext <- varnames[ext]        # additional variables
  if ( length(ext) == 0 ) ext <- FALSE
  
  regpar <- if ( regularize ) {
    control_reg(gamma = control_reg$gamma, alpha = control_reg$alpha, data = data, K = K)
  } else NULL 
  
  
  R <- length(varnames)
  P <- length(Xnames)
  Q <- R - P
  N <- nrow(data)
  
  # parallel-------------------------------------------------------
  parallel <- if ( is.logical(parallel) ) {
    if(parallel) start_parallel(parallel) else FALSE
  } else { start_parallel(parallel) }
  # on.exit( if(parallel) parallel::stopCluster(attr(parallel, "cluster")) )
  `%DO%` <- if(parallel) `%dopar%` else `%do%`
  #----------------------------------------------------------------
  
  # learning phase-------------------------------------------------
  learn0 <- learn
  # ---------------------------------------------------------------
  
  # initialize clustering variables-----------------------------------------
  if ( !is.null(start) ) {
    # start from a given set of variables from the learning phase
    tmp <- learn
    noStart <- setdiff(Xnames, start)
    j <- match( noStart, Xnames )
    tmp$mu <- tmp$mu[-j,, drop = FALSE]
    tmp$sigma <- tmp$sigma[-j,-j,, drop = FALSE]
    clus <- start
    learn <- tmp
  } else {
    # start from all the variables of the learning set
    clus <- Xnames
  }
  noclus <- setdiff(varnames, clus)
  
  clusMod <- damda(learn, data, H = H, control_em = control_em, control_reg = control_reg, regularize = regularize, verbose = FALSE)
  bicClusMod <- clusMod$bic
  if ( verbose ) cat(" Starting clustering set: \n", clus, "\n")
  #-------------------------------------------------------------------------
  
  # variable selection------------------------------------------------------
  CRIT <- TRUE
  toRem <- toAdd <- NULL
  zero <- (.Machine$double.eps)^(1/3)
  it <- 0
  
  tm <- Sys.time()
  while ( CRIT ) {
    # add step.......................................................
    #
    # do not add the variable just removed
    Noclus <- if ( !is.null(toRem) ) {
      if ( is.na(match(toRem, noclus)) ) noclus else noclus[ -match(toRem, noclus) ]
    } else noclus
    #
    if ( length(Noclus) != 0 ) {            # there are no variables to add
      out <- foreach( i = Noclus, .packages = c("mclust", "damda", "BMA"),
                      .combine = "rbind", .errorhandling = "pass",
                      .final = function(o) matrix(o, ncol = 2) ) %DO% {
                        set <- intersect( varnames, c(clus, i) )           # clustering set after adding variable i
                        bicReg <- bic_reg( y = data[,i], x = data[,clus] )
                        #
                        if ( i %in% Obs ) {
                          # add a variable from the observed set
                          tmp <- learn
                          vn <- intersect( colnames(learn0$sigma[,,1]), c(i, colnames(tmp$sigma[,,1])) )
                          tmp$mu <- learn0$mu[vn,, drop = FALSE]
                          tmp$sigma <- learn0$sigma[vn,vn,, drop = FALSE]
                          
                        } else {
                          # add a variable from the set of new ones
                          tmp <- learn
                        }
                        #
                        modClus <- damda(tmp, data[,set], H = H, regularize = regularize,
                                         control_em = control_em, control_reg = control_reg, verbose = FALSE)
                        return( c(bicReg$bic, modClus$bic) )
                      }
      rownames(out) <- Noclus
      bicNoClus <- bicClusMod + out[,1]
      bicAdd <- out[,2]
      bicDiff <- bicAdd - bicNoClus
      
      # no adaptive DA model was fitted
      if ( all( is.na(bicDiff) ) ) bicDiff <- -Inf
      
      # variable proposed for adding
      best <- which.max(bicDiff)
      toAdd <- names(bicDiff)[best]
      
      # add
      if ( verbose ) cat(" add ", toAdd)
      if ( bicDiff[best] > zero ) {
        added <- TRUE
        bicClusMod <- bicAdd[best]
        j <- match(toAdd, noclus)
        clus <- c(clus, toAdd)
        noclus <- noclus[-j]
        if ( verbose ) cat(" : accepted -- BIC difference =", bicDiff[best], "\n")
        #
        # update DA model if added an observed variable
        if ( toAdd %in% Obs ) {
          vn <- intersect( colnames(learn0$sigma[,,1]), c(toAdd, colnames(learn$sigma[,,1])) )
          learn$mu <- learn0$mu[vn,, drop = FALSE]
          learn$sigma <- learn0$sigma[vn,vn,, drop = FALSE]
        }
      } else {
        added <- FALSE
        toAdd <- NULL
        if ( verbose ) cat(" : rejected -- BIC difference =", bicDiff[best], "\n")
      }
      #
    } else { # if ( length(Noclus) != 0 )
      added <- FALSE
      toAdd <- NULL
    }
    #................................................................
    
    # removal step...................................................
    #
    # do not remove the variable just added
    Clus <- if ( !is.null(toAdd) ) {
      if ( is.na(match(toAdd, clus)) ) clus else clus[ -match(toAdd, clus) ]
    } else  clus
    #
    # do not remove the last variable from learning phase
    inlearn <- intersect(Xnames, Clus)
    #
    # need at least 2 variables from the learning phase
    if ( length(inlearn) < 3 ) break 
    
    out <- foreach( i = Clus, .packages = c("mclust", "damda", "BMA"),
                    .combine = "rbind", .errorhandling = "remove",
                    .final = function(o) matrix(o, ncol = 2) ) %DO% {
                      set <- setdiff(clus, i)                     # clustering set after removing variable i
                      bicReg <- bic_reg( y = data[,i], x = data[,set] )
                      #
                      # remove a variable from the observed set
                      if ( i %in% Obs ) {
                        tmp <- learn
                        vn <- colnames(tmp$sigma[,,1])
                        j <- match( i, vn )
                        tmp$mu <- tmp$mu[-j,, drop = FALSE]
                        tmp$sigma <- tmp$sigma[-j,-j,, drop = FALSE]
                        
                      } else {
                        # remove a variable from the set of new ones
                        tmp <- learn
                      }
                      #
                      modNoClus <- damda(tmp, data[,set], H = H, regularize = regularize,
                                         control_em = control_em, control_reg = control_reg, verbose = FALSE)
                      return( c(bicReg$bic, modNoClus$bic) )
                    }
    rownames(out) <- Clus
    bicRem <- rowSums(out)
    bicDiff <- bicRem - bicClusMod
    
    # no adaptive DA model was fitted
    if ( all( is.na(bicDiff) ) ) bicDiff <- -Inf
    
    # variable proposed for removal
    best <- which.max(bicDiff)
    toRem <- names(bicDiff)[best]
    
    # remove
    if ( verbose ) cat(" remove ", toRem)
    if ( bicDiff[best] > zero ) {
      if ( length(intersect(clus, Xnames)) == 1 ) {
        removed <- FALSE
        toRem <- NULL
        if ( verbose ) cat(": no, it's the last one observed in the learning phase!", "\n")
      } else {
        removed <- TRUE
        j <- match(toRem, clus)
        clus <- clus[-j]
        noclus <- c(noclus, toRem)
        bicClusMod <- out[best,2]
        if ( verbose ) cat(" : accepted -- BIC difference =", bicDiff[best], "\n")
        #
        # update DA model if removed an observed variable
        if ( toRem %in% Obs ) {
          vn <- colnames(learn$sigma[,,1])
          j <- match( toRem, vn )
          # learn$data <- learn$data[, setdiff(vn, toRem), drop = FALSE ]
          learn$mu <- learn$mu[-j,, drop = FALSE]
          learn$sigma <- learn$sigma[-j,-j,, drop = FALSE]
        }
      }
      #
    } else {
      removed <- FALSE
      toRem <- NULL
      if ( verbose ) cat(" : rejected, BIC difference =", bicDiff[best], "\n")
      if ( length(clus) == R ) break      # all variables are selected
    }
    #................................................................
    
    it <- it + 1
    
    # check
    CRIT <- ( (removed | added) & (length(clus) >= 2) & (it < maxit) )
  }
  #-------------------------------------------------------------------------
  
  # stop parallel computing
  if ( parallel ) parallel::stopCluster( attr(parallel, "cluster") )
  
  # fit model on selected variables
  mod <- damda(learn, data[,intersect(varnames, clus)], H = H, regularize = regularize,
               control_em = control_em, control_reg = control_reg, verbose = FALSE)
  
  df <- difftime( Sys.time(), tm, units = "secs" )
  
  return( list(variables = intersect(varnames, clus), model = mod, time = as.numeric(df)) )
}
