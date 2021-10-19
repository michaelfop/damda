\name{damda_varsel}
\alias{damda_varsel}

\title{
Inductive Variable Selection for Dimension-Adaptive Mixture Discriminant Analysis
}

\description{
The function implements a fast inductive variable selection approach for the Dimension-Adaptive Mixture Discriminant Analysis classifier. The method allows to find the optimal subset of variables having the most useful information at discriminating the classes in the test data. A greedy stepwise-forward search is implemented, which can also be run exploiting parallel computing functionalities.
}

\usage{
damda_varsel(learn, data,
             K = learn$K,
             H = 0:1,
             regularize = FALSE,
             start = NULL,
             control_em = damda::control_em(),
             control_reg = damda::control_reg(),
             parallel = FALSE,
             verbose = TRUE,
             maxit = 100)
}


\arguments{
  \item{learn}{
  A list containing a collection of class-specific parameters estimated in the training phase, a.k.a. the learning phase. The parameters typically are those corresponding to a Gaussian mixture discriminant analysis classifier. The list must include the following slots:
  \describe{
        \item{\code{pro}}{A vector containing the class mixing proportions (class proportions).}
        \item{\code{mu}}{The mean for each class, arranged column-wise, i.e. columns denote the classes.}
        \item{\code{sigma}}{An array containing the class-specific covariance matrices.}
        \item{\code{K}}{The number of classes observed in the training set.}
  }
}
  \item{data}{
  A matrix or data.frame containing the test data.
}
  \item{K}{
  The number of classes observed in the training data. No need to be specified if the list in argument \code{learn} already includes the number of classes in the training set.
}
  \item{H}{
  An integer vector specifying the numbers of extra classes for which the BIC is to be calculated. Default is to look from 0 to 1 extra classes in the test data.
}
  \item{regularize}{
  A logical argument indicating if Bayesian regularization should be performed. Default to \code{FALSE}.
}
  
  \item{start}{
  An optional vector containing the indexes of variables used to initialize the variable selection algorithm. If \code{NULL}, all the variables present in the learning phase are used as the initial set of variables useful for classification.
}
  
  \item{control_em}{
  A list of control parameters used in the EM algorithm for inductive model estimation; see also \code{\link{control_em}}.
}
  
  \item{control_reg}{
  A list of hyper parameters for Bayesian regularization. Only used when \code{regularization = TRUE}; see also \code{\link{control_reg}}.
}

  \item{parallel}{
 A logical argument indicating if parallel computation should be used for the variable selection stepwise search. If TRUE, all the available cores are used. The argument could also be set to a numeric integer value specifying the number of cores to be employed.
}

  \item{verbose}{
  If \code{TRUE} a progress bar will be shown.
}

 \item{maxit}{The maximum number of iterations (of addition and removal steps) the variable selection algorithm is allowed to run for.}
}


\details{
  The function implements an inductive variable selection procedure for the Dimension-Adaptive Mixture Discriminant Analysis (D-AMDA) classifier. A greedy forward-stepwise search is used to search the model space, where variables are added and removed in turn from the current set of variables useful for classification. 
  
To assess the discriminating power of a variable, the selection procedure computes at each iteration the BIC difference between a model where the variable is useful for classification against a model where the variable is uninformative or redundant. In addition to perform variable selection, the algorithm also returns the optimal number of hidden classes in the test data (if any).
}

\value{
A list including the following slots: 
\item{variables}{The set of selected relevant classification variables.}
\item{model}{An object of class \code{\link{damda}} containing the D-AMDA model fitted on the selected variables.}
\item{time}{The time taken to run the variable selection procedure.}
}

\references{
Fop, M., Mattei, P. A., Bouveyron, C., Murphy, T. B. (2021). Unobserved classes and extra variables in high-dimensional discriminant analysis. \emph{Advances in Data Analysis and Classification}, accepted.
}

