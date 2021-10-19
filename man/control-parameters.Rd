\name{control-parameters}
\alias{control_em}
\alias{control_reg}

\title{
Set Control Parameters
}
\description{
Set control parameters for the EM algorithm for inductive estimation and for Bayesian regularization of the class-specific covariance matrices
}
\usage{
control_em(tol = 1e-05, 
           iter = 1e03, 
           printmsg = FALSE, 
           init = c("mclust", "hc"))

control_reg(gamma = 0.01, 
            alpha = 0, ...)
}

\arguments{
  \item{tol}{Tolerance value for judging when convergence of the EM algorithm has been reached.}
  \item{iter}{Maximum number of iterations of the EM algorithm.}
  \item{printmsg}{Logical argument to print any warning/error messages during the estimation procedure. By default messages are suppressed (and typically user would use default setting).}
  \item{init}{Initialization of the EM algorithm. If \code{init = "mclust"}, the initial classification of the observations in the test data is obtained using model based clustering; if \code{init = "hc"}, agglomerative model-based hierarchical clustering is employed to obtain the initial classifiation. Default is \code{init = "mclust"}.}
  \item{gamma}{A scaling factor for regularization of the covariance matrices.}
  \item{alpha}{A scaling factor for regularization on the expected number of observations in each class (note: not tested).}
  \item{...}{Further arguments not meant to be modified by the user.}
}

\details{
Function \code{control_em} is employed to set standard parameters of the EM algorithm for inductive estimation.

Function \code{control_reg} is used to set hyperparameters for Bayesian regularization of the class covariance matrices. The hyperparameter \code{gamma} controls the amount of regularization, such that the determinant of the regularizing matrix is proportional to \eqn{\gamma/C}{gamma/C}, where C is the total number of classes present in the test data (those observed in the training set plus those additional classes which might be present in the test data). It is suggested to set \code{gamma} to \eqn{(log R)/N^2}{(log R)/N^2}, where R is the total number of variables observed in the test data, and N is the number of observations in the test data. The hyperparameter \code{alpha} controls the (minimum) expected number of observations in each class. Use of this hyperparameter has not been thoroughly tested and default is \code{alpha = 0}, which corresponds to no regularization on the expected number of observations in the classes. 
}

\value{
A list of parameters values.
}

\references{
Fop, M., Mattei, P. A., Bouveyron, C., Murphy, T. B. (2021). Unobserved classes and extra variables in high-dimensional discriminant analysis. \emph{Advances in Data Analysis and Classification}, accepted.}