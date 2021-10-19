\name{damda}
\alias{damda}
\alias{predict.damda}

\title{
Dimension-Adaptive Mixture Discriminant Analysis
}

\description{
Implements the Dimension-Adaptive Mixture Discriminant Analysis classifier for settings where the test data might include additional unknown classes and/or extra dimensions. The function performs automatic selection of the number of hidden classes using BIC.

}

\usage{
damda(learn, data, 
      K = learn$K, 
      H = 0:2, 
      regularize = FALSE, 
      control_em = damda::control_em(), 
      control_reg = damda::control_reg(), 
      verbose = TRUE)
      
\method{predict}{damda}(object, newdata, ...)
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
  An integer vector specifying the numbers of extra classes for which the BIC is to be calculated. Default is to look from 0 to 2 extra classes in the test data.
}
  \item{regularize}{
  A logical argument indicating if Bayesian regularization should be performed. Default to \code{FALSE}.
}
  \item{control_em}{
  A list of control parameters used in the EM algorithm for inductive model estimation; see also \code{\link{control_em}}.
}
  \item{control_reg}{
  A list of hyper parameters for Bayesian regularization. Only used when \code{regularization = TRUE}; see also \code{\link{control_reg}}.
}
  \item{verbose}{
  If \code{TRUE} a progress bar will be shown.
}

 \item{object}{An object of class \code{'damda'} resulting from a call to function \code{\link{damda}}.
}

 \item{newdata}{A data frame or matrix giving the data for which predictions need to be obtained. If missing, the data employed in the call to \code{\link{data}} are classified.
}

 \item{\dots}{Further arguments passed to or from other methods.
}

}

\details{
  The function implements Dimension-Adaptive Mixture Discriminant Analysis (D-AMDA) classifier for supervised learning in settings where the test data might include unknown classes and/or additional dimensions/variables. The model estimation procedure is based on an EM algorithm embedded in an inductive estimation framework.
  
  Note that only the test data are required in input, and no training data need to be provided. Indeed, the function requires the parameter estimates obtained during the training stage, and these must correspond to class-specific proportions, means, and covariance matrices. The training stage can potentially be performed using any type of classifier, as long as corresponding class-related parameters are provide. More in line with the proposed framework of D-AMDA, a mixture discriminant analysis classifier would be implemented in the training phase.
  
  Model selection in the context of adaptive classification corresponds to detection of hidden classes (if any), not previously observed during the training phase. To this purpose, the BIC is employed to select the optimal model.
}

\value{
An object of class \code{'damda'} containing the optimal D-AMDA classifier. The object in output is a list containing:
\item{learn}{The parameters learned during the training phase.}
\item{K}{Number of classes observed in the training data.}
\item{H}{Selected optimal number of hidden classes according to BIC.}
\item{parameters}{A list including the parameters of the training phase provided in input and those estimated during the discovery phase.}
\item{z}{A matrix whose \code{[i,k]}th entry is the probability that observation \code{i} of the test data belongs to the \code{k}th class.} 
\item{classification}{Predicted classification of the observations in the test set, corresponding to the maximum a posteriori of matrix \code{z}.}
\item{loglik}{Value of the maximized log-likelihood.}
\item{N}{Number of observations in the test data.}
\item{npar}{Number of estimated parameters.}
\item{obs}{A vector containing the indexes of the variables in the test set observed also in the training set.}
\item{ext}{A vector containing the indexes of the additional variables present in the test set but not observed in the training set.}
\item{bic}{Optimal BIC value.}
\item{BIC}{All BIC values.}
\item{data}{The test data matrix provided in input.}
}

\references{
Fop, M., Mattei, P. A., Bouveyron, C., Murphy, T. B. (2021). Unobserved classes and extra variables in high-dimensional discriminant analysis. \emph{Advances in Data Analysis and Classification}, accepted.
}

