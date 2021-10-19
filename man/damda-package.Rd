\name{damda-package}
\alias{damda-package}
\docType{package}

\title{Dimension-Adaptive Mixture Discriminant Analysis}

\description{Dimension-Adaptive Mixture Discriminant Analysis with variable selection for classification problems where the test data includes unobserved classes and is recorded on additional variables compared to the training data.}

\details{A package implementing an adaptive mixture discriminant analysis classifier for the situation where the test data might be recorded on extra variables and include classes not observed during the training phase. The package implements a fast inductive estimation framework and allows to perform variable selection with forward-greedy search to detect the variables most useful to discriminate between the classes.}

\section{How to cite this package}{
To cite \pkg{damda} in publications use:
\preformatted{

@article{fop_2021,
  author = {Fop, M. and Mattei, P. A. and Bouveyron, C. and Murphy, T. B.},
  journal = {Advances in Data Analysis and Classification},
  number = {Accepted},
  pages = {},
  title = {Unobserved classes and extra variables in high-dimensional discriminant analysis},
  volume = {},
  year = {2021}
}

}}

\author{
Michael Fop.

Maintainer: Michael Fop \email{michael.fop@ucd.ie}
}

\references{
Fop, M., Mattei, P. A., Bouveyron, C., Murphy, T. B. (2021). Unobserved classes and extra variables in high-dimensional discriminant analysis. \emph{Advances in Data Analysis and Classification}, accepted.
}

\keyword{package}
