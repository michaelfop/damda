# damda - Dimension-Adaptive Discriminant Analysis

An [R](https://www.r-project.org/) package implementing *Dimension-Adaptive Mixture Discriminant Analysis* for supervised classification in settings where the test data might include unknown classes and/or being recorded on additional dimensions/variables. The package also implements a variable selection procedure for detecting the most informative variables for classification in high-dimensional situations.

## Installation

You can install the package directly from GitHub:

```{r}
# install.packages("devtools")
devtools::install_github("michaelfop/damda")
```

## Usage

Check the package vignette for some examples
```{r}
vignette("damda")
```

## References

Fop, M., Mattei, P. A., Bouveyron, C., Murphy, T. B. (2021). Unobserved classes and extra variables in high-dimensional discriminant analysis. *Advances in Data Analysis and Classification*, accepted. [ArXiv preprint](https://arxiv.org/abs/2102.01982).
