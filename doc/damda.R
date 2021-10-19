## ----setup, include=FALSE---------------------------------------------------------------
# setwd("~/Dropbox/Work/PhD_Work/MbC_VariableSelectionReview/stat_surveys/review/R_code")
library(knitr)
opts_chunk$set(fig.align="center",
               fig.width=5, fig.height=4.5,
               dev.args=list(pointsize=8))
options("width" = 90)

## ---------------------------------------------------------------------------------------
library(damda)

## ---- message=FALSE---------------------------------------------------------------------
library(MBCbook)
data(wine27)

data <- scale(wine27[,1:27])
class <- wine27$Type

# set training and test set
set.seed(987321)

R <- ncol(data)
obs <- sort( sample(1:R, R*1/3) )
ext <- setdiff(1:R, obs)

n <- nrow(data)
n_test <- n/2    # test set larger than training set

sel_test <- sample(1:n, n_test)
data_test <- data[sel_test, ]
class_test <- class[sel_test]

rem <- which(class == "Barolo")
data_train <- data[-union(rem, sel_test), obs]
class_train <- factor( class[-union(rem, sel_test)] )

## ---- message=FALSE---------------------------------------------------------------------
library(mclust)
mc <- MclustDA(data_train, class_train, modelType = "EDDA", modelNames = "VVI")

# predict using only variables available in training set
mc_class <- predict.MclustDA(mc, data_test[,obs])$class  
table(class_test, mc_class)

## ---- message=FALSE---------------------------------------------------------------------
# learning phase

P <- ncol(data_train)
K <- nlevels(class_train) 

learn <- list(pro = mc$prop)
tmp <- lapply(mc$models, "[[", "parameters")
learn$mu <- sapply(tmp, "[[", "mean")
names <- colnames(data_train)
rownames(learn$mu) <- colnames(data_train)
learn$sigma <- array( unlist( lapply(tmp, function(x) x$variance$sigma) ), c(P,P, K),
                      dimnames = list(names, names) )
learn$K <- length(learn$pro)

## ---------------------------------------------------------------------------------------
fit <- damda(learn, data_test, H = 0:2, regularize = TRUE, verbose = FALSE)
fit$BIC
fit$H

table(class_test, fit$classification)   # nice!

## ----message=FALSE----------------------------------------------------------------------
require(clusterGeneration)

C <- 3
K <- C - 1

set.seed(123098)
gen_rand_clust <- genRandomClust(
  numClust = C,
  sepVal = 0.1,
  rangeN = c(1000, 2000),
  numNonNoisy = 4,
  numNoisy = 12,
  numReplicate = 1,
  outputDatFlag = FALSE, outputLogFlag = FALSE)

data <- gen_rand_clust$datList$test_1
class <- factor(gen_rand_clust$memList$test_1, labels = LETTERS[1:3])

## ---------------------------------------------------------------------------------------
data <- scale(data)

R <- ncol(data)
obs <- sort( sample(1:R, R/2) )
ext <- setdiff(1:R, obs)
P <- length(obs)

n <- nrow(data)
n_train <- 500
n_test <- 300

sel_test <- sample(1:n, n_test)
data_test <- data[sel_test, ]
class_test <- class[sel_test]

rem <- which(class == "C")
data_tmp <- data[-union(rem, sel_test), obs]
class_tmp <- factor( class[-union(rem, sel_test)] )

sel_train <- sample(1:nrow(data_tmp), n_train)
data_train <- data_tmp[sel_train, ]
class_train <- class_tmp[sel_train]

## ---- message=FALSE---------------------------------------------------------------------
mc <- MclustDA(data_train, class_train, modelType = "EDDA")
mc$models[[1]]$modelName      # selected covariance model

# predict using only variables available in training set
mc_class <- predict.MclustDA(mc, data_test[,obs])$class  
table(class_test, mc_class)


# learning phase
learn <- list(pro = mc$prop)
tmp <- lapply(mc$models, "[[", "parameters")
learn$mu <- sapply(tmp, "[[", "mean")
names <- colnames(data_train)
rownames(learn$mu) <- colnames(data_train)
learn$sigma <- array( unlist( lapply(tmp, function(x) x$variance$sigma) ), c(P,P, K),
                      dimnames = list(names, names) )
learn$K <- length(learn$pro)

## ---------------------------------------------------------------------------------------
fit <- damda(learn, data_test, H = 0:2, verbose = FALSE)
fit$H

table(class_test, fit$classification)

## ----message=FALSE----------------------------------------------------------------------
vsel <- damda_varsel(learn, data_test, H = 0:2)

vsel$variables             # inspect selected variables
gen_rand_clust$noisyList   # these are the "bad" variables

table(class_test, vsel$model$classification)   # nice!

## ---------------------------------------------------------------------------------------
sessionInfo()

