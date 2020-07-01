### tests for stochastic gradient ascent
### Test the following things:
###
### 1. basic 2-D SGA
###    SGA without function, only gradient
###    SGA neither function nor gradient
###    SGA in 1-D case
### 2. SGA w/momentum
### 3. SGA full batch
### 4. SGA, no gradient supplied
###    SGA, return numeric hessian, gradient provided
###    SGA, return numeric hessian, no gradient provided
###    SGA, printlevel 1, storeValues
###    SGA, NA as iterlim: should give informative error
###    SGA, storeValues but no fn (should fail)
###
### using highly unequally scaled data
###    SGA without gradient clipping (fails)
###    SGA with gradient clipping (works, although does not converge)

library(maxLik)
library(testthat)

## ---------- OLS 
## log-likelihood function(s):
## return log-likelihood on validation data
loglik <- function(beta, index) {
   e <- yValid - XValid %*% beta
   -crossprod(e)/length(y)
}
## gradlik: work on training data
gradlik <- function(beta, index) {
   e <- yTrain[index] - XTrain[index,,drop=FALSE] %*% beta
   g <- t(-2*t(XTrain[index,,drop=FALSE]) %*% e)
   -g/length(index)
}

### create random data
set.seed(1)
N <- 1000
x <- rnorm(N)
X <- cbind(1, x)
y <- 1 + x + rnorm(N)
## training-validation
iTrain <- sample(N, 0.8*N)
XTrain <- X[iTrain,,drop=FALSE]
XValid <- X[-iTrain,,drop=FALSE]
yTrain <- y[iTrain]
yValid <- y[-iTrain]
cat("Analytic solution (training data):\n")
start <- c(const=10, x=10)
b0 <- drop(solve(crossprod(XTrain)) %*% crossprod(XTrain, yTrain))
names(b0) <- names(start)
tol <- 1e-1  # coefficient tolerance

## ---------- 1. working example
res <- maxSGA(loglik, gradlik, start=start,
            control=list(printLevel=0, iterlim=200,
                         SG_batchSize=100, SG_learningRate=0.1,
                         storeValues=TRUE),
            nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)
                           # SGA usually ends with gradient not equal to 0 so we don't test that

## ---------- store parameters
res <- maxSGA(loglik, gradlik, start=start,
              control=list(printLevel=0, iterlim=20,
                           SG_batchSize=100, SG_learningRate=0.1,
                           storeParameters=TRUE),
              nObs=length(yTrain))
expect_equal(dim(storedParameters(res)), c(1 + nIter(res), 2))

## ---------- no function, only gradient
res <- maxSGA(grad=gradlik, start=start,
            control=list(printLevel=0, iterlim=10, SG_batchSize=100),
            nObs=length(yTrain))

## ---------- neither function nor gradient
try(res <- maxSGA(start=start,
                  control=list(printLevel=0, iterlim=10, SG_batchSize=100),
                  nObs=length(yTrain))
    )

## ---------- 1D case
N1 <- 1000
t <- rexp(N1, 2)
loglik1 <- function(theta, index) sum(log(theta) - theta*t[index])
gradlik1 <- function(theta, index) sum(1/theta - t[index])
res <- maxSGA(loglik1, gradlik1, start=1,
              control=list(iterlim=300, SG_batchSize=20), nObs=length(t))
expect_equal(coef(res), 1/mean(t), tolerance=0.2)
expect_null(hessian(res))

## ---------- 2. SGA with momentum
res <- maxSGA(loglik, gradlik, start=start,
            control=list(printLevel=0, iterlim=200,
                         SG_batchSize=100, SG_learningRate=0.1, SGA_momentum=0.9),
            nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)

## ---------- 3. full batch
res <- maxSGA(loglik, gradlik, start=start,
            control=list(printLevel=0, iterlim=200,
                         SG_batchSize=NULL, SG_learningRate=0.1),
            nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)

## ---------- 4. no gradient
res <- maxSGA(loglik, start=start,
              control=list(iterlim=1000, SG_learningRate=0.02), nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)

## ---------- return Hessian, gradient provided
res <- maxSGA(loglik, gradlik, start=start,
              control=list(iterlim=1000, SG_learningRate=0.02),
              nObs=length(yTrain),
              finalHessian=TRUE)
expect_equal(coef(res), b0, tolerance=tol)
expect_equal(dim(hessian(res)), c(2,2))

## ---------- return Hessian, no gradient
res <- maxSGA(loglik, start=start,
              control=list(iterlim=1000, SG_learningRate=0.02),
              nObs=length(yTrain),
              finalHessian=TRUE)
expect_equal(coef(res), b0, tolerance=tol)
expect_equal(dim(hessian(res)), c(2,2))

### ---------- SGA, printlevel 1, storeValues ----------
### it should just work
res <- maxSGA(loglik, gradlik, start=start,
              control=list(iterlim=2, storeValues=TRUE, printLevel=1),
              nObs=length(yTrain),
              finalHessian=TRUE)

### ---------- SGA, NA as iterlim ----------
### should give informative error
try(res <- maxSGA(loglik, gradlik, start=start,
                  control=list(iterlim=NA),
                  nObs=length(yTrain),
                  finalHessian=TRUE)
    )

### ---------- SGA, fn missing but storeValues=TRUE
### should give informative error
try(res <- maxSGA(grad=gradlik, start=start,
                  control=list(iterlim=10, storeValues=TRUE),
                  nObs=length(yTrain))
                  )

## ---------- gradient by observations
gradlikO <- function(beta, index) {
   e <- yTrain[index] - XTrain[index,,drop=FALSE] %*% beta
   g <- -2*drop(e)*XTrain[index,,drop=FALSE]
   -g/length(index)
}
res <- maxSGA(grad=gradlikO, start=start,
              control=list(printLevel=0, iterlim=100,
                           SG_batchSize=100),
              nObs=length(yTrain))
expect_equal(coef(res), b0, tolerance=tol)


### -------------------- create unequally scaled data
set.seed(1)
N <- 1000
x <- rnorm(N, sd=100)
XTrain <- cbind(1, x)
yTrain <- 1 + x + rnorm(N)
start <- c(const=10, x=10)

## ---------- no gradient clipping:
## should fail with informative "NA/Inf in gradient" message
try(res <- maxSGA(loglik, gradlik, start=start,
                  control=list(iterlim=100, SG_learningRate=0.5),
                  nObs=length(yTrain))
    )

## ---------- gradient clipping: should not fail
res <- maxSGA(loglik, gradlik, start=start,
              control=list(iterlim=100, SG_learningRate=0.5,
                           SG_clip=1e6),
              nObs=length(yTrain)
              )
