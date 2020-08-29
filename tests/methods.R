## Test methods.  Note: only test if methods work in terms of dim, length, etc,
## not in terms of values here
##
## ...
## * printing summary with max.columns, max.rows

library(maxLik)
require(testthat)
require(sandwich)
set.seed(0)

## Test standard methods for "lm"
x <- runif(20)
y <- x + rnorm(20)
m <- lm(y ~ x)
print(nObs(m))
print(stdEr(m))

## Test maxControl methods:
set.seed(9)
x <- rnorm(20, sd=2)
ll1 <- function(par) dnorm(x, mean=par, sd=1, log=TRUE)
ll2 <- function(par) dnorm(x, mean=par[1], sd=par[2], log=TRUE)
for(method in c("NR", "BFGS", "BFGSR")) {
   cat("-- method", method, "--\n")
   m <- maxLik(ll2, start=c(0, 2), method=method, control=list(iterlim=1))
   expect_equal(maxValue(m), -41.35, tolerance=0.01)
   expect_true(is.vector(gradient(m)), info="'gradient' returns a vector")
   expect_equal(length(gradient(m)), 2, info="'gradient(m)' is of length 2")
   expect_true(is.matrix(estfun(m)), info="'estfun' returns a matrix")
   expect_equal(dim(estfun(m)), c(20,2), info="'estfun(m)' is 20x2 matrix")
   cat("MaxControl structure:\n")
   show(maxControl(m))
}

## Test methods for non-likelihood optimization
hatf <- function(theta) exp(- theta %*% theta)
for(optimizer in c(maxNR, maxBFGSR, maxBFGS, maxNM, maxSANN, maxCG)) {
   name <- as.character(quote(optimizer))
   res <- optimizer(hatf, start=c(1,1))
   if(name %in% c("maxNR", "maxBFGS", "maxNM", "maxCG")) {
      expect_equal(coef(res), c(0,0), tol=1e-5,
                   info=paste0(name, ": result (0,0)"))
   }
   expect_equal(objectiveFn(res), hatf, info=paste0(name, ": objectiveFn correct"))
}

## Test maxLik vcov related methods
set.seed( 15 )
t <- rexp(20, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
a <- maxLik(loglik, start=1)
expect_equal(dim(vcov(a)), c(1,1), info="vcov 1D numeric correct")
expect_equal(length(stdEr(a)), 1, info="stdEr 1D numeric correct")
a <- maxLik(loglik, gradlik, hesslik, start=1)
expect_equal(dim(vcov(a)), c(1,1), info="vcov 1D analytic correct")
expect_equal(length(stdEr(a)), 1, info="stdEr 1D analytic correct")

### ---------- printing ----------
### ---------- max.columns, max.rows ----------
loglik <- function(beta) {
   e <- y - X %*% beta
   -crossprod(e)
}
gradlik <- function(beta) {
   e <- y - X %*% beta
   l <- crossprod(e)
   g <- t(-2*t(X) %*% e)
   -g
}
## linear regression with many columns
X <- matrix(rnorm(20*15), 20, 15)
beta <- rep(1, ncol(X))
y <- X %*% beta + rnorm(20, sd=0.3)
m <- maxNR(loglik, gradlik, start=rep(1, ncol(X)), iterlim=1)
print(summary(m, hessian=TRUE), max.rows=4, max.cols=2, digits=1)
                           # digits=1: very crude rounding but we want just to test the lines
