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
ll <- function(par) dnorm(x, mean=par[1], sd=par[2], log=TRUE)
for(method in c("NR", "BFGS", "BFGSR")) {
   cat("-- method", method, "--\n")
   m <- maxLik(ll, start=c(0, 2), method=method, control=list(iterlim=1))
   expect_equal(maxValue(m), -41.35, tolerance=0.01)
   expect_true(is.vector(gradient(m)), info="'gradient' returns a vector")
   expect_equal(length(gradient(m)), 2, info="'gradient(m)' is of length 2")
   expect_true(is.matrix(estfun(m)), info="'estfun' returns a matrix")
   expect_equal(dim(estfun(m)), c(20,2), info="'estfun(m)' is 20x2 matrix")
   cat("MaxControl structure:\n")
   show(maxControl(m))
}
