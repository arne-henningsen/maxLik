library(maxLik)
require(testthat)
set.seed(0)

## Test standard methods for "lm"
x <- runif(20)
y <- x + rnorm(20)
m <- lm(y ~ x)
print(nObs(m))
print(stdEr(m))

## Test maxControl methods:
set.seed(9)
x <- rnorm(20)
ll <- function(mu) dnorm(x, mean=mu, log=TRUE)
for(method in c("NR", "BFGS", "BFGSR")) {
   cat("-- method", method, "--\n")
   m <- maxLik(ll, start=0, method=method, control=list(iterlim=1))
   expect_equal(maxValue(m), -27.5157350664669, tolerance=0.001)
   cat("MaxControl structure:\n")
   show(maxControl(m))
}
