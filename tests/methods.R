library(maxLik)
set.seed(0)

## Test standard methods for "lm"
x <- runif(100)
y <- x + rnorm(100)
m <- lm(y ~ x)
print(nObs(m))
print(stdEr(m))
