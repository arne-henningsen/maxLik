library(maxLik)
set.seed(0)

## Test standard methods for "lm"
x <- runif(20)
y <- x + rnorm(20)
m <- lm(y ~ x)
print(nObs(m))
print(stdEr(m))

