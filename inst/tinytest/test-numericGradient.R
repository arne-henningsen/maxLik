### test numeric methods, including NA-s in results
### 
if(!requireNamespace("tinytest", quietly = TRUE)) {
   cat("These tests require 'tinytest' package\n")
   q("no")
}
library(tinytest)
suppressPackageStartupMessages(library(maxLik))

f <- function(x) {
   if(x[1] <= 0)
       return(NA)
                           # support of x[1] is (0, Inf)
   return(c(log(x[1]),x[2]))
}

## Should work well
ng <- matrix(c(5.108256, 0, 0, 1), 2, 2)
expect_equal(numericGradient(f, c(0.2,1), eps=0.1),
             ng, tolerance = 1e-4)

## Warning: not a scalar function
expect_warning(numericHessian(f, t0=c(0.2,1), eps=0.1))

## Warning: NA-s as x[1] too close to 0
expect_warning(numericGradient(f, c(0.01,1), eps=0.1))

