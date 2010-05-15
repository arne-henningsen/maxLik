## load the maxLik package
library( maxLik )

## fitting an exponential distribution by ML,
## e.g. estimation of an exponential duration model

# generate data
set.seed( 4 )
t <- rexp( 100, 2 )

# log-likelihood function, gradient, and Hessian
loglik <- function(theta) log(theta) - theta*t
loglikSum <- function(theta) sum( log(theta) - theta*t )
gradlik <- function(theta) 1/theta - t
gradlikSum <- function(theta) sum( 1/theta - t )
hesslik <- function(theta) -100/theta^2


## NR estimation
# Estimate with only function values
ml <- maxLik( loglik, start = 1 )
print( ml )
summary( ml )
nObs( ml )
print.default( ml )
# log-likelihood value summed over all observations
mlSum <- maxLik( loglikSum, start = 1 )
all.equal( mlSum[], ml[-11] )

# Estimate with analytic gradient
mlg <- maxLik( loglik, gradlik, start = 1 )
nObs( mlg )
all.equal( mlg, ml )
# gradient summed over all observations
mlgSum <- maxLik( loglikSum, gradlikSum, start = 1 )
all.equal( mlgSum[], mlg[-11] )

# Estimate with analytic gradient and Hessian
mlgh <- maxLik( loglik, gradlik, hesslik, start = 1 )
all.equal( mlgh, mlg )


## BHHH estimation
# Estimate with only function values
mlBhhh <- maxLik( loglik, start = 1, method = "BHHH" )
print( mlBhhh )
summary( mlBhhh )
nObs( mlBhhh )
all.equal( mlBhhh[ -c( 5, 6, 10 ) ], ml[ -c( 5, 6, 10 ) ] )

# Estimate with analytic gradient
mlgBhhh <- maxLik( loglik, gradlik, start = 1, method = "BHHH" )
nObs( mlgBhhh )
all.equal( mlgBhhh, mlBhhh )

# Estimate with analytic gradient and Hessian (unused during estimation)
mlghBhhh <- maxLik( loglik, gradlik, hesslik, start = 1, method = "BHHH" )
all.equal( mlghBhhh, mlgBhhh )

## BFGS estimation
# Estimate with only function values
mlBfgs <- maxLik( loglik, start = 1, method = "BFGS" )
print( mlBfgs )
summary( mlBfgs )
nObs( mlBfgs )
all.equal( mlBfgs[ -c( 5, 6, 9, 10, 11 ) ], ml[ -c( 5, 6, 9, 10 ) ] )
# log-likelihood value summed over all observations
mlSumBfgs <- maxLik( loglikSum, start = 1, method = "BFGS" )
all.equal( mlSumBfgs[], mlBfgs[-12] )

# Estimate with analytic gradient
mlgBfgs <- maxLik( loglik, gradlik, start = 1, method = "BFGS" )
nObs( mlgBfgs )
all.equal( mlgBfgs, mlBfgs )
# gradient summed over all observations
mlgSumBfgs <- maxLik( loglikSum, gradlikSum, start = 1, method = "BFGS" )
all.equal( mlgSumBfgs[], mlgBfgs[-12] )

# Estimate with analytic gradient and Hessian (unused during estimation)
mlghBfgs <- maxLik( loglik, gradlik, hesslik, start = 1, method = "BFGS" )
all.equal( mlghBfgs, mlgBfgs )

## NM estimation
# Estimate with only function values
mlNm <- maxLik( loglik, start = 1, method = "NM" )
print( mlNm )
summary( mlNm )
nObs( mlNm )
all.equal( mlNm[ -c( 5, 6, 9, 10, 11 ) ], ml[ -c( 5, 6, 9, 10 ) ] )

# Estimate with analytic gradient (unused during estimation)
mlgNm <- maxLik( loglik, gradlik, start = 1, method = "NM" )
nObs( mlgNm )
all.equal( mlgNm, mlNm )

# Estimate with analytic gradient and Hessian (both unused during estimation)
mlghNm <- maxLik( loglik, gradlik, hesslik, start = 1, method = "NM" )
all.equal( mlghNm, mlgNm )

## SANN estimation
# Estimate with only function values
mlSann <- maxLik( loglik, start = 1, method = "SANN" )
print( mlSann )
summary( mlSann )
nObs( mlSann )
all.equal( mlSann[ -c( 5, 6, 9, 10, 11 ) ], ml[ -c( 5, 6, 9, 10 ) ] )

# Estimate with analytic gradient (unused during estimation)
mlgSann <- maxLik( loglik, gradlik, start = 1, method = "SANN" )
nObs( mlgSann )
all.equal( mlgSann, mlSann )

# Estimate with analytic gradient and Hessian (both unused during estimation)
mlghSann <- maxLik( loglik, gradlik, hesslik, start = 1, method = "SANN" )
all.equal( mlghSann, mlgSann )


