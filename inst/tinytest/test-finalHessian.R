### Test the 'finalHessian' argument of optimization routines

## do not run unless 'NOT_CRAN' explicitly defined
## (Suggested by Sebastian Meyer and others)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    message("skipping final Hessian tests, not necessary on CRAN")
    q("no")
}
if(!requireNamespace("tinytest", quietly = TRUE)) {
   message("These tests require 'tinytest' package\n")
   q("no")
}
require(maxLik)
set.seed( 4 )

## log-likelihood function, gradient, and Hessian for 1-parameter case
## (exponential distribution)
##
## Individual likelihoods
ll1i <- function(theta) {
   if(!all(theta > 0))
       return(NA)
   log(theta) - theta*t
}
ll1 <- function(theta) sum( log(theta) - theta*t )
gr1i <- function(theta) 1/theta - t
## Aggregated likelihoods
gr1 <- function(theta) sum( 1/theta - t )
hs1 <- function(theta) -100/theta^2
t <- rexp( 100, 2 )

## the same functions for 2-variable case (normal distribution)
##
## Aggregated likelihoods
ll2 <- function( param ) {
   ## log likelihood function
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
}
gr2 <- function( param ) {
   ## function to calculate analytical gradients
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llGrad <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llGrad )
}

## Individual likelihoods
ll2i <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   return( llValues )
}
gr2i <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llGrads <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llGrads )
}

## analytical Hessians
hs2 <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llHess <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llHess )
}
x <- rnorm(100, 1, 2)


## NR
# Estimate with only function values (single parameter)
expect_silent(a <- maxLik(ll1, gr1, start = 1, method = "NR" ))
expect_equal(dim(hessian(a)), c(1,1))
expect_warning(b <- maxLik(ll1, gr1, start = 1, method = "NR", finalHessian="bhhh"))
                           # BHHH not possible
expect_null(hessian(b))
expect_silent(c <- maxLik(ll1i, gr1i, start = 1, method = "NR", finalHessian=FALSE))
expect_null(hessian(c))
expect_silent(d <- maxLik(ll1i, gr1i, start = 1, method = "NR", finalHessian="bhhh"))
expect_equal(dim(hessian(d)), c(1,1))
## (vector parameter)
expect_silent(a <- maxLik( ll2, gr2, start = c(0,1), method = "NR" ))
expect_equal(dim(hessian(a)), c(2,2))
expect_warning(b <- maxLik( ll2, gr2, start = c(0,1), method = "NR", finalHessian="bhhh"))
                           # BHHH not possible
expect_null(hessian(b))
expect_silent(c <- maxLik(ll2, gr2, start = c(0,1), method = "NR", finalHessian=FALSE))
expect_null(hessian(c))

## BFGSR
# Estimate with only function values (single parameter)
expect_silent(a <- maxLik(ll1, gr1, start = 1, method = "BFGSR"))
expect_equal(dim(hessian(a)), c(1,1))
expect_warning(b <- maxLik(ll1, gr1, start = 1, method = "BFGSR", finalHessian="bhhh"))
                           # should issue a warning as BHHH not possible
expect_null(hessian(b))
expect_silent(c <- maxLik( ll1i, gr1i, start = 1, method = "BFGSR", finalHessian=FALSE))
expect_null(hessian(c))
expect_silent(d <- maxLik(ll1i, gr1i, start = 1, method = "BFGSR", finalHessian="bhhh"))
expect_equal(dim(hessian(d)), c(1,1))
# Estimate with only function values (vector parameter)
expect_silent(a <- maxLik(ll2, gr2, start = c(0,1), method = "BFGSR" ))
expect_equal(dim(hessian(a)), c(2,2))
expect_warning(b <- maxLik(ll2, gr2, start = c(0,1), method = "BFGSR", finalHessian="bhhh"))
                           # should issue a warning as BHHH not possible
expect_null(hessian(b))
expect_silent(c <- maxLik(ll2, gr2, start = c(0,1), method = "BFGSR", finalHessian=FALSE))
expect_null(hessian(c))

### Nelder-Mead
## Individual observations only
expect_silent(b <- maxLik(ll2i, start = c(0,1), method = "NM", finalHessian="bhhh"))
expect_equal(dim(hessian(b)), c(2,2))
## Individual observations, summed gradient
expect_warning(b <- maxLik( ll2i, gr2, start = c(0,1), method = "NM", finalHessian="bhhh"))
                           # gradient is aggregated
                           # (yes, could do it based on individual likelihood and numeric gradient)
expect_null(hessian(b))
