# load the 'maxLik' package
library(maxLik)

## data to fit a normal distribution
# set seed for pseudo random numbers
set.seed( 123 )
# generate a variable from normally distributed random numbers
x <- rnorm( 100, 1, 2 )
xSaved <- x

## log likelihood function
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
}

## log likelihood function (individual observations)
llfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   return( llValues )
}

## function to calculate analytical gradients
gf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llGrad <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llGrad )
}

## function to calculate analytical gradients (individual observations)
gfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llGrads <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llGrads )
}

## function to calculate analytical Hessians
hf <- function( param ) {
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

# start values
startVal <- c( mu = 0, sigma = 1 )

## NR method
ml <- maxLik( llf, start = startVal )
print( ml )
summary( ml )
activePar( ml )
AIC( ml )
coef( ml )
condiNumber( ml )
hessian( ml )
logLik( ml )
maximType( ml )
nIter( ml )
nObs( ml )
nParam( ml )
returnCode( ml )
returnMessage( ml )
vcov( ml )
logLik( summary( ml ) )
mlInd <- maxLik( llfInd, start = startVal )
summary( mlInd )
all.equal( ml, mlInd )

# with analytical gradients
mlg <- maxLik( llf, gf, start = startVal )
summary( mlg )
all.equal( ml, mlg )
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
all.equal( mlInd, mlgInd )
all.equal( mlg, mlgInd )

# with analytical gradients and Hessians
mlgh <- maxLik( llf, gf, hf, start = startVal )
all.equal( mlg, mlgh )


## BHHH method
mlBHHH <- try( maxLik( llf, start = startVal, method = "BHHH" ) )
x <- xSaved[1]
try( maxLik( llfInd, start = startVal, method = "BHHH" ) )
x <- xSaved[1:2]
try( maxLik( llfInd, start = startVal, method = "BHHH" ) )
x <- xSaved
mlBHHH <- maxLik( llfInd, start = startVal, method = "BHHH" )
print( mlBHHH )
summary( mlBHHH )
activePar( mlBHHH )
AIC( mlBHHH )
coef( mlBHHH )
condiNumber( mlBHHH )
hessian( mlBHHH )
logLik( mlBHHH )
maximType( mlBHHH )
nIter( mlBHHH )
nObs( mlBHHH )
nParam( mlBHHH )
returnCode( mlBHHH )
returnMessage( mlBHHH )
vcov( mlBHHH )
logLik( summary( mlBHHH ) )
all.equal( ml, mlBHHH )

# with analytical gradients
mlgBHHH <- try( maxLik( llf, gf, start = startVal, method = "BHHH" ) )
mlgBHHH <- try( maxLik( llfInd, gf, start = startVal, method = "BHHH" ) )
x <- xSaved[1]
try( maxLik( llf, gfInd, start = startVal, method = "BHHH" ) )
try( maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ) )
x <- xSaved[1:2]
try( maxLik( llf, gfInd, start = startVal, method = "BHHH" ) )
try( maxLik( llfInd, gfInd, start = startVal, method = "BHHH" ) )
x <- xSaved
mlgBHHH <- maxLik( llfInd, gfInd, start = startVal, method = "BHHH" )
summary( mlgBHHH )
all.equal( mlBHHH, mlgBHHH )
all.equal( mlg, mlgBHHH )
mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlgBHHH2 )

# with unused Hessian
mlghBHHH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlghBHHH )


## BFGS method
mlBFGS <- maxLik( llf, start = startVal, method = "BFGS" )
print( mlBFGS )
summary( mlBFGS )
activePar( mlBFGS )
AIC( mlBFGS )
coef( mlBFGS )
condiNumber( mlBFGS )
hessian( mlBFGS )
logLik( mlBFGS )
maximType( mlBFGS )
nIter( mlBFGS )
nObs( mlBFGS )
nParam( mlBFGS )
returnCode( mlBFGS )
returnMessage( mlBFGS )
vcov( mlBFGS )
logLik( summary( mlBFGS ) )
all.equal( ml, mlBFGS )
mlIndBFGS <- maxLik( llfInd, start = startVal, method = "BFGS" )
summary( mlIndBFGS )
all.equal( mlBFGS, mlIndBFGS )
# with individual log likelihood values
mlIndBFGS <- maxLik( llfInd, start = startVal, method = "BFGS" )
summary( mlIndBFGS )
all.equal( mlBFGS, mlIndBFGS )

# with analytical gradients
mlgBFGS <- maxLik( llf, gf, start = startVal, method = "BFGS" )
summary( mlgBFGS )
all.equal( mlBFGS, mlgBFGS )
all.equal( mlg, mlgBFGS )
mlgIndBFGS <- maxLik( llfInd, gfInd, start = startVal, method = "BFGS" )
all.equal( mlgBFGS, mlgIndBFGS )

# with unused Hessian
mlghBFGS <- maxLik( llf, gf, hf, start = startVal, method = "BFGS" )
all.equal( mlgBFGS, mlghBFGS )


## NM method
mlNM <- maxLik( llf, start = startVal, method = "NM" )
print( mlNM )
summary( mlNM )
activePar( mlNM )
AIC( mlNM )
coef( mlNM )
condiNumber( mlNM )
hessian( mlNM )
logLik( mlNM )
maximType( mlNM )
nIter( mlNM )
nObs( mlNM )
nParam( mlNM )
returnCode( mlNM )
returnMessage( mlNM )
vcov( mlNM )
logLik( summary( mlNM ) )
all.equal( ml, mlNM )
# with individual log likelihood values
mlIndNM <- maxLik( llfInd, start = startVal, method = "NM" )
summary( mlIndNM )
all.equal( mlNM, mlIndNM )

# with unused analytical gradients
mlgNM <- maxLik( llf, gf, start = startVal, method = "NM" )
summary( mlgNM )
all.equal( mlNM, mlgNM )
# with individual log likelihood values and gradients
mlgIndNM <- maxLik( llfInd, gfInd, start = startVal, method = "NM" )
summary( mlgIndNM )
all.equal( mlgNM, mlgIndNM )

# with unused analytical gradients and Hessian
mlghNM <- maxLik( llf, gf, hf, start = startVal, method = "NM" )
all.equal( mlgNM, mlghNM )


## SANN method
set.seed( 345 )
mlSANN <- maxLik( llf, start = startVal, method = "SANN" )
print( mlSANN )
summary( mlSANN )
activePar( mlSANN )
AIC( mlSANN )
coef( mlSANN )
condiNumber( mlSANN )
hessian( mlSANN )
logLik( mlSANN )
maximType( mlSANN )
nIter( mlSANN )
nObs( mlSANN )
nParam( mlSANN )
returnCode( mlSANN )
returnMessage( mlSANN )
vcov( mlSANN )
logLik( summary( mlSANN ) )
all.equal( ml, mlSANN )
# with individual log likelihood values
mlIndSANN <- maxLik( llfInd, start = startVal, method = "SANN" )
summary( mlIndSANN )
all.equal( mlSANN, mlIndSANN )

# with unused analytical gradients
set.seed( 345 )
mlgSANN <- maxLik( llf, gf, start = startVal, method = "SANN" )
summary( mlgSANN )
all.equal( mlSANN, mlgSANN )
# with individual log likelihood values and gradients
mlgIndSANN <- maxLik( llfInd, gfInd, start = startVal, method = "SANN" )
summary( mlgIndSANN )
all.equal( mlgSANN, mlgIndSANN )

# with unused analytical gradients and Hessian
set.seed( 345 )
mlghSANN <- maxLik( llf, gf, hf, start = startVal, method = "SANN" )
all.equal( mlgSANN, mlghSANN )
