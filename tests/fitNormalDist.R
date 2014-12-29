# load the 'maxLik' package
library(maxLik)
options(digits = 4)
                           # just to avoid so many differences when comparing these output files
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
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
}

## log likelihood function (individual observations)
llfInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
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

## log likelihood function with gradients as attributes
llfGrad <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   return( llValues )
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

## log likelihood function with gradients and Hessian as attributes
llfGradHess <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   attributes( llValue )$gradient <- c( sum( ( x - mu ) / sigma^2 ),
      - N / sigma + sum( ( x - mu )^2 / sigma^3 ) )
   attributes( llValue )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValue )
}

## log likelihood function with gradients as attributes (individual observations)
llfGradHessInd <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   if(!(sigma > 0))
       return(NA)
                           # to avoid warnings in the output
   N <- length( x )
   llValues <- -0.5 * log( 2 * pi ) - log( sigma ) -
      0.5 * ( x - mu )^2 / sigma^2
   attributes( llValues )$gradient <- cbind( ( x - mu ) / sigma^2,
      - 1 / sigma + ( x - mu )^2 / sigma^3 )
   attributes( llValues )$hessian <- matrix( c(
      N * ( - 1 / sigma^2 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      sum( - 2 * ( x - mu ) / sigma^3 ),
      N / sigma^2 + sum( - 3 * ( x - mu )^2 / sigma^4 ) ),
      nrow = 2, ncol = 2 )
   return( llValues )
}


# start values
startVal <- c( mu = 0, sigma = 1 )

## NR method
ml <- maxLik( llf, start = startVal )
print( ml )
print( summary( ml ), digits = 2 )
activePar( ml )
AIC( ml )
coef( ml )
condiNumber( ml )
round( hessian( ml ), 1 )
logLik( ml )
maximType( ml )
nIter( ml )
try( nObs( ml ) )
nParam( ml )
returnCode( ml )
returnMessage( ml )
round( vcov( ml ), 3 )
logLik( summary( ml ) )
mlInd <- maxLik( llfInd, start = startVal )
print( summary( mlInd ), digits = 2 )
all.equal( ml[-c(3,4,5,6)], mlInd[ -c(3,4,5,6,11) ], tolerance = 1e-3 )
                           # 3  gradient, should be close to 0, but may vary enormously in relative terms
mlInd[[11]][sample(nrow(mlInd[[11]]), 10),]
                           # just print a sample of 10
nObs( mlInd )

# with analytical gradients
mlg <- maxLik( llf, gf, start = startVal )
print( summary( mlg ), digits = 2 )
all.equal( ml[-c(5,6)], mlg[-c(5,6)], tolerance = 1e-3 )
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
all.equal( mlInd, mlgInd, tolerance = 1e-3 )
all.equal( mlg[ ], mlgInd[ -11 ], tolerance = 1e-3 )
round( mlgInd[[ 11 ]], 3 )

# with analytical gradients as attribute
mlG <- maxLik( llfGrad, start = startVal )
all.equal( mlG, mlg, tolerance = 1e-3 )
all.equal( mlG$gradient, gf( coef( mlG ) ), check.attributes = FALSE,
   tolerance = 1e-3 )
mlGInd <- maxLik( llfGradInd, start = startVal )
all.equal( mlGInd, mlgInd, tolerance = 1e-3 )
all.equal( mlGInd$gradient, colSums( gfInd( coef( mlGInd ) ) ),
   check.attributes = FALSE, tolerance = 1e-3 )
all.equal( mlGInd$gradientObs, gfInd( coef( mlGInd ) ),
   check.attributes = FALSE, tolerance = 1e-3 )

# with analytical gradients as argument and attribute
mlgG <- maxLik( llfGrad, gf, start = startVal )
all.equal( mlgG, mlg, tolerance = 1e-3 )
all.equal( mlgG, mlG, tolerance = 1e-3 )

# with analytical gradients and Hessians
mlgh <- maxLik( llf, gf, hf, start = startVal )
all.equal( mlg, mlgh, tolerance = 1e-3 )

# with analytical gradients and Hessian as attribute
mlGH <- maxLik( llfGradHess, start = startVal )
all.equal( mlGH, mlgh, tolerance = 1e-3 )

# with analytical gradients and Hessian as argument and attribute
mlgGhH <- maxLik( llfGradHess, gf, hf, start = startVal )
all.equal( mlgGhH, mlgh, tolerance = 1e-3 )
all.equal( mlgGhH, mlGH, tolerance = 1e-3 )


## BHHH method
mlBHHH <- try( maxLik( llf, start = startVal, method = "BHHH" ) )
x <- xSaved[1]
try( maxLik( llfInd, start = startVal, method = "BHHH" ) )
x <- xSaved[1:2]
try( maxLik( llfInd, start = startVal, method = "BHHH" ) )
x <- xSaved
mlBHHH <- maxLik( llfInd, start = startVal, method = "BHHH" )
print( mlBHHH )
print( summary( mlBHHH ), digits = 2 )
activePar( mlBHHH )
AIC( mlBHHH )
coef( mlBHHH )
condiNumber( mlBHHH )
round( hessian( mlBHHH ), 1 )
logLik( mlBHHH )
maximType( mlBHHH )
nIter( mlBHHH )
nParam( mlBHHH )
returnCode( mlBHHH )
returnMessage( mlBHHH )
round( vcov( mlBHHH ), 3 )
logLik( summary( mlBHHH ) )
all.equal( ml[-c(4,5,6,9,10) ], mlBHHH[ -c(4,5,6,9,10,11) ], tolerance = 1e-3 )
round( mlBHHH[[ 11 ]], 3 )
nObs( mlBHHH )
# final Hessian = usual Hessian
mlBhhhH <- maxLik( llfInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlBhhhH[-4], mlBHHH[-4], tolerance = 1e-3 )
round( hessian( mlBhhhH ), 1 )
print( summary( mlBhhhH ) , digits = 2 )

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
print( summary( mlgBHHH ), digits = 2 )
all.equal( mlBHHH, mlgBHHH, tolerance = 1e-3 )
all.equal( mlg[-c(4,5,6,9,10)], mlgBHHH[-c(4,5,6,9,10,11)], tolerance = 1e-3 )
round( mlgBHHH[[ 11 ]], 3 )
mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlgBHHH2, tolerance = 1e-3 )
# final Hessian = usual Hessian
mlgBhhhH <- maxLik( llf, gfInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlgBhhhH, mlBhhhH, tolerance = 1e-3 )
all.equal( mlgBhhhH[-4], mlgBHHH[-4], tolerance = 1e-3 )
round( hessian( mlgBhhhH ), 1 )

# with analytical gradients as attribute
try( maxLik( llfGrad, start = startVal, method = "BHHH" ) )
x <- xSaved[1]
try( maxLik( llfGrad, start = startVal, method = "BHHH" ) )
try( maxLik( llfGradInd, start = startVal, method = "BHHH" ) )
x <- xSaved[1:2]
try( maxLik( llfGrad, start = startVal, method = "BHHH" ) )
try( maxLik( llfGradInd, start = startVal, method = "BHHH" ) )
x <- xSaved
mlGBHHH <- maxLik( llfGradInd, start = startVal, method = "BHHH" )
all.equal( mlGBHHH, mlgBHHH, tolerance = 1e-3 )
# final Hessian = usual Hessian
mlGBhhhH <- maxLik( llfGradInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlGBhhhH, mlgBhhhH, tolerance = 1e-3 )

# with analytical gradients as argument and attribute
mlgGBHHH <- maxLik( llfGradInd, gfInd, start = startVal, method = "BHHH" )
all.equal( mlgGBHHH, mlgBHHH, tolerance = 1e-3 )
all.equal( mlgGBHHH, mlGBHHH, tolerance = 1e-3 )

# with unused Hessian
mlghBHHH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlghBHHH, tolerance = 1e-3 )
# final Hessian = usual Hessian
mlghBhhhH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlghBhhhH[-4], mlghBHHH[-4], tolerance = 1e-3 )
all.equal( mlghBhhhH, mlgBhhhH, tolerance = 1e-3 )

# with unused Hessian as attribute
mlGHBHHH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH" )
all.equal( mlGHBHHH, mlghBHHH, tolerance = 1e-3 )
# final Hessian = usual Hessian
mlGHBhhhH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlGHBhhhH, mlghBhhhH, tolerance = 1e-3 )

# with analytical gradients and Hessian as argument and attribute
mlgGhHBHHH <- maxLik( llfGradHessInd, gfInd, hf, start = startVal, method = "BHHH" )
all.equal( mlgGhHBHHH, mlghBHHH, tolerance = 1e-3 )
all.equal( mlgGhHBHHH, mlGHBHHH, tolerance = 1e-3 )


### BFGS-YC method
mlBFGSYC <- maxLik( llf, start = startVal, method = "bfgsr" )
print( mlBFGSYC )
print( summary( mlBFGSYC ), digits = 2 )
activePar( mlBFGSYC )
AIC( mlBFGSYC )
coef( mlBFGSYC )
condiNumber( mlBFGSYC )
round( hessian( mlBFGSYC ), 1 )
logLik( mlBFGSYC )
maximType( mlBFGSYC )
nIter( mlBFGSYC )
try( nObs( mlBFGSYC ) )
nParam( mlBFGSYC )
returnCode( mlBFGSYC )
returnMessage( mlBFGSYC )
round( vcov( mlBFGSYC ), 3 )
logLik( summary( mlBFGSYC ) )
all.equal( ml[-c(3,4,5,6,9,10)], mlBFGSYC[-c(3,4,5,6,9,10)], tolerance = 1e-3 )
all.equal( ml[-c(5,6,9,10)], mlBFGSYC[-c(5,6,9,10)], tolerance = 1e-2 )
mlIndBFGSYC <- maxLik( llfInd, start = startVal, method = "BFGSR" )
print( summary( mlIndBFGSYC ), digits = 2 )
all.equal( mlBFGSYC[-c(3,4,9)], mlIndBFGSYC[ -c(3,4,9,11) ], tolerance = 1e-3 )
round( mlIndBFGSYC[[ 11 ]], 3 )
nObs( mlIndBFGSYC )

# with analytical gradients
mlgBFGSYC <- maxLik( llf, gf, start = startVal, method = "BFGSR" , print.level=1)
print( summary(mlgBFGSYC), digits = 2 )
all.equal( mlBFGSYC[-4], mlgBFGSYC[-4], tolerance = 1e-3 )
mlgIndBFGSYC <- maxLik( llfInd, gfInd, start = startVal,
   method = "BFGSR" )
all.equal( mlIndBFGSYC, mlgIndBFGSYC, tolerance = 1e-3 )
all.equal( mlgBFGSYC[ -c(3,9) ], mlgIndBFGSYC[ -c(3,9,11) ], tolerance = 1e-3 )
round( mlgIndBFGSYC[[ 11 ]], 3 )

# with analytical gradients as attribute
mlGBFGSYC <- maxLik( llfGrad, start = startVal, method = "BFGSR" , print.level=1)
all.equal( mlGBFGSYC, mlgBFGSYC, tolerance = 1e-3 )
mlGIndBFGSYC <- maxLik( llfGradInd, start = startVal, method = "BFGSR" )
all.equal( mlGIndBFGSYC, mlgIndBFGSYC, tolerance = 1e-3 )

# with analytical gradients as argument and attribute
mlgGBFGSYC <- maxLik( llfGrad, gf, start = startVal, method = "BFGSR" )
all.equal( mlgGBFGSYC, mlgBFGSYC, tolerance = 1e-3 )
all.equal( mlgGBFGSYC, mlGBFGSYC, tolerance = 1e-3 )

# with analytical gradients and Hessians
mlghBFGSYC <- maxLik( llf, gf, hf, start = startVal, method = "BFGSR" )
all.equal( mlgBFGSYC, mlghBFGSYC, tolerance = 1e-3 )

# with analytical gradients and Hessian as attribute
mlGHBFGSYC <- maxLik( llfGradHess, start = startVal, method = "BFGSR" )
all.equal( mlGHBFGSYC, mlghBFGSYC, tolerance = 1e-3 )

# with analytical gradients and Hessian as argument and attribute
mlgGhHBFGSYC <- maxLik( llfGradHess, gf, hf, start = startVal, method = "BFGSR" )
all.equal( mlgGhHBFGSYC, mlghBFGSYC, tolerance = 1e-3 )
all.equal( mlgGhHBFGSYC, mlGHBFGSYC, tolerance = 1e-3 )


## BFGS method
mlBFGS <- maxLik( llf, start = startVal, method = "BFGS" )
print( mlBFGS )
print( summary( mlBFGS ), digits = 2 )
activePar( mlBFGS )
AIC( mlBFGS )
coef( mlBFGS )
condiNumber( mlBFGS )
round( hessian( mlBFGS ), 1 )
logLik( mlBFGS )
maximType( mlBFGS )
nIter( mlBFGS )
nParam( mlBFGS )
returnCode( mlBFGS )
returnMessage( mlBFGS )
round( vcov( mlBFGS ), 3 )
logLik( summary( mlBFGS ) )
all.equal( ml[-c(4,5,6,9,10)], mlBFGS[-c(4,5,6,9,10,11)], tolerance = 1e-3 )
# with individual log likelihood values
mlIndBFGS <- maxLik( llfInd, start = startVal, method = "BFGS" )
print( summary( mlIndBFGS ), digits = 2 )
all.equal( mlBFGS[-4], mlIndBFGS[-c(4,12)], tolerance = 1e-3 )
mlIndBFGS[12]
nObs( mlIndBFGS )

# with analytical gradients
mlgBFGS <- maxLik( llf, gf, start = startVal, method = "BFGS" )
print( summary( mlgBFGS ), digits = 2 )
all.equal( mlBFGS[-4], mlgBFGS[-4], tolerance = 1e-3 )
all.equal( mlg[-c(5,6,9,10)], mlgBFGS[-c(5,6,9,10,11)], tolerance = 1e-3 )
mlgIndBFGS <- maxLik( llfInd, gfInd, start = startVal, method = "BFGS" )
all.equal( mlgBFGS[], mlgIndBFGS[-12], tolerance = 1e-3 )
mlgIndBFGS[12]

# with analytical gradients as attribute
mlGBFGS <- maxLik( llfGrad, start = startVal, method = "BFGS" )
all.equal( mlGBFGS, mlgBFGS, tolerance = 1e-3 )
mlGIndBFGS <- maxLik( llfGradInd, start = startVal, method = "BFGS" )
all.equal( mlGIndBFGS, mlgIndBFGS, tolerance = 1e-3 )

# with analytical gradients as argument and attribute
mlgGBFGS <- maxLik( llfGrad, gf, start = startVal, method = "BFGS" )
all.equal( mlgGBFGS, mlgBFGS, tolerance = 1e-3 )
all.equal( mlgGBFGS, mlGBFGS, tolerance = 1e-3 )

# with unused Hessian
mlghBFGS <- maxLik( llf, gf, hf, start = startVal, method = "BFGS" )
all.equal( mlgBFGS, mlghBFGS, tolerance = 1e-3 )

# with analytical gradients and Hessian as attribute
mlGHBFGS <- maxLik( llfGradHess, start = startVal, method = "BFGS" )
all.equal( mlGHBFGS, mlghBFGS, tolerance = 1e-3 )

# with analytical gradients and Hessian as argument and attribute
mlgGhHBFGS <- maxLik( llfGradHess, gf, hf, start = startVal, method = "BFGS" )
all.equal( mlgGhHBFGS, mlghBFGS, tolerance = 1e-3 )
all.equal( mlgGhHBFGS, mlGHBFGS, tolerance = 1e-3 )


## NM method
mlNM <- maxLik( llf, start = startVal, method = "NM" )
print( mlNM )
print( summary( mlNM ), digits = 2 )
activePar( mlNM )
AIC( mlNM )
coef( mlNM )
condiNumber( mlNM )
round( hessian( mlNM ), 1 )
logLik( mlNM )
maximType( mlNM )
nIter( mlNM )
nParam( mlNM )
returnCode( mlNM )
returnMessage( mlNM )
round( vcov( mlNM ), 3 )
logLik( summary( mlNM ) )
all.equal( ml[-c(3,4,5,6,9,10)], mlNM[-c(3,4,5,6,9,10,11)], tolerance = 1e-3 )
# with individual log likelihood values
mlIndNM <- maxLik( llfInd, start = startVal, method = "NM" )
print( summary( mlIndNM ), digits = 2 )
all.equal( mlNM[-4], mlIndNM[-c(4,12)], tolerance = 1e-3 )
mlIndNM[12]
nObs( mlIndNM )

# with unused analytical gradients
mlgNM <- maxLik( llf, gf, start = startVal, method = "NM" )
print( summary( mlgNM ), digits = 2 )
all.equal( mlNM[-4], mlgNM[-4], tolerance = 1e-3 )
# with individual log likelihood values and gradients
mlgIndNM <- maxLik( llfInd, gfInd, start = startVal, method = "NM" )
print( summary( mlgIndNM ), digits = 2 )
all.equal( mlgNM[], mlgIndNM[-12], tolerance = 1e-3 )
mlgIndNM[12]

# with (unused) analytical gradients as attribute
mlGNM <- maxLik( llfGrad, start = startVal, method = "NM" )
all.equal( mlGNM, mlgNM, tolerance = 1e-3 )
mlGIndNM <- maxLik( llfGradInd, start = startVal, method = "NM" )
all.equal( mlGIndNM, mlgIndNM, tolerance = 1e-3 )

# with analytical gradients as argument and attribute
mlgGNM <- maxLik( llfGrad, gf, start = startVal, method = "NM" )
all.equal( mlgGNM, mlgNM, tolerance = 1e-3 )
all.equal( mlgGNM, mlGNM, tolerance = 1e-3 )

# with unused analytical gradients and Hessian
mlghNM <- maxLik( llf, gf, hf, start = startVal, method = "NM" )
all.equal( mlgNM, mlghNM, tolerance = 1e-3 )


## SANN method
mlSANN <- maxLik( llf, start = startVal, method = "SANN" )
print( mlSANN )
print( summary( mlSANN ), digits = 2 )
activePar( mlSANN )
AIC( mlSANN )
coef( mlSANN )
condiNumber( mlSANN )
round( hessian( mlSANN ), 1 )
logLik( mlSANN )
maximType( mlSANN )
nIter( mlSANN )
nParam( mlSANN )
returnCode( mlSANN )
returnMessage( mlSANN )
round( vcov( mlSANN ), 3 )
logLik( summary( mlSANN ) )
all.equal( ml[-c(3,4,5,6,9,10)], mlSANN[-c(3,4,5,6,9,10,11)], tolerance = 1e-3 )
# with individual log likelihood values
mlIndSANN <- maxLik( llfInd, start = startVal, method = "SANN" )
print( summary( mlIndSANN ), digits = 2 )
all.equal( mlSANN[-4], mlIndSANN[-c(4,12)], tolerance = 1e-3 )
mlIndSANN[12]
nObs( mlIndSANN )

# with unused analytical gradients
mlgSANN <- maxLik( llf, gf, start = startVal, method = "SANN" )
print( summary( mlgSANN ), digits = 2 )
all.equal( mlSANN[-4], mlgSANN[-4], tolerance = 1e-3 )
# with individual log likelihood values and gradients
mlgIndSANN <- maxLik( llfInd, gfInd, start = startVal, method = "SANN" )
print( summary( mlgIndSANN ), digits = 2 )
all.equal( mlgSANN[], mlgIndSANN[-12], tolerance = 1e-3 )
mlgIndSANN[12]

# with unused analytical gradients and Hessian
mlghSANN <- maxLik( llf, gf, hf, start = startVal, method = "SANN" )
all.equal( mlgSANN, mlghSANN, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSANNCand <- maxLik( llf, start = startVal, method = "SANN",
   cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
print( summary( mlSANNCand ), digits = 2 )
all.equal( mlSANNCand[-c(3,4)], mlSANN[-c(3,4)], tolerance = 1e-2 )

############### with fixed parameters ###############
# start values
startValFix <- c( mu = 1, sigma = 1 )

# fix mu (the mean ) at its start value
isFixed <- c( TRUE, FALSE )

## NR method with fixed parameters
mlFix <- maxLik( llf, start = startValFix, activePar = !isFixed )
mlFix1 <- maxLik( llf, start = startValFix, activePar = 2 )
all.equal( mlFix, mlFix1, tolerance = 1e-3 )
mlFix2 <- maxLik( llf, start = startValFix, fixed = isFixed )
all.equal( mlFix, mlFix2, tolerance = 1e-3 )
mlFix3 <- maxLik( llf, start = startValFix, fixed = "mu" )
all.equal( mlFix, mlFix3, tolerance = 1e-3 )
mlFix4 <- maxLik( llf, start = startValFix, fixed = 1 )
all.equal( mlFix, mlFix4, tolerance = 1e-3 )
print( mlFix )
print( summary( mlFix ), digits = 2 )
activePar( mlFix )
AIC( mlFix )
coef( mlFix )
condiNumber( mlFix )
round( hessian( mlFix ), 1 )
logLik( mlFix )
maximType( mlFix )
nIter( mlFix )
nParam( mlFix )
returnCode( mlFix )
returnMessage( mlFix )
round( vcov( mlFix ), 3 )
logLik( summary( mlFix ) )
mlIndFix <- maxLik( llfInd, start = startValFix, activePar = !isFixed )
mlIndFix1 <- maxLik( llfInd, start = startValFix, activePar = 2 )
all.equal( mlIndFix, mlIndFix1, tolerance = 1e-3 )
mlIndFix2 <- maxLik( llfInd, start = startValFix, fixed = isFixed )
all.equal( mlIndFix, mlIndFix2, tolerance = 1e-3 )
mlIndFix3 <- maxLik( llfInd, start = startValFix, fixed = "mu" )
all.equal( mlIndFix, mlIndFix3, tolerance = 1e-3 )
mlIndFix4 <- maxLik( llfInd, start = startValFix, fixed = 1 )
all.equal( mlIndFix, mlIndFix4, tolerance = 1e-3 )
print( summary( mlIndFix ), digits = 2 )
all.equal( mlFix[ ], mlIndFix[ -11 ], tolerance = 1e-3 )
round( mlFix[[3]], 5 )
round( mlIndFix[[3]], 5 )
round( mlIndFix[[ 11 ]], 3 )
nObs( mlIndFix )

# with analytical gradients
mlgFix <- maxLik( llf, gf, start = startValFix, activePar = !isFixed )
mlgFix1 <- maxLik( llf, gf, start = startValFix, activePar = 2 )
all.equal( mlgFix, mlgFix1, tolerance = 1e-3 )
mlgFix2 <- maxLik( llf, gf, start = startValFix, fixed = isFixed )
all.equal( mlgFix, mlgFix2, tolerance = 1e-3 )
print( summary( mlgFix ), digits = 2 )
mlgIndFix <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed )
all.equal( mlIndFix, mlgIndFix, tolerance = 1e-3 )
all.equal( mlgFix[ ], mlgIndFix[ -11 ], tolerance = 1e-3 )
round( mlgIndFix[[ 11 ]], 3 )

# with analytical gradients and Hessians
mlghFix <- maxLik( llf, gf, hf, start = startValFix, activePar = !isFixed )
all.equal( mlgFix, mlghFix, tolerance = 1e-3 )
mlgFix[[4]]
mlghFix[[4]]

## BHHH method with fixed parameters
mlFixBHHH <- maxLik( llfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
mlFixBHHH1 <- maxLik( llfInd, start = startValFix, activePar = 2,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH1, tolerance = 1e-3 )
mlFixBHHH2 <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH2, tolerance = 1e-3 )
mlFixBHHH3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH3, tolerance = 1e-3 )
mlFixBHHH4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH4, tolerance = 1e-3 )
print( mlFixBHHH )
print( summary( mlFixBHHH ), digits = 2 )
activePar( mlFixBHHH )
AIC( mlFixBHHH )
coef( mlFixBHHH )
condiNumber( mlFixBHHH )
round( hessian( mlFixBHHH ), 1 )
logLik( mlFixBHHH )
maximType( mlFixBHHH )
nIter( mlFixBHHH )
nParam( mlFixBHHH )
returnCode( mlFixBHHH )
returnMessage( mlFixBHHH )
round( vcov( mlFixBHHH ), 3 )
logLik( summary( mlFixBHHH ) )
all.equal( mlFix[ -c( 4, 5, 6, 9, 10 ) ], mlFixBHHH[ -c( 4, 5, 6, 9, 10, 11 ) ],
   tolerance = 1e-3 )
round( mlFixBHHH[[ 11 ]], 3 )
nObs( mlFixBHHH )

# with analytical gradients
mlgFixBHHH <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
mlgFixBHHH1 <- maxLik( llfInd, gfInd, start = startValFix, activePar = 2,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH1, tolerance = 1e-3 )
mlgFixBHHH2 <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH2, tolerance = 1e-3 )
mlgFixBHHH3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH3, tolerance = 1e-3 )
mlgFixBHHH4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH4, tolerance = 1e-3 )
print( summary( mlgFixBHHH ), digits = 2 )
all.equal( mlFixBHHH, mlgFixBHHH, tolerance = 1e-3 )
mlgFixBHHH2 <- maxLik( llf, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH")
all.equal( mlgFixBHHH, mlgFixBHHH2, tolerance = 1e-3 )

# with unused Hessians
mlghFixBHHH <- maxLik( llfInd, gfInd, hf, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlghFixBHHH, tolerance = 1e-3 )

## BFGS method with fixed parameters
mlFixBfgs <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
mlFixBfgs3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlFixBfgs, mlFixBfgs3, tolerance = 1e-3 )
mlFixBfgs4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlFixBfgs, mlFixBfgs4, tolerance = 1e-3 )
print( mlFixBfgs )
print( summary( mlFixBfgs ), digits = 2 )
activePar( mlFixBfgs )
AIC( mlFixBfgs )
coef( mlFixBfgs )
condiNumber( mlFixBfgs )
round( hessian( mlFixBfgs ), 1 )
logLik( mlFixBfgs )
maximType( mlFixBfgs )
nIter( mlFixBfgs )
nParam( mlFixBfgs )
returnCode( mlFixBfgs )
returnMessage( mlFixBfgs )
round( vcov( mlFixBfgs ), 3 )
logLik( summary( mlFixBfgs ) )
all.equal( mlghFix[ -c( 5, 6, 9, 10 ) ], mlFixBfgs[ -c( 5, 6, 9, 10, 11 ) ],
   tolerance = 1e-3 )
mlIndFixBfgs <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "BFGS" )
all.equal( mlFixBfgs[-c(4,9)], mlIndFixBfgs[ -c(4,9,12) ], tolerance = 1e-3 )
print(formatC(mlIndFixBfgs$gradientObs, format="f", digits=4, width=7), quote=FALSE)
                           # print fradient, only 4 digits to avoid clutter in R CMD tests
mlIndFixBfgs3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlIndFixBfgs, mlIndFixBfgs3, tolerance = 1e-3 )
mlIndFixBfgs4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlIndFixBfgs, mlIndFixBfgs4, tolerance = 1e-3 )
nObs( mlIndFixBfgs )

# with analytical gradients
mlgFixBfgs <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
mlgFixBfgs3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlgFixBfgs, mlgFixBfgs3, tolerance = 1e-3 )
mlgFixBfgs4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlgFixBfgs, mlgFixBfgs4, tolerance = 1e-3 )
print( summary( mlgFixBfgs ), digits = 2 )
all.equal( mlFixBfgs[ -9 ], mlgFixBfgs[ -9 ], tolerance = 1e-3 )
mlgIndFixBfgs <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "BFGS")
all.equal( mlgFixBfgs[ ], mlgIndFixBfgs[ -12 ], tolerance = 1e-3 )
round( mlgIndFixBfgs[[ 12 ]], 3 )
mlgIndFixBfgs3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlgIndFixBfgs, mlgIndFixBfgs3, tolerance = 1e-3 )
mlgIndFixBfgs4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlgIndFixBfgs, mlgIndFixBfgs4, tolerance = 1e-3 )

# with unused Hessians
mlghFixBfgs <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
all.equal( mlgFixBfgs, mlghFixBfgs, tolerance = 1e-3 )
mlghFixBfgs3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlghFixBfgs, mlghFixBfgs3, tolerance = 1e-3 )
mlghFixBfgs4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlghFixBfgs, mlghFixBfgs4, tolerance = 1e-3 )

## NM method with fixed parameters
mlFixNm <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "NM" )
mlFixNm3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlFixNm, mlFixNm3, tolerance = 1e-3 )
mlFixNm4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlFixNm, mlFixNm4, tolerance = 1e-3 )
print( mlFixNm )
print( summary( mlFixNm ), digits = 2 )
activePar( mlFixNm )
AIC( mlFixNm )
coef( mlFixNm )
condiNumber( mlFixNm )
round( hessian( mlFixNm ), 1 )
logLik( mlFixNm )
maximType( mlFixNm )
nIter( mlFixNm )
nParam( mlFixNm )
returnCode( mlFixNm )
returnMessage( mlFixNm )
round( vcov( mlFixNm ), 3 )
logLik( summary( mlFixNm ) )
all.equal( mlFixBfgs[ -c(4,9,10) ], mlFixNm[ -c(4,9,10) ], tolerance = 1e-3 )
mlIndFixNm <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "NM" )
all.equal( mlFixNm[-4], mlIndFixNm[-c(4,12)], tolerance = 1e-3 )
round( mlIndFixNm[[ 12 ]], 3 )
mlIndFixNm3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlIndFixNm, mlIndFixNm3, tolerance = 1e-3 )
mlIndFixNm4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlIndFixNm, mlIndFixNm4, tolerance = 1e-3 )
nObs( mlIndFixNm )

# with analytical gradients
mlgFixNm <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "NM" )
mlgFixNm3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlgFixNm, mlgFixNm3, tolerance = 1e-3 )
mlgFixNm4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlgFixNm, mlgFixNm4, tolerance = 1e-3 )
print( summary( mlgFixNm ), digits = 2 )
all.equal( mlFixNm, mlgFixNm, tolerance = 1e-3 )
mlgIndFixNm <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "NM")
all.equal( mlgFixNm[ ], mlgIndFixNm[ -12 ], tolerance = 1e-3 )
round( mlgIndFixNm[[ 12 ]], 3 )

# with unused Hessians
mlghFixNm <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "NM" )
all.equal( mlgFixNm, mlghFixNm, tolerance = 1e-3 )
mlghFixNm3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlghFixNm, mlghFixNm3, tolerance = 1e-3 )
mlghFixNm4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlghFixNm, mlghFixNm4, tolerance = 1e-3 )

## SANN method with fixed parameters
mlFixSann <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "SANN" )
mlFixSann3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "SANN" )
all.equal( mlFixSann, mlFixSann3, tolerance = 1e-3 )
mlFixSann4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "SANN" )
all.equal( mlFixSann, mlFixSann4, tolerance = 1e-3 )
print( mlFixSann )
print( summary( mlFixSann ), digits = 2 )
activePar( mlFixSann )
AIC( mlFixSann )
coef( mlFixSann )
condiNumber( mlFixSann )
round( hessian( mlFixSann ), 1 )
logLik( mlFixSann )
maximType( mlFixSann )
nIter( mlFixSann )
nParam( mlFixSann )
returnCode( mlFixSann )
returnMessage( mlFixSann )
round( vcov( mlFixSann ), 3 )
logLik( summary( mlFixSann ) )
all.equal( mlFixBfgs[ -c(4,9,10) ], mlFixSann[ -c(4,9,10) ], 
   tolerance = 1e-3 )
mlIndFixSann <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "SANN" )
all.equal( mlFixSann[ ], mlIndFixSann[ -12 ], tolerance = 1e-2 )
round( mlIndFixSann[[ 12 ]], 3 )
nObs( mlIndFixSann )

# with analytical gradients
mlgFixSann <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "SANN" )
print( summary( mlgFixSann ), digits = 2 )
all.equal( mlFixSann[-4], mlgFixSann[-4], tolerance = 1e-3 )
mlgIndFixSann <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "SANN")
all.equal( mlgFixSann[ ], mlgIndFixSann[ -12 ], tolerance = 1e-3 )
round( mlgIndFixSann[[ 12 ]], 3 )

# with unused Hessians
mlghFixSann <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "SANN" )
all.equal( mlgFixSann, mlghFixSann, tolerance = 1e-3 )


############### with parameter constraints ###############
A <- matrix( -1, nrow = 1, ncol = 2 )


############### inequality constraints ###############
inEq <- list( ineqA = A, ineqB = 2.5 )

## NR method with inequality constraints
try( maxLik( llf, start = startVal, constraints = inEq, method = "NR" ) )

## BHHH method with inequality constraints
try( maxLik( llf, start = startVal, constraints = inEq, method = "BHHH" ) )

## BFGS method with inequality constraints
mlBfgsInEq <- maxLik( llf, start = startVal, constraints = inEq,
   method = "BFGS" )
print( mlBfgsInEq )
print( summary( mlBfgsInEq ), digits = 2 )
activePar( mlBfgsInEq )
AIC( mlBfgsInEq )
coef( mlBfgsInEq )
condiNumber( mlBfgsInEq )
round( hessian( mlBfgsInEq ), 1 )
logLik( mlBfgsInEq )
maximType( mlBfgsInEq )
nIter( mlBfgsInEq )
nParam( mlBfgsInEq )
returnCode( mlBfgsInEq )
returnMessage( mlBfgsInEq )
round( vcov( mlBfgsInEq ), 3 )
logLik( summary( mlBfgsInEq ) )
mlBfgsInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
print( summary( mlBfgsInEqInd ), digits = 2 )
all.equal( mlBfgsInEq[ ], mlBfgsInEqInd[ -12 ], tolerance = 1e-3 )
round( mlBfgsInEqInd[[ 12 ]], 3 )
nObs( mlBfgsInEqInd )

# with analytical gradients
mlgBfgsInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlBfgsInEq, mlgBfgsInEq, tolerance = 1e-3 )
mlgBfgsInEqInd <- maxLik( llfInd, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd[ -12 ], mlgBfgsInEq[ ], tolerance = 1e-3 )
round( mlgBfgsInEqInd[[ 12 ]], 3 )
mlgBfgsInEqInd2 <- maxLik( llf, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd, mlgBfgsInEqInd2, tolerance = 1e-3 )

# with unused Hessian
mlghBfgsInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEq, mlghBfgsInEq, tolerance = 1e-3 )

## NM method with inequality constraints
mlNmInEq <- maxLik( llf, start = startVal, constraints = inEq, method = "NM" )
print( mlNmInEq )
print( summary( mlNmInEq ), digits = 2 )
activePar( mlNmInEq )
AIC( mlNmInEq )
coef( mlNmInEq )
condiNumber( mlNmInEq )
round( hessian( mlNmInEq ), 1 )
logLik( mlNmInEq )
maximType( mlNmInEq )
nIter( mlNmInEq )
nParam( mlNmInEq )
returnCode( mlNmInEq )
returnMessage( mlNmInEq )
round( vcov( mlNmInEq ), 3 )
logLik( summary( mlNmInEq ) )
all.equal( mlBfgsInEq[-c(9,10,11)], mlNmInEq[-c(9,10,11)], tolerance = 1e-3 )
mlNmInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "NM" )
print( summary( mlNmInEqInd ), digits = 2 )
all.equal( mlNmInEq[-4], mlNmInEqInd[-c(4,12)], tolerance = 1e-3 )
round( mlNmInEqInd[[ 12 ]], 3 )
nObs( mlNmInEqInd )

# with unused analytical gradients
mlgNmInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlNmInEq, mlgNmInEq, tolerance = 1e-3 )

# with unused analytical gradients and Hessians
mlghNmInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlgNmInEq, mlghNmInEq, tolerance = 1e-3 )

## SANN method with inequality constraints
mlSannInEq <- maxLik( llf, start = startVal, constraints = inEq,
   method = "SANN" )
print( mlSannInEq )
print( summary( mlSannInEq ), digits = 2 )
activePar( mlSannInEq )
AIC( mlSannInEq )
coef( mlSannInEq )
condiNumber( mlSannInEq )
round( hessian( mlSannInEq ), 1 )
logLik( mlSannInEq )
maximType( mlSannInEq )
nIter( mlSannInEq )
nParam( mlSannInEq )
returnCode( mlSannInEq )
returnMessage( mlSannInEq )
round( vcov( mlSannInEq ), 3 )
logLik( summary( mlSannInEq ) )
all.equal( mlBfgsInEq[-c(2,3,4,9,10,11)], mlSannInEq[-c(2,3,4,9,10,11)], 
   tolerance = 1e-3 )
all.equal( mlBfgsInEq[-c(3,4,9,10,11)], mlSannInEq[-c(3,4,9,10,11)], 
   tolerance = 1e-2 )
# with unused analytical gradients
mlgSannInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "SANN" )
all.equal( mlSannInEq, mlgSannInEq, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSannInEqCand <- maxLik( llf, start = startVal, constraints = inEq,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
print( summary( mlSannInEqCand ), digits = 2 )
all.equal( mlSannInEqCand[-c(2,3,4)], mlSannInEq[-c(2,3,4)], tolerance = 1e-3 )
all.equal( mlSannInEqCand, mlSannInEq, tolerance = 1e-1 )

############### equality constraints ###############
eqCon <- list( eqA = A, eqB = 2.5 )

## NR method with equality constraints
mlCon <- maxLik( llf, start = startVal, constraints = eqCon )
print( mlCon )
print( summary( mlCon ), digits = 2 )
activePar( mlCon )
AIC( mlCon )
coef( mlCon )
condiNumber( mlCon )
round( hessian( mlCon ), 1 )
logLik( mlCon )
maximType( mlCon )
nIter( mlCon )
nParam( mlCon )
returnCode( mlCon )
returnMessage( mlCon )
round( vcov( mlCon ), 3 )
logLik( summary( mlCon ) )
mlConInd <- maxLik( llfInd, start = startVal, constraints = eqCon )
print( summary( mlConInd ), digits = 2 )
all.equal( mlCon[-4], mlConInd[-c(4,11)], tolerance = 1e-3 )
mlConInd[11]
nObs( mlConInd )

# with analytical gradients
mlgCon <- maxLik( llf, gf, start = startVal, constraints = eqCon )
print( summary( mlgCon ), digits = 2 )
all.equal( mlCon[ -c(2,3,4,5,6,7,9,11) ], mlgCon[ -c(2,3,4,5,6,7,9,11) ], 
   tolerance = 1e-3 )
all.equal( mlCon[ -c( 5, 6, 7, 9, 11 ) ], mlgCon[ -c( 5, 6, 7, 9, 11 ) ], 
   tolerance = 1e-1 )
mlgConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon )
all.equal( mlConInd[ -c(2,3,4,5,6,7,9,11,12) ], mlgConInd[ -c(2,3,4,5,6,7,9,11,12) ],
   tolerance = 1e-3 )
all.equal( mlConInd[ -c(5,6,7,9,12) ], mlgConInd[ -c(5,6,7,9,12) ],
   tolerance = 1e-1 )
all.equal( mlgCon[], mlgConInd[-11], tolerance = 1e-3 )
mlgConInd[11]

# with analytical gradients as attribute
mlGCon <- maxLik( llfGrad, start = startVal, constraints = eqCon )
all.equal( mlGCon, mlgCon, tolerance = 1e-3 )
all.equal( mlGCon[-c(2,3,4,5,6,7,9,11)], mlCon[-c(2,3,4,5,6,7,9,11)], 
   tolerance = 1e-3 )
all.equal( mlGCon[-c(5,6,7,9,11)], mlCon[-c(5,6,7,9,11)], 
   tolerance = 1e-1 )

# with analytical gradients and Hessians
mlghCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon )
all.equal( mlgCon, mlghCon, tolerance = 1e-3 )

# with analytical gradients and Hessians as attributes
mlGHCon <- maxLik( llfGradHess, start = startVal, constraints = eqCon )
all.equal( mlGHCon, mlghCon, tolerance = 1e-3 )
all.equal( mlGHCon[-c(2,3,4,5,6,7,9,11)], mlCon[-c(2,3,4,5,6,7,9,11)], 
   tolerance = 1e-3 )
all.equal( mlGHCon[-c(5,6,7,9,11)], mlCon[-c(5,6,7,9,11)], 
   tolerance = 1e-1 )


## BHHH method with equality constraints
mlBhhhCon <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
print( mlBhhhCon )
print( summary( mlBhhhCon ), digits = 2 )
activePar( mlBhhhCon )
AIC( mlBhhhCon )
coef( mlBhhhCon )
condiNumber( mlBhhhCon )
round( hessian( mlBhhhCon ), 1 )
logLik( mlBhhhCon )
maximType( mlBhhhCon )
nIter( mlBhhhCon )
nParam( mlBhhhCon )
returnCode( mlBhhhCon )
returnMessage( mlBhhhCon )
round( vcov( mlBhhhCon ), 3 )
logLik( summary( mlBhhhCon ) )
all.equal( mlCon[ -c( 5, 6, 7, 9, 10 ) ], mlBhhhCon[ -c( 5, 6, 7, 9, 10, 11 ) ],
   tolerance = 5e-3 )
mlBhhhCon[11]
nObs( mlBhhhCon )

# with analytical gradients
mlgBhhhCon <- maxLik( llf, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
print( summary( mlgBhhhCon ), digits = 2 )
all.equal( mlBhhhCon[-c(2,3,4,5,6,7,9,11,12)], mlgBhhhCon[-c(2,3,4,5,6,7,9,11,12)],
   tolerance = 1e-3 )
all.equal( mlBhhhCon[-c(5,6,7,9,12)], mlgBhhhCon[-c(5,6,7,9,12)],
   tolerance = 1e-1 )
mlgBhhhConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlgBhhhCon, mlgBhhhConInd, tolerance = 1e-3 )

# with analytical gradients as attribute
mlGBhhhCon <- maxLik( llfGradInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
print( summary( mlGBhhhCon ), digits = 2 )
all.equal( mlGBhhhCon, mlgBhhhCon, tolerance = 1e-3 )
all.equal( mlGBhhhCon[-c(2,3,4,5,6,7,9,11,12)], mlBhhhCon[-c(2,3,4,5,6,7,9,11,12)],
   tolerance = 1e-3 )
all.equal( mlGBhhhCon[-c(5,6,7,9,12)], mlBhhhCon[-c(5,6,7,9,12)],
   tolerance = 1e-1 )

# with analytical gradients and unused Hessians
mlghBhhhCon <- maxLik( llf, gfInd, hf, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlgBhhhCon, mlghBhhhCon, tolerance = 1e-3 )

# with analytical gradients and unused Hessians as attributes
mlGHBhhhCon <- maxLik( llfGradHessInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlGHBhhhCon, mlghBhhhCon, tolerance = 1e-3 )
all.equal( mlGHBhhhCon, mlGBhhhCon, tolerance = 1e-3 )


## BFGS method with equality constraints
mlBfgsCon <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "BFGS" )
print( mlBfgsCon )
print( summary( mlBfgsCon ), digits = 2 )
activePar( mlBfgsCon )
AIC( mlBfgsCon )
coef( mlBfgsCon )
condiNumber( mlBfgsCon )
round( hessian( mlBfgsCon ), 1 )
logLik( mlBfgsCon )
maximType( mlBfgsCon )
nIter( mlBfgsCon )
nParam( mlBfgsCon )
returnCode( mlBfgsCon )
returnMessage( mlBfgsCon )
round( vcov( mlBfgsCon ), 3 )
logLik( summary( mlBfgsCon ) )
all.equal( mlBfgsCon[ -c( 4, 5, 6, 9, 10 ) ], mlCon[ -c( 4, 5, 6, 9, 10 ) ],
   tolerance = 1e-3 )
mlBfgsConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
print( summary( mlBfgsConInd ), digits = 2 )
all.equal( mlBfgsCon[-c(4,9)], mlBfgsConInd[-c(4,9,12)], tolerance = 1e-3 )
mlBfgsConInd[12]
nObs( mlBfgsConInd )

# with analytical gradients
mlgBfgsCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "BFGS" )
print( summary( mlgBfgsCon ), digits = 2 )
all.equal( mlBfgsCon[-c(3,4,9,11)], mlgBfgsCon[-c(3,4,9,11)], tolerance = 1e-2 )
mlgBfgsConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon[], mlgBfgsConInd[-12], tolerance = 1e-3 )
mlgBfgsConInd[12]

# with analytical gradients and unused Hessians
mlghBfgsCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon, mlghBfgsCon, tolerance = 1e-3 )

## NM method with equality constraints
mlNmCon <- maxLik( llf, start = startVal, constraints = eqCon, method = "NM", SUMTTol=0)
print( mlNmCon )
print( summary( mlNmCon ), digits = 2 )
activePar( mlNmCon )
AIC( mlNmCon )
coef( mlNmCon )
condiNumber( mlNmCon )
round( hessian( mlNmCon ), 1 )
logLik( mlNmCon )
maximType( mlNmCon )
nIter( mlNmCon )
nParam( mlNmCon )
returnCode( mlNmCon )
returnMessage( mlNmCon )
round( vcov( mlNmCon ), 3 )
logLik( summary( mlNmCon ) )
all.equal( mlNmCon[ -c( 4, 5, 6, 9, 10 ) ], mlCon[ -c( 4, 5, 6, 9, 10 ) ],
   tolerance = 1e-3 )
mlNmConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
print( summary( mlNmConInd ), digits = 2 )
all.equal( mlNmCon[], mlNmConInd[-12], tolerance = 1e-3 )
mlNmConInd[12]
nObs( mlNmConInd )

# with unused analytical gradients
mlgNmCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlNmCon, mlgNmCon, tolerance = 1e-3 )
mlgNmConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlgNmCon[], mlgNmConInd[-12], tolerance = 1e-3 )
mlgNmConInd[12]

# with unused analytical gradients and Hessians
mlghNmCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlgNmCon, mlghNmCon, tolerance = 1e-3 )

## SANN method with equality constraints
mlSannCon <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", SUMTTol=0)
print( mlSannCon )
print( summary( mlSannCon ), digits = 2 )
activePar( mlSannCon )
AIC( mlSannCon )
coef( mlSannCon )
condiNumber( mlSannCon )
round( hessian( mlSannCon ), 1 )
logLik( mlSannCon )
maximType( mlSannCon )
nIter( mlSannCon )
nParam( mlSannCon )
returnCode( mlSannCon )
returnMessage( mlSannCon )
round( vcov( mlSannCon ), 3 )
logLik( summary( mlSannCon ) )
all.equal( mlSannCon[ -c(2,3,4,5,6,9,10,11) ], mlBfgsCon[ -c(2,3,4,5,6,9,10,11) ],
   tolerance = 1e-3 )
all.equal( mlSannCon[ -c(3,4,5,6,9,10,11) ], mlBfgsCon[ -c(3,4,5,6,9,10,11) ],
   tolerance = 1e-2 )

# with unused analytical gradients
mlgSannCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "SANN", SUMTTol=0)
all.equal( mlSannCon, mlgSannCon, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSannConCand <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
print( summary( mlSannConCand ), digits = 2 )
all.equal( mlSannConCand[-c(1,2,3,4,11)], mlSannCon[-c(1,2,3,4,11)], 
   tolerance = 1e-3 )
all.equal( mlSannConCand[-c(2,3,4,11)], mlSannCon[-c(2,3,4,11)], 
   tolerance = 1e-1 )


## test for method "estfun"
library( sandwich )
try( estfun( ml ) )
estfun( mlInd )[ 1:5, ]
estfun( mlgInd )[ 1:5, ]
estfun( mlBHHH )[ 1:5, ]
estfun( mlgBHHH )[ 1:5, ]
estfun( mlIndBFGS )[ 1:5, ]
estfun( mlgIndBFGS )[ 1:5, ]
estfun( mlIndNM )[ 1:5, ]
estfun( mlgIndNM )[ 1:5, ]
estfun( mlIndSANN )[ 1:5, ]
estfun( mlgIndSANN )[ 1:5, ]
estfun( mlIndFix )[ 1:5, ]
estfun( mlgIndFix )[ 1:5, ]
estfun( mlFixBHHH )[ 1:5, ]
estfun( mlgFixBHHH )[ 1:5, ]
estfun( mlIndFixBfgs )[ 1:5, ]
estfun( mlgIndFixBfgs )[ 1:5, ]
estfun( mlIndFixNm )[ 1:5, ]
estfun( mlgIndFixNm )[ 1:5, ]
estfun( mlIndFixSann )[ 1:5, ]
estfun( mlgIndFixSann )[ 1:5, ]
estfun( mlBfgsInEqInd )[ 1:5, ]
estfun( mlgBfgsInEqInd )[ 1:5, ]
estfun( mlNmInEqInd )[ 1:5, ]
estfun( mlConInd )[ 1:5, ]
estfun( mlgConInd )[ 1:5, ]
estfun( mlBhhhCon )[ 1:5, ]
estfun( mlgBhhhCon )[ 1:5, ]
estfun( mlBfgsConInd )[ 1:5, ]
estfun( mlgBfgsConInd )[ 1:5, ]
estfun( mlNmConInd )[ 1:5, ]
estfun( mlgNmConInd )[ 1:5, ]


## test for method "bread"
try( bread( ml ) )
round( bread( mlInd ), 2 )
round( bread( mlgInd ), 2 )
round( bread( mlBHHH ), 2 )
round( bread( mlgBHHH ), 2 )
round( bread( mlIndBFGS ), 2 )
round( bread( mlgIndBFGS ), 2 )
round( bread( mlIndNM ), 2 )
round( bread( mlgIndNM ), 2 )
round( bread( mlIndSANN ), 2 )
round( bread( mlgIndSANN ), 2 )
round( bread( mlIndFix ), 2 )
round( bread( mlgIndFix ), 2 )
round( bread( mlFixBHHH ), 2 )
round( bread( mlgFixBHHH ), 2 )
round( bread( mlIndFixBfgs ), 2 )
round( bread( mlgIndFixBfgs ), 2 )
round( bread( mlIndFixNm ), 2 )
round( bread( mlgIndFixNm ), 2 )
round( bread( mlIndFixSann ), 2 )
round( bread( mlgIndFixSann ), 2 )
round( bread( mlBfgsInEqInd ), 2 )
round( bread( mlgBfgsInEqInd ), 2 )
round( bread( mlNmInEqInd ), 2 )
round( bread( mlConInd ), 2 )
round( bread( mlgConInd ), 2 )
round( bread( mlBhhhCon ), 2 )
round( bread( mlgBhhhCon ), 2 )
round( bread( mlBfgsConInd ), 2 )
round( bread( mlgBfgsConInd ), 2 )
round( bread( mlNmConInd ), 2 )
round( bread( mlgNmConInd ), 2 )


## test for method "sandwich"
try( sandwich( ml ) )
printSandwich <- function( x ) {
   print( round( sandwich( x ), 2 ) )
   tmp <- all.equal( sandwich( x ), vcov( x ) )
   if( isTRUE( tmp ) ) {
      print( tmp )
   }
}
printSandwich( mlInd )
printSandwich( mlgInd )
printSandwich( mlBHHH )
printSandwich( mlgBHHH )
printSandwich( mlIndBFGS )
printSandwich( mlgIndBFGS )
printSandwich( mlIndNM )
printSandwich( mlgIndNM )
printSandwich( mlIndSANN )
printSandwich( mlgIndSANN )
printSandwich( mlIndFix )
printSandwich( mlgIndFix )
printSandwich( mlFixBHHH )
printSandwich( mlgFixBHHH )
printSandwich( mlIndFixBfgs )
printSandwich( mlgIndFixBfgs )
printSandwich( mlIndFixNm )
printSandwich( mlgIndFixNm )
printSandwich( mlIndFixSann )
printSandwich( mlgIndFixSann )
printSandwich( mlBfgsInEqInd )
printSandwich( mlgBfgsInEqInd )
printSandwich( mlNmInEqInd )
printSandwich( mlConInd )
printSandwich( mlgConInd )
printSandwich( mlBhhhCon )
printSandwich( mlgBhhhCon )
printSandwich( mlBfgsConInd )
printSandwich( mlgBfgsConInd )
printSandwich( mlNmConInd )
printSandwich( mlgNmConInd )
