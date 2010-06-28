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

## log likelihood function with gradients as attributes
llfGrad <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
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
summary( ml )
activePar( ml )
AIC( ml )
coef( ml )
condiNumber( ml )
hessian( ml )
logLik( ml )
maximType( ml )
nIter( ml )
try( nObs( ml ) )
nParam( ml )
returnCode( ml )
returnMessage( ml )
vcov( ml )
logLik( summary( ml ) )
mlInd <- maxLik( llfInd, start = startVal )
summary( mlInd )
all.equal( ml[ ], mlInd[ -11 ] )
mlInd[ 11 ]
nObs( mlInd )

# with analytical gradients
mlg <- maxLik( llf, gf, start = startVal )
summary( mlg )
all.equal( ml, mlg )
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
all.equal( mlInd, mlgInd )
all.equal( mlg[ ], mlgInd[ -11 ] )
mlgInd[ 11 ]

# with analytical gradients as attribute
mlG <- maxLik( llfGrad, start = startVal )
all.equal( mlG, mlg, tolerance = 1e-10 )
all.equal( mlG$gradient, gf( coef( mlG ) ), check.attributes = FALSE )
mlGInd <- maxLik( llfGradInd, start = startVal )
all.equal( mlGInd, mlgInd, tolerance = 1e-10 )
all.equal( mlGInd$gradient, colSums( gfInd( coef( mlGInd ) ) ),
   check.attributes = FALSE )
all.equal( mlGInd$gradientObs, gfInd( coef( mlGInd ) ),
   check.attributes = FALSE )

# with analytical gradients and Hessians
mlgh <- maxLik( llf, gf, hf, start = startVal )
all.equal( mlg, mlgh )

# with analytical gradients and Hessian as attribute
mlGH <- maxLik( llfGradHess, start = startVal )
all.equal( mlGH, mlgh, tolerance = 1e-10 )


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
nParam( mlBHHH )
returnCode( mlBHHH )
returnMessage( mlBHHH )
vcov( mlBHHH )
logLik( summary( mlBHHH ) )
all.equal( ml[ ], mlBHHH[ -11 ] )
mlBHHH[ 11 ]
nObs( mlBHHH )

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
all.equal( mlg[ ], mlgBHHH[ -11 ] )
mlgBHHH[ 11 ]
mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlgBHHH2 )

# with analytical gradients as attribute
mlGBHHH <- maxLik( llfGradInd, start = startVal, method = "BHHH" )
all.equal( mlGBHHH, mlgBHHH, tolerance = 1e-10 )

# with unused Hessian
mlghBHHH <- maxLik( llfInd, gfInd, hf, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlghBHHH )

# with unused Hessian as attribute
mlGHBHHH <- maxLik( llfGradHessInd, start = startVal, method = "BHHH" )
all.equal( mlGHBHHH, mlghBHHH, tolerance = 1e-10 )


### BFGS-YC method
mlgBFGSYC <- maxLik( llf, gf, start = startVal, method = "bfgsyc" , print.level=1)
print( mlgBFGSYC )
summary(mlgBFGSYC)
activePar( mlgBFGSYC )
AIC( mlgBFGSYC )
coef( mlgBFGSYC )
condiNumber( mlgBFGSYC )
hessian( mlgBFGSYC )
logLik( mlgBFGSYC )
maximType( mlgBFGSYC )
nIter( mlgBFGSYC )
nParam( mlgBFGSYC )
returnCode( mlgBFGSYC )
returnMessage( mlgBFGSYC )
vcov( mlgBFGSYC )
logLik( summary( mlgBFGSYC ) )
all.equal( mlgBFGSYC[ -c( 5, 6, 9, 10 ) ], mlg[ -c( 5, 6, 9, 10 ) ] )

# with analytical gradients as attribute
mlGBFGSYC <- maxLik( llfGrad, start = startVal, method = "bfgsyc" , print.level=1)
all.equal( mlGBFGSYC, mlgBFGSYC )


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
nParam( mlBFGS )
returnCode( mlBFGS )
returnMessage( mlBFGS )
vcov( mlBFGS )
logLik( summary( mlBFGS ) )
all.equal( ml, mlBFGS )
# with individual log likelihood values
mlIndBFGS <- maxLik( llfInd, start = startVal, method = "BFGS" )
summary( mlIndBFGS )
all.equal( mlBFGS[], mlIndBFGS[-12] )
mlIndBFGS[12]
nObs( mlIndBFGS )

# with analytical gradients
mlgBFGS <- maxLik( llf, gf, start = startVal, method = "BFGS" )
summary( mlgBFGS )
all.equal( mlBFGS, mlgBFGS )
all.equal( mlg, mlgBFGS )
mlgIndBFGS <- maxLik( llfInd, gfInd, start = startVal, method = "BFGS" )
all.equal( mlgBFGS[], mlgIndBFGS[-12] )
mlgIndBFGS[12]

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
nParam( mlNM )
returnCode( mlNM )
returnMessage( mlNM )
vcov( mlNM )
logLik( summary( mlNM ) )
all.equal( ml, mlNM )
# with individual log likelihood values
mlIndNM <- maxLik( llfInd, start = startVal, method = "NM" )
summary( mlIndNM )
all.equal( mlNM[], mlIndNM[-12] )
mlIndNM[12]
nObs( mlIndNM )

# with unused analytical gradients
mlgNM <- maxLik( llf, gf, start = startVal, method = "NM" )
summary( mlgNM )
all.equal( mlNM, mlgNM )
# with individual log likelihood values and gradients
mlgIndNM <- maxLik( llfInd, gfInd, start = startVal, method = "NM" )
summary( mlgIndNM )
all.equal( mlgNM[], mlgIndNM[-12] )
mlgIndNM[12]

# with unused analytical gradients and Hessian
mlghNM <- maxLik( llf, gf, hf, start = startVal, method = "NM" )
all.equal( mlgNM, mlghNM )


## SANN method
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
nParam( mlSANN )
returnCode( mlSANN )
returnMessage( mlSANN )
vcov( mlSANN )
logLik( summary( mlSANN ) )
all.equal( ml, mlSANN )
# with individual log likelihood values
mlIndSANN <- maxLik( llfInd, start = startVal, method = "SANN" )
summary( mlIndSANN )
all.equal( mlSANN[], mlIndSANN[-12] )
mlIndSANN[12]
nObs( mlIndSANN )

# with unused analytical gradients
mlgSANN <- maxLik( llf, gf, start = startVal, method = "SANN" )
summary( mlgSANN )
all.equal( mlSANN, mlgSANN )
# with individual log likelihood values and gradients
mlgIndSANN <- maxLik( llfInd, gfInd, start = startVal, method = "SANN" )
summary( mlgIndSANN )
all.equal( mlgSANN[], mlgIndSANN[-12] )
mlgIndSANN[12]

# with unused analytical gradients and Hessian
mlghSANN <- maxLik( llf, gf, hf, start = startVal, method = "SANN" )
all.equal( mlgSANN, mlghSANN )

# with a user-specified function to generate a new candidate point
mlSANNCand <- maxLik( llf, start = startVal, method = "SANN",
   cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSANNCand )
all.equal( mlSANNCand, mlSANN )

############### with fixed parameters ###############
# start values
startValFix <- c( mu = 1, sigma = 1 )

# fix mu (the mean ) at its start value
isFixed <- c( TRUE, FALSE )

## NR method with fixed parameters
mlFix <- maxLik( llf, start = startValFix, activePar = !isFixed )
mlFix1 <- maxLik( llf, start = startValFix, activePar = 2 )
all.equal( mlFix, mlFix1 )
mlFix2 <- maxLik( llf, start = startValFix, fixed = isFixed )
all.equal( mlFix, mlFix2 )
mlFix3 <- maxLik( llf, start = startValFix, fixed = "mu" )
all.equal( mlFix, mlFix3 )
mlFix4 <- maxLik( llf, start = startValFix, fixed = 1 )
all.equal( mlFix, mlFix4 )
print( mlFix )
summary( mlFix )
activePar( mlFix )
AIC( mlFix )
coef( mlFix )
condiNumber( mlFix )
hessian( mlFix )
logLik( mlFix )
maximType( mlFix )
nIter( mlFix )
nParam( mlFix )
returnCode( mlFix )
returnMessage( mlFix )
vcov( mlFix )
logLik( summary( mlFix ) )
mlIndFix <- maxLik( llfInd, start = startValFix, activePar = !isFixed )
mlIndFix1 <- maxLik( llfInd, start = startValFix, activePar = 2 )
all.equal( mlIndFix, mlIndFix1 )
mlIndFix2 <- maxLik( llfInd, start = startValFix, fixed = isFixed )
all.equal( mlIndFix, mlIndFix2 )
mlIndFix3 <- maxLik( llfInd, start = startValFix, fixed = "mu" )
all.equal( mlIndFix, mlIndFix3 )
mlIndFix4 <- maxLik( llfInd, start = startValFix, fixed = 1 )
all.equal( mlIndFix, mlIndFix4 )
summary( mlIndFix )
all.equal( mlFix[ ], mlIndFix[ -11 ] )
mlFix[[3]]
mlIndFix[[3]]
mlIndFix[ 11 ]
nObs( mlIndFix )

# with analytical gradients
mlgFix <- maxLik( llf, gf, start = startValFix, activePar = !isFixed )
mlgFix1 <- maxLik( llf, gf, start = startValFix, activePar = 2 )
all.equal( mlgFix, mlgFix1 )
mlgFix2 <- maxLik( llf, gf, start = startValFix, fixed = isFixed )
all.equal( mlgFix, mlgFix2 )
summary( mlgFix )
all.equal( mlFix, mlgFix )
mlFix[[3]]
mlgFix[[3]]
mlFix[[4]]
mlgFix[[4]]
mlgIndFix <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed )
all.equal( mlIndFix, mlgIndFix )
mlIndFix[[3]]
mlgIndFix[[3]]
mlIndFix[[4]]
mlgIndFix[[4]]
all.equal( mlgFix[ ], mlgIndFix[ -11 ] )
mlgIndFix[ 11 ]

# with analytical gradients and Hessians
mlghFix <- maxLik( llf, gf, hf, start = startValFix, activePar = !isFixed )
all.equal( mlgFix, mlghFix )
mlgFix[[4]]
mlghFix[[4]]

## BHHH method with fixed parameters
mlFixBHHH <- maxLik( llfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
mlFixBHHH1 <- maxLik( llfInd, start = startValFix, activePar = 2,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH1 )
mlFixBHHH2 <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH2 )
mlFixBHHH3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH3 )
mlFixBHHH4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "BHHH" )
all.equal( mlFixBHHH, mlFixBHHH4 )
print( mlFixBHHH )
summary( mlFixBHHH )
activePar( mlFixBHHH )
AIC( mlFixBHHH )
coef( mlFixBHHH )
condiNumber( mlFixBHHH )
hessian( mlFixBHHH )
logLik( mlFixBHHH )
maximType( mlFixBHHH )
nIter( mlFixBHHH )
nParam( mlFixBHHH )
returnCode( mlFixBHHH )
returnMessage( mlFixBHHH )
vcov( mlFixBHHH )
logLik( summary( mlFixBHHH ) )
all.equal( mlFix[ -c( 5, 6, 9, 10 ) ], mlFixBHHH[ -c( 5, 6, 9, 10, 11 ) ] )
mlFix[[ 3 ]]
mlFixBHHH[[ 3 ]]
mlFix[[ 4 ]]
mlFixBHHH[[ 4 ]]
mlFixBHHH[ 11 ]
nObs( mlFixBHHH )

# with analytical gradients
mlgFixBHHH <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
mlgFixBHHH1 <- maxLik( llfInd, gfInd, start = startValFix, activePar = 2,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH1 )
mlgFixBHHH2 <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH2 )
mlgFixBHHH3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH3 )
mlgFixBHHH4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlgFixBHHH4 )
summary( mlgFixBHHH )
all.equal( mlFixBHHH, mlgFixBHHH )
mlgFixBHHH2 <- maxLik( llf, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH")
all.equal( mlgFixBHHH, mlgFixBHHH2 )

# with unused Hessians
mlghFixBHHH <- maxLik( llfInd, gfInd, hf, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlghFixBHHH )

## BFGS method with fixed parameters
mlFixBfgs <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
mlFixBfgs3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlFixBfgs, mlFixBfgs3 )
mlFixBfgs4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlFixBfgs, mlFixBfgs4 )
print( mlFixBfgs )
summary( mlFixBfgs )
activePar( mlFixBfgs )
AIC( mlFixBfgs )
coef( mlFixBfgs )
condiNumber( mlFixBfgs )
hessian( mlFixBfgs )
logLik( mlFixBfgs )
maximType( mlFixBfgs )
nIter( mlFixBfgs )
nParam( mlFixBfgs )
returnCode( mlFixBfgs )
returnMessage( mlFixBfgs )
vcov( mlFixBfgs )
logLik( summary( mlFixBfgs ) )
all.equal( mlghFix[ -c( 5, 6, 9, 10 ) ], mlFixBfgs[ -c( 5, 6, 9, 10, 11 ) ] )
mlIndFixBfgs <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "BFGS" )
all.equal( mlFixBfgs[ -9 ], mlIndFixBfgs[ -c(9,12) ] )
mlIndFixBfgs[ 12 ]
mlIndFixBfgs3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlIndFixBfgs, mlIndFixBfgs3 )
mlIndFixBfgs4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlIndFixBfgs, mlIndFixBfgs4 )
nObs( mlIndFixBfgs )

# with analytical gradients
mlgFixBfgs <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
mlgFixBfgs3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlgFixBfgs, mlgFixBfgs3 )
mlgFixBfgs4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlgFixBfgs, mlgFixBfgs4 )
summary( mlgFixBfgs )
all.equal( mlFixBfgs[ -9 ], mlgFixBfgs[ -9 ] )
mlgIndFixBfgs <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "BFGS")
all.equal( mlgFixBfgs[ ], mlgIndFixBfgs[ -12 ] )
mlgIndFixBfgs[ 12 ]
mlgIndFixBfgs3 <- maxLik( llfInd, gfInd, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlgIndFixBfgs, mlgIndFixBfgs3 )
mlgIndFixBfgs4 <- maxLik( llfInd, gfInd, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlgIndFixBfgs, mlgIndFixBfgs4 )

# with unused Hessians
mlghFixBfgs <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "BFGS" )
all.equal( mlgFixBfgs, mlghFixBfgs )
mlghFixBfgs3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
   method = "BFGS" )
all.equal( mlghFixBfgs, mlghFixBfgs3 )
mlghFixBfgs4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
   method = "BFGS" )
all.equal( mlghFixBfgs, mlghFixBfgs4 )

## NM method with fixed parameters
mlFixNm <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "NM" )
mlFixNm3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlFixNm, mlFixNm3 )
mlFixNm4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlFixNm, mlFixNm4 )
print( mlFixNm )
summary( mlFixNm )
activePar( mlFixNm )
AIC( mlFixNm )
coef( mlFixNm )
condiNumber( mlFixNm )
hessian( mlFixNm )
logLik( mlFixNm )
maximType( mlFixNm )
nIter( mlFixNm )
nParam( mlFixNm )
returnCode( mlFixNm )
returnMessage( mlFixNm )
vcov( mlFixNm )
logLik( summary( mlFixNm ) )
all.equal( mlFixBfgs[ -c( 9, 10 ) ], mlFixNm[ -c( 9, 10 ) ] )
mlIndFixNm <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "NM" )
all.equal( mlFixNm[ ], mlIndFixNm[ -12 ] )
mlIndFixNm[ 12 ]
mlIndFixNm3 <- maxLik( llfInd, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlIndFixNm, mlIndFixNm3 )
mlIndFixNm4 <- maxLik( llfInd, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlIndFixNm, mlIndFixNm4 )
nObs( mlIndFixNm )

# with analytical gradients
mlgFixNm <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "NM" )
mlgFixNm3 <- maxLik( llf, gf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlgFixNm, mlgFixNm3 )
mlgFixNm4 <- maxLik( llf, gf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlgFixNm, mlgFixNm4 )
summary( mlgFixNm )
all.equal( mlFixNm, mlgFixNm )
mlgIndFixNm <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "NM")
all.equal( mlgFixNm[ ], mlgIndFixNm[ -12 ] )
mlgIndFixNm[ 12 ]

# with unused Hessians
mlghFixNm <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "NM" )
all.equal( mlgFixNm, mlghFixNm )
mlghFixNm3 <- maxLik( llf, gf, hf, start = startValFix, fixed = "mu",
   method = "NM" )
all.equal( mlghFixNm, mlghFixNm3 )
mlghFixNm4 <- maxLik( llf, gf, hf, start = startValFix, fixed = 1,
   method = "NM" )
all.equal( mlghFixNm, mlghFixNm4 )

## SANN method with fixed parameters
mlFixSann <- maxLik( llf, start = startValFix, fixed = isFixed,
   method = "SANN" )
mlFixSann3 <- maxLik( llf, start = startValFix, fixed = "mu",
   method = "SANN" )
all.equal( mlFixSann, mlFixSann3 )
mlFixSann4 <- maxLik( llf, start = startValFix, fixed = 1,
   method = "SANN" )
all.equal( mlFixSann, mlFixSann4 )
print( mlFixSann )
summary( mlFixSann )
activePar( mlFixSann )
AIC( mlFixSann )
coef( mlFixSann )
condiNumber( mlFixSann )
hessian( mlFixSann )
logLik( mlFixSann )
maximType( mlFixSann )
nIter( mlFixSann )
nParam( mlFixSann )
returnCode( mlFixSann )
returnMessage( mlFixSann )
vcov( mlFixSann )
logLik( summary( mlFixSann ) )
all.equal( mlFixBfgs[ -c( 9, 10 ) ], mlFixSann[ -c( 9, 10 ) ] )
mlIndFixSann <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "SANN" )
all.equal( mlFixSann[ ], mlIndFixSann[ -12 ] )
mlIndFixSann[ 12 ]
nObs( mlIndFixSann )

# with analytical gradients
mlgFixSann <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "SANN" )
summary( mlgFixSann )
all.equal( mlFixSann, mlgFixSann )
mlgIndFixSann <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "SANN")
all.equal( mlgFixSann[ ], mlgIndFixSann[ -12 ] )
mlgIndFixSann[ 12 ]

# with unused Hessians
mlghFixSann <- maxLik( llf, gf, hf, start = startValFix, fixed = isFixed,
   method = "SANN" )
all.equal( mlgFixSann, mlghFixSann )


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
summary( mlBfgsInEq )
activePar( mlBfgsInEq )
AIC( mlBfgsInEq )
coef( mlBfgsInEq )
condiNumber( mlBfgsInEq )
hessian( mlBfgsInEq )
logLik( mlBfgsInEq )
maximType( mlBfgsInEq )
nIter( mlBfgsInEq )
nParam( mlBfgsInEq )
returnCode( mlBfgsInEq )
returnMessage( mlBfgsInEq )
vcov( mlBfgsInEq )
logLik( summary( mlBfgsInEq ) )
mlBfgsInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
summary( mlBfgsInEqInd )
all.equal( mlBfgsInEq[ ], mlBfgsInEqInd[ -12 ] )
mlBfgsInEqInd[ 12 ]
nObs( mlBfgsInEqInd )

# with analytical gradients
mlgBfgsInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlBfgsInEq, mlgBfgsInEq )
mlgBfgsInEqInd <- maxLik( llfInd, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd[ -12 ], mlgBfgsInEq[ ] )
mlgBfgsInEqInd[ 12 ]
mlgBfgsInEqInd2 <- maxLik( llf, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd, mlgBfgsInEqInd2 )

# with unused Hessian
mlghBfgsInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEq, mlghBfgsInEq )

## NM method with inequality constraints
mlNmInEq <- maxLik( llf, start = startVal, constraints = inEq, method = "NM" )
print( mlNmInEq )
summary( mlNmInEq )
activePar( mlNmInEq )
AIC( mlNmInEq )
coef( mlNmInEq )
condiNumber( mlNmInEq )
hessian( mlNmInEq )
logLik( mlNmInEq )
maximType( mlNmInEq )
nIter( mlNmInEq )
nParam( mlNmInEq )
returnCode( mlNmInEq )
returnMessage( mlNmInEq )
vcov( mlNmInEq )
logLik( summary( mlNmInEq ) )
all.equal( mlBfgsInEq, mlNmInEq )
mlNmInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "NM" )
summary( mlNmInEqInd )
all.equal( mlNmInEq[ ], mlNmInEqInd[ -12 ] )
mlNmInEqInd[ 12 ]
nObs( mlNmInEqInd )

# with unused analytical gradients
mlgNmInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlNmInEq, mlgNmInEq )

# with unused analytical gradients and Hessians
mlghNmInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlgNmInEq, mlghNmInEq )

## SANN method with inequality constraints
mlSannInEq <- maxLik( llf, start = startVal, constraints = inEq,
   method = "SANN" )
print( mlSannInEq )
summary( mlSannInEq )
activePar( mlSannInEq )
AIC( mlSannInEq )
coef( mlSannInEq )
condiNumber( mlSannInEq )
hessian( mlSannInEq )
logLik( mlSannInEq )
maximType( mlSannInEq )
nIter( mlSannInEq )
nParam( mlSannInEq )
returnCode( mlSannInEq )
returnMessage( mlSannInEq )
vcov( mlSannInEq )
logLik( summary( mlSannInEq ) )
all.equal( mlBfgsInEq, mlSannInEq )

# with unused analytical gradients
mlgSannInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "SANN" )
all.equal( mlSannInEq, mlgSannInEq )

# with a user-specified function to generate a new candidate point
mlSannInEqCand <- maxLik( llf, start = startVal, constraints = inEq,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSannInEqCand )
all.equal( mlSannInEqCand, mlSannInEq )

############### equality constraints ###############
eqCon <- list( eqA = A, eqB = 2.5 )

## NR method with equality constraints
mlCon <- maxLik( llf, start = startVal, constraints = eqCon )
print( mlCon )
summary( mlCon )
activePar( mlCon )
AIC( mlCon )
coef( mlCon )
condiNumber( mlCon )
hessian( mlCon )
logLik( mlCon )
maximType( mlCon )
nIter( mlCon )
nParam( mlCon )
returnCode( mlCon )
returnMessage( mlCon )
vcov( mlCon )
logLik( summary( mlCon ) )
mlConInd <- maxLik( llfInd, start = startVal, constraints = eqCon )
summary( mlConInd )
all.equal( mlCon[], mlConInd[-12] )
mlConInd[12]
nObs( mlConInd )

# with analytical gradients
mlgCon <- maxLik( llf, gf, start = startVal, constraints = eqCon )
summary( mlgCon )
all.equal( mlCon[ -c( 5, 6, 7, 9 ) ], mlgCon[ -c( 5, 6, 7, 9 ) ] )
mlgConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon )
all.equal( mlConInd, mlgConInd )
all.equal( mlgCon[], mlgConInd[-12] )
mlgConInd[12]

# with analytical gradients and Hessians
mlghCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon )
all.equal( mlgCon, mlghCon )

## BHHH method with equality constraints
mlBhhhCon <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
print( mlBhhhCon )
summary( mlBhhhCon )
activePar( mlBhhhCon )
AIC( mlBhhhCon )
coef( mlBhhhCon )
condiNumber( mlBhhhCon )
hessian( mlBhhhCon )
logLik( mlBhhhCon )
maximType( mlBhhhCon )
nIter( mlBhhhCon )
nParam( mlBhhhCon )
returnCode( mlBhhhCon )
returnMessage( mlBhhhCon )
vcov( mlBhhhCon )
logLik( summary( mlBhhhCon ) )
all.equal( mlCon[ -c( 5, 6, 7, 9, 10 ) ], mlBhhhCon[ -c( 5, 6, 7, 9, 10, 12 ) ] )
mlBhhhCon[12]
nObs( mlBhhhCon )

# with analytical gradients
mlgBhhhCon <- maxLik( llf, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
summary( mlgBhhhCon )
all.equal( mlBhhhCon, mlgBhhhCon )
mlgBhhhConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlgBhhhCon, mlgBhhhConInd )

# with analytical gradients and unused Hessians
mlghBhhhCon <- maxLik( llf, gfInd, hf, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlgBhhhCon, mlghBhhhCon )

## BFGS method with equality constraints
mlBfgsCon <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "BFGS" )
print( mlBfgsCon )
summary( mlBfgsCon )
activePar( mlBfgsCon )
AIC( mlBfgsCon )
coef( mlBfgsCon )
condiNumber( mlBfgsCon )
hessian( mlBfgsCon )
logLik( mlBfgsCon )
maximType( mlBfgsCon )
nIter( mlBfgsCon )
nParam( mlBfgsCon )
returnCode( mlBfgsCon )
returnMessage( mlBfgsCon )
vcov( mlBfgsCon )
logLik( summary( mlBfgsCon ) )
all.equal( mlBfgsCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ] )
mlBfgsConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlBfgsConInd )
all.equal( mlBfgsCon[], mlBfgsConInd[-12] )
mlBfgsConInd[12]
nObs( mlBfgsConInd )

# with analytical gradients
mlgBfgsCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlgBfgsCon )
all.equal( mlBfgsCon, mlgBfgsCon )
mlgBfgsConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon[], mlgBfgsConInd[-12] )
mlgBfgsConInd[12]

# with analytical gradients and unused Hessians
mlghBfgsCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon, mlghBfgsCon )

## NM method with equality constraints
mlNmCon <- maxLik( llf, start = startVal, constraints = eqCon, method = "NM", SUMTTol=0)
print( mlNmCon )
summary( mlNmCon )
activePar( mlNmCon )
AIC( mlNmCon )
coef( mlNmCon )
condiNumber( mlNmCon )
hessian( mlNmCon )
logLik( mlNmCon )
maximType( mlNmCon )
nIter( mlNmCon )
nParam( mlNmCon )
returnCode( mlNmCon )
returnMessage( mlNmCon )
vcov( mlNmCon )
logLik( summary( mlNmCon ) )
all.equal( mlNmCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ] )
mlNmConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
summary( mlNmConInd )
all.equal( mlNmCon[], mlNmConInd[-12] )
mlNmConInd[12]
nObs( mlNmConInd )

# with unused analytical gradients
mlgNmCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlNmCon, mlgNmCon )
mlgNmConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlgNmCon[], mlgNmConInd[-12] )
mlgNmConInd[12]

# with unused analytical gradients and Hessians
mlghNmCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
all.equal( mlgNmCon, mlghNmCon )

## SANN method with equality constraints
mlSannCon <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", SUMTTol=0)
print( mlSannCon )
summary( mlSannCon )
activePar( mlSannCon )
AIC( mlSannCon )
coef( mlSannCon )
condiNumber( mlSannCon )
hessian( mlSannCon )
logLik( mlSannCon )
maximType( mlSannCon )
nIter( mlSannCon )
nParam( mlSannCon )
returnCode( mlSannCon )
returnMessage( mlSannCon )
vcov( mlSannCon )
logLik( summary( mlSannCon ) )
all.equal( mlSannCon[ -c( 5, 6, 9, 10 ) ], mlBfgsCon[ -c( 5, 6, 9, 10 ) ] )

# with unused analytical gradients
mlgSannCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "SANN", SUMTTol=0)
all.equal( mlSannCon, mlgSannCon )

# with a user-specified function to generate a new candidate point
mlSannConCand <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSannConCand )
all.equal( mlSannConCand, mlSannCon )


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
bread( mlInd )
bread( mlgInd )
bread( mlBHHH )
bread( mlgBHHH )
bread( mlIndBFGS )
bread( mlgIndBFGS )
bread( mlIndNM )
bread( mlgIndNM )
bread( mlIndSANN )
bread( mlgIndSANN )
bread( mlIndFix )
bread( mlgIndFix )
bread( mlFixBHHH )
bread( mlgFixBHHH )
bread( mlIndFixBfgs )
bread( mlgIndFixBfgs )
bread( mlIndFixNm )
bread( mlgIndFixNm )
bread( mlIndFixSann )
bread( mlgIndFixSann )
bread( mlBfgsInEqInd )
bread( mlgBfgsInEqInd )
bread( mlNmInEqInd )
bread( mlConInd )
bread( mlgConInd )
bread( mlBhhhCon )
bread( mlgBhhhCon )
bread( mlBfgsConInd )
bread( mlgBfgsConInd )
bread( mlNmConInd )
bread( mlgNmConInd )


## test for method "sandwich"
try( sandwich( ml ) )
printSandwich <- function( x ) {
   print( sandwich( x ) )
   print( all.equal( sandwich( x ), vcov( x ) ) )
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
