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
all.equal( ml[-3], mlInd[ c(-3,-11) ], tolerance = 1e-3 )
                           # 3  gradient, should be close to 0, but may vary enormously in relative terms
mlInd[[11]][sample(nrow(mlInd[[11]]), 10),]
                           # just print a sample of 10
nObs( mlInd )

# with analytical gradients
mlg <- maxLik( llf, gf, start = startVal )
summary( mlg )
all.equal( ml, mlg, tolerance = 1e-3 )
mlgInd <- maxLik( llfInd, gfInd, start = startVal )
all.equal( mlInd, mlgInd, tolerance = 1e-3 )
all.equal( mlg[ ], mlgInd[ -11 ], tolerance = 1e-3 )
mlgInd[ 11 ]

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
all.equal( ml[ ], mlBHHH[ -11 ], tolerance = 1e-3 )
mlBHHH[ 11 ]
nObs( mlBHHH )
# final Hessian = usual Hessian
mlBhhhH <- maxLik( llfInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlBhhhH, mlBHHH, tolerance = 1e-3 )
hessian( mlBhhhH ) 
summary( mlBhhhH ) 

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
all.equal( mlBHHH, mlgBHHH, tolerance = 1e-3 )
all.equal( mlg[ ], mlgBHHH[ -11 ], tolerance = 1e-3 )
mlgBHHH[ 11 ]
mlgBHHH2 <- maxLik( llf, gfInd, start = startVal, method = "BHHH" )
all.equal( mlgBHHH, mlgBHHH2, tolerance = 1e-3 )
# final Hessian = usual Hessian
mlgBhhhH <- maxLik( llf, gfInd, start = startVal, method = "BHHH", 
   finalHessian = TRUE )
all.equal( mlgBhhhH, mlBhhhH, tolerance = 1e-3 )
all.equal( mlgBhhhH, mlgBHHH, tolerance = 1e-3 )
hessian( mlgBhhhH ) 

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
all.equal( mlghBhhhH, mlghBHHH, tolerance = 1e-3 )
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
summary( mlBFGSYC )
activePar( mlBFGSYC )
AIC( mlBFGSYC )
coef( mlBFGSYC )
condiNumber( mlBFGSYC )
hessian( mlBFGSYC )
logLik( mlBFGSYC )
maximType( mlBFGSYC )
nIter( mlBFGSYC )
try( nObs( mlBFGSYC ) )
nParam( mlBFGSYC )
returnCode( mlBFGSYC )
returnMessage( mlBFGSYC )
vcov( mlBFGSYC )
logLik( summary( mlBFGSYC ) )
all.equal( ml[-c(5,6,9,10)], mlBFGSYC[-c(5,6,9,10)], tolerance = 1e-3 )
mlIndBFGSYC <- maxLik( llfInd, start = startVal, method = "BFGSR" )
summary( mlIndBFGSYC )
all.equal( mlBFGSYC[ -9 ], mlIndBFGSYC[ -c(9,11) ], tolerance = 1e-3 )
mlIndBFGSYC[ 11 ]
nObs( mlIndBFGSYC )

# with analytical gradients
mlgBFGSYC <- maxLik( llf, gf, start = startVal, method = "BFGSR" , print.level=1)
summary(mlgBFGSYC)
all.equal( mlBFGSYC, mlgBFGSYC, tolerance = 1e-3 )
mlgIndBFGSYC <- maxLik( llfInd, gfInd, start = startVal,
   method = "BFGSR" )
all.equal( mlIndBFGSYC, mlgIndBFGSYC, tolerance = 1e-3 )
all.equal( mlgBFGSYC[ -9 ], mlgIndBFGSYC[ -c(9,11) ], tolerance = 1e-3 )
mlgIndBFGSYC[ 11 ]

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
all.equal( ml, mlBFGS, tolerance = 1e-3 )
# with individual log likelihood values
mlIndBFGS <- maxLik( llfInd, start = startVal, method = "BFGS" )
summary( mlIndBFGS )
all.equal( mlBFGS[], mlIndBFGS[-12], tolerance = 1e-3 )
mlIndBFGS[12]
nObs( mlIndBFGS )

# with analytical gradients
mlgBFGS <- maxLik( llf, gf, start = startVal, method = "BFGS" )
summary( mlgBFGS )
all.equal( mlBFGS, mlgBFGS, tolerance = 1e-3 )
all.equal( mlg, mlgBFGS, tolerance = 1e-3 )
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
all.equal( ml, mlNM, tolerance = 1e-3 )
# with individual log likelihood values
mlIndNM <- maxLik( llfInd, start = startVal, method = "NM" )
summary( mlIndNM )
all.equal( mlNM[], mlIndNM[-12], tolerance = 1e-3 )
mlIndNM[12]
nObs( mlIndNM )

# with unused analytical gradients
mlgNM <- maxLik( llf, gf, start = startVal, method = "NM" )
summary( mlgNM )
all.equal( mlNM, mlgNM, tolerance = 1e-3 )
# with individual log likelihood values and gradients
mlgIndNM <- maxLik( llfInd, gfInd, start = startVal, method = "NM" )
summary( mlgIndNM )
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
all.equal( ml, mlSANN, tolerance = 1e-3 )
# with individual log likelihood values
mlIndSANN <- maxLik( llfInd, start = startVal, method = "SANN" )
summary( mlIndSANN )
all.equal( mlSANN[], mlIndSANN[-12], tolerance = 1e-3 )
mlIndSANN[12]
nObs( mlIndSANN )

# with unused analytical gradients
mlgSANN <- maxLik( llf, gf, start = startVal, method = "SANN" )
summary( mlgSANN )
all.equal( mlSANN, mlgSANN, tolerance = 1e-3 )
# with individual log likelihood values and gradients
mlgIndSANN <- maxLik( llfInd, gfInd, start = startVal, method = "SANN" )
summary( mlgIndSANN )
all.equal( mlgSANN[], mlgIndSANN[-12], tolerance = 1e-3 )
mlgIndSANN[12]

# with unused analytical gradients and Hessian
mlghSANN <- maxLik( llf, gf, hf, start = startVal, method = "SANN" )
all.equal( mlgSANN, mlghSANN, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSANNCand <- maxLik( llf, start = startVal, method = "SANN",
   cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSANNCand )
all.equal( mlSANNCand, mlSANN, tolerance = 1e-3 )

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
all.equal( mlIndFix, mlIndFix1, tolerance = 1e-3 )
mlIndFix2 <- maxLik( llfInd, start = startValFix, fixed = isFixed )
all.equal( mlIndFix, mlIndFix2, tolerance = 1e-3 )
mlIndFix3 <- maxLik( llfInd, start = startValFix, fixed = "mu" )
all.equal( mlIndFix, mlIndFix3, tolerance = 1e-3 )
mlIndFix4 <- maxLik( llfInd, start = startValFix, fixed = 1 )
all.equal( mlIndFix, mlIndFix4, tolerance = 1e-3 )
summary( mlIndFix )
all.equal( mlFix[ ], mlIndFix[ -11 ], tolerance = 1e-3 )
mlFix[[3]]
mlIndFix[[3]]
mlIndFix[ 11 ]
nObs( mlIndFix )

# with analytical gradients
mlgFix <- maxLik( llf, gf, start = startValFix, activePar = !isFixed )
mlgFix1 <- maxLik( llf, gf, start = startValFix, activePar = 2 )
all.equal( mlgFix, mlgFix1, tolerance = 1e-3 )
mlgFix2 <- maxLik( llf, gf, start = startValFix, fixed = isFixed )
all.equal( mlgFix, mlgFix2, tolerance = 1e-3 )
summary( mlgFix )
all.equal( mlFix, mlgFix, tolerance = 1e-3 )
mlFix[[3]]
mlgFix[[3]]
mlFix[[4]]
mlgFix[[4]]
mlgIndFix <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed )
all.equal( mlIndFix, mlgIndFix, tolerance = 1e-3 )
mlIndFix[[3]]
mlgIndFix[[3]]
mlIndFix[[4]]
mlgIndFix[[4]]
all.equal( mlgFix[ ], mlgIndFix[ -11 ], tolerance = 1e-3 )
mlgIndFix[ 11 ]

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
all.equal( mlFix[ -c( 5, 6, 9, 10 ) ], mlFixBHHH[ -c( 5, 6, 9, 10, 11 ) ],
   tolerance = 1e-3 )
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
summary( mlgFixBHHH )
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
all.equal( mlghFix[ -c( 5, 6, 9, 10 ) ], mlFixBfgs[ -c( 5, 6, 9, 10, 11 ) ],
   tolerance = 1e-3 )
mlIndFixBfgs <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "BFGS" )
all.equal( mlFixBfgs[ -9 ], mlIndFixBfgs[ -c(9,12) ], tolerance = 1e-3 )
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
summary( mlgFixBfgs )
all.equal( mlFixBfgs[ -9 ], mlgFixBfgs[ -9 ], tolerance = 1e-3 )
mlgIndFixBfgs <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "BFGS")
all.equal( mlgFixBfgs[ ], mlgIndFixBfgs[ -12 ], tolerance = 1e-3 )
mlgIndFixBfgs[ 12 ]
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
all.equal( mlFixBfgs[ -c( 9, 10 ) ], mlFixNm[ -c( 9, 10 ) ], tolerance = 1e-3 )
mlIndFixNm <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "NM" )
all.equal( mlFixNm[ ], mlIndFixNm[ -12 ], tolerance = 1e-3 )
mlIndFixNm[ 12 ]
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
summary( mlgFixNm )
all.equal( mlFixNm, mlgFixNm, tolerance = 1e-3 )
mlgIndFixNm <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "NM")
all.equal( mlgFixNm[ ], mlgIndFixNm[ -12 ], tolerance = 1e-3 )
mlgIndFixNm[ 12 ]

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
all.equal( mlFixBfgs[ -c( 9, 10 ) ], mlFixSann[ -c( 9, 10 ) ], 
   tolerance = 1e-3 )
mlIndFixSann <- maxLik( llfInd, start = startValFix, fixed = isFixed,
   method = "SANN" )
all.equal( mlFixSann[ ], mlIndFixSann[ -12 ], tolerance = 1e-3 )
mlIndFixSann[ 12 ]
nObs( mlIndFixSann )

# with analytical gradients
mlgFixSann <- maxLik( llf, gf, start = startValFix, fixed = isFixed,
   method = "SANN" )
summary( mlgFixSann )
all.equal( mlFixSann, mlgFixSann, tolerance = 1e-3 )
mlgIndFixSann <- maxLik( llfInd, gfInd, start = startValFix, fixed = isFixed,
   method = "SANN")
all.equal( mlgFixSann[ ], mlgIndFixSann[ -12 ], tolerance = 1e-3 )
mlgIndFixSann[ 12 ]

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
all.equal( mlBfgsInEq[ ], mlBfgsInEqInd[ -12 ], tolerance = 1e-3 )
mlBfgsInEqInd[ 12 ]
nObs( mlBfgsInEqInd )

# with analytical gradients
mlgBfgsInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlBfgsInEq, mlgBfgsInEq, tolerance = 1e-3 )
mlgBfgsInEqInd <- maxLik( llfInd, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd[ -12 ], mlgBfgsInEq[ ], tolerance = 1e-3 )
mlgBfgsInEqInd[ 12 ]
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
all.equal( mlBfgsInEq, mlNmInEq, tolerance = 1e-3 )
mlNmInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "NM" )
summary( mlNmInEqInd )
all.equal( mlNmInEq[ ], mlNmInEqInd[ -12 ], tolerance = 1e-3 )
mlNmInEqInd[ 12 ]
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
all.equal( mlBfgsInEq, mlSannInEq, tolerance = 1e-3 )

# with unused analytical gradients
mlgSannInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "SANN" )
all.equal( mlSannInEq, mlgSannInEq, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSannInEqCand <- maxLik( llf, start = startVal, constraints = inEq,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSannInEqCand )
all.equal( mlSannInEqCand, mlSannInEq, tolerance = 1e-3 )

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
all.equal( mlCon[], mlConInd[-11], tolerance = 1e-3 )
mlConInd[11]
nObs( mlConInd )

# with analytical gradients
mlgCon <- maxLik( llf, gf, start = startVal, constraints = eqCon )
summary( mlgCon )
all.equal( mlCon[ -c( 5, 6, 7, 9 ) ], mlgCon[ -c( 5, 6, 7, 9 ) ], 
   tolerance = 1e-3 )
mlgConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon )
all.equal( mlConInd, mlgConInd, tolerance = 1e-3 )
all.equal( mlgCon[], mlgConInd[-11], tolerance = 1e-3 )
mlgConInd[11]

# with analytical gradients as attribute
mlGCon <- maxLik( llfGrad, start = startVal, constraints = eqCon )
all.equal( mlGCon, mlgCon, tolerance = 1e-3 )
all.equal( mlGCon, mlCon, tolerance = 1e-3 )

# with analytical gradients and Hessians
mlghCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon )
all.equal( mlgCon, mlghCon, tolerance = 1e-3 )

# with analytical gradients and Hessians as attributes
mlGHCon <- maxLik( llfGradHess, start = startVal, constraints = eqCon )
all.equal( mlGHCon, mlghCon, tolerance = 1e-3 )
all.equal( mlGHCon, mlCon, tolerance = 1e-3 )


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
all.equal( mlCon[ -c( 5, 6, 7, 9, 10 ) ], mlBhhhCon[ -c( 5, 6, 7, 9, 10, 11 ) ],
   tolerance = 1e-3 )
mlBhhhCon[11]
nObs( mlBhhhCon )

# with analytical gradients
mlgBhhhCon <- maxLik( llf, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
summary( mlgBhhhCon )
all.equal( mlBhhhCon, mlgBhhhCon, tolerance = 1e-3 )
mlgBhhhConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
all.equal( mlgBhhhCon, mlgBhhhConInd, tolerance = 1e-3 )

# with analytical gradients as attribute
mlGBhhhCon <- maxLik( llfGradInd, start = startVal, constraints = eqCon,
   method = "BHHH" )
summary( mlGBhhhCon )
all.equal( mlGBhhhCon, mlgBhhhCon, tolerance = 1e-3 )
all.equal( mlGBhhhCon, mlBhhhCon, tolerance = 1e-3 )

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
all.equal( mlBfgsCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ],
   tolerance = 1e-3 )
mlBfgsConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlBfgsConInd )
all.equal( mlBfgsCon[], mlBfgsConInd[-12], tolerance = 1e-3 )
mlBfgsConInd[12]
nObs( mlBfgsConInd )

# with analytical gradients
mlgBfgsCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlgBfgsCon )
all.equal( mlBfgsCon, mlgBfgsCon, tolerance = 1e-3 )
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
all.equal( mlNmCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ],
   tolerance = 1e-3 )
mlNmConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "NM", SUMTTol=0)
summary( mlNmConInd )
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
all.equal( mlSannCon[ -c( 5, 6, 9, 10 ) ], mlBfgsCon[ -c( 5, 6, 9, 10 ) ],
   tolerance = 1e-3 )

# with unused analytical gradients
mlgSannCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "SANN", SUMTTol=0)
all.equal( mlSannCon, mlgSannCon, tolerance = 1e-3 )

# with a user-specified function to generate a new candidate point
mlSannConCand <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", cand = function(x)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSannConCand )
all.equal( mlSannConCand, mlSannCon, tolerance = 1e-3 )


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
