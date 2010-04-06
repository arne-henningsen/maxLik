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
mlgSANN <- maxLik( llf, gf, start = startVal, method = "SANN" )
summary( mlgSANN )
all.equal( mlSANN, mlgSANN )
# with individual log likelihood values and gradients
mlgIndSANN <- maxLik( llfInd, gfInd, start = startVal, method = "SANN" )
summary( mlgIndSANN )
all.equal( mlgSANN, mlgIndSANN )

# with unused analytical gradients and Hessian
mlghSANN <- maxLik( llf, gf, hf, start = startVal, method = "SANN" )
all.equal( mlgSANN, mlghSANN )

# with a user-specified function to generate a new candidate point
mlSANNCand <- maxLik( llf, start = startVal, method = "SANN",
   cand = function(x,fnOrig=NULL)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSANNCand )
all.equal( mlSANNCand, mlSANN )

############### with fixed parameters ###############
# start values
startValFix <- c( mu = 1, sigma = 1 )

# fix mu (the mean ) at its start value
isFixed <- c( TRUE, FALSE )

## NR method
mlFix <- maxLik( llf, start = startValFix, activePar = !isFixed )
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
nObs( mlFix )
nParam( mlFix )
returnCode( mlFix )
returnMessage( mlFix )
vcov( mlFix )
logLik( summary( mlFix ) )
mlIndFix <- maxLik( llfInd, start = startValFix, activePar = !isFixed )
summary( mlIndFix )
all.equal( mlFix, mlIndFix )
mlFix[[3]]
mlIndFix[[3]]

# with analytical gradients
mlgFix <- maxLik( llf, gf, start = startValFix, activePar = !isFixed )
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
all.equal( mlgFix, mlgIndFix )

# with analytical gradients and Hessians
mlghFix <- maxLik( llf, gf, hf, start = startValFix, activePar = !isFixed )
all.equal( mlgFix, mlghFix )
mlgFix[[4]]
mlghFix[[4]]

## BHHH method
mlFixBHHH <- maxLik( llfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
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
nObs( mlFixBHHH )
nParam( mlFixBHHH )
returnCode( mlFixBHHH )
returnMessage( mlFixBHHH )
vcov( mlFixBHHH )
logLik( summary( mlFixBHHH ) )
all.equal( mlFix[ -c( 5, 6, 9, 10 ) ], mlFixBHHH[ -c( 5, 6, 9, 10 ) ] )
mlFix[[ 3 ]]
mlFixBHHH[[ 3 ]]
mlFix[[ 4 ]]
mlFixBHHH[[ 4 ]]

# with analytical gradients
mlgFixBHHH <- maxLik( llfInd, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
summary( mlgFixBHHH )
all.equal( mlFixBHHH, mlgFixBHHH )
mlgFixBHHH2 <- maxLik( llf, gfInd, start = startValFix, activePar = !isFixed,
   method = "BHHH")
all.equal( mlgFixBHHH, mlgFixBHHH2 )

# with unused Hessians
mlghFixBHHH <- maxLik( llfInd, gfInd, hf, start = startValFix, activePar = !isFixed,
   method = "BHHH" )
all.equal( mlgFixBHHH, mlghFixBHHH )


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
nObs( mlBfgsInEq )
nParam( mlBfgsInEq )
returnCode( mlBfgsInEq )
returnMessage( mlBfgsInEq )
vcov( mlBfgsInEq )
logLik( summary( mlBfgsInEq ) )
mlBfgsInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
summary( mlBfgsInEqInd )
all.equal( mlBfgsInEq, mlBfgsInEqInd )

# with analytical gradients
mlgBfgsInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlBfgsInEq, mlgBfgsInEq )
mlgBfgsInEqInd <- maxLik( llfInd, gfInd, start = startVal, constraints = inEq,
   method = "BFGS" )
all.equal( mlgBfgsInEqInd, mlgBfgsInEq )
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
nObs( mlNmInEq )
nParam( mlNmInEq )
returnCode( mlNmInEq )
returnMessage( mlNmInEq )
vcov( mlNmInEq )
logLik( summary( mlNmInEq ) )
all.equal( mlBfgsInEq, mlNmInEq )
mlNmInEqInd <- maxLik( llfInd, start = startVal, constraints = inEq,
   method = "NM" )
summary( mlNmInEqInd )
all.equal( mlNmInEq, mlNmInEqInd )

# with unused analytical gradients
mlgNmInEq <- maxLik( llf, gf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlNmInEq, mlgNmInEq )

# with unused analytical gradients and Hessians
mlghNmInEq <- maxLik( llf, gf, hf, start = startVal, constraints = inEq,
   method = "NM" )
all.equal( mlgNmInEq, mlghNmInEq )

## SANN method with inequality constraints
try( maxLik( llf, start = startVal, constraints = inEq, method = "SANN" ) )


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
nObs( mlCon )
nParam( mlCon )
returnCode( mlCon )
returnMessage( mlCon )
vcov( mlCon )
logLik( summary( mlCon ) )
# mlConInd <- maxLik( llfInd, start = startVal, constraints = eqCon )
   ##### this optimization seems to end up in an infinite loop
# summary( mlConInd )
# all.equal( mlCon, mlConInd )

# with analytical gradients
mlgCon <- maxLik( llf, gf, start = startVal, constraints = eqCon )
summary( mlgCon )
all.equal( mlCon[ -c( 5, 6, 7, 9 ) ], mlgCon[ -c( 5, 6, 7, 9 ) ] )
mlgConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon )
# all.equal( mlConInd, mlgConInd )
all.equal( mlgCon, mlgConInd )

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
nObs( mlBhhhCon )
nParam( mlBhhhCon )
returnCode( mlBhhhCon )
returnMessage( mlBhhhCon )
vcov( mlBhhhCon )
logLik( summary( mlBhhhCon ) )
all.equal( mlCon[ -c( 5, 6, 7, 9, 10 ) ], mlBhhhCon[ -c( 5, 6, 7, 9, 10 ) ] )

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
nObs( mlBfgsCon )
nParam( mlBfgsCon )
returnCode( mlBfgsCon )
returnMessage( mlBfgsCon )
vcov( mlBfgsCon )
logLik( summary( mlBfgsCon ) )
all.equal( mlBfgsCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ] )
mlBfgsConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlBfgsConInd )
all.equal( mlBfgsCon, mlBfgsConInd )

# with analytical gradients
mlgBfgsCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "BFGS" )
summary( mlgBfgsCon )
all.equal( mlBfgsCon, mlgBfgsCon )
mlgBfgsConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon, mlgBfgsConInd )

# with analytical gradients and unused Hessians
mlghBfgsCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "BFGS" )
all.equal( mlgBfgsCon, mlghBfgsCon )

## NM method with equality constraints
mlNmCon <- maxLik( llf, start = startVal, constraints = eqCon, method = "NM" )
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
nObs( mlNmCon )
nParam( mlNmCon )
returnCode( mlNmCon )
returnMessage( mlNmCon )
vcov( mlNmCon )
logLik( summary( mlNmCon ) )
all.equal( mlNmCon[ -c( 5, 6, 9, 10 ) ], mlCon[ -c( 5, 6, 9, 10 ) ] )
mlNmConInd <- maxLik( llfInd, start = startVal, constraints = eqCon,
   method = "NM" )
summary( mlNmConInd )
all.equal( mlNmCon, mlNmConInd )

# with unused analytical gradients
mlgNmCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "NM" )
all.equal( mlNmCon, mlgNmCon )
mlgNmConInd <- maxLik( llfInd, gfInd, start = startVal, constraints = eqCon,
   method = "NM" )
all.equal( mlgNmCon, mlgNmConInd )

# with unused analytical gradients and Hessians
mlghNmCon <- maxLik( llf, gf, hf, start = startVal, constraints = eqCon,
   method = "NM" )
all.equal( mlgNmCon, mlghNmCon )

## SANN method with equality constraints
mlSannCon <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN" )
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
nObs( mlSannCon )
nParam( mlSannCon )
returnCode( mlSannCon )
returnMessage( mlSannCon )
vcov( mlSannCon )
logLik( summary( mlSannCon ) )
all.equal( mlSannCon[ -c( 5, 6, 9, 10 ) ], mlBfgsCon[ -c( 5, 6, 9, 10 ) ] )

# with unused analytical gradients
mlgSannCon <- maxLik( llf, gf, start = startVal, constraints = eqCon,
   method = "SANN" )
all.equal( mlSannCon, mlgSannCon )

# with a user-specified function to generate a new candidate point
mlSannConCand <- maxLik( llf, start = startVal, constraints = eqCon,
   method = "SANN", cand = function(x,fnOrig=NULL)c(rnorm(1,x[1]),rnorm(1,x[2])) )
summary( mlSannConCand )
all.equal( mlSannConCand, mlSannCon )

