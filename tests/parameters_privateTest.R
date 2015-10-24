
### Test battery for various optimization parameters for different optimizers.
### 
library(maxLik)
options(digits = 4)
                           # just to avoid so many differences when comparing these output files
## data to fit a normal distribution
set.seed( 123 )
# generate a variable from normally distributed random numbers
N <- 50
x <- rnorm(N, 1, 2 )

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

# start values
startVal <- c( mu = 0, sigma = 1 )

# 
ml <- maxLik( llf, start = startVal )
print(summary(ml))
## tol
mlTol <- maxLik( llf, start = startVal, tol=1)
print(summary(mlTol))
mlTolC <- maxLik(llf, start=startVal, control=list(tol=1))
print(all.equal(mlTol, mlTolC))
try(ml <- maxLik( llf, start = startVal, tol=-1))
try(ml <- maxLik( llf, start = startVal, tol=c(1,2)))
try(ml <- maxLik( llf, start = startVal, tol=TRUE))
try(ml <- maxLik( llf, start = startVal, control=list(tol=-1)))
try(ml <- maxLik( llf, start = startVal, control=list(tol=c(1,2))))
try(ml <- maxLik( llf, start = startVal, control=list(tol=TRUE)))
## reltol
mlRelTol <- maxLik( llf, start = startVal, reltol=1)
print(summary(mlRelTol))
mlRelTolC <- maxLik(llf, start=startVal, control=list(reltol=1))
print(all.equal(mlRelTol, mlRelTolC))
try(ml <- maxLik( llf, start = startVal, reltol=-1))
try(ml <- maxLik( llf, start = startVal, reltol=c(1,2)))
try(ml <- maxLik( llf, start = startVal, reltol=TRUE))
try(ml <- maxLik( llf, start = startVal, control=list(reltol=-1)))
try(ml <- maxLik( llf, start = startVal, control=list(reltol=c(1,2))))
try(ml <- maxLik( llf, start = startVal, control=list(reltol=TRUE)))
## gradtol
mlGradtol <- maxLik( llf, start = startVal, gradtol=1e-2)
print(summary(mlGradtol))
mlGradtolC <- maxLik(llf, start=startVal, control=list(gradtol=1e-2))
print(all.equal(mlGradtol, mlGradtolC))
try(ml <- maxLik( llf, start = startVal, gradtol=-1))
try(ml <- maxLik( llf, start = startVal, gradtol=c(1,2)))
try(ml <- maxLik( llf, start = startVal, gradtol=TRUE))
try(ml <- maxLik( llf, start = startVal, control=list(gradtol=-1)))
try(ml <- maxLik( llf, start = startVal, control=list(gradtol=c(1,2))))
try(ml <- maxLik( llf, start = startVal, control=list(gradtol=TRUE)))
## examples with steptol, lambdatol
## qac
mlMarq <- maxLik( llf, start = startVal, qac="marquardt")
print(summary(mlMarq))
mlMarqC <- maxLik(llf, start=startVal, control=list(qac="marquardt"))
print(all.equal(mlMarq, mlMarqC))
try(ml <- maxLik( llf, start = startVal, qac=-1))
try(ml <- maxLik( llf, start = startVal, qac=c("a", "b")))
try(ml <- maxLik( llf, start = startVal, qac=TRUE))
try(ml <- maxLik( llf, start = startVal, control=list(qac=-1)))
try(ml <- maxLik( llf, start = startVal, control=list(qac=c("a", "b"))))
try(ml <- maxLik( llf, start = startVal, control=list(qac=TRUE)))
mlMarqCl <- a <- maxLik(llf, start = startVal,
                        control=list(qac="marquardt", lambda0=1000, lambdaStep=4))
print(all.equal(coef(mlMarqCl), coef(mlMarq)))
## NM: alpha, beta, gamma
mlNM <- maxLik( llf, start = startVal, method="nm")
print(summary(mlNM))
mlNMAlpha <- maxLik(llf, start=startVal, method="nm", beta=0.8)
mlNMAlphaC <- maxLik(llf, start=startVal, method="nm", control=list(beta=0.8))
print(all.equal(mlNMAlpha, mlNMAlphaC))

## likelihood function with additional parameter
llf1 <- function( param, sigma ) {
   mu <- param
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
       0.5 * sum( ( x - mu )^2 / sigma^2 )
   return( llValue )
       }

## log-lik mixture
logLikMix <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll
}

## loglik mixture with additional parameter
logLikMixA <- function(param, rho) {
   mu1 <- param[1]
   mu2 <- param[2]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
   ll
}

## Test the following with all the main optimizers:
for(method in c("NR", "BFGS", "BFGSR")) {
   ## two parameters at the same time
   ## iterlim, printLevel
   cat("-- method", method, "--\n")
   N <- 100
   x <- rnorm(N, 1, 2 )
   startVal <- c(1,2)
   ml2 <- maxLik( llf, start=startVal, method=method, iterlim=1, printLevel=2)
   print(summary(ml2))
   ml2C <- maxLik(llf, start=startVal, method=method,
                  control=list(iterlim=1, printLevel=2))
   print(all.equal(ml2, ml2C))
   ## what about additional parameters for the loglik function?
   mls <- maxLik(llf1, start=0, method=method, sigma=1)
   print(coef(mls))
   mlsM <- maxLik(llf1, start=0, method=method, tol=1, sigma=1)
   mlsCM <- maxLik(llf1, start=0, method=method, control=list(tol=1), sigma=1)
   cat("Additional parameters to loglik: open == control()?\n")
   print(all.equal(mlsM, mlsCM))
   ## And what about unused parameters?
   cat("What about unused parameters?\n")
   try(maxLik(llf1, start=0, method=method, control=list(tol=1),
              sigma=1, unusedPar=2))
                           # error
   N <- 100
   ## Does this work with constraints?
   x <- c(rnorm(N, mean=-1), rnorm(N, mean=1))
   ## First test inequality constraints
   ## Inequality constraints: x + y + z < 0.5
   A <- matrix(c(-1, 0, 0,
                 0, -1, 0,
                 0, 0, 1), 3, 3, byrow=TRUE)
   B <- rep(0.5, 3)
   start <- c(0.4, 0, 0.9)
   ## analytic gradient
   cat("Inequality constraints, analytic gradient & Hessian\n")
   mix <- try(maxLik(logLikMix, 
                     start=start, method=method,
                     constraints=list(ineqA=A, ineqB=B)))
   if(!inherits(mix, "try-error")) {
      print(summary(mix))
   }
   mixGT <- try(maxLik(logLikMix, 
                       start=start, method=method,
                       constraints=list(ineqA=A, ineqB=B),
                       tol=1))
   if(!inherits(mixGT, "try-error")) {
      print(summary(mixGT))
   }
   mixGTC <- try(maxLik(logLikMix, 
                    start=start, method=method,
                    constraints=list(ineqA=A, ineqB=B),
                    control=list(tol=1)))
   if(!inherits(mixGTC, "try-error")) {
      print(all.equal(mixGT, mixGTC))
   }
   ## 2d inequality constraints: x + y < 0.5
   A2 <- matrix(c(-1, -1), 1, 2, byrow=TRUE)
   B2 <- 0.5
   start2 <- c(-0.5, 0.5)
   cat("Inequality constraints, additional parameters\n")
   mixA <- try(maxLik(logLikMixA, 
                      start=start2, method=method,
                      constraints=list(ineqA=A2, ineqB=B2),
                      tol=1,
                      rho=0.5))
   mixAC <- try(maxLik(logLikMixA, 
                       start=start2, method=method,
                       constraints=list(ineqA=A2, ineqB=B2),
                       control=list(tol=1),
                       rho=0.5))
   if(!inherits(mixA, "try-error") & !inherits(mixAC, "try-error")) {
      cat("Coefficients equal?\n")
      print(all.equal(coef(mixA), coef(mixAC)))
      cat("Hessians equal?\n")
      print(all.equal(hessian(mixA), hessian(mixAC)))
   }
}

### Test adding both default and user-specified parameters through control list
estimate <- function(control=NULL, ...) {
   return(maxLik(llf, start=c(1,1),
                 control=c(list(iterlim=100), control),
                 ...))
}
m <- estimate(control=list(iterlim=1), fixed=2)
show(maxControl(m))
                           # iterlim should be 1
print(coef(m))
                           # sigma should be 1.000
## Does print.level overwrite 'printLevel'?
m <- estimate(control=list(printLevel=2, print.level=1))
show(maxControl(m))

## Does open parameters override everything?
m <- estimate(control=list(printLevel=2, print.level=1), print.level=0)
show(maxControl(m))

### does both printLevel, print.level work for condiNumber?
condiNumber(hessian(m), print.level=0) # no output
condiNumber(hessian(m), printLevel=0)  # no output
condiNumber(hessian(m), printLevel=0, print.level=1) # no output

