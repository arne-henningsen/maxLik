
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

logLikMix <- function(param) {
   rho <- param[1]
   if(rho < 0 || rho > 1)
       return(NA)
   mu1 <- param[2]
   mu2 <- param[3]
   ll <- log(rho*dnorm(x - mu1) + (1 - rho)*dnorm(x - mu2))
#   ll <- sum(ll)
   ll
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
## NM: alpha, beta, gamma
mlNM <- maxLik( llf, start = startVal, method="nm")
print(summary(mlNM))
mlNMAlpha <- maxLik(llf, start=startVal, method="nm", beta=0.8)
mlNMAlphaC <- maxLik(llf, start=startVal, method="nm", control=list(beta=0.8))
print(all.equal(mlNMAlpha, mlNMAlphaC))

## two parameters at the same time
## iterlim, printLevel
ml2 <- maxLik( llf, start = startVal, method="nm", iterlim=1, printLevel=2)
print(summary(ml2))
ml2C <- maxLik(llf, start=startVal, control=list(iterlim=1, printLevel=2))
print(all.equal(ml2, ml2C))

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
mix <- maxLik(logLikMix, 
              start=start,
              constraints=list(ineqA=A, ineqB=B))
print(summary(mix))
mixGT <- maxLik(logLikMix, 
              start=start,
              constraints=list(ineqA=A, ineqB=B),
                gradtol=1e-2)
print(summary(mixGT))
mixGTC <- maxLik(logLikMix, 
                 start=start,
                 constraints=list(ineqA=A, ineqB=B),
                 control=list(gradtol=1e-2))
print(all.equal(mixGT, mixGTC))


