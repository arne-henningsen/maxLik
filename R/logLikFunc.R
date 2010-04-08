## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, ...) {

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   sum( fnOrig( theta, ... ) )
}
