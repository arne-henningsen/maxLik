## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, gradOrig, hessOrig, ...) {
   sum( fnOrig( theta, ... ) )
}
