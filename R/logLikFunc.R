## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, ...) {
   sum( fnOrig( theta, ... ) )
}
