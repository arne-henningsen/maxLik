## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, ...) {
   sum( fn( theta, ... ) )
}
