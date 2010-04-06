## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, func, ...) {
   sum( func( theta, ... ) )
}
