### Return #of observations for models
### using the standard generic 'nobs'
nobs.maxLik <- function(object, ...) {
   .nObs( object )
}

setMethod("nobs", "maxLik", nobs.maxLik)


## Return #of observations for models using our own
## generic 'nObs' (deprecated)
nObs.maxLik <- function(x, ...) {
   .nObs( x )
}
