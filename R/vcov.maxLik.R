### Return variance-covariance matrix.
### If Hessian does not exist, return matrix of NA-s
### 
vcov.maxLik <- function(object, eigentol=1e-12, ...) {
   return( .vcov( object, eigentol=eigentol ) )
}
