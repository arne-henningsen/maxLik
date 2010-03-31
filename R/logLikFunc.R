## We wrap the objective function in two versions:
## 1) sum over possible individual likelihoods
## 2) strip possible SUMT parameters and sum thereafter
## The former is for passing to the optimizer from within sumt
## The latter is necessary for passing '...' to the function
logLikFunc <- function(theta, ...) {
   sum( fn( theta, ... ) )
}
logLikFuncSumt <- function(theta, ...) {
   ## this wrapper makes a) single-valued function (in case of BHHH
   ## vector-valued); and b) strips the SUMT extra arguments
   f <- match.call()
   f[names(formals(sumt))] <- NULL
   f[[1]] <- as.name("logLikFunc")
   names(f)[2] <- ""
   f1 <- eval(f, sys.frame(sys.parent()))
   return( f1 )
}
