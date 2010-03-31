## We wrap the gradeint function in two versions:
## 1) sum over possible individual gradients
## 2) strip possible SUMT parameters and sum thereafter
## The former is for passing to the optimizer from withing sumt
## The latter is necessary for passing '...' to the function
logLikGrad <- function(theta, ...) {
   if(is.null(grad)) {
      g <- numericGradient(logLikFunc, theta, ...)
   } else {
      g <- grad(theta, ...)
   }
   if(!is.null(dim(g))) {
      g <- colSums(g)
   }
   names( g ) <- names( start )
   return( g )
}
logLikGradSumt <- function(theta, ...) {
   g <- match.call()
   g[names(formals(sumt))] <- NULL
   g[[1]] <- as.name("logLikGrad")
   names(g)[2] <- ""
   g <- eval(g, sys.frame(sys.parent()))
   return( g )
}
