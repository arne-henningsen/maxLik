## We wrap the gradeint function in two versions:
## 1) sum over possible individual gradients
## 2) strip possible SUMT parameters and sum thereafter
## The former is for passing to the optimizer from withing sumt
## The latter is necessary for passing '...' to the function
logLikGrad <- function(theta, ...) {
   if(!is.null(grad)) {
      g <- grad(theta, ...)
      if(!is.null(dim(g))) {
         if(nrow(g) > 1) {
            g <- colSums( g )
         }
      }
      names( g ) <- names( start )
      return( g )
   }
   g <- numericGradient(logLikFunc, theta, ...)
   if(!is.null(dim(g))) {
      return(colSums(g))
   } else {
      return(g)
   }
}
logLikGradSumt <- function(theta, ...) {
   if(!is.null(grad)) {
      g <- match.call()
      g[names(formals(sumt))] <- NULL
      g[[1]] <- as.name("grad")
      names(g)[2] <- ""
      g <- eval(g, sys.frame(sys.parent()))
      if(!is.null(dim(g))) {
         if(nrow(g) > 1) {
            g <- colSums( g )
         }
      }
      names( g ) <- names( start )
      return( g )
   }
   g <- match.call()
   g[names(formals(sumt))] <- NULL
   g[[1]] <- as.name("numericGradient")
   names(g)[2] <- "t0"
   g$f <- logLikFunc
   g <- eval(g, sys.frame(sys.parent()))
   if(!is.null(dim(g))) {
      return(colSums(g))
   } else {
      return(g)
   }
}
