## strip possible SUMT parameters and call the function thereafter
callWithoutSumt <- function(theta, fName, ...) {
   f <- match.call()
   f[names(formals(sumt))] <- NULL
   f[[1]] <- as.name(fName)
   names(f)[2] <- ""
   f[[3]] <- NULL
   f1 <- eval(f, sys.frame(sys.parent()))
   return( f1 )
}
