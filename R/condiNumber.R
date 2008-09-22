### condiNumber: print matrix' condition number adding columns one by one.
### In this way user may investigate the which columns cause problems with singularity

condiNumber <- function(x, ...)
    UseMethod("condiNumber")

condiNumber.default <- function(x, exact=FALSE, norm=FALSE, print.level=1, ...) {
   ## x:  a matrix, condition number of which are to be printed
   ## exact: whether the condition number have to be exact or approximated (see 'kappa')
   ## norm: whether to normalise the matrix' columns.
   ## print.level: whether to print the condition numbers while calculating.  Useful for interactive testing. 
   if(dim(x)[2] > dim(x)[1]) {
      warning(paste(dim(x)[1], "rows and", dim(x)[2], "columns, use transposed matrix"))
      x <- t(x)
   }
   cn <- numeric(ncol(x))
    if(norm) {
                                        # Now normalise column vectors
       x <- apply(x, 2, FUN=function(v) v/sqrt(sum(v*v)))
    }
    for(i in seq(length=ncol(x))) {
        m <- x[,1:i]
        cn[i] <- kappa(m, exact=exact)
        if(print.level > 0)
            cat(colnames(x)[i], "\t", cn[i], "\n")
    }
   names(cn) <- colnames(x)
   invisible(cn)
}

condiNumber.maxLik <- function(x, ...)
    condiNumber.default(hessian(x, ...))
