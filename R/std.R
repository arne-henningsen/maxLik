### methods for extracting standard errors from the models

stdEr <- function(x, ...)
    ## Extract standard deviations from models (as coefficients)
    UseMethod("stdEr")

stdEr.default <- function(x, ...) {
   if(!is.null(x$std))
       return(x$std)
   s <- sqrt(diag(vcov(x)))
   names(s) <- names(coef(x))
   s
}

stdEr.lm <- function(x, ...)
    sqrt(diag(vcov(x)))

