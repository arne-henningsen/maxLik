### methods for extracting standard errors from the models

std <- function(x, ...)
    ## Extract standard deviations from models (as coefficients)
    UseMethod("std")

std.default <- function(x, ...) {
   if(!is.null(x$std))
       return(x$std)
   s <- sqrt(diag(vcov(x)))
   names(s) <- names(coef(x))
   s
}

std.lm <- function(x, ...)
    sqrt(diag(vcov(x)))

