maximType <- function(x)
    UseMethod("maximType")

maximType.default <- function(x)
    .field(x, "maximType")

maximType.maxim <- function(x)
    .maximType(x)
