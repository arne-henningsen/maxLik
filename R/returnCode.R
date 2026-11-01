### Returns return code of maxim objects
### This signals the reason the iterations ended,
### such as successful convergence, or an error

returnCode <- function(x, ...)
    UseMethod("returnCode")

returnCode.default <- function(x, ...)
    .field(x, "returnCode")

returnCode.maxim <- function(x, ...)
    .returnCode(x)

returnCode.maxLik <- function(x, ...)
    .returnCode(x)
