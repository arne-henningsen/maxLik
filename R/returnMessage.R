
returnMessage <- function(x, ...)
    UseMethod("returnMessage")

returnMessage.default <- function(x, ...)
    .field(x, "returnMessage")

returnMessage.maxim <- function(x, ...)
    .returnMessage(x)

returnMessage.maxLik <- function(x, ...)
    .returnMessage(x)
