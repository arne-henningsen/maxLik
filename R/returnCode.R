### Returns return code of maxim objects
### This signals the reason the iterations ended,
### such as successful convergence, or an error

returnCode <- function(x, ...)
    UseMethod("returnCode")

## The default method reads the component 'returnCode' directly, not via
## the internal accessor '.returnCode()', which prefers the component
## 'code'.  Foreign objects are conventionally documented with the long
## name, and routing the default method through the accessor would change
## the precedence of 'code' versus 'returnCode' for such objects.
returnCode.default <- function(x, ...)
    .field(x, "returnCode")

returnCode.maxim <- function(x, ...)
    .returnCode(x)

returnCode.maxLik <- function(x, ...)
    .returnCode(x)
