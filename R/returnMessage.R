
returnMessage <- function(x, ...)
    UseMethod("returnMessage")

## The default method reads the component 'returnMessage' directly, not
## via the internal accessor '.returnMessage()', which prefers the
## component 'message'.  Foreign objects are conventionally documented
## with the long name, and routing the default method through the accessor
## would change the precedence of 'message' versus 'returnMessage' for
## such objects.
returnMessage.default <- function(x, ...)
    .field(x, "returnMessage")

returnMessage.maxim <- function(x, ...)
    .returnMessage(x)

returnMessage.maxLik <- function(x, ...)
    .returnMessage(x)
