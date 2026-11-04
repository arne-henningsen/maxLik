maximType <- function(x)
    UseMethod("maximType")

## The default method reads the component 'maximType' directly, not via
## the internal accessor '.maximType()', which prefers the component
## 'type'.  Foreign objects are conventionally documented with the long
## name, and routing the default method through the accessor would change
## the precedence of 'type' versus 'maximType' for such objects.
maximType.default <- function(x)
    .field(x, "maximType")

maximType.maxim <- function(x)
    .maximType(x)
